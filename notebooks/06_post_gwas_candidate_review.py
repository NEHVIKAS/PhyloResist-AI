"""Create a reproducible post-GWAS candidate-locus and LD review.

Run from any working directory with the project's ``phyloresist`` environment:

    python notebooks/06_post_gwas_candidate_review.py

Inputs are read from ``data/raw``, ``data/reference``, and the existing
mrMLM final-result files. The script writes derived reports under ``results``;
it never edits the raw CSVs, the mrMLM result files, or earlier scripts.

The annotation uses the paper's default +/-3 kb window. The small LD review
checks pairwise r-squared around the nearest reported QTN for each gene named
in the paper. It does not reproduce the authors' manual LD-boundary decisions
and must not be interpreted as proof that a nearby gene is causal.
"""

from __future__ import annotations

import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path, PurePosixPath
from typing import Any
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
REFERENCE_DIR = PROJECT_ROOT / "data" / "reference"
GWAS_DIR = PROJECT_ROOT / "results" / "gwas" / "mrmlm"
REVIEW_DIR = GWAS_DIR / "review"
REPORT_DIR = PROJECT_ROOT / "results" / "reports"

GENOTYPE_PATH = RAW_DIR / "genotype.csv"
ANNOTATION_PATH = REFERENCE_DIR / "full_annotation.xlsx"
PHENOTYPE_PATHS = {
    "propiconazole": RAW_DIR / "phenotype_propiconazole.csv",
    "tebuconazole": RAW_DIR / "phenotype_tebuconazole.csv",
}
FINAL_RESULT_PATHS = {
    fungicide: GWAS_DIR / fungicide / "1_Final result.csv"
    for fungicide in PHENOTYPE_PATHS
}

DEFAULT_WINDOW_BP = 3_000
LD_SCAN_RADIUS_BP = 100_000
HIGH_LD_R2 = 0.75
EXCEL_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"

# Genes explicitly discussed in the supplied study manuscript. This list is
# only used to focus the review; it does not define which loci are significant.
PAPER_CANDIDATES = {
    "FGRAMPH1_01G09621": "Multidrug-resistance ABC transporter; paper reports tebuconazole only.",
    "FGRAMPH1_01G13811": "Quinate transporter discussed by the paper.",
    "FGRAMPH1_01G16125": "MFS transporter discussed by the paper.",
    "FGRAMPH1_01G19393": "CTR copper transporter discussed by the paper.",
    "FGRAMPH1_01G22013": "Transport-related protein discussed by the paper.",
    "FGRAMPH1_01G25701": "Fusarin C cluster transporter discussed by the paper.",
    "FGRAMPH1_01G08375": "Nonribosomal peptide synthetase discussed by the paper.",
    "FGRAMPH1_01G20955": "Nonribosomal peptide synthetase discussed by the paper.",
    "FGRAMPH1_01G08383": "Polyketide synthase discussed by the paper.",
    "FGRAMPH1_01G21231": "HSP70 candidate discussed for both fungicides.",
}


def require_inputs() -> None:
    """Stop early with a clear list if a required input is absent."""
    paths = [GENOTYPE_PATH, ANNOTATION_PATH, *PHENOTYPE_PATHS.values(), *FINAL_RESULT_PATHS.values()]
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input file(s):\n" + "\n".join(str(p) for p in missing))


def normalize_scaffold(value: Any) -> str:
    """Map marker labels such as SHG970334.1 to workbook labels HG970334."""
    text = str(value).strip()
    if text.startswith("S"):
        text = text[1:]
    return re.sub(r"\.\d+$", "", text)


def parse_marker(marker: Any) -> tuple[str, int]:
    """Return normalized scaffold and coordinate from a GWAS marker ID."""
    text = str(marker).strip()
    try:
        prefix, position = text.rsplit("_", 1)
        scaffold = normalize_scaffold(prefix)
        coordinate = int(position)
    except (ValueError, TypeError) as error:
        raise ValueError(f"Could not parse marker ID {marker!r}; expected scaffold_position form") from error
    if not scaffold or coordinate < 1:
        raise ValueError(f"Invalid marker ID {marker!r}")
    return scaffold, coordinate


def normalize_id(value: Any) -> str:
    """Normalize isolate IDs without dropping meaningful leading zeros."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return re.sub(r"\.0+$", "", text)


def find_isolate_id_column(frame: pd.DataFrame) -> str:
    """Identify the phenotype ID field despite label capitalization/format."""
    labels = list(frame.columns)
    explicit = [
        column for column in labels
        if re.search(r"(^|[^a-z])(id|isolate|sample|strain)([^a-z]|$)", str(column).strip().casefold())
    ]
    if explicit:
        return explicit[0]
    phenotype_labels = [column for column in labels if "phenotype" in str(column).casefold()]
    return (phenotype_labels or labels[:1])[0]


def validate_isolate_ids(genotype: pd.DataFrame) -> list[str]:
    """Check that each phenotype file contains exactly the genotype isolates."""
    if genotype.shape[1] <= 3:
        raise ValueError("Expected three genotype metadata columns followed by isolate columns.")
    isolate_ids = [normalize_id(value) for value in genotype.columns[3:]]
    if not isolate_ids or any(not value for value in isolate_ids):
        raise ValueError("Genotype isolate headers contain blank IDs.")
    if len(isolate_ids) != len(set(isolate_ids)):
        raise ValueError("Genotype isolate headers contain duplicate IDs.")

    for fungicide, path in PHENOTYPE_PATHS.items():
        phenotype = pd.read_csv(path, dtype=str, keep_default_na=True)
        id_column = find_isolate_id_column(phenotype)
        phenotype_ids = [normalize_id(value) for value in phenotype[id_column]]
        if not phenotype_ids or any(not value for value in phenotype_ids):
            raise ValueError(f"Blank phenotype isolate ID in {path.name} (column {id_column!r}).")
        if len(phenotype_ids) != len(set(phenotype_ids)):
            raise ValueError(f"Duplicate phenotype isolate IDs in {path.name}.")
        if set(phenotype_ids) != set(isolate_ids):
            only_genotype = sorted(set(isolate_ids) - set(phenotype_ids))
            only_phenotype = sorted(set(phenotype_ids) - set(isolate_ids))
            raise ValueError(
                f"Isolate IDs do not match for {fungicide}; "
                f"genotype-only={only_genotype[:10]}, phenotype-only={only_phenotype[:10]}"
            )
        trait_columns = [column for column in phenotype.columns if str(column).strip().casefold() == "trait1"]
        if len(trait_columns) != 1:
            raise ValueError(f"Expected exactly one trait1 field in {path.name}.")
        trait_values = pd.to_numeric(phenotype[trait_columns[0]], errors="coerce")
        if trait_values.isna().any() or not np.isfinite(trait_values.to_numpy(dtype=float)).all():
            raise ValueError(f"trait1 contains missing, non-numeric, or non-finite values in {path.name}.")
    return isolate_ids


def _excel_column_index(cell_reference: str) -> int:
    """Convert an Excel cell reference (for example, A1 or AC4) to zero-based column."""
    letters = re.match(r"[A-Z]+", cell_reference.upper())
    if not letters:
        raise ValueError(f"Invalid Excel cell reference: {cell_reference!r}")
    index = 0
    for letter in letters.group(0):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def read_annotation_workbook(path: Path) -> pd.DataFrame:
    """Read the source workbook with Python's standard library.

    This avoids requiring an additional Excel-reading package in the user's
    conda environment. It reads the first worksheet named ``full annotation``.
    """
    with zipfile.ZipFile(path) as workbook:
        workbook_root = ET.fromstring(workbook.read("xl/workbook.xml"))
        relationship_root = ET.fromstring(workbook.read("xl/_rels/workbook.xml.rels"))
        relationships = {
            relation.attrib["Id"]: relation.attrib["Target"]
            for relation in relationship_root.findall(f"{{{PKG_REL_NS}}}Relationship")
        }

        sheet_path = None
        for sheet in workbook_root.findall(f".//{{{EXCEL_NS}}}sheet"):
            if sheet.attrib.get("name", "").strip().casefold() == "full annotation":
                relationship_id = sheet.attrib.get(f"{{{REL_NS}}}id")
                target = relationships.get(relationship_id or "")
                if target:
                    target_path = PurePosixPath(target)
                    sheet_path = str(target_path).lstrip("/") if target.startswith("/") else str(PurePosixPath("xl") / target_path)
                break
        if not sheet_path or sheet_path not in workbook.namelist():
            raise ValueError("Could not find the 'full annotation' worksheet in the reference workbook.")

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in workbook.namelist():
            shared_root = ET.fromstring(workbook.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(node.text or "" for node in item.iter(f"{{{EXCEL_NS}}}t"))
                for item in shared_root.findall(f"{{{EXCEL_NS}}}si")
            ]

        sheet_root = ET.fromstring(workbook.read(sheet_path))
        sheet_data = sheet_root.find(f"{{{EXCEL_NS}}}sheetData")
        if sheet_data is None:
            raise ValueError("The annotation worksheet does not contain sheet data.")

        rows: list[list[str]] = []
        for row_element in sheet_data.findall(f"{{{EXCEL_NS}}}row"):
            cells: dict[int, str] = {}
            for cell in row_element.findall(f"{{{EXCEL_NS}}}c"):
                column_index = _excel_column_index(cell.attrib.get("r", ""))
                cell_type = cell.attrib.get("t", "")
                if cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{EXCEL_NS}}}t"))
                else:
                    value_node = cell.find(f"{{{EXCEL_NS}}}v")
                    value = value_node.text if value_node is not None and value_node.text else ""
                    if cell_type == "s" and value:
                        value = shared_strings[int(value)]
                cells[column_index] = value.strip()
            if cells:
                row_width = max(cells) + 1
                values = [""] * row_width
                for index, value in cells.items():
                    values[index] = value
                rows.append(values)

    if not rows:
        raise ValueError("The annotation worksheet is empty.")
    header = [str(value).strip() for value in rows[0]]
    records = [row + [""] * max(0, len(header) - len(row)) for row in rows[1:]]
    annotation = pd.DataFrame(records, columns=header)
    required_columns = {"Gene ID", "Product Description", "Chr", "start", "stop"}
    missing = required_columns - set(annotation.columns)
    if missing:
        raise ValueError(f"Annotation worksheet is missing columns: {sorted(missing)}")

    annotation = annotation.rename(
        columns={
            "Gene ID": "gene_id",
            "Product Description": "product_description",
            "Chr": "scaffold",
            "start": "gene_start",
            "stop": "gene_stop",
        }
    )
    annotation["scaffold"] = annotation["scaffold"].map(normalize_scaffold)
    annotation["gene_start"] = pd.to_numeric(annotation["gene_start"], errors="coerce")
    annotation["gene_stop"] = pd.to_numeric(annotation["gene_stop"], errors="coerce")
    annotation = annotation.dropna(subset=["gene_id", "scaffold", "gene_start", "gene_stop"])
    annotation["gene_start"] = annotation["gene_start"].astype(int)
    annotation["gene_stop"] = annotation["gene_stop"].astype(int)
    annotation = annotation.loc[annotation["gene_id"].astype(str).str.strip() != ""].copy()
    return annotation.reset_index(drop=True)


def normalized_column_lookup(frame: pd.DataFrame) -> dict[str, str]:
    """Map case/space-normalized names to the actual result CSV headers."""
    return {re.sub(r"\s+", " ", str(column).strip()).casefold(): str(column) for column in frame.columns}


def require_column(lookup: dict[str, str], requested: str, file_path: Path) -> str:
    key = re.sub(r"\s+", " ", requested.strip()).casefold()
    if key not in lookup:
        raise ValueError(f"Required column {requested!r} not found in {file_path.name}; found {list(lookup.values())}")
    return lookup[key]


def read_gwas_results() -> tuple[dict[str, pd.DataFrame], dict[str, dict[str, dict[str, Any]]]]:
    """Read final results and aggregate multiple method calls per marker."""
    by_fungicide: dict[str, pd.DataFrame] = {}
    loci_by_fungicide: dict[str, dict[str, dict[str, Any]]] = {}
    required = ["RS#", "Method", "LOD score", "Chromosome", "Marker position (bp)", "QTN effect", "r2 (%)", "MAF"]

    for fungicide, path in FINAL_RESULT_PATHS.items():
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        lookup = normalized_column_lookup(frame)
        cols = {name: require_column(lookup, name, path) for name in required}
        parsed: list[dict[str, Any]] = []
        for _, row in frame.iterrows():
            marker = str(row[cols["RS#"]]).strip()
            scaffold, position = parse_marker(marker)
            parsed.append(
                {
                    "fungicide": fungicide,
                    "marker": marker,
                    "scaffold": scaffold,
                    "position": position,
                    "chromosome": str(row[cols["Chromosome"]]).strip(),
                    "method": str(row[cols["Method"]]).strip(),
                    "lod": float(row[cols["LOD score"]]),
                    "effect": str(row[cols["QTN effect"]]).strip(),
                    "r2_percent": str(row[cols["r2 (%)"]]).strip(),
                    "maf": str(row[cols["MAF"]]).strip(),
                }
            )
        parsed_frame = pd.DataFrame(parsed)
        by_fungicide[fungicide] = parsed_frame

        loci: dict[str, dict[str, Any]] = {}
        for marker, group in parsed_frame.groupby("marker", sort=False):
            top_row = group.sort_values("lod", ascending=False).iloc[0]
            loci[str(marker)] = {
                "fungicide": fungicide,
                "marker": str(marker),
                "scaffold": str(top_row["scaffold"]),
                "position": int(top_row["position"]),
                "chromosome": str(top_row["chromosome"]),
                "methods": sorted(set(group["method"].astype(str))),
                "method_call_count": int(len(group)),
                "top_lod": float(top_row["lod"]),
                "effect_at_top_lod": str(top_row["effect"]),
                "r2_percent_at_top_lod": str(top_row["r2_percent"]),
                "maf_at_top_lod": str(top_row["maf"]),
            }
        loci_by_fungicide[fungicide] = loci

    return by_fungicide, loci_by_fungicide


def validate_markers_against_genotype(
    genotype: pd.DataFrame, loci_by_fungicide: dict[str, dict[str, dict[str, Any]]]
) -> None:
    """Ensure reported marker IDs and coordinates exist in the raw genotype data."""
    marker_column = genotype.columns[0]
    position_column = next((column for column in genotype.columns[:3] if str(column).strip().casefold() == "pos"), None)
    if position_column is None:
        raise ValueError(f"Expected a 'pos' column among the first three genotype columns; found {list(genotype.columns[:3])}")
    marker_to_row = genotype.drop_duplicates(subset=[marker_column]).set_index(marker_column)
    issues: list[str] = []
    for loci in loci_by_fungicide.values():
        for locus in loci.values():
            marker = locus["marker"]
            if marker not in marker_to_row.index:
                issues.append(f"{marker}: not present in genotype rs# column")
                continue
            raw_position = pd.to_numeric(marker_to_row.loc[marker, position_column], errors="coerce")
            if pd.isna(raw_position) or int(raw_position) != locus["position"]:
                issues.append(f"{marker}: marker-ID position differs from genotype pos={raw_position}")
    if issues:
        raise ValueError("GWAS markers did not match genotype metadata:\n" + "\n".join(issues[:20]))


def annotate_loci(
    loci_by_fungicide: dict[str, dict[str, dict[str, Any]]], annotation: pd.DataFrame
) -> pd.DataFrame:
    """Return one row per marker-gene overlap within the default +/-3 kb window."""
    genes_by_scaffold = {
        scaffold: group.sort_values("gene_start")
        for scaffold, group in annotation.groupby("scaffold", sort=False)
    }
    shared_markers = set(loci_by_fungicide["propiconazole"]) & set(loci_by_fungicide["tebuconazole"])
    output: list[dict[str, Any]] = []

    for fungicide, loci in loci_by_fungicide.items():
        for locus in loci.values():
            left = max(1, locus["position"] - DEFAULT_WINDOW_BP)
            right = locus["position"] + DEFAULT_WINDOW_BP
            candidates = genes_by_scaffold.get(locus["scaffold"], pd.DataFrame())
            overlaps = candidates.loc[
                (candidates["gene_start"] <= right) & (candidates["gene_stop"] >= left)
            ] if not candidates.empty else candidates
            shared = locus["marker"] in shared_markers
            common = {
                "fungicide": fungicide,
                "marker": locus["marker"],
                "scaffold": locus["scaffold"],
                "position": locus["position"],
                "paper_chromosome": locus["chromosome"],
                "detected_by": ";".join(locus["methods"]),
                "method_call_count": locus["method_call_count"],
                "top_lod": locus["top_lod"],
                "effect_at_top_lod": locus["effect_at_top_lod"],
                "r2_percent_at_top_lod": locus["r2_percent_at_top_lod"],
                "maf_at_top_lod": locus["maf_at_top_lod"],
                "shared_exact_marker_across_fungicides": shared,
                "window_start": left,
                "window_end": right,
                "window_bp_each_side": DEFAULT_WINDOW_BP,
                "window_note": "Default +/-3 kb only; manually adjusted LD boundaries are not reconstructed.",
            }
            if overlaps.empty:
                output.append({**common, "gene_id": "", "product_description": "", "gene_start": "", "gene_stop": "", "distance_marker_to_gene_bp": ""})
                continue
            for _, gene in overlaps.iterrows():
                distance = 0 if gene["gene_start"] <= locus["position"] <= gene["gene_stop"] else min(
                    abs(locus["position"] - int(gene["gene_start"])), abs(locus["position"] - int(gene["gene_stop"]))
                )
                output.append(
                    {
                        **common,
                        "gene_id": str(gene["gene_id"]),
                        "product_description": str(gene["product_description"]),
                        "gene_start": int(gene["gene_start"]),
                        "gene_stop": int(gene["gene_stop"]),
                        "distance_marker_to_gene_bp": int(distance),
                    }
                )
    return pd.DataFrame(output)


def minor_allele_indicator(values: pd.Series) -> tuple[np.ndarray | None, list[str]]:
    """Encode one of two nucleotide alleles as 0/1; A/C/G/T labels stay nominal."""
    clean = values.astype("string").str.strip().str.upper().replace({"": pd.NA, "NA": pd.NA, "NAN": pd.NA})
    counts = clean.dropna().value_counts()
    alleles = sorted(str(value) for value in counts.index.tolist())
    if len(alleles) != 2:
        return None, alleles
    minor = str(counts.idxmin())
    indicator = (clean == minor).astype("Float64").to_numpy(dtype=float, na_value=np.nan)
    return indicator, alleles


def pairwise_r2(left: np.ndarray | None, right: np.ndarray | None) -> tuple[float | None, int]:
    """Calculate pairwise-complete r-squared between two biallelic markers."""
    if left is None or right is None:
        return None, 0
    valid = np.isfinite(left) & np.isfinite(right)
    n = int(valid.sum())
    if n < 3 or np.unique(left[valid]).size != 2 or np.unique(right[valid]).size != 2:
        return None, n
    # Compute Pearson r-squared directly. Avoid np.corrcoef here because some
    # Windows conda builds can crash inside its covariance helper.
    x_values = left[valid].tolist()
    y_values = right[valid].tolist()
    mean_x = sum(x_values) / n
    mean_y = sum(y_values) / n
    ss_x = sum((value - mean_x) ** 2 for value in x_values)
    ss_y = sum((value - mean_y) ** 2 for value in y_values)
    if ss_x == 0 or ss_y == 0:
        return None, n
    cross_product = sum((x_values[index] - mean_x) * (y_values[index] - mean_y) for index in range(n))
    r_squared = (cross_product * cross_product) / (ss_x * ss_y)
    return float(min(1.0, max(0.0, r_squared))), n


def build_ld_review(
    genotype: pd.DataFrame,
    isolate_ids: list[str],
    loci_by_fungicide: dict[str, dict[str, dict[str, Any]]],
    annotation: pd.DataFrame,
    annotated_loci: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Review local LD around the nearest reported QTN for each paper candidate."""
    marker_column = genotype.columns[0]
    marker_to_index = {str(marker).strip(): index for index, marker in enumerate(genotype[marker_column])}
    candidates_by_id = annotation.drop_duplicates(subset=["gene_id"]).set_index("gene_id")
    all_loci = [locus for loci in loci_by_fungicide.values() for locus in loci.values()]
    summary_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    priority_rows: list[dict[str, Any]] = []
    exact_shared = set(loci_by_fungicide["propiconazole"]) & set(loci_by_fungicide["tebuconazole"])

    for gene_id, paper_note in PAPER_CANDIDATES.items():
        if gene_id not in candidates_by_id.index:
            raise ValueError(f"Paper candidate gene {gene_id} was not found in full_annotation.xlsx")
        gene = candidates_by_id.loc[gene_id]
        scaffold = str(gene["scaffold"])
        start, stop = int(gene["gene_start"]), int(gene["gene_stop"])
        same_scaffold = [locus for locus in all_loci if locus["scaffold"] == scaffold]
        nearest = min(
            same_scaffold,
            key=lambda locus: max(start - locus["position"], 0, locus["position"] - stop),
            default=None,
        )
        if nearest is None:
            distance_to_gene = None
        else:
            distance_to_gene = max(start - nearest["position"], 0, nearest["position"] - stop)

        gene_hits = annotated_loci.loc[annotated_loci["gene_id"] == gene_id]
        traits = sorted(set(gene_hits["fungicide"].astype(str)))
        marker_labels = sorted(set(f"{row.fungicide}::{row.marker}" for row in gene_hits.itertuples()))
        methods = sorted(
            {
                method
                for value in gene_hits["detected_by"].astype(str)
                for method in value.split(";")
                if method
            }
        )
        hit_marker_ids = set(gene_hits["marker"].astype(str))
        has_shared_locus = bool(hit_marker_ids & exact_shared)

        if len(traits) == 2 and has_shared_locus:
            priority = "1 - shared marker near paper-named gene"
        elif len(traits) == 2:
            priority = "1 - candidate region observed for both fungicides"
        elif gene_hits.shape[0]:
            priority = "2 - paper-named gene within default window"
        elif distance_to_gene is not None and distance_to_gene <= LD_SCAN_RADIUS_BP:
            priority = "3 - outside default window; inspect LD interval"
        else:
            priority = "3 - outside default window; coordinate/LD review needed"

        closest_qtn_for_gene = None
        if nearest is not None and distance_to_gene is not None and distance_to_gene <= LD_SCAN_RADIUS_BP:
            # Prefer the fungicide for which the gene is reported in the paper
            # when tied; otherwise keep the first stable row from the input map.
            candidate_loci = [locus for locus in same_scaffold if locus["position"] == nearest["position"]]
            closest_qtn_for_gene = candidate_loci[0] if candidate_loci else nearest

        gene_interval_variants = 0
        max_gene_r2: float | None = None
        gene_high_ld_count = 0
        nearest_high_ld_distance: int | None = None
        nearest_high_ld_marker = ""
        spanning_high_ld_distance: int | None = None
        spanning_high_ld_marker = ""
        spanning_high_ld_r2: float | None = None
        scanned_count = 0

        if closest_qtn_for_gene is not None:
            qtn = closest_qtn_for_gene
            qtn_index = marker_to_index.get(qtn["marker"])
            if qtn_index is None:
                raise ValueError(f"Lead QTN {qtn['marker']} is absent from genotype calls.")
            qtn_values, qtn_alleles = minor_allele_indicator(genotype.loc[qtn_index, isolate_ids])
            qtn_marker_rows = genotype.loc[genotype[marker_column].astype(str).str.strip().str.startswith("S" + scaffold + ".")]
            local_rows: list[tuple[int, str, int]] = []
            for row_index, marker in qtn_marker_rows[marker_column].items():
                marker = str(marker).strip()
                try:
                    marker_scaffold, position = parse_marker(marker)
                except ValueError:
                    continue
                if marker_scaffold == scaffold and abs(position - qtn["position"]) <= LD_SCAN_RADIUS_BP:
                    local_rows.append((int(row_index), marker, position))
            scanned_count = len(local_rows)

            for row_index, neighbor_marker, position in local_rows:
                neighbor_values, neighbor_alleles = minor_allele_indicator(genotype.loc[row_index, isolate_ids])
                r2, pair_n = pairwise_r2(qtn_values, neighbor_values)
                overlaps_gene = start <= position <= stop
                if overlaps_gene and r2 is not None:
                    gene_interval_variants += 1
                    max_gene_r2 = r2 if max_gene_r2 is None else max(max_gene_r2, r2)
                    if r2 >= HIGH_LD_R2:
                        gene_high_ld_count += 1
                if neighbor_marker != qtn["marker"] and r2 is not None and r2 >= HIGH_LD_R2:
                    distance = abs(position - qtn["position"])
                    if nearest_high_ld_distance is None or distance < nearest_high_ld_distance:
                        nearest_high_ld_distance = distance
                        nearest_high_ld_marker = neighbor_marker
                    # A high-LD marker can support an expanded interval only
                    # when it lies beyond the gene in the QTN-to-gene direction.
                    spans_gene = (
                        qtn["position"] < start and position >= stop
                    ) or (
                        qtn["position"] > stop and position <= start
                    )
                    if spans_gene and (spanning_high_ld_distance is None or distance < spanning_high_ld_distance):
                        spanning_high_ld_distance = distance
                        spanning_high_ld_marker = neighbor_marker
                        spanning_high_ld_r2 = r2
                if r2 is not None:
                    detail_rows.append(
                        {
                            "gene_id": gene_id,
                            "fungicide": qtn["fungicide"],
                            "qtn_marker": qtn["marker"],
                            "qtn_position": qtn["position"],
                            "gene_scaffold": scaffold,
                            "gene_start": start,
                            "gene_stop": stop,
                            "neighbor_marker": neighbor_marker,
                            "neighbor_position": position,
                            "distance_from_qtn_bp": abs(position - qtn["position"]),
                            "neighbor_alleles": ";".join(neighbor_alleles),
                            "qtn_alleles": ";".join(qtn_alleles),
                            "r2": r2,
                            "pairwise_isolates": pair_n,
                            "overlaps_gene_interval": overlaps_gene,
                            "high_ld_marker_spans_gene_interval": bool(
                                r2 >= HIGH_LD_R2
                                and neighbor_marker != qtn["marker"]
                                and ((qtn["position"] < start and position >= stop)
                                     or (qtn["position"] > stop and position <= start))
                            ),
                        }
                    )

        summary_rows.append(
            {
                "gene_id": gene_id,
                "product_description": str(gene["product_description"]),
                "paper_context": paper_note,
                "gene_scaffold": scaffold,
                "gene_start": start,
                "gene_stop": stop,
                "local_default_window_fungicides": ";".join(traits),
                "local_default_window_markers": ";".join(marker_labels),
                "methods": ";".join(methods),
                "shared_exact_marker_in_default_window": has_shared_locus,
                "nearest_qtn_fungicide": closest_qtn_for_gene["fungicide"] if closest_qtn_for_gene else "",
                "nearest_qtn_marker": closest_qtn_for_gene["marker"] if closest_qtn_for_gene else "",
                "nearest_qtn_position": closest_qtn_for_gene["position"] if closest_qtn_for_gene else "",
                "distance_qtn_to_gene_interval_bp": distance_to_gene if nearest is not None else "",
                "variants_scanned_within_100kb": scanned_count,
                "genotype_variants_in_gene_interval_with_pairwise_r2": gene_interval_variants,
                "max_r2_with_variant_inside_gene": max_gene_r2 if max_gene_r2 is not None else "",
                "gene_interval_variants_r2_ge_0_75": gene_high_ld_count,
                "nearest_nonlead_marker_r2_ge_0_75_bp": nearest_high_ld_distance if nearest_high_ld_distance is not None else "",
                "nearest_nonlead_marker_r2_ge_0_75": nearest_high_ld_marker,
                "high_ld_marker_spanning_gene_interval": spanning_high_ld_marker,
                "spanning_marker_distance_from_qtn_bp": spanning_high_ld_distance if spanning_high_ld_distance is not None else "",
                "spanning_marker_r2": spanning_high_ld_r2 if spanning_high_ld_r2 is not None else "",
                "priority": priority,
                "interpretation_note": "Exploratory; does not reproduce manually adjusted LD windows or establish causality.",
            }
        )

    summary = pd.DataFrame(summary_rows)
    details = pd.DataFrame(detail_rows)
    prioritized = summary[
        [
            "priority", "gene_id", "product_description", "paper_context",
            "local_default_window_fungicides", "local_default_window_markers", "methods",
            "shared_exact_marker_in_default_window", "nearest_qtn_fungicide", "nearest_qtn_marker",
            "distance_qtn_to_gene_interval_bp", "nearest_nonlead_marker_r2_ge_0_75_bp",
            "nearest_nonlead_marker_r2_ge_0_75", "high_ld_marker_spanning_gene_interval",
            "spanning_marker_distance_from_qtn_bp", "spanning_marker_r2", "interpretation_note",
        ]
    ].copy()
    return summary, details, prioritized


def make_crosscheck_report(
    by_fungicide: dict[str, pd.DataFrame],
    loci_by_fungicide: dict[str, dict[str, dict[str, Any]]],
) -> str:
    """Summarize reproducible counts from final-result tables."""
    lines = [
        "# GWAS Results Cross-Check",
        "",
        "Generated from the repository's final mrMLM result CSVs by `06_post_gwas_candidate_review.py`.",
        "",
        "## Method-call counts and unique markers",
        "",
    ]
    for fungicide, frame in by_fungicide.items():
        lines.append(f"### {fungicide.capitalize()}")
        lines.append("")
        for method, group in frame.groupby("method", sort=True):
            lines.append(f"- {method}: {len(group)} calls")
        lines.append(f"- Unique marker IDs after cross-method deduplication: {len(loci_by_fungicide[fungicide])}")
        multi_method = sum(1 for locus in loci_by_fungicide[fungicide].values() if len(locus["methods"]) > 1)
        lines.append(f"- Markers detected by more than one method: {multi_method}")
        lines.append("")

    prop = set(loci_by_fungicide["propiconazole"])
    tebu = set(loci_by_fungicide["tebuconazole"])
    shared = prop & tebu
    distinct_total = len(prop | tebu)
    lines.extend(
        [
            "## Cross-fungicide comparison",
            "",
            f"- Shared exact marker IDs: {len(shared)}",
            f"- Distinct marker IDs across both fungicides: {distinct_total}",
            "- The supplied preprint's Results section reports 13 shared QTNs while its Abstract says 12; these outputs match the Results count.",
            "- These are mrMLM output counts, not a new significance test or external validation.",
            "",
        ]
    )
    lines.extend(["## Reported lead loci", ""])
    lead_checks = [
        ("propiconazole", "SHG970334.1_6640556", "mrMLM", "~10% variance in the Results text"),
        ("tebuconazole", "SHG970333.1_8322036", "pKWmEB", "3.06% variance in the Results text"),
    ]
    for fungicide, marker, method, paper_value in lead_checks:
        rows = by_fungicide[fungicide]
        match = rows.loc[(rows["marker"] == marker) & (rows["method"] == method)]
        if match.empty:
            lines.append(f"- {fungicide}: manuscript lead marker {marker} / {method} was not found in the final result table.")
        else:
            row = match.iloc[0]
            lines.append(
                f"- {fungicide}: {marker} on chromosome {row['chromosome']} has r²={row['r2_percent']}% "
                f"under {method}; the supplied Results section reports {paper_value}."
            )
    lines.extend(["", "## Phenotype scale and range", ""])
    lines.append(
        "The study author confirmed `trait1` is an isolate-level least-squares mean EC50 in mM. "
        "The GWAS uses the stored values without log transformation. Converted summaries below are for comparison only."
    )
    for fungicide, path in PHENOTYPE_PATHS.items():
        phenotype = pd.read_csv(path, dtype=str, keep_default_na=True)
        trait_name = next(column for column in phenotype.columns if str(column).strip().casefold() == "trait1")
        trait = pd.to_numeric(phenotype[trait_name], errors="coerce")
        mean_um = float(trait.mean()) * 1_000
        minimum_um = float(trait.min()) * 1_000
        maximum_um = float(trait.max()) * 1_000
        lines.append(
            f"- {fungicide}: n={len(trait)}, mean={mean_um:.2f} µM, range={minimum_um:.2f}–{maximum_um:.2f} µM."
        )
    lines.append("")
    return "\n".join(lines)


def write_candidate_report(
    summary: pd.DataFrame,
    loci_by_fungicide: dict[str, dict[str, dict[str, Any]]],
) -> str:
    """Write an interpretable gene-priority report from the generated tables."""
    prop_markers = set(loci_by_fungicide["propiconazole"])
    tebu_markers = set(loci_by_fungicide["tebuconazole"])
    shared_markers = prop_markers & tebu_markers
    lines = [
        "# Candidate Gene Prioritization",
        "",
        "This report is generated from the final GWAS result CSVs and `data/reference/full_annotation.xlsx`.",
        "It is a hypothesis shortlist, not evidence that a gene causes fungicide resistance.",
        "",
        "## Highest-priority cross-fungicide candidates",
        "",
    ]
    top = summary.loc[summary["priority"].str.startswith("1 -")]
    if top.empty:
        lines.append("No paper-mentioned candidate gene was found in default-window regions shared across traits.")
    for _, row in top.iterrows():
        lines.append(
            f"- **{row['gene_id']} ({row['product_description']})**: "
            f"{row['local_default_window_markers'] or 'nearby locus'}; "
            f"observed for {row['local_default_window_fungicides'] or 'both fungicides'}; "
            f"methods: {row['methods'] or 'not recorded'}."
        )
    lines.extend(["", "## Other paper-mentioned candidates", ""])
    for _, row in summary.loc[~summary["priority"].str.startswith("1 -")].iterrows():
        qtn_note = ""
        if str(row["nearest_qtn_marker"]):
            qtn_note = (
                f" Nearest QTN: {row['nearest_qtn_marker']} "
                f"({row['distance_qtn_to_gene_interval_bp']} bp from the gene interval)."
            )
        spanning_marker = str(row["high_ld_marker_spanning_gene_interval"])
        if spanning_marker:
            qtn_note += (
                f" A marker beyond the gene in the QTN-to-gene direction, {spanning_marker}, "
                f"has r²={float(row['spanning_marker_r2']):.3f}."
            )
        elif str(row["nearest_nonlead_marker_r2_ge_0_75"]):
            qtn_note += (
                f" The nearest non-lead marker with r²>=0.75 is "
                f"{row['nearest_nonlead_marker_r2_ge_0_75']} at "
                f"{row['nearest_nonlead_marker_r2_ge_0_75_bp']} bp, but it does not span the gene interval."
            )
        lines.append(
            f"- **{row['gene_id']} ({row['product_description']})** — {row['priority']}."
            f" {row['paper_context']}{qtn_note}"
        )
    lines.extend(
        [
            "",
            "## Boundaries",
            "",
            "Gene overlaps use the paper's default +/-3 kb window only. The authors manually adjusted some boundaries after reviewing local LD. The accompanying LD scan is exploratory, covers at most 100 kb around each nearest reported QTN, and does not reconstruct those manual decisions.",
            "",
            f"The final result files contain {len(prop_markers)} unique propiconazole markers and {len(tebu_markers)} unique tebuconazole markers, with {len(shared_markers)} exact markers shared across fungicides.",
            "",
        ]
    )
    return "\n".join(lines)


def build_candidate_evidence_matrix(
    summary: pd.DataFrame,
    loci_by_fungicide: dict[str, dict[str, dict[str, Any]]],
) -> pd.DataFrame:
    """Combine reported QTN statistics with gene proximity and local LD evidence."""
    rows: list[dict[str, Any]] = []
    for row in summary.to_dict(orient="records"):
        fungicide = str(row.get("nearest_qtn_fungicide", ""))
        marker = str(row.get("nearest_qtn_marker", ""))
        qtn = loci_by_fungicide.get(fungicide, {}).get(marker, {})
        scaffold, position = parse_marker(marker) if marker else ("", "")

        if row.get("shared_exact_marker_in_default_window"):
            evidence_class = "Exact marker shared across fungicides within default gene window"
        elif row.get("high_ld_marker_spanning_gene_interval"):
            evidence_class = "Nearby QTN plus a high-LD marker spanning the gene interval"
        elif row.get("local_default_window_markers"):
            evidence_class = "Nearby reported QTN within default gene window"
        else:
            evidence_class = "Candidate mentioned by study; gene-level link needs further evidence"

        rows.append(
            {
                "gene_id": row.get("gene_id", ""),
                "product_description": row.get("product_description", ""),
                "study_candidate_rationale": row.get("paper_context", ""),
                "evidence_class": evidence_class,
                "nearest_qtn_fungicide": fungicide,
                "nearest_qtn_marker": marker,
                "nearest_qtn_chromosome": qtn.get("chromosome", ""),
                "nearest_qtn_position_bp": position,
                "nearest_qtn_methods": ";".join(qtn.get("methods", [])),
                "nearest_qtn_effect_at_top_lod": qtn.get("effect_at_top_lod", ""),
                "nearest_qtn_top_lod": qtn.get("top_lod", ""),
                "nearest_qtn_r2_percent_at_top_lod": qtn.get("r2_percent_at_top_lod", ""),
                "nearest_qtn_maf_at_top_lod": qtn.get("maf_at_top_lod", ""),
                "distance_qtn_to_gene_interval_bp": row.get("distance_qtn_to_gene_interval_bp", ""),
                "reported_markers_within_default_window": row.get("local_default_window_markers", ""),
                "methods_within_default_window": row.get("methods", ""),
                "exact_marker_shared_across_fungicides": row.get("shared_exact_marker_in_default_window", False),
                "genotyped_variants_inside_gene_with_pairwise_r2": row.get("genotype_variants_in_gene_interval_with_pairwise_r2", ""),
                "maximum_qtn_r2_with_variant_inside_gene": row.get("max_r2_with_variant_inside_gene", ""),
                "nearest_nonlead_marker_r2_ge_0_75": row.get("nearest_nonlead_marker_r2_ge_0_75", ""),
                "nearest_nonlead_marker_distance_from_qtn_bp": row.get("nearest_nonlead_marker_r2_ge_0_75_bp", ""),
                "high_ld_marker_spanning_gene_interval": row.get("high_ld_marker_spanning_gene_interval", ""),
                "spanning_marker_r2": row.get("spanning_marker_r2", ""),
                "interpretation_limit": "Association/proximity/LD are prioritization evidence only; none establish that this gene causes fungicide sensitivity.",
            }
        )
    return pd.DataFrame(rows)


def write_candidate_evidence_report(evidence: pd.DataFrame) -> str:
    """Write a compact guide to the evidence matrix and its interpretation."""
    lines = [
        "# Candidate Evidence Matrix",
        "",
        "Generated from the final mrMLM result tables, gene annotation, and local LD scan by `06_post_gwas_candidate_review.py`.",
        "",
        "The CSV keeps reported association statistics, gene proximity, and genotype-based LD evidence in separate columns. `study_candidate_rationale` summarizes how the source study discussed the gene; it is not independent functional validation.",
        "",
        "The default gene window is ±3 kb. The exploratory LD scan extends to 100 kb around the nearest reported QTN and evaluates biallelic nucleotide markers. An exact shared marker means the same marker ID was reported for both fungicides within the default gene window.",
        "",
        "This is a prioritization aid, not a causal-gene claim. Candidate signals require confirmation from the final paper/supplementary material and independent biological or population-level validation.",
        "",
        f"Candidates summarized: {len(evidence)}",
        "",
        "Full table: `candidate_evidence_matrix.csv`.",
        "",
    ]
    return "\n".join(lines)


def write_ld_report(summary: pd.DataFrame) -> str:
    """Write a gene-by-gene explanation of local LD evidence."""
    lines = [
        "# Local LD Review",
        "",
        "Pairwise r-squared was calculated from biallelic nucleotide indicators in the 152 genotype columns. "
        "Nucleotide labels were treated as nominal alleles. The scan includes genotype markers on the same "
        "normalized scaffold within 100 kb of each paper-mentioned candidate gene's nearest reported QTN.",
        "",
    ]
    for _, row in summary.iterrows():
        lines.extend([f"## {row['gene_id']} — {row['product_description']}", ""])
        if not row["nearest_qtn_marker"]:
            lines.extend(["No reported QTN was found on the gene's scaffold within the scan radius.", ""])
            continue
        lines.append(
            f"Nearest {row['nearest_qtn_fungicide']} QTN: `{row['nearest_qtn_marker']}`; "
            f"distance to gene interval: {row['distance_qtn_to_gene_interval_bp']} bp."
        )
        lines.append(
            f"Biallelic genotype markers inside gene interval with calculable r²: "
            f"{row['genotype_variants_in_gene_interval_with_pairwise_r2']}."
        )
        spanning_marker = str(row["high_ld_marker_spanning_gene_interval"])
        if spanning_marker:
            lines.append(
                f"A non-lead marker `{spanning_marker}` lies beyond the gene interval in the QTN-to-gene "
                f"direction and has r²={float(row['spanning_marker_r2']):.3f} with the QTN "
                f"({row['spanning_marker_distance_from_qtn_bp']} bp from the QTN). The interval from the QTN "
                "to that marker spans the gene; this may explain a manually expanded candidate region."
            )
        else:
            nearest_ld = str(row["nearest_nonlead_marker_r2_ge_0_75"])
            if nearest_ld:
                lines.append(
                    f"The nearest non-lead marker with r²>=0.75 is `{nearest_ld}` "
                    f"({row['nearest_nonlead_marker_r2_ge_0_75_bp']} bp from the QTN), but it does not "
                    "extend across the gene in the QTN-to-gene direction. It does not directly support "
                    "including this gene in an LD-expanded window."
                )
            else:
                lines.append("No non-lead marker with r²>=0.75 was found within the scan radius.")
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "The scan is exploratory and is not a reconstruction of the authors' manually adjusted LD windows. "
            "It does not prove that any nearby gene is causal. Full pairwise values are in "
            "`candidate_gene_ld_variants.csv`; gene-level summaries are in `candidate_gene_ld_check.csv`.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    require_inputs()
    genotype = pd.read_csv(GENOTYPE_PATH, dtype=str, keep_default_na=False, low_memory=False)
    isolate_ids = validate_isolate_ids(genotype)
    annotation = read_annotation_workbook(ANNOTATION_PATH)
    by_fungicide, loci_by_fungicide = read_gwas_results()
    validate_markers_against_genotype(genotype, loci_by_fungicide)

    # The genotype matrix is keyed by rs#; the first three columns are metadata.
    # Check allele values are nucleotide calls before using them in biallelic LD.
    genotype_calls = genotype.iloc[:, 3:].astype("string").apply(lambda col: col.str.strip().str.upper())
    observed_calls = set(genotype_calls.stack().dropna().astype(str))
    unexpected = observed_calls - {"A", "C", "G", "T"}
    if unexpected:
        raise ValueError(f"Unexpected genotype code(s) prevent nucleotide-based LD: {sorted(unexpected)[:20]}")

    annotated_locus_rows = annotate_loci(loci_by_fungicide, annotation)
    ld_summary, ld_details, prioritization = build_ld_review(
        genotype, isolate_ids, loci_by_fungicide, annotation, annotated_locus_rows
    )
    evidence_matrix = build_candidate_evidence_matrix(ld_summary, loci_by_fungicide)

    GWAS_DIR.mkdir(parents=True, exist_ok=True)
    REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    annotated_locus_rows.to_csv(GWAS_DIR / "candidate_loci_annotation.csv", index=False)
    prioritization.to_csv(REVIEW_DIR / "candidate_gene_prioritization.csv", index=False)
    ld_summary.to_csv(REVIEW_DIR / "candidate_gene_ld_check.csv", index=False)
    ld_details.to_csv(REVIEW_DIR / "candidate_gene_ld_variants.csv", index=False)
    evidence_matrix.to_csv(REVIEW_DIR / "candidate_evidence_matrix.csv", index=False)
    (REVIEW_DIR / "candidate_evidence_matrix.md").write_text(
        write_candidate_evidence_report(evidence_matrix), encoding="utf-8"
    )
    (REVIEW_DIR / "candidate_gene_prioritization.md").write_text(
        write_candidate_report(ld_summary, loci_by_fungicide), encoding="utf-8"
    )
    (REVIEW_DIR / "candidate_gene_ld_review.md").write_text(
        write_ld_report(ld_summary), encoding="utf-8"
    )
    (REPORT_DIR / "gwas_results_crosscheck.md").write_text(
        make_crosscheck_report(by_fungicide, loci_by_fungicide), encoding="utf-8"
    )

    print("Post-GWAS candidate review completed.")
    print(f"Genotype matrix: {len(genotype):,} variants x {len(isolate_ids):,} isolates")
    print(f"Annotation genes read: {len(annotation):,}")
    print(f"Candidate annotation rows: {len(annotated_locus_rows):,}")
    print(f"Candidate-gene LD rows: {len(ld_details):,}")
    print(f"Candidate evidence rows: {len(evidence_matrix):,}")
    print(f"Candidate annotation: {GWAS_DIR / 'candidate_loci_annotation.csv'}")
    print(f"Review reports: {REVIEW_DIR}")
    print(f"Cross-check: {REPORT_DIR / 'gwas_results_crosscheck.md'}")
    print("All outputs are derived; raw data and original GWAS files were not modified.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Post-GWAS candidate review failed: {error}", file=sys.stderr)
        raise
