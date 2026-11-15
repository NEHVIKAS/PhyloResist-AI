"""Extract and validate isolate provenance from the study supplement.

Run from the repository root in the phyloresist environment:
    python notebooks/08_prepare_isolate_metadata.py

Reads the original workbook under data/reference and the isolate headers in the
raw genotype file. It writes a matched, source-traceable CSV and a match report.
The source workbook and the three raw CSVs are never modified.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime, timedelta
from pathlib import Path, PurePosixPath
from typing import Any
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_WORKBOOK = ROOT / "data" / "reference" / "fungicide_study_supplement_media-1.xlsx"
GENOTYPE_PATH = ROOT / "data" / "raw" / "genotype.csv"
OUTPUT_PATH = ROOT / "data" / "reference" / "isolate_metadata.csv"
REPORT_PATH = ROOT / "results" / "reports" / "isolate_metadata_match_report.md"

EXCEL_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
BUILTIN_DATE_FORMAT_IDS = {14, 15, 16, 17, 18, 19, 20, 21, 22, 45, 46, 47}


def normalize_id(value: Any) -> str:
    """Normalize numeric-like identifiers without changing leading zeroes."""
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0+", text):
        text = text.split(".", 1)[0]
    return text


def excel_column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference.upper())
    if not letters:
        raise ValueError(f"Invalid Excel cell reference: {reference!r}")
    index = 0
    for letter in letters.group(0):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def read_workbook_table(path: Path, wanted_sheet: str) -> tuple[list[dict[str, str]], int]:
    """Read an OpenXML worksheet with the standard library (no Excel add-on)."""
    with ZipFile(path) as workbook:
        workbook_root = ET.fromstring(workbook.read("xl/workbook.xml"))
        relationships_root = ET.fromstring(workbook.read("xl/_rels/workbook.xml.rels"))
        relationships = {
            relation.attrib["Id"]: relation.attrib["Target"]
            for relation in relationships_root.findall(f"{{{PKG_REL_NS}}}Relationship")
        }

        sheet_path = ""
        for sheet in workbook_root.findall(f".//{{{EXCEL_NS}}}sheet"):
            if sheet.attrib.get("name", "").strip().casefold() == wanted_sheet.casefold():
                target = relationships.get(sheet.attrib.get(f"{{{REL_NS}}}id", ""), "")
                if target:
                    sheet_path = (
                        target.lstrip("/") if target.startswith("/")
                        else str(PurePosixPath("xl") / PurePosixPath(target))
                    )
                break
        if not sheet_path or sheet_path not in workbook.namelist():
            raise ValueError(f"Worksheet {wanted_sheet!r} was not found in {path.name}.")

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in workbook.namelist():
            shared_root = ET.fromstring(workbook.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(node.text or "" for node in item.iter(f"{{{EXCEL_NS}}}t"))
                for item in shared_root.findall(f"{{{EXCEL_NS}}}si")
            ]

        # Resolve the cell style for the collection-date field. The source uses
        # Excel's built-in month/year number format for serial-date entries.
        style_num_formats: list[int] = []
        if "xl/styles.xml" in workbook.namelist():
            styles_root = ET.fromstring(workbook.read("xl/styles.xml"))
            xfs = styles_root.find(f"{{{EXCEL_NS}}}cellXfs")
            if xfs is not None:
                style_num_formats = [int(xf.attrib.get("numFmtId", "0")) for xf in xfs]

        sheet_root = ET.fromstring(workbook.read(sheet_path))
        sheet_data = sheet_root.find(f"{{{EXCEL_NS}}}sheetData")
        if sheet_data is None:
            raise ValueError(f"Worksheet {wanted_sheet!r} has no sheet data.")

        rows: list[tuple[int, dict[int, tuple[str, bool]]]] = []
        for row_element in sheet_data.findall(f"{{{EXCEL_NS}}}row"):
            cells: dict[int, tuple[str, bool]] = {}
            for cell in row_element.findall(f"{{{EXCEL_NS}}}c"):
                column = excel_column_index(cell.attrib.get("r", ""))
                cell_type = cell.attrib.get("t", "")
                value_node = cell.find(f"{{{EXCEL_NS}}}v")
                value = value_node.text if value_node is not None and value_node.text else ""
                if cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{EXCEL_NS}}}t"))
                elif cell_type == "s" and value:
                    value = shared_strings[int(value)]
                style_index = int(cell.attrib.get("s", "0"))
                num_format_id = style_num_formats[style_index] if style_index < len(style_num_formats) else 0
                cells[column] = (value.strip(), num_format_id in BUILTIN_DATE_FORMAT_IDS)
            if cells:
                rows.append((int(row_element.attrib.get("r", "0")), cells))

    if len(rows) < 3:
        raise ValueError(f"Worksheet {wanted_sheet!r} does not contain a title, header and data rows.")
    header_row = rows[1][1]
    headers = [header_row.get(column, ("", False))[0] for column in range(max(header_row) + 1)]
    if not headers or headers[0].strip().casefold() not in {"isolate", "isolate id", "isolate_id"}:
        raise ValueError(f"Unexpected header row in {wanted_sheet!r}: {headers}")

    records: list[dict[str, str]] = []
    for excel_row, cells in rows[2:]:
        record = {headers[column]: cells.get(column, ("", False))[0] for column in range(len(headers))}
        record["_excel_row"] = str(excel_row)
        # Record whether the source's date cell is explicitly formatted as a date.
        date_index = next((i for i, name in enumerate(headers) if name.strip().casefold() == "date of collection"), None)
        record["_collection_date_is_excel_date"] = (
            "true" if date_index is not None and cells.get(date_index, ("", False))[1] else "false"
        )
        records.append(record)
    return records, len(records)


def normalize_collection_date(raw: str, is_excel_date: str) -> tuple[str, str]:
    """Return a faithfully reported date label and its four-digit year."""
    text = str(raw).strip()
    if not text:
        return "", ""
    if is_excel_date.casefold() == "true":
        try:
            serial = float(text)
            date_value = datetime(1899, 12, 30) + timedelta(days=serial)
            # The workbook style is mmm-yy: retain month/year, not an invented day.
            date_label = date_value.strftime("%Y-%m")
            return date_label, date_value.strftime("%Y")
        except (ValueError, OverflowError):
            pass
    if re.fullmatch(r"\d{4}(?:\.0+)?", text):
        year = text.split(".", 1)[0]
        return year, year
    return text, text[:4] if re.match(r"\d{4}", text) else ""


def main() -> int:
    if not SOURCE_WORKBOOK.is_file():
        raise FileNotFoundError(f"Study supplement not found: {SOURCE_WORKBOOK}")
    if not GENOTYPE_PATH.is_file():
        raise FileNotFoundError(f"Genotype input not found: {GENOTYPE_PATH}")

    genotype_headers = pd.read_csv(GENOTYPE_PATH, dtype=str, nrows=0).columns.tolist()
    if len(genotype_headers) <= 3:
        raise ValueError("Expected three variant metadata columns followed by isolate IDs.")
    genotype_ids = [normalize_id(value) for value in genotype_headers[3:]]
    if len(genotype_ids) != len(set(genotype_ids)) or any(not value for value in genotype_ids):
        raise ValueError("Genotype isolate IDs are blank or duplicated.")

    records, source_row_count = read_workbook_table(SOURCE_WORKBOOK, "Supplemental Table 1")
    supplement = pd.DataFrame(records)
    isolate_column = next(column for column in supplement.columns if column.strip().casefold() in {"isolate", "isolate id", "isolate_id"})
    supplement["isolate_id"] = supplement[isolate_column].map(normalize_id)
    source_rows = supplement.loc[supplement["isolate_id"].ne("")].copy()
    source_ids = set(source_rows["isolate_id"])
    genotype_set = set(genotype_ids)
    missing_from_supplement = sorted(genotype_set - source_ids)
    supplement_only = sorted(source_ids - genotype_set)
    if missing_from_supplement:
        raise ValueError(
            "Some genotype isolate IDs are absent from the supplement: "
            + ", ".join(missing_from_supplement[:25])
        )

    matched_rows = source_rows.loc[source_rows["isolate_id"].isin(genotype_set)].copy()
    duplicates = matched_rows.loc[matched_rows["isolate_id"].duplicated(False), "isolate_id"].unique().tolist()
    if duplicates:
        raise ValueError(f"Matched isolate IDs are duplicated in the supplement: {duplicates[:20]}")
    by_id = matched_rows.set_index("isolate_id", drop=False)
    matched = by_id.loc[genotype_ids].reset_index(drop=True)
    column_lookup = {str(column).strip().casefold(): column for column in matched.columns}

    def field(label: str) -> pd.Series:
        actual = column_lookup.get(label.casefold())
        return matched[actual].astype(str).str.strip().replace({"nan": "", "None": ""}) if actual else pd.Series([""] * len(matched))

    date_reported, collection_year = zip(
        *(normalize_collection_date(value, is_date)
          for value, is_date in zip(field("Date of collection"), field("_collection_date_is_excel_date")))
    )
    output = pd.DataFrame(
        {
            "isolate_id": genotype_ids,
            "publication": field("publication"),
            "alternate_id": field("alternate ID"),
            "collection_state": field("State"),
            "collection_city_or_county": field("City/County"),
            "collection_date_as_reported": list(date_reported),
            "collection_year": list(collection_year),
            "host": field("host"),
            "chemotype": field("Chemotype"),
            "population_label": field("Population"),
            "state_validation_group": field("State"),
            "source_table": "Supplemental Table 1",
            "source_excel_row": field("_excel_row"),
        }
    )
    if output["isolate_id"].duplicated().any() or len(output) != len(genotype_ids):
        raise ValueError("Matched metadata output does not have exactly one row per genotype isolate.")

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT_PATH, index=False)

    state_counts = output["collection_state"].value_counts(dropna=False)
    population_counts = output["population_label"].value_counts(dropna=False)
    missing_counts = output.replace("", pd.NA).isna().sum()
    lines = [
        "# Isolate Metadata Match Report",
        "",
        "Source: `data/reference/fungicide_study_supplement_media-1.xlsx`, worksheet `Supplemental Table 1`.",
        "The uploaded source workbook is preserved; `isolate_metadata.csv` contains only rows matched to the current genotype IDs, in genotype-file order.",
        "",
        "## Match summary",
        "",
        f"- Genotype isolate IDs: {len(genotype_ids)}",
        f"- Worksheet rows containing stored cell data parsed: {source_row_count}",
        f"- Nonblank unique IDs in the supplement: {source_rows['isolate_id'].nunique()}",
        f"- Genotype isolates matched: {output['isolate_id'].nunique()} of {len(genotype_ids)}",
        f"- Supplement-only IDs not present in this genotype panel: {len(supplement_only)} ({', '.join(supplement_only) if supplement_only else 'none'})",
        f"- Genotype IDs missing from the supplement: {len(missing_from_supplement)}",
        "",
        "## Metadata coverage",
        "",
    ]
    lines.extend(f"- {column}: {len(output) - int(missing_counts[column])} of {len(output)} populated" for column in output.columns if column not in {"isolate_id", "source_table"})
    lines.extend(["", "## Population and geography", ""])
    lines.extend(f"- State {state}: n={count}" for state, count in state_counts.items())
    lines.extend(["", "Population labels:"])
    lines.extend(f"- {population}: n={count}" for population, count in population_counts.items())
    lines.extend(
        [
            "",
            "`Population` has few non-NA1 labels and is not suitable for a population-held-out model validation split in this panel. `state_validation_group` preserves the source State field for exploratory leave-one-state-out evaluation; state is a geographic grouping, not a synonym for an independent population. Some states have very small sample counts, so state-level scores need per-state sample counts and cautious interpretation.",
            "",
            "Collection dates formatted in the source as month/year are represented as YYYY-MM; entries supplied only as a year remain YYYY. No exact day is inferred from month/year-formatted Excel dates.",
            "",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")

    print("Isolate metadata extraction completed.")
    print(f"Matched isolates: {output['isolate_id'].nunique()} / {len(genotype_ids)}")
    print(f"Supplement-only isolate IDs: {', '.join(supplement_only) if supplement_only else 'none'}")
    print(f"Metadata CSV: {OUTPUT_PATH}")
    print(f"Match report: {REPORT_PATH}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Isolate metadata extraction failed: {error}", file=sys.stderr)
        raise
