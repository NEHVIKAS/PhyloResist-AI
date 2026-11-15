# Isolate Metadata Match Report

Source: `data/reference/fungicide_study_supplement_media-1.xlsx`, worksheet `Supplemental Table 1`.
The uploaded source workbook is preserved; `isolate_metadata.csv` contains only rows matched to the current genotype IDs, in genotype-file order.

## Match summary

- Genotype isolate IDs: 152
- Worksheet rows containing stored cell data parsed: 998
- Nonblank unique IDs in the supplement: 152
- Genotype isolates matched: 152 of 152
- Supplement-only IDs not present in this genotype panel: 0 (none)
- Genotype IDs missing from the supplement: 0

## Metadata coverage

- publication: 150 of 152 populated
- alternate_id: 152 of 152 populated
- collection_state: 152 of 152 populated
- collection_city_or_county: 150 of 152 populated
- collection_date_as_reported: 134 of 152 populated
- collection_year: 134 of 152 populated
- host: 150 of 152 populated
- chemotype: 152 of 152 populated
- population_label: 152 of 152 populated
- state_validation_group: 152 of 152 populated
- source_excel_row: 152 of 152 populated

## Population and geography

- State NY: n=73
- State MN: n=35
- State ND: n=33
- State MT: n=8
- State MI: n=1
- State VA: n=1
- State LA: n=1

Population labels:
- NA1: n=147
- unassigned: n=4
- NA2: n=1

`Population` has few non-NA1 labels and is not suitable for a population-held-out model validation split in this panel. `state_validation_group` preserves the source State field for exploratory leave-one-state-out evaluation; state is a geographic grouping, not a synonym for an independent population. Some states have very small sample counts, so state-level scores need per-state sample counts and cautious interpretation.

Collection dates formatted in the source as month/year are represented as YYYY-MM; entries supplied only as a year remain YYYY. No exact day is inferred from month/year-formatted Excel dates.
