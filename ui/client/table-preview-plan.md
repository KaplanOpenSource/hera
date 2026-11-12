# Plan: table preview (Parquet, CSV, Excel)

## Today

- `.parquet`, `.csv`, `.xlsx` and `.xls` files get no preview.
- Parquet: binary, typed columns. The header holds the row count and the columns.
- Parquet saved by Dask is a folder of part files, not one file.
- CSV: plain text. No types, no row count; both must be read from the file.
- Excel: one table per sheet. `.xlsx` is the new format, `.xls` the old one.
- Test Parquet files: `hera/doc/jupyter/User/data/lowfrequency/MASSADA.parquet`, `hera/doc/jupyter/User/data/experiment/data/Sonic.parquet`.
- No CSV or Excel file in the repo; tests make their own.
- In `requirements.txt`: `pandas`, `pyarrow`, `openpyxl` (reads `.xlsx`). Not there: `xlrd` (reads `.xls`).

## Phase 1

1. Server: one route, `GET /table?path=...&sheet=...`. The file ending picks the reader.
2. Server: read only the first `PREVIEW_ROWS` rows (start at 100). Never the whole file.
3. Server, Parquet: `pyarrow.dataset`, for a file or a folder alike; `count_rows()` and `head()`.
4. Server, CSV: `pandas.read_csv(nrows=PREVIEW_ROWS)`; types as pandas guesses them; row count by counting lines.
5. Server, Excel: `pandas.read_excel(sheet_name=..., nrows=PREVIEW_ROWS)`; first sheet if none given; row count from the sheet size.
6. Server: send JSON: `{ rowCount, columns: [{ name, type }], rows: [[...]], sheets }`. `sheets` only for Excel. Dates as ISO text, NaN as null.
7. Server: add `xlrd` to `requirements.txt` for `.xls` (ask first: the Docker image needs a rebuild).
8. Client: `ServerApi.table(path, sheet?)` with `getJson`.
9. Client: `isTableResource`: the resource ends in `.parquet`, `.csv`, `.xlsx` or `.xls`.
10. Client: `TablePreview` shows "N rows, M columns", then the rows in a MUI `DataGrid`.
11. Client: the column header shows the name, the type on hover.
12. Client: for Excel, MUI `Tabs` above the grid, one per sheet.
13. Client: add tables to `hasPreview` and `PreviewPanel`.
14. Add `@mui/x-data-grid` (free MIT version) to `package.json`.
15. Measure on a big CSV: how long the line count takes.

## Later

- More rows on scroll: `offset` and `limit` on `/table`.
- CSV with `;` or tabs: let pandas guess the separator.
- Skip the CSV line count if it is too slow on big files; show "100+ rows".

## Tests

- Server, each format: right row count, columns, types and first rows on a small file made in the test.
- Server: a Parquet folder of several parts reads like one file.
- Server: never more than `PREVIEW_ROWS` rows.
- Server: dates and NaN come out as ISO text and null.
- Server, Excel: sheet names listed; `sheet` picks the sheet.
- Server: a missing path gives a 404; an unknown ending gives a 400.
- Client: each ending gets a preview; the counts and rows reach the grid.
- Client: Excel sheet tabs switch the sheet.
- Client: a failed call shows no grid.
- Then the full `TEST_UI.md` checklist.
