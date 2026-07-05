from pathlib import Path
import pandas as pd
import re

folder = Path(".")  # folder containing the CSV files

dose_map = {
    "control": 0,
    "05": 0.5,
    "1": 1,
    "5": 5,
}

all_data = []

for file in folder.glob("*.csv"):
    stem = file.stem
    parts = stem.split("_")

    # Expected: 2A_caco2_05_nabu.csv
    if len(parts) < 3:
        print(f"Skipping {file.name}: unexpected filename")
        continue

    well = parts[0]
    cell_line = parts[1]
    dose_text = parts[2]

    match = re.match(r"(\d+)([A-Za-z]+)", well)
    if not match:
        print(f"Skipping {file.name}: could not parse well")
        continue

    well_column = int(match.group(1))
    well_row = match.group(2)

    if dose_text not in dose_map:
        print(f"Skipping {file.name}: unknown dose '{dose_text}'")
        continue

    df = pd.read_csv(file)

    if "Area" not in df.columns:
        print(f"Skipping {file.name}: no Area column")
        continue

    out = df[["Area"]].copy()
    out["well_column"] = well_column
    out["well_row"] = well_row
    out["cell_line"] = cell_line
    out["dose"] = dose_map[dose_text]
    out["source_file"] = file.name

    all_data.append(out)

combined = pd.concat(all_data, ignore_index=True)

combined.to_csv("combined_cell_areas.csv", index=False)

print("Saved combined_cell_areas.csv")
