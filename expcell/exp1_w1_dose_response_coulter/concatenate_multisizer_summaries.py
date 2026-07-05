from pathlib import Path
import pandas as pd

folder = Path(".")

rows = []

for file in folder.glob("*.#m4."):

    with open(file, "r", errors="ignore") as f:
        lines = f.readlines()

    stats = {}
    in_stats = False

    for line in lines:
        line = line.strip()

        if line == "[SizeStats]":
            in_stats = True
            continue

        if in_stats and line.startswith("[") and line != "[SizeStats]":
            break

        if in_stats and "=" in line:
            key, value = line.split("=", 1)
            stats[key.strip()] = value.strip()

    rows.append({
        "file": file.name,
        "mean_diameter": stats.get("Mean"),
        "median_diameter_D50": stats.get("Median"),
        "mode_diameter": stats.get("MOAOAode"),
        "cell_count": stats.get("SampleSize"),
        "D10": stats.get("D10"),
        "D90": stats.get("D90")
    })

df = pd.DataFrame(rows)

for col in df.columns:
    if col != "file":
        df[col] = pd.to_numeric(df[col], errors="coerce")

df.to_csv("multisizer_summary.csv", index=False)
