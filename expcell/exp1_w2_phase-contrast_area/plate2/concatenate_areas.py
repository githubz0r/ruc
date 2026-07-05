import pandas as pd
import glob
import os

# Find all CSV files in the current directory
csv_files = glob.glob("*.csv")

# Read and combine them
df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)

# Extract well letter and number
df["well_letter"] = df["well"].str[0]
df["well_number"] = df["well"].str[1:].astype(int)

# NaBu concentration
df["mM_NaBu"] = df["well_letter"].map({
    "A": 0,
    "B": 0,
    "C": 1,
    "D": 1
})

# Resveratrol concentration
resv_map = {
    1: 0,
    2: 2.5,
    3: 5,
    4: 10,
    5: 20,
    6: 40
}

df["uM_Resveratrol"] = df["well_number"].map(resv_map)

# Arrange columns
df = df[
    [
        "cell_line",
        "ncell",
        "well",
        "mM_NaBu",
        "uM_Resveratrol",
        "area",
    ]
]

# Sort by well if desired
df = df.sort_values(["well"])

# Save
df.to_csv("combined_cell_area.csv", index=False)

print(f"Combined {len(csv_files)} files.")
print(f"Saved {len(df)} rows to combined_cell_area.csv")
