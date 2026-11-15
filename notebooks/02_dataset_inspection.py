
import pandas as pd
from pathlib import Path

# Locate the raw data folder
DATA_DIR = Path("data/raw")

# Find all CSV files
files = list(DATA_DIR.glob("*.csv"))

print("\nPHYLORESIST AI - DATASET INSPECTION")
print("=" * 60)

if not files:
    print("No CSV files found!")
    print("Check that your files are inside data/raw/")
else:
    for file in files:

        print("\n" + "=" * 60)
        print("FILE:", file.name)

        # Read the CSV file
        df = pd.read_csv(file)

        # Dataset dimensions
        print("\nDataset dimensions:")
        print("Rows:", df.shape[0])
        print("Columns:", df.shape[1])

        # Column names
        print("\nColumn names:")
        print(df.columns.tolist())

        # First five rows
        print("\nFirst five rows:")
        print(df.head().to_string())

        # Data types
        print("\nData types:")
        print(df.dtypes)

        # Missing values
        print("\nMissing values:")
        print(df.isnull().sum())

        # Duplicate rows
        print("\nDuplicate rows:")
        print(df.duplicated().sum())

print("\nINSPECTION COMPLETED")
