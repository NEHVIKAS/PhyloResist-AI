import pandas as pd
import numpy as np

print("Welcome to PhyloResist AI!")

# Example fungicide sensitivity data
data = {
    "Isolate": ["FG001", "FG002", "FG003", "FG004"],
    "Tebuconazole_EC50": [0.15, 0.32, 0.85, 1.20],
    "Propiconazole_EC50": [0.21, 0.45, 0.91, 1.50]
}

df = pd.DataFrame(data)

print("\nExperimental dataset:")
print(df)

# Log transformation
df["Log10_Tebuconazole"] = np.log10(
    df["Tebuconazole_EC50"]
)

print("\nLog-transformed data:")
print(df)

print("\nProject environment working!")
