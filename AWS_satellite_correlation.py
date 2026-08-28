from pathlib import Path

import AWS_data
import pandas as pd
import interpolation_altimetry_AWS
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

from paths import FIGURES_DIR as _BASE_FIGURES_DIR


# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Correlation plotting ─────────────────────────────────────────────────
"""
def subplot_correlation(STATION, fig_name, list_satellites, variables, colors) -- builds a 6x6 grid of satellite correlation subplots (incomplete/broken: missing a colon on the def line).
"""

# ── Anomaly helpers ───────────────────────────────────────────────────────
"""
def gps_alt_anomaly(aws_path) -- loads an AWS CSV and adds a column with the GPS altitude anomaly relative to its first valid value.

def z_surf_combined_anomaly(aws) -- adds a column with the combined surface elevation anomaly relative to its first valid value.
"""

# ── Example script (module-level, not a function) ────────────────────────
"""
Module-level code below (outside any function) loads the Copernicus Climate Data Store series matched to the UPE_L station, computes the surface elevation anomaly, computes the Pearson correlation between the AWS and satellite anomalies, and plots a 1:1 scatter of the two.
"""


FIGURES_DIR = _BASE_FIGURES_DIR
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def subplot_correlation(
        STATION,
        fig_name,
        list_satellites,
        variables,
        colors,
)
    
    OUTPUT_FILE = FIGURES_DIR / (fig_name + ".png")
    NCOLS = 6
    NROWS = 6

    fig, axes = plt.subplots(
        NROWS,
        NCOLS,
        sharey=True,
        figsize=(14, NROWS * 3.5),
        constrained_layout=True,
    )




upe_l = interpolation_altimetry_AWS.satellite_on_aws(
    AWS_data.STATION_UPE_L["UPE_L"]["file"], "Copernicus_Climate_Data_Store", "dh"
)


def gps_alt_anomaly(aws_path):
    aws = pd.read_csv(aws_path, parse_dates=["time"])
    aws["gps_alt_anomaly"] = aws["gps_alt"] - aws["gps_alt"].dropna().iloc[0]
    return aws


def z_surf_combined_anomaly(aws):
    aws["z_surf_combined_anomaly"] = (
        aws["z_surf_combined"] - aws["z_surf_combined"].dropna().iloc[0]
    )
    return aws


upe_l_ano = z_surf_combined_anomaly(upe_l)
upe_l_ano = upe_l_ano[["z_surf_combined_anomaly", "dh"]].dropna()

r_pearson, p_pearson = pearsonr(upe_l_ano["z_surf_combined_anomaly"], upe_l_ano["dh"])
print(f"Pearson : r = {r_pearson:.3f}, p = {p_pearson:.3g}")

lims = [
    min(upe_l_ano["z_surf_combined_anomaly"].min(), upe_l_ano["dh"].min()),
    max(upe_l_ano["z_surf_combined_anomaly"].max(), upe_l_ano["dh"].max()),
]

plt.plot(
    lims, lims,
    linestyle='--',
    linewidth=1,
    color="black",
    label="1:1"
)

plt.scatter(
    upe_l_ano["z_surf_combined_anomaly"], upe_l_ano["dh"],
    color="violet",
    s=15,
    alpha=1
)
plt.grid(
    True,
    linestyle="--",
    alpha=0.6
)
plt.title("UPE_L")
plt.xlabel("AWS SEC (m)")
plt.ylabel("Copernicus SEC (m)")
plt.legend()
plt.show()