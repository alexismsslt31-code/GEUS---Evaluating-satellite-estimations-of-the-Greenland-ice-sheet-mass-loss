from paths import AWS_DATA_DIR, SATELLITE_DATA_DIR

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Regional groupings (monthly data, file + daily_data sub-entry) ────────
"""
STATIONS_NW, STATIONS_NE, STATIONS_SE, STATIONS_SW -- every PROMICE/GC-Net station grouped by geographic quadrant of the ice sheet, each entry pointing to both its monthly CSV ("file") and its daily CSV ("daily_data"/"file").
"""

# ── Ablation vs accumulation groupings ─────────────────────────────────────
"""
STATIONS_ABLATION_month, STATIONS_ACCUMULATION, STATIONS_ABLATION_day, STATIONS_ACCUMULATION_day, STATIONS_accumulation (duplicate of STATIONS_ACCUMULATION, redefined further below), STATIONS_ablation -- stations grouped by ablation-zone vs accumulation-zone membership instead of by quadrant, at monthly and daily cadence.
"""

# ── Single-variable selections ─────────────────────────────────────────────
"""
STATIONS_TAS, STATIONS_accumulation, STATIONS_daily, STATIONS_TAS_lat_lon -- narrower dicts keyed by "STATION variable" (e.g. "SDL z_surf_combined") rather than by station alone, each entry carrying an explicit "col" pointing at one AWS column -- used where a plotting call needs to pick a single variable per station rather than the whole record.
"""

# ── Per-region "currently available data" subsets ──────────────────────────
"""
STATIONS_SE_AVAILABLE_month, STATIONS_NE_AVAILABLE_month, STATIONS_NO_AVAILABLE_month, STATIONS_SO_AVAILABLE_month, STATIONS_SE_AVAILABLE_day, STATIONS_NE_AVAILABLE_day, STATIONS_NO_AVAILABLE_day, STATIONS_SO_AVAILABLE_day -- trimmed-down, hand-curated per-quadrant subsets (some stations commented out) reflecting which station records were usable at the time these were built, at monthly and daily cadence.
"""

# ── Small test / single-station selections ──────────────────────────────────
"""
STATIONS_test_ablation, STATIONS_test_accumulation, STATION_TUN, STATION_SWC, STATION_EGP, STATION_NAE, STATION_NSE, STATION_SDL, STATION_SDM, STATION_HUM, STATION_CEN, STATION_NEM, STATION_NAU, STATION_DY2, STATION_KAN_U, STATION_KPC_L, STATION_KPC_U, STATION_TAS_L, STATION_TAS_A, STATION_TAS_U, STATION_QAS_L, STATION_QAS_M, STATION_QAS_A, STATION_QAS_U, STATION_NUK_L, STATION_NUK_U, STATION__NUK_U, STATION_KAN_L, STATION_KAN_M, STATION_KAN_T, STATION_JAR, STATION_UPE_L, STATION_UPE_U, STATION_THU_L -- one-station dicts, used to quickly plot/inspect a single station without building a full regional dict.
"""

# ── Ad hoc working selections ────────────────────────────────────────────────
"""
STATION_use, STATION_sample -- the station subset currently being worked on/plotted in other modules (contents change as the analysis focus moves between stations).
"""


STATIONS_NW = {
    "THU_L": {
        "file": AWS_DATA_DIR / "month/THU_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/THU_L_day.csv"},
    },
    "THU_L2": {
        "file": AWS_DATA_DIR / "month/THU_L2_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/THU_L2_day.csv"},
    },
    "NEM": {
        "file": AWS_DATA_DIR / "month/NEM_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NEM_day.csv"},
    },
    "THU_U": {
        "file": AWS_DATA_DIR / "month/THU_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/THU_U_day.csv"},
    },
    "HUM": {
        "file": AWS_DATA_DIR / "month/HUM_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/HUM_day.csv"},
    },
    "CEN": {
        "file": AWS_DATA_DIR / "month/CEN_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/CEN_day.csv"},
    },
    "UPE_L": {
        "file": AWS_DATA_DIR / "month/UPE_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/UPE_L_day.csv"},
    },
    "UPE_U": {
        "file": AWS_DATA_DIR / "month/UPE_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/UPE_U_day.csv"},
    },
    "NAU": {
        "file": AWS_DATA_DIR / "month/NAU_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NAU_day.csv"},
    },
}

STATIONS_NE = {
    "KPC_L": {
        "file": AWS_DATA_DIR / "month/KPC_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/KPC_L_day.csv"},
    },
    "KPC_U": {
        "file": AWS_DATA_DIR / "month/KPC_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/KPC_U_day.csv"},
    },
    "TUN": {
        "file": AWS_DATA_DIR / "month/TUN_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/TUN_day.csv"},
    },
    "ZAC_U": {
        "file": AWS_DATA_DIR / "month/ZAC_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/ZAC_U_day.csv"},
    },
    "ZAC_L": {
        "file": AWS_DATA_DIR / "month/ZAC_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/ZAC_L_day.csv"},
    },
    "ZAC_A": {
        "file": AWS_DATA_DIR / "month/ZAC_A_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/ZAC_A_day.csv"},
    },
    "EGP": {
        "file": AWS_DATA_DIR / "month/EGP_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/EGP_day.csv"},
    },
    "NAE": {
        "file": AWS_DATA_DIR / "month/NAE_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NAE_day.csv"},
    },
    "SCO_L": {
        "file": AWS_DATA_DIR / "month/SCO_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/SCO_L_day.csv"},
    },
    "SCO_U": {
        "file": AWS_DATA_DIR / "month/SCO_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/SCO_U_day.csv"},
    },
    "FRE": {
        "file": AWS_DATA_DIR / "month/FRE_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/FRE_day.csv"},
    },
}

STATIONS_SE = {
    "NSE": {
        "file": AWS_DATA_DIR / "month/NSE_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NSE_day.csv"},
    },
    "SDM": {
        "file": AWS_DATA_DIR / "month/SDM_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/SDM_day.csv"},
    },
    "MIT": {
        "file": AWS_DATA_DIR / "month/MIT_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/MIT_day.csv"},
    },
    "TAS_U": {
        "file": AWS_DATA_DIR / "month/TAS_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/TAS_U_day.csv"},
    },
    "TAS_L": {
        "file": AWS_DATA_DIR / "month/TAS_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/TAS_L_day.csv"},
    },
    "TAS_A": {
        "file": AWS_DATA_DIR / "month/TAS_A_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/TAS_A_day.csv"},
    },
    "SDL": {
        "file": AWS_DATA_DIR / "month/SDL_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/SDL_day.csv"},
    },
}

STATIONS_SW = {
    "DY2": {
        "file": AWS_DATA_DIR / "month/NSE_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/DY2_day.csv"},
    },
    "KAN_U": {
        "file": AWS_DATA_DIR / "month/KAN_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/KAN_U_day.csv"},
    },
    "KAN_L": {
        "file": AWS_DATA_DIR / "month/KAN_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/KAN_L_day.csv"},
    },
    "KAN_B": {
        "file": AWS_DATA_DIR / "month/KAN_B_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/KAN_B_day.csv"},
    },
    "KAN_T": {
        "file": AWS_DATA_DIR / "month/KAN_T_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/KAN_T_day.csv"},
    },
    "KAN_M": {
        "file": AWS_DATA_DIR / "month/KAN_M_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/KAN_M_day.csv"},
    },
    "NUK_B": {
        "file": AWS_DATA_DIR / "month/NUK_B_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NUK_B_day.csv"},
    },
    "NUK_U": {
        "file": AWS_DATA_DIR / "month/NUK_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NUK_U_day.csv"},
    },
    "NUK_L": {
        "file": AWS_DATA_DIR / "month/NUK_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NUK_L_day.csv"},
    },
    "NUK_K": {
        "file": AWS_DATA_DIR / "month/NUK_K_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/NUK_K_day.csv"},
    },
    "QAS_A": {
        "file": AWS_DATA_DIR / "month/QAS_A_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/QAS_A_day.csv"},
    },
    "QAS_U": {
        "file": AWS_DATA_DIR / "month/QAS_U_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/QAS_U_day.csv"},
    },
    "QAS_L": {
        "file": AWS_DATA_DIR / "month/QAS_L_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/QAS_L_day.csv"},
    },
    "QAS_M": {
        "file": AWS_DATA_DIR / "month/QAS_M_month.csv",
        "daily_data": {"file": AWS_DATA_DIR / "day/QAS_M_day.csv"},
    },
}

STATIONS_ABLATION_month = {
    "RED_L": {"file": AWS_DATA_DIR / "month/RED_L_month.csv", "col": "gps_alt"},
    "THU_L": {"file": AWS_DATA_DIR / "month/THU_L_month.csv", "col": "gps_alt"},
    "THU_L2": {"file": AWS_DATA_DIR / "month/THU_L2_month.csv", "col": "gps_alt"},
    "UPE_L": {"file": AWS_DATA_DIR / "month/UPE_L_month.csv", "col": "gps_alt"},
    "UPE_U": {"file": AWS_DATA_DIR / "month/UPE_U_month.csv", "col": "gps_alt"},
    "WEG_L": {"file": AWS_DATA_DIR / "month/WEG_L_month.csv", "col": "gps_alt"},
    "WEG_B": {"file": AWS_DATA_DIR / "month/WEG_B_month.csv", "col": "gps_alt"},
    "LYN_T": {"file": AWS_DATA_DIR / "month/LYN_T_month.csv", "col": "gps_alt"},
    "LYN_L": {"file": AWS_DATA_DIR / "month/LYN_L_month.csv", "col": "gps_alt"},
    "JAR": {"file": AWS_DATA_DIR / "month/JAR_month.csv", "col": "gps_alt"},
    "SWC": {"file": AWS_DATA_DIR / "month/SWC_month.csv", "col": "gps_alt"},
    "KAN_B": {"file": AWS_DATA_DIR / "month/KAN_B_month.csv", "col": "gps_alt"},
    "KAN_T": {"file": AWS_DATA_DIR / "month/KAN_T_month.csv", "col": "gps_alt"},
    "KAN_L": {"file": AWS_DATA_DIR / "month/KAN_L_month.csv", "col": "gps_alt"},
    "KAN_M": {"file": AWS_DATA_DIR / "month/KAN_M_month.csv", "col": "gps_alt"},
    "NUK_K": {"file": AWS_DATA_DIR / "month/NUK_K_month.csv", "col": "gps_alt"},
    "NUK_B": {"file": AWS_DATA_DIR / "month/NUK_B_month.csv", "col": "gps_alt"},
    "NUK_L": {"file": AWS_DATA_DIR / "month/NUK_L_month.csv", "col": "gps_alt"},
    "NUK_U": {"file": AWS_DATA_DIR / "month/NUK_U_month.csv", "col": "gps_alt"},
    "NUK_N": {"file": AWS_DATA_DIR / "month/NUK_N_month.csv", "col": "gps_alt"},
    "QAS_L": {"file": AWS_DATA_DIR / "month/QAS_L_month.csv", "col": "gps_alt"},
    "QAS_M": {"file": AWS_DATA_DIR / "month/QAS_M_month.csv", "col": "gps_alt"},
    "QAS_U": {"file": AWS_DATA_DIR / "month/QAS_U_month.csv", "col": "gps_alt"},
    "QAS_A": {"file": AWS_DATA_DIR / "month/QAS_A_month.csv", "col": "gps_alt"},
    "TAS_L": {"file": AWS_DATA_DIR / "month/TAS_L_month.csv", "col": "gps_alt"},
    "TAS_U": {"file": AWS_DATA_DIR / "month/TAS_U_month.csv", "col": "gps_alt"},
    "TAS_A": {"file": AWS_DATA_DIR / "month/TAS_A_month.csv", "col": "gps_alt"},
    "SER_B": {"file": AWS_DATA_DIR / "month/SER_B_month.csv", "col": "gps_alt"},
    "MIT": {"file": AWS_DATA_DIR / "month/MIT_month.csv", "col": "gps_alt"},
    "SCO_L": {"file": AWS_DATA_DIR / "month/SCO_L_month.csv", "col": "gps_alt"},
    "SCO_U": {"file": AWS_DATA_DIR / "month/SCO_U_month.csv", "col": "gps_alt"},
    "FRE": {"file": AWS_DATA_DIR / "month/FRE_month.csv", "col": "gps_alt"},
    "ZAC_L": {"file": AWS_DATA_DIR / "month/ZAC_L_month.csv", "col": "gps_alt"},
    "ZAC_U": {"file": AWS_DATA_DIR / "month/ZAC_U_month.csv", "col": "gps_alt"},
    "KPC_L": {"file": AWS_DATA_DIR / "month/KPC_L_month.csv", "col": "gps_alt"},
    "KPC_U": {"file": AWS_DATA_DIR / "month/KPC_U_month.csv", "col": "gps_alt"},
}

STATIONS_ACCUMULATION = {
    "CEN": {
        "file": AWS_DATA_DIR / "month/CEN_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_CEN.csv"
        },
    },
    "HUM": {
        "file": AWS_DATA_DIR / "month/HUM_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_HUM.csv"
        },
    },
    "NEM": {
        "file": AWS_DATA_DIR / "month/NEM_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_NEM.csv"
        },
    },
    "TUN": {
        "file": AWS_DATA_DIR / "month/TUN_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_TUN.csv"
        },
    },
    "EGP": {
        "file": AWS_DATA_DIR / "month/EGP_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_EGP.csv"
        },
    },
    "NAE": {
        "file": AWS_DATA_DIR / "month/NAE_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_NAE.csv"
        },
    },
    "ZAC_A": {
        "file": AWS_DATA_DIR / "month/ZAC_A_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR
            / "Khan_interpolation_AWS/satellite_data_ZAC_A.csv"
        },
    },
    "NAU": {
        "file": AWS_DATA_DIR / "month/NAU_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_NAU.csv"
        },
    },
    "CP1": {
        "file": AWS_DATA_DIR / "month/CP1_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_CP1.csv"
        },
    },
    "KAN_U": {
        "file": AWS_DATA_DIR / "month/KAN_U_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR
            / "Khan_interpolation_AWS/satellite_data_KAN_U.csv"
        },
    },
    "DY2": {
        "file": AWS_DATA_DIR / "month/DY2_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_DY2.csv"
        },
    },
    "NSE": {
        "file": AWS_DATA_DIR / "month/NSE_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_NSE.csv"
        },
    },
    "SDL": {
        "file": AWS_DATA_DIR / "month/SDL_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_SDL.csv"
        },
    },
    "SDM": {
        "file": AWS_DATA_DIR / "month/SDM_month.csv",
        "Khan": {
            "file": SATELLITE_DATA_DIR / "Khan_interpolation_AWS/satellite_data_SDM.csv"
        },
    },
}


STATIONS_ABLATION_day = {
    "RED_L": {"file": AWS_DATA_DIR / "day/RED_L_day.csv", "col": "gps_alt"},
    "THU_L": {"file": AWS_DATA_DIR / "day/THU_L_day.csv", "col": "gps_alt"},
    "THU_L2": {"file": AWS_DATA_DIR / "day/THU_L2_day.csv", "col": "gps_alt"},
    "UPE_L": {"file": AWS_DATA_DIR / "day/UPE_L_day.csv", "col": "gps_alt"},
    "UPE_U": {"file": AWS_DATA_DIR / "day/UPE_U_day.csv", "col": "gps_alt"},
    "WEG_L": {"file": AWS_DATA_DIR / "day/WEG_L_day.csv", "col": "gps_alt"},
    "WEG_B": {"file": AWS_DATA_DIR / "day/WEG_B_day.csv", "col": "gps_alt"},
    "LYN_T": {"file": AWS_DATA_DIR / "day/LYN_T_day.csv", "col": "gps_alt"},
    "LYN_L": {"file": AWS_DATA_DIR / "day/LYN_L_day.csv", "col": "gps_alt"},
    "JAR": {"file": AWS_DATA_DIR / "day/JAR_day.csv", "col": "gps_alt"},
    "SWC": {"file": AWS_DATA_DIR / "day/SWC_day.csv", "col": "gps_alt"},
    "KAN_B": {"file": AWS_DATA_DIR / "day/KAN_B_day.csv", "col": "gps_alt"},
    "KAN_T": {"file": AWS_DATA_DIR / "day/KAN_T_day.csv", "col": "gps_alt"},
    "KAN_L": {"file": AWS_DATA_DIR / "day/KAN_L_day.csv", "col": "gps_alt"},
    "KAN_M": {"file": AWS_DATA_DIR / "day/KAN_M_day.csv", "col": "gps_alt"},
    "NUK_K": {"file": AWS_DATA_DIR / "day/NUK_K_day.csv", "col": "gps_alt"},
    "NUK_B": {"file": AWS_DATA_DIR / "day/NUK_B_day.csv", "col": "gps_alt"},
    "NUK_L": {"file": AWS_DATA_DIR / "day/NUK_L_day.csv", "col": "gps_alt"},
    "NUK_U": {"file": AWS_DATA_DIR / "day/NUK_U_day.csv", "col": "gps_alt"},
    "NUK_N": {"file": AWS_DATA_DIR / "day/NUK_N_day.csv", "col": "gps_alt"},
    "QAS_L": {"file": AWS_DATA_DIR / "day/QAS_L_day.csv", "col": "gps_alt"},
    "QAS_M": {"file": AWS_DATA_DIR / "day/QAS_M_day.csv", "col": "gps_alt"},
    "QAS_U": {"file": AWS_DATA_DIR / "day/QAS_U_day.csv", "col": "gps_alt"},
    "QAS_A": {"file": AWS_DATA_DIR / "day/QAS_A_day.csv", "col": "gps_alt"},
    "TAS_L": {"file": AWS_DATA_DIR / "day/TAS_L_day.csv", "col": "gps_alt"},
    "TAS_U": {"file": AWS_DATA_DIR / "day/TAS_U_day.csv", "col": "gps_alt"},
    "TAS_A": {"file": AWS_DATA_DIR / "day/TAS_A_day.csv", "col": "gps_alt"},
    "SER_B": {"file": AWS_DATA_DIR / "day/SER_B_day.csv", "col": "gps_alt"},
    "MIT": {"file": AWS_DATA_DIR / "day/MIT_day.csv", "col": "gps_alt"},
    "SCO_L": {"file": AWS_DATA_DIR / "day/SCO_L_day.csv", "col": "gps_alt"},
    "SCO_U": {"file": AWS_DATA_DIR / "day/SCO_U_day.csv", "col": "gps_alt"},
    "FRE": {"file": AWS_DATA_DIR / "day/FRE_day.csv", "col": "gps_alt"},
    "ZAC_L": {"file": AWS_DATA_DIR / "day/ZAC_L_day.csv", "col": "gps_alt"},
    "ZAC_U": {"file": AWS_DATA_DIR / "day/ZAC_U_day.csv", "col": "gps_alt"},
    "KPC_L": {"file": AWS_DATA_DIR / "day/KPC_L_day.csv", "col": "gps_alt"},
    "KPC_U": {"file": AWS_DATA_DIR / "day/KPC_U_day.csv", "col": "gps_alt"},
}

STATIONS_ACCUMULATION_day = {
    "CEN": {"file": AWS_DATA_DIR / "day/CEN_day.csv", "col": "gps_alt"},
    "HUM": {"file": AWS_DATA_DIR / "day/HUM_day.csv", "col": "gps_alt"},
    "NEM": {"file": AWS_DATA_DIR / "day/NEM_day.csv", "col": "gps_alt"},
    "TUN": {"file": AWS_DATA_DIR / "day/TUN_day.csv", "col": "gps_alt"},
    "EGP": {"file": AWS_DATA_DIR / "day/EGP_day.csv", "col": "gps_alt"},
    "NAE": {"file": AWS_DATA_DIR / "day/NAE_day.csv", "col": "gps_alt"},
    "ZAC_A": {"file": AWS_DATA_DIR / "day/ZAC_A_day.csv", "col": "gps_alt"},
    "NAU": {"file": AWS_DATA_DIR / "day/NAU_day.csv", "col": "gps_alt"},
    "CP1": {"file": AWS_DATA_DIR / "day/CP1_day.csv", "col": "gps_alt"},
    "KAN_U": {"file": AWS_DATA_DIR / "day/KAN_U_day.csv", "col": "gps_alt"},
    "DY2": {"file": AWS_DATA_DIR / "day/DY2_day.csv", "col": "gps_alt"},
    "NSE": {"file": AWS_DATA_DIR / "day/NSE_day.csv", "col": "gps_alt"},
    "SDL": {"file": AWS_DATA_DIR / "day/SDL_day.csv", "col": "gps_alt"},
    "SDM": {"file": AWS_DATA_DIR / "day/SDM_day.csv", "col": "gps_alt"},
}

STATIONS_TAS = {
    "TAS_L z_surf_combined": {
        "file": AWS_DATA_DIR / "month/TAS_L_month.csv",
        "col": "z_surf_combined",
    },
    "TAS_L gps_alt": {"file": AWS_DATA_DIR / "month/TAS_L_month.csv", "col": "gps_alt"},
    "TAS_U z_surf_combined": {
        "file": AWS_DATA_DIR / "month/TAS_U_month.csv",
        "col": "z_surf_combined",
    },
    "TAS_U gps_alt": {"file": AWS_DATA_DIR / "month/TAS_U_month.csv", "col": "gps_alt"},
    "TAS_A z_surf_combined": {
        "file": AWS_DATA_DIR / "month/TAS_A_month.csv",
        "col": "z_surf_combined",
    },
    "TAS_A gps_alt": {"file": AWS_DATA_DIR / "month/TAS_A_month.csv", "col": "gps_alt"},
}

STATIONS_accumulation = {
    "NSE z_surf_combined": {
        "file": AWS_DATA_DIR / "month/NSE_month.csv",
        "col": "z_surf_combined",
    },
    "NSE gps_alt": {"file": AWS_DATA_DIR / "month/NSE_month.csv", "col": "gps_alt"},
    "SDM z_surf_combined": {
        "file": AWS_DATA_DIR / "month/SDM_month.csv",
        "col": "z_surf_combined",
    },
    "SDM gps_alt": {"file": AWS_DATA_DIR / "month/SDM_month.csv", "col": "gps_alt"},
    "SDL z_surf_combined": {
        "file": AWS_DATA_DIR / "month/SDL_month.csv",
        "col": "z_surf_combined",
    },
    "SDL gps_alt": {"file": AWS_DATA_DIR / "month/SDL_month.csv", "col": "gps_alt"},
}

STATIONS_daily = {
    "NSE z_surf_combined": {
        "file": AWS_DATA_DIR / "day/NSE_day.csv",
        "col": "z_surf_combined",
    },
    "NSE gps_alt": {"file": AWS_DATA_DIR / "day/NSE_day.csv", "col": "gps_alt"},
    "TAS_L z_surf_combined": {
        "file": AWS_DATA_DIR / "day/TAS_L_day.csv",
        "col": "z_surf_combined",
    },
    "TAS_L gps_alt": {"file": AWS_DATA_DIR / "day/TAS_L_day.csv", "col": "gps_alt"},
}

STATIONS_TAS_lat_lon = {
    "TAS_L latitude": {
        "file": AWS_DATA_DIR / "month/TAS_L_month.csv",
        "col": "gps_lat",
    },
    "TAS_L longitude": {
        "file": AWS_DATA_DIR / "month/TAS_L_month.csv",
        "col": "gps_lon",
    },
    "TAS_U latitude": {
        "file": AWS_DATA_DIR / "month/TAS_U_month.csv",
        "col": "gps_lat",
    },
    "TAS_U longitude": {
        "file": AWS_DATA_DIR / "month/TAS_U_month.csv",
        "col": "gps_lon",
    },
    "TAS_A latitude": {
        "file": AWS_DATA_DIR / "month/TAS_A_month.csv",
        "col": "gps_lat",
    },
    "TAS_A longitude": {
        "file": AWS_DATA_DIR / "month/TAS_A_month.csv",
        "col": "gps_lon",
    },
}

STATIONS_SE_AVAILABLE_month = {
    "TAS_A": {"file": AWS_DATA_DIR / "month/TAS_A_month.csv"},
    "TAS_U": {"file": AWS_DATA_DIR / "month/TAS_U_month.csv"},
    "TAS_L": {"file": AWS_DATA_DIR / "month/TAS_L_month.csv"},
    #     "NSE": {"file": AWS_DATA_DIR/ "month/NSE_month.csv"},
    #     "SDL": {"file": AWS_DATA_DIR/ "month/SDL_month.csv"},
    #     "SDM": {"file": AWS_DATA_DIR/ "month/SDM_month.csv"}
}

STATIONS_NE_AVAILABLE_month = {
    #     "TUN": {"file": AWS_DATA_DIR/ "month/TUN_month.csv"},
    #     "EGP": {"file": AWS_DATA_DIR/ "month/EGP_month.csv"},
    #     "NAE": {"file": AWS_DATA_DIR/ "month/NAE_month.csv"},
    #     "ZAC_A": {"file": AWS_DATA_DIR/ "month/ZAC_A_month.csv"},
}

STATIONS_NO_AVAILABLE_month = {
    # "HUM": {"file": AWS_DATA_DIR/ "month/HUM_month.csv"},
    # "CEN": {"file": AWS_DATA_DIR/ "month/CEN_month.csv"},
    # "NEM": {"file": AWS_DATA_DIR/ "month/NEM_month.csv"},
    # "NAU": {"file": AWS_DATA_DIR/ "month/NAU_month.csv"},
    "THU_L": {"file": AWS_DATA_DIR / "month/THU_L_month.csv"},
    "UPE_L": {"file": AWS_DATA_DIR / "month/UPE_L_month.csv"},
    "UPE_U": {"file": AWS_DATA_DIR / "month/UPE_U_month.csv"},
}

STATIONS_SO_AVAILABLE_month = {
    # "DY2": {"file": AWS_DATA_DIR/ "month/DY2_month.csv"},
    "SWC": {"file": AWS_DATA_DIR / "month/SWC_month.csv"},
    "KAN_L": {"file": AWS_DATA_DIR / "month/KAN_L_month.csv"},
    "KAN_M": {"file": AWS_DATA_DIR / "month/KAN_M_month.csv"},
    "KAN_U": {"file": AWS_DATA_DIR / "month/KAN_U_month.csv"},
    "NUK_U": {"file": AWS_DATA_DIR / "month/NUK_U_month.csv"},
    "QAS_L": {"file": AWS_DATA_DIR / "month/QAS_L_month.csv"},
    "QAS_M": {"file": AWS_DATA_DIR / "month/QAS_M_month.csv"},
    "QAS_U": {"file": AWS_DATA_DIR / "month/QAS_U_month.csv"},
}

# _______________________________________________________________________________________________________________________

STATIONS_SE_AVAILABLE_day = {
    "TAS_A": {"file": AWS_DATA_DIR / "day/TAS_A_day.csv"},
    "TAS_U": {"file": AWS_DATA_DIR / "day/TAS_U_day.csv"},
    "TAS_L": {"file": AWS_DATA_DIR / "day/TAS_L_day.csv"},
    # "NSE": {"file": AWS_DATA_DIR/ "day/NSE_day.csv"},
    # "SDL": {"file": AWS_DATA_DIR/ "day/SDL_day.csv"},
    # "SDM": {"file": AWS_DATA_DIR/ "day/SDM_day.csv"}
}

STATIONS_NE_AVAILABLE_day = {
    "TUN": {"file": AWS_DATA_DIR / "day/TUN_day.csv"},
    "EGP": {"file": AWS_DATA_DIR / "day/EGP_day.csv"},
    "NAE": {"file": AWS_DATA_DIR / "day/NAE_day.csv"},
    "ZAC_A": {"file": AWS_DATA_DIR / "day/ZAC_A_day.csv"},
}

STATIONS_NO_AVAILABLE_day = {
    # "HUM": {"file": AWS_DATA_DIR/ "day/HUM_day.csv"},
    # "CEN": {"file": AWS_DATA_DIR/ "day/CEN_day.csv"},
    # "NEM": {"file": AWS_DATA_DIR/ "day/NEM_day.csv"},
    # "NAU": {"file": AWS_DATA_DIR/ "day/NAU_day.csv"},
    "THU_L": {"file": AWS_DATA_DIR / "day/THU_L_day.csv"},
    "UPE_L": {"file": AWS_DATA_DIR / "day/UPE_L_day.csv"},
    "UPE_U": {"file": AWS_DATA_DIR / "day/UPE_U_day.csv"},
}

STATIONS_SO_AVAILABLE_day = {
    # "DY2": {"file": AWS_DATA_DIR/ "day/DY2_day.csv"},
    "SWC": {"file": AWS_DATA_DIR / "day/SWC_day.csv"},
    "KAN_L": {"file": AWS_DATA_DIR / "day/KAN_L_day.csv"},
    "KAN_M": {"file": AWS_DATA_DIR / "day/KAN_M_day.csv"},
    "KAN_U": {"file": AWS_DATA_DIR / "day/KAN_U_day.csv"},
    "NUK_U": {"file": AWS_DATA_DIR / "day/NUK_U_day.csv"},
    "QAS_L": {"file": AWS_DATA_DIR / "day/QAS_L_day.csv"},
    "QAS_M": {"file": AWS_DATA_DIR / "day/QAS_M_day.csv"},
    "QAS_U": {"file": AWS_DATA_DIR / "day/QAS_U_day.csv"},
}


STATIONS_test_ablation = {
    "QAS_U": {"file": AWS_DATA_DIR / "month/QAS_U_month.csv"},
    "NUK_U": {"file": AWS_DATA_DIR / "month/NUK_U_month.csv"},
}

STATIONS_test_accumulation = {
    "NAE": {"file": AWS_DATA_DIR / "month/NAE_month.csv"},
    "TAS_A": {"file": AWS_DATA_DIR / "month/TAS_A_month.csv"},
}

STATIONS_accumulation = {
    "TUN": {
        "file": AWS_DATA_DIR / "month/TUN_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/TUN_hour.csv",
    },
    "EGP": {
        "file": AWS_DATA_DIR / "month/EGP_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/EGP_hour.csv",
    },
    "NAE": {
        "file": AWS_DATA_DIR / "month/NAE_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NAE_hour.csv",
    },
    "NSE": {
        "file": AWS_DATA_DIR / "month/NSE_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NSE_hour.csv",
    },
    "SDL": {
        "file": AWS_DATA_DIR / "month/SDL_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/SDL_hour.csv",
    },
    "SDM": {
        "file": AWS_DATA_DIR / "month/SDM_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/SDM_hour.csv",
    },
    "HUM": {
        "file": AWS_DATA_DIR / "month/HUM_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/HUM_hour.csv",
    },
    "CEN": {
        "file": AWS_DATA_DIR / "month/CEN_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/CEN_hour.csv",
    },
    "NEM": {
        "file": AWS_DATA_DIR / "month/NEM_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NEM_hour.csv",
    },
    "NAU": {
        "file": AWS_DATA_DIR / "month/NAU_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NAU_hour.csv",
    },
    "DY2": {
        "file": AWS_DATA_DIR / "month/DY2_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/DY2_hour.csv",
    },
    "KAN_U": {
        "file": AWS_DATA_DIR / "month/KAN_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_U_hour.csv",
    },
}


STATIONS_ablation = {
    "KPC_L": {
        "file": AWS_DATA_DIR / "month/KPC_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KPC_L_hour.csv",
    },
    "KPC_U": {
        "file": AWS_DATA_DIR / "month/KPC_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KPC_U_hour.csv",
    },
    "TAS_L": {
        "file": AWS_DATA_DIR / "month/TAS_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/TAS_L_hour.csv",
    },
    "TAS_A": {
        "file": AWS_DATA_DIR / "month/TAS_A_month.csv",
        "hourly_data": AWS_DATA_DIR / "day/TAS_A_day.csv",
    },
    "TAS_U": {
        "file": AWS_DATA_DIR / "month/TAS_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/TAS_U_hour.csv",
    },
    "QAS_L": {
        "file": AWS_DATA_DIR / "month/QAS_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_L_hour.csv",
    },
    "QAS_M": {
        "file": AWS_DATA_DIR / "month/QAS_M_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_M_hour.csv",
    },
    "QAS_A": {
        "file": AWS_DATA_DIR / "month/QAS_A_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_A_hour.csv",
    },
    "QAS_U": {
        "file": AWS_DATA_DIR / "month/QAS_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_U_hour.csv",
    },
    "NUK_L": {
        "file": AWS_DATA_DIR / "month/NUK_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NUK_L_hour.csv",
    },
    "NUK_U": {
        "file": AWS_DATA_DIR / "month/NUK_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NUK_U_hour.csv",
    },
    "KAN_L": {
        "file": AWS_DATA_DIR / "month/KAN_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_L_hour.csv",
    },
    "KAN_M": {
        "file": AWS_DATA_DIR / "month/KAN_M_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_M_hour.csv",
    },

    "KAN_T": {
        "file": AWS_DATA_DIR / "month/KAN_T_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_T_hour.csv",
    },
    "SWC": {
        "file": AWS_DATA_DIR / "month/SWC_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/SWC_hour.csv",
    },
    "JAR": {
        "file": AWS_DATA_DIR / "month/JAR_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/JAR_hour.csv",
    },
    "UPE_L": {
        "file": AWS_DATA_DIR / "month/UPE_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/UPE_L_hour.csv",
    },
    "UPE_U": {
        "file": AWS_DATA_DIR / "month/UPE_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/UPE_U_hour.csv",
    },
    "THU_L": {
        "file": AWS_DATA_DIR / "month/THU_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/THU_L_hour.csv",
    },
}


STATION_TUN = {
    "TUN": {
        "file": AWS_DATA_DIR / "month/TUN_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/TUN_hour.csv",
    }
}

STATION_SWC = {
    "SWC": {
        "file": AWS_DATA_DIR / "month/SWC_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/SWC_hour.csv",
    }
}

STATION_EGP = {
    "EGP": {
        "file": AWS_DATA_DIR / "month/EGP_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/EGP_hour.csv",
    }
}

STATION_NAE = {
    "NAE": {
        "file": AWS_DATA_DIR / "month/NAE_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NAE_hour.csv",
    }
}

STATION_NSE = {
    "NSE": {
        "file": AWS_DATA_DIR / "month/NSE_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NSE_hour.csv",
    }
}

STATION_SDL = {
    "SDL": {
        "file": AWS_DATA_DIR / "month/SDL_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/SDL_hour.csv",
    }
}

STATION_SDM = {
    "SDM": {
        "file": AWS_DATA_DIR / "month/SDM_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/SDM_hour.csv",
    }
}

STATION_HUM = {
    "HUM": {
        "file": AWS_DATA_DIR / "month/HUM_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/HUM_hour.csv",
    }
}

STATION_CEN = {
    "CEN": {
        "file": AWS_DATA_DIR / "month/CEN_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/CEN_hour.csv",
    }
}

STATION_NEM = {
    "NEM": {
        "file": AWS_DATA_DIR / "month/NEM_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NEM_hour.csv",
    }
}

STATION_NAU = {
    "NAU": {
        "file": AWS_DATA_DIR / "month/NAU_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NAU_hour.csv",
    }
}

STATION_DY2 = {
    "DY2": {
        "file": AWS_DATA_DIR / "month/DY2_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/DY2_hour.csv",
    }
}

STATION_KAN_U = {
    "KAN_U": {
        "file": AWS_DATA_DIR / "month/KAN_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_U_hour.csv",
    }
}

STATION_KPC_L = {
    "KPC_L": {
        "file": AWS_DATA_DIR / "month/KPC_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KPC_L_hour.csv",
    }
}

STATION_KPC_U = {
    "KPC_U": {
        "file": AWS_DATA_DIR / "month/KPC_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KPC_U_hour.csv",
    }
}

STATION_TAS_L = {
    "TAS_L": {
        "file": AWS_DATA_DIR / "month/TAS_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/TAS_L_hour.csv",
    }
}

STATION_TAS_A = {
    "TAS_A": {
        "file": AWS_DATA_DIR / "month/TAS_A_month.csv",
        "hourly_data": AWS_DATA_DIR / "day/TAS_A_day.csv",
    }
}

STATION_TAS_U = {
    "TAS_U": {
        "file": AWS_DATA_DIR / "month/TAS_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/TAS_U_hour.csv",
    }
}

STATION_QAS_L = {
    "QAS_L": {
        "file": AWS_DATA_DIR / "month/QAS_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_L_hour.csv",
    }
}

STATION_QAS_M = {
    "QAS_M": {
        "file": AWS_DATA_DIR / "month/QAS_M_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_M_hour.csv",
    }
}

STATION_QAS_A = {
    "QAS_A": {
        "file": AWS_DATA_DIR / "month/QAS_A_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_A_hour.csv",
    }
}

STATION_QAS_U = {
    "QAS_U": {
        "file": AWS_DATA_DIR / "month/QAS_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/QAS_U_hour.csv",
    }
}

STATION_NUK_L = {
    "NUK_L": {
        "file": AWS_DATA_DIR / "month/NUK_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NUK_L_hour.csv",
    }
}

STATION_NUK_U = {
    "NUK_U": {
        "file": AWS_DATA_DIR / "month/NUK_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NUK_U_hour.csv",
    }
}

STATION__NUK_U = {
    "NUK_U": {
        "file": AWS_DATA_DIR / "month/NUK_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/NUK_U_hour.csv",
    }
}

STATION_KAN_L = {
    "KAN_L": {
        "file": AWS_DATA_DIR / "month/KAN_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_L_hour.csv",
    }
}

STATION_KAN_M = {
    "KAN_M": {
        "file": AWS_DATA_DIR / "month/KAN_M_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_M_hour.csv",
    }
}

STATION_KAN_T = {
    "KAN_T": {
        "file": AWS_DATA_DIR / "month/KAN_T_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/KAN_T_hour.csv",
    }
}

STATION_SWC = {
    "SWC": {
        "file": AWS_DATA_DIR / "month/SWC_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/SWC_hour.csv",
    }
}

STATION_JAR = {
    "JAR": {
        "file": AWS_DATA_DIR / "month/JAR_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/JAR_hour.csv",
    }
}

STATION_UPE_L = {
    "UPE_L": {
        "file": AWS_DATA_DIR / "month/UPE_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/UPE_L_hour.csv",
    }
}

STATION_UPE_U = {
    "UPE_U": {
        "file": AWS_DATA_DIR / "month/UPE_U_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/UPE_U_hour.csv",
    }
}

STATION_THU_L = {
    "THU_L": {
        "file": AWS_DATA_DIR / "month/THU_L_month.csv",
        "hourly_data": AWS_DATA_DIR / "hour/THU_L_hour.csv",
    }
}

STATION_use = {
    "SDL": {
            "file": AWS_DATA_DIR / "month/SDL_month.csv",
            "hourly_data": AWS_DATA_DIR / "hour/SDL_hour.csv",
        },
    "NSE": {
            "file": AWS_DATA_DIR / "month/NSE_month.csv",
            "hourly_data": AWS_DATA_DIR / "hour/NSE_hour.csv",
        },
}

STATION_sample = {
    "EGP": {
                "file": AWS_DATA_DIR / "month/EGP_month.csv",
                "hourly_data": AWS_DATA_DIR / "hour/EGP_hour.csv",
            },
    "KAN_U": {
                "file": AWS_DATA_DIR / "month/KAN_U_month.csv",
                "hourly_data": AWS_DATA_DIR / "hour/KAN_U_hour.csv",
            },
    "CEN": {
                "file": AWS_DATA_DIR / "month/CEN_month.csv",
                "hourly_data": AWS_DATA_DIR / "hour/CEN_hour.csv",
            },
    "ZAC_L": {
                "file": AWS_DATA_DIR / "month/ZAC_L_month.csv",
                "hourly_data": AWS_DATA_DIR / "hour/ZAC_L_hour.csv",
            },
    "THU_L": {
                "file": AWS_DATA_DIR / "month/THU_L_month.csv",
                "hourly_data": AWS_DATA_DIR / "hour/THU_L_hour.csv",
            },
    "SWC": {
                "file": AWS_DATA_DIR / "month/SWC_month.csv",
                "hourly_data": AWS_DATA_DIR / "hour/SWC_hour.csv",
            },
}