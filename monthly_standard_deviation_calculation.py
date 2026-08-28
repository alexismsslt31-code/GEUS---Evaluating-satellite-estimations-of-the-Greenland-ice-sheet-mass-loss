import numpy as np
import pandas as pd

# _________________________________________________________________________________________________________________________

#                                                      File organisation
# _________________________________________________________________________________________________________________________


# ── Monthly standard deviation from daily data ─────────────────────────────
"""
def monthly_std(daily_path, columns, min_valid_fraction) -- computes, for a daily AWS CSV, the monthly standard deviation of one or more columns, grouped by real calendar year-month (not just month-of-year), with an optional minimum-valid-fraction filter below which a month's std is set to NaN.
"""


def monthly_std(daily_path, columns=("z_surf_combined", "gps_alt"),
                 min_valid_fraction=0.0):
    """Calcule l'écart-type mensuel de plusieurs séries journalières issues
    d'un fichier AWS, groupées par année-mois réel (contrairement à un simple
    regroupement par mois du calendrier qui fusionnerait janvier 2015 et
    janvier 2016, par exemple).

    input :
        daily_path : str ou Path, chemin vers le CSV journalier (doit contenir
            une colonne 'time' ainsi que les colonnes listées dans `columns`)
        columns : tuple[str], noms des colonnes à traiter
            (par défaut : 'z_surf_combined' et 'gps_alt')
        min_valid_fraction : float entre 0 et 1, fraction minimale de valeurs
            non-NaN requise dans un mois pour calculer un std (sinon NaN).
            0.0 = pas de filtrage (comportement par défaut).

    output :
        pd.DataFrame indexé par mois (Timestamp, début de mois), avec une
        colonne d'écart-type par variable demandée (ex: 'z_surf_combined',
        'gps_alt'). Les NaN sont ignorés dans le calcul (np.nanstd), sauf
        si le mois est vide ou insuffisamment rempli.
    """
    daily_df = pd.read_csv(daily_path, parse_dates=["time"])
    daily_df = daily_df.set_index("time")
    daily_df = daily_df.sort_index()

    def _std_or_nan(group):
        n_valid = group.notna().sum()
        n_total = len(group)
        if n_total == 0:
            return np.nan
        if min_valid_fraction > 0.0 and (n_valid / n_total) < min_valid_fraction:
            return np.nan
        return np.nanstd(group.values)

    results = {}
    for col in columns:
        if col not in daily_df.columns:
            raise KeyError(f"Colonne '{col}' absente de {daily_path}")
        results[col] = daily_df[col].resample("MS").apply(_std_or_nan)

    monthly = pd.DataFrame(results)
    return monthly