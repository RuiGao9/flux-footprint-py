# flux_footprint/core.py
from .kljun_ffp import FFP, FFP_climatology

def run_kljun_point(zm, z0, ustar, ol, sigmav, wind_dir, **kwargs):
    # preprocess the input parameters if necessary
    return FFP(zm=zm, z0=z0, ustar=ustar, ol=ol, sigmav=sigmav, wind_dir=wind_dir, **kwargs)

def run_kljun_climatology(df, site_params, **kwargs):
    # convert the DataFrame to the list format required by the Kljun model
    return FFP_climatology(
        zm=site_params['zm'],
        z0=df['z0'].tolist(),
        ustar=df['ustar'].tolist(),
        # other parameters can be added here as needed
        **kwargs
    )