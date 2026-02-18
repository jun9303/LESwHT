import numpy as np

def funcbody(x, y, z, t):
    """
    Function for an immersed body.
    Returns <= 0 inside body, > 0 outside body.
    Logic preserved from funcbody.f90 (Channel Flow).
    """
    # Channel flow logic: Walls at Y >= 0.5 and Y <= -0.5
    if abs(y) >= 0.5:
        return -1.0
    else:
        return 1.0

def rotate_zaxis(xr, yr, zr, theta):
    """Rotate X and Y coordinate to counterclockwise direction around Z axis."""
    pi = np.pi
    rad = theta * pi / 180.0
    
    xtempo = xr
    ytempo = yr
    
    xr_new = np.cos(rad) * xtempo - np.sin(rad) * ytempo
    yr_new = np.sin(rad) * xtempo + np.cos(rad) * ytempo
    zr_new = zr
    
    return xr_new, yr_new, zr_new

def rotate_yaxis(xr, yr, zr, theta):
    """Rotate X and Z coordinate to counterclockwise direction around Y axis."""
    pi = np.pi
    rad = theta * pi / 180.0
    
    xtempo = xr
    ztempo = zr
    
    xr_new = np.cos(rad) * xtempo - np.sin(rad) * ztempo
    zr_new = -np.sin(rad) * xtempo - np.cos(rad) * ztempo
    yr_new = yr
    
    return xr_new, yr_new, zr_new

def rotate_xaxis(xr, yr, zr, theta):
    """Rotate Y and Z coordinate to counterclockwise direction around X axis."""
    pi = np.pi
    rad = theta * pi / 180.0
    
    ytempo = yr
    ztempo = zr
    
    zr_new = -np.cos(rad) * ztempo + np.sin(rad) * ytempo
    yr_new = np.sin(rad) * ztempo + np.cos(rad) * ytempo
    xr_new = xr
    
    return xr_new, yr_new, zr_new