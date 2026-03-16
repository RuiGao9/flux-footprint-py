import numpy as np
import sys
import numbers
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from scipy import signal as sg
from numpy import ma

# ===============================================================================
# Defining exception types (merged from Kljun's FFP and FFP_climatology)
# ===============================================================================
exTypes = {'message': 'Message',
           'alert': 'Alert',
           'error': 'Error',
           'fatal': 'Fatal error'}

exceptions = [
    {'code': 1, 'type': exTypes['fatal'], 'msg': 'At least one required parameter is missing.'},
    {'code': 2, 'type': exTypes['error'], 'msg': 'zm (measurement height) must be larger than zero.'},
    {'code': 3, 'type': exTypes['error'], 'msg': 'z0 (roughness length) must be larger than zero.'},
    {'code': 4, 'type': exTypes['error'], 'msg': 'h (BPL height) must be larger than 10 m.'},
    {'code': 5, 'type': exTypes['error'], 'msg': 'zm must be smaller than h (PBL height).'},
    {'code': 6, 'type': exTypes['alert'], 'msg': 'zm should be above roughness sub-layer (12.5*z0).'},
    {'code': 7, 'type': exTypes['error'], 'msg': 'zm/ol must be equal or larger than -15.5.'},
    {'code': 8, 'type': exTypes['error'], 'msg': 'sigmav must be larger than zero.'},
    {'code': 9, 'type': exTypes['error'], 'msg': 'ustar must be >=0.1.'},
    {'code': 10, 'type': exTypes['error'], 'msg': 'wind_dir must be >=0 and <=360.'},
    {'code': 11, 'type': exTypes['fatal'], 'msg': 'Passed data arrays don\'t all have the same length.'},
    {'code': 12, 'type': exTypes['fatal'], 'msg': 'No valid zm passed.'},
    {'code': 13, 'type': exTypes['alert'], 'msg': 'Using z0, ignoring umean if passed.'},
    {'code': 14, 'type': exTypes['alert'], 'msg': 'No valid z0 passed, using umean.'},
    {'code': 15, 'type': exTypes['fatal'], 'msg': 'No valid z0 or umean array passed.'},
    {'code': 16, 'type': exTypes['error'], 'msg': 'At least one required input is invalid. Skipping current footprint.'},
    {'code': 17, 'type': exTypes['alert'], 'msg': 'Only one value of zm passed. Using it for all footprints.'},
    {'code': 18, 'type': exTypes['fatal'], 'msg': 'rs must be a number or a list of numbers.'},
    {'code': 19, 'type': exTypes['alert'], 'msg': 'rs value(s) larger than 90% were found and eliminated.'},
    {'code': 20, 'type': exTypes['error'], 'msg': 'zm must be above roughness sub-layer (12.5*z0).'},
]

def raise_ffp_exception(code, verbosity=2):
    '''Showing abnormal behavior or raising exceptions based on code and verbosity level.'''
    ex = [it for it in exceptions if it['code'] == code][0]
    string = ex['type'] + '(' + str(ex['code']).zfill(4) + '):\n ' + ex['msg']
    if verbosity > 0:
        if ex['type'] == exTypes['fatal']:
            raise Exception(string + '\n Execution aborted.')
        elif verbosity > 1:
            print(string + '\n Execution continues.')

# ===============================================================================
# Functions helping for footprint calculations (FFP and FFP_climatology shared)
# ===============================================================================

def check_ffp_inputs(ustar, sigmav, h, ol, wind_dir, zm, z0, umean, rslayer, verbosity):
    if zm <= 0.: raise_ffp_exception(2, verbosity); return False
    if z0 is not None and umean is None and z0 <= 0.: raise_ffp_exception(3, verbosity); return False
    if h <= 10.: raise_ffp_exception(4, verbosity); return False
    if zm > h: raise_ffp_exception(5, verbosity); return False
    if z0 is not None and umean is None and zm <= 12.5*z0:
        if rslayer == 1: raise_ffp_exception(6, verbosity)
        else: raise_ffp_exception(20, verbosity); return False
    if float(zm)/ol <= -15.5: raise_ffp_exception(7, verbosity); return False
    if sigmav <= 0: raise_ffp_exception(8, verbosity); return False
    if ustar <= 0.1: raise_ffp_exception(9, verbosity); return False
    if wind_dir is not None and (wind_dir > 360 or wind_dir < 0): raise_ffp_exception(10, verbosity); return False
    return True

def get_contour_levels(f, dx, dy, rs=None):
    if rs is None: rs = list(np.linspace(0.10, 0.90, 9))
    if isinstance(rs, (int, float)): rs = [rs]
    pclevs = np.full(len(rs), np.nan)
    ars = np.full(len(rs), np.nan)
    sf = np.sort(f, axis=None)[::-1]
    msf = ma.masked_array(sf, mask=(np.isnan(sf) | np.isinf(sf)))
    csf = msf.cumsum().filled(np.nan) * dx * dy
    for ix, r in enumerate(rs):
        dcsf = np.abs(csf - r)
        pclevs[ix] = sf[np.nanargmin(dcsf)]
        ars[ix] = csf[np.nanargmin(dcsf)]
    return [(round(r, 3), ar, pclev) for r, ar, pclev in zip(rs, ars, pclevs)]

def get_contour_vertices(x, y, f, lev):
    cs = plt.contour(x, y, f, [lev])
    plt.close()
    if not cs.allsegs[0]: return [None, None]
    segs = cs.allsegs[0][0]
    xr, yr = segs[:, 0], segs[:, 1]
    if x.min() >= min(xr) or max(xr) >= x.max() or y.min() >= min(yr) or max(yr) >= y.max():
        return [None, None]
    return [xr, yr]

def plot_footprint(x_2d, y_2d, fs, clevs=None, show_heatmap=True, colormap=None, line_width=0.5):
    if not isinstance(fs, list): fs = [fs]
    if colormap is None: colormap = cm.jet
    fig, ax = plt.subplots(figsize=(10, 8))
    if clevs is not None:
        clevs = [c for c in clevs[::-1] if c is not None]
        for f in fs:
            ax.contour(x_2d, y_2d, f, clevs, colors='w', linewidths=line_width)
    if show_heatmap:
        im = ax.imshow(fs[0], cmap=colormap, extent=(x_2d.min(), x_2d.max(), y_2d.min(), y_2d.max()), origin='lower', aspect=1)
        plt.colorbar(im, shrink=1.0, format='%.3e')
    plt.xlabel('x [m]'); plt.ylabel('y [m]')
    plt.show()
    return fig, ax

# ===============================================================================
# Core function 1: FFP (calculate footprint for a single record)
# ===============================================================================
def FFP(zm=None, z0=None, umean=None, h=None, ol=None, sigmav=None, ustar=None, 
        wind_dir=None, rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], rslayer=0,
        nx=1000, crop=False, fig=False, **kwargs):
    """Single record FFP calculation"""
    #===========================================================================
    # Get kwargs
    show_heatmap = kwargs.get('show_heatmap', True)

    #===========================================================================
    ## Input check
    flag_err = 0
        
    ## Check existence of required input pars
    if None in [zm, h, ol, sigmav, ustar] or (z0 is None and umean is None):
        raise_ffp_exception(1)

    # Define rslayer if not passed
    if rslayer == None: rslayer == 0

    # Define crop if not passed
    if crop == None: crop == 0

    # Define fig if not passed
    if fig == None: fig == 0

    # Check passed values
    if zm <= 0.: raise_ffp_exception(2)
    if z0 is not None and umean is None and z0 <= 0.: raise_ffp_exception(3)
    if h <= 10.: raise_ffp_exception(4)
    if zm > h: raise_ffp_exception(5)        
    if z0 is not None and umean is None and zm <= 12.5*z0:
        if rslayer is 1: raise_ffp_exception(6)
        else: raise_ffp_exception(12)
    if float(zm)/ol <= -15.5: raise_ffp_exception(7)
    if sigmav <= 0: raise_ffp_exception(8)
    if ustar <= 0.1: raise_ffp_exception(9)
    if wind_dir is not None:
        if wind_dir> 360 or wind_dir < 0: raise_ffp_exception(10)
    if nx < 600: raise_ffp_exception(11)

    # Resolve ambiguity if both z0 and umean are passed (defaults to using z0)
    if None not in [z0, umean]: raise_ffp_exception(13)

    #===========================================================================
    # Handle rs
    if rs is not None:

        # Check that rs is a list, otherwise make it a list
        if isinstance(rs, numbers.Number): 
            if 0.9 < rs <= 1 or 90 < rs <= 100: rs = 0.9
            rs = [rs]
        if not isinstance(rs, list): raise_ffp_exception(14)

        # If rs is passed as percentages, normalize to fractions of one
        if np.max(rs) >= 1: rs = [x/100. for x in rs]

        # Eliminate any values beyond 0.9 (90%) and inform user
        if np.max(rs) > 0.9:
            raise_ffp_exception(15)
            rs = [item for item in rs if item <= 0.9]

        # Sort levels in ascending order
        rs = list(np.sort(rs))


    #===========================================================================
    # Model parameters
    a = 1.4524
    b = -1.9914
    c = 1.4622
    d = 0.1359
    ac = 2.17 
    bc = 1.66
    cc = 20.0

    xstar_end = 30
    oln = 5000 #limit to L for neutral scaling
    k = 0.4 #von Karman

    #===========================================================================
    # Scaled X* for crosswind integrated footprint
    xstar_ci_param = np.linspace(d, xstar_end, nx+2)
    xstar_ci_param = xstar_ci_param[1:]

    # Crosswind integrated scaled F* 
    fstar_ci_param = a * (xstar_ci_param-d)**b * np.exp(-c/ (xstar_ci_param-d))
    ind_notnan     = ~np.isnan(fstar_ci_param)
    fstar_ci_param = fstar_ci_param[ind_notnan]
    xstar_ci_param = xstar_ci_param[ind_notnan]

    # Scaled sig_y*
    sigystar_param = ac * np.sqrt(bc * xstar_ci_param**2 / (1 + cc * xstar_ci_param))

    #===========================================================================
    # Real scale x and f_ci
    if z0 is not None:
        # Use z0
        if ol <= 0 or ol >= oln:
            xx  = (1 - 19.0 * zm/ol)**0.25
            psi_f = np.log((1 + xx**2) / 2.) + 2. * np.log((1 + xx) / 2.) - 2. * np.arctan(xx) + np.pi/2
        elif ol > 0 and ol < oln:
            psi_f = -5.3 * zm / ol

        x = xstar_ci_param * zm / (1. - (zm / h)) * (np.log(zm / z0) - psi_f)
        if np.log(zm / z0) - psi_f > 0:
            x_ci = x
            f_ci = fstar_ci_param / zm * (1. - (zm / h)) / (np.log(zm / z0) - psi_f)
        else:
            x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d = None
            flag_err = 1
    else:
        # Use umean if z0 not available
        x = xstar_ci_param * zm / (1. - zm / h) * (umean / ustar * k)
        if umean / ustar > 0:
            x_ci = x
            f_ci = fstar_ci_param / zm * (1. - zm / h) / (umean / ustar * k)
        else:
            x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d = None
            flag_err = 1
                        
    #Maximum location of influence (peak location)
    xstarmax = -c / b + d
    if z0 is not None:
        x_ci_max = xstarmax * zm / (1. - (zm / h)) * (np.log(zm / z0) - psi_f)
    else:
        x_ci_max = xstarmax * zm / (1. - (zm / h)) * (umean / ustar * k)

    #Real scale sig_y
    if abs(ol) > oln:
        ol = -1E6
    if ol <= 0:   #convective
        scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.80
    elif ol > 0:  #stable
        scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.55
    if scale_const > 1:
        scale_const = 1.0
    sigy = sigystar_param / scale_const * zm * sigmav / ustar
    sigy[sigy < 0] = np.nan

    #Real scale f(x,y)
    dx = x_ci[2] - x_ci[1]
    y_pos = np.arange(0, (len(x_ci) / 2.) * dx * 1.5, dx)
    #f_pos = np.full((len(f_ci), len(y_pos)), np.nan)
    f_pos = np.empty((len(f_ci), len(y_pos)))
    f_pos[:] = np.nan
    for ix in range(len(f_ci)):
        f_pos[ix,:] = f_ci[ix] * 1 / (np.sqrt(2 * np.pi) * sigy[ix]) * np.exp(-y_pos**2 / ( 2 * sigy[ix]**2))

    #Complete footprint for negative y (symmetrical)
    y_neg = - np.fliplr(y_pos[None, :])[0]
    f_neg = np.fliplr(f_pos)
    y = np.concatenate((y_neg[0:-1], y_pos))
    f = np.concatenate((f_neg[:, :-1].T, f_pos.T)).T

    # #Matrices for output
    # x_2d = np.tile(x[:,None], (1,len(y)))
    # y_2d = np.tile(y.T,(len(x),1))
    # f_2d = f
    # 1. 首先生成基础网格 (你之前检查过的那段)
    # Matrices for output
    x_2d = np.tile(x[:,None], (1,len(y)))
    y_2d = np.tile(y.T,(len(x),1))
    f_2d = f        

    #===========================================================================
    # Derive footprint ellipsoid incorporating R% of the flux, if requested,
    # starting at peak value.
    dy = dx
    if rs is not None:
        clevs = get_contour_levels(f_2d, dx, dy, rs)
        frs = [item[2] for item in clevs]
        xrs = []
        yrs = []
        for ix, fr in enumerate(frs):
            xr,yr = get_contour_vertices(x_2d, y_2d, f_2d, fr)
            if xr is None:
                frs[ix] = None
            xrs.append(xr)
            yrs.append(yr)
    else:
        if crop:
            rs_dummy = 0.8 #crop to 80%
            clevs = get_contour_levels(f_2d, dx, dy, rs_dummy)
            xrs = []
            yrs = []
            xrs,yrs = get_contour_vertices(x_2d, y_2d, f_2d, clevs[0][2])

    # 2. 紧接着进行旋转 (替换掉原有的旋转代码)
    #===========================================================================
    # 旋转逻辑 (将足迹对准地理风向)
    if wind_dir is not None:
        # 将角度转为弧度
        wind_dir_rad = wind_dir * np.pi / 180.
        
        # 计算极坐标
        dist = np.sqrt(x_2d**2 + y_2d**2)
        angle = np.arctan2(y_2d, x_2d)
        
        # 重新映射坐标 (标准地理/数学坐标映射)
        # 注意：这里的 np.pi/2 偏移量取决于你的风向定义是否 0度为北
        x_2d = dist * np.cos(wind_dir_rad + angle + np.pi/2) 
        y_2d = dist * np.sin(wind_dir_rad + angle + np.pi/2)

        # 同样需要旋转白线(等值线)的顶点坐标，否则白线和阴影会错位
        if rs is not None:
            for ix, r in enumerate(rs):
                xr_lev = np.array([px for px in xrs[ix] if px is not None])    
                yr_lev = np.array([py for py in yrs[ix] if py is not None])    
                dist_r = np.sqrt(xr_lev**2 + yr_lev**2)
                angle_r = np.arctan2(yr_lev, xr_lev)
                
                # 使用相同的逻辑旋转白线顶点
                xrs[ix] = list(dist_r * np.cos(wind_dir_rad + angle_r + np.pi/2))
                yrs[ix] = list(dist_r * np.sin(wind_dir_rad + angle_r + np.pi/2))
                
    #===========================================================================
    # Crop domain and footprint to the largest rs value
    if crop:
        xrs_crop = [x for x in xrs if x is not None]
        yrs_crop = [x for x in yrs if x is not None]
        if rs is not None:
            dminx = np.floor(min(xrs_crop[-1]))
            dmaxx = np.ceil(max(xrs_crop[-1]))
            dminy = np.floor(min(yrs_crop[-1]))
            dmaxy = np.ceil(max(yrs_crop[-1]))
        else:
            dminx = np.floor(min(xrs_crop))
            dmaxx = np.ceil(max(xrs_crop))
            dminy = np.floor(min(yrs_crop))
            dmaxy = np.ceil(max(yrs_crop))
        jrange = np.where((y_2d[0] >= dminy) & (y_2d[0] <= dmaxy))[0]
        jrange = np.concatenate(([jrange[0]-1], jrange, [jrange[-1]+1]))
        jrange = jrange[np.where((jrange>=0) & (jrange<=y_2d.shape[0]-1))[0]]
        irange = np.where((x_2d[:,0] >= dminx) & (x_2d[:,0] <= dmaxx))[0]
        irange = np.concatenate(([irange[0]-1], irange, [irange[-1]+1]))
        irange = irange[np.where((irange>=0) & (irange<=x_2d.shape[1]-1))[0]]
        jrange = [[it] for it in jrange]
        x_2d = x_2d[irange,jrange]
        y_2d = y_2d[irange,jrange]
        f_2d = f_2d[irange,jrange]

    #===========================================================================
    #Rotate 3d footprint if requested
    if wind_dir is not None:			
        wind_dir = wind_dir * np.pi / 180.
        dist = np.sqrt(x_2d**2 + y_2d**2)
        angle = np.arctan2(y_2d, x_2d)
        x_2d = dist * np.sin(wind_dir - angle)
        y_2d = dist * np.cos(wind_dir - angle)

        if rs is not None:
            for ix, r in enumerate(rs):
                xr_lev = np.array([x for x in xrs[ix] if x is not None])    
                yr_lev = np.array([x for x in yrs[ix] if x is not None])    
                dist = np.sqrt(xr_lev**2 + yr_lev**2)
                angle = np.arctan2(yr_lev,xr_lev)
                xr = dist * np.sin(wind_dir - angle)
                yr = dist * np.cos(wind_dir - angle)
                xrs[ix] = list(xr) 
                yrs[ix] = list(yr) 

    #===========================================================================
    # Plot footprint
    if fig:
        fig_out,ax = plot_footprint(x_2d=x_2d, y_2d=y_2d, fs=f_2d,
                                    show_heatmap=show_heatmap,clevs=frs)
        
    #===========================================================================
    # Fill output structure
    if rs is not None:
        return {'x_ci_max': x_ci_max, 'x_ci': x_ci, 'f_ci': f_ci,
                'x_2d': x_2d, 'y_2d': y_2d, 'f_2d': f_2d,
                'rs': rs, 'fr': frs, 'xr': xrs, 'yr': yrs, 'flag_err':flag_err}
    else:
        return {'x_ci_max': x_ci_max, 'x_ci': x_ci, 'f_ci': f_ci,
                'x_2d': x_2d, 'y_2d': y_2d, 'f_2d': f_2d, 'flag_err':flag_err}

# ===============================================================================
# Core function 2: FFP_climatology (batch processing/climatology)
# ===============================================================================
def FFP_climatology(zm=None, z0=None, umean=None, h=None, ol=None, sigmav=None, ustar=None,
                    wind_dir=None, domain=None, dx=None, dy=None, nx=None, ny=None, 
                    rs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], rslayer=0,
                    smooth_data=1, crop=False, pulse=None, verbosity=2, fig=False, **kwargs):
    """Calculation of cumulative FFP for long sequences of data, with options for smoothing and plotting."""
    
    #===========================================================================
    # Get kwargs
    show_heatmap = kwargs.get('show_heatmap', True)


    #===========================================================================
    # Input check
    flag_err = 0
        
    # Check existence of required input pars
    if None in [zm, h, ol, sigmav, ustar] or (z0 is None and umean is None):
        raise_ffp_exception(1, verbosity)

    # Convert all input items to lists
    if not isinstance(zm, list): zm = [zm]
    if not isinstance(h, list): h = [h]
    if not isinstance(ol, list): ol = [ol]
    if not isinstance(sigmav, list): sigmav = [sigmav]
    if not isinstance(ustar, list): ustar = [ustar]
    if not isinstance(wind_dir, list): wind_dir = [wind_dir]
    if not isinstance(z0, list): z0 = [z0]
    if not isinstance(umean, list): umean = [umean]

    # Check that all lists have same length, if not raise an error and exit
    ts_len = len(ustar)
    if any(len(lst) != ts_len for lst in [sigmav, wind_dir, h, ol]):
        # at least one list has a different length, exit with error message
        raise_ffp_exception(11, verbosity)

    # Special treatment for zm, which is allowed to have length 1 for any
    # length >= 1 of all other parameters
    if all(val is None for val in zm): raise_ffp_exception(12, verbosity)
    if len(zm) == 1:
        raise_ffp_exception(17, verbosity)
        zm = [zm[0] for i in range(ts_len)]

    # Resolve ambiguity if both z0 and umean are passed (defaults to using z0)
    # If at least one value of z0 is passed, use z0 (by setting umean to None)
    if not all(val is None for val in z0):
        raise_ffp_exception(13, verbosity)
        umean = [None for i in range(ts_len)]
        # If only one value of z0 was passed, use that value for all footprints
        if len(z0) == 1: z0 = [z0[0] for i in range(ts_len)]
    elif len(umean) == ts_len and not all(val is None for val in umean):
        raise_ffp_exception(14, verbosity)
        z0 = [None for i in range(ts_len)]
    else:
        raise_ffp_exception(15, verbosity)

    # Rename lists as now the function expects time series of inputs
    ustars, sigmavs, hs, ols, wind_dirs, zms, z0s, umeans = \
            ustar, sigmav, h, ol, wind_dir, zm, z0, umean

    #===========================================================================
    # Handle rs
    if rs is not None:

        # Check that rs is a list, otherwise make it a list
        if isinstance(rs, numbers.Number): 
            if 0.9 < rs <= 1 or 90 < rs <= 100: rs = 0.9
            rs = [rs]
        if not isinstance(rs, list): raise_ffp_exception(18, verbosity)

        # If rs is passed as percentages, normalize to fractions of one
        if np.max(rs) >= 1: rs = [x/100. for x in rs]

        # Eliminate any values beyond 0.9 (90%) and inform user
        if np.max(rs) > 0.9:
            raise_ffp_exception(19, verbosity)
            rs = [item for item in rs if item <= 0.9]

        # Sort levels in ascending order
        rs = list(np.sort(rs))

    #===========================================================================
    # Define computational domain
    # Check passed values and make some smart assumptions
    if isinstance(dx, numbers.Number) and dy is None: dy = dx
    if isinstance(dy, numbers.Number) and dx is None: dx = dy
    if not all(isinstance(item, numbers.Number) for item in [dx, dy]): dx = dy = None
    if isinstance(nx, int) and ny is None: ny = nx
    if isinstance(ny, int) and nx is None: nx = ny
    if not all(isinstance(item, int) for item in [nx, ny]): nx = ny = None
    if not isinstance(domain, list) or len(domain) != 4: domain = None

    if all(item is None for item in [dx, nx, domain]):
        # If nothing is passed, default domain is a square of 2 Km size centered
        # at the tower with pizel size of 2 meters (hence a 1000x1000 grid)
        domain = [-1000., 1000., -1000., 1000.]
        dx = dy = 2.
        nx = ny = 1000
    elif domain is not None:
        # If domain is passed, it takes the precendence over anything else
        if dx is not None:
            # If dx/dy is passed, takes precendence over nx/ny
            nx = int((domain[1]-domain[0]) / dx)
            ny = int((domain[3]-domain[2]) / dy)
        else:
            # If dx/dy is not passed, use nx/ny (set to 1000 if not passed)
            if nx is None: nx = ny = 1000
            # If dx/dy is not passed, use nx/ny
            dx = (domain[1]-domain[0]) / float(nx)
            dy = (domain[3]-domain[2]) / float(ny)
    elif dx is not None and nx is not None:
        # If domain is not passed but dx/dy and nx/ny are, define domain
        domain = [-nx*dx/2, nx*dx/2, -ny*dy/2, ny*dy/2]
    elif dx is not None:
        # If domain is not passed but dx/dy is, define domain and nx/ny
        domain = [-1000, 1000, -1000, 1000]
        nx = int((domain[1]-domain[0]) / dx)
        ny = int((domain[3]-domain[2]) / dy)
    elif nx is not None:
        # If domain and dx/dy are not passed but nx/ny is, define domain and dx/dy
        domain = [-1000, 1000, -1000, 1000]
        dx = (domain[1]-domain[0]) / float(nx)
        dy = (domain[3]-domain[2]) / float(nx)

    # Put domain into more convenient vars
    xmin, xmax, ymin, ymax = domain

    # Define rslayer if not passed
    if rslayer == None: rslayer == 0

    # Define smooth_data if not passed
    if smooth_data == None: smooth_data == 1

    # Define crop if not passed
    if crop == None: crop == 0

    # Define pulse if not passed
    if pulse == None:
        if ts_len <= 20:
            pulse = 1
        else:
            pulse = int(ts_len / 20)

    # Define fig if not passed
    if fig == None: fig == 0

    #===========================================================================
    # Model parameters
    a = 1.4524
    b = -1.9914
    c = 1.4622
    d = 0.1359
    ac = 2.17
    bc = 1.66
    cc = 20.0
        
    oln = 5000 #limit to L for neutral scaling
    k = 0.4 #von Karman

    #===========================================================================
    # Define physical domain in cartesian and polar coordinates
    # Cartesian coordinates
    x = np.linspace(xmin, xmax, nx + 1)
    y = np.linspace(ymin, ymax, ny + 1)
    x_2d, y_2d = np.meshgrid(x, y)

    # Polar coordinates
    # Set theta such that North is pointing upwards and angles increase clockwise
    rho = np.sqrt(x_2d**2 + y_2d**2)
    theta = np.arctan2(x_2d, y_2d)

    # initialize raster for footprint climatology
    fclim_2d = np.zeros(x_2d.shape)

    #===========================================================================
    # Loop on time series

    # Initialize logic array valids to those 'timestamps' for which all inputs are
    # at least present (but not necessarily phisically plausible)
    valids = [True if not any([val is None for val in vals]) else False \
              for vals in zip(ustars, sigmavs, hs, ols, wind_dirs, zms)]

    if verbosity > 1: print ('')
    for ix, (ustar, sigmav, h, ol, wind_dir, zm, z0, umean) \
            in enumerate(zip(ustars, sigmavs, hs, ols, wind_dirs, zms, z0s, umeans)):

        # Counter
        if verbosity > 1 and ix % pulse == 0:
            print ('Calculating footprint ', ix+1, ' of ', ts_len)

        valids[ix] = check_ffp_inputs(ustar, sigmav, h, ol, wind_dir, zm, z0, umean, rslayer, verbosity)

        # If inputs are not valid, skip current footprint
        if not valids[ix]:
            raise_ffp_exception(16, verbosity)
        else:
            #===========================================================================
            # Rotate coordinates into wind direction
            if wind_dir is not None:
                rotated_theta = theta - wind_dir * np.pi / 180.

            #===========================================================================
            # Create real scale crosswind integrated footprint and dummy for
            # rotated scaled footprint
            fstar_ci_dummy = np.zeros(x_2d.shape)
            f_ci_dummy = np.zeros(x_2d.shape)
            xstar_ci_dummy = np.zeros(x_2d.shape)
            px = np.ones(x_2d.shape)
            if z0 is not None:
                # Use z0
                if ol <= 0 or ol >= oln:
                    xx = (1 - 19.0 * zm/ol)**0.25
                    psi_f = (np.log((1 + xx**2) / 2.) + 2. * np.log((1 + xx) / 2.) - 2. * np.arctan(xx) + np.pi/2)
                elif ol > 0 and ol < oln:
                    psi_f = -5.3 * zm / ol
                if (np.log(zm / z0)-psi_f)>0:
                    xstar_ci_dummy = (rho * np.cos(rotated_theta) / zm * (1. - (zm / h)) / (np.log(zm / z0) - psi_f))
                    px = np.where(xstar_ci_dummy > d)
                    fstar_ci_dummy[px] = a * (xstar_ci_dummy[px] - d)**b * np.exp(-c / (xstar_ci_dummy[px] - d))
                    f_ci_dummy[px] = (fstar_ci_dummy[px] / zm * (1. - (zm / h)) / (np.log(zm / z0) - psi_f))
                else:
                    flag_err = 3
                    valids[ix] = 0
            else:
                # Use umean if z0 not available
                xstar_ci_dummy = (rho * np.cos(rotated_theta) / zm * (1. - (zm / h)) / (umean / ustar * k))
                px = np.where(xstar_ci_dummy > d)
                fstar_ci_dummy[px] = a * (xstar_ci_dummy[px] - d)**b * np.exp(-c / (xstar_ci_dummy[px] - d))
                f_ci_dummy[px] = (fstar_ci_dummy[px] / zm * (1. - (zm / h)) / (umean / ustar * k))

            #===========================================================================
            # Calculate dummy for scaled sig_y* and real scale sig_y
            sigystar_dummy = np.zeros(x_2d.shape)
            sigystar_dummy[px] = (ac * np.sqrt(bc * np.abs(xstar_ci_dummy[px])**2 / (1 +
                                  cc * np.abs(xstar_ci_dummy[px]))))

            if abs(ol) > oln:
                ol = -1E6
            if ol <= 0:   #convective
                scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.80
            elif ol > 0:  #stable
                scale_const = 1E-5 * abs(zm / ol)**(-1) + 0.55
            if scale_const > 1:
                scale_const = 1.0

            sigy_dummy = np.zeros(x_2d.shape)
            sigy_dummy[px] = (sigystar_dummy[px] / scale_const * zm * sigmav / ustar)
            sigy_dummy[sigy_dummy < 0] = np.nan

            #===========================================================================
            # Calculate real scale f(x,y)
            f_2d = np.zeros(x_2d.shape)
            f_2d[px] = (f_ci_dummy[px] / (np.sqrt(2 * np.pi) * sigy_dummy[px]) *
                        np.exp(-(rho[px] * np.sin(rotated_theta[px]))**2 / ( 2. * sigy_dummy[px]**2)))

            #===========================================================================
            # Add to footprint climatology raster
            fclim_2d = fclim_2d + f_2d;

    #===========================================================================
    # Continue if at least one valid footprint was calculated
    n = sum(valids)
    vs = None
    clevs = None
    if n==0:
        print ("No footprint calculated")
        flag_err = 1
    else:

        #===========================================================================
        # Normalize and smooth footprint climatology
        fclim_2d = fclim_2d / n;

        if smooth_data is not None:
            skernel  = np.matrix('0.05 0.1 0.05; 0.1 0.4 0.1; 0.05 0.1 0.05')
            fclim_2d = sg.convolve2d(fclim_2d,skernel,mode='same');
            fclim_2d = sg.convolve2d(fclim_2d,skernel,mode='same');

        #===========================================================================
        # Derive footprint ellipsoid incorporating R% of the flux, if requested,
        # starting at peak value.
        if rs is not None:
            clevs = get_contour_levels(fclim_2d, dx, dy, rs)
            frs = [item[2] for item in clevs]
            xrs = []
            yrs = []
            for ix, fr in enumerate(frs):
                xr,yr = get_contour_vertices(x_2d, y_2d, fclim_2d, fr)
                if xr is None:
                    frs[ix]  = None
                    flag_err = 2
                xrs.append(xr)
                yrs.append(yr)
        else:
            if crop:
                rs_dummy = 0.8 #crop to 80%
                clevs = get_contour_levels(fclim_2d, dx, dy, rs_dummy)
                xrs = []
                yrs = []
                xrs,yrs = get_contour_vertices(x_2d, y_2d, fclim_2d, clevs[0][2])

        #===========================================================================
        # Crop domain and footprint to the largest rs value
        if crop:
            xrs_crop = [x for x in xrs if x is not None]
            yrs_crop = [x for x in yrs if x is not None]
            if rs is not None:
                dminx = np.floor(min(xrs_crop[-1]))
                dmaxx = np.ceil(max(xrs_crop[-1]))
                dminy = np.floor(min(yrs_crop[-1]))
                dmaxy = np.ceil(max(yrs_crop[-1]))
            else:
                dminx = np.floor(min(xrs_crop))
                dmaxx = np.ceil(max(xrs_crop))
                dminy = np.floor(min(yrs_crop))
                dmaxy = np.ceil(max(yrs_crop))
                
            if dminy>=ymin and dmaxy<=ymax:
                jrange = np.where((y_2d[:,0] >= dminy) & (y_2d[:,0] <= dmaxy))[0]
                jrange = np.concatenate(([jrange[0]-1], jrange, [jrange[-1]+1]))
                jrange = jrange[np.where((jrange>=0) & (jrange<=y_2d.shape[0]))[0]]
            else:
                jrange = np.linspace(0, 1, y_2d.shape[0]-1)
                        
            if dminx>=xmin and dmaxx<=xmax:
                irange = np.where((x_2d[0,:] >= dminx) & (x_2d[0,:] <= dmaxx))[0]
                irange = np.concatenate(([irange[0]-1], irange, [irange[-1]+1]))
                irange = irange[np.where((irange>=0) & (irange<=x_2d.shape[1]))[0]]
            else:
                irange = np.linspace(0, 1, x_2d.shape[1]-1)

            jrange = [[it] for it in jrange]
            x_2d = x_2d[jrange,irange]
            y_2d = y_2d[jrange,irange]
            fclim_2d = fclim_2d[jrange,irange]

            
        #===========================================================================
        # Plot footprint
        if fig:
            fig_out,ax = plot_footprint(x_2d=x_2d, y_2d=y_2d, fs=fclim_2d,
                                        show_heatmap=show_heatmap,clevs=frs)

            
    #===========================================================================
    # Fill output structure
    if rs is not None:
        return {'x_2d': x_2d, 'y_2d': y_2d, 'fclim_2d': fclim_2d,
                'rs': rs, 'fr': frs, 'xr': xrs, 'yr': yrs, 'n':n, 'flag_err':flag_err}
    else:
        return {'x_2d': x_2d, 'y_2d': y_2d, 'fclim_2d': fclim_2d,
                'n':n, 'flag_err':flag_err}
    pass