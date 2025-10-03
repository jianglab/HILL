import numpy as np
import pandas as pd
import os, pathlib
from joblib import Memory

username = os.getenv("USER")
memory = Memory(location=f"/tmp/{username}_joblib_cache", verbose=0)

from shiny import reactive
from shiny.express import ui, render, module, expressify
import plotly.graph_objects as go
import helicon
from itertools import product

from finufft import nufft2d2
from numba import set_num_threads, prange

from scipy.ndimage import map_coordinates
import scipy.fft
import scipy.fftpack as fp
from scipy.spatial.transform import Rotation as R
from scipy.ndimage import affine_transform, map_coordinates
from scipy.signal import correlate
from scipy.interpolate import splrep, splev
from scipy.interpolate import RegularGridInterpolator
from scipy.special import jnp_zeros
from scipy.optimize import minimize, fmin

from streamlit_bokeh import streamlit_bokeh
from bokeh.events import MouseMove, MouseEnter, DoubleTap
from bokeh.io import export_png
from bokeh.layouts import gridplot, column, layout
from bokeh.models import Button, ColumnDataSource, CustomJS, Label, LinearColorMapper, Slider, Span, Spinner
from bokeh.models.tools import CrosshairTool, HoverTool
from bokeh.plotting import figure

def create_image_figure(image, dx, dy, title="", title_location="below", plot_width=None, plot_height=None, x_axis_label='x', y_axis_label='y', tooltips=None, show_axis=True, show_toolbar=True, crosshair_color="white", aspect_ratio=None):
    #from bokeh.plotting import figure
    h, w = image.shape
    if aspect_ratio is None:
        if plot_width and plot_height:
            aspect_ratio = plot_width/plot_height
        else:
            aspect_ratio = w*dx/(h*dy)
    tools = 'box_zoom,crosshair,pan,reset,save,wheel_zoom,tap'
    fig = figure(title_location=title_location, 
        frame_width=plot_width, frame_height=plot_height, 
        x_axis_label=x_axis_label, y_axis_label=y_axis_label,
        x_range=(-w//2*dx, (w//2-1)*dx), y_range=(-h//2*dy, (h//2-1)*dy), 
        tools=tools, aspect_ratio=aspect_ratio)
    fig.grid.visible = False
    if title:
        fig.title.text=title
        fig.title.align = "center"
        fig.title.text_font_size = "18px"
        fig.title.text_font_style = "normal"
    if not show_axis: fig.axis.visible = False
    if not show_toolbar: fig.toolbar_location = None

    source_data = ColumnDataSource(data=dict(image=[image], x=[-w//2*dx], y=[-h//2*dy], dw=[w*dx], dh=[h*dy]))
    #from bokeh.models import LinearColorMapper
    color_mapper = LinearColorMapper(palette='Greys256')    # Greys256, Viridis256
    image = fig.image(source=source_data, image='image', color_mapper=color_mapper,
                x='x', y='y', dw='dw', dh='dh'
            )

    # add hover tool only for the image
    #from bokeh.models.tools import HoverTool, CrosshairTool
    if not tooltips:
        tooltips = [("x", "$x Å"), ('y', '$y Å'), ('val', '@image')]
    image_hover = HoverTool(renderers=[image], tooltips=tooltips)
    fig.add_tools(image_hover)
    fig.hover[0].attachment="vertical"
    crosshair = [t for t in fig.tools if isinstance(t, CrosshairTool)]
    if crosshair: 
        for ch in crosshair: ch.line_color = crosshair_color
    return fig, source_data

def create_layerline_image_figure(data, cutoff_res_x, cutoff_res_y, helical_radius, tilt, phase=None, fft_top_only=False, pseudo_color=True, const_image_color="", title="", yaxis_visible=True, tooltips=None, hover_phase=False):
    ny, nx = data.shape
    dsy = 1/(ny//2*cutoff_res_y)
    dsx = 1/(nx//2*cutoff_res_x)
    x_range = (-(nx//2+0.5)*dsx, (nx//2-0.5)*dsx)
    if fft_top_only:
        y_range = (-(ny//2 * 0.01)*dsy, (ny//2-0.5)*dsy)
    else:
        y_range = (-(ny//2+0.5)*dsy, (ny//2-0.5)*dsy)

    bessel = bessel_n_image(ny, nx, cutoff_res_x, cutoff_res_y, helical_radius, tilt).astype(np.int16)

    tools = 'box_zoom,pan,reset,save,wheel_zoom'
    fig = figure(title_location="below", frame_width=nx, frame_height=ny, 
        x_axis_label=None, y_axis_label=None, x_range=x_range, y_range=y_range, tools=tools, active_drag="box_zoom")
    fig.grid.visible = False
    fig.title.text = title
    fig.title.align = "center"
    fig.title.text_font_size = "20px"
    fig.yaxis.visible = yaxis_visible   # leaving yaxis on will make the crosshair x-position out of sync with other figures

    source_data = ColumnDataSource(data=dict(image=[data.astype(np.float16)], x=[-nx//2*dsx], y=[-ny//2*dsy], dw=[nx*dsx], dh=[ny*dsy], bessel=[bessel]))
    if phase is not None: source_data.add(data=[np.fmod(np.rad2deg(phase)+360, 360).astype(np.float16)], name="phase")
    if const_image_color:
        palette = (const_image_color,)
    else:
        palette = 'Viridis256' if pseudo_color else 'Greys256'
    color_mapper = LinearColorMapper(palette=palette)    # Greys256, Viridis256
    image = fig.image(source=source_data, image='image', color_mapper=color_mapper, x='x', y='y', dw='dw', dh='dh')
    if tooltips is None:
        tooltips = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Val', '@image')]
    #if phase is not None: tooltips.append(("Phase", "@phase °"))
    if hover_phase: tooltips.append(("Phase", "@phase °"))
    image_hover = HoverTool(renderers=[image], tooltips=tooltips, attachment="vertical")
    fig.add_tools(image_hover)

    # avoid the need for embedding resr/resy/resx image -> smaller fig object and less data to transfer
    mousemove_callback_code = """
    var x = cb_obj.x
    var y = cb_obj.y
    var resr = Math.round((1./Math.sqrt(x*x + y*y) + Number.EPSILON) * 100) / 100
    var resy = Math.abs(Math.round((1./y + Number.EPSILON) * 100) / 100)
    var resx = Math.abs(Math.round((1./x + Number.EPSILON) * 100) / 100)
    hover.tooltips[0][1] = resr.toString() + " Å"
    hover.tooltips[1][1] = resy.toString() + " Å"
    hover.tooltips[2][1] = resx.toString() + " Å"
    """
    mousemove_callback = CustomJS(args={"hover":fig.hover[0]}, code=mousemove_callback_code)
    fig.js_on_event(MouseMove, mousemove_callback)
    
    return fig, source_data

def resize_rescale_power_spectra(data, nyquist_res, cutoff_res=None, output_size=None, log=True, low_pass_fraction=0, high_pass_fraction=0, norm=1):
    #from scipy.ndimage import map_coordinates
    ny, nx = data.shape
    ony, onx = output_size
    res_y, res_x = cutoff_res
    Y, X = np.meshgrid(np.arange(ony, dtype=np.float32)-(ony//2+0.5), np.arange(onx, dtype=np.float32)-(onx//2+0.5), indexing='ij')
    Y = Y/(ony//2+0.5) * nyquist_res/res_y * ny//2 + ny//2+0.5
    X = X/(onx//2+0.5) * nyquist_res/res_x * nx//2 + nx//2+0.5
    pwr = map_coordinates(data, (Y.flatten(), X.flatten()), order=3, mode='constant').reshape(Y.shape)
    if log: pwr = np.log1p(np.abs(pwr))
    if 0<low_pass_fraction<1 or 0<high_pass_fraction<1:
        pwr = low_high_pass_filter(pwr, low_pass_fraction=low_pass_fraction, high_pass_fraction=high_pass_fraction)
    if norm: pwr = normalize(pwr, percentile=(0, 100))
    return pwr

def compute_power_spectra(data, apix, cutoff_res=None, output_size=None, log=True, low_pass_fraction=0, high_pass_fraction=0):
    fft = fft_rescale(data, apix=apix, cutoff_res=cutoff_res, output_size=output_size)
    fft = np.fft.fftshift(fft)  # shift fourier origin from corner to center

    if log: pwr = np.log1p(np.abs(fft))
    else: pwr = np.abs(fft)
    if 0<low_pass_fraction<1 or 0<high_pass_fraction<1:
        pwr = low_high_pass_filter(pwr, low_pass_fraction=low_pass_fraction, high_pass_fraction=high_pass_fraction)
    pwr = normalize(pwr, percentile=(0, 100))

    phase = np.angle(fft, deg=False)
    return pwr, phase

def fft_rescale(image, apix=1.0, cutoff_res=None, output_size=None):
    if cutoff_res:
        res_limit_y, cutoff_res_x = cutoff_res
    else:
        res_limit_y, cutoff_res_x = 2*apix, 2*apix
    if output_size:
        ony, onx = output_size
    else:
        ony, onx = image.shape
    freq_y = np.fft.fftfreq(ony) * 2*apix/res_limit_y
    freq_x = np.fft.fftfreq(onx) * 2*apix/cutoff_res_x
    Y, X = np.meshgrid(freq_y, freq_x, indexing='ij')
    Y = (2*np.pi * Y).flatten(order='C')
    X = (2*np.pi * X).flatten(order='C')

    #from finufft import nufft2d2
    fft = nufft2d2(x=Y, y=X, f=image.astype(np.complex128), eps=1e-6)
    fft = fft.reshape((ony, onx))

    # phase shifts for real-space shifts by half of the image box in both directions
    phase_shift = np.ones(fft.shape)
    phase_shift[1::2, :] *= -1
    phase_shift[:, 1::2] *= -1
    fft *= phase_shift
    # now fft has the same layout and phase origin (i.e. np.fft.ifft2(fft) would obtain original image)
    return fft

def auto_correlation(data, sqrt=True, high_pass_fraction=0):
    #from scipy.signal import correlate2d
    fft = np.fft.rfft2(data)
    product = fft*np.conj(fft)
    if sqrt: product = np.sqrt(product)
    if 0<high_pass_fraction<=1:
        ny, nx = product.shape
        Y, X = np.meshgrid(np.arange(-ny//2, ny//2, dtype=float), np.arange(-nx//2, nx//2, dtype=float), indexing='ij')
        Y /= ny//2
        X /= nx//2
        f2 = np.log(2)/(high_pass_fraction**2)
        filter = 1.0 - np.exp(- f2 * Y**2) # Y-direction only
        product *= np.fft.fftshift(filter)
    corr = np.fft.fftshift(np.fft.irfft2(product))
    corr /= np.max(corr)
    return corr

def low_high_pass_filter(data, low_pass_fraction=0, high_pass_fraction=0):
    fft = np.fft.fft2(data)
    ny, nx = fft.shape
    Y, X = np.meshgrid(np.arange(ny, dtype=np.float32)-ny//2, np.arange(nx, dtype=np.float32)-nx//2, indexing='ij')
    Y /= ny//2
    X /= nx//2
    if 0<low_pass_fraction<1:
        f2 = np.log(2)/(low_pass_fraction**2)
        filter_lp = np.exp(- f2 * (X**2+Y**2))
        fft *= np.fft.fftshift(filter_lp)
    if 0<high_pass_fraction<1:
        f2 = np.log(2)/(high_pass_fraction**2)
        filter_hp = 1.0 - np.exp(- f2 * (X**2+Y**2))
        fft *= np.fft.fftshift(filter_hp)
    ret = np.abs(np.fft.ifft2(fft))
    return ret

def normalize(data, percentile=(0, 100)):
    p0, p1 = percentile
    vmin, vmax = sorted(np.percentile(data, (p0, p1)))
    data2 = (data-vmin)/(vmax-vmin)
    return data2

def bessel_1st_peak_positions(n_max:int = 100):
    #import numpy as np
    ret = np.zeros(n_max+1, dtype=np.float32)
    #from scipy.special import jnp_zeros
    for i in range(1, n_max+1):
        ret[i] = jnp_zeros(i, 1)[0]
    return ret

def bessel_n_image(ny, nx, nyquist_res_x, nyquist_res_y, radius, tilt):
    #import numpy as np
    table = bessel_1st_peak_positions()
    
    if tilt:
        dsx = 1./(nyquist_res_x*nx//2)
        dsy = 1./(nyquist_res_x*ny//2)
        Y, X = np.meshgrid(np.arange(ny, dtype=np.float32)-ny//2, np.arange(nx, dtype=np.float32)-nx//2, indexing='ij')
        Y = 2*np.pi * np.abs(Y)*dsy * radius
        X = 2*np.pi * np.abs(X)*dsx * radius
        Y /= np.cos(np.deg2rad(tilt))
        X = np.hypot(X, Y*np.sin(np.deg2rad(tilt)))
        X = np.expand_dims(X.flatten(), axis=-1)
        indices = np.abs(table - X).argmin(axis=-1)
        return np.reshape(indices, (ny, nx)).astype(np.int16)
    else:
        ds = 1./(nyquist_res_x*nx//2)
        xs = 2*np.pi * np.abs(np.arange(nx)-nx//2)*ds * radius
        xs = np.expand_dims(xs, axis=-1)
        indices = np.abs(table - xs).argmin(axis=-1)
        return np.tile(indices, (ny, 1)).astype(np.int16)

def simulate_helix(twist, rise, csym, helical_radius, ball_radius, ny, nx, apix, tilt=0, az0=None):
    def simulate_projection(centers, sigma, ny, nx, apix):
        sigma2 = sigma*sigma
        d = np.zeros((ny, nx))
        Y, X = np.meshgrid(np.arange(0, ny, dtype=np.float32)-ny//2, np.arange(0, nx, dtype=np.float32)-nx//2, indexing='ij')
        X *= apix
        Y *= apix
        for ci in range(len(centers)):
            yc, xc = centers[ci]
            x = X-xc
            y = Y-yc
            d += np.exp(-(x*x+y*y)/sigma2)
        return d
    def helical_unit_positions(twist, rise, csym, radius, height, tilt=0, az0=0):
        imax = int(height/rise)
        i0 = -imax
        i1 = imax
        
        centers = np.zeros(((2*imax+1)*csym, 3), dtype=np.float32)
        for i in range(i0, i1+1):
            z = rise*i
            for si in range(csym):
                angle = np.deg2rad(twist*i + si*360./csym + az0 + 90)   # start from +y axis
                x = np.cos(angle) * radius
                y = np.sin(angle) * radius
                centers[i*csym+si, 0] = x
                centers[i*csym+si, 1] = y
                centers[i*csym+si, 2] = z
        if tilt:
            #from scipy.spatial.transform import Rotation as R
            rot = R.from_euler('x', tilt, degrees=True)
            centers = rot.apply(centers)
        centers = centers[:, [2, 0]]    # project along y
        return centers
    if az0 is None: az0 = np.random.uniform(0, 360)
    centers = helical_unit_positions(twist, rise, csym, helical_radius, height=ny*apix, tilt=tilt, az0=az0)
    projection = simulate_projection(centers, ball_radius, ny, nx, apix)
    return projection

@helicon.cache(
    cache_dir=str(helicon.cache_dir / "hill"), expires_after=7, verbose=0
)  # 7 days
def symmetrize_transform_map(
    data,
    apix,
    twist_degree,
    rise_angstrom,
    csym=1,
    fraction=1.0,
    new_size=None,
    new_apix=None,
    axial_rotation=0,
    tilt=0,
):
    if new_apix > apix:
        data_work = helicon.low_high_pass_filter(
            data, low_pass_fraction=apix / new_apix
        )
    else:
        data_work = data
    m = helicon.apply_helical_symmetry(
        data=data_work,
        apix=apix,
        twist_degree=twist_degree,
        rise_angstrom=rise_angstrom,
        csym=csym,
        new_size=new_size,
        new_apix=new_apix,
        fraction=fraction,
        cpu=helicon.available_cpu(),
    )
    if axial_rotation or tilt:
        m = helicon.transform_map(m, rot=axial_rotation, tilt=tilt)
    return m


def compute_layer_line_positions(twist, rise, csym, radius, tilt, cutoff_res, m_max=-1):
    if cutoff_res<=0:
        return []
    table = bessel_1st_peak_positions()/(2*np.pi*radius)

    if m_max<1:
        m_max = int(np.floor(np.abs(rise/cutoff_res)))+3
    m = list(range(-m_max, m_max+1))
    m.sort(key=lambda x: (abs(x), x))   # 0, -1, 1, -2, 2, ...
    
    smax = 1./cutoff_res

    tf = 1./np.cos(np.deg2rad(tilt))
    tf2 = np.sin(np.deg2rad(tilt))
    m_groups = {} # one group per m order
    for mi in range(len(m)):
        d = {}
        sy0 = m[mi] / rise

        # first peak positions of each layer line
        p = twist2pitch(twist, rise)
        ds_p = 1/p
        ll_i_top = int(np.abs(smax - sy0)/ds_p) * 2
        ll_i_bottom = -int(np.abs(-smax - sy0)/ds_p) * 2
        ll_i = np.array([i for i in range(ll_i_bottom, ll_i_top+1) if not i%csym], dtype=np.int32)
        sy = sy0 + ll_i * ds_p
        sx = table[np.clip(np.abs(ll_i), 0, len(table)-1)]
        if tilt:
            sy = np.array(sy, dtype=np.float32) * tf
            sx = np.sqrt(np.power(np.array(sx, dtype=np.float32), 2) - np.power(sy*tf2, 2))
            sx[np.isnan(sx)] = 1e-6
        px  = list(sx) + list(-sx)
        py  = list(sy) + list(sy)
        n = list(ll_i) + list(ll_i)
        d["LL"] = (px, py, n)
        d["m"] = m

        m_groups[m[mi]] = d
    return m_groups

def compute_phase_difference_across_meridian(phase):
    # https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html
    phase_diff = phase * 0
    phase_diff[..., 1:] = phase[..., 1:] - phase[..., 1:][..., ::-1]
    phase_diff = np.rad2deg(np.arccos(np.cos(phase_diff)))   # set the range to [0, 180]. 0 -> even order, 180 - odd order
    return phase_diff

@helicon.cache(expires_after=7, cache_dir=helicon.cache_dir / "hill", verbose=0)
def transform_2d_filament(data, angle, dx, dy, negate, apix):
    if negate:
        data = -data

    if (angle or dx or dy):
        data= rotate_shift_image(
            data, 
            angle = -angle, 
            post_shift = (dy / apix, dx / apix), 
            order = 1
        )
    return data

def rotate_shift_image(data, angle=0, pre_shift=(0, 0), post_shift=(0, 0), rotation_center=None, order=1):
    # pre_shift/rotation_center/post_shift: [y, x]
    if angle==0 and pre_shift==[0,0] and post_shift==[0,0]: return data*1.0
    ny, nx = data.shape
    if rotation_center is None:
        rotation_center = np.array((ny//2, nx//2), dtype=np.float32)
    ang = np.deg2rad(angle)
    m = np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]], dtype=np.float32)
    pre_dy, pre_dx = pre_shift    
    post_dy, post_dx = post_shift

    offset = -np.dot(m, np.array([post_dy, post_dx], dtype=np.float32).T) # post_rotation shift
    offset += np.array(rotation_center, dtype=np.float32).T - np.dot(m, np.array(rotation_center, dtype=np.float32).T)  # rotation around the specified center
    offset += -np.array([pre_dy, pre_dx], dtype=np.float32).T     # pre-rotation shift

    #from scipy.ndimage import affine_transform
    ret = affine_transform(data, matrix=m, offset=offset, order=order, mode='constant')
    return ret

@helicon.cache(expires_after=7, cache_dir=helicon.cache_dir / "hill", verbose=0)
def mask_2d_filament(data, mask_radius, apix, mask_len_fraction):
    _, nx = data.shape
    fraction_x = mask_radius/(nx//2*apix)
    tapering_image = generate_tapering_filter(image_size=data.shape, fraction_start=[mask_len_fraction, fraction_x], fraction_slope=(1.0-mask_len_fraction)/2.)
    return data*tapering_image

def get_class2d_from_uploaded_file(fileobj):
    import os, tempfile

    orignal_filename = fileobj.name
    print("file obj name: ", orignal_filename)
    # add back!! os.listdir("C:\\Users\\anika\\AppData\\Local\\Temp\\")
    suffix = os.path.splitext(orignal_filename)[-1]
    with tempfile.NamedTemporaryFile(suffix=suffix) as temp:
        temp.write(fileobj.read())
        print("TN: ", temp.name)
        return get_class2d_from_file(temp.name)

def generate_tapering_filter(image_size, fraction_start=[0, 0], fraction_slope=0.1):
    ny, nx = image_size
    fy, fx = fraction_start
    if not (0<fy<1 or 0<fx<1): return np.ones((ny, nx))
    Y, X = np.meshgrid(np.arange(0, ny, dtype=np.float32)-ny//2, np.arange(0, nx, dtype=np.float32)-nx//2, indexing='ij')
    filter = np.ones_like(Y)
    if 0<fy<1:
        Y = np.abs(Y / (ny//2))
        inner = Y<fy
        outer = Y>fy+fraction_slope
        Y = (Y-fy)/fraction_slope
        Y = (1. + np.cos(Y*np.pi))/2.0
        Y[inner]=1
        Y[outer]=0
        filter *= Y
    if 0<fx<1:
        X = np.abs(X / (nx//2))
        inner = X<fx
        outer = X>fx+fraction_slope
        X = (X-fx)/fraction_slope
        X = (1. + np.cos(X*np.pi))/2.0
        X[inner]=1
        X[outer]=0
        filter *= X
    return filter

def estimate_radial_range(data, thresh_ratio=0.1):
    proj_y = np.sum(data, axis=0)
    n = len(proj_y)
    background = np.mean(proj_y[[0,1,2,-3,-2,-1]])
    thresh = (proj_y.max() - background) * thresh_ratio + background
    indices = np.nonzero(proj_y<thresh)[0]
    try:
        xmin = np.max(indices[indices<np.argmax(proj_y[:n//2])])
    except:
        xmin = 0
    try:
        xmax = np.min(indices[indices>np.argmax(proj_y[n//2:])+n//2])
    except:
        xmax = n-1
    mask_radius = max(abs(n//2-xmin), abs(xmax-n//2))
    proj_y -= thresh
    proj_y[proj_y<0] = 0
    def fitRadialProfile(x, radProfile):
        a, b, w, rcore, rmax= x  # y = a*(sqrt(rmax^2-x^2)+(w-1)*sqrt(rcore^2-x^2))+b
        try:
            n = len(radProfile)
            x = np.abs(np.arange(n, dtype=float)-n/2)
            yshell = radProfile * 0
            mask = x<=abs(rmax)
            yshell[mask] = np.sqrt(rmax*rmax - x[mask]*x[mask])
            ycore = radProfile * 0
            mask = x<=abs(rcore)
            ycore[mask] = np.sqrt(rcore*rcore - x[mask]*x[mask])
            y = a*(yshell+(w-1)*ycore)+b
            score = np.linalg.norm(y-radProfile)
        except:
            score = 1e10
        return score
    #from scipy.optimize import minimize
    #from itertools import product
    bounds = ((0, None), (None, None), (0, None), (0, mask_radius), (0, mask_radius))
    vals_a = (1, 2, 4, 8)
    vals_w = (0, 0.5)
    vals_rcore = (0, mask_radius/2)
    results = []
    for val_a, val_w, val_rcore in product(vals_a, vals_w, vals_rcore):
        x0 = (val_a, 0, val_w, val_rcore, mask_radius)
        res = minimize(fitRadialProfile, x0, args=(proj_y,), method='Nelder-Mead', bounds=bounds, tol=1e-6)
        a, b, w, rcore, rmax = res.x
        results.append((res.fun, w, rcore, rmax, val_a, val_w, val_rcore))
    result = sorted(results)[0]
    w, rcore, rmax = result[1:4]
    rmean = 0.5 * (rmax*rmax+(w-1)*rcore*rcore) / (rmax+(w-1)*rcore)
    return float(rmean), float(mask_radius)    # pixel

@memory.cache
def get_class2d_from_url(url):

    print("start function")
    url_final = get_direct_url(url)  # convert cloud drive indirect url to direct url
    print("url final: ", url_final)
    fileobj = download_file_from_url(url_final)
    print("fileobj: ", fileobj)
    if fileobj is None:
        raise ValueError(
            f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview"
        )
    print("before data")
    data, apix, crs = get_class2d_from_file(fileobj.name)
    print("data: ", data)
    return data, apix, crs


def get_class2d_from_file(classFile):
    import mrcfile
    #print("before mrcfile")
    #classFile = "C:\\Users\\anika\\Downloads\\run_it020_classes.mrcs"
    try: 
        open(classFile, "r")
    except Exception as e:
        print("open exception: ", e)

    #print("after opening")
    with mrcfile.open(classFile) as mrc:
        #print("opens mrc file")
        map_crs = [int(mrc.header.mapc), int(mrc.header.mapr), int(mrc.header.maps)]
        apix = float(mrc.voxel_size.x)
        data = mrc.data
    #print("got data: ", data)
    return data, round(apix, 4), map_crs

def download_file_from_url(url):
    import tempfile
    import requests
    import os

    if pathlib.Path(url).is_file():
        return open(url, "rb")
    try:
        filesize = get_file_size(url)
        local_filename = url.split("/")[-1]
        suffix = "." + local_filename
        fileobj = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
        with requests.get(url) as r:
            r.raise_for_status()  # Check for request success
            fileobj.write(r.content)
        return fileobj
    except requests.exceptions.RequestException as e:
        print(e)
        return None


def get_direct_url(url):
    import re

    if url.startswith("https://drive.google.com/file/d/"):
        hash = url.split("/")[5]
        return f"https://drive.google.com/uc?export=download&id={hash}"
    elif url.startswith("https://app.box.com/s/"):
        hash = url.split("/")[-1]
        return f"https://app.box.com/shared/static/{hash}"
    elif url.startswith("https://www.dropbox.com"):
        if url.find("dl=1") != -1:
            return url
        elif url.find("dl=0") != -1:
            return url.replace("dl=0", "dl=1")
        else:
            return url + "?dl=1"
    elif url.find("sharepoint.com") != -1 and url.find("guestaccess.aspx") != -1:
        return url.replace("guestaccess.aspx", "download.aspx")
    elif url.startswith("https://1drv.ms"):
        import base64

        data_bytes64 = base64.b64encode(bytes(url, "utf-8"))
        data_bytes64_String = (
            data_bytes64.decode("utf-8").replace("/", "_").replace("+", "-").rstrip("=")
        )
        return (
            f"https://api.onedrive.com/v1.0/shares/u!{data_bytes64_String}/root/content"
        )
    else:
        return url


def get_file_size(url):
    import requests

    response = requests.head(url)
    if "Content-Length" in response.headers:
        file_size = int(response.headers["Content-Length"])
        return file_size
    else:
        return None

def extract_emdb_id(url):
    import re
    pattern = r'EMD-(\d+)'
    match = re.search(pattern, url)
    if match:
        return f"EMD-{match.group(1)}"
    return None


def extract_url(input_string):
    # Check if the input string contains a URL in square brackets with parentheses format
    import re
    
    # Pattern to match [text](url) format
    pattern = r'\[.*?\]\((.*?)\)'
    
    match = re.search(pattern, input_string)
    if match:
        return match.group(1)
    else:
        return None


def get_images_from_url(url):
    url_final = helicon.get_direct_url(url)  # convert cloud drive indirect url to direct url
    fileobj = helicon.download_file_from_url(url_final)
    if fileobj is None:
        raise ValueError(
            f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview"
        )
    data, apix = get_images_from_file(fileobj.name)
    return data, apix


def get_images_from_file(imageFile):
    import mrcfile

    with mrcfile.open(imageFile) as mrc:
        apix = float(mrc.voxel_size.x)
        data = mrc.data
    return data, round(apix, 4)

def get_amyloid_n_sub_1_symmetry(twist, rise, max_n=10):
    ret = 1
    for n in range(max_n, 1, -1):
        if not (4.5 < rise * n < 5):
            continue
        if abs(360 - abs(twist * n)) > 90:
            continue
        ret = n
        break
    return ret 

def get_one_map_xyz_projects(map_info): # length_z, map_projection_xyz_choices):
    label = map_info.label
    try:
        data, apix = map_info.get_data()
    except Exception as e:
        if map_info.filename:
            msg = f"Failed to obtain uploaded map {label}"
        elif map_info.url:
            msg = f"Failed to download the map from {map_info.url}"
        elif map_info.emd_id:
            msg = f"Failed to download the map from EMDB for {map_info.emd_id}"
        raise ValueError(msg)
    
    images = []
    image_labels = []
    # if 'z' in map_projection_xyz_choices:
    #     rise = map_info.rise
    #     if rise>0:
    #         rise *=  get_amyloid_n_sub_1_symmetry(twist=map_info.twist, rise=map_info.rise)
    #         images += [helicon.crop_center_z(data, n=max(1, int(0.5 + length_z * rise / apix))).sum(axis=0)]
    #     else:
    #         images += [data.sum(axis=0)]
    #     image_labels += [label + ':Z']
    images += [data.sum(axis=1)]
    image_labels += [label + ':Y']
    # if 'x' in map_projection_xyz_choices:
    #     images += [data.sum(axis=2)]
    #     image_labels += [label + ':X']
    print("data, apix: ", data, apix)
    print("sum: ", np.sum(data))
        
    return images, image_labels, apix

def plot_2d_class(micrograph, title, apix, plot_height=None, plot_width=None):

    fig = go.FigureWidget()

    h, w = micrograph.shape

    fig.add_trace(
        go.Heatmap(
            name="image",
            z=micrograph,
            x=np.arange(w) * apix,
            y=np.arange(h) * apix,
            colorscale="gray",
            showscale=False,  # Removed colorbar
            #reversescale=True, # Greys colorscale is white to black by default??? why?
            hoverongaps=False,
            hovertemplate=(
                "x: %{x:.1f} Å<br>"
                + "y: %{y:.1f} Å<br>"
                + "val: %{z:.1f}<br>"
                + "<extra></extra>"
            ),
        )
    )

    layout_params = {
        "title": {
            "text": title,
            "x": 0.5,
            "y": 0.95,
            "xanchor": "center",
            "font": {"size": 14},
        },
        "xaxis": {"visible": False, "range": [0, w * apix]},
        "yaxis": {
            "visible": False,
            "range": [0, h * apix],
            "scaleanchor": "x",
            "autorange": "reversed",
        },
        "plot_bgcolor": "white",
        "showlegend": False,
        "width": plot_width,
        "height": plot_height,
    }

    if plot_width or plot_height:
        if plot_width:
            layout_params["width"] = plot_width
        if plot_height:
            layout_params["height"] = plot_height
    else:
        layout_params["autosize"] = True

    layout_params["margin"] = dict(l=0, r=0, t=50, b=0)

    layout_params["modebar"] = {
        "remove": [],
        "bgcolor": "rgba(255, 255, 255, 0.7)",
    }

    fig.update_layout(**layout_params)

    return fig

def symmetrize_project_align_one_map(map_info, image_query, image_query_label, image_query_apix, rescale_apix, length_xy_factor, match_sf, angle_range, scale_range):
    if abs(map_info.twist) < 1e-3:
        return map_info, None
    
    try:
        data, apix = map_info.get_data()
    except:
        return map_info, None

    twist = map_info.twist
    rise = map_info.rise
    csym = map_info.csym
    label = map_info.label
    
    nz, ny, nx = data.shape
    if rescale_apix:
        image_ny, image_nx = image_query.shape
        new_apix = image_query_apix
        twist_work = helicon.set_to_periodic_range(twist, min=-180, max=180)
        if abs(twist_work)<90:
            pitch = 360/abs(twist_work) * rise 
        elif abs(twist_work)<180:
            pitch = 360/(180-abs(twist_work)) * rise
        else:
            pitch = image_nx * new_apix
        length = int(pitch / new_apix + image_nx * length_xy_factor)//2*2
        new_size = (length, image_ny, image_ny)

        data_work = helicon.low_high_pass_filter(data, low_pass_fraction=apix/new_apix)
    else:
        new_apix = apix
        new_size = (nz, ny, nx)
        data_work = data

    fraction = 5 * rise / (nz * apix)
    
    data_sym = helicon.apply_helical_symmetry(
        data = data_work,
        apix = apix,
        twist_degree = twist,
        rise_angstrom = rise,
        csym = csym,
        fraction = fraction,
        new_size = new_size,
        new_apix = new_apix,
        cpu = helicon.available_cpu()
    )
    proj = data_sym.sum(axis=2).T
        
    flip, scale, rotation_angle, shift_cartesian, similarity_score, aligned_image_moving = helicon.align_images(image_moving=image_query, image_ref=proj, scale_range=scale_range, angle_range=angle_range, check_polarity=True, check_flip=True, return_aligned_moving_image=True) 

    if match_sf:
        mask = aligned_image_moving > 0
        proj = helicon.match_structural_factors(data=proj, apix=new_apix, data_target=aligned_image_moving, apix_target=new_apix, mask=mask)

    return map_info, (flip, scale, rotation_angle, shift_cartesian, similarity_score, aligned_image_moving, image_query_label, proj, label)


def set_to_periodic_range(v, min=-180, max=180):
    if min <= v <= max: return v
    #from math import fmod
    tmp = fmod(v-min, max-min)
    if tmp>=0: tmp+=min
    else: tmp+=max
    return tmp

def twist2pitch(twist, rise):
    if twist:
        return 360. * rise/abs(twist)
    else:
        return rise

def pitch2twist(pitch, rise):
    if pitch>rise:
        return set_to_periodic_range(360. * rise/pitch)
    else:
        return 0.

# def auto_vertical_center(data, n_theta=180):
#   #from skimage.transform import radon
#   #from scipy.signal import correlate
  
#   data_work = np.clip(data, 0, None)
  
#   theta = np.linspace(start=0., stop=180., num=n_theta, endpoint=False)
#   #import warnings
#   with warnings.catch_warnings(): # ignore outside of circle warnings
#     warnings.simplefilter('ignore')
#     sinogram = radon(data_work, theta=theta)
#   sinogram += sinogram[::-1, :]
#   y = np.std(sinogram, axis=0)
#   theta_best = -theta[np.argmax(y)]

#   rotated_data = rotate_shift_image(data_work, angle=theta_best)
#   # now find best vertical shift
#   yproj = np.sum(rotated_data, axis=0)
#   yproj_xflip = yproj*1.0
#   yproj_xflip[1:] = yproj[1:][::-1]
#   corr = correlate(yproj, yproj_xflip, mode='same')
#   shift_best = -(np.argmax(corr) - len(corr)//2)/2

#   # refine to sub-degree, sub-pixel level
#   def score_rotation_shift(x):
#     theta, shift_x = x
#     data_tmp=rotate_shift_image(data_work, angle=theta, post_shift=(0, shift_x), order=1)
#     xproj = np.sum(data_tmp, axis=0)[1:]
#     xproj += xproj[::-1]
#     score = -np.std(xproj)
#     return score
#   #from scipy.optimize import fmin
#   res = fmin(score_rotation_shift, x0=(theta_best, shift_best), xtol=1e-2, disp=0)
#   theta_best, shift_best = res
#   return set_to_periodic_range(theta_best), shift_best


#@st.cache_data(persist='disk', show_spinner=False)
def gen_filament_template(length, diameter, angle=0, center_offset=(0, 0), image_size=(1024, 1024), apix=1.0, order=5):
    ny, nx = image_size
    y = (np.arange(0, ny) - ny//2)*apix
    x = (np.arange(0, nx) - nx//2)*apix
    Y, X = np.meshgrid(y, x, indexing='ij')
    # flattop gaussian: order>2
    d = np.exp( -np.log(2)*(np.abs(np.power((Y)/(length/2), order))+np.abs(np.power((X)/(diameter/2), order))) )
    if angle!=0:
        #from skimage import transform
        d = transform.rotate(image=d, angle=angle, center=(nx//2, ny//2))
    if center_offset!=(0, 0):
        #from skimage import transform
        xform = transform.EuclideanTransform(
            translation = (center_offset[0]/apix, center_offset[1]/apix)
        )
        d = transform.warp(d, xform.inverse)
    return d

#@st.cache_data(persist='disk', show_spinner=False)
def pad_to_size(array, ny, nx):
    h, w = array.shape
    a = (ny - h) // 2
    aa = ny - a - h
    b = (nx - w) // 2
    bb = nx - b - w
    return np.pad(array, pad_width=((a, aa), (b, bb)), mode='constant')

#@st.cache_data(persist='disk', show_spinner=False)
def filament_transform_fft(image, filament_template, angle_step):
    #import scipy.fft
    #import skimage.transform
    ny, nx = image.shape
    fny, fnx = filament_template.shape
    if ny != fny or nx != fnx:
        pad = True
    else:
        pad = False
    angles = np.arange(0, 180, angle_step)
    res_cc = np.zeros(shape=(len(angles), ny, nx), dtype=np.float32)
    res_ang = np.zeros(shape=(len(angles), ny, nx), dtype=np.float32)
    image_fft = scipy.fft.rfft2(image)
    for ai, angle in enumerate(angles):
        template = transform.rotate(image=filament_template, angle=angle, center=(fnx//2, fny//2))
        if pad:
            template = pad_to_size(template, ny, nx)
        template_fft = np.conj(scipy.fft.rfft2(template))
        res_cc[ai] = scipy.fft.fftshift(scipy.fft.irfft2(image_fft * template_fft))
    fft1d = scipy.fft.fft(res_cc, axis=0)
    fft1d_abs = np.abs(fft1d)
    ret_amp2f = fft1d_abs[1, :, :]/np.sum(fft1d_abs, axis=0)
    ret_ang = np.rad2deg(np.angle(fft1d)[1, :, :])
    ret_ang[ret_ang<0] += 360
    return (ret_amp2f, ret_ang)

#@st.cache_data(persist='disk', show_spinner=False)
def sample_axis_dots(data, apix, nx, ny, r_filament_pixel, l_template_pixel, da, num_samples, lp_x, lp_y):
    # fill in potential black backgrounds with helical boxer
    #data_slice_median=np.median(data)
    #for i in range(ny):
    #    for j in range(nx):
    #        if data[i,j]==0:
    #            data[i,j]=data_slice_median
    #        else:
    #            break
    #    for j in range(nx):
    #        if data[i,-(j+1)]==0:
    #            data[i,-(j+1)]=data_slice_median
    #        else:
    #            break
    
    if nx < 2 * 2 * r_filament_pixel:
        data = pad_to_size(data, 2 * 2 * r_filament_pixel,ny)
        nx = 2 * 2 * r_filament_pixel


    # apply low pass filter
    data_fft=fp.fftshift(fp.fft2(data))

    #kernel = Gaussian2DKernel(lp_x,lp_y,0,x_size=nx,y_size=ny).array
    kernel = gen_filament_template(length=lp_y, diameter=lp_x, image_size=(ny, nx), apix=apix, order=2)
    max_k=np.max(kernel)
    min_k=np.min(kernel)
    kernel=(kernel-min_k)/(max_k-min_k)
    kernel_shape=np.shape(kernel)

    data_fft_filtered=np.multiply(data_fft,kernel)

    data_filtered=fp.ifft2(fp.ifftshift(data_fft_filtered)).real

    # normalize
    #data_filtered=(data_filtered-np.mean(data_filtered))/np.std(data_filtered)
    vmin = data_filtered.min()
    vmax = data_filtered.max()
    data_filtered = (vmax-data_filtered)/(vmax-vmin)
    
    diameter = 2*r_filament_pixel*apix
    length = l_template_pixel

    template_size = round(max(diameter, length)/apix*1.2)//2*2

    filament_template = gen_filament_template(length=length, diameter=diameter, image_size=(np.min([template_size,ny]), np.min([template_size,nx])), apix=apix, order=2)
    filament_transform_method = filament_transform_fft
    cc, ang = filament_transform_method(image=data_filtered, filament_template=filament_template, angle_step=3)
    cc_vmin = cc.min()
    cc_vmax = cc.max()
    cc = (cc_vmax-cc)/(cc_vmax-cc_vmin)
    cc_template = np.repeat([np.mean(cc,axis=0)],repeats=length,axis=0)
    cc, ang = filament_transform_method(image=cc, filament_template=cc_template, angle_step=da)

    ####################################################################
    # center point detection
    dots=np.zeros(np.shape(data))
    centers=cc.argmax(axis=1)

    xs=[]
    ys=[]
    row_offset=int(ny/num_samples/2)
    #num_samples=10
    for i in range(row_offset,ny,int(ny/num_samples)):
        xs.append(centers[i])
        ys.append(i)
    #xs.append(centers[-1])
    #ys.append(ny-1)


    xs=np.array(xs)
    ys=np.array(ys)

    return xs, ys

def create_fit_spline_figure(data,xs,ys,new_xs,apix):
    h, w = data.shape
    xs = (np.array(xs)-w//2)*apix
    new_xs = (np.array(new_xs)-w//2)*apix
    ys = (np.array(ys)-h//2)*apix

    aspect_ratio = w/h
    tools = 'box_zoom,crosshair,pan,reset,save,wheel_zoom'
    fig = figure(x_range=(-w//2*apix, (w//2-1)*apix), y_range=(-h//2*apix, (h//2-1)*apix), 
        tools=tools, aspect_ratio=aspect_ratio)
    fig.grid.visible = False
    fig.axis.visible = False
    fig.toolbar_location = None

    source_data = ColumnDataSource(data=dict(image=[data], x=[-w//2*apix], y=[-h//2*apix], dw=[w*apix], dh=[h*apix]))
    color_mapper = LinearColorMapper(palette='Greys256')    # Greys256, Viridis256
    data = fig.image(source=source_data, image='image', color_mapper=color_mapper,
                x='x', y='y', dw='dw', dh='dh'
            )

    # add hover tool only for the image
    #from bokeh.models.tools import HoverTool, CrosshairTool
    tooltips = [("x", "$xÅ"), ('y', '$yÅ'), ('val', '@image')]
    image_hover = HoverTool(renderers=[data], tooltips=tooltips)
    fig.add_tools(image_hover)
    fig.hover[0].attachment="vertical"
       
    # Plot the line and points
    source_points = ColumnDataSource(data=dict(xs=xs, ys=ys))
    source_line = ColumnDataSource(data=dict(new_xs=new_xs, ys=ys))
    fig.line('new_xs', 'ys', source=source_line, line_color="red")
    fig.circle('xs', 'ys', source=source_points, color="red", size=5)
    
    return fig

#@st.cache_data(persist='disk', show_spinner=False)
def fit_spline(_disp_col,data,xs,ys,apix,display=False):
    # fit spline
    ny,nx=data.shape
    tck = splrep(ys,xs,s=20)
    new_xs = splev(ys,tck)
    
    return new_xs,tck

# test the straightening part by forcing using the center of a straight filament: works
# TODO: check the nans in the output images
#@st.cache_data(persist='disk', show_spinner=False)
def filament_straighten(_disp_col,data,tck,new_xs,ys,r_filament_pixel_display,apix):
    ny,nx=data.shape
    # resample pixels
    #st.info(tck)
    y0=0
    x0=splev(y0,tck)
    for i in range(ny):
        dxdy=splev(y0,tck,der=1)
        orthog_dxdy=-(1.0/dxdy)
        tangent_x0y0=lambda y: dxdy*y + (x0-dxdy*y0)
        normal_x0y0=lambda y: orthog_dxdy*y + (x0-orthog_dxdy*y0)
        rev_normal_x0y0=lambda x: (x+orthog_dxdy*y0-x0)/orthog_dxdy
        #new_row_xs=np.arange(-int(nx/2),int(nx/2),1).T*np.abs(orthog_dxdy)/np.sqrt(1+orthog_dxdy*orthog_dxdy)+x0
        new_row_xs = np.arange(-int(r_filament_pixel_display), int(r_filament_pixel_display), 1).T * np.abs(orthog_dxdy) / np.sqrt(
            1 + orthog_dxdy * orthog_dxdy) + x0
        new_row_ys=rev_normal_x0y0(new_row_xs)
        y0=y0+np.sqrt((1-dxdy*dxdy))
        #st.info(dxdy)
        x0=splev(y0,tck)

    # interpolate resampled pixles
    x_coord=np.arange(0,nx,1)
    y_coord=np.arange(0,ny,1)
    interpol=RegularGridInterpolator((x_coord,y_coord),np.transpose(data),bounds_error=False,fill_value=0)

    nx = 2*int(r_filament_pixel_display)

    new_im=np.zeros((ny,nx))
    y_init=0
    x_init=splev(y_init,tck)
    curr_y=y_init
    curr_x=x_init

    for row in range(0,ny):
        dxdy=splev(curr_y,tck,der=1)
        orthog_dxdy=-(1.0/dxdy)
        #tangent_x0y0=lambda y: dxdy*y + (curr_x-dxdy*curr_y)
        #normal_x0y0=lambda y: orthog_dxdy*y + (curr_x-orthog_dxdy*curr_y)
        
        # TODO: check 2536: RuntimeWarning: invalid value encountered in multiply (also 1671,1672 lines). Potentially related to the nans in the output
        rev_normal_x0y0=lambda x: (x+orthog_dxdy*curr_y-curr_x)/orthog_dxdy
        #new_row_xs=np.arange(-int(nx/2),int(nx/2),1).T*np.abs(orthog_dxdy)/np.sqrt(1+orthog_dxdy*orthog_dxdy)+curr_x
        new_row_xs = np.arange(-int(r_filament_pixel_display), int(r_filament_pixel_display), 1).T * np.abs(orthog_dxdy) / np.sqrt(
            1 + orthog_dxdy * orthog_dxdy) + curr_x
        new_row_ys=rev_normal_x0y0(new_row_xs)
        new_row_coords=np.vstack([new_row_xs,new_row_ys]).T
        new_row=interpol(new_row_coords)
        new_im[row,:]=new_row

        curr_y=curr_y+np.sqrt((1-dxdy*dxdy))
        curr_x=splev(curr_y,tck)

    ## fill the zeros on the edge again
    #data_slice_mean=np.median(data)
    #for i in range(ny):
    #    for j in range(nx):
    #        if new_im[i,j]==0:
    #            new_im[i,j]=data_slice_mean
    #        else:
    #            break
    #    for j in range(nx):
    #        if new_im[i,-(j+1)]==0:
    #            new_im[i,-(j+1)]=data_slice_mean
    #        else:
    #            break
    
    # with _disp_col:
    #    st.write("Straightened image:")
    #    fig = create_image_figure(new_im, apix, apix, plot_width=nx, plot_height=ny, x_axis_label=None, y_axis_label=None, tooltips=None, show_axis=False, show_toolbar=False, crosshair_color="white")
    #    st.bokeh_chart(fig, use_container_width=True)

    return new_im