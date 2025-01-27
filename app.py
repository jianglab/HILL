import numpy as np
import pandas as pd
import plotly.express as px
from shiny import App, Inputs, Outputs, Session, render, ui
from shiny.ui import div, HTML

import plotly.express as px
from shinywidgets import render_plotly

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import ImgData
from shiny import session

from shinywidgets import output_widget, render_widget, render_plotly
import pandas as pd
from shiny import reactive, req
import compute
import image_trace
import shiny_test
import plotly.graph_objects as go

# may be unnecessary
from skimage.transform import radon # resize_local_mean
import mrcfile
from shiny.types import FileInfo
import pathlib
import matplotlib.pyplot as plt

# from streamlit imports
import argparse, base64, gc, io, os, pathlib, random, socket, stat, tempfile, urllib, warnings
from getpass import getuser
from itertools import product
from math import fmod
from os import getpid
from urllib.parse import parse_qs

# TODO: delete bokeh graphs and reimplement
from bokeh.events import MouseMove, MouseEnter, DoubleTap
from bokeh.io import export_png
from bokeh.layouts import gridplot, column, layout
from bokeh.models import Button, ColumnDataSource, CustomJS, Label, LinearColorMapper, Slider, Span, Spinner
from bokeh.models.tools import CrosshairTool, HoverTool
from bokeh.plotting import figure

from finufft import nufft2d2

from numba import jit, set_num_threads, prange

import pandas as pd
from PIL import Image
from psutil import virtual_memory, Process

import qrcode

import scipy.fft
import scipy.fftpack as fp
from scipy.spatial.transform import Rotation as R
from scipy.ndimage import affine_transform, map_coordinates
from scipy.signal import correlate
from scipy.interpolate import splrep, splev
from scipy.interpolate import RegularGridInterpolator
from scipy.special import jnp_zeros
from scipy.optimize import minimize, fmin

from skimage.io import imread
from skimage import transform



# TODO: add helicon shiny here (install -r)

def compute_layer_line_positions(twist, rise, csym, radius, tilt, cutoff_res, m_max):
    # TODO Complete this function
    return {i: f"Group {i}" for i in range(m_max + 1)}


def import_with_auto_install(packages, scope=locals()):
    import importlib, site, subprocess, sys
    if isinstance(packages, str):
        packages = [packages]
    for package in packages:
        if ":" in package:
            package_import_name, package_pip_name = package.split(":")
        else:
            package_import_name, package_pip_name = package, package
        try:
            scope[package_import_name] = importlib.import_module(package_import_name)
        except ImportError:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', package_pip_name])
            importlib.reload(site)
            scope[package_import_name] = importlib.import_module(package_import_name)

#ny, nx = 200, 200 # Get the dimensions of the loaded image
nx = reactive.value(0)  # default value
ny = reactive.value(0)  # default value
nz = reactive.value(1)    # default value
MAX_SIZE = 500 # temp value


urls = {
    "empiar-10940_job010": (
        "https://tinyurl.com/y5tq9fqa",
    )
}
url_key = "empiar-10940_job010"

# may not be needed
data_all = reactive.value(None)
apix_reactive = reactive.value(2)
image_size = reactive.value(0)
displayed_class_images = reactive.value([])
displayed_class_labels = reactive.value([])
initial_selected_image_indices = reactive.value([])
selected_images = reactive.value([])
selected_image_labels = reactive.value([])

first_point = reactive.Value(None)
second_point = reactive.Value(None)

apix_auto = reactive.value(0)
negate_auto = reactive.value(0)
angle_auto = reactive.value(0)
dx_auto = reactive.value(0)
dy_auto = reactive.value(0)
mask_radius_auto = reactive.value(0)
mask_len_percent_auto = reactive.value(0)
radius_auto = reactive.value(0.0)

apix_value = reactive.value(0.0)
angle_value = reactive.value(0.0)
dx_value = reactive.value(0.0)
dy_value = reactive.value(0.0)
mask_radius_value = reactive.value(0.0)
mask_len_value = reactive.value(0.0)

value = reactive.value(200)
helical_radius = reactive.value(1)
m_max_auto_reactive = reactive.value(3)
twist = reactive.value(0.0)
min_pitch = reactive.value(0.0)
pitch = reactive.value(0.0)
rise = reactive.value(0.0)
min_rise = reactive.value(0.0)
max_rise = reactive.value(0.0)

is_3d_reactive = reactive.value(None)
data = reactive.value(None)


# my code
app_ui = ui.page_fluid(
    ui.tags.style("""
        
        .scrollable-container {
            display: flex;
            flex-direction: row;
            overflow-x: auto;
            white-space: nowrap;  /* Prevent wrapping */
        }
        .scrollable-container > .shiny-column {
            min-width: 250px;  /* Minimum width for each column */
            flex-shrink: 0;    /* Prevent columns from shrinking */
        }
        .button-container .btn {
            font-size: 1em;  /* Adjust font size */
            padding: 5px 10px; /* Adjust padding */
        }
        .wrap-text {
            word-wrap: break-word;  /* Allow text to wrap within the column */
            overflow-wrap: break-word; /* For better compatibility */
            white-space: normal;  /* Ensure normal white space handling */
        }        
        
        
    """),
    ui.layout_sidebar(
        ui.sidebar(
            ui.accordion(
                ui.accordion_panel(
                    "README",
                    ui.p("This Web app considers a biological helical structure as the product of a continous helix and a set of parallel planes. Based on the covolution theory, "
                         "the Fourier Transform (FT) of a helical structure would be the convolution of the FT of the continous helix and the FT of the planes.  \nThe FT of a continous helix consists "
                         "of equally spaced layer planes (3D) or layerlines (2D projection) that can be described by Bessel functions of increasing orders (0, ±1, ±2, ...) from the Fourier origin "
                         "(i.e. equator). The spacing between the layer planes/lines is determined by the helical pitch (i.e. the shift along the helical axis for a 360° turn of the helix). If the structure"
                          " has additional cyclic symmetry (for example, C6) around the helical axis, only the layer plane/line orders of integer multiplier of the symmetry (e.g. 0, ±6, ±12, ...) are visible. "
                          "The primary peaks of the layer lines in the power spectra form a pattern similar to a X symbol.  \nThe FT of the parallel planes consists of equally spaced points along the helical "
                          "axis (i.e. meridian) with the spacing being determined by the helical rise.  \nThe convolution of these two components (X-shaped pattern of layer lines and points along the meridian) "
                          "generates the layer line patterns seen in the power spectra of the projection images of helical structures. The helical indexing task is thus to identify the helical rise, pitch"
                           " (or twist), and cyclic symmetry that would predict a layer line pattern to explain the observed layer lines in the power spectra. This Web app allows you to interactively change "
                           "the helical parameters and superimpose the predicted layer liines on the power spectra to complete the helical indexing task.  \n  \nPS: power spectra; PD: phase differences across "
                           "the meridian; YP: Y-axis power spectra profile; LL: layer lines; m: indices of the X-patterns along the meridian; Jn: Bessel order"),
                    value="readme_panel"
                ),
                id="sidebar_accordion",
                open=False
            ),
            ui.accordion(
                ui.accordion_panel(
                    "Input Mode",
                    ui.input_radio_buttons(
                        "input_mode_params",
                        "How to obtain the input image/map:",
                        {"1": "upload", "2": "url", "3": "emd-xxxxx"},
                    ),
                    value="input_mode"
                )
            ),
            ui.output_ui("sidebar_text"),
            class_="custom-sidebar",
        ),
        ui.row(
            ui.column(12,
                ui.div(
                    ui.h2("HILL: Helical Indexing using Layer Lines"),
                    style="text-align: center; margin-bottom: 20px;"
                )
            )
        ),
        ui.div(
            ui.div(
                ui.column(2,
                    ui.div(
                        ui.output_ui("col_one"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 100%;",
                        class_="wrap-text"
                    ),
                    class_="button-container"
                ),
                ui.column(1,
                    ui.div(
                        ui.output_ui("col_two"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 100%;",
                        class_="wrap-text"
                    ),
                    class_="button-container"
                ),
                ui.column(3,
                    ui.div(
                        ui.output_ui("col_three"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 15%;",
                        class_="wrap-text"
                    ),
                    class_="button-container"
                ),
                ui.column(3,
                    ui.div(
                        ui.output_ui("col_four"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 15%;",
                        class_="wrap-text"
                    ),
                    class_="button-container"
                ),
                ui.column(3,
                    ui.div(
                        ui.output_ui("col_five"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 15%;",
                        class_="wrap-text"
                    ),
                    class_="button-container"
                ),
                class_="scrollable-container"
            ),
        ),
        ui.output_ui("dynamic_plot"),
        ui.row(
            ui.column(12,
                ui.div(
                    ui.markdown("*Developed by the [Jiang Lab@Purdue University](https://jiang.bio.purdue.edu). Report problems to [HILL@GitHub](https://github.com/jianglab/hill/issues)*"),
                )
            )
        ),
        shiny_test.setup_adjustable_sidebar(width="10vw"),
    ),
)




def server(input, output, session):

    def update_dimensions(image):
        if image is not None:
            if len(image.shape) == 2:
                nx.set(image.shape[0])
                ny.set(image.shape[1])
                nz.set(1)
                #ny, nx = image.shape
                #nz = 1
            elif len(image.shape) == 3:
                nx.set(image.shape[0])
                ny.set(image.shape[1])
                nz.set(image.shape[2])
                #nz, ny, nx = image.shape

    @output
    @render.text
    def m_max_auto():
        rise = input.rise()
        cutoff_res_y = input.cutoff_res_y()
        return int(np.floor(np.abs(rise/cutoff_res_y))) + 3

    @reactive.Calc
    def m_groups():
        return compute_layer_line_positions(
            twist=input.twist(),
            rise=input.rise(),
            csym=input.csym(),
            radius=input.helical_radius(),
            tilt=input.tilt(),
            cutoff_res=input.cutoff_res_y(),
            m_max=input.m_max()
        )

    @output
    @render.ui
    def show_choices():
        m_max = input.m_max()
        checkboxes = [
            ui.input_checkbox(
                f"m_{i}",
                label=str(i),
                value=i in [0, 1],
            )
            for i in range(m_max, -1, -1)
        ]
        return ui.div(*checkboxes)

    @output
    @render.ui
    def col_one():
        return ui.TagList(
                ui.input_action_button("col1_saveTwist", "Save twist/rise↶", class_="save-btn-success"),
                ui.input_numeric('csym', 'csym', value=6, min=1, step=1),
                ui.input_numeric('filament_diameter', 'Filament/tube diameter (Å)', value=value(), min=1.0, max=1000.0, step=10.),
                ui.input_numeric("out_of_plane_tilt", "Out-of-plane tilt (°)", value=0.0, min=-90.0, max=90.0, step=1.0),
                ui.input_numeric('res_limit_x', 'Resolution limit - X (Å)', value=3*apix_reactive(), min=2*apix_reactive(), step=1),
                ui.input_numeric("res_limit_y", "Resolution limit - Y (Å)", value=2*apix_reactive(), min=2*apix_reactive(), step=1.0),
                ui.accordion(
                    ui.accordion_panel(
                        "Additional settings",
                        ui.input_checkbox("fft_top_only", "Only display the top half of FFT", value=False),
                        ui.input_checkbox("log_amp", "Log(amplitude)", value=True),
                        ui.input_text("const_image_color", "Flatten the PS/PD image in this color", value="", placeholder="white"),
                        ui.input_text("ll_colors", 'Layerline colors', value="lime cyan violet salmon silver"),
                        ui.input_numeric("hp_fraction", 'Fourier high-pass (%)', value=0.4 * 100, min=0.0, max=100.0, step=0.1),
                        ui.input_numeric("lp_fraction", 'Fourier low-pass (%)', value=0.0 * 100, min=0.0, max=100.0, step=10.0),
                        ui.input_numeric("pnx", 'FFT X-dim size (pixels)', value=512, min=min(nx(), 128), step=2),
                        ui.input_numeric("pny", 'FFT Y-dim size (pixels)', value=1024, min=min(ny(), 512), step=2),
                        value="add_settings"
                    ),
                ),
                ui.accordion(
                    ui.accordion_panel(
                        "Simulation",
                        ui.input_numeric("ball_radius", 'Gaussian radius (Å)', value=0.0, min=0.0, max=helical_radius(), step=5.0),
                        value="simulation"
                    ),
                ),
                ui.input_checkbox("share_url", "Show/Reload sharable URL", value=False),
            )

    @output
    @render.ui
    def col_two():
        return ui.TagList(
            ui.input_action_button("col2_saveTwist", "Save twist/rise◀", class_="save-btn-success"),
            ui.h5("Display:"),
            ui.input_checkbox("PS", "PS", value=True),
            ui.input_checkbox("YP", "YP", value=True),
            ui.input_checkbox("Phase", "Phase", value=False),
            ui.input_checkbox("PD", "PD", value=True),
            ui.input_checkbox("Color", "Color", value=False),
            ui.input_checkbox("LL", "LL", value=True),
            ui.input_checkbox("LLText", "LLText", value=True),
            ui.h5("m:"),
            ui.input_numeric("m_max", "Max=", value=m_max_auto_reactive(), min=1, step=1),
            ui.output_ui("show_choices")
        )

    @output
    @render.ui
    def col_three():
        return ui.TagList(
            ui.input_numeric("spinner_twist", "Twist (°)", value=twist(), min=-180.0, max=180.0, step=1.0, width="300px"),
            ui.input_slider("slider_twist", "Twist (°)", min=-180.0, max=180.0, value=twist(), step=0.01, width="300px"),
        )
        
    @output
    @render.ui
    def col_four():
        min_pitch_v = abs(rise())
        min_pitch.set(min_pitch_v)
        return ui.TagList(
            ui.input_numeric("spinner_pitch", "Pitch (Å)", value=max(min_pitch(), pitch()), min=min_pitch(), step=1.0, width="300px"),
            ui.input_slider("slider_pitch", "Pitch (Å)", min=pitch()/2.0, max=pitch()*2.0, value=pitch(), step=pitch()*0.002, width="300px"),
        )
        
    @output
    @render.ui
    def col_five():
        return ui.TagList(
            ui.input_numeric("spinner_rise", "Rise (Å)", value=rise(), min=min_rise(), max=max_rise(), step=1.0, width="300px"),
            ui.input_slider("slider_rise", "Rise (Å)", min=rise()/2.0, max=min(max_rise(), rise()*2.0), value=rise(), step=min(max_rise(), rise()*2.0)*0.001, width="300px")
        )


    # code for the sidebar
    # when adding the code for sharable url, you can get rid of the url variables since it's only for uploading images and parsing

    @output
    @render.ui
    def conditional_3D(): 
        is_3d = input.is_3d()

        if input.is_3d():
            return ui.TagList( 
                ui.accordion(
                        ui.accordion_panel(
                            ui.p("Generate 2-D projection from the 3-D map"),
                            ui.input_checkbox("apply_helical_sym", "Apply helical symmetry", value=0),
                            ui.output_ui("helical_sym_output"),
                            ui.input_numeric("az", "Rotation around the helical axis (°):", min=0.0, max=360., value=0.0, step=1.0),
                            ui.input_numeric("tilt", "Tilt (°):", min=-180.0, max=180., value=0.0, step=1.0),
                            ui.input_numeric("noise", "Add noise (σ):", min=0.0, value=0.0, step=0.5),
                            value="2D_projection"
                        )
                ),
                ui.accordion(
                        ui.accordion_panel(
                            ui.p("Image Parameters"),
                            ui.input_radio_buttons("input_type", "Input is:", choices=["image", "PS", "PD"], inline=True),
                            ui.output_ui("image_para_values"),
                            value="image_parameters_1"
                        )
                    )
                )
        else:
            return ui.TagList(
                ui.accordion(
                    ui.accordion_panel(
                        ui.div(
                            ui.p("Choose an image"),
                            ui.output_ui("dynamic_image_select_url"),
                            style="max-height: 500px; overflow-y: auto;"
                        ),
                        value="image_selection"
                    ),
                ),
                ui.accordion(
                    ui.accordion_panel(
                        ui.p("Image Parameters"),
                        ui.input_action_button("skip_image", "Skip this image"),
                        ui.input_radio_buttons("input_type", "Input is:", choices=["image", "PS", "PD"], inline=True),
                        ui.output_ui("image_para_values"),
                        value="image_parameters_2"
                    )
                ),
            )
        

    # updating apix value
    @reactive.Effect
    @reactive.event(input.input_type, input.apix)
    def set_apix():
        input_type = input.input_type()
        apix_value = input.apix()
        
        if input_type in ["PS", "PD"]:
            apix_reactive.set(apix_value * 0.5)
            return apix_value * 0.5
        
        apix_reactive.set(apix_value)
        return apix_value

    @output
    @render.ui
    @reactive.event(input.is_3d, angle_auto, apix_reactive, dx_auto, nx, ny, mask_radius_auto, mask_len_percent_auto)
    def image_para_values():
        # TODO: ISSUE HAS SOMETHING TO DO WITH THIS FUNCTION
        # print("1: ", mask_radius_auto()*apix_reactive())
        # print("2: ", nx()/2*apix_reactive())
        # print("nx, ny: ", nx(), ny())
        # print("apix: ", apix_reactive())
        # print("mask radius: ", mask_radius_auto())


        value = input.input_type()

        if value in ['image']:
            return ui.TagList(
                            ui.input_numeric("apix", "Pixel size (Å/pixel)", value=apix_auto(), min=0.1, max=30.0, step=0.01),
                            ui.output_ui("add_transpose"),
                            ui.input_checkbox("negate", "Invert the image contrast", value=negate_auto()),
                            ui.input_checkbox("straightening", "Straighten the filament", value=False),
                            ui.input_numeric("angle", "Rotate (°)", value=-angle_auto(), min=-180.0, max=180.0, step=1.0),
                            ui.input_numeric("dx", "Shift along X-dim (Å)", value=dx_auto()*apix_reactive(), min=-nx()*apix_reactive(), max=nx()*apix_reactive(), step=1.0),
                            ui.input_numeric("dy", "Shift along Y-dim (Å)", value=0.0, min=-ny()*apix_reactive(), max=ny()*apix_reactive(), step=1.0),
                            ui.input_numeric("mask_radius", "Mask radius (Å)", value=min(mask_radius_auto()*apix_reactive(), nx()/2*apix_reactive()), min=1.0, max=nx()/2*apix_reactive(), step=1.0),
                            ui.input_numeric("mask_len", "Mask length (%)", value=mask_len_percent_auto(), min=10.0, max=100.0, step=1.0),
            )
        elif value in ['PS', 'PD']:
            return ui.TagList(
                            ui.input_numeric("apix", "Nyquist res (Å)", value=2*apix_auto(), min=0.1, max=30.0, step=0.01),
                            ui.output_ui("add_transpose"),
                            ui.input_checkbox("negate", "Invert the image contrast", value=negate_auto()),
                            ui.input_numeric("angle", "Rotate (°)", value=-angle_auto(), min=-180.0, max=180.0, step=1.0),
                            ui.input_numeric("dx", "Shift along X-dim (Å)", value=dx_auto()*apix_reactive(), min=-nx()*apix_reactive(), max=nx()*apix_reactive(), step=1.0),
                            ui.input_numeric("dy", "Shift along Y-dim (Å)", value=0.0, min=-ny()*apix_reactive(), max=ny()*apix_reactive(), step=1.0),
            )

    @output
    @render.ui
    @reactive.event(input.is_3d)
    def add_transpose():
        ####### CONTINUE FROM HERE!!
        transpose_auto = input.input_mode_params() not in [2, 3] and nx() > ny()
        print("is_3d: ", input.is_3d())
        if input.is_3d():
            return ui.TagList(
                    ui.input_checkbox("transpose", "Transpose the image", value=transpose_auto),
            )

    # updating apix_auto value
    @reactive.Effect
    @reactive.event(input.input_mode_params, input.url_params)
    def set_apix_auto_url():
        mode = input.input_mode_params()
        # {"1": "upload", "2": "url", "3": "emd-xxxxx"}

        if mode == "1": # upload
            pass

            fileobj = input.upload_classes()
            data_all_v, map_crs_auto_v, apix_auto_v = get_2d_image_from_uploaded_file(fileobj)
            apix_auto.set(apix_auto_v)

            # is_3d_auto = guess_if_3d(filename=fileobj.name, data=data_all())
            # is_3d_reactive.set(is_3d_auto)
        elif mode == "2": # url
            image_url = input.url_params()
            data_all_v, map_crs_auto_v, apix_auto_v = get_2d_image_from_url(image_url)
            apix_auto.set(apix_auto_v)
        else: # emd
            # data_all, map_crs_auto, apix_auto = get_emdb_map(emd_id)
            pass
        
    @reactive.Effect
    @reactive.event(input.input_mode_params, input.upload_classes)
    def set_apix_auto_upload():
        mode = input.input_mode_params()
        if mode == "1" and input.upload_classes(): # upload
            fileobj = input.upload_classes()
            file_info = fileobj[0]['datapath']
            data_all_v, map_crs_auto_v, apix_auto_v = get_2d_image_from_uploaded_file(input.upload_classes()[0])
            apix_auto.set(apix_auto_v)

    # TODO: make set_apix_auto for emd selection

    @reactive.Effect
    @reactive.event(input.input_type, data, selected_images)
    def set_angle_auto():
        mode = input.input_type() 

        if mode in ["PS", "PD"]:
            angle_auto.set(0.0)
            dx_auto.set(0.0)
        elif is_3d_reactive():
            angle_auto.set(0.0)
            dx_auto.set(0.0)
        elif selected_images() and len(selected_images()) > 0 and data() is not None:
            ang_auto_v, dx_auto_v = auto_vertical_center(selected_images()[0])# data())
            angle_auto.set(ang_auto_v)
            dx_auto.set(dx_auto_v)

        if selected_images() and len(selected_images()) > 0:
            aspect_ratio = float(nx() / ny())
        
        if input.straightening() and aspect_ratio < 1 and selected_images() and len(selected_images()) > 0:
            aspect_ratio = float(nx() / ny())
            angle_auto.set(0.0)

        # print("angle auto: ", angle_auto())


    @reactive.Effect
    @reactive.event(input.input_type, selected_images, data)
    def set_mask_radius_auto(): 
        # radius_auto, mask_radius_auto = estimate_radial_range(data, thresh_ratio=0.1)
        input_type = input.input_type()
        if input_type in ["image"] and selected_images() and len(selected_images()) > 0 and data() is not None:
            radius_auto_v, mask_radius_auto_v = estimate_radial_range(data(), thresh_ratio=0.1)
            # print("mask radius", mask_radius_auto_v)
            mask_radius_auto.set(mask_radius_auto_v)
            radius_auto.set(radius_auto_v)

        mask_len_percent_auto_v = 90.0
        if input.straightening():
            mask_len_percent_auto_v = 100.0
        mask_len_percent_auto.set(mask_len_percent_auto_v)

    
    @reactive.Effect
    @reactive.event(input.angle, input.dx, input.dy, input.apix, selected_images)
    def set_data():
        # Initialize data_v
        data_v = None
        
        # Check if there are selected images
        if selected_images() and len(selected_images()) > 0:
            data_v = selected_images()[0]
            # print("Initial data_v retrieved.")

        # Apply transformations conditionally
        if input.is_3d():
            # print("Processing as 3D input.")
            if input.transpose():
                # print("Applying transpose.")
                data_v = data_v.T
            
        if input.negate():
            # print("Inverting image contrast.")
            data_v = -data_v

        # Handle rotation and shift
        angle_value = input.angle()
        dx_value = input.dx()
        dy_value = input.dy()
        apix_value = input.apix()
        
        # print(f"Transform parameters - Angle: {angle_value}, Dx: {dx_value}, Dy: {dy_value}, Apix: {apix_value}")

        if (angle_value or dx_value or dy_value) and data_v is not None:
            # print("Applying rotation and shift.")
            data_v = rotate_shift_image(
                data_v, 
                angle=-angle_value, 
                post_shift=(dy_value / apix_value, dx_value / apix_value), 
                order=1
            )
            # print("Rotation and shift applied.")

        # Update the reactive `data` object
        if data_v is not None:
            data.set(data_v)
            # print("Updated data reactive with transformed data_v.")
   
    # code added from the Helical Pitch Image Selection:
    @output
    @render.ui
    def sidebar_text():
        selection = input.input_mode_params()
        
        # Create the is_3d checkbox with current dimensions
        is_3d_checkbox = ui.input_checkbox(
            "is_3d",
            f"The input ({nx()}x{ny()}x{nz()}) is a 3D map",
            value=(nz() > 1)  # Automatically set based on nz dimension
        )
        
        if selection == "1":  # Upload
            return ui.TagList(
                ui.p("Upload a mrc or mrcs file"),
                ui.input_file(
                    "upload_classes",
                    "Upload the class averages in MRC format (.mrcs, .mrc)",
                    accept=[".mrcs", ".mrc"],
                    placeholder="mrcs or mrc file",
                ),
                ui.input_action_button("run", label="Run", style="width: 100%;"),
                ui.input_checkbox("is_3d", "The input is a 3D map", value=False),
                # ui.input_checkbox("ignore_blank", "Ignore blank classes", value=True),
                ui.output_ui("conditional_3D"),
                output_widget("display_micrograph"),
                output_widget("transformed_display_micrograph"),
                output_widget("plot_graph"),
                output_widget("acf_plot"),
            )
        elif selection == "2":  # URL
            return ui.TagList(
                ui.input_text(
                    "url_params",
                    "Input a url of 2D image(s) or a 3D map:",
                    value=urls[url_key][0],
                ),
                ui.input_action_button("run", label="Run", style="width: 100%;"),
                ui.input_checkbox("is_3d", "The input is a 3D map", value=False),
                # ui.input_checkbox("ignore_blank", "Ignore blank classes", value=True),
                ui.output_ui("conditional_3D"),
                output_widget("display_micrograph"),
                output_widget("transformed_display_micrograph"),
                output_widget("plot_graph"),
                output_widget("acf_plot"),

            )
        elif selection == "3":  # EMD-xxxxx
            return ui.TagList(
                    ui.p("You have selected to use an EMD file. Please enter the EMD accession number (e.g., EMD-xxxxx)."),
                    ui.input_text("input_emd", "Input an EMDB ID (emd-xxxxx):", value="emd-10499"),
                    ui.p("EMD-10499 | resolution=3.9Å\ntwwist=-31.44° | rise=6.721Å | c2"), # temp values
                    ui.accordion(
                        ui.accordion_panel(
                            ui.p("Generate 2-D projection from the 3-D map"),
                            ui.input_checkbox("apply_helical_sym", "Apply helical symmetry", value=0),
                            ui.input_numeric("az", "Rotation around the helical axis (°):", min=0.0, max=360., value=0.0, step=1.0),
                            ui.input_numeric("tilt", "Tilt (°):", min=-180.0, max=180., value=0.0, step=1.0),
                            ui.input_numeric("noise", "Add noise (σ):", min=0.0, value=0.0, step=0.5),
                            value="2D_projection_emd"
                        )
                    ), 
                    ui.accordion(
                            ui.accordion_panel(
                                ui.p("Image Parameters"),
                                ui.input_radio_buttons("input_type", "Input is:", choices=["image", "PS", "PD"], inline=True),
                                ui.input_numeric("apix", "Pixel size (Å/pixel)", value=apix_auto(), min=0.1, max=30.0, step=0.01),
                                ui.input_checkbox("transpose", "Transpose the image", value=negate_auto()),
                                ui.input_checkbox("negate", "Invert the image contrast", value=negate_auto()),
                                ui.input_numeric("angle", "Rotate (°)", value=-angle_auto(), min=-180.0, max=180.0, step=1.0),
                                ui.input_numeric("dx", "Shift along X-dim (Å)", value=dx_auto()*apix_reactive(), min=-nx()*apix_reactive(), max=nx()*apix_reactive(), step=1.0),
                                ui.input_numeric("dy", "Shift along Y-dim (Å)", value=0.0, min=-ny()*apix_reactive(), max=ny()*apix_reactive(), step=1.0),
                                ui.input_numeric("mask_radius", "Mask radius (Å)", value=min(mask_radius_auto()*apix_reactive(), nx()/2*apix_reactive()), min=1.0, max=nx()/2*apix_reactive(), step=1.0),
                                ui.input_numeric("mask_len", "ask length (%)", value=mask_len_percent_auto(), min=10.0, max=100.0, step=1.0),
                                value="image_parameters_emd"
                            )
                        )             
                )
            # try to do random input, otherwise it's ok to leave out
        else:
            return ui.p("Please select an option to proceed.")

    @output
    @render.ui
    @reactive.event(input.apply_helical_sym)
    def helical_sym_output():
        value = input.apply_helical_sym()
        print("here")
        if value:
            print("in val")
            return ui.TagList(
                    # HOW TO GET THE FOLLOWING SESSION STATE VALUES:
                    # ui.input_numeric("twist_ahs", "Twist (°):", value=???, min=-180.0, max=180.0, step=1.0),
                    # ui.input_numeric("rise_ahs", "Rise (Å):", value=???, min=0.0, step=1.0),
                    # ui.input_numeric("csym_ahs", "Csym:", value=???, min=1, step=1),
                    ui.input_numeric("apix_map", "Current map pixel size (Å):", value=apix_auto(), min=0.0, step=1.0),
                    ui.output_ui("update_apix_map_vals"),
                    ui.hr(),
            )

    @output
    @render.ui
    @reactive.event(input.apix_map)
    def update_apix_map_vals():
        return ui.TagList(
            ui.input_numeric("apix_ahs", "New map pixel size (Å):", value=input.apix_map(), min=0.0, step=1.0),
            # ui.input_numeric(
            #     "fraction_ahs",
            #     "Center fraction (0-1):",
            #     value=1.0,
            #     min=input.rise_ahs()/(nz()*input.apix_map()),
            #     max=1.0,
            #     step=0.1,
            # ),
            # ui.input_numeric(
            #     "length_ahs",
            #     "Box length (Å):",
            #     value=input.apix_map()*max(nz(),nx()),
            #     min=input.rise_ahs(),
            #     step=1.0,
            # ),
            ui.input_numeric("width_ahs", "Box width (Å):", value=input.apix_map()*nx(), min=0.0, step=1.0),
        )


    @reactive.effect
    @reactive.event(data_all) # , input.ignore_blank)
    def get_displayed_class_images():
        if data_all() is not None:
            data = data_all()
            n = len(data)
            images = [data[i] for i in range(n)]
            image_size.set(max(images[0].shape))

            update_dimensions(images[0])

            # if input.ignore_blank():
            #     included = []
            #     included_images = []
            #     for i in range(n):
            #         image = images[i]
            #         if np.max(image) > np.min(image):
            #             included.append(i)
            #             included_images.append(image)
            #     images = included_images
            # else:
            #     included = list(range(n))

            included = []
            included_images = []
            for i in range(n):
                image = images[i]
                if np.max(image) > np.min(image):
                    included.append(i)
                    included_images.append(image)
            images = included_images
            
            image_labels = [f"Class {i+1}" for i in included]
            displayed_class_labels.set(image_labels)
            displayed_class_images.set(images)
            
            # Reset selected images when display changes
            selected_images.set([])
            selected_image_labels.set([])
    
    @output
    @render.ui
    def display_selected_images():
        imgs = selected_images()
        labels = selected_image_labels()
        if len(imgs) > 0 and len(labels) > 0:
            return shiny_test.image_gallery(
                id="display_selected_image",
                label="Selected classe(s):",
                images=imgs,
                image_labels=labels,
                enable_selection=False,
                allow_multiple_selection=False,
            )
        return ui.div("No images selected")
    
    
    @output
    @render.ui
    def dynamic_image_select_upload():
        if len(displayed_class_images()) > 0:
            return shiny_test.image_gallery(
                id="select_classes",
                label="Select class(es):",
                images=displayed_class_images(),
                image_labels=displayed_class_labels(),
                image_size=128,
                initial_selected_indices=initial_selected_image_indices(),
                enable_selection=True,
                allow_multiple_selection=False,
            )
        return ui.div("No images available")
    
    @output
    @render.ui
    def dynamic_image_select_url():
        if len(displayed_class_images()) > 0:
            return shiny_test.image_gallery(
                id="select_classes",
                label="Select class(es):",
                images=displayed_class_images(),
                image_labels=displayed_class_labels(),
                image_size=128,
                initial_selected_indices=initial_selected_image_indices(),
                enable_selection=True,
                allow_multiple_selection=False,
            )
        return ui.div("No images available")

    @reactive.effect
    @reactive.event(input.select_classes)
    def update_selected_images():
        if input.select_classes() is not None:
            selected = [displayed_class_images()[i] for i in input.select_classes()]
            selected_images.set(
                [displayed_class_images()[i] for i in input.select_classes()]
            )
            selected_image_labels.set(
                [displayed_class_labels()[i] for i in input.select_classes()]
            )
            if selected:
                update_dimensions(selected[0])
            

    @reactive.Effect
    @reactive.event(input.run)
    def get_class2d():
        if input.input_mode_params() == "1":
            fileinfo = input.upload_classes()
            if fileinfo is not None:
                class_file = fileinfo[0]["datapath"]
                try:
                    data, apix = compute.get_class2d_from_file(class_file)
                    #apix_reactive.set(apix)
                    update_dimensions(data[0])
                    #nx = data.shape[-1]
                    image_size.set(nx)
                    data_all.set(data)
                    image_size.set(nx)
                except Exception as e:
                    print(e)
                    modal = ui.modal(
                        f"Failed to read the uploaded 2D class average images from {fileinfo[0]['name']}",
                        title="File upload error",
                        easy_close=True,
                        footer=None
                    )
                    ui.modal_show(modal)
        
        elif input.input_mode_params() == "2" and input.url_params():
            url = input.url_params()
            try:
                data, apix = compute.get_class2d_from_url(url)
                #apix_reactive.set(apix)
                update_dimensions(data[0])
                #nx = data.shape[-1]
                image_size.set(nx())
                data_all.set(data)
                image_size.set(nx())
            except Exception as e:
                print(e)
                modal = ui.modal(
                    f"Failed to download 2D class average images from {url}",
                    title="URL download error",
                    easy_close=True,
                    footer=None
                )
                ui.modal_show(modal)

    # functions for outputting graphs
    @render_plotly
    @reactive.event(selected_images)
    def display_micrograph():
        images = selected_images()
        if not images:  # Check for None or empty list
            return px.scatter(title="No image selected")

        h, w = selected_images()[0].shape[:2]
        nx.set(h)
        ny.set(w)

        fig = image_trace.plot_micrograph(
            micrograph=255 - selected_images()[0],
            title=f"Original image ({nx()}x{ny()})",
            apix=apix_reactive(),
        )

        def plot_micrograph_on_click(trace, points, selector):
            if selector.shift == True:
                first_point.set((points.xs[0], points.ys[0]))
            else:
                first_point.set(None)
                second_point.set(None)

        def plot_micrograph_on_hover(trace, points, selector):
            if first_point() is None:
                second_point.set(None)
                return

            if selector.shift == True:
                second_point.set((points.xs[0], points.ys[0]))

        for data in fig.data:
            if data.name == "image":
                data.on_click(plot_micrograph_on_click)
                data.on_hover(plot_micrograph_on_hover)

        return fig
    

    @render_plotly
    @reactive.event(selected_images, input.negate)
    def transformed_display_micrograph():
        images = selected_images()
        if not images:
            return px.scatter(title="No image selected")

        h, w = selected_images()[0].shape[:2]
        nx.set(h)
        ny.set(w)

        image_to_display = (
            selected_images()[0] if input.negate() else 255 - selected_images()[0]
        )

        fig = image_trace.plot_micrograph(
            micrograph= image_to_display, #255 - selected_images()[0],
            title=f"Transformed image ({nx()}x{ny()})",
            apix=apix_reactive(),
        )

        def plot_micrograph_on_click(trace, points, selector):
            if selector.shift == True:
                first_point.set((points.xs[0], points.ys[0]))
            else:
                first_point.set(None)
                second_point.set(None)

        def plot_micrograph_on_hover(trace, points, selector):
            if first_point() is None:
                second_point.set(None)
                return

            if selector.shift == True:
                second_point.set((points.xs[0], points.ys[0]))

        for data in fig.data:
            if data.name == "image":
                data.on_click(plot_micrograph_on_click)
                data.on_hover(plot_micrograph_on_hover)

        return fig
    

    # @reactive.Calc
    # def mask_parameters():
    #     images = selected_images()
    #     if not images:
    #         return px.scatter(title="No image selected")
        
    #     data = selected_images()[0]

    #     straightening = input.straightening()
    #     radius_auto, mask_radius_auto_v = estimate_radial_range(data)
    #     mask_radius_auto.set(mask_radius_auto_v)
    #     mask_radius = input.mask_radius()
    #     mask_len_percent_auto_v = 100.0 if straightening else 90.0
    #     mask_len_percent_auto.set(mask_len_percent_auto_v)
    #     mask_len_fraction = input.mask_len() / 100.0
    #     return mask_radius, mask_len_fraction

    @render_plotly
    @reactive.event(selected_images)
    def plot_graph():
        images = selected_images()
        if not images:
            #print("No selected images!")
            return px.scatter(title="No image selected")
        
        data = selected_images()[0]
        mask_radius = mask_radius_auto()
        mask_len_fraction = input.mask_len() / 100.0
        
        x = np.arange(-nx() // 2, nx() // 2) * apix_reactive()
        ymax = np.max(data, axis=0)
        ymean = np.mean(data, axis=0)
        
        fig = go.Figure()
        
        # Add max line
        fig.add_trace(go.Scatter(x=x, y=ymax, mode='lines', name='max', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=-x, y=ymax, mode='lines', name='max flipped', line=dict(color='red', dash='dash')))
        
        # Add mean line
        fig.add_trace(go.Scatter(x=x, y=ymean, mode='lines', name='mean', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=-x, y=ymean, mode='lines', name='mean flipped', line=dict(color='blue', dash='dash')))
        
        # Add mask radius spans
        fig.add_shape(type="line", x0=-mask_radius, x1=-mask_radius, y0=0, y1=np.max(ymax), 
                      line=dict(color="green", dash="dash"), name='mask_min')
        fig.add_shape(type="line", x0=mask_radius, x1=mask_radius, y0=0, y1=np.max(ymax), 
                      line=dict(color="green", dash="dash"), name='mask_max')

        fig.update_layout(
            xaxis_title="x (Å)",
            yaxis_title="Pixel value",
            showlegend=True,
            legend=dict(x=0.9, y=1.1),
            height=400,
            #width=250
            autosize=True,
            margin=dict(l=10, r=10, t=50, b=10),
            plot_bgcolor="white",
            title=dict(
                text="x (Å) vs. Pixel value",
                x=0.5,
                y=0.95,
                xanchor="center",
                font=dict(size=14)
            )
        )


        return fig

    @reactive.Calc
    def acf_data():
        images = selected_images()
        if not images:
            #print("No selected images!")
            return None, None, None

        data = images[0]

        acf = auto_correlation(data, sqrt=True, high_pass_fraction=0.1)

        if acf is None or np.isnan(acf).any():
            print("Invalid auto-correlation data!")
            return None, None, None

        ny = acf.shape[0]
        y = np.arange(-ny // 2, ny // 2) * apix_reactive() 
        xmax = np.max(acf, axis=1)  

        return xmax, y, acf


    @render_plotly
    @reactive.event(selected_images)
    def acf_plot():
        xmax, y, acf = acf_data()

        if xmax is None or y is None:
            return px.scatter(title="No valid data available")

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=xmax, y=y, mode='lines', name="ACF", line=dict(color='red')
        ))

        fig.update_layout(
            xaxis_title="Auto-correlation",
            yaxis_title="Axial Shift (Å)",
            showlegend=True,
            height=400,
            #width=250,
            hovermode="x unified",
            autosize=True,
            margin=dict(l=10, r=10, t=50, b=10),
            plot_bgcolor="white",
            title=dict(
                text="Auto-correlation vs. Axial Shift (Å)",
                x=0.5,
                y=0.95,
                xanchor="center",
                font=dict(size=14)
            )
        )

        return fig



app = App(app_ui, server)


###########################

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

def compute_layer_line_positions(twist, rise, csym, radius, tilt, cutoff_res, m_max=-1):
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
        cutoff_res_y, cutoff_res_x = cutoff_res
    else:
        cutoff_res_y, cutoff_res_x = 2*apix, 2*apix
    if output_size:
        ony, onx = output_size
    else:
        ony, onx = image.shape
    freq_y = np.fft.fftfreq(ony) * 2*apix/cutoff_res_y
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

def auto_vertical_center(data, n_theta=180):
  #from skimage.transform import radon
  #from scipy.signal import correlate
  
  data_work = np.clip(data, 0, None)
  
  theta = np.linspace(start=0., stop=180., num=n_theta, endpoint=False)
  #import warnings
  with warnings.catch_warnings(): # ignore outside of circle warnings
    warnings.simplefilter('ignore')
    sinogram = radon(data_work, theta=theta)
  sinogram += sinogram[::-1, :]
  y = np.std(sinogram, axis=0)
  theta_best = -theta[np.argmax(y)]

  rotated_data = rotate_shift_image(data_work, angle=theta_best)
  # now find best vertical shift
  yproj = np.sum(rotated_data, axis=0)
  yproj_xflip = yproj*1.0
  yproj_xflip[1:] = yproj[1:][::-1]
  corr = correlate(yproj, yproj_xflip, mode='same')
  shift_best = -(np.argmax(corr) - len(corr)//2)/2

  # refine to sub-degree, sub-pixel level
  def score_rotation_shift(x):
    theta, shift_x = x
    data_tmp=rotate_shift_image(data_work, angle=theta, post_shift=(0, shift_x), order=1)
    xproj = np.sum(data_tmp, axis=0)[1:]
    xproj += xproj[::-1]
    score = -np.std(xproj)
    return score
  #from scipy.optimize import fmin
  res = fmin(score_rotation_shift, x0=(theta_best, shift_best), xtol=1e-2, disp=0)
  theta_best, shift_best = res
  return set_to_periodic_range(theta_best), shift_best

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

def generate_projection(data, az=0, tilt=0, noise=0, output_size=None):
    #from scipy.spatial.transform import Rotation as R
    #from scipy.ndimage import affine_transform
    # note the convention change
    # xyz in scipy is zyx in cryoEM maps
    if az or tilt:
        rot = R.from_euler('zx', [tilt, az], degrees=True)  # order: right to left
        m = rot.as_matrix()
        nx, ny, nz = data.shape
        bcenter = np.array((nx//2, ny//2, nz//2), dtype=np.float32)
        offset = bcenter.T - np.dot(m, bcenter.T)
        data_work = affine_transform(data, matrix=m, offset=offset, mode='nearest')
    else:
        data_work = data
    ret = data_work.sum(axis=1)   # integrate along y-axis
    if output_size is not None:
        ony, onx = output_size
        ny, nx = ret.shape
        if ony!=ny or onx!=nx:
            top_bottom_mean = np.mean(ret[(0,-1),:])
            ret2 = np.zeros((ony, onx), dtype=np.float32) + top_bottom_mean
            y0 = ony//2 - ny//2
            x0 = onx//2 - nx//2
            ret2[y0:y0+ny, x0:x0+nx] = ret
            ret = ret2 
    ret = normalize(ret)
    if noise:
        ret += np.random.normal(loc=0.0, scale=noise*np.std(data[data!=0]), size=ret.shape)
    return ret

def apply_helical_symmetry(data, apix, twist_degree, rise_angstrom, csym=1, fraction=1.0, new_size=None, new_apix=None, cpu=1):
  if new_apix is None: new_apix = apix  
  nz0, ny0, nx0 = data.shape
  if new_size != data.shape:
    nz1, ny1, nx1 = new_size
    nz2, ny2, nx2 = max(nz0, nz1), max(ny0, ny1), max(nx0, nx1)
    data_work = np.zeros((nz2, ny2, nx2), dtype=np.float32)
  else:
    data_work = np.zeros((nz0, ny0, nx0), dtype=np.float32)

  nz, ny, nx = data_work.shape
  w = np.zeros((nz, ny, nx), dtype=np.float32)

  hsym_max = max(1, int(nz*new_apix/rise_angstrom))
  hsyms = range(-hsym_max, hsym_max+1)
  csyms = range(csym)

  mask = (data!=0)*1
  z_nonzeros = np.nonzero(mask)[0]
  z0 = np.min(z_nonzeros)
  z1 = np.max(z_nonzeros)
  z0 = max(z0, nz0//2-int(nz0*fraction+0.5)//2)
  z1 = min(nz0-1, min(z1, nz0//2+int(nz0*fraction+0.5)//2))

  set_num_threads(cpu)

  for hi in hsyms:
    for k in prange(nz):
      k2 = ((k-nz//2)*new_apix + hi * rise_angstrom)/apix + nz0//2
      if k2 < z0 or k2 >= z1: continue
      k2_floor, k2_ceil = int(np.floor(k2)), int(np.ceil(k2))
      wk = k2 - k2_floor

      for ci in csyms:
        rot = np.deg2rad(twist_degree * hi + 360*ci/csym)
        m = np.array([
              [ np.cos(rot),  np.sin(rot)],
              [-np.sin(rot),  np.cos(rot)]
            ])
        for j in prange(ny):
          for i in prange(nx):
            j2 = (m[0,0]*(j-ny//2) + m[0,1]*(i-nx/2))*new_apix/apix + ny0//2
            i2 = (m[1,0]*(j-ny//2) + m[1,1]*(i-nx/2))*new_apix/apix + nx0//2

            j2_floor, j2_ceil = int(np.floor(j2)), int(np.ceil(j2))
            i2_floor, i2_ceil = int(np.floor(i2)), int(np.ceil(i2))
            if j2_floor<0 or j2_floor>=ny0-1: continue
            if i2_floor<0 or i2_floor>=nx0-1: continue

            wj = j2 - j2_floor
            wi = i2 - i2_floor

            data_work[k, j, i] += (
                (1 - wk) * (1 - wj) * (1 - wi) * data[k2_floor, j2_floor, i2_floor] +
                (1 - wk) * (1 - wj) * wi * data[k2_floor, j2_floor, i2_ceil] +
                (1 - wk) * wj * (1 - wi) * data[k2_floor, j2_ceil, i2_floor] +
                (1 - wk) * wj * wi * data[k2_floor, j2_ceil, i2_ceil] +
                wk * (1 - wj) * (1 - wi) * data[k2_ceil, j2_floor, i2_floor] +
                wk * (1 - wj) * wi * data[k2_ceil, j2_floor, i2_ceil] +
                wk * wj * (1 - wi) * data[k2_ceil, j2_ceil, i2_floor] +
                wk * wj * wi * data[k2_ceil, j2_ceil, i2_ceil]
            )
            w[k, j, i] += 1.0
  mask = w>0
  data_work = np.where(mask, data_work / w, data_work)
  if data_work.shape != new_size:
    nz1, ny1, nx1 = new_size
    data_work = data_work[nz//2-nz1//2:nz//2+nz1//2, ny//2-ny1//2:ny//2+ny1//2, nx//2-nx1//2:nx//2+nx1//2]
  return data_work

def normalize(data, percentile=(0, 100)):
    p0, p1 = percentile
    vmin, vmax = sorted(np.percentile(data, (p0, p1)))
    data2 = (data-vmin)/(vmax-vmin)
    return data2

def nonzero_images(data, thresh_ratio=1e-3):
    assert(len(data.shape) == 3)
    sigmas = np.std(data, axis=(1,2))
    thresh = sigmas.max() * thresh_ratio
    nonzeros = np.where(sigmas>thresh)[0]
    if len(nonzeros)>0: 
        return nonzeros
    else:
        None

def guess_if_is_phase_differences_across_meridian(data, err=30):
    if np.any(data[:, 0]):
        return False
    if not (data.min()==0 and (0<=180-data.max()<err)):
        return False
    sym_diff = data[:, 1:] - data[:, 1:][:, ::-1]
    if np.any(sym_diff):
        return False
    return True

def guess_if_is_power_spectra(data, thresh=15):
    median = np.median(data)
    max = np.max(data)
    sigma = np.std(data)
    if (max-median)>thresh*sigma: return True
    else: return False

def guess_if_is_positive_contrast(data):
    y_proj = np.sum(data, axis=0)
    mean_edge = np.mean(y_proj[[0,1,2,-3,-2,-1]])
    if np.max(y_proj)-mean_edge > abs(np.min(y_proj)-mean_edge): return True
    else: return False

def guess_if_3d(filename, data=None):
    if filename.endswith(".mrcs"): return False
    if filename.startswith("cryosparc") and filename.endswith("_class_averages.mrc"): return False    # cryosparc_P*_J*_*_class_averages.mrc
    if data is None: return None
    if len(data.shape)<3: return False
    if len(data.shape)>3: return None
    nz, ny, nx = data.shape
    if nz==1: return False
    if nz==ny and nz==nx: return True
    if ny==nx and nz in [50, 100, 200]: return False
    return None

def get_2d_image_from_uploaded_file(fileobj):
    #import os, tempfile
    #original_filename = fileobj.name
    # data_all = mrcfile.open(fileobj, permissive=True).data
    # map_crs_auto = None
    # apix_auto = None

    # if data_all is not None:
    #     apix_auto = 1.0  # or whatever default value you want to use
    #     map_crs_auto = data_all.copy()  # if you need a copy of the data
    # original_filename = fileobj.name
    if fileobj is None:
        return 0, 0, 0
    
    original_filename = fileobj['name']
    suffix = os.path.splitext(original_filename)[-1]
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as temp:
        with open(fileobj['datapath'], 'rb') as src:
            temp.write(src.read())  # Copy content to temporary file
        temp.flush()
    data, map_crs, apix = get_2d_image_from_file(temp.name)

    return data.astype(np.float32), map_crs, apix

    # return data_all, map_crs_auto, apix_auto

def get_emdb_ids():
    try:
        import_with_auto_install(["pandas"])
        #import pandas as pd
        entries_all = pd.read_csv('https://www.ebi.ac.uk/emdb/api/search/current_status:"REL"?rows=1000000&wt=csv&download=true&fl=emdb_id,structure_determination_method,resolution,image_reconstruction_helical_delta_z_value,image_reconstruction_helical_delta_phi_value,image_reconstruction_helical_axial_symmetry_details')
        entries_all["emdb_id"] = entries_all["emdb_id"].str.split("-", expand=True).iloc[:, 1].astype(str)
        emdb_ids_all = list(entries_all["emdb_id"])
        methods = list(entries_all["structure_determination_method"])
        resolutions = list(entries_all["resolution"])
        emdb_helical = entries_all[entries_all["structure_determination_method"]=="helical"].rename(columns={"image_reconstruction_helical_delta_z_value": "rise", "image_reconstruction_helical_delta_phi_value": "twist", "image_reconstruction_helical_axial_symmetry_details": "csym"}).reset_index()
        emdb_ids_helical = list(emdb_helical["emdb_id"])
    except:
        emdb_ids_all = []
        methods = []
        resolutions = []
        emdb_helical = None
        emdb_ids_helical = []
        #st.warning("WARNING: failed to obtain the list of EMDB entries")
    return emdb_ids_all, methods, resolutions, emdb_helical, emdb_ids_helical

def get_emdb_helical_parameters(emd_id):
    emdb_helical, emdb_ids_helical = get_emdb_ids()[-2:]
    if emdb_helical is not None and emd_id in emdb_ids_helical:
        row_index = emdb_helical.index[emdb_helical["emdb_id"] == emd_id].tolist()[0]
        row = emdb_helical.iloc[row_index]
        ret = {}
        try:
            csym = int(row.csym[1:])
        except: 
            csym = 1
            ret["csym_known"] = False
        ret.update({"resolution":row.resolution, "twist":float(row.twist), "rise":float(row.rise), "csym":csym})
    else:
        ret = {}
    return ret

def get_emdb_map_url(emd_id: str):
    #server = "https://files.wwpdb.org/pub"    # Rutgers University, USA
    server = "https://ftp.ebi.ac.uk/pub/databases" # European Bioinformatics Institute, England
    #server = "http://ftp.pdbj.org/pub" # Osaka University, Japan
    url = f"{server}/emdb/structures/EMD-{emd_id}/map/emd_{emd_id}.map.gz"
    return url

def get_emdb_map(emd_id: str):
    url = get_emdb_map_url(emd_id)
    fileobj = download_file_from_url(url)
    if fileobj is None:
        #st.error(f"ERROR: EMD-{emd_id} map ({url}) could not be downloaded")
        #st.stop()
        @reactive.Effect
        def show_error():
            ui.output_text("error_message", f"ERROR: EMD-{emd_id} map ({url}) could not be downloaded")
        raise Exception(f"EMD-{emd_id} map could not be downloaded")
    #import mrcfile
    with mrcfile.open(fileobj.name, mode='r') as mrc:
        map_crs = [int(mrc.header.mapc), int(mrc.header.mapr), int(mrc.header.maps)]
        vmin, vmax = np.min(mrc.data), np.max(mrc.data)
        data = ((mrc.data - vmin) / (vmax - vmin))
        apix = mrc.voxel_size.x.item()
    return data.astype(np.float32), map_crs, apix

def get_2d_image_from_url(url):
    url_final = get_direct_url(url)    # convert cloud drive indirect url to direct url
    fileobj = download_file_from_url(url_final)
    if fileobj is None:
        #st.error(f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview")
        #st.stop()
        @reactive.Effect
        def show_error():
            ui.output_text(
                "error_message",
                (
                    f"ERROR: {url} could not be downloaded. If this URL points to a cloud drive file, "
                    "make sure the link is a direct download link instead of a link for preview."
                ),
            )
        # Raise an exception to stop further processing
        raise Exception("File download failed")
    data = get_2d_image_from_file(fileobj.name)
    return data

def get_2d_image_from_file(filename):
    try:
        with mrcfile.open(filename, permissive=True) as mrc:
            map_crs = [int(mrc.header.mapc), int(mrc.header.mapr), int(mrc.header.maps)]
            data = mrc.data.astype(np.float32)
            apix = mrc.voxel_size.x.item()
    except Exception as e:
        print(f"Error reading MRC file: {e}")
        try:
            from skimage.io import imread
            data = imread(filename, as_gray=True) * 1.0
            data = data[::-1, :]
            apix = 1.0
            map_crs = [1, 2, 3]
        except Exception as e:
            print(f"Error reading file as image: {e}")
            raise ValueError("Unsupported file format or corrupted file.")
    return data, map_crs, apix

def download_file_from_url(url):
    import tempfile
    import requests
    try:
        filesize = get_file_size(url)
        local_filename = url.split('/')[-1]
        suffix = '.' + local_filename
        fileobj = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
        msg = f'Downloading {url}'
        if filesize is not None:
            msg += f" ({filesize/2**20:.1f} MB)"
        """with st.spinner(msg):
            with requests.get(url) as r:
                r.raise_for_status()  # Check for request success
                fileobj.write(r.content)"""
        @reactive.Effect
        def show_spinner():
            ui.output_text("spinner_message", msg)

        with requests.get(url) as r:
            r.raise_for_status()  # Check for request success
            fileobj.write(r.content)

        return fileobj
    except Exception as e:
        return None

def get_file_size(url):
    import requests
    response = requests.head(url)
    if 'Content-Length' in response.headers:
        file_size = int(response.headers['Content-Length'])
        return file_size
    else:
        return None

def change_mrc_map_crs_order(data, current_order, target_order=[1, 2, 3]):
    if current_order == target_order: return data
    map_crs_to_np_axes = {1:2, 2:1, 3:0}
    current_np_axes_order = [map_crs_to_np_axes[int(i)] for i in current_order]
    target_np_axes_order = [map_crs_to_np_axes[int(i)] for i in target_order]
    import numpy as np
    ret = np.moveaxis(data, current_np_axes_order, target_np_axes_order)
    return ret

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

def encode_numpy(img, hflip=False, vflip=False):
    if img.dtype != np.dtype('uint8'):
        vmin, vmax = img.min(), img.max()
        tmp = (255*(img-vmin)/(vmax-vmin)).astype(np.uint8)
    else:
        tmp = img
    if hflip:
        tmp = tmp[:, ::-1]
    if vflip:
        tmp = tmp[::-1, :]
    #import io, base64
    #from PIL import Image
    pil_img = Image.fromarray(tmp)
    buffer = io.BytesIO()
    pil_img.save(buffer, format="JPEG")
    encoded = base64.b64encode(buffer.getvalue()).decode()
    return f"data:image/jpeg;base64, {encoded}"

class Data:
    def __init__(self, twist, rise, csym, diameter, rotate=None, dx=None, apix_or_nqyuist=None, url=None, input_type="image"):
        self.input_type = input_type
        self.twist = twist
        self.rise = rise
        self.csym = csym
        self.diameter = diameter
        self.rotate = rotate
        self.dx = dx
        if self.input_type in ["PS", "PD"]:
            self.nyquist = apix_or_nqyuist
        else:
            self.apix = apix_or_nqyuist
        self.url = url

data_examples = [
    Data(twist=29.40, rise=21.92, csym=6, diameter=138, url="https://tinyurl.com/y5tq9fqa"),
    #Data(twist=36.0, rise=3.4, csym=1, diameter=20, dx=5, input_type="PS", apix_or_nqyuist=2.5, url="https://upload.wikimedia.org/wikipedia/en/b/b2/Photo_51_x-ray_diffraction_image.jpg")
]

def clear_twist_rise_csym_in_session_state():
    #if "twist" in st.session_state: del st.session_state["twist"]
    #if "rise" in st.session_state: del st.session_state["rise"]
    #if "csym" in st.session_state: del st.session_state["csym"]
    return

int_types = {'apply_helical_sym_0':0, 'apply_helical_sym_1':0, 'csym':1, 'csym_ahs_0':1, 'csym_ahs_1':1, 'do_random_embid_0':0, 'do_random_embid_1':0, 'fft_top_only':0, 'image_index_0':0, 'image_index_1':0, 'input_mode_0':1, 'input_mode_1':1, 'is_3d_0':0, 'is_3d_1':0, 'm_0':1, 'm_1':1, 'm_max':3, 'negate_0':0, 'negate_1':0, 'pnx':512, 'pny':1024, 'show_LL':1, 'show_LL_text':1, 'show_phase_diff':1, 'show_pwr':1, 'show_yprofile':1, 'transpose_0':0, 'transpose_1':0, 'share_url':0, 'show_qr':0, 'useplotsize':0}
float_types = {'angle_0':0, 'angle_1':0, 'apix_0':0, 'apix_1':0, 'apix_ahs_0':0, 'apix_ahs_1':0, 'apix_map_0':0, 'apix_map_1':0, 'apix_nyquist_0':0, 'apix_nyquist_1':0, 'az_0':0, 'az_1':0, 'ball_radius':0, 'cutoff_res_x':0, 'cutoff_res_y':0, 'diameter':0, 'dx_0':0, 'dx_1':0, 'dy_0':0, 'dy_1':0, 'fraction_ahs_0':0, 'fraction_ahs_1':0, 'length_ahs_0':0, 'length_ahs_1':0, 'mask_len_0':90, 'mask_len_1':90, 'mask_radius_0':0, 'mask_radius_1':0, 'noise_0':0, 'noise_1':0, 'resolution':0, 'rise':0, 'rise_ahs_0':0, 'rise_ahs_1':0, 'simuaz':0, 'simunoise':0, 'tilt':0, 'tilt_0':0, 'tilt_1':0, 'twist':0, 'twist_ahs_0':0, 'twist_ahs_1':0, 'width_ahs_0':0, 'width_ahs_1':1}
other_types = {'const_image_color':'', 'emd_id_0':'', 'emd_id_1':'', 'input_type_0':'image', 'input_type_1':'image', 'll_colors':'lime cyan violet salmon silver', 'url_0':'', 'url_0':''}

def get_direct_url(url):
    #import re
    if url.startswith("https://drive.google.com/file/d/"):
        hash = url.split("/")[5]
        return f"https://drive.google.com/uc?export=download&id={hash}"
    elif url.startswith("https://app.box.com/s/"):
        hash = url.split("/")[-1]
        return f"https://app.box.com/shared/static/{hash}"
    elif url.startswith("https://www.dropbox.com"):
        if url.find("dl=1")!=-1: return url
        elif url.find("dl=0")!=-1: return url.replace("dl=0", "dl=1")
        else: return url+"?dl=1"
    elif url.find("sharepoint.com")!=-1 and url.find("guestaccess.aspx")!=-1:
        return url.replace("guestaccess.aspx", "download.aspx")
    elif url.startswith("https://1drv.ms"):
        #import base64
        data_bytes64 = base64.b64encode(bytes(url, 'utf-8'))
        data_bytes64_String = data_bytes64.decode('utf-8').replace('/','_').replace('+','-').rstrip("=")
        return f"https://api.onedrive.com/v1.0/shares/u!{data_bytes64_String}/root/content"
    else:
        return url

def set_to_periodic_range(v, min=-180, max=180):
    if min <= v <= max: return v
    #from math import fmod
    tmp = fmod(v-min, max-min)
    if tmp>=0: tmp+=min
    else: tmp+=max
    return tmp

def mem_info():
    import_with_auto_install(["psutil"])
    #from psutil import virtual_memory
    mem = virtual_memory()
    mb = pow(2, 20)
    return (mem.total/mb, mem.available/mb, mem.used/mb, mem.percent)

def mem_quota():
    fqdn = get_hostname()
    if fqdn.find("heroku")!=-1:
        return 512  # MB
    username = get_username()
    if username.find("appuser")!=-1:    # streamlit share
        return 1024  # MB
    available_mem = mem_info()[1]
    return available_mem

def mem_used():
    import_with_auto_install(["psutil"])
    #from psutil import Process
    #from os import getpid
    mem = Process(getpid()).memory_info().rss / 1024**2   # MB
    return mem

def host_uptime():
    import_with_auto_install(["uptime"])
    from uptime import uptime
    t = uptime()
    if t is None: t = 0
    return t

def get_username():
    #from getpass import getuser
    return getuser()

def get_hostname():
    #import socket
    fqdn = socket.getfqdn()
    return fqdn

def is_hosted(return_host=False):
    hosted = False
    host = ""
    fqdn = get_hostname()
    if fqdn.find("heroku")!=-1:
        hosted = True
        host = "heroku"
    username = get_username()
    if username.find("appuser")!=-1:
        hosted = True
        host = "streamlit"
    if not host:
        host = "localhost"
    if return_host:
        return hosted, host
    else:
        return hosted

def read_mrc_data(mrc):
    # read mrc data
    mrc_data = mrcfile.open(mrc, 'r+')
    mrc_data.set_image_stack

    v_size=mrc_data.voxel_size
    nx=mrc_data.header['nx']
    ny=mrc_data.header['ny']
    nz=mrc_data.header['nz']
    apix=v_size['x']

    data=mrc_data.data
    return data, nx, ny

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
    
    # for debugging interpolation
    # st.write(xs)
    # xs = np.array(xs) * 0 + nx//2
    
    tck = splrep(ys,xs,s=20)

    new_xs = splev(ys,tck)
    
    if display:
        def render_spline():
            ui.output_text("spline_fit_msg", "Fitted spline:")
            p = create_fit_spline_figure(data, xs, ys, new_xs, apix)
            ui.output_plot("spline_fit_plot", p)

        reactive.Effect(render_spline)
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






