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

from plotly.graph_objects import Figure
from shinywidgets import output_widget, render_widget, render_plotly


import pandas as pd
from shiny import reactive, req
import compute
import image_trace
import shiny_test
import plotly.graph_objects as go

# may be unnecessary
from skimage.transform import resize_local_mean
import mrcfile
from shiny.types import FileInfo
import pathlib
import mrcfile
import matplotlib.pyplot as plt


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

ny, nx = 200, 200  # Get the dimensions of the loaded image
MAX_SIZE = 500 # temp value

urls = {
    "empiar-10940_job010": (
        "https://tinyurl.com/y5tq9fqa",
    )
}
url_key = "empiar-10940_job010"

# may not be needed
"""params_orig = reactive.value(None)
params_work = reactive.value(None)
apix_micrograph_auto = reactive.value(0)
apix_micrograph_auto_source = reactive.value(None)
apix_micrograph = reactive.value(0)
apix_particle = reactive.value(0)

data_all = reactive.value(None)
abundance = reactive.value([])
apix_class = reactive.value(0)
image_size = reactive.value(0)

displayed_class_ids = reactive.value([])
displayed_class_images = reactive.value([])
displayed_class_labels = reactive.value([])

min_len = reactive.value(0)
selected_image_indices = reactive.value([])
selected_images = reactive.value([])
selected_image_labels = reactive.value([])

selected_helices = reactive.value(([], [], 0))
retained_helices_by_length = reactive.value([])
pair_distances = reactive.value([]) """

data_all = reactive.value(None)
apix_reactive = reactive.value(1)
image_size = reactive.value(0)
displayed_class_images = reactive.value([])
displayed_class_labels = reactive.value([])
initial_selected_image_indices = reactive.value([])
selected_images = reactive.value([])
selected_image_labels = reactive.value([])


first_point = reactive.Value(None)
second_point = reactive.Value(None)

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
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 100%;",
                        class_="wrap-text"
                    ),
                    class_="button-container"
                ),
                ui.column(3,
                    ui.div(
                        ui.output_ui("col_four"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 100%;",
                        class_="wrap-text"
                    ),
                    class_="button-container"
                ),
                ui.column(3,
                    ui.div(
                        ui.output_ui("col_five"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px; height: 100%;",
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
        style="display: flex; flex-direction: column; height: 100%;"
    ),
    
)

def server(input, output, session):

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
        value = 200 # temp temp val
        apix = 1 # temp val
        nx = 5 # temp val
        ny = 5 #temp val
        helical_radius = 1 #temp val
        return ui.TagList(
                ui.input_action_button("col1_saveTwist", "Save twist/rise↶", class_="save-btn-success"),
                ui.input_numeric('csym', 'csym', value=6, min=1, step=1),
                ui.input_numeric('filament_diameter', 'Filament/tube diameter (Å)', value=value, min=1.0, max=1000.0, step=10.),
                ui.input_numeric("out_of_plane_tilt", "Out-of-plane tilt (°)", value=0.0, min=-90.0, max=90.0, step=1.0),
                ui.input_numeric('res_limit_x', 'Resolution limit - X (Å)', value=3*apix, min=2*apix, step=1),
                ui.input_numeric("res_limit_y", "Resolution limit - Y (Å)", value=2*apix, min=2*apix, step=1.0),
                ui.accordion(
                    ui.accordion_panel(
                        "Additional settings",
                        ui.input_checkbox("fft_top_only", "Only display the top half of FFT", value=False),
                        ui.input_checkbox("log_amp", "Log(amplitude)", value=True),
                        ui.input_text("const_image_color", "Flatten the PS/PD image in this color", value="", placeholder="white"),
                        ui.input_text("ll_colors", 'Layerline colors', value="lime cyan violet salmon silver"),
                        ui.input_numeric("hp_fraction", 'Fourier high-pass (%)', value=0.4 * 100, min=0.0, max=100.0, step=0.1),
                        ui.input_numeric("lp_fraction", 'Fourier low-pass (%)', value=0.0 * 100, min=0.0, max=100.0, step=10.0),
                        ui.input_numeric("pnx", 'FFT X-dim size (pixels)', value=512, min=min(nx, 128), step=2),
                        ui.input_numeric("pny", 'FFT Y-dim size (pixels)', value=1024, min=min(ny, 512), step=2),
                        value="add_settings"
                    ),
                ),
                ui.accordion(
                    ui.accordion_panel(
                        "Simulation",
                        ui.input_numeric("ball_radius", 'Gaussian radius (Å)', value=0.0, min=0.0, max=helical_radius, step=5.0),
                        value="simulation"
                    ),
                ),
                ui.input_checkbox("share_url", "Show/Reload sharable URL", value=False),
            )

    @output
    @render.ui
    def col_two():
        m_max_auto = 3 # temp val
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
            ui.input_numeric("m_max", "Max=", value=m_max_auto, min=1, step=1),
            ui.output_ui("show_choices")
        )

    @output
    @render.ui
    def col_three():
        twist = 0.0 # temp val
        return ui.TagList(
            ui.input_numeric("spinner_twist", "Twist (°)", value=twist, min=-180.0, max=180.0, step=1.0, width="300px"),
            ui.input_slider("slider_twist", "Twist (°)", min=-180.0, max=180.0, value=twist, step=0.01, width="300px"),
        )
        
    @output
    @render.ui
    def col_four():
        min_pitch = 0.0 # temp val
        pitch = 0.0 # temp val
        return ui.TagList(
            # TODO fix this line: ui.input_numeric("spinner_pitch", "Pitch (Å)", value=max(min_pitch, session.get("pitch", min_pitch)), min=min_pitch, step=1.0, width="300px"),
            ui.input_slider("slider_pitch", "Pitch (Å)", min=pitch/2.0, max=pitch*2.0, value=pitch, step=pitch*0.002, width="300px"),
        )
        
    @output
    @render.ui
    def col_five():
        min_rise = 0.0 # temp val
        max_rise = 0.0 # temp val
        rise = 0.0 # temp val
        return ui.TagList(
            ui.input_numeric("spinner_rise", "Rise (Å)", value=rise, min=min_rise, max=max_rise, step=1.0, width="300px"),
            ui.input_slider("slider_rise", "Rise (Å)", min=rise/2.0, max=min(max_rise, rise*2.0), value=rise, step=min(max_rise, rise*2.0)*0.001, width="300px")
        )

    #url = reactive.Value(urls[url_key][0])

    # code for the sidebar
    


    # when adding the code for sharable url, you can get rid of the url variables since it's only for uploading images and parsing

    @output
    @render.ui
    def conditional_3D():
        is_3d = input.is_3d()

        apix_auto = 0 # temp value
        negate_auto = 0 # temp value
        angle_auto = 1 # temp value
        dx_auto = 1 # temp value
        apix = 1 # temp value
        mask_radius_auto = 1 # temp value
        mask_len_percent_auto = 1 # temp value

        if is_3d:
            return ui.TagList( 
                ui.accordion(
                        ui.accordion_panel(
                            ui.p("Generate 2-D projection from the 3-D map"),
                            ui.input_checkbox("apply_helical_sym", "Apply helical symmetry", value=0),
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
                            ui.input_numeric("apix", "Pixel size (Å/pixel)", value=apix_auto, min=0.1, max=30.0, step=0.01),
                            ui.input_checkbox("transpose", "Transpose the image", value=negate_auto),
                            ui.input_checkbox("negate", "Invert the image contrast", value=negate_auto),
                            ui.input_numeric("angle", "Rotate (°)", value=-angle_auto, min=-180.0, max=180.0, step=1.0),
                            ui.input_numeric("dx", "Shift along X-dim (Å)", value=dx_auto*apix, min=-nx*apix, max=nx*apix, step=1.0),
                            ui.input_numeric("dy", "Shift along Y-dim (Å)", value=0.0, min=-ny*apix, max=ny*apix, step=1.0),
                            ui.input_numeric("mask_radius", "Mask radius (Å)", value=min(mask_radius_auto*apix, nx/2*apix), min=1.0, max=nx/2*apix, step=1.0),
                            ui.input_numeric("mask_len", "Mask length (%)", value=mask_len_percent_auto, min=10.0, max=100.0, step=1.0),
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
                        ui.input_numeric("apix", "Pixel size (Å/pixel)", value=apix_auto, min=0.1, max=30.0, step=0.01),
                        ui.input_checkbox("negate", "Invert the image contrast", value=negate_auto),
                        ui.input_checkbox("straightening", "Straighten the filament", value=False),
                        ui.input_numeric("angle", "Rotate (°)", value=-angle_auto, min=-180.0, max=180.0, step=1.0),
                        ui.input_numeric("dx", "Shift along X-dim (Å)", value=dx_auto*apix, min=-nx*apix, max=nx*apix, step=1.0),
                        ui.input_numeric("dy", "Shift along Y-dim (Å)", value=0.0, min=-ny*apix, max=ny*apix, step=1.0),
                        ui.input_numeric("mask_radius", "Mask radius (Å)", value=min(mask_radius_auto*apix, nx/2*apix), min=1.0, max=nx/2*apix, step=1.0),
                        ui.input_numeric("mask_len", "Mask length (%)", value=mask_len_percent_auto, min=10.0, max=100.0, step=1.0),
                        value="image_parameters_2"
                    )
                ),
            )

    # code added from the Helical Pitch Image Selection:
    @output
    @render.ui
    def sidebar_text():
        selection = input.input_mode_params()
        
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
                ui.input_checkbox("ignore_blank", "Ignore blank classes", value=True),
                ui.accordion(
                    ui.accordion_panel(
                        "Dynamic Image Selector",
                        ui.div(
                            ui.output_ui("dynamic_image_select_upload"),
                            style="max-height: 500px; overflow-y: auto;"
                        )
                    ),
                    id="acc_single", 
                    multiple=False,
                ),
                #ui.output_ui("dynamic_image_select_upload")
            )
        elif selection == "2":  # URL
            return ui.TagList(
                ui.input_text(
                    "url_params",
                    "Input a url of 2D image(s) or a 3D map:",
                    value=urls[url_key][0],
                ),
                ui.input_action_button("run", label="Run", style="width: 100%;"),
                ui.input_checkbox("is_3d", "The input ({nx}x{ny}x{nz}) is a 3D map", value=False),
                ui.input_checkbox("ignore_blank", "Ignore blank classes", value=True),
                ui.output_ui("conditional_3D"),
                #ui.output_ui("display_selected_images"),
                output_widget("display_micrograph"),
                output_widget("transformed_display_micrograph"),
                output_widget("plot_graph"),
                output_widget("acf_plot"),

            )
        elif selection == "3":  # EMD-xxxxx
            apix_auto = 0 # temp value
            negate_auto = 0 # temp value
            angle_auto = 1 # temp value
            dx_auto = 1 # temp value
            apix = 1 # temp value
            mask_radius_auto = 1 # temp value
            mask_len_percent_auto = 1 # temp value

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
                                ui.input_numeric("apix", "Pixel size (Å/pixel)", value=apix_auto, min=0.1, max=30.0, step=0.01),
                                ui.input_checkbox("transpose", "Transpose the image", value=negate_auto),
                                ui.input_checkbox("negate", "Invert the image contrast", value=negate_auto),
                                ui.input_numeric("angle", "Rotate (°)", value=-angle_auto, min=-180.0, max=180.0, step=1.0),
                                ui.input_numeric("dx", "Shift along X-dim (Å)", value=dx_auto*apix, min=-nx*apix, max=nx*apix, step=1.0),
                                ui.input_numeric("dy", "Shift along Y-dim (Å)", value=0.0, min=-ny*apix, max=ny*apix, step=1.0),
                                ui.input_numeric("mask_radius", "Mask radius (Å)", value=min(mask_radius_auto*apix, nx/2*apix), min=1.0, max=nx/2*apix, step=1.0),
                                ui.input_numeric("mask_len", "Mask length (%)", value=mask_len_percent_auto, min=10.0, max=100.0, step=1.0),
                                value="image_parameters_emd"
                            )
                        )             
                )
            # try to do random input, otherwise it's ok to leave out
        else:
            return ui.p("Please select an option to proceed.")

    @reactive.effect
    @reactive.event(data_all, input.ignore_blank)
    def get_displayed_class_images():
        if data_all() is not None:
            data = data_all()
            n = len(data)
            images = [data[i] for i in range(n)]
            image_size.set(max(images[0].shape))

            if input.ignore_blank():
                included = []
                included_images = []
                for i in range(n):
                    image = images[i]
                    if np.max(image) > np.min(image):
                        included.append(i)
                        included_images.append(image)
                images = included_images
            else:
                included = list(range(n))
            
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
            selected_images.set(
                [displayed_class_images()[i] for i in input.select_classes()]
            )
            selected_image_labels.set(
                [displayed_class_labels()[i] for i in input.select_classes()]
            )

    @reactive.Effect
    @reactive.event(input.run)
    def get_class2d():
        if input.input_mode_params() == "1":
            fileinfo = input.upload_classes()
            if fileinfo is not None:
                class_file = fileinfo[0]["datapath"]
                try:
                    data, apix = compute.get_class2d_from_file(class_file)
                    nx = data.shape[-1]
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
                nx = data.shape[-1]
                data_all.set(data)
                image_size.set(nx)
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
        #req(selected_images() is not None)
        #print("Selected image: ", selected_images())
        images = selected_images()
        if not images:  # Check for None or empty list
            #print("No selected images!")
            return px.scatter(title="No image selected")

        #print("Selected image: ", selected_images())

        nx, ny = selected_images()[0].shape[:2]

        fig = image_trace.plot_micrograph(
            micrograph=255 - selected_images()[0],
            title=f"Original image ({nx}x{ny})",
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
            #print("No selected images!")
            return px.scatter(title="No image selected")

        #print("Selected image: ", selected_images())

        nx, ny = selected_images()[0].shape[:2]

        image_to_display = (
            selected_images()[0] if input.negate() else 255 - selected_images()[0]
        )

        fig = image_trace.plot_micrograph(
            micrograph= image_to_display, #255 - selected_images()[0],
            title=f"Transformed image ({nx}x{ny})",
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
    

    @reactive.Calc
    def mask_parameters():
        images = selected_images()
        if not images:
            #print("No selected images!")
            return px.scatter(title="No image selected")
        
        data = selected_images()[0]

        straightening = input.straightening()
        radius_auto, mask_radius_auto = estimate_radial_range(data)
        mask_radius = input.mask_radius()
        mask_len_percent_auto = 100.0 if straightening else 90.0
        mask_len_fraction = input.mask_len() / 100.0
        return mask_radius, mask_len_fraction

    @render_plotly
    @reactive.event(selected_images)
    def plot_graph():
        images = selected_images()
        if not images:
            #print("No selected images!")
            return px.scatter(title="No image selected")
        
        data = selected_images()[0]
        mask_radius, mask_len_fraction = mask_parameters()
        
        # Calculate data for the plot
        x = np.arange(-nx // 2, nx // 2) * apix_reactive()
        ymax = np.max(data, axis=0)
        ymean = np.mean(data, axis=0)
        
        # Create Plotly figure
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

        # Update layout
        fig.update_layout(
            xaxis_title="x (Å)",
            yaxis_title="Pixel value",
            showlegend=True,
            legend=dict(x=0.9, y=1.1),
            height=400
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
            hovermode="x unified"
        )

        return fig



app = App(app_ui, server)



# Placeholder function for estimating radial range
def estimate_radial_range(data, thresh_ratio=0.1):
    # Simulated radius_auto and mask_radius_auto
    radius_auto = 10  # Example value
    mask_radius_auto = 20  # Example value
    return radius_auto, mask_radius_auto

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






