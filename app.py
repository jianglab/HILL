import numpy as np
import pandas as pd

from shinywidgets import render_plotly

from shiny import App, reactive, render, ui

from shinywidgets import output_widget, render_plotly
import pandas as pd
from shiny import reactive, req
import compute
import util
import plotly.graph_objects as go

from finufft import nufft2d2

from numba import set_num_threads, prange

import pandas as pd

from scipy.spatial.transform import Rotation as R
from scipy.ndimage import map_coordinates

import helicon

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
nx = reactive.value(256)  # default value
ny = reactive.value(256)  # default value
nz = reactive.value(256)    # default value
MAX_SIZE = 500 # temp value


urls = {
    "empiar-10940_job010": (
        "https://tinyurl.com/y5tq9fqa",
    )
}
url_key = "empiar-10940_job010"

# may not be needed
input_data = reactive.value(None)
apix_from_file = reactive.value(1.0)
data_all_2d = reactive.value([])
data_all_2d_labels = reactive.value([])
initial_selected_image_indices = reactive.value([])
selected_images = reactive.value([])
selected_image_labels = reactive.value([])
image_display_size = reactive.value(128)

first_point = reactive.Value(None)
second_point = reactive.Value(None)

apix_auto = reactive.value(0.0)
negate_auto = reactive.value(0.0)
angle_auto = reactive.value(0.0)
dx_auto = reactive.value(0.0)
dy_auto = reactive.value(0.0)
mask_radius_auto = reactive.value(0.0)
mask_len_percent_auto = reactive.value(0.0)
radius_auto = reactive.value(0.0)
auto_transformation_estimated = reactive.value(False)

apix_value = reactive.value(0.0)
angle_value = reactive.value(0.0)
dx_value = reactive.value(0.0)
dy_value = reactive.value(0.0)
mask_radius_value = reactive.value(0.0)
mask_len_value = reactive.value(0.0)

helical_radius = reactive.value(1)
m_max_auto_reactive = reactive.value(3)
twist = reactive.value(29.40)
min_pitch = reactive.value(0.0)
pitch = reactive.value(268.41)
rise = reactive.value(21.92)
min_rise = reactive.value(2.0)
max_rise = reactive.value(268.41)

is_3d_reactive = reactive.value(None)
data = reactive.value(None)
data_2d_transformed = reactive.value(None)
data_2d_transformed_masked = reactive.value(None)

transform_params = reactive.Value(None)
change_image = reactive.value(0.0)
prev_dx_val = reactive.Value(-1.0)
curr_dx_val = reactive.Value(1.0)

key_emd_id = reactive.Value(None)
emd_id = reactive.value("")
csym = reactive.Value(None)
max_map_size = reactive.Value(None)
max_map_dim = reactive.Value(None)
stop_map_size = reactive.Value(None)


emdb_df_original = reactive.value(None)
emdb_df = reactive.value(None)
emd_url = reactive.value("")
msg_reactive = reactive.value("")

maps = reactive.value([])
map_xyz_projections = reactive.value([])
map_xyz_projection_title = reactive.value("Map XYZ projections:")
map_xyz_projection_labels = reactive.value([])
map_xyz_projection_display_size = reactive.value(128)
recalculate_emd_params = reactive.Value(False)

ps_data = reactive.value(None)
pd_data = reactive.value(None)

user_inputs = reactive.value({
    'apix': None,
    'angle': None,
    'dx': None,
    'dy': None,
    'mask_radius': None
})


app_ui = ui.page_fluid(
    ui.tags.style("""       
        .scrollable-sidebar {
            max-height: 100vh; /* Adjust height as needed */
            overflow-y: auto; /* Enables vertical scrolling */
            border: 1px solid #ccc; /* Optional: for visual clarity */
            padding: 10px;
        }
        .scrollable-container {
            max-height: 100vh; 
            overflow-y: auto;  
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
                        selected="2",
                    ),
                    value="input_mode"
                )
            ),
            ui.output_ui("conditional_input_uis"),
            ui.output_ui("conditional_3d_transformation_uis"),
            #ui.input_action_button("run", label="Run"),
            ui.output_ui("input_display_uis"),
            output_widget("display_selected_image"),
            ui.accordion(
                ui.accordion_panel(
                    ui.p("Image Parameters"),
                    ui.input_radio_buttons("input_type", "Input is:", choices=["image", "PS", "PD"], inline=True),
                    ui.output_ui("img_2d_transformation_uis"),
                    value="image_parameters_emd"
                )
            ),
            output_widget("display_transformed_data"),
            output_widget("plot_radial_profile"),
            output_widget("acf_plot"),
            #class_="custom-sidebar",
            class_="scrollable-sidebar",
            width = "20vw"            
        ),
        ui.div(
            ui.row(
                ui.column(12,
                    ui.div(
                        ui.h2("HILL: Helical Indexing using Layer Lines"),
                        style="text-align: center; margin-bottom: 20px;"
                    )
                )
            ),
            ui.div(
                ui.layout_columns(
                    ui.div(
                        ui.input_numeric('csym', 'csym', value=6, min=1, step=1, update_on="blur"),
                        ui.input_numeric('filament_diameter', 'Filament/tube diameter (Å)', value=100*2, min=1.0, max=1000.0, step=10, update_on="blur"),
                        ui.input_numeric("out_of_plane_tilt", "Out-of-plane tilt (°)", value=0.0, min=-90.0, max=90.0, step=1.0, update_on="blur"),
                        ui.input_numeric('res_limit_x', 'Resolution limit - X (Å)', value=round(3*1.0,4), min=2*1.0, step=1.0, update_on="blur"),
                        ui.input_numeric("res_limit_y", "Resolution limit - Y (Å)", value=round(2*1.0,4), min=2*1.0, step=1.0, update_on="blur"),
                        ui.accordion(
                            ui.accordion_panel(
                                "Additional settings",
                                ui.input_checkbox("fft_top_only", "Only display the top half of FFT", value=False),
                                ui.input_checkbox("log_amp", "Log(amplitude)", value=True),
                                ui.input_text("const_image_color", "Flatten the PS/PD image in this color", value="", placeholder="white"),
                                ui.input_text("ll_colors", 'Layerline colors', value="lime cyan violet salmon silver"),
                                ui.input_numeric("hp_fraction", 'Fourier high-pass (%)', value=0.4 * 100, min=0.0, max=100.0, step=0.1, update_on="blur"),
                                ui.input_numeric("lp_fraction", 'Fourier low-pass (%)', value=0.0 * 100, min=0.0, max=100.0, step=10.0, update_on="blur"),
                                ui.input_numeric("pnx", 'FFT X-dim size (pixels)', value=512, min=128, step=2, update_on="blur"),
                                ui.input_numeric("pny", 'FFT Y-dim size (pixels)', value=1024, min=512, step=2, update_on="blur"),
                                value="add_settings"
                            ),
                        ),
                        ui.accordion(
                            ui.accordion_panel(
                                "Simulation",
                                ui.input_numeric("ball_radius", 'Gaussian radius (Å)', value=0.0, min=0.0, max=100.0, step=5.0, update_on="blur"),
                                value="simulation"
                            ),
                        ),
                        ui.input_checkbox("share_url", "Show/Reload sharable URL", value=False),
                        #ui.output_ui("col_one"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px;",
                        class_="wrap-text"
                    ),
                    ui.div(
                        ui.output_ui("col_two"),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px;",
                        class_="wrap-text"
                    ),
                    ui.div(
                        ui.layout_columns(
                            ui.output_ui("helical_param_inputs"),
                            output_widget("main_plots"),
                            col_widths=(12,12),
                            #row_heights=(1,4)
                        ),
                        style="text-align: left; border: 1px solid #ddd; padding: 10px;",
                    ),
                    col_widths=(2,2,8)
                ),
                class_="scrollable-container"
            ),
        ),
        ui.row(
            ui.column(12,
                ui.div(
                    ui.markdown("*Developed by the [Jiang Lab@Penn State](https://jianglab.science.psu.edu). Report problems to [HILL@GitHub](https://github.com/jianglab/hill/issues)*"),
                )
            )
        ),
        helicon.shiny.setup_ajdustable_sidebar(width="20vw"),
    ),
)


def server(input, output, session):

    @output
    @render.text
    def m_max_auto():
        rise = input.rise()
        res_limit_y = input.res_limit_y()
        return int(np.floor(np.abs(rise/res_limit_y))) + 3

    @reactive.calc
    def m_groups():
        req(input.filament_diameter()>0)
        req(input.res_limit_y()>0)
        print("Updating mgroups")
        return compute.compute_layer_line_positions(
            twist=input.spinner_twist(),
            rise=input.spinner_rise(),
            csym=input.csym(),
            radius=input.filament_diameter()/2,
            tilt=input.out_of_plane_tilt(),
            cutoff_res=input.res_limit_y(),
            m_max=input.m_max()
        )

    @reactive.effect
    @reactive.event(input.spinner_twist)
    def _():
        if twist() != input.spinner_twist() and input.use_twist_pitch() == "Twist":
            twist.set(input.spinner_twist())
            pitch.set(round(360/np.abs(twist()) * rise(),2))
            max_rise.set(pitch())

    @reactive.effect
    @reactive.event(input.slider_twist)
    def _():
        if twist() != input.slider_twist() and input.use_twist_pitch() == "Twist":
            twist.set(input.slider_twist())
            pitch.set(round(360/np.abs(twist()) * rise(),2))
            max_rise.set(pitch())

    @reactive.effect
    @reactive.event(input.spinner_pitch)
    def _():
        if pitch() != input.spinner_pitch() and input.use_twist_pitch() == "Pitch":
            pitch.set(input.spinner_pitch())
            twist.set(round(360/np.abs(pitch()) * rise(),2))
            max_rise.set(pitch())

    @reactive.effect
    @reactive.event(input.slider_pitch) 
    def _():
        if pitch() != input.slider_pitch() and input.use_twist_pitch() == "Pitch":
            pitch.set(input.slider_pitch())
            twist.set(round(360/np.abs(pitch()) * rise(),2))
            max_rise.set(pitch())

    @reactive.effect
    @reactive.event(input.spinner_rise)
    def _():
        if rise() != input.spinner_rise():
            rise.set(input.spinner_rise())
            if input.use_twist_pitch() == "Twist":
                pitch.set(round(360/np.abs(twist()) * input.spinner_rise(),2))
            else:
                twist.set(round(360/np.abs(pitch()) * input.spinner_rise(),2))
            min_pitch.set(rise())

    @reactive.effect
    @reactive.event(input.slider_rise)
    def _():
        if rise() != input.slider_rise():
            rise.set(input.slider_rise())
            if input.use_twist_pitch() == "Twist":
                pitch.set(round(360/np.abs(twist()) * input.slider_rise(),2))
            else:
                twist.set(round(360/np.abs(pitch()) * input.slider_rise(),2))
            min_pitch.set(rise())
    
    @reactive.effect
    @reactive.event(input.apix)
    def set_min_rise():
        min_rise.set(round(2*input.apix(),2))

    @output
    @render.ui
    def show_choices():
        m_max = input.m_max()
        ms = dict()
        for m in range(m_max, -m_max-1, -1):
            ms[m] = str(m)
        checkboxes = ui.input_checkbox_group("ms","m:",ms, selected=[-1, 0, 1])
        return checkboxes

    @output
    @render.ui
    def col_one():
        return [
                ui.input_numeric('csym', 'csym', value=6, min=1, step=1, update_on="blur"),
                ui.input_numeric('filament_diameter', 'Filament/tube diameter (Å)', value=radius_auto()*2, min=1.0, max=1000.0, step=10, update_on="blur"),
                ui.input_numeric("out_of_plane_tilt", "Out-of-plane tilt (°)", value=0.0, min=-90.0, max=90.0, step=1.0, update_on="blur"),
                ui.input_numeric('res_limit_x', 'Resolution limit - X (Å)', value=round(3*apix_from_file(),4), min=2*apix_from_file(), step=1, update_on="blur"),
                ui.input_numeric("res_limit_y", "Resolution limit - Y (Å)", value=round(2*apix_from_file(),4), min=2*apix_from_file(), step=1.0, update_on="blur"),
                ui.accordion(
                    ui.accordion_panel(
                        "Additional settings",
                        ui.input_checkbox("fft_top_only", "Only display the top half of FFT", value=False),
                        ui.input_checkbox("log_amp", "Log(amplitude)", value=True),
                        ui.input_text("const_image_color", "Flatten the PS/PD image in this color", value="", placeholder="white"),
                        ui.input_text("ll_colors", 'Layerline colors', value="lime cyan violet salmon silver"),
                        ui.input_numeric("hp_fraction", 'Fourier high-pass (%)', value=0.4 * 100, min=0.0, max=100.0, step=0.1, update_on="blur"),
                        ui.input_numeric("lp_fraction", 'Fourier low-pass (%)', value=0.0 * 100, min=0.0, max=100.0, step=10.0, update_on="blur"),
                        ui.input_numeric("pnx", 'FFT X-dim size (pixels)', value=512, min=min(nx(),128), step=2, update_on="blur"),
                        ui.input_numeric("pny", 'FFT Y-dim size (pixels)', value=1024, min=min(ny(), 512), step=2, update_on="blur"),
                        value="add_settings"
                    ),
                ),
                ui.accordion(
                    ui.accordion_panel(
                        "Simulation",
                        ui.input_numeric("ball_radius", 'Gaussian radius (Å)', value=0.0, min=0.0, max=helical_radius(), step=5.0, update_on="blur"),
                        value="simulation"
                    ),
                ),
                ui.input_checkbox("share_url", "Show/Reload sharable URL", value=False),
        ]

    @output
    @render.ui
    def col_two():
        return [
            ui.input_action_button("copy_twist", "Copy twist/rise◀"),
            ui.input_radio_buttons("use_twist_pitch", "Use twist/pitch:", choices=["Twist", "Pitch"], selected="Twist"),
            ui.h5("Display:"),
            ui.input_checkbox("PS", "PS", value=True),
            ui.input_checkbox("YP", "YP", value=True),
            ui.input_checkbox("Phase", "Phase", value=False),
            ui.input_checkbox("PD", "PD", value=True),
            ui.input_checkbox("Color", "Color", value=True),
            ui.input_checkbox("LL", "LL", value=True),
            ui.input_checkbox("LLText", "LLText", value=True),
            ui.h5("m:"),
            ui.input_numeric("m_max", "Max=", value=m_max_auto_reactive(), min=1, step=1, update_on="blur"),
            ui.output_ui("show_choices")
        ]

    @output
    @render.ui
    def col_three():
        return [
            ui.input_numeric("spinner_twist", "Twist (°)", value=twist(), min=-180.0, max=180.0, step=1.0, width="300px", update_on="blur"),
            ui.input_slider("slider_twist", "Twist (°)", min=-180.0, max=180.0, value=twist(), step=0.01, width="300px"),
        ]
        
    @output
    @render.ui
    def col_four():
        return [
            ui.input_numeric("spinner_pitch", "Pitch (Å)", value=max(min_pitch(), pitch()), min=min_pitch(), step=1.0, width="300px", update_on="blur"),
            ui.input_slider("slider_pitch", "Pitch (Å)", min=pitch()/2.0, max=pitch()*2.0, value=pitch(), step=pitch()*0.002, width="300px"),
        ]
        
    @output
    @render.ui
    def col_five():
        return [
            ui.input_numeric("spinner_rise", "Rise (Å)", value=rise(), min=min_rise(), max=max_rise(), step=1.0, width="300px", update_on="blur"),
            ui.input_slider("slider_rise", "Rise (Å)", min=min_rise(), max=min(max_rise(), rise()*2.0), value=rise(), step=min(max_rise(), rise()*2.0)*0.001, width="300px")
        ]
    
    @output
    @render.ui
    def helical_param_inputs():
        return ui.layout_columns(
            ui.input_numeric("spinner_twist", "Twist (°)", value=twist(), min=-180.0, max=180.0, step=1.0, width="300px", update_on="blur"),
            ui.input_numeric("spinner_pitch", "Pitch (Å)", value=max(min_pitch(), pitch()), min=min_pitch(), step=1.0, width="300px", update_on="blur"),
            ui.input_numeric("spinner_rise", "Rise (Å)", value=rise(), min=min_rise(), max=max_rise(), step=1.0, width="300px", update_on="blur"),
            ui.input_slider("slider_twist", "Twist (°)", min=-180.0, max=180.0, value=twist(), step=0.01, width="300px"),
            ui.input_slider("slider_pitch", "Pitch (Å)", min=pitch()/2.0, max=pitch()*2.0, value=pitch(), step=round(pitch()*0.002,2), width="300px"),
            ui.input_slider("slider_rise", "Rise (Å)", min=min_rise(), max=min(max_rise(), rise()*2.0), value=rise(), step=round(min(max_rise(), rise()*2.0)*0.001,2), width="300px"),
            col_widths=(4,4,4,4,4,4)
        )
    
    main_plot = reactive.value(go.FigureWidget().set_subplots(rows=1,cols=3, subplot_titles=("", "", ""), column_widths=[0.4,0.2,0.4], vertical_spacing=0,horizontal_spacing=0))
    
    @render_plotly
    def main_plots():
        return main_plot()
    
    @reactive.effect
    @reactive.event(ps_data)
    def update_main_plots():
        req(input.res_limit_y()>0)
        req(input.res_limit_x()>0)
        req(ps_data() is not None)
        print('computing ll position')
        print(input.res_limit_y())
        fig = main_plot()
        data_ps = ps_data()
        yprofile=np.mean(data_ps,axis=1)
        yprofile/=yprofile.max()
        data_pd = pd_data()
        
        print("entered main plots")
        print(input.ms())

        if len(fig.data) > 0:
            print("update main plots")
            fig.update_traces(z=data_ps,selector=dict(name="ps"))
            fig.update_traces(x=yprofile,selector=dict(name="yprofile"))
            fig.update_traces(z=data_pd,selector=dict(name="pd"))
        else:
            print("create main plots")
            ny, nx = data_ps.shape
            dsy = 1/(ny//2*input.res_limit_y())
            dsx = 1/(nx//2*input.res_limit_x())
            x_range = (-(nx//2+0.5)*dsx, (nx//2-0.5)*dsx)
            if input.fft_top_only():
                y_range = (-(ny//2 * 0.01)*dsy, (ny//2-0.5)*dsy)
            else:
                y_range = (-(ny//2+0.5)*dsy, (ny//2-0.5)*dsy)
            
            if input.Color():
                plot_color = "viridis"
            else:
                plot_color = "gray"

            xs = np.linspace(x_range[0],x_range[1],input.pnx())
            ys = np.linspace(y_range[0],y_range[1],input.pny())
            print("xs:", x_range)
            print("ys:", y_range)

            fig.add_trace(
                go.Heatmap(
                    name="ps",
                    z=data_ps,
                    x=xs,
                    y=ys,
                    xaxis="x1",
                    yaxis="y1",
                    colorscale=plot_color,
                    showscale=False,
                    ),
                row=1,col=1
            )
            fig.add_trace(
                go.Line(
                    name="yprofile",
                    x=yprofile,
                    y=ys,
                    xaxis="x2",
                    yaxis="y1",
                    ),
                row=1,col=2
            )
            fig.add_trace(
                go.Heatmap(
                    name="pd",
                    z=data_pd,
                    x=xs,
                    y=ys,
                    xaxis="x3",
                    yaxis="y1",
                    colorscale=plot_color,
                    showscale=False,
                    ),
                row=1,col=3
            )

            fig.update_xaxes(showticklabels=False, range=[x_range[0],x_range[1]])
            fig.update_yaxes(showticklabels=False, range=[y_range[0],y_range[1]])
            fig.update_xaxes(range=[0,1], row=1,col=2)
            
            fig.update_xaxes(showspikes=True, spikemode="across")
            fig.update_yaxes(showspikes=True, spikemode="across")
            fig.update_layout(
                hoversubplots="axis",
                yaxis_title="",
                hovermode="y unified",
                showlegend=False,
                height=input.pny(),
                autosize=True,
                coloraxis_showscale=False,
            )
            fig.layout.xaxis1.title="Power Spectrum"
            fig.layout.xaxis2.title="Y-profile"
            fig.layout.xaxis3.title="Phase Difference"
            fig.update_traces(yaxis='y1')

            # curr_m_groups=compute_layer_line_positions(
            #     twist=input.spinner_twist(),
            #     rise=input.spinner_rise(),
            #     csym=input.csym(),
            #     radius=input.filament_diameter()/2,
            #     tilt=input.out_of_plane_tilt(),
            #     cutoff_res=input.res_limit_y(),
            #     m_max=input.m_max()
            # )

            curr_m_groups = m_groups()
            #print("Current m_groups:")
            #print(curr_m_groups[0]["LL"])
            
            curr_ll_colors=input.ll_colors()
            if input.LL() and max(curr_m_groups[0]["LL"][0])>0:
                x, y, n = curr_m_groups[0]["LL"]
                tmp_x = np.sort(np.unique(x))
                width = np.mean(tmp_x[1:]-tmp_x[:-1])
                height = width/5
                shapes = []
                for mi, m in enumerate(curr_m_groups.keys()):
                    if str(m) not in input.ms(): continue
                    xs, ys, bessel_order = curr_m_groups[m]["LL"]
                    tags = [m, bessel_order]
                    color = curr_ll_colors.split()[abs(m)%len(curr_ll_colors)]
                    fills = [n%2*1.0 for n in bessel_order]
                    if input.LLText():
                        texts = [str(int(n)) for n in bessel_order]
                        for (x,y,text) in zip(xs,ys,texts):
                            shapes += [
                                dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, label=dict(text=text, font=dict(color=color)), opacity=0, line_color=color, name='ll', xref="x1",yref="y1"),
                                dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, label=dict(text=text, font=dict(color=color)), opacity=0, line_color=color, name='ll', xref="x3",yref="y1"),
                            ]
                    else:
                        for (x,y,fill) in zip(xs,ys,fills):
                            if fill:
                                shapes += [
                                    dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, fillcolor=color, name='ll', xref="x1",yref="y1"),
                                    dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, fillcolor=color, name='ll', xref="x3",yref="y1"),
                                ]
                            else:
                                shapes += [
                                    dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, name='ll', xref="x1",yref="y1"),
                                    dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, name='ll', xref="x3",yref="y1"),
                                ]
            
                fig.layout.shapes = shapes
        #return fig

    @reactive.effect
    @reactive.event(m_groups, input.ms, input.LL, input.LLText, input.ll_colors)
    def update_layer_lines():
        print("Updating layer lines")
        print("input.LL():", input.LL())
        print("input.LLText():", input.LLText())
        fig = main_plot()
        if input.LL():
            if len(fig.data) > 0:
                curr_m_groups = m_groups()
                #print("Current m_groups:")
                #print(curr_m_groups[0]["LL"])
                
                curr_ll_colors=input.ll_colors()
                if max(curr_m_groups[0]["LL"][0])>0:
                    x, y, n = curr_m_groups[0]["LL"]
                    tmp_x = np.sort(np.unique(x))
                    width = np.mean(tmp_x[1:]-tmp_x[:-1])
                    height = width/5
                    shapes = []
                    for mi, m in enumerate(curr_m_groups.keys()):
                        if str(m) not in input.ms(): continue
                        xs, ys, bessel_order = curr_m_groups[m]["LL"]
                        tags = [m, bessel_order]
                        color = curr_ll_colors.split()[abs(m)%len(curr_ll_colors)]
                        fills = [n%2*1.0 for n in bessel_order]
                        if input.LLText():
                            texts = [str(int(n)) for n in bessel_order]                
                            for (x,y,text) in zip(xs,ys,texts):
                                shapes += [
                                    dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, label=dict(text=text, font=dict(color=color)), opacity=0, line_color=color, fillcolor=color, name='ll', xref="x1",yref="y1"),
                                    dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, label=dict(text=text, font=dict(color=color)), opacity=0, line_color=color, fillcolor=color, name='ll', xref="x3",yref="y1"),
                                ]
                        else:
                            for (x,y,fill) in zip(xs,ys,fills):
                                if fill:
                                    shapes += [
                                        dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, fillcolor=color, name='ll', xref="x1",yref="y1"),
                                        dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, fillcolor=color, name='ll', xref="x3",yref="y1"),
                                    ]
                                else:
                                    shapes += [
                                        dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, name='ll', xref="x1",yref="y1"),
                                        dict(type='circle',x0=x-width/2, y0=y-height/2, x1=x+width/2, y1=y+height/2, line_color=color, name='ll', xref="x3",yref="y1"),
                                    ]
                
                    fig.layout.shapes = shapes
        else:
            print("Removing layer lines")
            #fig.update_layout(shapes=[])
            fig.layout.shapes = []

    @output
    @render.ui
    def conditional_3d_transformation_uis(): 
        if input.is_3d():
            return [
                ui.accordion(
                    ui.accordion_panel(
                        ui.p("Generate 2-D projection from the 3-D map"),
                        ui.input_checkbox("apply_helical_sym", "Apply helical symmetry", value=0),
                        ui.output_ui("apply_helical_sym_uis"),
                        ui.input_numeric("az", "Rotation around the helical axis (°):", min=0.0, max=360., value=0.0, step=1.0, update_on="blur"),
                        ui.input_numeric("tilt", "Tilt (°):", min=-180.0, max=180., value=0.0, step=1.0, update_on="blur"),
                        ui.input_numeric("gauss_noise_std", "Add Gaussian noise (σ):", min=0.0, value=0.0, step=0.5, update_on="blur"),
                        ui.input_action_button("generate_2d_projection", "Generate 2D projection"),
                        value="2D_projection",
                        open=False
                    )
                ),
            ]
        else:
            return None

    @output
    @render.ui
    def img_2d_transformation_uis():
        value = input.input_type()
        if value in ['image']:
            return [
                ui.input_numeric(
                    "apix", 
                    "Pixel size (Å/pixel)", 
                    #value=apix_from_file(),
                    value=1.0,
                    min=0.1, 
                    max=30.0, 
                    step=0.01,
                    update_on="blur"
                ),
                #ui.output_ui("add_transpose"),
                ui.input_checkbox("negate", "Invert the image contrast", value=negate_auto()),
                ui.input_checkbox("straightening", "Straighten the filament", value=False),
                ui.input_numeric(
                    "angle", 
                    "Rotate (°)", 
                    value=angle_auto(),
                    min=-180.0, 
                    max=180.0, 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "dx", 
                    "Shift along X-dim (Å)", 
                    value=dx_auto(),
                    min=-500.0, 
                    max=500.0, 
                    # min=-nx()*apix_from_file(), 
                    # max=nx()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "dy", 
                    "Shift along Y-dim (Å)", 
                    value=dy_auto(),
                    min=-500.0, 
                    max=500.0, 
                    # min=-ny()*apix_from_file(), 
                    # max=ny()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "mask_radius", 
                    "Mask radius (Å)", 
                    value=mask_radius_auto(),
                    min=1.0, 
                    max=500.0, 
                    # max=nx()/2*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "mask_len", 
                    "Mask length (%)", 
                    value=mask_len_percent_auto(), 
                    min=10.0, 
                    max=100.0, 
                    step=1.0,
                    update_on="blur"
                ),
            ]
        elif value in ['PS', 'PD']:
            return [
                ui.input_numeric(
                    "apix", 
                    "Nyquist res (Å)", 
                    #value=2*apix_from_file(),
                    value=1.0,
                    min=0.1, 
                    max=30.0, 
                    step=0.01,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "angle", 
                    "Rotate (°)", 
                    value=angle_auto(),
                    min=-180.0, 
                    max=180.0, 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "dx", 
                    "Shift along X-dim (Å)", 
                    value=dx_auto(),
                    min=-nx()*apix_from_file(), 
                    max=nx()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "dy", 
                    "Shift along Y-dim (Å)", 
                    value=dy_auto(),
                    min=-ny()*apix_from_file(), 
                    max=ny()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
            ]

    @render.ui
    @reactive.event(input.is_3d, input.input_mode_params)
    def add_transpose():
        transpose_auto = input.input_mode_params() not in [2, 3] and nx() > ny()
        if not input.is_3d():
            return ui.input_checkbox("transpose", "Transpose the image", value=transpose_auto),

    # @reactive.effect
    # @reactive.event(selected_images)
    # def set_transformation_parameters_auto():
    #     req(selected_images())
    #     # if input.input_mode_params() == "3":
    #     #     angle_auto.set(0.0)
    #     #     dx_auto.set(0.0)
    #     #     dy_auto.set(0.0)
    #     #     return
        
    #     mode = input.input_type() 

    #     if mode in ["PS", "PD"]:
    #         # angle_auto.set(0.0)
    #         # dx_auto.set(0.0)
    #         # dy_auto.set(0.0)
    #         ui.update_numeric('angle', value=0.0)
    #         ui.update_numeric('dx', value=0.0)
    #         ui.update_numeric('dy', value=0.0)
    #     else:
    #         print("Estimating auto transformation parameters...")
    #         ang_auto_v, dy_auto_v, diameter = helicon.estimate_helix_rotation_center_diameter(selected_images()[0], estimate_center=True, estimate_rotation=True)
    #         # angle_auto.set(round(ang_auto_v+90,2))
    #         # dy_auto.set(round(dy_auto_v*input.apix(),2))
    #         # mask_radius_auto.set(round(diameter/2*input.apix(),2))
    #         print("Estimated angle:", ang_auto_v+90)
    #         print("input angle before:", input.angle())
    #         ui.update_numeric('angle', value=round(ang_auto_v+90,2))
    #         print("input angle after:", input.angle())
    #         ui.update_numeric('dy', value=round(dy_auto_v*input.apix(),2))
    #         ui.update_numeric('mask_radius', value=round(diameter/2*input.apix(),2))
    #         mask_len_percent_auto_default = 90.0
    #         if input.straightening():
    #             mask_len_percent_auto_default = 100.0
    #         #mask_len_percent_auto.set(mask_len_percent_auto_default)
    #         ui.update_numeric('mask_len', value=mask_len_percent_auto_default)
    #     print(selected_images())
    #     print("Auto transformation parameters set: angle=", input.angle())

    @reactive.effect
    @reactive.event(selected_images, input.apix, input.angle, input.dx, input.dy, input.mask_radius)
    def set_data_2d_transformed():
        req(input.input_type() in ["image"])
        req(len(selected_images())>0)
        req(input.apix()!=0)
        if (input.angle() or input.dx() or input.dy() or input.negate() or input.mask_radius()):
            data_2d_transformed.set(compute.mask_2d_filament(compute.transform_2d_filament(selected_images()[0],input.angle(),input.dx(),input.dy(),input.negate(),input.apix()),input.mask_radius(), input.apix(), input.mask_len()/100.0))
    
    @reactive.effect
    @reactive.event(data_2d_transformed)
    def set_mask_radius_auto():
        req(data_2d_transformed() is not None) 
        input_type = input.input_type()
        if input_type in ["image"]:
            radius_auto_v, _ = compute.estimate_radial_range(data_2d_transformed(), thresh_ratio=0.1)
            #radius_auto.set(round(radius_auto_v*input.apix(),2))
            ui.update_numeric('filament_diameter', value=round(radius_auto_v*input.apix(),2))
    
    @reactive.effect
    def get_ps_pd_data():
        req(data_2d_transformed() is not None)
        req(input.apix()>0)
        req(input.res_limit_y()>0)
        req(input.res_limit_x()>0)
        pwr, phase = compute_power_spectra(data_2d_transformed(), apix=input.apix(), cutoff_res=(input.res_limit_y(), input.res_limit_x()), 
                output_size=(input.pny(), input.pnx()), log=True, low_pass_fraction=0, high_pass_fraction=0.4*100)
        ps_data.set(pwr)
        pd_data.set(compute.compute_phase_difference_across_meridian(phase))
    
    @output
    @render.ui
    def conditional_input_uis():
        selection = input.input_mode_params()
        if selection == "1":  # Upload
            return [
                ui.p("Upload a mrc or mrcs file"),
                ui.input_file(
                    "img_file_upload",
                    "Upload the image file (.mrcs, .mrc, .map)",
                    accept=[".mrcs", ".mrc", ".map",".map.gz"],
                    placeholder="image file",
                ),
                ui.input_checkbox("is_3d", "The input is a 3D map", value=False),
                # ui.input_checkbox("ignore_blank", "Ignore blank classes", value=True),
            ]
        elif selection == "2":  # URL
            return [
                ui.input_text(
                    "img_file_url",
                    "Input a url of 2D image(s) or a 3D map:",
                    value=urls[url_key][0],
                ),
                ui.input_checkbox("is_3d", "The input is a 3D map", value=False),
            ]
        elif selection == "3":  # EMD-xxxxx
            return [
                    ui.input_text("img_file_emd_id", "Input an EMDB ID (emd-xxxxx):", value="emd-10499"),
                    #ui.output_ui("get_link_emd"),
                    #ui.p("You have selected to use an EMD file. Please enter the EMD accession number (e.g., EMD-xxxxx)."),
                    ui.input_action_button("select_random_emdb", "Select a random EMDB ID", value=False),
                    ui.input_checkbox("is_3d", "The input is a 3D map", value=True),
            ]
        else:
            return ui.p("Please select an option to proceed.")
    
    # @reactive.event(input_data, input.is_3d)
    # def get_2d_data():
    #     if input.is_3d():
    #         pass
    
    @render.ui
    @reactive.event(input_data, input.is_3d)   
    def input_display_uis():
        return ui.div(
                    ui.output_ui("input_image_selection_gallery"),
                    style="max-height: 500px; overflow-y: auto;"
                )

    @reactive.effect
    def max_map_size_set():
        if util.is_hosted():
            max_map_size_t  = util.mem_quota()/2    # MB
            max_map_dim_t   = int(pow(max_map_size*pow(2, 20)/4, 1./3.)//10*10)    # pixels in any dimension
            stop_map_size_t = util.mem_quota()*0.75 # MB

            max_map_size.set(max_map_size_t)
            max_map_dim.set(max_map_dim_t)
            stop_map_size.set(stop_map_size_t)
        else:
            max_map_size_t = -1   # no limit
            max_map_dim_t  = -1
            max_map_size.set(max_map_size_t)
            max_map_dim.set(max_map_dim_t)
        if max_map_size()>0:
            warning_map_size = f"Due to the resource limit, the maximal map size should be {max_map_dim}x{max_map_dim}x{max_map_dim} voxels or less to avoid crashing the server process"
    
    
    
    
    #@output
    @render.ui
    @reactive.event (input.input_mode_params, input.img_file_emd_id)
    def get_link_emd():
        req(input.input_mode_params() == "3")
        req(input.img_file_emd_id())
        emdb = helicon.dataset.EMDB()
        emd_ids = emdb.helical_structure_ids()
        print(emd_ids)
        return len(emd_ids)

        # emdb_ids_all, methods, resolutions, _, emdb_ids_helical = util.get_emdb_ids()
        # url = "https://www.ebi.ac.uk/emdb/search/*%20AND%20structure_determination_method:%22helical%22?rows=10&sort=release_date%20desc"
        # if emdb_ids_helical is None:
        #     return
        # msg = ""

        # selected_id = str(10499) # str(random.choice(arr))
        # # key_emd_id.set('emd-' + selected_id) # random.choice(emdb_ids_helical))

        # id = input.img_file_emd_id().lower().split("emd-")[-1]
        
        # if id != selected_id:
        #     emd_id.set(id)
        # elif not key_emd_id():
        #     selected_id = str(10499) # str(random.choice(arr))
        #     key_emd_id.set('emd-' + selected_id)
        #     emd_id.set(key_emd_id().lower().split("emd-")[-1])

        #     msg = f'[EMD-{emd_id()}](https://www.ebi.ac.uk/emdb/entry/EMD-{emd_id()})'
        #     emd_url.set(msg)
        #     resolution = resolutions[emdb_ids_all.index(emd_id())]
        #     msg += f' | resolution={resolution}Å'
        #     params = util.get_emdb_helical_parameters(emd_id())

        #     if params and ("twist" in params and "rise" in params and "csym" in params):           
        #         msg += f"  \ntwist={params['twist']}° | rise={params['rise']}Å"
        #         if params.get("csym_known", True): msg += f" | c{params['csym']}"
        #         else: msg += " | csym n/a"
        #         rise.set(params['rise'])
        #         twist.set(params['twist'])
        #         pitch.set(compute.twist2pitch(twist=twist(), rise=rise()))
        #         csym.set(params['csym'])
        #     else:
        #         msg +=  "  \n*helical params not available*"
        #     msg_reactive.set(msg + f'\n[All {len(emdb_ids_helical)} helical structures in EMDB]({url})')
    
    @reactive.effect
    @reactive.event(input.select_random_emdb)
    def select_random_emdb_id():
        req(input.input_mode_params() == "3")
        rand_id = '0013'
        ui.update_text('img_file_emd_id', value='emd-' + rand_id)
    
    @output
    @render.ui
    @reactive.event(msg_reactive)
    def emdb_info_uis():
        req(input.input_mode_params() == "3")
        return ui.markdown(msg_reactive())

    # @reactive.effect
    # @reactive.event(input.input_mode_params, input.img_file_emd_id)
    # def get_map_from_emd_id():
    #     req(input.input_mode_params() == "3")
    #     try:
    #         map_info = compute.MapInfo(
    #             emd_id = input.img_file_emd_id(),
    #         )
    #         input_data.set(map_info.data)
    #         apix_from_file.set(map_info.apix)
    #     except Exception as e:
    #         print("Error in get_map_from_emd_d:", e)

    # @reactive.effect
    # @reactive.event(maps)
    # def get_map_xyz_projections():
    #     try:
    #         if not maps() or len(maps()) == 0:
    #             print("Maps is empty, cannot generate projections")
    #             return
                
    #         map_xyz_projections.set([])
    #         images = []
    #         image_labels = []
    #         apix_val = None
            
    #         map_xyz_projection_title.set(f"Map Y projections:")
    #         # print("projection err 1")
    #         with ui.Progress(min=0, max=len(maps())) as p:
    #             p.set(message="Generating x/yz/ projections", detail="This may take a while ...")
    #             # print("projection err 2")
    #             for mi, m in enumerate(maps()):
    #                 # print("M: : ", m)
    #                 p.set(mi, message=f"{mi+1}/{len(maps())}: x/y/z projecting {m.label}")
    #                 # print("projection err 3")
    #                 try:
    #                     tmp_images, tmp_image_labels, apix = compute.get_one_map_xyz_projects(map_info=m)
    #                     # print("projection err 3.5")
    #                     images += tmp_images
    #                     image_labels += tmp_image_labels

    #                     if apix_val is None:
    #                         apix_val = apix
                        
    #                     # print(f"Generated {len(tmp_images)} projections for map {m.label}")
    #                 except Exception as e:
    #                     print(f"Error generating projections for map {m.label}:", e)
    #             # print("projection err 4")
    #             if images:
    #                 map_xyz_projection_labels.set(image_labels)
    #                 map_xyz_projections.set(images)
    #                 # print(np.shape(images))
    #                 # print(f"Set {len(images)} projections in total")
                
    #             selected_images.set(map_xyz_projections())
                
    #             # Add this part to set apix value for EMD files
    #             if apix_val is not None:
    #                 apix_from_file.set(apix_val)
    #                 apix_auto.set(apix_val)
    #                 # Reset user inputs to force using auto values
    #                 user_inputs.set({
    #                     'apix': None,
    #                     'angle': None,
    #                     'dx': None,
    #                     'dy': None,
    #                     'mask_radius': None
    #                 })
    #                 # Force recalculation of auto parameters for EMD
    #                 recalculate_emd_params.set(True)
    #             # print("projection err 5")
    #     except Exception as e:
    #         print("Error in get_map_xyz_projections:", e)




    @output
    @render.ui
    @reactive.event(input.apply_helical_sym)
    def apply_helical_sym_uis():
        req(input.apply_helical_sym())
        return [
                ui.input_numeric("twist_ahs", "Twist (°):", value=twist(), min=-180.0, max=180.0, step=1.0),
                ui.input_numeric("rise_ahs", "Rise (Å):", value=rise(), min=0.0, step=1.0),
                ui.input_numeric("csym_ahs", "Csym:", value=input.csym(), min=1, step=1),
                ui.input_numeric("apix_map", "Current map pixel size (Å):", value=apix_from_file(), min=0.0, step=1.0),
                ui.input_numeric("apix_ahs", "New map pixel size (Å):", value=apix_from_file(), min=0.0, step=1.0),
                ui.input_numeric("fraction_ahs","Center fraction (0-1):",value=1.0,min=0,max=1.0,step=0.1),
                ui.input_numeric("length_ahs","Box length (Pixel):",value=max(nx(),nz()),min=0,step=1),
                ui.input_numeric("width_ahs", "Box width (Pixel):", value=max(nx(),nz()), min=0, step=1),
                ui.hr(),
        ]

    # @output
    # @render.ui
    # @reactive.event(input.apix_map)
    # def update_apix_map_vals():
    #     if nz()*input.apix_map() == 0:
    #         return
        
    #     return [
    #         ui.input_numeric("apix_ahs", "New map pixel size (Å):", value=input.apix_map(), min=0.0, step=1.0),
    #         ui.input_numeric(
    #             "fraction_ahs",
    #             "Center fraction (0-1):",
    #             value=1.0,
    #             min=input.rise_ahs()/(nz()*input.apix_map()),
    #             max=1.0,
    #             step=0.1,
    #         ),
    #         ui.input_numeric(
    #             "length_ahs",
    #             "Box length (Å):",
    #             value=input.apix_map()*max(nz(),nx()),
    #             min=input.rise_ahs(),
    #             step=1.0,
    #         ),
    #         ui.input_numeric("width_ahs", "Box width (Å):", value=input.apix_map()*nx(), min=0.0, step=1.0),
    #     ]


    @reactive.effect
    @reactive.event(input_data, input.is_3d) # , input.ignore_blank)
    def get_data_all_2d():
        req(input_data() is not None)
        req(not input.is_3d())
        data = input_data()
        n = len(data)
        images = [data[i] for i in range(n)]

        included = []
        included_images = []
        for i in range(n):
            image = images[i]
            if np.max(image) > np.min(image):
                included.append(i)
                included_images.append(image)
        images = included_images
        
        image_labels = [f"{i+1}" for i in included]
        data_all_2d_labels.set(image_labels)
        data_all_2d.set(images)
        
        # Reset selected images when display changes
        selected_images.set([])
        selected_image_labels.set([])
    
    @reactive.effect
    @reactive.event(input_data, input.is_3d) # , input.ignore_blank)
    def get_data_all_2d_from_3d():
        req(input_data() is not None)
        req(not input.is_3d())
        data = input_data()
        n = len(data)
        images = [data[i] for i in range(n)]

        included = []
        included_images = []
        for i in range(n):
            image = images[i]
            if np.max(image) > np.min(image):
                included.append(i)
                included_images.append(image)
        images = included_images
        
        image_labels = [f"{i+1}" for i in included]
        data_all_2d_labels.set(image_labels)
        data_all_2d.set(images)
        
        # Reset selected images when display changes
        selected_images.set([])
        selected_image_labels.set([])

    @reactive.effect
    @reactive.event(input.generate_2d_projection)
    def update_all_images_from_3d_input_data():
        req(input_data() is not None)
        req(input.is_3d)
        data_3d = input_data()
        if input.apply_helical_sym():
            print(f"Applying helical symmetry: twist={input.twist_ahs()}°, rise={input.rise_ahs()}Å, csym={input.csym_ahs()}, apix from {input.apix_map()}Å to {input.apix_ahs()}Å")
            m = compute.symmetrize_transform_map(
                data=data_3d,
                apix=input.apix_map(),
                twist_degree=input.twist_ahs(),
                rise_angstrom=input.rise_ahs(),
                csym=input.csym_ahs(),
                new_size=(input.length_ahs(),input.width_ahs(), input.width_ahs()),
                new_apix=input.apix_ahs(),
                axial_rotation=input.az(),
                tilt=input.tilt(),
            )
            print(np.shape(data_3d),np.shape(m))
            ui.update_numeric('apix', value=input.apix_ahs())
        else:
            print(f"Applying rotation {input.az()}° and tilt {input.tilt()}° to the map")
            m = helicon.transform_map(data_3d, rot=input.az(), tilt=input.tilt())

        proj = np.transpose(m.sum(axis=-1))[:, ::-1]
        proj = proj[np.newaxis, :, :]

        def add_noise(image, noise, thres=1e-3):
            sigma = np.std(image[image > thres])  # ignore background pixels
            image += np.random.normal(scale=sigma * noise, size=image.shape)
            return image

        if input.gauss_noise_std() > 0:
            proj[0, :, :] = add_noise(proj[0, :, :], input.gauss_noise_std())
        proj[0, :, :] = helicon.transform_image(proj[0, :, :], rotation=90.0)

        data_all_2d.set(proj)
        data_all_2d_labels.set(["Projection"])
        selected_images.set(proj)
        selected_image_labels.set(["Projection"])
        data_2d_transformed.set(proj[0])

    #@output
    @render.ui
    def input_image_selection_gallery():
        if len(data_all_2d()) > 0:
            _, ny, nx = np.shape(data_all_2d())
            return helicon.shiny.image_gallery(
                id="select_image",
                label=reactive.value("Select a image:"),
                images=data_all_2d,
                image_labels=data_all_2d_labels,
                image_size=image_display_size,
                initial_selected_indices=initial_selected_image_indices,
                enable_selection=True,
                allow_multiple_selection=False,
            )
        return ui.div("No images available")

    @reactive.effect
    @reactive.event(input.select_image)
    def update_selected_images():
        if input.select_image() is not None:
            #selected = [data_all_2d()[i] for i in input.select_image()]
            selected_images.set(
                [data_all_2d()[i] for i in input.select_image()]
            )
            selected_image_labels.set(
                [data_all_2d_labels()[i] for i in input.select_image()]
            )

            mode = input.input_type() 

            if mode in ["PS", "PD"]:
                # angle_auto.set(0.0)
                # dx_auto.set(0.0)
                # dy_auto.set(0.0)
                ui.update_numeric('angle', value=0.0)
                ui.update_numeric('dx', value=0.0)
                ui.update_numeric('dy', value=0.0)
            else:
                print("Estimating auto transformation parameters...")
                ang_auto_v, dy_auto_v, diameter = helicon.estimate_helix_rotation_center_diameter(selected_images()[0], estimate_center=True, estimate_rotation=True)
                # angle_auto.set(round(ang_auto_v+90,2))
                # dy_auto.set(round(dy_auto_v*input.apix(),2))
                # mask_radius_auto.set(round(diameter/2*input.apix(),2))
                print("Estimated angle:", ang_auto_v+90)
                print("input angle before:", input.angle())
                ui.update_numeric('angle', value=round(ang_auto_v+90,2))
                print("input angle after:", input.angle())
                ui.update_numeric('dy', value=round(dy_auto_v*input.apix(),2))
                ui.update_numeric('mask_radius', value=round(diameter/2*input.apix(),2))
                mask_len_percent_auto_default = 90.0
                if input.straightening():
                    mask_len_percent_auto_default = 100.0
                #mask_len_percent_auto.set(mask_len_percent_auto_default)
                ui.update_numeric('mask_len', value=mask_len_percent_auto_default)
            print(selected_images())
            print("Auto transformation parameters set: angle=", input.angle())            
            
    @reactive.effect
    @reactive.event(input.img_file_upload)
    def get_images_from_input_2d_upload():
        if input.input_mode_params() == "1":
            print("Uploading image file...")
            fileinfo = input.img_file_upload()
            if fileinfo is not None:
                class_file = fileinfo[0]["datapath"]
                try:
                    data, apix, crs = compute.get_class2d_from_file(class_file)
                    print(crs)
                    nz_file, ny_file, nx_file = data.shape
                    apix_from_file.set(apix)
                    update_with_apix_from_file(apix)
                    #update_dimensions(data)
                    input_data.set(data)
                    nx.set(nx_file)
                    ny.set(ny_file)
                    nz.set(nz_file)
                    update_with_ny_nx(ny_file, nx_file)
                except Exception as e:
                    print(e)
                    modal = ui.modal(
                        f"Failed to read the uploaded images from {fileinfo[0]['name']}",
                        title="File upload error",
                        easy_close=True,
                        footer=None
                    )
                    ui.modal_show(modal)
    
    
    @reactive.effect
    @reactive.event(input.img_file_url)
    def get_images_from_input_2d_url():
        if input.input_mode_params() == "2":
            print("Downloading image file from URL...")
            url = input.img_file_url()
            try:
                data, apix, crs = compute.get_class2d_from_url(url)
                print(crs)
                nz_file, ny_file, nx_file = data.shape
                apix_from_file.set(apix)
                update_with_apix_from_file(apix)
                #update_dimensions(data)
                input_data.set(data)
                nx.set(nx_file)
                ny.set(ny_file)
                nz.set(nz_file)
                update_with_ny_nx(ny_file, nx_file)
            except Exception as e:
                print(e)
                modal = ui.modal(
                    f"Failed to download images from {url}",
                    title="URL download error",
                    easy_close=True,
                    footer=None
                )
                ui.modal_show(modal)

    @reactive.effect
    @reactive.event(input.img_file_emd_id)
    def get_images_from_input_emd():
        return

    def update_with_apix_from_file(apix):
        ui.update_numeric('apix', value=apix)
        ui.update_numeric('res_limit_x', value=round(3*apix_from_file(),4), min=2*apix_from_file())
        ui.update_numeric('res_limit_y', value=round(2*apix_from_file(),4), min=2*apix_from_file())

    # @reactive.Effect
    # @reactive.event(input.is_3d)
    # def get_data_2d():
    #     if input.is_3d():
    #         pass
    #     else:
    #         data_2d = selected_images()[0]
    
    # functions for outputting graphs
    @render_plotly
    @reactive.event(selected_images, input.apix)
    def display_selected_image():
        req(len(selected_images())>0)
        print("update display_selected_image:"+str(input.apix()))

        h, w = selected_images()[0].shape[:2]
        # nx.set(h)
        # ny.set(w)
        # update_with_ny_nx(w, h)

        img = selected_images()[0]

        images = img
        print("origin image to display: ", images.shape)
        print("min image to display: ", np.min(images))
        print("max image to display: ", np.max(images))

        # Plot the micrograph
        fig = compute.plot_2d_class(
            micrograph=images,
            title=f"Original image ({nx()}x{ny()})",
            apix=input.apix(),
        )

        return fig
    
    def update_with_ny_nx(ny_val, nx_val):
        ui.update_numeric('pnx', min=min(ny_val, 128))
        ui.update_numeric('pny', min=min(nx_val, 512))
    
    @render_plotly
    @reactive.event(data_2d_transformed)
    def display_transformed_data():
        req(data_2d_transformed() is not None)
        images = data_2d_transformed() # selected_images()
        if data_2d_transformed() is None:
            return None # px.scatter(title="No image selected")

        h, w = data_2d_transformed().shape[:2]
        #nx.set(h)
        #ny.set(w)

        images = data_2d_transformed()

        image_to_display = (
            0-images if input.negate() else images #selected_images()[0]
        )

        print("transformed image to display: ", image_to_display.shape)
        print("average image to display: ", np.average(image_to_display))
        print("apix_from_file:", apix_from_file())

        fig = compute.plot_2d_class(
            micrograph= image_to_display, # image_to_display, #255 - selected_images()[0],
            title=f"Transformed image ({nx()}x{ny()})",
            apix= input.apix(),
        )

        return fig
    
    @render_plotly
    @reactive.event(data_2d_transformed)
    def plot_radial_profile():
        req(data_2d_transformed() is not None)
        
        data = data_2d_transformed()
        mask_radius = input.mask_radius()
        
        x = np.arange(-nx() // 2, nx() // 2) * input.apix()
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
            #height=400,
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
        if len(images)<=0:
            #print("No selected images!")
            return None, None, None

        data = images[0]

        acf = auto_correlation(data, sqrt=True, high_pass_fraction=0.1)

        if acf is None or np.isnan(acf).any():
            print("Invalid auto-correlation data!")
            return None, None, None

        ny = acf.shape[0]
        y = np.arange(-ny // 2, ny // 2) * input.apix() 
        xmax = np.max(acf, axis=1)  

        return xmax, y, acf


    @render_plotly
    @reactive.event(data_2d_transformed)
    def acf_plot():
        req(data_2d_transformed() is not None)
        xmax, y, acf = acf_data()

        if xmax is None or y is None:
            return None # px.scatter(title="No valid data available")

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=xmax, y=y, mode='lines', name="ACF", line=dict(color='red')
        ))

        fig.update_layout(
            xaxis_title="Auto-correlation",
            yaxis_title="Axial Shift (Å)",
            showlegend=True,
            #height=400,
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

