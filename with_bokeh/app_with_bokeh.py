import numpy as np
import pandas as pd

from shiny import App, reactive, render, ui

from shinywidgets import output_widget, render_widget, render_plotly, bokeh_dependency
import pandas as pd
from shiny import reactive, req
import compute
import util
import plotly.graph_objects as go

from bokeh.plotting import figure
from bokeh.layouts import gridplot, column, layout
from bokeh.events import MouseMove, MouseEnter, DoubleTap
from bokeh.models import Button, ColumnDataSource, CustomJS, Legend, Span, Spinner, Slider
from bokeh.models.tools import CrosshairTool, HoverTool

from jupyter_bokeh import BokehModel

import pandas as pd

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
nx_curr_img = reactive.value(256)  # default value
ny_curr_img = reactive.value(256)  # default value
nz_curr_img = reactive.value(256)    # default value
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

markers_straighten_x = reactive.value([])
markers_straighten_y = reactive.value([])

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
    bokeh_dependency(),
    ui.tags.style("""       
        .main-page-container {
            max-height: 100vh; 
            overflow: auto;  
            //display: flex;
            //flex-direction: row;
            //white-space: nowrap;  /* Prevent wrapping */
        }
        .scrollable-sidebar {
            //height: auto;
            max-height: 150vh; /* Adjust height as needed */
            overflow-y: auto; /* Enables vertical scrolling */
            border: 1px solid #ccc; /* Optional: for visual clarity */
            padding: 10px;
        }
        .scrollable-container {
            max-height: 100%; 
            display: inline-block;
            //display: inline-flex; /* this or inline-block */
            //flex-direction: column;
            white-space: nowrap;  /* Prevent wrapping */
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
        .footer {
            // https://stackoverflow.com/questions/67763901/footer-position-in-shiny
            position: absolute;
            bottom: 0;
        }
    """),
    ui.layout_sidebar(
        ui.sidebar(
            ui.navset_pill(
                ui.nav_panel("Inputs",
                    ui.div(
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
                                    inline=True
                                ),
                                value="input_mode"
                            )
                        ),
                        ui.output_ui("conditional_input_uis"),
                        ui.output_ui("conditional_3d_transformation_uis"),
                        #ui.input_action_button("run", label="Run"),
                        ui.output_ui("input_display_uis"),
                        output_widget("display_selected_image",height="auto"), # shinywidgets.css has flex: 1 1 400px which makes a gap after the widget
                        ui.accordion(
                            ui.accordion_panel(
                                ui.p("Image Parameters"),
                                ui.input_radio_buttons("input_type", "Input is:", choices=["Image", "PS", "PD"], inline=True),
                                ui.output_ui("img_2d_transformation_uis"),
                                value="image_parameters_emd"
                            )
                        ),
                        output_widget("display_transformed_data",height="auto"),
                        output_widget("plot_radial_profile",height="auto"),
                        output_widget("acf_plot",height="auto"),
                    )
                ),
                ui.nav_panel("Parameters",
                    ui.layout_columns(
                        ui.input_radio_buttons("use_twist_pitch", "Use twist/pitch:", choices=["Twist", "Pitch"], selected="Twist", inline=True),
                        ui.input_numeric('csym', 'csym', value=6, min=1, step=1, update_on="blur"),
                        ui.input_numeric('filament_diameter', 'Filament/tube diameter (Å)', value=69.0*2, min=1.0, max=1000.0, step=10.0, update_on="blur"),
                        ui.input_numeric("out_of_plane_tilt", "Out-of-plane tilt (°)", value=0.0, min=-90.0, max=90.0, step=1.0, update_on="blur"),
                        ui.input_numeric('res_limit_x', 'Resolution limit - X (Å)', value=round(3*1.0,4), min=2*1.0, step=1.0, update_on="blur"),
                        ui.input_numeric("res_limit_y", "Resolution limit - Y (Å)", value=round(2*1.0,4), min=2*1.0, step=1.0, update_on="blur"),
                        ui.input_checkbox("fft_top_only", "Only display the top half of FFT", value=False),

                        ui.input_checkbox("log_amp", "Log(amplitude)", value=True),
                        ui.input_text("const_image_color", "Flatten the PS/PD image in this color", value="", placeholder="white"),
                        ui.input_text("ll_colors", 'Layerline colors', value="lime cyan violet salmon silver"),
                        ui.input_numeric("hp_fraction", 'Fourier high-pass (%)', value=0.40, min=0.0, max=100.0, step=0.1, update_on="blur"),
                        ui.input_numeric("lp_fraction", 'Fourier low-pass (%)', value=0.00, min=0.0, max=100.0, step=0.1, update_on="blur"),
                        ui.input_numeric("pnx", 'FFT X-dim size (pixels)', value=512, min=128, step=2, update_on="blur"),
                        ui.input_numeric("pny", 'FFT Y-dim size (pixels)', value=1024, min=512, step=2, update_on="blur"),
                        col_widths=6,
                        style="align-items: flex-end;"
                    ),

                    #ui.input_checkbox("share_url", "Show/Reload sharable URL", value=False),
                    ui.div(
                        #ui.input_action_button("copy_twist", "Copy twist/rise◀"),
                        ui.h5("Display:"),
                        ui.layout_columns(
                            ui.input_checkbox("PS", "PS", value=True),
                            ui.input_checkbox("YP", "YP", value=True),
                            ui.input_checkbox("Phase", "Phase", value=False),
                            ui.input_checkbox("PD", "PD", value=True),
                            ui.input_checkbox("Color", "Color", value=True),
                            ui.input_checkbox("LL", "LL", value=True),
                            ui.input_checkbox("LLText", "LLText", value=True),
                            col_widths=6,
                            style="align-items: flex-end;"
                        ),
                        #ui.h5("m:"),
                        ui.input_numeric("m_max", "Max=", value=3, min=1, step=1, update_on="blur"),
                        ui.output_ui("show_m_choices"),
                        style="display: inline; text-align: left;",
                        class_="wrap-text"
                    ),
                ),
                ui.nav_menu("Others",
                    ui.nav_panel("Filament Straightening",
                        ui.div(
                            ui.input_numeric("num_markers_straighten", "Number of markers", value=2, min=2, step=1, update_on="blur"),
                            ui.input_numeric("template_diameter_straighten", "Template diameter (Å)", value=50, min=1, step=1, update_on="blur"),
                            ui.input_numeric("template_length_straighten", "Template length (Å)", value=100, min=1, step=1, update_on="blur"),
                            ui.input_numeric("angular_step_straighten", "Angular step (°)", value=1.0, min=0.1, step=0.1, update_on="blur"),
                            ui.input_numeric("lp_angst_straighten", "Low-pass filter (Å)", value=10.0, min=2.0, step=1.0, update_on="blur"),
                            ui.input_numeric("hp_angst_straighten", "High-pass filter (Å)", value=200.0, min=2.0, step=1.0, update_on="blur"),
                            ui.input_numeric("output_width_straighten", "Output width (pixels)", value=256, min=16, step=1, update_on="blur"),
                            ui.input_numeric("output_height_straighten", "Output height (pixels)", value=256, min=16, step=1, update_on="blur"),
                            style="text-align: left; border: 1px solid #ddd; padding: 10px;",
                            class_="wrap-text"
                        ),                    
                    ),
                    ui.nav_panel("Simulation",
                        ui.input_numeric("ball_radius_sim", 'Gaussian radius (Å)', value=0.0, min=0.0, max=100.0, step=5.0, update_on="blur"),
                    ),
                ),
                id="sidebar_navset"
            ),
            class_="scrollable-sidebar",
            width = "20vw"
        ),
        ui.navset_hidden(
            ui.nav_panel(None,
                ui.div(
                    ui.row(
                        ui.column(12,
                            ui.div(
                                ui.h2("HILL: Helical Indexing using Layer Lines"),
                                style="text-align: center; margin-bottom: 20px;"
                            )
                        )
                    ),
                    ui.row(
                        ui.column(12,
                            ui.div(
                                output_widget("main_plots", height="auto"),
                                style="text-align: center; margin-bottom: 0px;"
                            )
                        )
                    ),
                    ui.input_numeric("rise", "Rise (Å)", value=20.0, min=2.0, step=0.1, update_on="blur"),
                    ui.row(
                        ui.column(12,
                            ui.div(
                                ui.markdown("*Developed by the [Jiang Lab@Penn State](https://jianglab.science.psu.edu). Report problems to [HILL@GitHub](https://github.com/jianglab/hill/issues)*"),
                                class_="footer"
                            )
                        )
                    ),
                    #style="overflow-x: auto; max-width: 100%"
                    class_="scrollable-container"
                ),
                value="main_content_tab",
            ),
            ui.nav_panel(None,
                ui.output_ui("filament_straightening_uis"),
                value="straighten_filament_tab"
            ),
            id="main_tabs"
        ),
        helicon.shiny.setup_ajdustable_sidebar(width="20vw"),
    ),
    class_="main-page-container",
)


def server(input, output, session):
    
    @reactive.effect
    @reactive.event(input.rise)
    def test_input_rise_update():
        print("input.rise updated to:", input.rise())
    @reactive.effect
    @reactive.event(rise)
    def test_rise_update():
        print("rise updated to:", rise())

    @reactive.effect
    @reactive.event(input.sidebar_navset)
    def main_page_tab_switch():
        print(input.sidebar_navset())
        if input.sidebar_navset() == "Filament Straightening":
            ui.update_navs("main_tabs", selected="straighten_filament_tab")
            print(input.apix())
        else:
            ui.update_navs("main_tabs", selected="main_content_tab")

    @output
    @render.text
    def m_max_auto():
        return int(np.floor(np.abs(rise()/input.res_limit_y()))) + 3

    @reactive.calc
    def m_groups():
        req(input.filament_diameter()>0)
        req(input.res_limit_y()>0)
        print("Updating mgroups")
        return compute.compute_layer_line_positions(
            twist=twist(),
            rise=rise(),
            csym=input.csym(),
            radius=input.filament_diameter()/2,
            tilt=input.out_of_plane_tilt(),
            cutoff_res=input.res_limit_y(),
            m_max=input.m_max()
        )
    
    @reactive.effect
    @reactive.event(input.apix)
    def set_min_rise():
        min_rise.set(round(2*input.apix(),2))

    @output(suspend_when_hidden=False)
    @render.ui
    def show_m_choices():
        m_max = input.m_max()
        ms = dict()
        for m in range(m_max, -m_max-1, -1):
            ms[m] = str(m)
        checkboxes = ui.input_checkbox_group("ms","m:",ms, selected=[-1, 0, 1])
        return checkboxes
    
    striaghten_fig = reactive.value(go.FigureWidget())
    
    @render_plotly
    @reactive.event(striaghten_fig)
    def figure_get_markers_straighten():
        return striaghten_fig()


    @output
    @render.ui
    def filament_straightening_uis():
        return ui.div(
            output_widget("figure_get_markers_straighten", width="auto"),
            ui.input_action_button("reset_markers", "Clear markers"),
            ui.input_action_button("straighten_run", "Straighten"),
            style="text-align: center; border: 1px solid #ddd; padding: 10px;",
            class_="wrap-text"
        )
        #return "Temp Text"

    @reactive.effect
    @reactive.event(input.reset_markers)
    def clear_markers():
        striaghten_fig().layout.shapes=[]
        markers_straighten_x.set([])
        markers_straighten_y.set([])
    
    @reactive.effect
    @reactive.event(input.straighten_run)
    def run_straightening():
        print(markers_straighten_x(), markers_straighten_y())
        sorted_indices = np.argsort(markers_straighten_y())
        xs = np.array(markers_straighten_x())[sorted_indices]/input.apix()
        ys = np.array(markers_straighten_y())[sorted_indices]/input.apix()
        print(xs,ys)
        new_xs, tck = compute.fit_spline(striaghten_fig(), data_2d_transformed(), xs, ys, input.apix(), display=False)
        data_2d_transformed.set(np.nan_to_num(compute.filament_straighten(None,data_2d_transformed(), tck, new_xs, ys, input.output_width_straighten()//2,input.apix()), 0.0))
    
    @render_plotly
    @reactive.event(data_2d_transformed, input.straightening)
    def update_figure_get_markers_straighten():
        req(input.straightening())
        req(data_2d_transformed() is not None)
        img = data_2d_transformed()
        h, w = img.shape
        apix = input.apix()
        marker_r_angst = input.template_diameter_straighten()/2.0 * apix

        fig = striaghten_fig()
        if len(fig.data) > 0:
            print("update straighten plot")
            fig.update_traces(z=data_2d_transformed(),selector=dict(name="image_to_straighten"))
        else:
            print("generating straighten plot")
            fig.add_trace(
                go.Heatmap(
                    name="image_to_straighten",
                    z=img,
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
            fig.add_trace(
                go.Line(
                    name="fitted_spline",
                    x=[],
                    y=[],
                ),
            )
            fig.update_layout(width=w, height=h)
            fig.update_xaxes(showticklabels=False)
            fig.update_yaxes(showticklabels=False, autorange="reversed")

            def update_markers(trace,points,state):
                click_x = points.xs[0]
                click_y = points.ys[0]
                markers_straighten_x.set(markers_straighten_x()+[click_x])
                markers_straighten_y.set(markers_straighten_y()+[click_y])
                fig.layout.shapes = list(fig.layout.shapes) + [dict(type='circle',x0=click_x-marker_r_angst, y0=click_y-marker_r_angst, x1=click_x+marker_r_angst, y1=click_y+marker_r_angst, line_color='red', name='markers')]

            
            fig.data[0].on_click(update_markers)
            markers_straighten_x.set([])
            markers_straighten_y.set([])
        return fig

    main_plot = reactive.value(None)
    
    @render_widget
    @reactive.event(main_plot)
    def main_plots():
        req(main_plot() is not None)
        return main_plot()

    #@render_widget
    @reactive.effect
    @reactive.event(ps_data)
    def update_main_plots():
        req(input.res_limit_y()>0)
        req(input.res_limit_x()>0)
        req(ps_data() is not None)
        print('entering main_plots')
        # print("input.PS():", input.PS())
        # print("input.max():", input.m_max())
        # print("input.ms():", input.ms())

        figs = []
        figs_image = []
        figs_width = 0
        pwr_work = ps_data()
        phase_work = pd_data()
        phase_diff_work = pd_data()
        if input.PS():
            tooltips = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Amp', '@image')]
            fig = compute.create_layerline_image_figure(
                pwr_work, 
                input.res_limit_x(), 
                input.res_limit_y(), 
                input.filament_diameter()/2.0, 
                input.out_of_plane_tilt(), 
                phase=phase_work if input.Phase() else None, #TODO: add phase data
                fft_top_only=input.fft_top_only(), 
                pseudo_color=input.Color(), 
                const_image_color=input.const_image_color(), 
                title="Power Spectra", 
                yaxis_visible=False, 
                tooltips=tooltips
            )
            figs.append(fig)
            figs_image.append(fig)
            figs_width += pwr_work.shape[-1]

        if input.YP():
            ny, nx = pwr_work.shape
            dsy = 1/(ny//2*input.res_limit_y())
            y=np.arange(-ny//2, ny//2)*dsy
            yinv = y*1.0
            yinv[yinv==0] = 1e-10
            yinv = 1/np.abs(yinv)
            yprofile = np.mean(pwr_work, axis=1)
            yprofile /= yprofile.max()
            source_data = ColumnDataSource(data=dict(yprofile=yprofile, y=y, resy=yinv))
            tools = 'box_zoom,hover,pan,reset,save,wheel_zoom'
            tooltips = [('Res y', '@resy Å'), ('Amp', '$x')]
            fig = figure(frame_width=nx//2, frame_height=figs[-1].frame_height, y_range=figs[-1].y_range, y_axis_location = "right", title=None, tools=tools, tooltips=tooltips)
            fig.line(source=source_data, x='yprofile', y='y', line_width=2, color='blue')
            fig.yaxis.visible = False
            fig.hover[0].attachment="vertical"
            figs.append(fig)
            figs_width += nx//2

        if input.PD():
            tooltips = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Phase Diff', '@image °')]
            fig = compute.create_layerline_image_figure(
                phase_diff_work, 
                input.res_limit_x(), 
                input.res_limit_y(), 
                input.filament_diameter()/2, 
                input.out_of_plane_tilt(), 
                phase=phase_work if input.Phase() else None, 
                fft_top_only=input.fft_top_only(), 
                pseudo_color=True, 
                const_image_color=input.const_image_color(), 
                title="Phase Diff Across Meridian", 
                yaxis_visible=False, 
                tooltips=tooltips
            )
            figs.append(fig)
            figs_image.append(fig)
            figs_width += phase_diff_work.shape[-1]   
        
        figs[-1].yaxis.fixed_location = figs[-1].x_range.end
        figs[-1].yaxis.visible = True
        crosshair_width = Span(dimension="width", line_color="red")
        crosshair_height = Span(dimension="height", line_color="red")
        add_linked_crosshair_tool(figs, overlay=(crosshair_width))
        add_linked_crosshair_tool(figs_image, overlay=(crosshair_width, crosshair_height))

        fig_ellipses = []
        curr_m_groups = m_groups()
        curr_ll_colors=input.ll_colors()
        if figs_image and input.LL():
            if max(curr_m_groups[0]["LL"][0])>0:
                x, y, n = curr_m_groups[0]["LL"]
                tmp_x = np.sort(np.unique(x))
                width = np.mean(tmp_x[1:]-tmp_x[:-1])
                height = width/5
                for mi, m in enumerate(curr_m_groups.keys()):
                    if str(m) not in input.ms(): continue
                    x, y, bessel_order = curr_m_groups[m]["LL"]
                    if input.LLText():
                        texts = [str(int(n)) for n in bessel_order]
                    tags = [m, bessel_order]
                    color = curr_ll_colors[abs(m)%len(curr_ll_colors)]
                    #bessel_colors = ["cyan","greenyellow"]
                    ellipse_alpha = [n%2*1.0 for n in bessel_order]
                    for f in figs_image:
                        if input.LLText():
                            text_labels = f.text(x, y, y_offset=2.0, text=texts, text_color=color, text_baseline="middle", text_align="center")
                            text_labels.tags = tags
                            text_labels.level = "overlay"
                            print("text labels:", text_labels)
                            fig_ellipses.append(text_labels)
                        else:
                            ellipses = f.ellipse(x, y, width=width, height=height, line_color=color, fill_color=color, fill_alpha=ellipse_alpha, line_width=1.0)
                            ellipses.tags = tags
                            fig_ellipses.append(ellipses)
            # else:
            #     st.warning(f"No off-equator layer lines to draw for Pitch={pitch:.2f} Csym={csym} combinations. Consider increasing Pitch or reducing Csym")
        
        #from bokeh.models import CustomJS
        #from bokeh.events import MouseEnter
        title_js = CustomJS(args=dict(title="HILL: Helical Indexing using Layer Lines"), code="""
            document.title=title
        """)
        figs[0].js_on_event(MouseEnter, title_js)

        if fig_ellipses:
            slider_width = figs_width//3 if len(figs)>1 else figs_width
            #from bokeh.models import Slider, CustomJS
            spinner_twist = Spinner(title='Twist (°)', low=-180.0, high=180.0, step=1.0, value=twist(), format="0.00", width=slider_width)
            spinner_pitch = Spinner(title='Pitch (Å)', low=min_pitch(), step=1.0, value=max(min_pitch(), pitch()), format="0.00", width=slider_width)
            spinner_rise = Spinner(title='Rise (Å)', low=min_rise(), high=max_rise(), step=1.0, value=rise(), format="0.00", width=slider_width)

            slider_twist = Slider(start=-180, end=180, value=twist(), step=0.01, title="Twist (°)", width=slider_width)
            slider_pitch = Slider(start=pitch()/2, end=pitch()*2.0, value=pitch(), step=pitch()*0.002, title="Pitch (Å)", width=slider_width)
            slider_rise = Slider(start=rise()/2, end=min(max_rise(), rise()*2.0), value=rise(), step=min(max_rise(), rise()*2.0)*0.001, title="Rise (Å)", width=slider_width)
            
            callback_rise_code = """
                slider_rise.value = spinner_rise.value
                Shiny.setInputValue("rise", slider_rise.value, {priority: 'event'})
            """
            callback_pitch_code = """
                slider_pitch.value = spinner_pitch.value
            """
            callback_twist_code = """
                slider_twist.value = spinner_twist.value
            """
            callback_twist = CustomJS(args=dict(spinner_twist=spinner_twist, spinner_pitch=spinner_pitch, spinner_rise=spinner_rise,slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise), code=callback_twist_code)
            callback_pitch = CustomJS(args=dict(spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise,slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise), code=callback_pitch_code)
            callback_rise = CustomJS(args=dict(spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise,slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise), code=callback_rise_code)

            spinner_twist.js_on_change('value', callback_twist)
            spinner_pitch.js_on_change('value', callback_pitch)
            spinner_rise.js_on_change('value', callback_rise)

            callback_rise_code = """
                var twist_sign = 1.
                if (slider_twist.value < 0) {
                    twist_sign = -1.
                }
                var slider_twist_to_update = twist_sign * 360/(slider_pitch.value/slider_rise.value)
                if (slider_twist_to_update != slider_twist.value) {
                    slider_twist.value = slider_twist_to_update
                }
                if (spinner_rise.value != slider_rise.value) {
                    spinner_rise.value = slider_rise.value
                }                  
                var pitch_inv = 1./slider_pitch.value
                var rise_inv = 1./slider_rise.value
                for (var fi = 0; fi < fig_ellipses.length; fi++) {
                    var ellipses = fig_ellipses[fi]
                    const m = ellipses.tags[0]
                    const ns = ellipses.tags[1]
                    var y = ellipses.data_source.data.y
                    for (var i = 0; i < ns.length; i++) {
                        const n = ns[i]
                        y[i] = m * rise_inv + n * pitch_inv
                    }
                    ellipses.data_source.change.emit()
                }
            """
            callback_pitch_code = """
                var twist_sign = 1.
                if (slider_twist.value < 0) {
                    twist_sign = -1.
                }
                var slider_twist_to_update = twist_sign * 360/(slider_pitch.value/slider_rise.value)
                if (slider_twist_to_update != slider_twist.value) {
                    slider_twist.value = slider_twist_to_update
                }
                if (spinner_pitch.value != slider_pitch.value) {
                    spinner_pitch.value = slider_pitch.value
                } 
                var pitch_inv = 1./slider_pitch.value
                var rise_inv = 1./slider_rise.value
                for (var fi = 0; fi < fig_ellipses.length; fi++) {
                    var ellipses = fig_ellipses[fi]
                    const m = ellipses.tags[0]
                    const ns = ellipses.tags[1]
                    var y = ellipses.data_source.data.y
                    for (var i = 0; i < ns.length; i++) {
                        const n = ns[i]
                        y[i] = m * rise_inv + n * pitch_inv
                    }
                    ellipses.data_source.change.emit()
                }
            """
            callback_twist_code = """
                var slider_pitch_to_update = Math.abs(360/slider_twist.value * slider_rise.value)                   
                if (slider_pitch_to_update != slider_pitch.value) {
                    slider_pitch.value = slider_pitch_to_update
                }
                if (spinner_twist.value != slider_twist.value) {
                    spinner_twist.value = slider_twist.value
                } 
            """
            callback_rise = CustomJS(args=dict(fig_ellipses=fig_ellipses, slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_rise_code)
            callback_pitch = CustomJS(args=dict(fig_ellipses=fig_ellipses, slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_pitch_code)
            callback_twist = CustomJS(args=dict(slider_twist=slider_twist, slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_twist_code)
            slider_twist.js_on_change('value', callback_twist)
            slider_pitch.js_on_change('value', callback_pitch)
            slider_rise.js_on_change('value', callback_rise)

            callback_code = """
                let url = new URL(document.location)
                let params = url.searchParams
                params.set("twist", Math.round(slider_twist.value*100.)/100.)
                params.set("rise", Math.round(slider_rise.value*100.)/100.)
                //document.location = url.href
                history.replaceState({}, document.title, url.href)
                if (reload) {
                    var class_names = ["css-1x8cf1d edgvbvh10"]
                    // <button kind="secondary" class="css-1x8cf1d edgvbvh10">
                    console.log(class_names)
                    var i
                    for (i=0; i<class_names.length; i++) {
                        console.log(i, class_names[i])
                        let reload_buttons = document.getElementsByClassName(class_names[i])
                        console.log(reload_buttons)
                        if (reload_buttons.length>0) {
                            reload_buttons[reload_buttons.length-1].click()
                            break
                        }
                    }
                }
            """
            reload = input.input_mode_params() in [1, 2, 3] #and input_mode2 in [None, 1, 2, 3]
            callback = CustomJS(args=dict(slider_twist=slider_twist, slider_pitch=slider_pitch, slider_rise=slider_rise, reload=reload), code=callback_code)
            slider_twist.js_on_change('value_throttled', callback)
            slider_pitch.js_on_change('value_throttled', callback)
            slider_rise.js_on_change('value_throttled', callback)

            callback_code = """
                let url = new URL(document.location)
                let params = url.searchParams
                params.set("twist", Math.round(spinner_twist.value*100.)/100.)
                params.set("rise", Math.round(spinner_rise.value*100.)/100.)
                //document.location = url.href
                history.replaceState({}, document.title, url.href)
                if (reload) {
                    var class_names = ["css-1x8cf1d edgvbvh10"]
                    // <button kind="secondary" class="css-1x8cf1d edgvbvh10">
                    console.log(class_names)
                    var i
                    for (i=0; i<class_names.length; i++) {
                        console.log(i, class_names[i])
                        let reload_buttons = document.getElementsByClassName(class_names[i])
                        console.log(reload_buttons)
                        if (reload_buttons.length>0) {
                            reload_buttons[reload_buttons.length-1].click()
                            break
                        }
                    }
                }
            """
            callback = CustomJS(args=dict(spinner_twist=spinner_twist, spinner_pitch=spinner_pitch, spinner_rise=spinner_rise, reload=reload), code=callback_code)                
            spinner_twist.js_on_change('value_throttled', callback)
            spinner_pitch.js_on_change('value_throttled', callback)
            spinner_rise.js_on_change('value_throttled', callback)
                            
            #spinner_rise.on_change('value_throttled', test_callback)
            
            #callback_code = """
            #    document.dispatchEvent(
            #        new CustomEvent("RiseUpdateEvent", {detail: {rise: spinner_rise.value}})
            #    )
            #"""
            #button_update_param = Button(label="Save parameters", button_type="success")
            #callback_button_update = CustomJS(args=dict(spinner_rise=spinner_rise), code=callback_code)
            #button_update_param.js_on_event('button_click', callback_button_update)

            if len(figs)==1:
                #from bokeh.layouts import column
                figs[0].toolbar_location="right"
                figs_grid = column(children=[[spinner_twist, spinner_pitch, spinner_rise],[slider_twist, slider_pitch, slider_rise], figs[0]])
                override_height = input.pny()+180
            else:
                #from bokeh.layouts import layout
                figs_row = gridplot(children=[figs], toolbar_location='right')
                figs_grid = layout(children=[[spinner_twist, spinner_pitch, spinner_rise],[slider_twist, slider_pitch, slider_rise], figs_row])
                override_height = input.pny()+120
        else:
            figs_grid = gridplot(children=[figs], toolbar_location='right')
            override_height = input.pny()+120
        
        print('end of main_plots')
        print(figs_grid)
        main_plot.set(figs_grid)
        #return figs_grid

    def add_linked_crosshair_tool(figures, overlay=()):
        # https://docs.bokeh.org/en/latest/docs/user_guide/interaction/linking.html
        # create a linked crosshair tool among the figures
        crosshair = CrosshairTool(overlay=overlay)
        for fig in figures:
            fig.add_tools(crosshair)

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
                        ui.input_numeric("az", "Rotation around the helical axis (°):", min=0.0, max=360.0, value=0.0, step=1.0, update_on="blur"),
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
        if value in ['Image']:
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
                    min=-nx_curr_img()*apix_from_file(), 
                    max=nx_curr_img()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "dy", 
                    "Shift along Y-dim (Å)", 
                    value=dy_auto(),
                    min=-ny_curr_img()*apix_from_file(), 
                    max=ny_curr_img()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
            ]

    @render.ui
    @reactive.event(input.is_3d, input.input_type)
    def add_transpose():
        transpose_auto = (input.input_type() in ["Image"]) and (nx_curr_img() > ny_curr_img())
        if not input.is_3d():
            return ui.input_checkbox("transpose", "Transpose the image", value=transpose_auto),

    @reactive.effect
    @reactive.event(selected_images, input.apix, input.angle, input.dx, input.dy, input.mask_radius, input.negate, input.mask_len)
    def set_data_2d_transformed():
        req(input.input_type() in ["Image"])
        req(len(selected_images())>0)
        req(input.apix()!=0)
        if (input.angle() or input.dx() or input.dy() or input.negate() or input.mask_radius()):
            data_2d_transformed.set(compute.mask_2d_filament(compute.transform_2d_filament(selected_images()[0],input.angle(),input.dx(),input.dy(),input.negate(),input.apix()),input.mask_radius(), input.apix(), input.mask_len()/100.0))
    
    @reactive.effect
    @reactive.event(data_2d_transformed)
    def set_mask_radius_auto():
        req(data_2d_transformed() is not None) 
        input_type = input.input_type()
        if input_type in ["Image"]:
            radius_auto_v, _ = compute.estimate_radial_range(data_2d_transformed(), thresh_ratio=0.1)
            #radius_auto.set(round(radius_auto_v*input.apix(),2))
            ui.update_numeric('filament_diameter', value=round(radius_auto_v*input.apix()*2,2))
    
    @reactive.effect
    def get_ps_pd_data():
        req(data_2d_transformed() is not None)
        req(input.apix()>0)
        req(input.res_limit_y()>0)
        req(input.res_limit_x()>0)
        pwr, phase = compute.compute_power_spectra(data_2d_transformed(), apix=input.apix(), cutoff_res=(input.res_limit_y(), input.res_limit_x()), 
                output_size=(input.pny(), input.pnx()), log=input.log_amp(), low_pass_fraction=input.lp_fraction()/100.0, high_pass_fraction=input.hp_fraction()/100.0)
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
                ui.input_numeric("length_ahs","Box length (Pixel):",value=max(nx_curr_img(),nz_curr_img()),min=0,step=1),
                ui.input_numeric("width_ahs", "Box width (Pixel):", value=max(nx_curr_img(),nz_curr_img()), min=0, step=1),
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
                    nx_curr_img.set(nx_file)
                    ny_curr_img.set(ny_file)
                    nz_curr_img.set(nz_file)
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
                nx_curr_img.set(nx_file)
                ny_curr_img.set(ny_file)
                nz_curr_img.set(nz_file)
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
    
    # functions for outputting graphs
    @render_widget
    @reactive.event(selected_images, input.apix)
    def display_selected_image():
        req(len(selected_images())>0)
        print("update display_selected_image:"+str(input.apix()))

        #h, w = selected_images()[0].shape[:2]

        img = selected_images()[0]

        print("origin image to display: ", img.shape)
        print("min image to display: ", np.min(img))
        print("max image to display: ", np.max(img))

        fig = compute.create_image_figure(
            img, 
            input.apix(), 
            input.apix(), 
            title=f"Original image ({nx_curr_img()}x{ny_curr_img()})", 
            title_location="below", 
            plot_width=None, 
            plot_height=None, 
            x_axis_label=None, 
            y_axis_label=None, 
            tooltips=None, 
            show_axis=False, 
            show_toolbar=False, 
            crosshair_color="white"
        )

        return fig
    
    def update_with_ny_nx(ny_val, nx_val):
        ui.update_numeric('pnx', min=min(ny_val, 128))
        ui.update_numeric('pny', min=min(nx_val, 512))
    
    @render_widget
    @reactive.event(data_2d_transformed)
    def display_transformed_data():
        req(data_2d_transformed() is not None)
        
        img = data_2d_transformed()
        #h, w = data_2d_transformed().shape[:2]

        print("transformed image to display: ", img.shape)
        print("average image to display: ", np.average(img))
        print("apix_from_file:", apix_from_file())

        fig = compute.create_image_figure(
            img, 
            input.apix(), 
            input.apix(), 
            title=f"Transformed image ({nx_curr_img()}x{ny_curr_img()})", 
            title_location="below", 
            plot_width=None, 
            plot_height=None, 
            x_axis_label=None, 
            y_axis_label=None, 
            tooltips=None, 
            show_axis=False, 
            show_toolbar=False, 
            crosshair_color="white"
        )

        return fig
    
    @render_widget
    @reactive.event(data_2d_transformed)
    def plot_radial_profile():
        req(data_2d_transformed() is not None)
        
        data = data_2d_transformed()
        mask_radius = input.mask_radius()
        
        x = np.arange(-nx_curr_img() // 2, nx_curr_img() // 2) * input.apix()
        ymax = np.max(data, axis=0)
        ymean = np.mean(data, axis=0)

        tools = 'box_zoom,crosshair,hover,pan,reset,save,wheel_zoom'
        tooltips = [("X", "@x{0.0}Å")]
        p = figure(x_axis_label="x (Å)", y_axis_label="pixel value", frame_height=200, tools=tools, tooltips=tooltips)
        p.line(x, ymax, line_width=2, color='red', legend_label="max")
        p.line(-x, ymax, line_width=2, color='red', line_dash="dashed", legend_label="max flipped")
        p.line(x, ymean, line_width=2, color='blue', legend_label="mean")
        p.line(-x, ymean, line_width=2, color='blue', line_dash="dashed", legend_label="mean flipped")
        rmin_span = Span(location=-mask_radius, dimension='height', line_color='green', line_dash='dashed', line_width=3)
        rmax_span = Span(location=mask_radius, dimension='height', line_color='green', line_dash='dashed', line_width=3)
        p.add_layout(rmin_span)
        p.add_layout(rmax_span)
        p.yaxis.visible = False
        p.legend.visible = False
        p.legend.location = "top_right"
        p.legend.click_policy="hide"
        toggle_legend_js_x = CustomJS(args=dict(leg=p.legend[0]), code="""
            if (leg.visible) {
                leg.visible = false
                }
            else {
                leg.visible = true
            }
        """)
        p.js_on_event(DoubleTap, toggle_legend_js_x)

        return p

    @reactive.Calc
    def acf_data():
        req(input.input_type() in ["Image"])
        req(len(selected_images())>=0)
        print("updating acf")
        images = selected_images()
        # if len(images)<=0:
        #     #print("No selected images!")
        #     return None, None, None

        data = images[0]

        acf = compute.auto_correlation(data, sqrt=True, high_pass_fraction=0.1)

        if acf is None or np.isnan(acf).any():
            print("Invalid auto-correlation data!")
            return None

        return acf


    @render_widget
    @reactive.event(data_2d_transformed)
    def acf_plot():
        req(data_2d_transformed() is not None)
        print("Entering acf plots")
        acf = acf_data()
        print("acf data obtained")
        ny = acf.shape[0]
        y = np.arange(-ny//2, ny//2)*input.apix()
        xmax = np.max(acf, axis=1)

        tools = 'box_zoom,crosshair,hover,pan,reset,save,wheel_zoom'
        tooltips = [("Axial Shift", "@y{0.0}Å")]
        p = figure(x_axis_label="Auto-correlation", y_axis_label="Axial Shift (Å)", frame_height=ny, tools=tools, tooltips=tooltips)
        p.line(xmax, y, line_width=2, color='red', legend_label="ACF")
        print("acf figure created")
        # if 0:
        #     xmean = np.mean(acf, axis=1)
        #     p.line(xmean, y, line_width=2, color='blue', legend_label="mean")
        p.hover[0].attachment = "above"
        #print("acf legends:")
        #print(p.legend.__dir__())
        p.legend.visible = False
        p.legend.location = "top_right"
        p.legend.click_policy="hide"
        toggle_legend_js_y = CustomJS(args=dict(leg=p.legend[0]), code="""
            if (leg.visible) {
                leg.visible = false
                }
            else {
                leg.visible = true
            }
        """)
        p.js_on_event(DoubleTap, toggle_legend_js_y)

        return p

app = App(app_ui, server)


###########################

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

