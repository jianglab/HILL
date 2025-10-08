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
from bokeh.events import MouseMove, MouseEnter, DoubleTap, MouseLeave, Tap
from bokeh.models import Button, ColumnDataSource, CustomJS, Legend, Span, Spinner, Slider
from bokeh.models.tools import CrosshairTool, HoverTool, PointDrawTool

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
        .inline_input_box label{ display: table-cell; text-align: left; vertical-align: middle; } 
        .inline_input_box input{ width: 4em; margin-left: 1em;} 
        .inline_input_box .form-group{display: table-row;}
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
                        output_widget("plot_acf", width="auto"), # sometimes cropped with fixed size. why?
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
                        ui.input_text("const_image_color", "Flatten the PS/PD image in this color", value="", placeholder="white black", update_on="blur"),
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
                            ui.input_checkbox("PS", "Power spectra", value=True),
                            ui.input_checkbox("YP", "Y Profile", value=True),
                            ui.input_checkbox("Phase", "Phase hover", value=False),
                            ui.input_checkbox("PD", "Phase difference", value=True),
                            ui.input_checkbox("Color", "Color", value=True),
                            ui.input_checkbox("LL", "Layer line", value=True),
                            ui.input_checkbox("LLText", "Layer line text", value=True),
                            col_widths=6,
                            style="align-items: flex-end;"
                        ),
                        #ui.h5("m:"),
                        ui.div(
                            ui.input_numeric("m_max", "Max M =", value=3, min=1, step=1, update_on="blur", width="4em"),
                            class_="inline_input_box"
                        ),
                        #ui.output_ui("show_m_choices"),
                        ui.input_checkbox_group("ms","m:", ["3","2","1","0","-1","-2","-3"], selected=["-1", "0", "1"]),
                        style="display: inline; text-align: left;",
                        class_="wrap-text"
                    ),
                    ui.div(
                        ui.input_numeric("twist", "Twist (°)", value=29.40, min=-180.0, max=180.0, update_on="blur"),
                        ui.input_numeric("pitch", "Pitch (Å)", value=268.41, min=1.0, update_on="blur"),
                        ui.input_numeric("rise", "Rise (Å)", value=21.92, min=1.0, update_on="blur"),
                        hidden=True # only use for receiving values from javascript
                    )
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
                ui.head_content(ui.tags.script("""
                    function roundToTwoDecimals(num) {
                        return +(Math.round(num + "e+2") + "e-2");
                    }
                """)),
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
        # adjustable sideba width now have official support in shiny 1.5.0?
        #helicon.shiny.setup_ajdustable_sidebar(width="20vw"),
    ),
    class_="main-page-container",
)

#TODO: Change input type doesn't load initial image/change to url input doesn't take input

def server(input, output, session):
    # used to be global variables, which causes syncing issues across sessions
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
    initial_selected_image_indices = reactive.value([0])
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
    #data_2d_transformed_masked = reactive.value(None)

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

    pwr_data = reactive.value(None)
    phase_data = reactive.value(None)
    pd_data = reactive.value(None)

    user_inputs = reactive.value({
        'apix': None,
        'angle': None,
        'dx': None,
        'dy': None,
        'mask_radius': None
    })

    ################################################

    # main plots initialization
    figs = []
    figs_image = []
    figs_width = 0
    init_cutoff_res_x = 7.03
    init_cutoff_res_y = 4.69
    init_nx = 256
    init_ny = 256
    init_pnx = 512
    init_pny = 1024
    init_apix = 2.3438
    init_helical_radius = 77.6
    init_img_data = np.ones((init_ny,init_nx))
    pwr_work = np.zeros((init_pny, init_pnx))#np.random.uniform(low=0, high=1.0, size=(init_pny,init_pnx))
    pwr_work[0][0]=1.0
    phase_work = np.zeros((init_pny, init_pnx))
    phase_diff_work = np.zeros((init_pny, init_pnx))
    phase_diff_work[0][0]=1.0
    span_width = Span(dimension="width", line_color="red")
    span_height = Span(dimension="height", line_color="red")

    # PS
    tooltips_ps = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Amp', '@image')]
    fig_ps, data_source_ps = compute.create_layerline_image_figure(
        data=pwr_work, 
        cutoff_res_x=init_cutoff_res_x, 
        cutoff_res_y=init_cutoff_res_y, 
        helical_radius=init_helical_radius, 
        tilt=0.0,
        phase=phase_work,
        fft_top_only=False, 
        pseudo_color=True, 
        const_image_color="", 
        title="Power Spectra", 
        yaxis_visible=False, 
        tooltips=tooltips_ps
    )
    
    figs.append(fig_ps)
    figs_image.append(fig_ps)
    figs[-1].add_tools(CrosshairTool(overlay=(span_width, span_height)))
    figs_width += pwr_work.shape[-1]

    # YP
    ny, nx = (init_pny, init_pnx)
    dsy = 1/(ny//2*init_cutoff_res_y)
    y=np.arange(-ny//2, ny//2)*dsy
    yinv = y*1.0
    yinv[yinv==0] = 1e-10
    yinv = 1/np.abs(yinv)
    yprofile = np.mean(pwr_work, axis=1)
    yprofile /= yprofile.max()
    data_source_yp = ColumnDataSource(data=dict(yprofile=yprofile, y=y, resy=yinv))
    tools_yp = 'box_zoom,hover,pan,reset,save,wheel_zoom'
    tooltips_yp = [('Res y', '@resy Å'), ('Amp', '$x')]
    fig_yp = figure(frame_width=nx//2, frame_height=figs[-1].frame_height, y_range=figs[-1].y_range, y_axis_location = "right", title=None, tools=tools_yp, tooltips=tooltips_yp)
    fig_yp.line(source=data_source_yp, x='yprofile', y='y', line_width=2, color='blue')
    fig_yp.yaxis.visible = False
    fig_yp.hover[0].attachment="vertical"
    figs.append(fig_yp)
    figs[-1].add_tools(CrosshairTool(overlay=(span_width)))
    figs_width += nx//2

    # PD
    tooltips_pd = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Phase Diff', '@image °')]
    fig_pd, data_source_pd = compute.create_layerline_image_figure(
        data=phase_diff_work, 
        cutoff_res_x=init_cutoff_res_x, 
        cutoff_res_y=init_cutoff_res_y, 
        helical_radius=init_helical_radius, 
        tilt=0.0, 
        phase=phase_work, 
        fft_top_only=False, 
        pseudo_color=True, 
        const_image_color="", 
        title="Phase Diff Across Meridian", 
        yaxis_visible=False, 
        tooltips=tooltips_pd
    )
    figs.append(fig_pd)
    figs_image.append(fig_pd)
    figs_width += phase_diff_work.shape[-1]       

    # add synced crosshairs
    figs[-1].yaxis.fixed_location = figs[-1].x_range.end
    figs[-1].yaxis.visible = True
    figs[-1].add_tools(CrosshairTool(overlay=(span_width, span_height)))
    #add_linked_crosshair_tool(figs, overlay=(crosshair_width))
    #add_linked_crosshair_tool(figs_image, overlay=(crosshair_width, crosshair_height))

    # add layer lines
    curr_twist = 29.40
    curr_rise = 21.92
    curr_pitch = 268.41
    fig_ellipses = []
    curr_m_groups = compute.compute_layer_line_positions(
            twist=29.40,
            rise=21.92,
            csym=6,
            radius=100.0,
            tilt=0.0,
            cutoff_res=init_cutoff_res_y,
            m_max=3
        )
    curr_ll_colors="lime cyan violet salmon silver".split()
    curr_ms=["-1","0","1"]
    curr_ll_text = True
    curr_ll = True
    if figs_image and curr_ll:
        if max(curr_m_groups[0]["LL"][0])>0:
            x, y, n = curr_m_groups[0]["LL"]
            tmp_x = np.sort(np.unique(x))
            width = np.mean(tmp_x[1:]-tmp_x[:-1])
            height = width/5
            for mi, m in enumerate(curr_m_groups.keys()):
                print("processing m:", m)
                #if str(m) not in curr_ms: continue
                x, y, bessel_order = curr_m_groups[m]["LL"]
                if curr_ll_text:
                    texts = [str(int(n)) for n in bessel_order]
                tags = [m, bessel_order]
                color = curr_ll_colors[abs(m)%len(curr_ll_colors)]
                ellipse_alpha = [n%2*1.0 for n in bessel_order]
                for f in figs_image:
                    if curr_ll_text:
                        text_labels = f.text(x, y, y_offset=0.0, text=texts, text_color=color, text_baseline="middle", text_align="center", visible=str(m) in curr_ms) # y offset = 2.0
                        text_labels.tags = tags
                        #text_labels.level = "overlay" # don't add, it makes texts go out of plots
                        #print("text labels:", text_labels)
                        fig_ellipses.append(text_labels)
                    else:
                        ellipses = f.ellipse(x, y, width=width, height=height, line_color=color, fill_color=color, fill_alpha=ellipse_alpha, line_width=1.0, visible=str(m) in curr_ms)
                        ellipses.tags = tags
                        fig_ellipses.append(ellipses)
    
    # add helical parameter controls
    slider_width = figs_width//3 if len(figs)>1 else figs_width
    #from bokeh.models import Slider, CustomJS
    spinner_twist = Spinner(title='Twist (°)', low=-360.0, high=360.0, step=1.0, value=curr_twist, format="0.00", width=slider_width)
    spinner_pitch = Spinner(title='Pitch (Å)', low=1.0, step=1.0, value=curr_pitch, format="0.00", width=slider_width)
    spinner_rise = Spinner(title='Rise (Å)', low=1.0, step=1.0, value=curr_rise, format="0.00", width=slider_width)

    slider_twist = Slider(start=-180.0, end=180.0, value=curr_twist, step=0.01, title="Twist (°)", width=slider_width)
    slider_pitch = Slider(start=curr_pitch/2, end=curr_pitch*2.0, value=curr_pitch, step=curr_pitch*0.002, title="Pitch (Å)", width=slider_width)
    slider_rise = Slider(start=curr_rise/2, end=min(curr_pitch, curr_rise*2.0), value=curr_rise, step=min(curr_pitch, curr_rise*2.0)*0.001, title="Rise (Å)", width=slider_width)
    
    #https://stackoverflow.com/questions/11832914/how-to-round-to-at-most-2-decimal-places-if-necessary
    callback_spinner_rise_change_code = """
        Shiny.setInputValue("rise", slider_rise.value, {priority: 'event'})
        slider_rise.value = spinner_rise.value
    """
    callback_spinner_pitch_change_code = """
        Shiny.setInputValue("pitch", slider_pitch.value, {priority: 'event'})
        slider_pitch.value = spinner_pitch.value
    """
    callback_spinner_twist_change_code = """
        Shiny.setInputValue("twist", slider_twist.value, {priority: 'event'})
        slider_twist.value = spinner_twist.value
    """
    callback_spinner_twist_change = CustomJS(args=dict(spinner_twist=spinner_twist, spinner_pitch=spinner_pitch, spinner_rise=spinner_rise,slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise), code=callback_spinner_twist_change_code)
    callback_spinner_pitch_change = CustomJS(args=dict(spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise,slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise), code=callback_spinner_pitch_change_code)
    callback_spinner_rise_change = CustomJS(args=dict(spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise,slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise), code=callback_spinner_rise_change_code)

    spinner_twist.js_on_change('value', callback_spinner_twist_change)
    spinner_pitch.js_on_change('value', callback_spinner_pitch_change)
    spinner_rise.js_on_change('value', callback_spinner_rise_change)

    callback_slider_rise_change_code = """
        var twist_sign = 1.
        if (slider_twist.value < 0) {
            twist_sign = -1.
        }
        console.log("rise_cb_use_twist_pitch:", $("input[type='radio'][name='use_twist_pitch']:checked").val())
        if ($("input[type='radio'][name='use_twist_pitch']:checked").val() == "Pitch") {
            var slider_twist_to_update = twist_sign * 360/(slider_pitch.value/slider_rise.value)
            if (slider_twist_to_update != slider_twist.value) {
                slider_twist.value = slider_twist_to_update
            }
        }
        else {
            var slider_pitch_to_update = Math.abs(360/slider_twist.value * slider_rise.value)
            if (slider_pitch_to_update != slider_pitch.value) {
                slider_pitch.value = slider_pitch_to_update
            }
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
    callback_slider_pitch_change_code = """
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
    callback_slider_twist_change_code = """
        var slider_pitch_to_update = Math.abs(360/slider_twist.value * slider_rise.value)                   
        if (slider_pitch_to_update != slider_pitch.value) {
            slider_pitch.value = slider_pitch_to_update
        }
        if (spinner_twist.value != slider_twist.value) {
            spinner_twist.value = slider_twist.value
        } 
    """
    callback_slider_rise_change = CustomJS(args=dict(fig_ellipses=fig_ellipses, slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_slider_rise_change_code)
    callback_slider_pitch_change = CustomJS(args=dict(fig_ellipses=fig_ellipses, slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_slider_pitch_change_code)
    callback_slider_twist_change = CustomJS(args=dict(slider_twist=slider_twist, slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_slider_twist_change_code)
    slider_twist.js_on_change('value', callback_slider_twist_change)
    slider_pitch.js_on_change('value', callback_slider_pitch_change)
    slider_rise.js_on_change('value', callback_slider_rise_change)

    callback_spinner_rise_change_throttled_code = """
        slider_rise.start = spinner_rise.value/2.0
        slider_rise.end = Math.min(spinner_rise.value*2.0, spinner_pitch.value)
        slider_rise.step = slider_rise.end*0.001
        spinner_pitch.low = slider_rise.value
    """
    callback_spinner_pitch_change_throttled_code = """
        slider_pitch.start = Math.max(spinner_pitch.value/2.0, spinner_rise.value)
        slider_pitch.end = Math.min(spinner_pitch.value*2.0, 10000.0)
        slider_pitch.step = slider_pitch.end*0.001
        spinner_rise.high = slider_pitch.value
    """
    callback_spinner_rise_change_throttled = CustomJS(args=dict(spinner_rise=spinner_rise, slider_rise=slider_rise, spinner_pitch=spinner_pitch), code=callback_spinner_rise_change_throttled_code)
    callback_spinner_pitch_change_throttled = CustomJS(args=dict(spinner_pitch=spinner_pitch, slider_pitch=slider_pitch, spinner_rise=spinner_rise), code=callback_spinner_pitch_change_throttled_code)
    spinner_rise.js_on_change('value_throttled', callback_spinner_rise_change_throttled)
    spinner_pitch.js_on_change('value_throttled', callback_spinner_pitch_change_throttled)

    figs_row = gridplot(children=[figs], toolbar_location='right')
    figs_grid = layout(children=[[spinner_twist, spinner_pitch, spinner_rise],[slider_twist, slider_pitch, slider_rise], figs_row])
    main_plot_widget = BokehModel(figs_grid)

    @reactive.effect
    @reactive.event(input.twist)
    def test_input_rise_update():
        print("input.twist updated to:", input.twist())
    @reactive.effect
    @reactive.event(input.pitch)
    def test_input_rise_update():
        print("input.pitch updated to:", input.pitch())
    @reactive.effect
    @reactive.event(input.rise)
    def test_input_rise_update():
        print("input.rise updated to:", input.rise())

    @reactive.effect
    @reactive.event(input.sidebar_navset)
    def main_page_tab_switch():
        print(input.sidebar_navset())
        if input.sidebar_navset() == "Filament Straightening":
            ui.update_navset("main_tabs", selected="straighten_filament_tab")
            print(input.apix())
        else:
            ui.update_navset("main_tabs", selected="main_content_tab")

    # @output
    # @render.text
    # def m_max_auto():
    #     return int(np.floor(np.abs(input.rise()/input.res_limit_y()))) + 3
    
    @reactive.effect
    @reactive.event(input.m_max, ignore_init=True)
    def update_m_choices():
        print("updating ms groups")
        prev_selected_ms = input.ms()
        ms = [ str(m) for m in range(input.m_max(), -input.m_max()-1, -1)]
        ui.update_checkbox_group("ms", choices=ms, selected=prev_selected_ms)

    @reactive.effect
    @reactive.event(input.m_max, input.filament_diameter, input.csym, input.out_of_plane_tilt, input.res_limit_y, input.ll_colors, input.LL, input.LLText, input.pnx, input.pny, ignore_init=True)
    def update_layerlines_all():
        print("updating plot sizes") # TODO: plot size updates only responses instantly here. Add a similar trigger for the frame size update
        fig_ps.frame_height = input.pny()
        fig_ps.frame_width = input.pnx()
        fig_yp.frame_height = input.pny()
        fig_yp.frame_width = input.pnx()//2
        fig_pd.frame_height = input.pny()
        fig_pd.frame_width = input.pnx()
        for e in fig_ellipses:
            for f in figs_image:
                if e in f.renderers:
                    f.renderers.remove(e)
        fig_ellipses.clear()
        if input.LL():
            curr_m_groups = compute.compute_layer_line_positions(
                twist=input.twist(),
                rise=input.rise(),
                csym=input.csym(),
                radius=input.filament_diameter()/2,
                tilt=input.out_of_plane_tilt(),
                cutoff_res=input.res_limit_y(),
                m_max=input.m_max()
            )
            curr_ll_colors = input.ll_colors().split()
            if max(curr_m_groups[0]["LL"][0])>0:
                x, y, n = curr_m_groups[0]["LL"]
                tmp_x = np.sort(np.unique(x))
                width = np.mean(tmp_x[1:]-tmp_x[:-1])
                height = width/5
                for mi, m in enumerate(curr_m_groups.keys()):
                    print("reprocessing m:", m)
                    #if str(m) not in curr_ms: continue
                    x, y, bessel_order = curr_m_groups[m]["LL"]
                    if input.LLText():
                        texts = [str(int(n)) for n in bessel_order]
                    tags = [m, bessel_order]
                    color = curr_ll_colors[abs(m)%len(curr_ll_colors)]
                    ellipse_alpha = [n%2*1.0 for n in bessel_order]
                    for f in figs_image:
                        if input.LLText():
                            text_labels = f.text(x, y, y_offset=0.0, text=texts, text_color=color, text_baseline="middle", text_align="center", visible=str(m) in input.ms()) # y offset = 2.0
                            text_labels.tags = tags
                            fig_ellipses.append(text_labels)
                        else:
                            ellipses = f.ellipse(x, y, width=width, height=height, line_color=color, fill_color=color, fill_alpha=ellipse_alpha, line_width=1.0, visible=str(m) in input.ms())
                            ellipses.tags = tags
                            fig_ellipses.append(ellipses)
        callback_slider_rise_change = CustomJS(args=dict(fig_ellipses=fig_ellipses, slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_slider_rise_change_code)
        callback_slider_pitch_change = CustomJS(args=dict(fig_ellipses=fig_ellipses, slider_twist=slider_twist,slider_pitch=slider_pitch, slider_rise=slider_rise, spinner_twist=spinner_twist,spinner_pitch=spinner_pitch, spinner_rise=spinner_rise), code=callback_slider_pitch_change_code)
        slider_pitch.js_on_change('value', callback_slider_pitch_change)
        slider_rise.js_on_change('value', callback_slider_rise_change)
    
    @reactive.effect
    @reactive.event(input.const_image_color, input.Color)
    def update_color_palette():
        from bokeh.models.glyphs import Image as ImageGlyph
        from bokeh.palettes import Viridis256, Greys256
        if input.const_image_color() != "":
            new_palette = tuple(input.const_image_color().split())
        else:
            if input.Color():
                new_palette = Viridis256
            else:
                new_palette = Greys256
        for f in figs_image:
            for rend in f.renderers:
                if getattr(rend, "glyph", None) and isinstance(rend.glyph, ImageGlyph):
                    rend.glyph.color_mapper.palette = new_palette
    
    # @reactive.effect
    # @reactive.event(input.pny, input.pnx)
    # def update_plot_frame_sizes():
    #     print("updating plot sizes")
    #     fig_ps.frame_height = input.pny()
    #     fig_ps.frame_width = input.pnx()
    #     fig_yp.frame_height = input.pny()
    #     fig_yp.frame_width = input.pnx()//2
    #     fig_pd.frame_height = input.pny()
    #     fig_pd.frame_width = input.pnx()
    
    @reactive.effect
    @reactive.event(input.Phase)
    def toggle_phase_hover():
        ps_h = next(t for t in fig_ps.tools if isinstance(t, HoverTool))
        pd_h = next(t for t in fig_pd.tools if isinstance(t, HoverTool))
        if input.Phase():
            ps_h.tooltips = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Amp', '@image'), ("Phase", "@phase °")]
            pd_h.tooltips = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Phase Diff', '@image °'), ("Phase", "@phase °")]
        else:
            ps_h.tooltips = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Amp', '@image')]
            pd_h.tooltips = [("Res r", "Å"), ('Res y', 'Å'), ('Res x', 'Å'), ('Jn', '@bessel'), ('Phase Diff', '@image °')]            
    
    @reactive.effect
    @reactive.event(input.filament_diameter, input.out_of_plane_tilt, input.res_limit_x, input.res_limit_y, input.fft_top_only, input.pnx, input.pny, ignore_init=True)
    def update_figure_bessel():
        # #ny, nx = figs[0].renderers[0].data_source.data['image'][0].shape
        # ny, nx = np.shape(data_2d_transformed())
        # bessel = compute.bessel_n_image(ny, nx, input.res_limit_x(), input.res_limit_y(), input.filament_diameter()/2.0, input.out_of_plane_tilt()).astype(np.int16)
        # data_source_ps.data = {
        #     **data_source_ps.data,
        #     "bessel": [bessel],
        # }
        # data_source_pd.data = {
        #     **data_source_pd.data,
        #     "bessel": [bessel],
        # }
        print("update_figure_bessel")
        ny = input.pny()
        nx = input.pnx()
        dsy = 1/(ny//2*input.res_limit_y())
        dsx = 1/(nx//2*input.res_limit_x())
        x_range = (-(nx//2+0.5)*dsx, (nx//2-0.5)*dsx)
        if input.fft_top_only():
            y_range = (-(ny//2 * 0.01)*dsy, (ny//2-0.5)*dsy)
        else:
            y_range = (-(ny//2+0.5)*dsy, (ny//2-0.5)*dsy)

        bessel = compute.bessel_n_image(ny, nx, input.res_limit_x(), input.res_limit_y(), input.filament_diameter()/2.0, input.out_of_plane_tilt()).astype(np.int16)

        # tools = 'box_zoom,pan,reset,save,wheel_zoom'
        # fig = figure(title_location="below", frame_width=nx, frame_height=ny, 
        #     x_axis_label=None, y_axis_label=None, x_range=x_range, y_range=y_range, tools=tools)
        # fig.grid.visible = False
        # fig.title.text = title
        # fig.title.align = "center"
        # fig.title.text_font_size = "20px"
        # fig.yaxis.visible = yaxis_visible   # leaving yaxis on will make the crosshair x-position out of sync with other figures

        data_source_ps.data = {
            **data_source_ps.data,
            "x":[-nx//2*dsx], 
            "y":[-ny//2*dsy], 
            "dw":[nx*dsx], 
            "dh":[ny*dsy], 
            "bessel":[bessel]
        }
        data_source_pd.data = {
            **data_source_pd.data,
            "x":[-nx//2*dsx], 
            "y":[-ny//2*dsy], 
            "dw":[nx*dsx], 
            "dh":[ny*dsy], 
            "bessel":[bessel]
        }

        print("set x_range:", x_range, "y_range:", y_range)
        for f in figs:
            f.y_range.start = y_range[0]
            f.y_range.end = y_range[1]
        for f in figs_image:
            f.x_range.start = x_range[0]
            f.x_range.end = x_range[1]

    @output
    @render.ui
    def filament_straightening_uis():
        return ui.div(
            output_widget("figure_get_markers_straighten", width="50%", height="auto"),
            ui.input_action_button("reset_markers", "Clear markers"),
            ui.input_action_button("straighten_run", "Straighten"),
            style="text-align: center; border: 1px solid #ddd; padding: 10px;",
            class_="wrap-text"
        )

    fig_straighten, data_source_straighten = compute.create_image_figure(
        init_img_data, 
        init_apix, 
        init_apix, 
        title=f"Original image ({init_nx}x{init_ny})", 
        title_location="below", 
        plot_width=None, 
        plot_height=None, 
        x_axis_label=None, 
        y_axis_label=None, 
        tooltips=None, 
        show_axis=False, 
        show_toolbar=True, 
        crosshair_color="white",
    )
    # straightening_draw_tool = PointDrawTool(renderers=fig_straighten.renderers, empty_value="black")
    # fig_straighten.add_tools(straightening_draw_tool)
    # fig_straighten.toolbar.active_tap = straightening_draw_tool
    source = ColumnDataSource({
        'x': [], 'y': []
    })

    straighten_marker_renderer = fig_straighten.scatter(x='x', y='y', source=source, color='red', size=20)

    draw_tool = PointDrawTool(renderers=[straighten_marker_renderer], default_overrides={"color":"red", "size":20})
    fig_straighten.add_tools(draw_tool)
    fig_straighten.toolbar.active_tap = draw_tool
    fig_straighten.js_on_event(Tap, CustomJS(code="console.log('straightening figure clicked')"))


    @render_widget
    def figure_get_markers_straighten():
        return fig_straighten

    @reactive.effect
    @reactive.event(input.reset_markers)
    def clear_markers():
        for m in fig_straighten:
            fig_straighten.renderers.remove(m)
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
        new_xs, tck = compute.fit_spline(data_2d_transformed(), xs, ys, input.apix(), display=False)
        data_2d_transformed.set(np.nan_to_num(compute.filament_straighten(None,data_2d_transformed(), tck, new_xs, ys, input.output_width_straighten()//2,input.apix()), 0.0))
    
    @reactive.effect
    @reactive.event(input.sidebar_navset)
    def update_figure_get_markers_straighten():
        req(input.sidebar_navset() == "Filament Straightening")
        req(len(selected_images())>0)
        print("update straightening figure")
        
        curr_img = selected_images()[0]

        h, w = curr_img.shape
        dx = input.apix()
        dy = input.apix()
        aspect_ratio = w*dx/(h*dy)
        fig_straighten.x_range.start = -w//2*dx
        fig_straighten.x_range.end = (w//2-1)*dx
        fig_straighten.y_range.start = (h//2-1)*dy
        fig_straighten.y_range.end = -h//2*dy
        fig_straighten.aspect_ratio = aspect_ratio
        fig_straighten.title.text = f"Original image ({nx_curr_img()}x{ny_curr_img()})"

        data_source_straighten.data = {
            **data_source_straighten.data,
            "image":[curr_img], 
            "x":[-w//2*dx], 
            "y":[-h//2*dy], 
            "dw":[w*dx], 
            "dh":[h*dy]
        }
        # img = data_2d_transformed()
        # h, w = img.shape
        # apix = input.apix()
        # marker_r_angst = input.template_diameter_straighten()/2.0 * apix

        # fig = striaghten_fig()
        # if len(fig.data) > 0:
        #     print("update straighten plot")
        #     fig.update_traces(z=data_2d_transformed(),selector=dict(name="image_to_straighten"))
        # else:
        #     print("generating straighten plot")
        #     fig.add_trace(
        #         go.Heatmap(
        #             name="image_to_straighten",
        #             z=img,
        #             x=np.arange(w) * apix,
        #             y=np.arange(h) * apix,
        #             colorscale="gray",
        #             showscale=False,  # Removed colorbar
        #             #reversescale=True, # Greys colorscale is white to black by default??? why?
        #             hoverongaps=False,
        #             hovertemplate=(
        #                 "x: %{x:.1f} Å<br>"
        #                 + "y: %{y:.1f} Å<br>"
        #                 + "val: %{z:.1f}<br>"
        #                 + "<extra></extra>"
        #             ),
        #         )
        #     )
        #     fig.add_trace(
        #         go.Line(
        #             name="fitted_spline",
        #             x=[],
        #             y=[],
        #         ),
        #     )
        #     fig.update_layout(width=w, height=h)
        #     fig.update_xaxes(showticklabels=False)
        #     fig.update_yaxes(showticklabels=False, autorange="reversed")

        #     def update_markers(trace,points,state):
        #         click_x = points.xs[0]
        #         click_y = points.ys[0]
        #         markers_straighten_x.set(markers_straighten_x()+[click_x])
        #         markers_straighten_y.set(markers_straighten_y()+[click_y])
        #         fig.layout.shapes = list(fig.layout.shapes) + [dict(type='circle',x0=click_x-marker_r_angst, y0=click_y-marker_r_angst, x1=click_x+marker_r_angst, y1=click_y+marker_r_angst, line_color='red', name='markers')]

            
        #     fig.data[0].on_click(update_markers)
        #     markers_straighten_x.set([])
        #     markers_straighten_y.set([])
        # return fig
    
    @render_widget
    def main_plots():
        return main_plot_widget

    @reactive.effect
    @reactive.event(input.PS, input.YP, input.PD, ignore_init=True)
    def toggle_plot_visibility():
        if input.PS():
            figs[0].visible = True
        else:
            figs[0].visible = False    
        if input.YP():
            figs[1].visible = True
        else:
            figs[1].visible = False
        if input.PD():
            figs[2].visible = True
        else:
            figs[2].visible = False
    
    @reactive.effect
    @reactive.event(input.ms)
    def toggle_layerline_visibility():
        #req(m_groups() is not None)
        selected_ms =[str(m) for m in input.ms()]
        print("toggling m visibility:", selected_ms)
        for ellipses in fig_ellipses:
            if str(ellipses.tags[0]) in selected_ms:
                ellipses.visible = True
            else:
                ellipses.visible = False

    #@render_widget
    @reactive.effect
    @reactive.event(pwr_data)
    def update_figure_pwr_yp_data():
        req(pwr_data() is not None)
        print("updating pwr data")
        pwr_work = pwr_data()
        phase_work = phase_data()
        data_source_ps.data = {
            **data_source_ps.data,
            "image": [np.asarray(pwr_work)],
            "phase": [np.fmod(np.rad2deg(phase_work)+360, 360).astype(np.float16)]
        }
        # dsy = 1/(input.pny()//2*input.res_limit_y())
        # if input.fft_top_only():
        #     y_range = (-(input.pny()//2 * 0.01)*dsy, (input.pny()//2-0.5)*dsy)
        # else:
        #     y_range = (-(input.pny()//2+0.5)*dsy, (input.pny()//2-0.5)*dsy)

        # fig_yp.y_range.start = y_range[0]
        # fig_yp.y_range.end = y_range[1]
        ny = input.pny()
        dsy = 1/(ny//2*input.res_limit_y())        
        y=np.arange(-ny//2, ny//2)*dsy
        yinv = y*1.0
        yinv[yinv==0] = 1e-10
        yinv = 1/np.abs(yinv)
        yprofile = np.mean(pwr_work, axis=1)
        yprofile /= yprofile.max()
        data_source_yp.data = {
            **data_source_yp.data,
            "yprofile": yprofile,
            "y":y,
            "resy":yinv
        }
    @reactive.effect
    @reactive.event(pd_data)
    def update_figure_pd_data():
        req(pd_data() is not None)
        print("updating pd data")
        pd_work = pd_data()
        phase_work = phase_data()
        data_source_pd.data = {
            **data_source_pd.data,
            "image": [np.asarray(pd_work)],
            "phase": [np.fmod(np.rad2deg(phase_work)+360, 360).astype(np.float16)]
        }

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
                    value=0.0,
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
                    value=0.0,
                    min=-180.0, 
                    max=180.0, 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "dx", 
                    "Shift along X-dim (Å)", 
                    value=0.0,
                    min=-nx_curr_img()*apix_from_file(), 
                    max=nx_curr_img()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
                ui.input_numeric(
                    "dy", 
                    "Shift along Y-dim (Å)", 
                    value=0.0,
                    min=-ny_curr_img()*apix_from_file(), 
                    max=ny_curr_img()*apix_from_file(), 
                    step=1.0,
                    update_on="blur"
                ),
            ]

    # @render.ui
    # @reactive.event(input.is_3d, input.input_type)
    # def add_transpose():
    #     transpose_auto = (input.input_type() in ["Image"]) and (nx_curr_img() > ny_curr_img())
    #     if not input.is_3d():
    #         return ui.input_checkbox("transpose", "Transpose the image", value=transpose_auto),

    @reactive.effect
    @reactive.event(input.apix, input.angle, input.dx, input.dy, input.mask_radius, input.negate, input.mask_len) #TODO: temp solution: remove event on select_image or auto_transformation_estimated. Need to check potential pitfalls
    def set_data_2d_transformed():
        req(input.input_type() in ["Image"])
        req(len(selected_images())>0)
        req(input.apix()!=0)
        #req(auto_transformation_estimated())
        temp_transformed = selected_images()[0]
        if (input.angle() or input.dx() or input.dy() or input.negate()):
            temp_transformed = compute.transform_2d_filament(temp_transformed,input.angle(),input.dx(),input.dy(),input.negate(),input.apix())
            if (input.mask_radius()>0.0 and input.mask_len()>0.0 and input.mask_len()<100.0):
                temp_transformed = compute.mask_2d_filament(temp_transformed,input.mask_radius(), input.apix(), input.mask_len()/100.0)
        data_2d_transformed.set(temp_transformed)
    
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
        phase_data.set(phase) # update order important? maybe think of a better implementation to update figures
        pwr_data.set(pwr)
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
                    update_on="blur",
                ),
                ui.input_checkbox("is_3d", "The input is a 3D map", value=False),
            ]
        elif selection == "3":  # EMD-xxxxx
            return [
                    ui.input_text("img_file_emd_id", "Input an EMDB ID (emd-xxxxx):", value="emd-10499", update_on="blur"),
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
        req(not input.is_3d())
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
        emdb = helicon.dataset.EMDB()
        ids = emdb.amyloid_atlas_ids()
        import random

        emdb_id = f"EMD-{random.choice(ids)}"
        ui.update_text("img_file_emd_id", value=emdb_id)

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

    @output
    @render.ui
    @reactive.event(input.apply_helical_sym)
    def apply_helical_sym_uis():
        req(input.apply_helical_sym())
        return [
                ui.input_numeric("twist_ahs", "Twist (°):", value=input.twist(), min=-180.0, max=180.0, step=1.0),
                ui.input_numeric("rise_ahs", "Rise (Å):", value=input.rise(), min=0.0, step=1.0),
                ui.input_numeric("csym_ahs", "Csym:", value=input.csym(), min=1, step=1),
                ui.input_numeric("apix_map", "Current map pixel size (Å):", value=apix_from_file(), min=0.0, step=1.0),
                ui.input_numeric("apix_ahs", "New map pixel size (Å):", value=apix_from_file(), min=0.0, step=1.0),
                ui.input_numeric("fraction_ahs","Center fraction (0-1):",value=1.0,min=0,max=1.0,step=0.1),
                ui.input_numeric("length_ahs","Box length (Pixel):",value=max(nx_curr_img(),nz_curr_img()),min=0,step=1),
                ui.input_numeric("width_ahs", "Box width (Pixel):", value=max(nx_curr_img(),nz_curr_img()), min=0, step=1),
                ui.hr(),
        ]

    @reactive.effect
    @reactive.event(input_data, input.is_3d) # , input.ignore_blank)
    def get_data_all_2d_from_input_data():
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

    @reactive.effect(priority=10)
    @reactive.event(input.is_3d)
    def clear_transformation_for_3d():
        if input.is_3d():
            ui.update_numeric('angle', value=0.0)
            ui.update_numeric('dx', value=0.0)
            ui.update_numeric('dy', value=0.0)
            ui.update_numeric('mask_radius', value=0.0)
            ui.update_numeric('mask_len', value=100.0)

    @reactive.effect
    @reactive.event(input.generate_2d_projection, input_data)
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
        
        #auto_transformation_estimated.set(False)
        data_all_2d.set(proj)
        data_all_2d_labels.set(["Projection"])
        selected_images.set(proj)
        selected_image_labels.set(["Projection"])
        data_2d_transformed.set(proj[0])
        # ui.update_numeric('angle', value=0.0)
        # ui.update_numeric('dx', value=0.0)
        # ui.update_numeric('dy', value=0.0)
        # ui.update_numeric('mask_radius', value=0.0)
        # ui.update_numeric('mask_len', value=100.0)
        #auto_transformation_estimated.set(True)

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
            #auto_transformation_estimated.set(False)
            selected_images.set(
                [data_all_2d()[i] for i in input.select_image()]
            )
            selected_image_labels.set(
                [data_all_2d_labels()[i] for i in input.select_image()]
            )

            mode = input.input_type() 

            if mode in ["PS", "PD"] or input.is_3d():
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
                ui.update_numeric('dx', value=0.0)
                ui.update_numeric('dy', value=round(dy_auto_v*input.apix(),2))
                ui.update_numeric('mask_radius', value=round(diameter/2*input.apix(),2))
                mask_len_percent_auto_default = 90.0
                if input.straightening():
                    mask_len_percent_auto_default = 100.0
                #mask_len_percent_auto.set(mask_len_percent_auto_default)
                ui.update_numeric('mask_len', value=mask_len_percent_auto_default)
                #auto_transformation_estimated.set(True)
            print(selected_images())
            print("Auto transformation parameters set: angle=", input.angle())            
            
    @reactive.effect
    @reactive.event(input.img_file_upload)
    def get_images_from_input_upload():
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
    @reactive.event(input.img_file_url, input.input_mode_params)
    def get_images_from_input_url():
        req(input.input_mode_params()=="2")
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

    #############################################################
    @reactive.effect
    @reactive.event(input.input_mode_params, input.img_file_emd_id)
    def get_images_from_input_emdb():
        req(input.input_mode_params() == "3")
        emdb_id = input.img_file_emd_id()
        req(len(emdb_id) > 0)
        try:
            data, apix = util.get_images_from_emdb(emdb_id=emdb_id)
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
            import traceback

            traceback.print_exc()
            print(e)
            data, apix = None, 1
            m = ui.modal(
                f"failed to obtain {emdb_id} map from EMDB",
                title="File download error",
                easy_close=True,
                footer=None,
            )
            ui.modal_show(m)
            return

    def update_with_apix_from_file(apix):
        ui.update_numeric('apix', value=apix)
        ui.update_numeric('res_limit_x', value=round(3*apix_from_file(),4), min=2*apix_from_file())
        ui.update_numeric('res_limit_y', value=round(2*apix_from_file(),4), min=2*apix_from_file())

    fig_orig_img, data_source_orig_img = compute.create_image_figure(
        init_img_data, 
        init_apix, 
        init_apix, 
        title=f"Original image ({init_nx}x{init_ny})", 
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

    @render_widget
    def display_selected_image():
        return fig_orig_img

    # functions for outputting graphs
    @reactive.effect
    @reactive.event(selected_images, input.apix)
    def update_display_selected_image():
        req(len(selected_images())>0)
        print("update display_selected_image:"+str(input.apix()))
        
        curr_img = selected_images()[0]

        h, w = curr_img.shape
        dx = input.apix()
        dy = input.apix()
        aspect_ratio = w*dx/(h*dy)
        fig_orig_img.x_range.start = -w//2*dx
        fig_orig_img.x_range.end = (w//2-1)*dx
        fig_orig_img.y_range.start = (h//2-1)*dy # flipped display
        fig_orig_img.y_range.end = -h//2*dy
        fig_orig_img.aspect_ratio = aspect_ratio
        fig_orig_img.title.text = f"Original image ({nx_curr_img()}x{ny_curr_img()})"

        data_source_orig_img.data = {
            **data_source_orig_img.data,
            "image":[curr_img], 
            "x":[-w//2*dx], 
            "y":[-h//2*dy], 
            "dw":[w*dx], 
            "dh":[h*dy]
        }

    def update_with_ny_nx(ny_val, nx_val):
        ui.update_numeric('pnx', min=min(ny_val, 128))
        ui.update_numeric('pny', min=min(nx_val, 512))

    fig_transformed_img, data_source_transformed_img = compute.create_image_figure(
        init_img_data, 
        init_apix, 
        init_apix, 
        title=f"Transformed image ({init_nx}x{init_ny})", 
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

    @render_widget
    def display_transformed_data():
        return fig_transformed_img

    # functions for outputting graphs
    @reactive.effect
    @reactive.event(data_2d_transformed)
    def update_display_transformed_data():
        req(data_2d_transformed() is not None)
        
        curr_img = data_2d_transformed()

        h, w = curr_img.shape
        dx = input.apix()
        dy = input.apix()
        aspect_ratio = w*dx/(h*dy)
        fig_transformed_img.x_range.start = -w//2*dx
        fig_transformed_img.x_range.end = (w//2-1)*dx
        fig_transformed_img.y_range.start = (h//2-1)*dy # flipped display
        fig_transformed_img.y_range.end = -h//2*dy
        fig_transformed_img.aspect_ratio = aspect_ratio
        fig_transformed_img.title.text = f"Transformed image ({nx_curr_img()}x{ny_curr_img()})"

        data_source_transformed_img.data = {
            **data_source_transformed_img.data,
            "image":[curr_img], 
            "x":[-w//2*dx], 
            "y":[-h//2*dy], 
            "dw":[w*dx], 
            "dh":[h*dy]
        }

    #############################
    # radial profile init and update

    init_mask_radius = 100.0
    
    radial_profile_x = np.arange(-init_nx // 2, init_nx // 2) * init_apix
    radial_profile_ymax = np.max(init_img_data, axis=0)
    radial_profile_ymean = np.mean(init_img_data, axis=0)

    radial_profile_tools = 'box_zoom,crosshair,hover,pan,reset,save,wheel_zoom'
    radial_profile_tooltips = [("X", "@x{0.0}Å")]
    radial_profile_p = figure(x_axis_label="x (Å)", y_axis_label="pixel value", frame_height=200, tools=radial_profile_tools, tooltips=radial_profile_tooltips)
    radial_profile_line_max = radial_profile_p.line(radial_profile_x, radial_profile_ymax, line_width=2, color='red', legend_label="max")
    radial_profile_line_max_flipped = radial_profile_p.line(-radial_profile_x, radial_profile_ymax, line_width=2, color='red', line_dash="dashed", legend_label="max flipped")
    radial_profile_line_mean = radial_profile_p.line(radial_profile_x, radial_profile_ymean, line_width=2, color='blue', legend_label="mean")
    radial_profile_line_mean_flipped = radial_profile_p.line(-radial_profile_x, radial_profile_ymean, line_width=2, color='blue', line_dash="dashed", legend_label="mean flipped")
    radial_profile_rmin_span = Span(location=-init_mask_radius, dimension='height', line_color='green', line_dash='dashed', line_width=3)
    radial_profile_rmax_span = Span(location=init_mask_radius, dimension='height', line_color='green', line_dash='dashed', line_width=3)
    radial_profile_p.add_layout(radial_profile_rmin_span)
    radial_profile_p.add_layout(radial_profile_rmax_span)
    radial_profile_p.yaxis.visible = False
    radial_profile_p.legend.visible = False
    radial_profile_p.legend.location = "top_right"
    radial_profile_p.legend.click_policy="hide"
    toggle_legend_js_x = CustomJS(args=dict(leg=radial_profile_p.legend[0]), code="""
        if (leg.visible) {
            leg.visible = false
            }
        else {
            leg.visible = true
        }
    """)
    radial_profile_p.js_on_event(DoubleTap, toggle_legend_js_x)


    @render_widget
    @reactive.event(input.input_type)
    def plot_radial_profile():
        req(input.input_type() in ["Image"])
        return radial_profile_p
    
    @reactive.effect
    @reactive.event(selected_images, input.apix, input.angle, input.dx, input.dy, input.mask_radius, input.negate, input.mask_len)
    def update_radial_profile_data():
        req(input.input_type() in ["Image"])
        req(len(selected_images())>0)
        req(input.apix()!=0)
        if (input.angle() or input.dx() or input.dy() or input.negate()):
            unmasked_data = compute.transform_2d_filament(selected_images()[0],input.angle(),input.dx(),input.dy(),input.negate(),input.apix())
        else:
            unmasked_data = selected_images()[0]
        
        curr_x = np.arange(-nx_curr_img() // 2, nx_curr_img() // 2) * input.apix()
        curr_ymax = np.max(unmasked_data, axis=0)
        curr_ymean = np.mean(unmasked_data, axis=0)

        radial_profile_line_max.data_source.data['x'] = curr_x
        radial_profile_line_max.data_source.data['y'] = curr_ymax        
        radial_profile_line_max_flipped.data_source.data['x'] = -curr_x
        radial_profile_line_max_flipped.data_source.data['y'] = curr_ymax   
        radial_profile_line_mean.data_source.data['x'] = curr_x
        radial_profile_line_mean.data_source.data['y'] = curr_ymean  
        radial_profile_line_mean_flipped.data_source.data['x'] = -curr_x
        radial_profile_line_mean_flipped.data_source.data['y'] = curr_ymean
        radial_profile_rmin_span.location = - input.mask_radius()
        radial_profile_rmax_span.location = input.mask_radius()
    
    #############################
    # acf init and update
    
    init_acf = compute.auto_correlation(init_img_data, sqrt=True, high_pass_fraction=0.1)
    init_apix = 2.3438
    acf_ny = init_ny
    acf_y = np.arange(-acf_ny//2, acf_ny//2)*init_apix
    acf_xmax = np.max(init_acf, axis=1)

    acf_tools = 'box_zoom,crosshair,hover,pan,reset,save,wheel_zoom'
    acf_tooltips = [("Axial Shift", "@y{0.0}Å")]
    acf_p = figure(x_axis_label="Auto-correlation", y_axis_label="Axial Shift (Å)", frame_height=acf_ny,y_range=[acf_ny//2*init_apix, -acf_ny//2*init_apix], sizing_mode="scale_both", tools=acf_tools, tooltips=acf_tooltips)
    acf_line = acf_p.line(acf_xmax, acf_y, line_width=2, color='red', legend_label="ACF")
    print("acf figure created")
    # if 0:
    #     xmean = np.mean(acf, axis=1)
    #     p.line(xmean, y, line_width=2, color='blue', legend_label="mean")
    acf_p.hover[0].attachment = "above"
    #print("acf legends:")
    #print(p.legend.__dir__())
    acf_p.legend.visible = False
    acf_p.legend.location = "top_right"
    acf_p.legend.click_policy="hide"
    toggle_legend_js_y = CustomJS(args=dict(leg=acf_p.legend[0]), code="""
        if (leg.visible) {
            leg.visible = false
            }
        else {
            leg.visible = true
        }
    """)
    acf_p.js_on_event(DoubleTap, toggle_legend_js_y)

    @render_widget
    @reactive.event(input.input_type)
    def plot_acf():
        req(input.input_type() in ["Image"])
        return acf_p
    
    @reactive.effect
    @reactive.event(data_2d_transformed)
    def update_acf_data():
        req(input.input_type() in ["Image"])
        req(data_2d_transformed() is not None)
        print("updating acf")

        ny, nx = data_2d_transformed().shape
        curr_acf = compute.auto_correlation(data_2d_transformed(), sqrt=True, high_pass_fraction=0.1)

        if curr_acf is None or np.isnan(curr_acf).any():
            print("Invalid auto-correlation data!")
            return

        acf_p.y_range.start = ny//2*input.apix()
        acf_p.y_range.end = -ny//2*input.apix()
        acf_p.frame_height = ny
        acf_line.data_source.data['x'] = np.max(curr_acf, axis=1)
        acf_line.data_source.data['y'] = np.arange(-ny//2, ny//2)*input.apix()
    
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

