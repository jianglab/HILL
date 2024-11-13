import numpy as np
from shiny import App, ui, reactive, render
import shiny_test
import compute

data_all = reactive.value(None)
image_size = reactive.value(0)
displayed_class_images = reactive.value([])
displayed_class_labels = reactive.value([])
initial_selected_image_indices = reactive.value([])
selected_images = reactive.value([])
selected_image_labels = reactive.value([])

default_url = "https://ftp.ebi.ac.uk/empiar/world_availability/10940/data/EMPIAR/Class2D/job010/run_it020_classes.mrcs"


app_ui = ui.page_fluid(
    ui.head_content(ui.tags.title("Class2D Viewer")),
    ui.tags.style(
        """
        * { font-size: 10pt; padding:0; border: 0; margin: 0; }
        aside {--_padding-icon: 10px;}
        """
    ),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_radio_buttons(
                "input_mode_classes",
                "How to obtain the class average images:",
                choices=["upload", "url"],
                selected="upload",
                inline=True,
            ),
            ui.panel_conditional(
                "input.input_mode_classes == 'upload'",
                ui.input_file(
                    "upload_classes",
                    "Upload the class averages in MRC format (.mrcs, .mrc)",
                    accept=[".mrcs", ".mrc"],
                    placeholder="mrcs or mrc file",
                )
            ),
            ui.panel_conditional(
                "input.input_mode_classes == 'url'",
                ui.input_text(
                    "url_classes",
                    "Download URL for MRC file",
                    value=default_url,
                )
            ),
            ui.input_action_button("run", label="Run", style="width: 100%;"),
            ui.output_ui("dynamic_image_select"),
        ),
        ui.column(
            12,
            ui.h1("Class2D Viewer", style="font-weight: bold;"),
            ui.output_ui("display_selected_images"),
            ui.accordion(
                ui.accordion_panel(
                    "Additional Parameters",
                    ui.input_checkbox("ignore_blank", "Ignore blank classes", value=True),
                    value="additional_parameters"
                ),
            )
        )
    ),
)

def server(input, output, session):
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
                enable_selection=True,
                allow_multiple_selection=True,
            )
        return ui.div("No images selected")
    
    @output
    @render.ui
    def dynamic_image_select():
        if len(displayed_class_images()) > 0:
            return shiny_test.image_gallery(
                id="select_classes",
                label="Select class(es):",
                images=displayed_class_images(),
                image_labels=displayed_class_labels(),
                image_size=128,
                initial_selected_indices=initial_selected_image_indices(),
                enable_selection=True,
                allow_multiple_selection=True,
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
        if input.input_mode_classes() == "upload":
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
        
        elif input.input_mode_classes() == "url" and input.url_classes():
            url = input.url_classes()
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

app = App(app_ui, server)