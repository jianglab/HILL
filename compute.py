import numpy as np
import pandas as pd
import os, pathlib

from joblib import Memory

username = os.getenv("USER")
memory = Memory(location=f"/tmp/{username}_joblib_cache", verbose=0)

from shiny import reactive
from shiny.express import ui, render, module, expressify


def compute_pair_distances(helices, lengths=None, target_total_count=-1):
    if lengths is not None:
        sorted_indices = (np.argsort(lengths))[::-1]
    else:
        sorted_indices = range(len(helices))
    min_len = 0
    dists_same_class = []
    for i in sorted_indices:
        _, segments_all_classes = helices[i]
        class_ids = np.unique(segments_all_classes["rlnClassNumber"])
        for ci in class_ids:
            mask = segments_all_classes["rlnClassNumber"] == ci
            segments = segments_all_classes.loc[mask, :]
            pos_along_helix = segments["rlnHelicalTrackLengthAngst"].values.astype(
                float
            )
            psi = segments["rlnAnglePsi"].values.astype(float)

            distances = np.abs(pos_along_helix[:, None] - pos_along_helix)
            distances = np.triu(distances)

            # Calculate pairwise distances only for segments with the same polarity
            mask = np.abs((psi[:, None] - psi + 180) % 360 - 180) < 90
            distances = distances[mask]
            dists_same_class.extend(
                distances[distances > 0]
            )  # Exclude zero distances (self-distances)
        if (
            lengths is not None
            and target_total_count > 0
            and len(dists_same_class) > target_total_count
        ):
            min_len = lengths[i]
            break
    if not dists_same_class:
        return [], 0
    else:
        return np.sort(dists_same_class), min_len


def select_helices_by_length(helices, lengths, min_len, max_len):
    min_len = 0 if min_len is None else min_len
    max_len = -1 if min_len is None else max_len

    helices_retained = []
    n_ptcls = 0
    for gi, (gn, g) in enumerate(helices):
        cond = max_len <= 0 and min_len <= lengths[gi]
        cond = cond or (
            max_len > 0 and (max_len > min_len and (min_len <= lengths[gi] < max_len))
        )
        if cond:
            n_ptcls += len(g)
            helices_retained.append((gn, g))
    return helices_retained, n_ptcls


def get_filament_length(helices, particle_box_length=0):
    filement_lengths = []
    for gn, g in helices:
        track_lengths = g["rlnHelicalTrackLengthAngst"].astype(float).values
        length = track_lengths.max() - track_lengths.min() + particle_box_length
        filement_lengths.append(length)
    return filement_lengths


def select_classes(params, class_indices):
    class_indices_tmp = np.array(class_indices) + 1
    mask = params["rlnClassNumber"].astype(int).isin(class_indices_tmp)
    particles = params.loc[mask, :]
    helices = list(particles.groupby(["rlnMicrographName", "rlnHelicalTubeID"]))
    return helices


def get_class_abundance(params, nClass):
    abundance = np.zeros(nClass, dtype=int)
    for gn, g in params.groupby("rlnClassNumber"):
        abundance[int(gn) - 1] = len(g)
    return abundance


def get_number_helices_classes(params):
    nHelices = len(list(params.groupby(["rlnMicrographName", "rlnHelicalTubeID"])))
    nClasses = len(params["rlnClassNumber"].unique())
    return nHelices, nClasses


def get_pixel_size(
    data,
    attrs=[
        "micrograph_blob/psize_A",
        "rlnMicrographPixelSize",
        "rlnMicrographOriginalPixelSize",
        "blob/psize_A",
        "rlnImagePixelSize",
    ],
    return_source=False,
):
    try:
        sources = [data.attrs["optics"]]
    except:
        sources = []
    sources += [data]
    for source in sources:
        for attr in attrs:
            if attr in source:
                if attr in ["rlnImageName", "rlnMicrographName"]:
                    import mrcfile, pathlib

                    folder = pathlib.Path(data["starFile"].iloc[0])
                    if folder.is_symlink():
                        folder = folder.readlink()
                    folder = folder.resolve().parent
                    filename = source[attr].iloc[0].split("@")[-1]
                    filename = str((folder / "../.." / filename).resolve())
                    with mrcfile.open(filename, header_only=True) as mrc:
                        apix = float(mrc.voxel_size.x)
                else:
                    apix = float(source[attr].iloc[0])
                if return_source:
                    return apix, attr
                else:
                    return apix
    return None


def assign_segment_id(data, inter_segment_distance):
    assert "rlnHelicalTrackLengthAngst" in data
    tmp = (
        data.loc[:, "rlnHelicalTrackLengthAngst"].astype(float) / inter_segment_distance
    )
    err = (tmp - tmp.round()).abs()
    if np.sum(np.where(err > 0.1)) > 0:
        print(
            f"WARNING: it appears that the helical segments were extracted with different inter-segment distances"
        )
    helical_segment_id = tmp.round().astype(int)
    return helical_segment_id


def estimate_inter_segment_distance(data):
    # data must have been sorted by micrograph, rlnHelicalTubeID, and rlnHelicalTrackLengthAngst
    helices = data.groupby(["rlnMicrographName", "rlnHelicalTubeID"], sort=False)

    import numpy as np

    dists_all = []
    for _, particles in helices:
        if len(particles) < 2:
            continue
        dists = np.sort(particles["rlnHelicalTrackLengthAngst"].astype(float).values)
        dists = dists[1:] - dists[:-1]
        dists_all.append(dists)
    dists_all = np.hstack(dists_all)
    dist_seg = np.median(dists_all)  # Angstrom
    return dist_seg


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
    data = get_class2d_from_file(fileobj.name)
    print("data: ", data)
    return data


def get_class2d_from_file(classFile):
    import mrcfile
    print("before mrcfile")
    #classFile = "C:\\Users\\anika\\Downloads\\run_it020_classes.mrcs"
    try: 
        open(classFile, "r")
    except Exception as e:
        print("open exception: ", e)

    print("after opening")
    with mrcfile.open(classFile) as mrc:
        print("opens mrc file")
        apix = float(mrc.voxel_size.x)
        data = mrc.data
    print("got data: ", data)
    return data, round(apix, 4)


@memory.cache
def get_class2d_params_from_url(url, url_cs_pass_through=None):
    url_final = get_direct_url(url)  # convert cloud drive indirect url to direct url
    #print("1 download file failed", url_final)
    fileobj = download_file_from_url(url_final)
    #print("2 file object", fileobj.name)
    if fileobj is None:
        raise ValueError(
            f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview"
        )
    
    if url_cs_pass_through is None:
        #print("url pass through is none")
        data = get_class2d_params_from_file(fileobj.name)
        #print("data: ", data)
        return data
    

    url_final_cs_pass_through = get_direct_url(
        url_cs_pass_through
    )  # convert cloud drive indirect url to direct url
    
    fileobj_cs_pass_through = download_file_from_url(url_final_cs_pass_through)
    if fileobj_cs_pass_through is None:
        raise ValueError(
            f"ERROR: {url_cs_pass_through} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview"
        )
    
    #print("fileobj: ", fileobj.name)
    data = get_class2d_params_from_file(fileobj.name, fileobj_cs_pass_through.name)
    return data


def get_class2d_params_from_file(params_file, cryosparc_pass_through_file=None):
    #print("start function")
    if params_file.endswith(".star"):
        #print("before params")
        
        #print("file obj name: ", params_file)
        #os.listdir("C:\\Users\\anika\\AppData\\Local\\Temp\\")
        params = star_to_dataframe(params_file)
        #print("end params, ", params)
    elif params_file.endswith(".cs"):
        assert cryosparc_pass_through_file is not None
        params = cs_to_dataframe(params_file, cryosparc_pass_through_file)
    required_attrs = np.unique(
        "rlnImageName rlnHelicalTubeID rlnHelicalTrackLengthAngst rlnClassNumber rlnAnglePsi".split()
    )
    missing_attrs = [attr for attr in required_attrs if attr not in params]
    if missing_attrs:
        raise ValueError(f"ERROR: parameters {missing_attrs} are not available")
    return params


def star_to_dataframe(starFile):
    import starfile
    import os 
    #print("file obj name: ", starFile)
    #temp = os.listdir("C:\\Users\\anika\\AppData\\Local\\Temp\\")
    #print(temp)
    d = starfile.read(starFile, always_dict=True)
    assert (
        "optics" in d and "particles" in d
    ), f"ERROR: {starFile} has {' '.join(d.keys())} but optics and particles are expected"
    data = d["particles"]
    data.attrs["optics"] = d["optics"]
    data.attrs["starFile"] = starFile
    return data


def cs_to_dataframe(cs_file, cs_pass_through_file):
    cs = np.load(cs_file)
    df_cs = pd.DataFrame.from_records(cs.tolist(), columns=cs.dtype.names)
    cs_passthrough = np.load(cs_pass_through_file)
    df_cs_passthrough = pd.DataFrame.from_records(
        cs_passthrough.tolist(), columns=cs_passthrough.dtype.names
    )
    data = pd.concat([df_cs, df_cs_passthrough], axis=1)
    data = data.loc[:, ~data.columns.duplicated()]
    # rlnImageName rlnHelicalTubeID rlnHelicalTrackLengthAngst rlnCoordinateX rlnCoordinateY rlnClassNumber rlnAnglePsi
    ret = pd.DataFrame()
    if "blob/idx" in data and "blob/path" in data:
        ret["rlnImageName"] = (
            (data["blob/idx"].astype(int) + 1).map("{:06d}".format)
            + "@"
            + data["blob/path"].str.decode("utf-8")
        )
    if "blob/psize_A" in data:
        ret["rlnImagePixelSize"] = data["blob/psize_A"]
        ret["blob/psize_A"] = data["blob/psize_A"]
    if "micrograph_blob/path" in data:
        ret["rlnMicrographName"] = data["micrograph_blob/path"]
    if "micrograph_blob/psize_A" in data:
        ret["rlnMicrographPixelSize"] = data["micrograph_blob/psize_A"]
        ret["micrograph_blob/psize_A"] = data["micrograph_blob/psize_A"]
    if "location/micrograph_path" in data:
        ret["rlnMicrographName"] = data["location/micrograph_path"]
    if (
        "location/center_x_frac" in data
        and "location/center_y_frac" in data
        and "location/micrograph_shape" in data
    ):
        locations = pd.DataFrame(data["location/micrograph_shape"].tolist())
        my = locations.iloc[:, 0]
        mx = locations.iloc[:, 1]
        ret["rlnCoordinateX"] = (
            (data["location/center_x_frac"] * mx).astype(float).round(2)
        )
        ret["rlnCoordinateY"] = (
            (data["location/center_y_frac"] * my).astype(float).round(2)
        )
    if "filament/filament_uid" in data:
        if "blob/path" in data:
            if data["filament/filament_uid"].min() > 1000:
                micrographs = data.groupby(["blob/path"])
                for _, m in micrographs:
                    mapping = {
                        v: i + 1
                        for i, v in enumerate(
                            sorted(m["filament/filament_uid"].unique())
                        )
                    }
                    ret.loc[m.index, "rlnHelicalTubeID"] = m[
                        "filament/filament_uid"
                    ].map(mapping)
            else:
                ret.loc[:, "rlnHelicalTubeID"] = data["filament/filament_uid"].astype(
                    int
                )

            if "filament/position_A" in data:
                filaments = data.groupby(["blob/path", "filament/filament_uid"])
                for _, f in filaments:
                    val = f["filament/position_A"].astype(np.float32).values
                    val -= np.min(val)
                    ret.loc[f.index, "rlnHelicalTrackLengthAngst"] = val.round(2)
        else:
            mapping = {
                v: i + 1
                for i, v in enumerate(sorted(data["filament/filament_uid"].unique()))
            }
            ret.loc[:, "rlnHelicalTubeID"] = data["filament/filament_uid"].map(mapping)
    if "filament/filament_pose" in data:
        ret.loc[:, "rlnAnglePsi"] = np.round(
            -np.rad2deg(data["filament/filament_pose"]), 1
        )
    # 2D class assignments
    if "alignments2D/class" in data:
        ret["rlnClassNumber"] = data["alignments2D/class"].astype(int) + 1
    if "alignments2D/shift" in data:
        shifts = pd.DataFrame(data["alignments2D/shift"].tolist()).round(2)
        ret["rlnOriginX"] = -shifts.iloc[:, 0]
        ret["rlnOriginY"] = -shifts.iloc[:, 1]
    if "alignments2D/pose" in data:
        ret["rlnAnglePsi"] = -np.rad2deg(data["alignments2D/pose"]).round(2)
    return ret


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


def plot_histogram(
    data,
    title,
    xlabel,
    ylabel,
    max_pair_dist=None,
    bins=50,
    log_y=True,
    show_pitch_twist={},
    multi_crosshair=False,
    fig=None,
):
    import plotly.graph_objects as go

    if max_pair_dist is not None and max_pair_dist > 0:
        data = [d for d in data if d <= max_pair_dist]

    hist, edges = np.histogram(data, bins=bins)
    hist_linear = hist
    if log_y:
        hist = np.log10(1 + hist)

    center = (edges[:-1] + edges[1:]) / 2

    hover_text = []
    for i, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        hover_info = f"{xlabel.replace(' (Å)', '')}: {center[i]:.0f} ({left:.0f}-{right:.0f})Å<br>{ylabel}: {hist_linear[i]}"
        if show_pitch_twist:
            rise = show_pitch_twist["rise"]
            csyms = show_pitch_twist["csyms"]
            for csym in csyms:
                twist = 360 / (center[i] * csym / rise)
                hover_info += f"<br>Twist for C{csym}: {twist:.2f}°"
        hover_text.append(hover_info)

    if fig:
        fig.data[0].x = center
        fig.data[0].y = hist
        fig.data[0].text = hover_text
        fig.layout.title.text = title
    else:
        fig = go.FigureWidget()

        histogram = go.Bar(
            x=center,
            y=hist,
            name="Histogram",
            marker_color="blue",
            hoverinfo="none",
        )

        fig.add_trace(histogram)

        fig.data[0].text = hover_text
        fig.data[0].hoverinfo = "text"
        fig.update_layout(
            template="plotly_white",
            title_text=title,
            title_x=0.5,
            title_font=dict(size=12),
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            autosize=True,
            hovermode="closest",
            hoverlabel=dict(bgcolor="white", font_size=12),
        )

        if multi_crosshair:
            for i in range(20):
                fig.add_vline(
                    x=0,
                    line_width=3 if i == 0 else 2,
                    line_dash="solid" if i == 0 else "dash",
                    line_color="green",
                    visible=False,
                )

            def update_vline(trace, points, state):
                if points.point_inds:
                    hover_x = points.xs[0]
                    with fig.batch_update():
                        for i, vline in enumerate(fig.layout.shapes):
                            x = hover_x * (i + 1)
                            vline.x0 = x
                            vline.x1 = x
                            if x <= fig.data[0].x.max():
                                vline.visible = True
                            else:
                                vline.visible = False

            fig.data[0].on_hover(update_vline)

    return fig




import pathlib
import helicon

def extract_emdb_id(url):
    import re
    pattern = r'EMD-(\d+)'
    match = re.search(pattern, url)
    if match:
        return f"EMD-{match.group(1)}"
    return None

class MapInfo:
    def __init__(self, data=None, filename=None, url=None, emd_id=None, label="", apix=None, twist=None, rise=None, csym=1):
        non_nones = [p for p in [data, filename, url, emd_id] if p is not None]
        if len(non_nones)>1:
            raise ValueError(f"MapInfo(): only one of these parameters can be set: data, filename, url, emd_id")
        elif len(non_nones)<1:
            raise ValueError(f"MapInfo(): one of these parameters must be set: data, filename, url, emd_id")
        self.data = data
        self.filename = filename
        self.url = url
        self.emd_id = emd_id
        self.label = label
        self.apix = apix
        self.twist = twist
        self.rise = rise
        self.csym = csym

    def __repr__(self):
        return (f"MapInfo(label={self.label}, emd_id={self.emd_id}, "
                f"twist={self.twist}, rise={self.rise}, csym={self.csym}, "
                f"apix={self.apix})")
        
    def get_data(self):
        if self.data is not None:
            return self.data, self.apix
        if isinstance(self.filename, str) and len(self.filename) and pathlib.Path(self.filename).exists():
            self.data, self.apix = get_images_from_file(self.filename)
            return self.data, self.apix
        if isinstance(self.url, str) and len(self.url):
            url = extract_url(self.url)
            print("url: ", url)
            self.data, self.apix = get_images_from_url(url) # self.url)
            return self.data, self.apix
        if isinstance(self.emd_id, str) and len(self.emd_id):
            emdb = helicon.dataset.EMDB()
            self.data, self.apix = emdb(self.emd_id)
            return self.data, self.apix
        raise ValueError(f"MapInfo.get_data(): failed to obtain data")


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





