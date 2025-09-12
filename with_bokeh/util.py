from typing import Optional
from pathlib import Path

import shiny
from shiny import reactive
from shiny import ui, module, render
from shiny.express import expressify

import socket
from getpass import getuser
from psutil import virtual_memory, Process
from os import getpid

import helicon
import os
import tempfile
import numpy as np
import mrcfile
import pandas as pd

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

def get_2d_image_from_uploaded_file(fileobj):
    #import os, tempfile
    #original_filename = fileobj.name
    # input_data = mrcfile.open(fileobj, permissive=True).data
    # map_crs_auto = None
    # apix_auto = None

    # if input_data is not None:
    #     apix_auto = 1.0  # or whatever default value you want to use
    #     map_crs_auto = input_data.copy()  # if you need a copy of the data
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

    # return input_data, map_crs_auto, apix_auto

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
    # print("emdb helical: ", emdb_helical)
    # print("emdb_ids_helical: ", emdb_ids_helical)

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
        @reactive.effect
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
    url_final = helicon.get_direct_url(url)    # convert cloud drive indirect url to direct url
    fileobj = download_file_from_url(url_final)
    if fileobj is None:
        #st.error(f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview")
        #st.stop()
        @reactive.effect
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
        filesize = helicon.get_file_size(url)
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

# def change_mrc_map_crs_order(data, current_order, target_order=[1, 2, 3]):
#     if current_order == target_order: return data
#     map_crs_to_np_axes = {1:2, 2:1, 3:0}
#     current_np_axes_order = [map_crs_to_np_axes[int(i)] for i in current_order]
#     target_np_axes_order = [map_crs_to_np_axes[int(i)] for i in target_order]
#     import numpy as np
#     ret = np.moveaxis(data, current_np_axes_order, target_np_axes_order)
#     return ret

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

int_types = {'apply_helical_sym_0':0, 'apply_helical_sym_1':0, 'csym':1, 'csym_ahs_0':1, 'csym_ahs_1':1, 'do_random_embid_0':0, 'do_random_embid_1':0, 'fft_top_only':0, 'image_index_0':0, 'image_index_1':0, 'input_mode_0':1, 'input_mode_1':1, 'is_3d_0':0, 'is_3d_1':0, 'm_0':1, 'm_1':1, 'm_max':3, 'negate_0':0, 'negate_1':0, 'pnx':512, 'pny':1024, 'show_LL':1, 'show_LL_text':1, 'show_phase_diff':1, 'show_pwr':1, 'show_yprofile':1, 'transpose_0':0, 'transpose_1':0, 'share_url':0, 'show_qr':0, 'useplotsize':0}
float_types = {'angle_0':0, 'angle_1':0, 'apix_0':0, 'apix_1':0, 'apix_ahs_0':0, 'apix_ahs_1':0, 'apix_map_0':0, 'apix_map_1':0, 'apix_nyquist_0':0, 'apix_nyquist_1':0, 'az_0':0, 'az_1':0, 'ball_radius':0, 'cutoff_res_x':0, 'res_limit_y':0, 'diameter':0, 'dx_0':0, 'dx_1':0, 'dy_0':0, 'dy_1':0, 'fraction_ahs_0':0, 'fraction_ahs_1':0, 'length_ahs_0':0, 'length_ahs_1':0, 'mask_len_0':90, 'mask_len_1':90, 'mask_radius_0':0, 'mask_radius_1':0, 'noise_0':0, 'noise_1':0, 'resolution':0, 'rise':0, 'rise_ahs_0':0, 'rise_ahs_1':0, 'simuaz':0, 'simunoise':0, 'tilt':0, 'tilt_0':0, 'tilt_1':0, 'twist':0, 'twist_ahs_0':0, 'twist_ahs_1':0, 'width_ahs_0':0, 'width_ahs_1':1}
other_types = {'const_image_color':'', 'emd_id_0':'', 'emd_id_1':'', 'input_type_0':'image', 'input_type_1':'image', 'll_colors':'lime cyan violet salmon silver', 'url_0':'', 'url_0':''}


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

# def get_file_size(url):
#     import requests
#     response = requests.head(url)
#     if 'Content-Length' in response.headers:
#         file_size = int(response.headers['Content-Length'])
#         return file_size
#     else:
#         return None
