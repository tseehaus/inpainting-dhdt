import os
import subprocess
import argparse

from Code.inpainting_shearlet import *


SETTINGS = {
    5: {'n_scales': 5, 
        'patch_size': 512, 
        'stride' : 384,
        'downsampling_factor' : 1,
        'apply_mask': True,
       },
    6: {'n_scales': 6, 
        'patch_size': 1024, 
        'stride' : 768,
        'downsampling_factor' : 1,
        'apply_mask': True,
       },
    7: {'n_scales': 6, 
        'patch_size': 2048, 
        'stride' : 1536,
        'downsampling_factor' : 2,
        'apply_mask': False,
       },
}

DEFAULTS = {
    'input_file': os.path.join("center_setup", "input_data", "on_ice_dh_center.tif"),
    'inv_mask_file': os.path.join("center_setup", "input_data", "voids-center_inv.tif"),
    'ice_mask_file': os.path.join("center_setup", "input_data", "ice_mask.tif"),
    'output_folder': os.path.join("outputs"),
    'workers': 1,
    'setting': 5,
}

def check_params(input_file, inv_mask_file, ice_mask_file, output_folder, setting):
    if not os.path.isfile(input_file) or not input_file.endswith('.tif'):
        raise ValueError(f"Invalid input_file: {input_file}\nPlease make sure to pass an existing .tif file")
    if not os.path.isfile(inv_mask_file) or not inv_mask_file.endswith('.tif'):
        raise ValueError("Invalid inv_mask_file,  {inv_mask_file}\nPlease make sure to pass an existing .tif file")
    if not os.path.isfile(ice_mask_file) or not ice_mask_file.endswith('.tif'):
        raise ValueError("Invalid ice_mask_file,  {ice_mask_file}\nPlease make sure to pass an existing .tif file")
    if not os.path.isdir(output_folder):
        raise ValueError("Invalid output_folder, {output_folder}")
    if not setting in SETTINGS.keys():
        raise ValueError(f"Invalid setting: please make sure that it is either of: {list(SETTINGS.keys())}")

if __name__ == "__main__":
    default_output_dir = "inpaint_out/SL/"
    default_input_file = "input_data/on_ice_dh_center_masked.tif"
    default_workers = 1
    default_setting = 5
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file',
                        type=str,
                        default=DEFAULTS['input_file'],
                        help=f"Path to the input geotif-file to be inpainted. \n(default: {DEFAULTS['input_file']})")
    parser.add_argument('--inv_mask_file',
                        type=str,                  
                        default=DEFAULTS['inv_mask_file'],
                        help=f"Path to the inverse mask geotif-file. In the mask file, non voids pixel are encoded as \"nan\"-values. \n(default:{DEFAULTS['inv_mask_file']})")
    parser.add_argument('--ice_mask_file',
                        type=str,
                        default=DEFAULTS['ice_mask_file'],
                        help=f"Path to the ice mask geotif-file. In the ice mask file, ice pixel are encoded as \"nan\"-values.  \n(default:{default_input_file})")
    parser.add_argument('--output_folder',
                        type=str,
                        default=DEFAULTS['output_folder'],
                        help=f"Path to the output folder to store the results. \n(default:{DEFAULTS['output_folder']})")
    parser.add_argument('--workers',
                        type=int,
                        default=DEFAULTS['workers'],
                        help=f"Maximum number of working processes that will be used to processes the patches in parallel. Type: int. \n(default:{DEFAULTS['workers']})")

    parser.add_argument('--setting',
                        type=int,
                        default=DEFAULTS['setting'],
                        help=f"The parameter setting for the inpainting method. Type: int. (5 for SL5, 6 for SL6, 7 for SL7). \n(default:{DEFAULTS['setting']})")
    

    args = parser.parse_args()
    input_file = str(args.input_file)
    inv_mask_file = str(args.inv_mask_file)
    ice_mask_file = str(args.ice_mask_file)
    output_folder = str(args.output_folder)
    max_workers = int(args.workers)
    setting = int(args.setting)
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    check_params(input_file, inv_mask_file, ice_mask_file, output_folder, setting)
    if max_workers<1:
        max_workers=DEFAULTS['workers']
        print(f"Warning: invalid max_workers argument. Default will be set to {DEFAULTS['workers']}")
    
    iterations = 400
    alpha = 0.001
    output_file = os.path.join(output_folder, f'IP_SL_{setting}.tif')
    kwargs = dict(max_workers=max_workers, **SETTINGS[setting])

    perform_patchwise_inpainting(input_file, inv_mask_file, ice_mask_file, output_file, iterations, alpha, **kwargs)
    subprocess.run(f'Rscript Code/geocode_shearlet.R --inpainted_file {output_file} --ice_mask_file {ice_mask_file}', shell=True)
    
    
 

        

    
