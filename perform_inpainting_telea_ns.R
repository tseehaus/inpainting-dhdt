source(file.path("Code", "inpainting_telea_ns.R"))
library(argparse)

DEFAULTS <- list( input_file = file.path("center_setup", "input_data", "on_ice_dh_center.tif"),
                  mask_file = file.path("center_setup", "input_data", "voids-center_inv.tif"),
                  ice_mask_file = file.path("center_setup", "input_data", "ice_mask.tif"),
                  output_folder = "outputs" )

check_params <- function(input_file, mask_file, ice_mask_file, output_folder) {
    if ( !file.exists(input_file) ||  !endsWith(input_file, '.tif') ) {
        stop(sprintf("Invalid input_file: %s \nPlease make sure to pass an existing .tif file", input_file))
    }
    if ( !file.exists(mask_file) ||  !endsWith(mask_file, '.tif') ) {
        stop(sprintf("Invalid mask_file: %s \nPlease make sure to pass an existing .tif file", mask_file))
    }
    if ( !file.exists(ice_mask_file) ||  !endsWith(ice_mask_file, '.tif') ) {
        stop(sprintf("Invalid ice_mask_file: %s \nPlease make sure to pass an existing .tif file", ice_mask_file))
    }
    if ( !dir.exists(output_folder) ) {
        stop(sprintf("Invalid output folder: %s ", output_folder))
    }
}

parser <- ArgumentParser(description='Inpaint elevation change maps')
parser$add_argument('--output_folder', 
                     help=sprintf("Path to the output folder to store the results. (default : %s)", DEFAULTS$output_folder),
                     default=DEFAULTS$output_folder)
parser$add_argument('--input_file', 
                    help=sprintf("Path to the input geotif-file to be inpainted. (default : %s)",
                                  DEFAULTS$input_file),
                    default=DEFAULTS$input_file)
parser$add_argument('--mask_file',  
                    help=sprintf('Path to the mask geotif-file. In the mask file, voids pixel are encoded as \"nan\"-values. (default : %s)', DEFAULTS$mask_file),
                    default=DEFAULTS$mask_file)
parser$add_argument('--ice_mask_file',  
                     help=sprintf('Path to the ice mask geotif-file. In the ice mask file, ice pixel are encoded as \"nan\"-values.  (default : %s)', DEFAULTS$ice_mask_file),
                    default=DEFAULTS$ice_mask_file)

args <- parser$parse_args()
output_folder <- args$output_folder
input_file <- args$input_file
mask_file <- args$mask_file
ice_mask_file <- args$ice_mask_file

if (!dir.exists(output_folder)) {dir.create(output_folder)}
check_params(input_file, mask_file, ice_mask_file, output_folder)

ip_modes <- c("Telea","NS")
ip_rads <- c(2,5,8,10,15,20)


print("start")
perform_inpainting(output_folder, input_file, mask_file, ice_mask_file, ip_modes, ip_rads)
