source(file.path("Code", "stat_analysis_local.R"))
library(argparse)

DEFAULTS <- list( inpaint_folder = file.path("center_setup", "inpaint_out"),
                  shape_folder = file.path("center_setup", "void_shapes"),
                  input_file = file.path("center_setup", "input_data", "on_ice_dh_center.tif"),
                  mask_file = file.path("center_setup", "input_data", "voids-center_inv.tif"),
                  ice_mask_file = file.path("center_setup", "input_data", "ice_mask.tif"),
                  output_folder = "outputs",
                  setup = "center" )

check_inputs <- function(shape_folder, inpaint_folder, input_file, mask_file, ice_mask_file, setup){
    if (!dir.exists(shape_folder)) {stop(paste0("invalid shape directory: ", shape_folder))}
    if (!dir.exists(inpaint_folder)) {stop(paste0("invalid inpainting directory: ", inpaint_folder))}
    if (!file.exists(input_file) || !endsWith(input_file, ".tif")) {stop(paste0("invalid input file: ", input_file))}
    if (!file.exists(mask_file) || !endsWith(mask_file, ".tif")) {stop(paste0("invalid mask file: ", mask_file))}
    if (!file.exists(ice_mask_file) || !endsWith(ice_mask_file, ".tif")) {stop(paste0("invalid ice mask file: ", ice_mask_file))}
    if (!(setup=="center") && !(setup=="juneau")) {stop(paste0("invalid setup: ", setup, ". expected: center or juneau"))}
}

parser <- ArgumentParser(description='Inpaint elevation change maps')
parser$add_argument('--inpaint_folder', 
                    help=sprintf("Directory where theinpainted files are stored. (default : %s)",
                                  DEFAULTS$inpaint_folder),
                     default=DEFAULTS$inpaint_folder)
parser$add_argument('--shape_folder', 
                     help=sprintf("Directory where the shape files are stored. (default : %s)",
                                  DEFAULTS$shape_folder),
                     default=DEFAULTS$shape_folder)
parser$add_argument('--input_file',  
                    help=sprintf("Path to the input geotif-file. (default : %s)",
                                  DEFAULTS$input_file),
                    default=DEFAULTS$input_file)
parser$add_argument('--mask_file',  
                    help=sprintf('Path to the mask geotif-file. In the mask file, voids pixel are encoded as \"nan\"-values. (default : %s)', DEFAULTS$mask_file),
                    default=DEFAULTS$mask_file)
parser$add_argument('--ice_mask_file',  
                     help=sprintf('Path to the ice mask geotif-file. In the ice mask file, ice pixel are encoded as \"nan\"-values.  (default : %s)', DEFAULTS$ice_mask_file),
                    default=DEFAULTS$ice_mask_file)
parser$add_argument('--output_folder', 
                     help=sprintf("Path to the output folder to store the results. (default : %s)", DEFAULTS$output_folder),
                     default=DEFAULTS$output_folder)
parser$add_argument('--setup', 
                     help=sprintf("Setup either 'center' or 'juneau'. (default : %s)", DEFAULTS$setup),
                     default=DEFAULTS$setup)

args <- parser$parse_args()
inpaint_folder <- args$inpaint_folder
shape_folder <- args$shape_folder
setup <- tolower(args$setup)
input_file <- args$input_file
mask_file <- args$mask_file
ice_mask_file <- args$ice_mask_file
output_folder <- args$output_folder


check_inputs(shape_folder, inpaint_folder, input_file, mask_file, ice_mask_file, setup)
shape_files <- get_shape_files(shape_folder, setup)
if (length(shape_files)==0) { stop(paste0("could not find shapefiles in: ", shape_dir)) }

if(!dir.exists(output_folder)){ dir.create(output_folder) }

print("perform statistical analysis (local)")
perform_statistical_analysis(input_file, mask_file, ice_mask_file, inpaint_folder, output_folder, shape_files, setup)
