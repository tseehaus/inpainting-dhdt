source("Code", "stat_analysis_global.R"))
library(argparse)


DEFAULTS <- list( inpaint_folder = file.path("correlation_setup", "inpaint_out"),
                  input_file = file.path("correlation_setup", "on_ice_ifsar_srtm_2012_2013_dh.tif"),
                  mask_file = file.path("correlation_setup", "corr-masks-inv", "50_corr_mask_inv.tif"),
                  ice_mask_file = file.path("correlation_setup", "ice_mask.tif"),
                  output_folder = "outputs"
                )


check_inputs <- function(inpaint_folder, input_file, mask_file, ice_mask_file){
    if (!dir.exists(inpaint_folder)) {stop(paste0("invalid inpainting directory: ", inpaint_folder))}
    if (!file.exists(input_file) || !endsWith(input_file, ".tif")) {stop(paste0("invalid input file: ", input_file))}
    if (!file.exists(mask_file) || !endsWith(mask_file, ".tif")) {stop(paste0("invalid input file: ", mask_file))}
    if (!file.exists(ice_mask_file) || !endsWith(ice_mask_file, ".tif")) {stop(paste0("invalid mask file: ", ice_mask_file))}
}

parser <- ArgumentParser(description='Inpaint elevation change maps')
parser$add_argument('--inpaint_folder', 
                    help=sprintf("Directory where the inpainted files are stored. (default : %s)",
                                  DEFAULTS$inpaint_folder),
                     default=DEFAULTS$inpaint_folder)
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

args <- parser$parse_args()
output_folder <- args$output_folder
inpaint_folder <- args$inpaint_folder
input_file <- args$input_file
mask_file <- args$mask_file
ice_mask_file <- args$ice_mask_file

check_inputs(inpaint_folder, input_file, mask_file, ice_mask_file)

if(!dir.exists(output_folder)){ dir.create(output_folder) }

csv_file_residuals <- file.path(output_folder, "residuals.csv")
csv_file_stats <- file.path(output_folder, "table3.csv")
png_file_plot <- file.path(output_folder, "figure6.png")

print("perform statistical analysis (global)")
perform_statistical_analysis_global(inpaint_folder, input_file, mask_file, csv_file_residuals, csv_file_stats, png_file_plot)
