require(raster)
library(argparse)



temp_dir <-(".temp")
if (!dir.exists(temp_dir)) {dir.create(temp_dir)}
rasterOptions(tmpdir=temp_dir, datatype = 'FLT4S',overwrite = T)

parser <- ArgumentParser(description='Add geocoding and ice mask to output of shearlet inpainting')
parser$add_argument('--inpainted_file', 
                     help='the output of the shearlet inpainting')
parser$add_argument('--ice_mask_file',  
                    help='the ice mask, i.e., a geotif file that specifies the ice locations (every NA pixel represents ice)')

args <- parser$parse_args()
ice_mask_file <- args$ice_mask_file
inpainted_file <- args$inpainted_file


ice_mask <- raster(ice_mask_file)
inpainted <- raster(inpainted_file)

crs(inpainted) <- crs(ice_mask)
extent(inpainted) <- extent(ice_mask)

inpainted <- inpainted*ice_mask

writeRaster(inpainted,inpainted_file,overwrite=T)


