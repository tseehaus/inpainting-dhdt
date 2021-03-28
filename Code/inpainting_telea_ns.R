require(raster)
require(rgeos)
library(Rvision)
library("stringr")

temp_dir <-(".temp")
if (!dir.exists(temp_dir)) {dir.create(temp_dir)}
rasterOptions(tmpdir=temp_dir, datatype = 'FLT4S',overwrite = T)


load_input_rasters <- function(input_raster_file, mask_file, ice_mask_file) {
    l = list()
    dem <- raster(input_raster_file)
    mask <- raster(mask_file)
    ice_mask <- raster(ice_mask_file)
    l[c("dem", "mask", "ice_mask")] = list(dem, mask, ice_mask)
    return (l)
}

apply_mask <- function (dem, mask) {
    cat ("masking DEMs using masks \n")
    l = list()
    all_gl<-na.omit(getValues(dem))
    dem_m<-mask(dem,mask,inverse=T)
    voids<-length(na.omit(getValues(dem_m)))
    voids_perc<-1-voids/(length(all_gl))
    l[c("all_gl", "dem_m", "voids", "voids_perc")] = list(all_gl, dem_m, voids, voids_perc)
    return(l)
}

log_masked_image <- function(input_file, dem_m) {
    l = list()
    root_path <- dirname(input_file)
    filename <- str_replace( basename(input_file), ".tif", "")
    dem_m_file <- paste0( file.path(root_path, filename), "_masked.tif")
    writeRaster(dem_m, dem_m_file)
}

rescale_to_int16_range <- function (dem_m) {
    l = list()
    qs<-c(-200,200) # remove strong outliers
    dem_m[dem_m>=qs[2]]<-qs[2]
    dem_m[dem_m<=qs[1]]<-qs[1]
    sc=(2^16)/(qs[2]-qs[1]) # scaling to 16 bit
    dem_m_sc<-(dem_m-qs[1])*sc 
    min<-qs[1]
    mIP<-dem_m
    mIP[!is.na(mIP)]<-0 ### create mask with voids are "white/255"
    mIP[is.na(mIP)]<-255
    l[c("dem_m", "dem_m_sc", "mIP", "min", "sc")] = list(dem_m, dem_m_sc, mIP, min, sc)
    return(l)
}
        
log_and_load_as_image <- function (input_file, dem_m_sc, mIP){
    l = list()
    root_path <- dirname(input_file)
    filename <- str_replace( basename(input_file), ".tif", "")
    dem_m_sc_file <- paste0( file.path(root_path, filename), "_masked_SC.tif")
    mIP_file <- file.path(root_path, "void_IP_mask.tif")
    
    writeRaster(dem_m_sc, dem_m_sc_file, datatype="INT2U" )
    writeRaster(mIP, mIP_file, datatype="INT1U" )

    # load dem as image object
    l[c("dem_m_sc_im", "m_im")] <- list( image(dem_m_sc_file),  image(mIP_file) )
    return(l)
}
    
log_stats <- function(all_gl, voids, voids_perc, sc, min){
    write.table(cbind( (length(all_gl)),
                        voids_perc,
                        voids*30*30/1000000,
                        sc,
                        min),
                "dem_masked_void_stats.txt",
                append=F, 
                quote=F,
                col.names = c("n-pixels","voids_%","Area_km²","scale","min_val"),
                row.names = F)
}

apply_inpainting_methods <- function(inpaint_dir, dem_m_sc_im, m_im, ice_mask, dem_m_sc, sc, min_sc, ip_mode, ip_rad){
    cat("starting In-painting \n")
    j<-1
    for (j in 1:length(ip_mode)){ # loop IP modes
        folder = file.path(inpaint_dir, ip_mode[j])
        if (! file_test('-d', folder)){
            dir.create(folder, recursive=T)
        }
        k<-1
        for (k in 1:length(ip_rad)){                        # loop IP radius
            cat ("inpainting using ", ip_mode[j], " and radius:", ip_rad[k],"\n")
            
            in_im<-inpaint(dem_m_sc_im, m_im, ip_rad[k], ip_mode[j])
            in_im_m<-as.matrix(in_im) # convert to matrix
            in_im_mr<-in_im_m[c(nrow(in_im_m):1),] #flip matix
            values(dem_m_sc)<-in_im_mr  # paste inpainted values in masked raster
            dem_m_sc<-dem_m_sc*ice_mask # mask to glacier areas
            dem_m<-dem_m_sc/sc+min_sc

            filename = file.path(folder, sprintf("IP_%s_%02d", ip_mode[j], ip_rad[k]) )
            png_file = paste0(filename, ".png")
            tif_file = paste0(filename, ".tif")
            
            write.Image(in_im, png_file)
            writeRaster(dem_m, tif_file)
        }
    }
}

perform_inpainting  <- function(output_dir, input_file, mask_file, ice_mask_file, ip_modes, ip_rads) {
    l = list()
    l[c("dem", "mask", "ice_mask")] <- load_input_rasters(input_file, mask_file, ice_mask_file)
    l[c("all_gl", "dem_m", "voids", "voids_perc")] <- apply_mask(l$dem, l$mask)
    log_masked_image(input_file, l$dem_m)
    l[c("dem_m", "dem_m_sc", "mIP", "min", "sc")] <- rescale_to_int16_range(l$dem_m)
    l[c("dem_m_sc_im", "m_im")] <- log_and_load_as_image(input_file, l$dem_m_sc, l$mIP)
    log_stats(l$all_gl, l$voids, l$voids_perc, l$sc, l$min)
    apply_inpainting_methods(output_dir, l$dem_m_sc_im, l$m_im, l$ice_mask, l$dem_m_sc, l$sc, l$min, ip_modes, ip_rads)
}
