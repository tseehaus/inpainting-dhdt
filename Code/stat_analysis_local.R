require(rgdal)
require(raster)
require(rgeos)
require(tools)

temp_dir <-(".temp")
if (!dir.exists(temp_dir)) {dir.create(temp_dir)}
rasterOptions(tmpdir=temp_dir, datatype = 'FLT4S',overwrite = T)

load_input_rasters <- function(dem_file, mask_file, ice_mask_file){
    dem <- raster(dem_file)
    mask_ <- raster(mask_file)
    ice_mask <- raster(ice_mask_file)
    l <- list(dem=dem, mask=mask_, ice_mask=ice_mask)
    return( l )
}

load_shapes <- function(shape_files){
    shapes <- list()
    shape_names <- names(shape_files)
    for (i in seq_along(shape_files) ){
        shape = shapefile(shape_files[[i]])
        shapes[shape_names[i]] = shape
        }
    return( shapes )
}

get_ip_types  <- function(IP_file){
    split_filename <- unlist(strsplit(file_path_sans_ext(basename(IP_file)), split="_"))
    df <- data.frame("IP_typ" = split_filename[2], 
                     "IP_styp" = split_filename[3])
    return( df )
}

prepare_data <- function(IP_dem, dem, filtered_dem, kernel, mask_, custom_mask){
    IP_dem_masked <- IP_dem * mask_ 
    offset <- IP_dem_masked - dem 
    relative_offset <- offset / dem * 100 # relative offset (%)
    
    offset <- mask(offset, custom_mask)
    relative_offset <- mask(relative_offset, custom_mask)

    filtered_offset <- focal(offset, w=kernel, na.rm=T, pad=T, padValue=NA) 
    filtered_offset <- mask(filtered_offset, offset)
    filtered_dem <- mask(filtered_dem, relative_offset)
    dem <- mask(dem, relative_offset)
    relative_filtered_offset <- filtered_offset / filtered_dem * 100
    
    offset <- na.omit(getValues(offset))
    relative_offset <- na.omit(getValues(relative_offset))
    filtered_offset <- na.omit(getValues(filtered_offset))
    relative_filtered_offset <- na.omit(getValues(relative_filtered_offset))
    dem <- na.omit(getValues(dem))
    filtered_dem <- na.omit(getValues(filtered_dem))
    
    l = list(offset=offset, relative_offset=relative_offset, filtered_offset=filtered_offset, 
             relative_filtered_offset=relative_filtered_offset, dem=dem, filtered_dem= filtered_dem)
    return( l )
}

calculate_stats <- function(offset, relative_offset, offset_filtered, relative_filtered_offset,  dem, dem_filtered, res){
    stats = data.frame(
        "area" = length(offset)*res^2/10^6,
        "mean_offset" = mean(offset),
        "median_offset" = median(offset),
        "min_offset" = min(offset),
        "max_offset" = max(offset),
        "RMSE" = sqrt(sum(offset^2)/length(offset)),
        "SD_offset" = sd(offset),
        "Vol_offset" = sum(offset)*res^2/10^9,
        "mad_offset" = mad(offset),
        "AAE" = mean(abs(offset)), 
        "AME" = median(abs(offset)),
        "ASD" = sd(abs(offset)), 
        "relative_AAE" = mean(abs(relative_offset)),
        "relative_AME" = median(abs(relative_offset)),
        "relative_ASD" = sd(abs(relative_offset)),
        "mean_org" = mean(dem),
        "mean_org_filtered" = mean(dem_filtered),
        "mean_offset_filtered" = mean(offset_filtered), 
        "AAE_filtered" = mean(abs(offset_filtered)), 
        "AME_filtered" = median(abs(offset_filtered)),
        "ASD_filtered" = sd(abs(offset_filtered)), 
        "relative_AAE_filtered" = mean(abs(relative_filtered_offset)), 
        "relative_AME_filtered" = median(abs(relative_filtered_offset)), 
        "relative_ASD_filtered" = sd(abs(relative_filtered_offset))
    )
    return( stats )
}

calculate_stats_table <- function(IP_dems, dem, filtered_dem, kernel, mask_, custom_mask, res){
    tab <- data.frame()
    for (j in 1:length(IP_dems)){
        #cat ("Gen. stats (IP area only) of ",basename(IP_dems[j]),"\n")
        IP_file <- IP_dems[j]
        IP_dem <- raster(IP_file)
        l <- prepare_data(IP_dem, dem, filtered_dem, kernel, mask_, custom_mask) 
        stats <- calculate_stats(l$offset, l$relative_offset, l$filtered_offset, l$relative_filtered_offset,  l$dem, l$filtered_dem, res)
        ip_types <- get_ip_types(IP_file)
        row <- merge(ip_types, stats)
        tab <- rbind(tab, row)
    }
    return( tab )
}

create_tables <- function(dem_file, mask_file, ice_mask_file, shape_files, inpaint_dir, output_dir, setup){
    rasters <- load_input_rasters(dem_file, mask_file, ice_mask_file)
    IP_dems<-dir(inpaint_dir, pattern=".tif$",full.names = T,recursive = T)
    resolution<-mean(res(rasters$dem))
    
    kernel <- focalWeight(rasters$dem, 300, "Gauss") 
    filtered_dem <- focal(rasters$dem, w=kernel,na.rm=T, pad=T, padValue=NA) * rasters$ice_mask
    shape_names <- names(shape_files)
    tables = list()
    for (i in 1:length(shape_files)){
        shape_mask <- shapefile(shape_files[[i]])
        custom_mask <- spTransform(shape_mask, CRS=proj4string(rasters$dem))
        tab <- calculate_stats_table(IP_dems, rasters$dem, filtered_dem, kernel, rasters$mask, custom_mask, resolution)
        tab_file <- file.path(output_dir, paste0("Stats_summary_masked_", setup, "_", shape_names[i], ".csv") )
        write.csv(tab, tab_file, row.names=F,na="")
        tables <-append(tables, list(tab))
    }
    return (tables)
}

prepare_plot_data <- function(data){
    colors <- c("red", "blue", "orange", "black")
    itypes <-c ("NS", "Telea", "SL", "CL")

    data$color <- character(nrow(data))
    for (i in seq_along(itypes)){
        data$color[which(data$IP_typ == itypes[i])] <- colors[i]
    }
    data$mean_p <- data$mean_offset / data$mean_org*100

    print("Successfull")
    return (data)
}

create_plot <- function(data, column_name, title, title_pos, pch, text){
    print("Create Plot")

    vmin <- min(0, min(data[[column_name]]))
    vmax <- max(0, max(data[[column_name]]))
    
    plot(x=data[[column_name]],
         y=1:nrow(data),
         names.arg=data$IP_styp,
         yaxt="n",
         xlab="",
         ylab="",
         space=0.1,
         col=data$color,
         pch=pch,
         cex=2,
         xlim=c(vmin, vmax),
         main="") 
    
    mtext(text, side=1, line=2.5, cex=2.5)
    legend(title_pos, title, text.font=2,cex=1.2)
    abline(v=0)
    axis(2,at=1:nrow(data), labels=paste0(data$IP_typ, "-", data$IP_styp),las=2)
    print("Successfull")
}

create_plots <- function(data_list, output_dir){
    png(filename=file.path(output_dir,"Plot_mdh_r_AAE_f.png"),width=30, height=20, units="in", res=300)
    op = par(mfrow = c(3,2), cex=2, mar=c(4,5,0,5)+0.1)
    plot_titles<-c("Center - Circle", "Center - Strip", "Center - Terminus")    
    for (i in seq_along(plot_titles)){
        data <- data_list[[i]]
        data <- prepare_plot_data(data)
        legend_pos = "topleft"
        create_plot(data, "mean_p", plot_titles[i], legend_pos, pch=19, text=expression(bar(dh^r)("%")) )
        #create_plot(data, "mean_offset", plot_titles[i], legend_pos, pch=19, text=expression(bar(dh^r)("%")) )           
        if (i!=1){
            legend_pos = "topright"
        }
        create_plot(data, "AAE_filtered", plot_titles[i], legend_pos, pch=15, text=expression(AAE^{f}*(m)) )
    }
    dev.off()
}

get_shape_files <- function(shape_dir, setup){
    if (setup == "center") {
        shape_files <- c(Circular = file.path(shape_dir, "circ-crop-center.shp"),
                            Strip = file.path(shape_dir, "strip-crop-center.shp"),
                            Terminus = file.path(shape_dir, "term-crop-center.shp"),
                            All = file.path(shape_dir, "all-crop-center.shp") )
    }
    else if (setup == "juneau"){
        shape_files <- c(Circular = file.path(shape_dir, "circ-crop-juneau.shp"),
                            Strip = file.path(shape_dir, "strip-crop-juneau.shp"),
                            All = file.path(shape_dir, "all-crop-juneau.shp") )
    }
    else {
        return(c())
    }
    indices = c()
    for (i in length(shape_files):1){
        if (file.exists(shape_files[[i]]) ) {
            indices <- c(indices, i)
        }
    }
    shape_files = c(shape_files[indices])
    return(shape_files)
}

perform_statistical_analysis <- function(dem_file, mask_file, ice_mask_file, inpaint_dir, output_dir, shape_files, setup){
    tables <- create_tables(dem_file, mask_file, ice_mask_file, shape_files, inpaint_dir, output_dir, setup)
    create_plots(tables, output_dir)
    
}
