require(raster)
require(rgdal)

temp_dir <-(".temp")
if (!dir.exists(temp_dir)) {dir.create(temp_dir)}
rasterOptions(tmpdir=temp_dir, datatype = 'FLT4S',overwrite = T)

load_rasters_files <- function(dem_file, mask_file, ice_mask_file=FALSE){
    dem <- raster(dem_file)
    mask_ <- raster(mask_file)
    if (ice_mask_file != FALSE) {
        ice_mask <- raster(ice_mask_file)
        rasters <- list(dem=dem, mask=mask_, ice_mask=ice_mask)
    } else {
        rasters <- list(dem=dem, mask=mask_)
    } 
    return(rasters)
}

sort_inpainted_files <-function(inpainted_dem_files, inpaint_type, inpaint_subtype, type_sort_order){
    sorting_indices <- c()
    for (z in 1:length(type_sort_order)){ ## sort accordint to "typ_order"
        idx<-which(inpaint_type %in% type_sort_order[z])
        sorting_indices <- c(sorting_indices, idx)
        }
    inpainted_dem_files <- inpainted_dem_files[sorting_indices]
    inpaint_subtype <- inpaint_subtype[sorting_indices]
    inpaint_type <- inpaint_type[sorting_indices]
    inpainted <- list(inpainted_dem_files=inpainted_dem_files, inpaint_subtype=inpaint_subtype, inpaint_type=inpaint_type)
    return(inpainted)
}

extract_inpainted_files <- function(inpaint_dir){
    inpainted_dem_files <-dir(inpaint_dir, pattern="Correlation.tif$",full.names = T,recursive = T) 
    inpaint_type <- matrix(unlist(strsplit(basename(inpainted_dem_files),split="_")),nrow=length(inpainted_dem_files),byrow=T)[,2] 
    inpaint_subtype <- matrix(unlist(strsplit(basename(inpainted_dem_files),split="_")),nrow=length(inpainted_dem_files),byrow=T)[,3]

    type_sort_order<-c("NS","Telea","SL","CL")
    inpainted <- sort_inpainted_files(inpainted_dem_files, inpaint_type, inpaint_subtype, type_sort_order)
    return(inpainted)
}

sample_elevation_differences <- function(inpainted_dem_files, dem, mask, corr_len){
    residuals <- list()
    for (i in 1:length(inpainted_dem_files)){
        # cat ("Gen. stats of ",basename(inpainted_dem_files[i]),"\n")
        dem_inpainted <- raster(inpainted_dem_files[i])
        residual  <- dem_inpainted *mask - dem # elevation difference on inpainted areas

        residual<-sampleRegular(residual, ncell(dem)/(corr_len[i]^2), asRaster=F) # regular sampling considering the correlation length
        residual<-na.omit(residual) # exclude NA values
        residuals[[i]]<-residual # add to list
    } 
    max_length <- max(unlist(lapply(residuals,length))) 
    residuals <- data.frame(lapply(residuals,"length<-",max_length)) # convert list to data frame
    return(residuals)
}

write_residuals_to_csv <- function(residuals, inpaint_type, inpaint_subtype, csv_file){
    colnames(residuals)<-paste(c("dh_IP(m)"), inpaint_type, inpaint_subtype, sep="_") # set col. names
    write.csv(residuals, csv_file, row.names = F, na="", quote=F) 
    }

calc_stats <- function(residuals, alpha=0.05){
    mean_offsets <- unlist(lapply(residuals,mean,na.rm=T)) 
    std_devs <- unlist(lapply(residuals, sd,na.rm=T)) 
    num_samples <- unlist(lapply(residuals,length)) 
    conf_interval_mean_offsets <- matrix(,ncol=2,nrow=length(residuals))
    conf_interval_std_devs <- matrix(,ncol=2,nrow=length(residuals))
    for (i in 1:length(residuals)){
      ttest<-t.test(residuals[[i]], conf.level=1-alpha/18^2)
      conf_interval_mean_offsets[i,]<-ttest$conf.int  
    }
    for (i in 1:length(residuals)){
        sd_up  <- sqrt((num_samples[i]-1)*std_devs[i]^2/qchisq(((alpha/18^2)/2), num_samples[i]-1, lower.tail = T))
        sd_low <- sqrt((num_samples[i]-1)*std_devs[i]^2/qchisq(((alpha/18^2)/2), num_samples[i]-1, lower.tail = F))
        conf_interval_std_devs[i,] <- cbind(sd_low,sd_up)
    }
    stats <- list(mean_offsets=mean_offsets, 
                  std_devs=std_devs, 
                  conf_interval_mean_offsets=conf_interval_mean_offsets, 
                  conf_interval_std_devs=conf_interval_std_devs,
                  num_samples=num_samples)
    return(stats)
}

create_ci_plots <- function(png_file, mean_offsets, conf_interval_mean_offsets, std_devs, conf_interval_std_devs, inpaint_type, inpaint_subtype){
    png(filename=png_file,width=16, height=10, units="in", res=300)
    p_colors<-c(rep("red",6),rep("blue",6),rep("orange",3),rep("black",3)) # plot colors
    op = par(mfrow = c(2,1),cex=1)
    yticks <- rev(1:length(p_colors))
    ytick_labels <- paste0(inpaint_type,"-",inpaint_subtype)
    
    par(mar=c(4,5,4,0.5)+0.1)
    plot(mean_offsets,yticks,xlim=c(min(conf_interval_mean_offsets[,1]),max(conf_interval_mean_offsets[,2])), yaxt="n",xaxt="n",ylab="",xlab=expression(bar(dh)    ("m")),pch=16,col=p_colors) #plot mean
    arrows(conf_interval_mean_offsets[,1], yticks, conf_interval_mean_offsets[,2], yticks, length=0.05, angle=90, code=3,col=p_colors) # conf. intervals
    abline(v=0)
    axis(1,at=seq(-1,1,0.05),labels=round(seq(-1,1,0.05),2),las=2,las=1)  #x axis labels
    axis(2,at=yticks, labels=ytick_labels, las=2) # yaxis labels

    par(mar=c(4,5,1,0.5)+0.1)
    plot(std_devs,yticks,xlim=c(0,max(conf_interval_std_devs)), yaxt="n",xaxt="n",ylab="",xlab=expression(sigma[dh]),pch=16,col=p_colors) 
    arrows(conf_interval_std_devs[,1], yticks, conf_interval_std_devs[,2], yticks, length=0.05, angle=90, code=3,col = p_colors) # pot conf.   interval
    abline(v=0)
    axis(1,at=seq(0,15,1),labels=round(seq(0,15,1),0),las=2,las=1) # x axis labels
    axis(2,at=yticks, labels=ytick_labels, las=2) # y axis labels
    dev.off()
    }
    
prepare_table_data <- function(inpaint_type, inpaint_subtype, corr_len, stats){
    table_data <- cbind(inpaint_type, 
                       inpaint_subtype,
                       round(stats$mean_offsets, 4),
                       round(stats$conf_interval_mean_offsets, 4),
                       round(stats$std_devs, 4),
                       round(stats$conf_interval_std_devs, 4),
                       corr_len)
    colnames(table_data)<-c("IP_mode",
                            "IP_para",
                            "mean",
                            "mean_conf_min",
                            "mean_conf_max",
                            "SD",
                            "SD_conf_min",
                            "SD_conf_max",
                            "cor_length")
    return(table_data)   
}

write_stats_to_csv <- function(csv_file, inpaint_type, inpaint_subtype, corr_len, stats){
    table_data <- prepare_table_data(inpaint_type, inpaint_subtype, corr_len, stats)
    write.csv(table_data, csv_file, row.names = F, na="",quote=F)
}


perform_statistical_analysis_global <- function(inpaint_dir, dem_file, mask_file, csv_file_residuals, csv_file_stats, png_file_plot){
    print("Loading files ...")
    rasters <- load_rasters_files(dem_file, mask_file)    
    inpainted <- extract_inpainted_files(inpaint_dir)
    correlation_lengths <- c(5,5,5,5,5,5, # NS
                             5,6,6,6,8,8, # Telea
                             4,4,5, # SL
                             4,10,10) # classical
    print("Sample elevation differences ...")
    residuals <- sample_elevation_differences(inpainted$inpainted_dem_files, rasters$dem, rasters$mask, correlation_lengths)
    write_residuals_to_csv(residuals, inpainted$inpaint_type, inpainted$inpaint_subtype, csv_file_residuals)
    print("Calculating statistics ...")
    stats <- calc_stats(residuals)
    write_stats_to_csv(csv_file_stats, inpainted$inpaint_type, inpainted$inpaint_subtype, correlation_lengths, stats)
    print("Creating plots ...")
    create_ci_plots(png_file_plot, stats$mean_offsets, stats$conf_interval_mean_offsets, stats$std_devs, stats$conf_interval_std_devs, inpainted$inpaint_type, inpainted$inpaint_subtype)
    print("Finished")
    }