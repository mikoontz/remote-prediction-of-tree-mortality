# This script: 
# 1) Loads the monthly climate data summaries Derek produced ("features/working-files/climate_data_matrix_jepson_PPN+SMC_central+south.Rdata")
# 2) Computes annual summaries (climate indices) from these. 
# 3) cbinds these to the EVI data matrix 

library(viridis)
library(raster)

load("features/working-files/climate_data_matrix_jepson_PPN+SMC_central+south.Rdata")

temp_mat <- as.data.frame(clim_mat[,219:434])
ppt_mat <- as.data.frame(clim_mat[,3:218])

# Calculate simple indices
# Annual precip (water year), and summer mean temp, winter mean temp, and annual mean temp
years <- 1999:2016
n.years <- length(years)
n.cells <- nrow(clim_mat)
ppt_wy = summer_temp = winter_temp = mean_temp = {}
for (i in 2:n.years) {
    ppt_sub <- cbind(ppt_mat[,grep(years[i-1], names(ppt_mat))][10:12], ppt_mat[,grep(years[i], names(ppt_mat))][1:9])
    tmp_sub <- cbind(temp_mat[,grep(years[i-1], names(temp_mat))][12],temp_mat[,grep(years[i], names(temp_mat))][1:11]) 
    ppt_wy <- cbind(ppt_wy, apply(ppt_sub, 1, sum))
    mean_temp <- cbind(mean_temp, apply(temp_mat, 1, mean))
    summer_temp <- cbind(summer_temp, apply(temp_mat[, 7:9], 1, mean))
    winter_temp <- cbind(winter_temp, apply(temp_mat[, 1:3], 1, mean))
}
ppt_wy = as.data.frame(ppt_wy)
mean_temp = as.data.frame(mean_temp)
winter_temp = as.data.frame(winter_temp)
summer_temp = as.data.frame(summer_temp)
for (i in 1:17) {
  names(ppt_wy)[i] <- paste("ppt_wy_", i+1999, sep="")
  names(mean_temp)[i] <- paste("mean_temp_", i+1999, sep="")
  names(winter_temp)[i] <- paste("winter_temp_", i+1999, sep="")
  names(summer_temp)[i] <- paste("summer_temp_", i+1999, sep="")
}

clim_data <- cbind(clim_mat[,c("x", "y")], ppt_wy, mean_temp, winter_temp, summer_temp)
clim_data$cell_number <- rownames(clim_data)
head(clim_data)


# Save to working files folder
write.csv(clim_data, "features/working-files/climate_data_summaries_jepson_PPN+SMC_central+south.csv")

# Merge with EVI summary data for exploratory check
load("features/working-files/evi_summary_PPN+SMC_jepson_central+south.Rdata")
evi_template = raster("features/sierra-nevada-250m-evi-template.tif")

evi_summary_clim <-  merge(evi_summary, clim_data,  by="cell_number")
head(evi_summary_clim)
cor(evi_summary_clim[,c(2:16, 27, 44, 61, 78)])
heatmap(cor(evi_summary_clim[,c(2:16, 27, 44, 61, 78)]))

# Quick EDA using 2008 climate data 
cols_to_scale <- c("evi_mean", "mean_temp_2008", "ppt_wy_2008", "wet_dry_diff")
for (i in 1:length(cols_to_scale)) evi_summary_clim[,cols_to_scale[i]] <- scale(evi_summary_clim[,cols_to_scale[i]])
summary(lm(sqrt(mort)~evi_mean + mean_temp_2008, data=evi_summary_clim))
summary(lm(sqrt(mort)~evi_mean * mean_temp_2008 * ppt_wy_2008 * wet_dry_diff, data=evi_summary_clim))

# Check by plotting
clim_sp = SpatialPointsDataFrame(evi_summary_clim, coords=evi_summary_clim[,c("x", "y")], proj4string = CRS(projection(evi_template)))
coldata = clim_sp$ppt_wy_2001
cols = round(coldata-min(coldata))
palette(viridis(max(cols)))
plot(clim_sp, col = cols)

coldata = clim_sp$summer_temp_2000
cols = round(coldata-min(coldata))
palette(heat.colors(max(cols)))
plot(clim_sp, col = cols)

coldata = clim_sp$wet_dry_diff
cols = round(coldata-min(coldata))
palette(viridis(max(cols)))
plot(clim_sp, col = cols)

m <- lm(sqrt(mort)~evi_mean * summer_temp_2008 * ppt_wy_2008 * wet_dry_diff, data=evi_summary_clim)
plot(clim_sp$mort~predict(m), pch=".")
