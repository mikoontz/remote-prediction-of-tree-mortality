# This script: 
# 1) Loads the monthly climate data summaries Derek produced ("features/working-files/climate_data_matrix_jepson_PPN+SMC_central+south.Rdata")
# 2) Computes annual summaries (climate indices) from these. 
# 3) cbinds these to the EVI data matrix 

load("features/working-files/climate_data_matrix_jepson_PPN+SMC_central+south.Rdata")

temp_mat <- as.data.frame(clim_mat[,207:410])
ppt_mat <- as.data.frame(clim_mat[,3:206])
#temp_mat$month <- ppt_mat$month <- month <- rep(1:12, 17)
#temp_mat$year <- ppt_mat$year <- rep(2000:2016, rep(12, 17))
#temp_mat$season <- ppt_mat$season <- ifelse(month%in%c(12,1,2),"Winter", ifelse(month%in%c(3,4,5),"Spring", ifelse(month%in%c(6,7,8), "Summer",ifelse(month%in%c(9,10,11),"Fall", NA))))

# Calculate indices
years <- 2000:2016
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
for (i in 1:16) {
  names(ppt_wy)[i] <- paste("ppt_wy_", i+2000, sep="")
  names(mean_temp)[i] <- paste("mean_temp_", i+2000, sep="")
  names(winter_temp)[i] <- paste("winter_temp_", i+2000, sep="")
  names(summer_temp)[i] <- paste("summer_temp_", i+2000, sep="")
}

clim_data <- cbind(clim_mat[,c("x", "y")], ppt_wy, mean_temp, winter_temp, summer_temp)
head(clim_data)

coldata = clim_data$ppt_wy_2001
cols = round(coldata-min(coldata))
palette(viridis(max(cols)))
plot(y~x, clim_data, col = cols)
     