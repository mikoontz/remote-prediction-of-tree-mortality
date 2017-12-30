# This script: 
# 1) Loads the monthly climate data summaries Derek produced ("features/working-files/climate_data_matrix_jepson_PPN+SMC_central+south.Rdata")
# 2) Computes annual summaries (climate indices) from these. 
# 3) cbinds these to the EVI data matrix 

load("features/working-files/climate_data_matrix_jepson_PPN+SMC_central+south.Rdata")

temp_mat <- as.data.frame(t(clim_mat[,207:410]))
ppt_mat <- as.data.frame(t(clim_mat[,3:206]))
temp_mat$month <- ppt_mat$month <- month <- rep(1:12, 17)
temp_mat$year <- ppt_mat$year <- rep(2000:2016, rep(12, 17))
temp_mat$season <- ppt_mat$season <- ifelse(month%in%c(12,1,2),"Winter", ifelse(month%in%c(3,4,5),"Spring", ifelse(month%in%c(6,7,8), "Summer",ifelse(month%in%c(9,10,11),"Fall", NA))))

# Calculate indices
years <- 2000:2016
n.years <- length(years)
n.cells <- nrow(clim_mat)
ppt_wy <- matrix(0, nrow=n.cells, ncol=n.years)
for (i in 1:n.cells) {
    ppt_sub <- matrix(ppt_mat[,i], ncol=12, byrow=T)
    ppt_wy <- apply(ppt_sub, 2, f<-function(x){return(sum())})
  }
}
