library(Rhabit)
library(terra)
library(ggplot2)
library(viridis)
library(gridExtra)





hbfull <- brick(here("case_study/aleut_habitat.grd"), values = TRUE)

covlist0 <- list(
  bathy  = hbfull$bathy,
  slope  = hbfull$slope,
  d2site = hbfull$d2site,
  d2shelf = hbfull$d2shelf
)

# Convert raster extent/projection from m to km
for (i in seq_along(covlist0)) {
  extent(covlist0[[i]]) <- extent(
    c(
      xmin(covlist0[[i]]), xmax(covlist0[[i]]),
      ymin(covlist0[[i]]), ymax(covlist0[[i]])
    ) / 1000
  )
  
  projection(covlist0[[i]]) <- gsub(
    "units=m", "units=km", projection(covlist0[[i]])
  )
}

# Resample all covariates to the same grid
ncov <- length(covlist0)
for (i in 2:ncov) {
  covlist0[[i]] <- resample(covlist0[[i]], covlist0[[1]])
}

# Crop covariates to area of interest
border <- 30
lim <- c(min(xy[,1])-border,max(xy[,1])+border,min(xy[,2])-border,max(xy[,2])+border)
covlist0 <- lapply(covlist0, crop, y=extent(lim))

covlist <- lapply(covlist0, raster_to_rhabit)







covnames <- c("Bathymetry", "Slope", "Distance to SSL sites", "Distance to shelf")
covplot <- list()
for(i in 1:ncov) {
  mycov <- covlist0[[i]]
  covmap <- data.frame(coordinates(mycov),val=values(mycov))
  covplot[[i]] <- ggplot(covmap,aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=bquote(paste("c"[.(i)]))) +
    ggtitle(covnames[i]) + xlab("Easting (km)") + ylab("Northing (km)") +
    theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
          legend.title = element_text(size=15), legend.text = element_text(size=12),
          legend.key.height=unit(2,"line"), title = element_text(size=12))
}
covplot[[ncov+1]] <- ncov/2 # "ncol" in grid.arrange call
names(covplot)[ncov+1] <- "ncol"

pdf(file = "SSLcovs.pdf", width=10, height=6)
do.call("grid.arrange",covplot)
dev.off()



# Plot locations
ggplot(tracks, aes(x/1000, y/1000, col=factor(ID), group=ID)) + geom_point(size=0.5) + geom_path() +
  coord_equal(xlim=lim[1:2], ylim=lim[3:4]) + xlab("Easting (km)") + ylab("Northing (km)") +
  scale_color_manual(values=c("firebrick","seagreen","royalblue"), name="ID")


do.call(grid.arrange, c(covplot[1:4], ncol = 2))

















