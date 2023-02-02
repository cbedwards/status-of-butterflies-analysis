# Read in spatial data ---------------------------------------------------------
library(sp)

# Read in species range polygon
range <- rgdal::readOGR("euphydryas_phaeton.kml")

# Alternatively, use the habitat suitability model instead of the polygon
hsm <- raster::raster("euphydryas_phaeton.tiff")


# Read in species observations
species <- read.csv("cleaned-updated-baltimore checkerspot.csv")

# Assign coordinates to the species object so it become a spatial data frame
coordinates(species) <- ~ lon + lat
sp::proj4string(species) <- sp::proj4string(range)


# Polygon approach -------------------------------------------------------------

raster::plot(species)
raster::plot(range, add=T, border="red", lwd=2)

# Which lat/longs are contained inside the polygon boundaries?
overlaid <- sp::over(species, as(range, "SpatialPolygons"))
within.range <- which(!is.na(overlaid))
outside.range <- which(is.na(overlaid))

# Just for visualization
raster::plot(species, col=c("red", "black")[factor(is.na(overlaid))])

# Then you can drop points based on in/out of range


# Raster approach --------------------------------------------------------------

# The nice thing with the new range map product is that it includes an estimate
# of habitat suitability, so we could also use a more continuous approach for
# deciding which points to keep

site.suitability <- raster::extract(hsm, species)

# What do use for a cutoff?
hist(site.suitability)

x <- raster::extract(hsm, range)
hist(unlist(x))

# Maybe .1? .2?
good.sites <- which(site.suitability>.1)
probably.not <- which(site.suitability<=.1)

raster::plot(species, col=c("red", "black")[factor(site.suitability<=.1)])
