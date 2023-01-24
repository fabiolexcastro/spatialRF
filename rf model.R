


# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(rnaturalearthdata, rnaturalearth, SSDM, raster, terra, sf, fs, glue, tidyverse, rgbif, geodata)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
spce <- 'Acrocomia mexicana' # Coyol Alimenticio. Los frutos llamados “coyol” son comestibles. Medicinal. Es utilizada en medicina tradicional en forma de extractos acuosos para el tratamiento de diabetes mellitus.
occr <- occ_data(scientificName = spce, limit = 2e5, hasCoordinate = T, hasGeospatialIssue = F)

# Selecting the dataframe from the object
occr <- occr[[2]]
colnames(occr)
unique(occr$country) %>% sort()
occr <- filter(occr, country == 'Mexico')

# Shapefiles
wrld <- ne_countries(returnclass = 'sf', scale = 50)

plot(st_geometry(wrld))
plot(mex1)
points(occr$decimalLongitude, occr$decimalLatitude, pch = 16, col = 'red')

mex1 <- geodata::gadm(country = 'MEX', level = 1, path = 'tmpr')
bioc <- geodata::worldclim_country(country = 'MEX', var = 'bioc', path = 'tmpr')
bioc <- terra::crop(bioc, mex1) %>% terra::mask(., mex1)
names(bioc) <- glue('bioc{1:19}')

occr <- dplyr::select(occr, x = decimalLongitude, y = decimalLatitude)
vles <- terra::extract(bioc, occr[,c('x', 'y')])
occr <- cbind(occr[,c('x', 'y')], vles[,-1])
occr <- as_tibble(occr)
occr <- mutate(occr, pb = 1)

# To generate background --------------------------------------------------
cell <- terra::extract(bioc[[1]], occr[,1:2], cells = T)$cell
duplicated(cell)
mask <- bioc[[1]] * 0
mask[cell] <- NA
back <- terra::as.data.frame(mask, xy = T) %>% as_tibble()
back <- sample_n(back, size = nrow(occr) * 2, replace = FALSE)
colnames(back)[3] <- 'pb'
back <- mutate(back, pb = 0)
back <- cbind(back, terra::extract(bioc, back[,c(1, 2)])[,-1])
back <- as_tibble(back)

# Join  -------------------------------------------------------------------
tble <- rbind(occr, back)

# Random forest -----------------------------------------------------------
bioc <- stack(bioc)
tble <- as.data.frame(tble)
sdrf <- modelling(algorithm = 'RF', Env = bioc, Occurrences = tble, Pcol = 'pb', Xcol = 'x', cv.parm = c(0.75, 0.25), Ycol = 'y', metric = 'TSS', select.metric = 'AUC')

plot(sdrf@projection)
plot(sdrf@binary)
sdrf@parameters

sdrf@name
sdrf@variable.importance
as.numeric(sdrf@variable.importance) %>% sum()

rstr <- sdrf@projection
rstr <- terra::rast(rstr)
rslt <- terra::as.data.frame(rstr, xy = T) %>% as_tibble()

# To make the map ---------------------------------------------------------


