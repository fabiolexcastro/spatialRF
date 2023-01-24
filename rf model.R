

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(rnaturalearthdata, rnaturalearth, cptcity, SSDM, ggspatial, raster, terra, sf, fs, glue, tidyverse, rgbif, geodata)

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
windowsFonts(georg = windowsFont('Georgia'))

gmap <- ggplot() + 
  geom_tile(data = rslt, aes(x = x, y = y, fill = Projection)) +
  scale_fill_gradientn(colors = cpt(pal = 'imagej_gyr_centre', n = 10, rev = TRUE)) +
  geom_sf(data = wrld, fill = NA, col = 'grey40', lwd = 0.2) + 
  geom_sf(data = st_as_sf(mex1), fill = NA, col = 'grey40', lwd = 0.3) + 
  coord_sf(xlim = ext(mex1)[1:2], ylim = ext(mex1)[3:4]) + 
  labs(x = 'Lon', y = 'Lat', fill = 'Puntaje de idoneidad') + 
  ggtitle(label = 'Idoneidad para la especie Coyol en el país de México', subtitle = 'Modelo Random Forest') +
  theme_bw() + 
  theme(text = element_text(family = 'georg', color = 'grey50'), 
        legend.position = 'bottom', 
        plot.title = element_text(hjust = 0.5, face = 'bold', color = 'grey30'),
        plot.subtitle = element_text(hjust = 0.5, face = 'bold', color = 'grey30'),
        # legend.key.width = unit(3, 'line'),
        panel.border = element_rect(color = 'grey80')) +
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) + 
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'georg', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'georg', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 

ggsave(plot = gmap, filename = 'png/mapa_rf.png', units = 'in', width = 9, height = 7, dpi = 300)



