library(rase)
library(ape)
library(rgeos)
library(rgdal)
library(maptools)
library(raster)
library(ks)
library(coda)
library(tmap)
library(RColorBrewer)
library(sp)
library(sf)
library(dplyr)


#########################################
#              Starting Vars            #
#########################################

setwd("~/Documentos/Omar/Projects/Amazonas_2020/2020/raseDatosRAmazonia/R/")

sp <- "melipona"
dist <- "log"
rep <- 50000

taxon <- paste0(sp,"-",dist,"-prunned",rep,".tre",sep="")
taxon2 <- paste0(sp,"/",sp,"AllSp",sep="")

rutaTrees <- "../data/trees/"
rutaDist <- "../data/shp/"

arbol <- read.tree( paste0(rutaTrees,taxon,sep="") )
datosS <- paste0(rutaDist,taxon2,".shp",sep="")

a <- rgdal::readOGR(datosS)

########################################
#               RASE p1                #
########################################

distributionData <- shape.to.rase(a)
#distributionData <- distributionData[-which(a$Sp%in%paste0(arbol$tip.label,".shp")==F)] #Solo stenodermatinae & cebidae
poly <- name.poly(distributionData, arbol,  poly.names = arbol$tip.label)

rasemcmc <- readRDS("../results/rase/melipona/melipona-log50000rase.RDS")


#######################################
#               Maps                  #
#######################################

coast <- readOGR("/media/omar/HD710 PRO/Omar/QGIS/CoastLine/ne_50m_coastline.shp")
CP <-   as(extent(-90, -30, -35, 14), "SpatialPolygons")

proj4string(CP) <- CRS(proj4string(coast))
co <-  gIntersection(coast, CP, byid=TRUE)

ras <- raster::raster("/media/omar/HD710 PRO/Omar/QGIS/alt/wc2.1_5m_elev.tif")
proj4string(CP) <- CRS(proj4string(co))
ra <- crop(ras, CP)
ra <- mask(ra, mask = CP)
proj4string(ra) <- CRS(proj4string(co))

pebas17 <- readOGR("../../shp/Pebas_Miocene_17Ma.shp")
pebas16 <- readOGR("../../shp/Pebas_Miocene_16Ma.shp")
pebas13 <- readOGR("../../shp/Pebas_Miocene_13-14Ma.shp")
pebas12 <- readOGR("../../shp/Pebas_Miocene_12Ma.shp")
pebas11 <- readOGR("../../shp/Pebas_Miocene_11Ma.shp")
pebas10 <- readOGR("../../shp/Pebas_Miocene_10Ma.shp")
amazonas8 <- readOGR("../../shp/Amazonas_A_8Ma.shp")
amazonas5 <- readOGR("../../shp/Amazonas_A_5Ma.shp")
andes33 <- crop(ra,CP) %>% mask(mask = readOGR("../../shp/Andes_60-33_Ma.shp"))
andes23 <- crop(ra,CP) %>% mask(mask = readOGR("../../shp/Andes_33-23_Ma.shp"))
andes10 <- crop(ra,CP) %>% mask(mask = readOGR("../../shp/Andes_23-10_Ma.shp"))
andes7 <- crop(ra,CP) %>% mask(mask = readOGR("../../shp/Andes_10-7_Ma.shp"))
andes3 <- crop(ra,CP) %>% mask(mask = readOGR("../../shp/Andes_7-3_Ma.shp"))
andes <- crop(ra,CP) %>% mask(mask = readOGR("../../shp/Andes_3-1_Ma.shp"))

an33 <- as(readOGR("../../shp/Andes_60-33_Ma.shp"),"SpatialLines")
an23 <- as(readOGR("../../shp/Andes_33-23_Ma.shp"),"SpatialLines")
an10 <- as(readOGR("../../shp/Andes_23-10_Ma.shp"),"SpatialLines")
an7 <- as(readOGR("../../shp/Andes_10-7_Ma.shp"),"SpatialLines")
an3 <- as(readOGR("../../shp/Andes_7-3_Ma.shp"),"SpatialLines")
an <- as(readOGR("../../shp/Andes_3-1_Ma.shp"),"SpatialLines")

pixs = 1e5

addalpha = function(colors, alpha=1.0) {
  r = col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] = alpha*255
  r = r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}


col.pals <-  c('BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys',
               'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral','Accent','Dark2','Paired','Pastel1','Pastel2','Set1','Set2','Set3','Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','Oranges', 'BuPu', 'PuRd', 'Greys')

palRas <- brewer.pal(9,"Greys")

nombre <- taxon


########################################
#           Slice of time              #
########################################

valor <- max(adephylo::distRoot(arbol)) - 1
limite <- floor(valor)
corte  <- 1
slice <- 1

while (corte < limite) {
  corte <- corte + slice
  slice_results <- rase.slice(arbol, slice = corte, res = rasemcmc, poly, niter = 1000)
  rsl <- slice_results

#plot(ra, col = brewer.pal(9, 'Greys'), maxpixels = pixs, legend = FALSE, las = 1, cex.axis = 1.4, col.axis = 'grey50', main = paste0(nombre,".",corte-1," My"))
#plot(co, col = 'grey', lwd = 1.5, add = TRUE)


###################################################################


  if((corte-1)%in%c(63:33)){
    map <- tm_shape(andes33)+
             tm_raster(palette = palRas, legend.show = F, alpha = .5)+
           tm_shape(an33)+
             tm_lines(col="white")
  }  
  if((corte-1)%in%c(32:23)){
    map <- tm_shape(andes23)+
             tm_raster(palette = palRas, legend.show = F, alpha = .5)+
           tm_shape(an23)+
             tm_lines(col="white")
  }
  if((corte-1)%in%c(22:10)){
    map <- tm_shape(andes10)+
             tm_raster(palette = palRas, legend.show = F, alpha = .5)+
           tm_shape(an10)+
             tm_lines(col="white")
  }
  if((corte-1)%in%c(9:7)){
    map <- tm_shape(andes7)+
             tm_raster(palette = palRas, legend.show = F, alpha = .5)+
           tm_shape(an7)+
             tm_lines(col="white")
  }
  if((corte-1)%in%c(6:3)){
    map <- tm_shape(andes3)+
             tm_raster(palette = palRas, legend.show = F, alpha = .5)+
           tm_shape(an3)+
             tm_lines(col="white")
  }
  if((corte-1)%in%c(2:1)){
    map <- tm_shape(andes)+
             tm_raster(palette = palRas, legend.show = F, alpha = .5)+
           tm_shape(an)+
             tm_lines(col="white")
  }
         
  map <- map + tm_grid(n.x = 5, n.y = 5, alpha = 0.25)+
         tm_shape(co)+
            tm_lines()+
         tm_layout(legend.position = c("right", "top"), 
                   title= paste0(corte-1," ","Ma"), 
                   title.position = c('right', 'top'),title.size = 3)

  if((corte-1) %in% c(17,15)){
    map <- map+tm_shape(pebas17)+
               tm_polygons(col="lightblue",border.col = NULL, alpha = 0.5)
  }
  
  if((corte-1)==16){
    map <- map+tm_shape(pebas16)+
              tm_polygons(col="lightblue",border.col = NULL, alpha = 0.5)
  }
  
  if((corte-1)%in%c(14,13)){
    map <- map+tm_shape(pebas13)+
               tm_polygons(col="lightblue",border.col = NULL, alpha = 0.5)
  }
  
  if((corte-1)==12){
    map <- map+tm_shape(pebas12)+
               tm_polygons(col="lightblue",border.col = NULL, alpha = 0.5)
  }

  if((corte-1)==11){
    map <- map+tm_shape(pebas11)+
              tm_polygons(col="lightblue",border.col = NULL, alpha = 0.5)
  }

  if((corte-1)%in%c(10,9)){
    map <- map+tm_shape(pebas10)+
               tm_polygons(col="lightblue",border.col = NULL, alpha = 0.5)
 }

  if((corte-1)%in%c(8,7,6)){
    map <- map+tm_shape(amazonas8)+
               tm_lines(col="skyblue",border.col = NULL, alpha = 0.75)
  }
  
  if((corte-1)%in%c(5:1)){
     map <- map+tm_shape(amazonas5)+
                tm_lines(col="skyblue",border.col = NULL, alpha = 0.75)
  }

  ccL <- list()

  for (i in 1:(ncol(rsl)/2)) {
  
    # make `x` & `y` object
    df = data.frame(rsl[,i], rsl[,i+(ncol(rsl)/2)])
  
    # estimate Highest Posterior Interval
    hh = Hpi(df, binned = TRUE)*1
  
    # make 2D dimensional smoothing
    dd = kde(df, H = hh)
  
    # create contour lines for polygons
    cc = contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]],
                      z = dd$estimate, levels = seq(0, max(dd$estimate),
                                                    length = 20))
    # set color
    col1 <- addalpha(brewer.pal(3, sample(col.pals,1) ), 0.4)
  
    # plot three stacked polygons for each level of the HPI
    #polygon(cc[[16]]$x, cc[[16]]$y, col = col1[1], border = 'grey60', lwd = 0.3)
    #polygon(cc[[17]]$x, cc[[17]]$y, col = col1[2], border = 'grey60', lwd = 0.3)
    #polygon(cc[[18]]$x, cc[[18]]$y, col = col1[3], border = 'grey60', lwd = 0.3)

    cc2 <- st_as_sf(ContourLines2SLDF(list(cc[[length(cc)-1]],
                                           cc[[length(cc)-3]],
                                           cc[[length(cc)-4]]) ) ) %>% st_polygonize()
    
    if(length(cc2$level)==2){cc2 <- st_as_sf(ContourLines2SLDF(list(cc[[length(cc)-1]],
                                                              cc[[length(cc)-sample(1:3,1)]],
                                                              cc[[length(cc)-sample(4:6,1)]]) ) ) %>% st_polygonize()}
    if(length(cc2$level)==2){cc2 <- st_as_sf(ContourLines2SLDF(list(cc[[length(cc)-1]],
                                                              cc[[length(cc)-sample(1:3,1)]],
                                                              cc[[length(cc)-sample(4:6,1)]]) ) ) %>% st_polygonize()}
    if(length(cc2$level)==2){cc2 <- st_as_sf(ContourLines2SLDF(list(cc[[length(cc)-1]],
                                                              cc[[length(cc)-sample(1:3,1)]],
                                                              cc[[length(cc)-sample(4:6,1)]]) ) ) %>% st_polygonize()}
    if(length(cc2$level)==2){cc2 <- st_as_sf(ContourLines2SLDF(list(cc[[length(cc)-1]],
                                                                    cc[[length(cc)-sample(1:3,1)]],
                                                                    cc[[length(cc)-sample(4:6,1)]]) ) ) %>% st_polygonize()}
    if(length(cc2$level)==2){cc2 <- st_as_sf(ContourLines2SLDF(list(cc[[length(cc)-1]],
                                                                    cc[[length(cc)-sample(1:3,1)]],
                                                                    cc[[length(cc)-sample(4:6,1)]]) ) ) %>% st_polygonize()}
    if(length(cc2$level)==2){cc2 <- st_as_sf(ContourLines2SLDF(list(cc[[length(cc)-1]],
                                                                    cc[[length(cc)-sample(1:3,1)]],
                                                                    cc[[length(cc)-sample(4:6,1)]]) ) ) %>% st_polygonize()}
    
    cc2 <- as(cc2,"Spatial")
    cc2$id <- c("1","2","3")
  
    #ccL[[i]] <- cc2
  
    #plot(cc2, add=T, col=col1)
    map <- map+tm_shape(cc2)+
               tm_polygons("id", style="fixed", palette=col1,breaks=c("1","2","3"), border.col = "grey", alpha=.5, legend.show=F)
  
    #map+tm_shape(cc2)+
    #tm_polygons("id", style="fixed", palette=col1,breaks=c("1","2","3"), border.col = "grey", alpha=.5)
  
  }

  tmap_save(tm = map, filename = paste0("../results/rase/melipona/",taxon,"_",corte-1,"_Ma",".pdf"),dpi = 300)
   
  try(system(paste0("pdftocairo -png"," ", "-png"," ","../results/rase/melipona/",taxon,"_",corte-1,"_Ma",".pdf"," ",
                    "../results/rase/melipona/",taxon,"_",corte-1,"_Ma",".png")))
  try(system("rm ../results/rase/melipona/*.pdf"))

}
