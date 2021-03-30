#Codigo control Omar

#########################################
#              Libraries                #
#########################################

set.seed(1331)

library(rase)
library(ape)
library(rgeos)
library(maptools)
library(raster)
library(coda)
library(spatstat.geom)


#########################################
#              Starting Vars            #
#########################################

setwd("~/Documentos/Omar/Projects/Amazonas_2020/2020/raseDatosRAmazonia/R/")

sp <- "melipona"
dist <- "log"
rep <- 50000

taxon <- paste0(sp,"-",dist,".tre",sep="")
taxon2 <- paste0(sp,"/",sp,"AllSp",sep="")

rutaTrees <- "../data/trees/"
rutaDist <- "../data/shp/"

arbol <- read.nexus( paste0(rutaTrees,taxon,sep="") )
datosS <- paste0(rutaDist,taxon2,".shp",sep="")

a <- rgdal::readOGR(datosS)

#########################################
#              Pretty Data              #
#########################################

newArbol <- drop.tip( arbol, which(!(arbol$tip.label%in%gsub(".shp","",a$Sp))), subtree=FALSE )
newArbol <- multi2di(newArbol, random=TRUE) #Chao politomias al azar

blMin <- min(newArbol$edge.length[newArbol$edge.length !=0]) /100
newArbol$edge.length[newArbol$edge.length == 0] <- blMin

nodo <- length(newArbol$tip.label)+1
root(newArbol,node=nodo) #Enraizar 
is.rooted(newArbol)

plot(newArbol)


write.tree(newArbol,file=paste0(rutaTrees,sp,"-",dist,"-prunned",rep,".tre",sep=""))

#solo para stenodermatinae
#a <- a[-13,]

#########################################
#                 RASE                  #
#########################################


distributionData <- shape.to.rase(a)

poly <- name.poly(distributionData, newArbol,  poly.names = newArbol$tip.label)

rase_results <-    rase(newArbol, poly , niter=rep)

rasemcmc <- coda::mcmc(rase_results)

#plot(rasemcmc)

#########################################
#              Saving RDS               #
#########################################

archivo <- paste0("../results/rase/",sp,"/",sp,"-",dist,rep,"rase.RDS",sep="")
saveRDS(rasemcmc,file=archivo)
