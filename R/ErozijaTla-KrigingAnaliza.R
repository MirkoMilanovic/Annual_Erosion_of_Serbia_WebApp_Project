# podešavanje radnog direktorijuma
setwd("C:/Users/Mirko/Desktop/DRZ. KART. PROJ/R")

# podešavanje polaznih parametara
options(prompt="> ", continue="+ ", digits=8, width=70, show.signif.stars=T)
# ucitavanje paketa
library(sp)
library(rgdal)
library(gstat)
library(rgeos)
library(foreign)
library(ggplot2)
library(landsat)

# ucitavanje tabele datih podataka za padavine
podaci <- read.table(file="podaci.txt", header=TRUE, sep="\t", dec=".", na.strings=c("NA", "-", "?"))

# histogram zavisnosti broja stanica od registrovane kolicine padavina u toku 2016. godine
hist(podaci$PADAVINE, xlab="Suma padavine u 2016. god. (mm)", col="yellow", main="Zavisnost broja klimatoloških stanica od registrovane kolicine padavina na njima u 2016. godini")

# kopiranje datih podataka i transformacija koordinata
podaci1 = podaci
coordinates(podaci1) <- c("DUZINA", "SIRINA")
proj4string(podaci1) <- CRS("+init=epsg:4326")
podaci1 <- spTransform(podaci1, CRS("+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 
                                                           +y_0=0 +ellps=bessel 
                                                           +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232   						 +units=m +no_defs"))


# bubble plot registrovane kolicine padavina na klimatološkim stanicama
bubble(subset(podaci1,!is.na(PADAVINE)), scales=list(draw=T), "PADAVINE",  col="red", pch=15, maxsize=2)

# ucitavanje rastera DEM Republike Srbije i za šumski pokrivac
dem <- readGDAL("dem250.tif")
#dem@proj4string
sume <- readGDAL("sume250.tif")

plot(dem)
plot(sume)

# racunanje nagiba na osnovu podataka DEM Srbije, export karte nagiba u tiff formatu
dem.slopeasp <- slopeasp(dem)
nagib <- dem.slopeasp$slope
plot(nagib)

writeGDAL(nagib, "C:/Users/Mirko/Desktop/DRZ. KART. PROJ/R/nagib250.tif", drivername="GTiff")
 
# prostorno povezivanje podataka, dodavanje podataka u baznu sliku
nagib@data[is.na(nagib@data)] <- 0
sume@data[is.na(sume@data)] <- 0
 
s = over(nagib, sume)

names(s)[names(s)=="band1"] <- "s"

nagib@data <- cbind(s, nagib@data)




# prostorno povezivanje podataka (tacke i grid)
podaci11 = over(podaci1, nagib)
podaci1@data <- cbind(podaci11, podaci1@data)

# pozicije klimatoloških stanica na podlozi DEM Pepublike Srbije prikazane paletom boja za kopno
plot(dem, main="Pozicije klimatoloških stanica", col=terrain.colors(20))
points(subset(podaci1, !is.na(PADAVINE)), col ="red" , pch = 17)
  

# aproksimaciona kriva zavisnosti nadmorske visine stanica od registrovanih padavina na njima
scatter.smooth(podaci1$VISINA, podaci1$PADAVINE, col="blue", xlab="Nadmorska visina (m)", ylab="Padavine (mm)")


# brisanje NA vrednosti za PADAVINE nepotrebnih podataka
drops <- c("VISINA")
podaci1 <- podaci1[,!(names(podaci1) %in% drops)]
podaci1 <- podaci1[!(is.na(podaci1$PADAVINE)),]


# variogram za blok kriging
v = variogram(PADAVINE~band1-s, podaci1)
m <- vgm(nugget=0, model="Exp", range=30000, psill=30000)
fitv = fit.variogram(v, m)
plot(v,m)


# universal blok kriging
EROZIJA.block <- krige(PADAVINE ~band1-s, podaci1, nagib, block = c(40,40), model=fitv)
spplot(EROZIJA.block["var1.pred"], main = "Blok kriging - predikcija")
spplot(EROZIJA.block["var1.var"],  main = "Blok kriging - varijansa")

# standardna greška blok kriging metode 
mean(sqrt(EROZIJA.block$var1.var), na.rm = TRUE)

# extraktovanje rastera predikcije
writeGDAL(EROZIJA.block["var1.pred"], "C:/Users/Mirko/Desktop/DRZ. KART. PROJ/R/erozija250.tif", drivername="GTiff")


# kopiranje ulaznih podataka, transformacija koordinata i prostorno preklapanje podataka
podaci2 = podaci
coordinates(podaci2) <- c("DUZINA", "SIRINA")
proj4string(podaci2) <- CRS("+init=epsg:4326")
podaci2 <- spTransform(podaci2, CRS("+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 
                                                           +y_0=0 +ellps=bessel 
                                                           +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232   					 +units=m +no_defs"))

block <- over(podaci2, EROZIJA.block)

podaci2 <- podaci2[,!(names(podaci2) %in% drops)]
View(podaci2)

# dodavanje vrednosti u tabelu sa stanicama
podaci2$block <- block$var1.pred

podaci2 <- podaci2[,!(names(podaci2) %in% drops)]
View(podaci2)

# bubble plot za predvidenu eroziju tla
bubble(podaci2, "block", scales=list(draw=T), col="green", pch=18, maxsize=2)


# formiranje i ekstrakcija tabele za eroziju za lokacije stanica i lokacije od interesa
podaci2$H <- podaci$VISINA
podaci2$EROZIJA <- podaci2$block
podaci2 <- as.data.frame(podaci2)
podaci2 <- podaci2[c("STANICA", "DUZINA", "SIRINA", "H", "EROZIJA")]
names(podaci2)[names(podaci2)=="DUZINA"] <- "Y"
names(podaci2)[names(podaci2)=="SIRINA"] <- "X"
View(podaci2)
write.table(podaci2, "podaciNovo.txt", sep="\t", row.names = FALSE)

