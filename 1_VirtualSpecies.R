library(dplyr)
library(ggplot2)
library(virtualspecies)
library(raster)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(stars)

source("SpatialProba.R")
source("rescale.R")
outdir <- "outputs20230824/figures/"
outsim <- "outputs20230824"

# Worldclim <- raster::getData('worldclim', var='bio', res = 10) 
BioData <- geodata::worldclim_global(var = c("bio"),
                                     path = "/your/home/dir",
                                     res = 10)
envData <- terra::crop(BioData, terra::ext(-12, 25, 36, 60))
names(envData) <- paste0("bio", seq_along(names(envData)))
envData <- envData[[c("bio1", "bio12") ]]
rm(BioData)

## 1. Virtual species ----
### 1.1 Virtual tree species ----
myCoeff <- c(intercept= 1, bio1 = -1, bio12= 0.01)
treeSuitab <- round(SpatialProba(coefs = myCoeff, env.rast = envData, quadr_term = NULL,  marginalPlots=FALSE),3)
plot(treeSuitab, main = "Virtual tree \n climatic suitability")

# create virtualspecies Pkg like object
my.first.species <- list()
my.first.species$suitab.raster <- raster::raster(treeSuitab)

# compute marginal effects
env_mat.tmp <- cbind.data.frame(as.data.frame(envData[[c("bio1", "bio12")]]))
myList_tree <- list(bio1 = data.frame(name="bio1 (°C)", 
                                    env=env_mat.tmp$bio1,
                                    PA=plogis(as.matrix(cbind.data.frame(intercept=1, 
                                                                         bio1=env_mat.tmp$bio1, 
                                                                         bio12=mean(env_mat.tmp$bio12, na.rm=TRUE)))%*%unname(myCoeff))),  
                  bio12 = data.frame(name="bio12 (mm)", 
                                     env=env_mat.tmp$bio12,  
                                     PA= plogis(as.matrix(cbind.data.frame(intercept=1, 
                                                                           bio1=mean(env_mat.tmp$bio1, na.rm=TRUE), 
                                                                           bio12=env_mat.tmp$bio12))%*%unname(myCoeff)))
)
myList_tree <- do.call(rbind, myList_tree)

### 1.2 Virtual herb species ----
#### 1.2.1 Germination rate ----
# Development of germination rate function depending on the virtual tree suitability
fake_tree_suitab <- seq(0,1, by = 0.01)  
a <- 1          
g.r <- 3 
shade.germ <- a*log(g.r*fake_tree_suitab+1)
shade.germ <- rescale(shade.germ, 
                      x.min = min(shade.germ), 
                      x.max = max(shade.germ), 
                      new.min = 0, new.max = 1)
plot(fake_tree_suitab, 
     shade.germ, 
     type="l", 
     xlab = "Virtual tree climatic suitability", 
     ylab = "Germination rate")

# Compute a germination probability raster
shade.germ.r <- a*exp(g.r*my.first.species$suitab.raster) 
shade.germ.r <- round(rescale(shade.germ.r, 
                              x.min = 1, 
                              x.max = max(terra::values(shade.germ.r), na.rm = TRUE), 
                              new.min = 0, 
                              new.max = 1), 3
)

par(mfrow=c(1,2))
plot(my.first.species$suitab.raster, main = "Virtual tree \n climatic suitability")
plot(shade.germ.r, main = "Germination rate of the virtual herb species")
par(mfrow=c(1,1))

#### 1.2.2 Virtual "herbaceous plant" species ----
myCoeff2 <- c(intercept= 1, bio1 = -0.85, bio12= 0.015)
herbSuitab <- SpatialProba(coefs = myCoeff2, env.rast =envData, quadr_term = NULL,  marginalPlots=FALSE)
plot(herbSuitab)

# create virtualspecies Pkg like object
my.second.species <- list()
my.second.species$suitab.raster <- herbSuitab
plot(my.second.species$suitab.raster, main = "Virtual herbaceous \n climatic suitability")

# compute marginal effects
env_mat.tmp <- cbind.data.frame(as.data.frame(envData[[c("bio1", "bio12")]]))
myList_herb <- list(bio1 = data.frame(name="bio1 (°C)", 
                                    env=env_mat.tmp$bio1,
                                    PA=plogis(as.matrix(cbind.data.frame(intercept=1, 
                                                                         bio1=env_mat.tmp$bio1, 
                                                                         bio12=mean(env_mat.tmp$bio12, na.rm=TRUE)))%*%unname(myCoeff2))),  
                  bio12 = data.frame(name="bio12 (mm)", 
                                     env=env_mat.tmp$bio12,  
                                     PA= plogis(as.matrix(cbind.data.frame(intercept=1, 
                                                                           bio1=mean(env_mat.tmp$bio1, na.rm=TRUE), 
                                                                           bio12=env_mat.tmp$bio12))%*%unname(myCoeff2)))
)

myList_herb <- do.call(rbind, myList_herb)

# intersect the climatic aspect of the herb species with the germination rate
tmp <- round(shade.germ.r* raster::raster(my.second.species$suitab.raster),3)

par(mfrow=c(2,2))
plot(my.first.species$suitab.raster, main = "Virtual tree climatic suitability")
plot(my.second.species$suitab.raster, main = "Virtual herb climatic suitability")
plot(shade.germ.r, main = "Germination probability of the \n herbaceous virtual species")
plot(tmp, main = "Virtual herb habitat suitability")
par(mfrow=c(1,1))

# assing the virtual herb species habitat suitability restricted by the germination rate to the virtualspecies like object
my.second.species$suitab.raster <- tmp

## 2. Response curves for both virtual species ----
myList_herb$species <- "Herb"
myList_tree$species <- "Tree"
cols<-c("Tree"="darkorange", "Herb"="darkgreen")

bio1_resp <- myList_herb %>% 
    bind_rows(myList_tree) %>% 
    as_tibble() %>% 
    mutate(species=factor(species, levels=c("Tree", "Herb"))) %>% 
    filter(name=="bio1 (°C)")  %>% 
    ggplot2::ggplot(ggplot2::aes(x = env, y = PA, col=species)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = cols)+
    ggplot2::xlim(-5,20)+
    ggplot2::labs(x="BIO1 (°C)", y="Suitability", col="Virtual species") +
    ggplot2::theme_classic()+
    theme(aspect.ratio = 1, 
          legend.background=element_blank(),
          panel.grid = element_blank(),
          legend.position = 'bottom',
          text = element_text(size=16), 
          strip.text = element_text(size=16),
          legend.text = element_text(size=16,angle = 0), 
          legend.title = element_text(size=16),
          legend.key.size = unit(1.5, 'cm'))

bio12_resp <- myList_herb %>% 
    bind_rows(myList_tree) %>% 
    as_tibble() %>% 
    mutate(species=factor(species, levels=c("Tree", "Herb"))) %>% 
    filter(name=="bio12 (mm)")  %>% 
    ggplot2::ggplot(ggplot2::aes(x = env, y = PA, col=species)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = cols)+
    ggplot2::xlim(0,1500)+
    ggplot2::labs(x="BIO12 (mm)", y="Suitability", col="Virtual species") +
    ggplot2::theme_classic()+
    theme(aspect.ratio = 1, 
          legend.background=element_blank(),
          panel.grid = element_blank(),
          legend.position = 'bottom',
          text = element_text(size=16), 
          strip.text = element_text(size=16),
          legend.text = element_text(size=16,angle = 0), 
          legend.title = element_text(size=16),
          legend.key.size = unit(1.5, 'cm'))

germRate_resp <- cbind.data.frame(treeSuitab=fake_tree_suitab,
                                  germ=round(shade.germ,3),
                                  species="Herb") %>% 
  ggplot2::ggplot(ggplot2::aes(x = treeSuitab, y = germ, col=species)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(values = cols)+
  ggplot2::labs(x="Virtual tree suitability", y="Germination rate", col="Virtual species") +
  ggplot2::theme_classic()+
  theme(aspect.ratio = 1, 
        legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), 
        legend.title = element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

pp<-ggpubr::ggarrange(bio1_resp, bio12_resp, germRate_resp,
                      ncol=3, nrow=1, 
                      labels=c('A', 'B','C'), 
                      common.legend = TRUE, legend="bottom")

pp
outname <- paste(outdir, "respCurves_", Sys.Date(),".pdf", sep="")
ggsave(pp, filename = outname, width = 16, height = 8, device='pdf', dpi=320)

# get species range
myList_herb %>% 
  bind_rows(myList_tree) %>% 
  as_tibble() %>%   
  filter(PA >0.025 & PA <0.975) %>% 
  dplyr::select(name, env, species) %>% 
  group_by(name, species) %>% 
  summarise(range=round(range(env))) %>% 
  tidyr::drop_na() 

## 3. Habitat suitability spatial outputs ----
World <- ne_countries(scale = "medium", returnclass = "sf")
eu <- subset(World, continent =="Europe")

p1 <- ggplot(data = World) +
        geom_sf()+
        geom_stars(data = stars::st_as_stars(envData$bio1) )+
        labs(x="Longitude",y="Latitude", fill="BIO1: Mean annual temperature (°C)")+
        scale_fill_viridis(option="inferno",  limits=c(-5,20),
                           oob = scales::squish, na.value = "transparent")+
        theme_void()+
        theme(plot.title = element_text(size=14,face = 'bold'),
              legend.background=element_blank(),
              strip.text.x = element_text(size=12,face = 'bold'),
              legend.position = 'bottom',
              legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
              legend.key.size = unit(2, 'cm'))+
        guides(fill = guide_colourbar(title.position="top", 
                                      title.hjust = 0.5, 
                                      barwidth = 20, barheight = 0.8),
               size = guide_legend(title.position="top", title.hjust = 0.5))+
        coord_sf(xlim = c(-12, 25), ylim = c(36, 60), expand = FALSE)

p2 <- ggplot(data = World) +
        geom_sf()+
        geom_stars(data = stars::st_as_stars(envData$bio12) )+
        labs(x="Longitude",y="Latitude", fill="BIO12: Mean annual precipitation (mm)")+
        scale_fill_viridis(option="cividis",  direction= -1, limits=c(0,2500),
                           oob = scales::squish, na.value = "transparent")+
        theme_void()+
        theme(plot.title = element_text(size=14,face = 'bold'),
              legend.background=element_blank(),
              strip.text.x = element_text(size=12,face = 'bold'),
              legend.position = 'bottom',
              legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
              legend.key.size = unit(2, 'cm'))+
        guides(fill = guide_colourbar(title.position="top", 
                                      title.hjust = 0.5, 
                                      barwidth = 20, barheight = 0.8),
               size = guide_legend(title.position="top", title.hjust = 0.5))+
        coord_sf(xlim = c(-12, 25), ylim = c(36, 60), expand = FALSE)

p3 <- ggplot(data = World) +
        geom_sf()+
        geom_stars(data = stars::st_as_stars(shade.germ.r) )+
        labs(x="Longitude",y="Latitude", fill="Virtual herb germination rate", color= "Observations")+
        scale_fill_viridis(option="rocket",  limits = c(0, 1),  oob = scales::squish, na.value = "transparent")+
        theme_void()+
        theme(plot.title = element_text(size=14,face = 'bold'),
              legend.background=element_blank(),
              strip.text.x = element_text(size=12,face = 'bold'),
              legend.position = 'bottom',
              legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
              legend.key.size = unit(2, 'cm'))+
        guides(fill = guide_colourbar(title.position="top", 
                                      title.hjust = 0.5, 
                                      barwidth = 20, barheight = 0.8),
               size = guide_legend(title.position="top", title.hjust = 0.5))+
        coord_sf(xlim = c(-12, 25), ylim = c(36, 60), expand = FALSE)

p4 <- ggplot(data = World) +
        geom_sf()+
        geom_stars(data = stars::st_as_stars(my.first.species$suitab.raster) )+
        labs(x="Longitude",y="Latitude", fill="Virtual tree suitability", color= "Observations")+
        scale_fill_viridis(option="viridis",   oob = scales::squish, na.value = "transparent")+
        theme_void()+
        theme(plot.title = element_text(size=14,face = 'bold'),
              legend.background=element_blank(),
              strip.text.x = element_text(size=12,face = 'bold'),
              legend.position = 'bottom',
              legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
              legend.key.size = unit(2, 'cm'))+
        guides(fill = guide_colourbar(title.position="top", 
                                      title.hjust = 0.5, 
                                      barwidth = 20, barheight = 0.8),
               size = guide_legend(title.position="top", title.hjust = 0.5))+
        coord_sf(xlim = c(-12, 25), ylim = c(36, 60), expand = FALSE)

p5 <- ggplot(data = World) +
        geom_sf()+
        geom_stars(data = stars::st_as_stars(my.second.species$suitab.raster) )+
        labs(x="Longitude",y="Latitude", fill="Virtual herb suitability", color= "Observations")+
        scale_fill_viridis(option="viridis",  limits = c(0, 1),  oob = scales::squish, na.value = "transparent")+
        theme_void()+
        theme(plot.title = element_text(size=14,face = 'bold'),
              legend.background=element_blank(),
              legend.position = 'bottom',
              legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
              legend.key.size = unit(2, 'cm'))+
        guides(fill = guide_colourbar(title.position="top", 
                                      title.hjust = 0.5, 
                                      barwidth = 20, barheight = 0.8),
               size = guide_legend(title.position="top", title.hjust = 0.5))+
        coord_sf(xlim = c(-12, 25), ylim = c(36, 60), expand = FALSE)

loc <- sf::st_as_sf(raster::sampleRandom(raster(envData$bio1), size =  500, sp = TRUE))
p6 <- ggplot(data = World) +
        geom_sf(col=NA)+
        geom_sf(data=loc, aes(col="Sampling locations"), alpha=0.5, size=0.8)+
        labs(x="Longitude",y="Latitude", colour="")+
        scale_color_manual(values="black")+
        theme_void()+
        theme(legend.position = "bottom",  text = element_text(size=14), 
              legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),)+
        guides(fill = guide_colourbar(title.position="top", 
                                      title.hjust = 0.5, barwidth = 20, barheight = 0.8),
               size = guide_legend(title.position="top", title.hjust = 0.5), 
               color = guide_legend(override.aes = list(size = 3)))+
        coord_sf(xlim = c(-12, 25), ylim = c(36, 60), expand = FALSE)

#assemble plots
library(cowplot)
pp <- plot_grid(p1, p2, p3, p4, p5, p6, nrow=2,  labels = "auto")
pp
outname <- paste0(outdir, "bios_VS_", Sys.Date(), ".pdf")
ggsave(pp, filename = outname, width = 16, height = 8, device='pdf', dpi=320)


## 4.  Convert virtual species suitability indices to presence-absence ----
new.pres.tree <- virtualspecies::convertToPA(x = my.first.species$suitab.raster,
                                             beta="random",
                                             alpha = -0.05, plot = FALSE,
                                             species.prevalence = 0.4)

# Sensitivity analysis on the herb species geographic prevalence
geog.prev <- c(0.25, 0.5, 0.75)
new.pres.herb <- lapply(geog.prev, function(y){convertToPA(x = my.second.species$suitab.raster, beta="random", alpha = -0.05, plot = FALSE, species.prevalence = y)})
names(new.pres.herb) <- paste0("geog.prev_", geog.prev)

par(mfrow=c(2,2))
plot(new.pres.tree$pa.raster, main = "Virtual tree species \n (geog. prev = 0.5)", legend = FALSE)
plot(new.pres.herb$geog.prev_0.25$pa.raster, main = "Virtual herbaceous species  \n (geog. prev = 0.25)")
plot(new.pres.herb$geog.prev_0.5$pa.raster, main = "Virtual herbaceous species  \n (geog. prev = 0.5)")
plot(new.pres.herb$geog.prev_0.75$pa.raster, main = "Virtual herbaceous species  \n (geog. prev = 0.75)")
par(mfrow=c(1,1))

### 4.1 Presence-absence spatial outputs ----
mySta <- raster::stack(new.pres.tree$pa.raster, 
                     new.pres.herb$geog.prev_0.25$pa.raster, 
                     new.pres.herb$geog.prev_0.5$pa.raster, 
                     new.pres.herb$geog.prev_0.75$pa.raster)

names(mySta)<- c("A", "B", "C", "D")
mySta <- stars::st_as_stars(mySta)

p <-ggplot(data = World) +
  geom_sf()+
  geom_stars(data = mySta)+
  labs(fill="")+
  scale_fill_manual(values =viridis::viridis(3),
                    na.value = "transparent", labels = c("Absence", "Presence"), 
                    na.translate = F)+
  theme_void()+
  facet_wrap(~band, ncol = 2, 
             labeller = labeller(band = c("A" = "Tree \n (geog. prev. = 0.40)",
                                          "B" = "Herb \n (geog. prev. = 0.25)", 
                                          "C" = "Herb \n (geog. prev. = 0.50)", 
                                          "D" = "Herb \n (geog. prev. = 0.75)")))+ 
  theme(plot.title = element_text(size=14,face = 'bold'),
        legend.background=element_blank(),
        legend.position = 'bottom',
        strip.text = element_text(size = 14),
        legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(1, 'cm'))+
  coord_sf(xlim = c(-12, 25), ylim = c(36, 60), expand = FALSE)

p

outname <- paste0(outdir, "VS_geogPrev", Sys.Date(), ".pdf")
ggsave(p, filename = outname, width = 10, height = 8, device='pdf', dpi=320)


## 5. Training dataset preparation ----
{set.seed(123)
  rastTemplate <- raster::raster(envData$bio1)
  rastTemplate[!is.na(rastTemplate)] <- 1
  loc <- raster::sampleRandom(rastTemplate, size =  500, sp = TRUE)}
plot(rastTemplate)
plot(loc, add=TRUE)
names(loc) <- "train"
crs(loc) <- crs(new.pres.tree$pa.raster)

# Stack variabels
pred.stack<-lapply(1:length(new.pres.herb), function(x){raster::stack(raster::raster(envData$bio1), 
                                                                      raster::raster(envData$bio12),
                                                                      shade.germ.r, 
                                                                      new.pres.tree$pa.raster,
                                                                      new.pres.herb[[x]]$pa.raster,
                                                                      new.pres.herb$geog.prev_0.25$suitab.raster)})

names(pred.stack[[1]])<-c("bio1", "bio12", "herb.germ", "tree.pa",  "herb.pa", "obs.herb.suit")
names(pred.stack[[2]])<-c("bio1", "bio12", "herb.germ", "tree.pa",  "herb.pa", "obs.herb.suit")
names(pred.stack[[3]])<-c("bio1", "bio12", "herb.germ", "tree.pa",  "herb.pa", "obs.herb.suit")


outname <- paste0(outsim, "/stacked_predictors_", Sys.Date(), ".RDS")
saveRDS(pred.stack, outname)

myGrid<-raster::raster(ext=extent(rastTemplate), res=c(8,9))
{set.seed(123)
  values(myGrid)<-sample( 1:ncell(myGrid), size=length(1:ncell(myGrid)))}
plot(myGrid)
plot(rastTemplate)
plot(raster::rasterToPolygons(myGrid), add=T)
loc$train<-raster::extract(myGrid, loc)

outname <- paste0(outsim, "/crossValGrid_", Sys.Date(), ".RDS")
saveRDS(myGrid, outname)

outname <- paste0(outsim, "/samplinLoc_", Sys.Date(), ".RDS")
saveRDS(loc, outname)