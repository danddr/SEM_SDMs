library(dplyr)
library(piecewiseSEM)
library(semEff)
library(raster)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(stars)
library(ggplot2)

outdir <- "/outputs20230824/figures/"
outsim <- "/outputs20230824/"
geog.prev <- c(0.25, 0.5, 0.75)

#---- 1. load and process spatial predictors ----
pred.stack<-readRDS("/outputs20230824/stacked_predictors_2023-08-24.RDS")
# plot(raster::crop(pred.stack[[2]]$obs.herb.suit, raster::extent(c(6, 12, 48, 52))))

# crop predictors to subset of the AOI
spatExt <- raster::extent(c(6, 12, 48, 52))
pred.stack.sub<-lapply(pred.stack, function(x){raster::as.data.frame(raster::crop(x, spatExt), xy=TRUE)})
pred.stack.sub

# store observed herb suitability
obs.herb.suit <-raster::crop(pred.stack[[1]]$obs.herb.suit, spatExt)

# scale predictors
pred.stack.sub <- lapply(pred.stack.sub, function(x){data.frame(x[,1:2], scale(x[,3:5], center = TRUE, scale = TRUE), x[, 6:8])})
names(pred.stack.sub)<-geog.prev

#---- 2. load SEM and GLM model coefficients ----
SEM.list <- list.files(outsim, pattern="semEffects", full.names = TRUE)
SEM.graph <- list.files(outsim, pattern="vs.semList", full.names = TRUE)
GLM.list <- list.files(outsim, pattern="glmEffects", full.names = TRUE)

GLM.list <- list("0.25" = stringr::str_subset(GLM.list, "0.25_"),
                 "0.5" = stringr::str_subset(GLM.list, "0.5_"), 
                 "0.75"= stringr::str_subset(GLM.list, "0.75_")
                  )

# storing objects
SEMprev0.25 <- list()
SEMprev0.5 <- list()
SEMprev0.75 <- list()

GLMprev0.25 <- list()
GLMprev0.5 <- list()
GLMprev0.75 <- list()

#---- 3. predictions ----
for(y in 1:length(SEM.list)){
  # y=1
  sem.eff <- readRDS(SEM.list[y])
  vs.semList <- readRDS(SEM.graph[y])
    for(i in 1:length(geog.prev)){
      # i=1
      message(paste("SEM and GLM spatial predictions: Processing fold", y, "with geographic prevalence",  geog.prev[i]))
      # Extract the (total) effects of each variable to use for prediction
      tot <- getTotEff(sem.eff[[i]], "herb.pa")
      tot.b <- getTotEff(sem.eff[[i]], "herb.pa", type = "boot")
      
      # SEM prediction
      mod <- vs.semList[[i]]$herb.pa
      tmp.df <- mod$data
      r.SEMpred <- pred.stack.sub[[i]][, c("x", "y")]
      newPredDat <- pred.stack.sub[[i]][,c("herb.germ", "tree.pa",  "bio1", "bio12"  )]
      r.SEMpred$herb.pred <- semEff::predEff(mod=mod, 
                                  newdata =  newPredDat, #tmp.df 
                                  effects = tot,
                                  # eff.boot = tot.b,
                                  type = "response")
      r.SEMpred <- raster::rasterFromXYZ(r.SEMpred)
      names(r.SEMpred)<-"SEM"
 
      # GLM pred
      myglm <- readRDS(GLM.list[[i]][y])
      r.GLMpred <- pred.stack.sub[[i]][, c("x", "y")]
      r.GLMpred$herb.pred <- predict(myglm, newdata=pred.stack.sub[[i]], type="response")
      r.GLMpred <- raster::rasterFromXYZ(r.GLMpred)
      names(r.GLMpred) <- "GLM"
      
      if(i == 1){
        SEMprev0.25[[y]] <- r.SEMpred
        GLMprev0.25[[y]] <- r.GLMpred
      }
      
      if(i == 2){
        SEMprev0.5[[y]] <- r.SEMpred
        GLMprev0.5[[y]] <- r.GLMpred
      }
      
      if(i ==3){
        SEMprev0.75[[y]] <- r.SEMpred
        GLMprev0.75[[y]] <- r.GLMpred
      }

    }
  
}

SEMprev0.25 <- do.call(raster::stack, SEMprev0.25)
SEMprev0.25 <- raster::calc(SEMprev0.25, fun=median)

SEMprev0.5 <- do.call(raster::stack, SEMprev0.5)
SEMprev0.5 <- raster::calc(SEMprev0.5, fun=median)

SEMprev0.75 <- do.call(raster::stack, SEMprev0.75)
SEMprev0.75 <- raster::calc(SEMprev0.75, fun=median)

SEMstack<-raster::stack(SEMprev0.25, SEMprev0.5, SEMprev0.75)
names(SEMstack) <- paste0("SEM_", geog.prev)

GLMprev0.25 <- do.call(raster::stack, GLMprev0.25)
GLMprev0.25 <- raster::calc(GLMprev0.25, fun=median)

GLMprev0.5 <- do.call(raster::stack, GLMprev0.5)
GLMprev0.5 <- raster::calc(GLMprev0.5, fun=median)

GLMprev0.75 <- do.call(raster::stack, GLMprev0.75)
GLMprev0.75 <- raster::calc(GLMprev0.75, fun=median)

GLMstack <- raster::stack(GLMprev0.25, GLMprev0.5, GLMprev0.75)
names(GLMstack) <- paste0("GLM_", geog.prev)

myStack <- raster::stack(SEMstack, 
                     GLMstack)


#---- 4. Outputs plots ----
World <- ne_countries(scale = "medium", returnclass = "sf")
eu <- subset(World, continent =="Europe")
bbox <- data.frame(lon=c(6, 12), 
           lat=c(48,52)) %>% 
  st_as_sf(coords = c("lon", "lat"), 
           crs = 4326) %>% 
  st_bbox() %>% 
  st_as_sfc()

#observed suitab
p1 <- ggplot(data = World) +
        geom_sf()+
        geom_stars(data = stars::st_as_stars(pred.stack[[1]]$obs.herb.suit) )+
        geom_sf(data=bbox, fill=NA, col="coral", linewidth=1.3)+
        labs(x="Longitude",y="Latitude", fill="Observed herb suitability")+
        scale_fill_viridis(option="viridis",  limits=c(0,1),
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


p1
outname<-paste0(outdir, "HerbSuit_obs", Sys.Date(), ".pdf")
ggsave(p1, filename = outname, width = 16, height = 8, device='pdf', dpi=320)


#predicted suitab
names(myStack) <- c("A", "B", "C", "D", "E", "F")
mySta <- stars::st_as_stars(myStack)
names(mySta)


p2 <- ggplot() +
    geom_stars(data = mySta)+
  labs(fill="Predicted herb suitability")+
  scale_fill_viridis(option="viridis",  limits = c(0, 1),  oob = scales::squish, na.value = "transparent")+
  theme_void()+
  facet_wrap(~band, ncol = 3, 
             labeller = labeller(band = c("A" = "SEM \n (geog. prev. = 0.25)",
                                          "B" = "SEM \n (geog. prev. = 0.50)", 
                                          "C" = "SEM \n (geog. prev. = 0.75)", 
                                          "D" = "GLM \n (geog. prev. = 0.25)",
                                          "E" = "GLM \n (geog. prev. = 0.50)", 
                                          "F" = "GLM \n (geog. prev. = 0.75)")))+ 
  theme(plot.title = element_text(size=14,face = 'bold'),
        legend.background=element_blank(),
        legend.position = 'bottom',
        strip.text = element_text(size = 14),
        legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(1, 'cm'))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  coord_cartesian()

p2
outname <- paste0(outdir, "HerbSuit_pred", Sys.Date(), ".pdf")
ggsave(p2, filename = outname, width = 8, height = 8, device='pdf', dpi=320)


#---- 5. Accuracy outputs ----
library(Metrics)
# compute cor and rmse with the observed distribution
sp.preds <- na.omit(raster::as.data.frame(stack(myStack, obs.herb.suit)))
rmse.distribution <- apply(sp.preds[, 1:6], 2, function(x){round(Metrics::rmse(sp.preds$obs.herb.suit, x),3)})
rmse.distribution <- data.frame(ModelType=attributes(rmse.distribution)$names, 
                                rmse.distribution=rmse.distribution)


names(rmse.distribution) <- c("Model", "RMSE")
rmse.distribution$Model <- rep(c("SEM", "GLM"), each= 3)
rmse.distribution$geog.prev <- rep(geog.prev,2 )

knitr::kable(rmse.distribution, align="c")

outname<-paste0(outdir, "gof_Geo_metrics_", Sys.Date(), ".csv")
write.csv(rmse.distribution, outname)

#spatial rmse
rmseOut<-list()
for(i in 1:nlayers(myStack)){
  rmseOut[[i]]<-sqrt(mean((obs.herb.suit - myStack[[i]])^2, na.rm = TRUE))
  }
rmseOut<-do.call(stack, rmseOut)
names(rmseOut) <- c("A", "B", "C", "D", "E", "F")

round(cellStats(rmseOut, median),3)

rmse_spatial <- ggplot() +
  geom_stars(data = stars::st_as_stars(rmseOut))+
  labs(fill="RMSE")+
  scale_fill_viridis(option="plasma",  
                     limits = c(0, 1),
                     direction = 1, 
                     oob = scales::squish, na.value = "transparent")+
  theme_void()+
  facet_wrap(~band, ncol = 3, 
             labeller = labeller(band = c("A" = "SEM \n (geog. prev. = 0.25)",
                                          "B" = "SEM \n (geog. prev. = 0.50)", 
                                          "C" = "SEM \n (geog. prev. = 0.75)", 
                                          "D" = "GLM \n (geog. prev. = 0.25)",
                                          "E" = "GLM \n (geog. prev. = 0.50)", 
                                          "F" = "GLM \n (geog. prev. = 0.75)")))+ 
  theme(plot.title = element_text(size=14,face = 'bold'),
        legend.background=element_blank(),
        legend.position = 'bottom',
        strip.text = element_text(size = 14),
        legend.text = element_text(size=12,angle = 0), legend.title = element_text(size=14),
        legend.key.size = unit(1, 'cm'))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  coord_cartesian()

rmse_spatial
outname <- paste0(outdir, "rmse_spatial", Sys.Date(), ".pdf")
ggsave(rmse_spatial, filename = outname, width = 8, height = 8, device='pdf', dpi=320)

