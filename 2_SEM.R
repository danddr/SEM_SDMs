library(dismo)
library(piecewiseSEM)
library(semEff)
library(pROC)
library(ecospat)
library(Metrics)

outdir <- "/outputs20230824/figures/"
outsim <- "/outputs20230824/"

geog.prev <- c(0.25, 0.5, 0.75)
loc <- readRDS("/outputs20230824/samplinLoc_2023-08-24.RDS")
myGrid <- readRDS("/outputs20230824/crossValGrid_2023-08-24.RDS")
pred.stack <- readRDS("/outputs20230824/stacked_predictors_2023-08-24.RDS")

train.perc <- 0.8
myFolds <- 10

for (f in 1:myFolds) {
  # f=1
  message("Processing fold", f)
  grid2sample <- sample(values(myGrid), train.perc*length(values(myGrid)))
  randpts <- loc
  randpts$train <- ifelse(randpts$train %in% grid2sample , 1, 0)

  ### 1. Extract variables and scaling them 
  dat <- lapply(pred.stack, function(x){cbind.data.frame(train=randpts$train,  raster::extract(x, randpts, df = TRUE , cellnumbers=FALSE, exact=TRUE))})
  dat.sc <- lapply(dat, function(x){data.frame(train=x[,1], scale(x[,3:5], center = TRUE, scale = TRUE), x[, 6:8])})
  names(dat.sc) <- geog.prev
  
  ### 2. Assembling the SEM
  vs.semList <- list()
  for(i in 1:length(dat.sc)) {
    # i=1
    tmp.df <- subset(dat.sc[[i]], train==1) # select training dataset
    vs.sem <- list(tree.pa = glm(tree.pa ~ bio1 + bio12, family = binomial,
                                 data = tmp.df),
                   herb.germ = glm(herb.germ ~ tree.pa + bio1 + bio12,
                                   data = tmp.df),
                   herb.pa = glm(herb.pa ~ herb.germ + bio1 + bio12 , family = binomial, 
                                 data = tmp.df)               
    )
    vs.semList[[i]] <- vs.sem
  }
  names(vs.semList) <- paste0("geog.prev", geog.prev)
  
  # Network visualization
  # plot(piecewiseSEM::as.psem(vs.semList$geog.prev0.75, Class = "psem"))
  outname <- paste0(outsim, "vs.semList", "_percTrain",train.perc, "_fold", f,"_", Sys.Date(), ".RDS" )
  saveRDS(vs.semList, outname)
  
  ### 3. Investigate SEM quality
  # Basis set = the smallest possible set of independence claims from a graph
  dSep.tab <- list()
  for(i in 1:length(geog.prev)){
    # i=1
    mytmp <- piecewiseSEM::dSep(piecewiseSEM::as.psem(vs.semList[[i]]))  
    mytmp$fisherC <-fisherC(mytmp)$Fisher.C
    mytmp$fisherC.df <-fisherC(mytmp)$df
    mytmp$fisherC.P.val <-fisherC(mytmp)$P.Value
    mytmp$geog.prev <- geog.prev[i]
    dSep.tab[[i]] <- mytmp
  }
  dSep.tab <- do.call(rbind.data.frame, dSep.tab)
  outname <- paste0(outsim, "dSeptab", "_percTrain",train.perc, "_fold", f,"_", Sys.Date(), ".RDS" )
  saveRDS(dSep.tab, outname)
  
  ### 4. Store summary
  sem.sum <- list()
  for(i in 1:length(geog.prev)){
    # i=1
    mytmp <- piecewiseSEM::as.psem(vs.semList[[i]], Class = "psem")  
    mytmp <- summary(mytmp, .progressBar = F, conserve = TRUE)
    sem.sum[[i]] <- mytmp
  }
  names(sem.sum) <- geog.prev
  # sem.sum$`0.25`$coefficients
  outname <- paste0(outsim, "semSummary", "_percTrain",train.perc, "_fold", f,"_", Sys.Date(), ".RDS" )
  saveRDS(sem.sum, outname)
  
  ### 5. Calculate the effects
  # To calculate the effects, we first bootstrap and save the standardized direct effects. 
  # Provide a sufficient number of bootstraps samples, since a small number may produce NAs (see [this source](https://stats.stackexchange.com/questions/37918/why-is-the-error-estimated-adjustment-a-is-na-generated-from-r-boot-package)).
  sem.eff <- list()
  for(i in 1:length(geog.prev)){
    # i=1
    message(paste("Effects: Processing geographic prevalence",  geog.prev[i]))
    mytmp <- piecewiseSEM::as.psem(vs.semList[[i]], Class = "psem")  
    mytmp <- bootEff(mytmp,  R = 1000,  seed = 13, parallel = "multicore", cl = 5)
    # Using the bootstrapped estimates we can now calculate direct, indirect, and total effects:
    mytmp <- semEff(mytmp)
    sem.eff[[i]] <- mytmp
  }
  names(sem.eff) <- geog.prev
  outname <- paste0(outsim, "semEffects", "_percTrain",train.perc, "_fold", f,"_", Sys.Date(),".RDS" )
  saveRDS(sem.eff, outname)
  
  ### 6. SEM predictions on the testing dataset
  pred.dat <- list() 
  for(i in 1:length(geog.prev)){
    # i=1
    tmp <- subset(dat.sc[[i]], train==0)
    tmp$geog.prev <- geog.prev[i]
    pred.dat[[i]]<-tmp
  }  
  names(pred.dat)<-geog.prev  
  
  for(i in 1:length(geog.prev)){
    # i=1
    message(paste("SEM predictions: Processing geographic prevalence",  geog.prev[i]))
    # Extract the (total) effects of each variable to use for prediction
    tot <- getTotEff(sem.eff[[i]], "herb.pa")
    tot.b <- getTotEff(sem.eff[[i]], "herb.pa", type = "boot")
    
    # SEM prediction
    mod <- vs.semList[[i]]$herb.pa
    test_simple_pred <- predEff(mod, 
                                newdata =  pred.dat[[i]], 
                                effects = tot,
                                # eff.boot = tot.b,
                                type = "response")
    pred.dat[[i]]$pred.herb.suit <- test_simple_pred #test_simple_pred$fit
    # pred.dat[[i]]$pred.herb.suit.se <- test_simple_pred$se.fit
    # pred.dat[[i]]$pred.herb.suit.ciLow <- test_simple_pred$ci.lower
    # pred.dat[[i]]$pred.herb.suit.ciUp <- test_simple_pred$ci.upper
  }
  
  ### 7. GLM predictions on the testing dataset
  for(i in 1:length(geog.prev)){
    # i=1
    message(paste("GLM predictions: Processing geographic prevalence",  geog.prev[i]))
    tmp.df <- subset(dat.sc[[i]], train==1) # select training dataset
    myglm <- glm(herb.pa~  herb.germ + bio1 + bio12 , data = tmp.df, family = binomial) 
    glm.pred <- predict(myglm, newdata=pred.dat[[i]], type="response")
    pred.dat[[i]]$pred.herb.suit.GLM <- glm.pred
    pred.dat[[i]]$R2.GLM <- round(modEvA::RsqGLM(myglm, plot=FALSE)$Nagelkerke,3)
    
    outname<-paste0(outsim, "glmEffects_", "prev", geog.prev[i], "_percTrain",train.perc, "_fold", f,"_", Sys.Date(), ".RDS" )
    saveRDS(myglm, outname)
    
  }
  
  outname<-paste0(outsim, "semPred", "_percTrain",train.perc, "_fold", f,"_", Sys.Date(), ".RDS" )
  saveRDS(pred.dat, outname)
  
  ### 8 Evaluate the predictions
  myOut<-list()
  for(i in 1:length(geog.prev)){
    # i=1
    tmp<-pred.dat[[i]][, c("geog.prev","herb.pa"  , "obs.herb.suit",  "pred.herb.suit", "pred.herb.suit.GLM", "R2.GLM")] #"pred.herb.suit.GLM.simple", "R2.GLM.simple"
    
    pROC_obj_SEM <- roc(tmp$herb.pa, tmp$pred.herb.suit,
                        smoothed = TRUE,
                        # arguments for ci
                        ci = TRUE, ci.alpha = 0.9, stratified = FALSE,
                        # arguments for plot
                        plot = FALSE, auc.polygon = TRUE, max.auc.polygon = TRUE,
                        grid = TRUE, print.auc = TRUE, show.thres = TRUE)
    
    sem.out<-rbind(data.frame("Model" = "SEM",
                              "geog.prev"=unique(tmp$geog.prev),  
                              "AUC" = round(pROC::auc(pROC_obj_SEM),3), 
                              round(coords(pROC_obj_SEM, "best",  transpose = FALSE,best.method="closest.topleft"),3),
                              "CBI" = ecospat.boyce(fit = tmp$pred.herb.suit,
                                                    obs = tmp[which(tmp$herb.pa==1),"pred.herb.suit"],
                                                    nclass = 0,
                                                    window.w = "default", res = 100, PEplot = FALSE)$cor,
                              "TSS" = round(ecospat.max.tss(Pred = tmp$pred.herb.suit,
                                                            Sp.occ = tmp$herb.pa)$max.TSS,3),
                              "RMSE" = round(Metrics::rmse(tmp$obs.herb.suit, tmp$pred.herb.suit),3),
                              "logLoss" = round(Metrics::logLoss(tmp$obs.herb.suit, tmp$pred.herb.suit),3),
                              "R2" = sem.sum[[i]]$R2$R.squared[3])
    )
    
    pROC_obj_GLM <- roc(tmp$herb.pa, tmp$pred.herb.suit.GLM,
                        smoothed = TRUE,
                        # arguments for ci
                        ci = TRUE, ci.alpha = 0.9, stratified = FALSE,
                        # arguments for plot
                        plot = FALSE, auc.polygon = TRUE, max.auc.polygon = TRUE,
                        grid = TRUE, print.auc = TRUE, show.thres = TRUE)
    
    glm.out<-rbind(data.frame("Model" = "GLM",
                              "geog.prev"= unique(tmp$geog.prev),  
                              "AUC" = round(pROC::auc(pROC_obj_GLM),3), 
                              round(coords(pROC_obj_GLM, "best",  transpose = FALSE,best.method="closest.topleft"),3),
                              "CBI" = ecospat.boyce(fit = tmp$pred.herb.suit.GLM,
                                                    obs = tmp[which(tmp$herb.pa==1),"pred.herb.suit.GLM"],
                                                    nclass = 0,
                                                    window.w = "default", res = 100, PEplot = FALSE)$cor,
                              "TSS" = round(ecospat.max.tss(Pred = tmp$pred.herb.suit.GLM,
                                                            Sp.occ = tmp$herb.pa)$max.TSS,3),
                              "RMSE" = round(Metrics::rmse(tmp$obs.herb.suit, tmp$pred.herb.suit.GLM),3),
                              "logLoss" = round(Metrics::logLoss(tmp$obs.herb.suit, tmp$pred.herb.suit.GLM),3),
                              "R2" = unique(tmp$R2.GLM))
    )
  
    myOut[[i]]<-rbind.data.frame(sem.out, glm.out) # rf.out, glm.out.simple
  }
  
  myOut<-do.call(rbind.data.frame, myOut)
  myOut$percTrain<-train.perc
  myOut$fold<-f
  
  outname<-paste0(outsim, "GOF", "_percTrain",train.perc, "_fold", f,"_", Sys.Date(), ".RDS" )
  saveRDS(myOut, outname)
}
