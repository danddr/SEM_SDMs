library(ggplot2)
library(dplyr)
outdir <- "/outputs20230824/figures/"
outsim <- "/outputs20230824/"


##---- 1 violin plot of the the GOFs ----
myGOF <- list.files(path = outsim, pattern = "GOF", full.names = TRUE)
myGOF <- lapply(myGOF, readRDS)
myGOF <- do.call(bind_rows, myGOF)
myGOF

p <-   myGOF %>%
    dplyr::select( -fold, -threshold) %>% 
    tidyr::pivot_longer(!c(Model, geog.prev, percTrain))  %>% 
    mutate(geog.prev=as.factor(geog.prev), 
           value=ifelse(value<0, 0, value),
           Model=factor(Model, levels = c("SEM", "GLM"))) %>% 
    dplyr::filter(name %in% c("AUC", "sensitivity", "specificity", "R2", "RMSE","TSS" )) %>%
    dplyr::mutate(name=dplyr::recode(name, "sensitivity"="Sensitivity", "specificity"="Specificity"),
                  name=factor(name, levels = c("R2", "RMSE",  "AUC", "Sensitivity", "Specificity", "TSS"))) %>%
    ggplot(aes(geog.prev, value, col=Model))+
    geom_violin()+
    stat_summary(fun = median, geom = "point",
                 size = 2, position=position_dodge(width = 0.9)) +
    ylim(0,1)+
    scale_color_manual(breaks = c( "SEM", "GLM"),
                       values=c(  "#0072B2",  "#FB8622" ))+
    facet_grid(~name)+
    labs(x="Geographic prevalence", y="Predictive accuracy metrics value", color="Model")+
    theme_bw()+
    theme(legend.background=element_blank(),
          panel.grid = element_blank(),
          legend.position = 'bottom',
          text = element_text(size=16), 
          strip.text = element_text(size=16),
          legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=14),
          legend.key.size = unit(1, 'cm'))

p
outname <- paste(outdir, "Herb_gof_", Sys.Date(),".png", sep="")
ggsave(p, filename = outname, width = 16, height = 8, device='png', dpi=320)


##---- 2 violin plot of the models' coefficients ----
### 2.1 process SEM effects ----
myList <- list.files(path = outsim, pattern = "semEff", full.names = TRUE)
myPercTrain <- 0.8

semEffOut<-list()
for(y in 1:length(myPercTrain)){
  # y=1
  myOut<-myList[grepl(paste0("percTrain", myPercTrain[y]), myList)]
  myOut<-lapply(myOut, readRDS)
  
  tmp<-list()
  for(i in 1:length(myOut)){
    # i=1
    tmp[[i]]<-rbind.data.frame(cbind.data.frame(
      Effects=c(rep('Direct',3), rep('Indirect',3), rep('Total',4), rep('Mediators',2)), 
      Predictors= c("bio1", "bio12", "herb.germ", "bio1", "bio12", "tree.pa", "bio1", "bio12", "tree.pa", "herb.germ",  "tree.pa", "herb.germ"),
      Estimate=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.25`$Summary$herb.pa$Effect))), 
      lCI=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.25`$Summary$herb.pa$`Lower CI`))),
      uCI=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.25`$Summary$herb.pa$`Upper CI`))), 
      geog.prev="0.25", 
      fold=i, 
      Model="SEM", 
      percTrain=myPercTrain[y]), 
      
      cbind.data.frame(
        Effects=c(rep('Direct',3), rep('Indirect',3), rep('Total',4), rep('Mediators',2)), 
        Predictors= c("bio1", "bio12", "herb.germ", "bio1", "bio12", "tree.pa", "bio1", "bio12", "tree.pa", "herb.germ",  "tree.pa", "herb.germ"),
        Estimate=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.5`$Summary$herb.pa$Effect))), 
        lCI=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.5`$Summary$herb.pa$`Lower CI`))),
        uCI=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.5`$Summary$herb.pa$`Upper CI`))), 
        geog.prev="0.5", 
        fold=i, 
        Model="SEM", 
        percTrain=myPercTrain[y]),
      
      cbind.data.frame(
        Effects=c(rep('Direct',3), rep('Indirect',3), rep('Total',4), rep('Mediators',2)), 
        Predictors= c("bio1", "bio12", "herb.germ", "bio1", "bio12", "tree.pa", "bio1", "bio12", "tree.pa", "herb.germ",  "tree.pa", "herb.germ"),
        Estimate=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.75`$Summary$herb.pa$Effect))), 
        lCI=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.75`$Summary$herb.pa$`Lower CI`))),
        uCI=as.numeric(na.omit(as.numeric(myOut[[i]]$`0.75`$Summary$herb.pa$`Upper CI`))),
        geog.prev="0.75", 
        fold=i, 
        Model="SEM",
        percTrain=myPercTrain[y])
    )
    
  }
  
  semEffOut[[y]]<-do.call(rbind.data.frame, tmp)
  
}

length(semEffOut)
str(semEffOut)
head(semEffOut)
semEffOut<-do.call(rbind.data.frame, semEffOut)

### 2.2 process GLM effects ----
myList <- list.files(path = outsim, pattern = "prev", full.names = TRUE)
myPercTrain<-c("0.8")
geog.prev<-c("0.25", "0.5", "0.75")
glmEffOut<-list()

for(y in 1:length(myPercTrain)){
  # y=1
  message("Processing percentage training", myPercTrain[y])
  mysubList<-myList[grepl(paste0("percTrain", myPercTrain[y]), myList)]
  tmp<-list()
  for(k in 1:length(geog.prev)){
    # k=2
    message("Processing geographic prevalence", geog.prev[k])
    myOut<-mysubList[grepl(paste0("prev", geog.prev[k]) , mysubList)] 
    myOut<-lapply(myOut, readRDS)
    
    tmp2<-list()
    for(i in 1:length(myOut)){
      # i=1
      tmp2[[i]]<- cbind.data.frame(
        Effects=c(rep('Total',4)), 
        Predictors= c("Intercept",  "herb.germ", "bio1", "bio12"),
        Estimate=myOut[[i]]$coefficients, 
        lCI=confint(myOut[[i]])[,1],
        uCI=confint(myOut[[i]])[,2], 
        geog.prev=geog.prev[k], 
        fold=i, 
        Model="GLM", 
        percTrain=myPercTrain[y]) 
    }
    
    tmp[[k]]<-do.call(rbind.data.frame, tmp2)
  }
  
  glmEffOut[[y]]<-do.call(rbind.data.frame, tmp)
  
}
glmEffOut<-do.call(rbind.data.frame, glmEffOut)
glmEffOut$percTrain<-as.numeric(glmEffOut$percTrain)

### 2.3 Plotting ----
coeffOut <-  bind_rows(glmEffOut,
                      subset(semEffOut, Effects=="Direct")) %>%  
  # dplyr::select(-fold) %>% 
  tidyr::pivot_longer(!c(geog.prev,  Model, percTrain, Effects,  Predictors, fold)) %>% 
  mutate(geog.prev=paste("Geographic prevalence", geog.prev)  
  ) %>% 
  filter( name=="Estimate",
         Predictors %in% c("bio1", "bio12", "herb.germ")) %>%   
  mutate(Predictors=recode(Predictors,  "herb.germ" = "Herb germination rate", 
                           "bio1" = "BIO1",  
                           "bio12" = "BIO12")) 


coeffOut %>% 
  group_by(Model, geog.prev, Predictors) %>% 
  summarise(med =median(value))

p <- coeffOut %>% 
        mutate(Predictors=as.factor(Predictors),
               Predictors=forcats::fct_relevel(Predictors, c("Herb germination rate", "BIO12", "BIO1")),
               Model=as.factor(Model),
               Model=factor(Model, levels=c("SEM", "GLM"))) %>% 
        ggplot(aes(Predictors, value, col=Model)) +
        geom_boxplot()+
        geom_hline(yintercept = 0, lty="dashed", col="darkgrey")+
        scale_color_manual(breaks = c( "SEM", "GLM" ),
                           values=c("#0072B2", "#D55E00"   ))+
        labs(x="", y="Coefficients")+
        # ylim(0,2)+
        facet_grid(~geog.prev, scales = "free")+
        coord_flip()+
        theme_bw()+
        theme(legend.background=element_blank(),
              panel.grid = element_blank(),
              legend.position = 'bottom',
              text = element_text(size=16), 
              strip.text = element_text(size=16),
              legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=16),
              legend.key.size = unit(1, 'cm'))

p
outname <- paste(outdir, "Herb_coeff_", Sys.Date(),".png", sep="")
ggsave(p, filename = outname, width = 16, height = 8, device='png', dpi=320)

