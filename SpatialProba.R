library(roxygen2)

#'function to generate the spatial layer with true probability of presece of a VS
#' @param coefs: a named vector of regression parameters. Names must match those of the env layers (except for intercept, and quadratic terms).
#' params for quadratic terms must have prefix 'quadr_' (followed by the name of the variable, e.g. quadr_bio1)
#' @param env.rast: a SpatRaster objects with environmental layers to generate spatial layer of probabilities
#' @param quadr_term: named vector with names of coefs for which a quadratic term is specified (without prefix 'quadr_')
#' @param marginalPlots: logical, if TRUE returns marginal plot

# myCoeff<-c(intercept = 1, bio1 = 1.6, quad_bio1= -0.10, bio12= 0.0008)
# envData.t<-terra::rast(envData)

# coefs = Fagus.coefs; env.rast =BioData;  quadr_term = "bio1"; marginalPlots=TRUE
# coefs = myCoeff; env.rast =envData; quadr_term = NULL; marginalPlots=TRUE

SpatialProba <- function(coefs=NULL, env.rast=NULL, quadr_term = NULL, marginalPlots=TRUE) {
  #check if names(coefs) is null
  if(isTRUE(is.null(names(coefs)))) stop("coefs must be a named vector")
  #check if the input env.rast is a SpatRaster
  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
  } 
  #check that if quadr_term is not null a 'quad_' name is in coefs
  if(isTRUE(!is.null(quadr_term))) {
    if(!isTRUE(any(grepl("quad", x = names(coefs))))) stop("quadr term specified, but no term with prefix 'quad' found in coefs")
  } 
  #get names of predictors (excluding intercept and names of quadratic terms, if any)
  if(isTRUE(!is.null(quadr_term))) {
    coefs_names <- names(coefs[!grepl("quad|intercept", x = names(coefs))]) 
  } else {
    coefs_names <- names(coefs[!grepl("intercept", x = names(coefs))]) 
  }
  #check if all coefs_names are in names(env.rast)
  if(isTRUE(!all(coefs_names %in% names(env.rast)))) stop("not all names in coefs are found in env.rast")
  #extract intercept
  mod_intrcpt <- coefs[["intercept"]]
  #subset env.rast if needed
  if(isTRUE(length(coefs_names) < terra::nlyr(env.rast))) env.rast <- env.rast[[coefs_names]]
  #coerce env.rast to a data.frame
  env_df <- terra::as.data.frame(env.rast, na.rm = TRUE, xy = TRUE)
  #get coords
  env_coords <- as.matrix(env_df[c("x", "y")])
  #get rid of coords
  env_df <- env_df[-c(1, 2)]
  #add quadr_term(s) (if any)
  if(!is.null(quadr_term)) {
    message(paste("Adding quadratic term for:", paste(quadr_term, collapse = " ")))
    quad_df <- data.frame(lapply(quadr_term, function(coef_nm) (env_df[[coef_nm]])^2))
    colnames(quad_df) <- paste0("quad_", quadr_term)
    env_df <- cbind(env_df, quad_df)
  }
  #coerce env_df to a matrix
  env_mat <- as.matrix(cbind("intercept" = 1, env_df))
  #re-order cols in env_mat to match coefs names
  env_mat <- env_mat[, names(coefs)]
  #check
  if(isTRUE(!all(colnames(env_mat) == names(coefs)))) stop("Names do not match between environmental matrix and coef vector")
  #compute probabilities (logit scale)
  proba_link <- env_mat%*%unname(coefs)
  #back-transform to proba (response) scale
  proba_resp <- plogis(proba_link)
  #rasterize to get the spatial layer
  spatial_proba <- terra::rasterize(x = env_coords, y = env.rast[[1]], values = proba_resp)
  names(spatial_proba) <- "TrueProba"
  
  #compute marginal effects: assuming only additive models now 
  if(isTRUE(marginalPlots)) {
    plotOut <- list()
    if(isTRUE(!is.null(quadr_term))) { ##BIG CHECK
      #compute marginal effect for the quadratic term
      env_mat.tmp<-env_mat
      #fix the non quadratic predictors taking the mean
      env_mat.tmp[ ,coefs_names[!coefs_names %in% c(quadr_term, paste0("quad_", quadr_term))]] <- mean( env_mat.tmp[ ,coefs_names[!coefs_names %in% c(quadr_term, paste0("quad_", quadr_term))]])
      proba_link <- env_mat.tmp%*%unname(coefs)
      proba_resp <- plogis(proba_link)
      marg.eff<-data.frame(quadr_term=env_mat.tmp[,quadr_term], PA=proba_resp)
      colnames(marg.eff)<-c(quadr_term, "PA")
      
      # plot
      p <- ggplot2::ggplot(data = marg.eff, ggplot2::aes(x = marg.eff[, 1], y = marg.eff[, 2])) +
        ggplot2::geom_line() +
        ggplot2::labs(x=quadr_term, y="PA") +
        ggplot2::theme_classic()
      
      plotOut[[quadr_term]]<-p 
    }
  
    #compute marginal effect for the others term
    env_mat.tmp <- env_mat
    
    #fix the other predictors taking the mean
    if(isTRUE(!is.null(quadr_term))) {
      env_mat.tmp <- cbind.data.frame(env_mat.tmp[ , which(!colnames(env_mat.tmp) %in% c(quadr_term, paste0("quad_", quadr_term)))],  
                     data.frame(t(colMeans(env_mat.tmp[ , c(quadr_term, paste0("quad_", quadr_term))]))))
    } else {
      
      env_mat.tmp <- cbind.data.frame(env_mat.tmp[ , which(!colnames(env_mat.tmp) %in% c(quadr_term, paste0("quad_", quadr_term)))],  
                                      data.frame(t(colMeans(env_mat.tmp[ , c(quadr_term, paste0("quad_", quadr_term))]))))
    }
    
    #re-order cols in env_mat to match coefs names
    env_mat.tmp <- as.matrix(env_mat.tmp[, names(coefs)])
    proba_link <- env_mat.tmp%*%unname(coefs)
    proba_resp <- plogis(proba_link)
    myname <- names(coefs[which(!colnames(env_mat.tmp) %in% c(quadr_term, paste0("quad_", quadr_term), "intercept"))])
    marg.eff<-data.frame(x=env_mat.tmp[, myname], PA=proba_resp)
    colnames(marg.eff)<-c(myname, "PA")
    # plot
    p <-ggplot2::ggplot(data = marg.eff, ggplot2::aes(x = marg.eff[, 1], y = marg.eff[, 2])) +
      ggplot2::geom_line() +
      ggplot2::labs(x=myname, y="PA") +
      ggplot2::theme_classic()
    plotOut[[myname]] <- p
    
    plotOut<-ggpubr::ggarrange(plotlist=plotOut)
    
    return(list(rast=spatial_proba, margEff = plotOut))

    } else {
    return(spatial_proba)
  }
  

  
}

# debugonce(SpatialProba)
# SpatialProba(coefs = Fagus.coefs, env.rast =BioData,  quadr_term = "bio1", marginalPlots=TRUE)
# 
# myCoeff<-c(intercept = 5, bio1 = 0.7, bio12= -0.2)
# plot(SpatialProba(coefs = myCoeff, env.rast = BioData))
