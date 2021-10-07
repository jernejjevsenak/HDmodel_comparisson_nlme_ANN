# This function is based on lmfor R package 
# (https://cran.r-project.org/web/packages/lmfor/index.html)
# All credits go to the authors of lmfor 

fit_nlme_two_levels <- function(d, h, plot, plot_PCA, modelName = "naslund", 
                                level = 1, na.omit = TRUE,
                                start=NA, bh=1.3, control=list(),random=NA,
                                SubModels=NA, vfstart=0){
  
if (!is.na(list(random))) nranp <- 1 # number of random parameters
  
  
  varf<-as.numeric(TRUE)
  
  if (na.omit) {
    sel<-is.na(h*d)
    d<-d[!sel]
    plot<-plot[!sel]
    h<-h[!sel]
  }
  
  dmean<-tapply(d,plot,mean)
  dmean<-dmean[match(plot,names(dmean))]
  
  dsd<-tapply(d,plot,sd)
  dsd<-dsd[match(plot,names(dsd))]
  dstd<-(d-dmean)/dsd

  if (varf==1)  w<-d  else if (varf==2) w<-pmax(1,dstd+3)
  
    # Nonlinear model fitting
    if (is.na(start[1])) {
      start<-eval(call(paste("startHD",modelName,sep=""),d=d,h=h,bh=bh))
      if (!is.na(SubModels)[1]) {
        start2<-c()
        for (i in 1:length(SubModels)) {
          start2<-c(start2,start[i])
        }
        start<-start2
      }
    }
    
        if (modelName == "curtis"){
          
          mod<-nlme(h ~ HDcurtis(d, a, b, bh = 1.3),
                    fixed = list(a~1, b ~1),
                    random = a + b ~ 1|plot_PCA/plot,
                    start=start,
                    weights=varPower(vfstart,~w), 
                    control=list(maxIter = 500, msmaxIter = 500))

        } else if (modelName == "naslund"){
          
          mod<-nlme(h ~ HDnaslund(d, a, b, bh = 1.3),
                    fixed = list(a~1, b ~1),
                    random = a + b ~ 1|plot_PCA/plot,
                    start=start,
                    weights=varPower(vfstart,~w), 
                    control=list(maxIter = 500, msmaxIter = 500))
          
        } else if (modelName == "wykoff"){
          
          mod<-nlme(h ~ HDwykoff(d, a, b, bh = 1.3),
                    fixed = list(a~1, b ~1),
                    random = a + b ~ 1|plot_PCA/plot,
                    start=start,
                    weights=varPower(vfstart,~w), 
                    control=list(maxIter = 500, msmaxIter = 500))
          
        } else if (modelName == "michailoff"){
          
          mod<-nlme(h ~ HDmichailoff(d, a, b, bh = 1.3),
                    fixed = list(a~1, b ~1),
                    random = a + b ~ 1|plot_PCA/plot,
                    start=start,
                    weights=varPower(vfstart,~w), 
                    control=list(maxIter = 500, msmaxIter = 500))
          
        } else {
          
          stop("Wrong model!")
          
        }
    
  mod
}




