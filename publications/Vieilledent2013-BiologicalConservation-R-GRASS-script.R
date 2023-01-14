################################################################################
#///////////////////////////////////////////////////////////////////////////////
#
# R/GRASS script
#
# Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
# June, 3rd 2013
#
# Vieilledent G., Cornu C., Cuni-Sanchez A., Leong Pock-Tsy J.-M. and
# Danthu P. 2013. Vulnerability of baobab species to climate change
# and effectiveness of the protected area network in Madagascar:
# towards new conservation priorities. Biological Conservation.
#
#///////////////////////////////////////////////////////////////////////////////
################################################################################

# Library
library(sp)
library(spgrass6)

# Including other mapsets
system("g.mapsets addmapset=climpres")
system("g.mapsets addmapset=clim2050")
system("g.mapsets addmapset=clim2080")

# Specifying the region resolution
system("g.region vect=Geol res=00:01")

# Show metadata
str(gmeta6())

# Examples using grass commands with function system()
system("g.list vect")
system("g.list rast")

#===================================================
# PCA to determine non-correlated climatic variables
#===================================================

#==============================================================
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (bio2/bio7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (bio5-bio6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month 
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter
#==============================================================

#== Creating a grid of vector point
system("v.mkgrid -p map=DataPCA grid=820,435 --o")
system("v.out.ascii DataPCA > DataPCA_temp.txt")
#== y coordinate must be inverted in the DataPCA_temp.txt  file, we use awk...
system("awk 'BEGIN{FS=\"|\";OFS=\"|\"} {x[NR]=$1;y[NR]=$2} END{i=1; while(i<=NR) {print x[i],y[NR+1-i]; i++}}' DataPCA_temp.txt > DataPCA.txt")
system("rm DataPCA_temp.txt")

#== To delete points which are unnecessary...
# ...we create a raster called Maskout indicating in(1) and out(0) raster points
system("r.mapcalc 'Maskout=if(isnull(bio1),0,1)'")
system("r.out.ascii Maskout@ghvi null=-9999 > Maskout.txt")
system("awk 'NR>=7{for (i=1;i<=(NF);i++) {print $i}}' Maskout.txt > Temp.txt")
system("awk '{OFS=\"|\"; getline x < \"Temp.txt\"; print $0,x}' DataPCA.txt > Results.txt")
system("mv Results.txt DataPCA.txt")

#== We add species data to DataPCA dataframe
names.species <- c("Ag","Ap","As")
for (s in names.species) {
  system(paste("r.out.ascii ",s,"_pres_abs null=-9999 > Maskout.txt",sep=""))
  system("awk 'NR>=7{for (i=1;i<=(NF);i++) {print $i}}' Maskout.txt > Temp.txt")
  system("awk '{OFS=\"|\"; getline x < \"Temp.txt\"; print $0,x}' DataPCA.txt > Results.txt")
  system("mv Results.txt DataPCA.txt")
}

#== We add all climatic variables to DataPCA dataframe
names.bio <- paste("bio",c(1:19),sep="")
for (v in names.bio) {
  system(paste("r.out.ascii ",v," null=-9999 > Maskout.txt",sep=""))
  system("awk 'NR>=7{for (i=1;i<=(NF);i++) {print $i}}' Maskout.txt > Temp.txt")
  system("awk '{OFS=\"|\"; getline x < \"Temp.txt\"; print $0,x}' DataPCA.txt > Results.txt")
  system("mv Results.txt DataPCA.txt")
}

#== We add the geologic variable to DataPCA dataframe
system("v.to.rast Geol output=Geol_rast col=RECLASS_ID --o") # Geol is transform into a raster
system("r.out.ascii Geol_rast@ghvi null=-9999 > Maskout.txt")
system("awk 'NR>=7{for (i=1;i<=(NF);i++) {print $i}}' Maskout.txt > Temp.txt")
system("awk '{OFS=\"|\"; getline x < \"Temp.txt\"; print $0,x}' DataPCA.txt > Results.txt")
system("mv Results.txt DataPCA.txt")

#== Cleaning DataPCA
DataPCA <- read.table(file="DataPCA.txt",header=FALSE,sep="|")
names(DataPCA) <- c("Long","Lat","Maskout",names.species,names.bio,"geol")
DataPCAup <- DataPCA[DataPCA$Maskout==1,-c(3)]
DataPCAup[DataPCAup==-9999] <- NA
Vect.Row.notisna <- !is.na(apply(DataPCAup,1,mean)) # Identifying row with the presence of NA
DataPCAup <- DataPCAup[Vect.Row.notisna,]
DataPCAup$ID <- c(1:length(DataPCAup[,1]))

#== Summary DataPCA
summary(DataPCAup)

#== Principal Component Analysis
library(ade4)

PCA <- dudi.pca(df=DataPCAup[,c(6:24)],center=TRUE,scale=TRUE,scannf=FALSE,nf=6)
save(PCA,file="PCA")
PCA$eig 
PC <- PCA$li # The row coordinates, i.e. the principal components

#== Plot
pdf(file="PCA.pdf",width=8,height=8)
par(mfrow=c(2,2),mar=c(2,2,2,0))
#=
x <- barplot(PCA$eig,ylim=c(0,8),main="Eigen values") # To choose the number of axes to keep
axis(1,at=x[c(1:6,seq(8,18,2))],labels=c(1:6,seq(8,18,2)))
#=
s.corcircle(PCA$co,xax=1,yax=2) # Correlation circle to interprete the axis
text(x=c(1.1),y=c(0),labels=c("Axis 1"),srt=c(-90),cex=1.2)
text(x=c(0),y=c(1.1),labels=c("Axis 2"),srt=c(0),cex=1.2)
#=
s.corcircle(PCA$co,xax=3,yax=4) # Correlation circle to interprete the axis
text(x=c(1.1),y=c(0),labels=c("Axis 3"),srt=c(-90),cex=1.2)
text(x=c(0),y=c(1.1),labels=c("Axis 4"),srt=c(0),cex=1.2)
#=
s.corcircle(PCA$co,xax=5,yax=6) # Correlation circle to interprete the axis
text(x=c(1.1),y=c(0),labels=c("Axis 5"),srt=c(-90),cex=1.2)
text(x=c(0),y=c(1.1),labels=c("Axis 6"),srt=c(0),cex=1.2)
dev.off()

Data.present <- cbind(DataPCAup,PC)

#== Supplementary individuals (see: http://pbil.univ-lyon1.fr/R/pdf/qr8.pdf)
# To do later with the future values of the bioclim variables for each cell

#===================================================
# Data-set for each species
#===================================================

Ag.data <- Data.present[,c(26,2,1,3,25,27:32)]
Ap.data <- Data.present[,c(26,2,1,4,25,27:32)]
As.data <- Data.present[,c(26,2,1,5,25,27:32)]
write.table(Ag.data,file="Ag.data.txt",sep="\t",row.names=FALSE,quote=FALSE)
write.table(Ap.data,file="Ap.data.txt",sep="\t",row.names=FALSE,quote=FALSE)
write.table(As.data,file="As.data.txt",sep="\t",row.names=FALSE,quote=FALSE)
write.table(Data.present,file="Data.present.txt",sep="\t",row.names=FALSE,quote=FALSE)

#===================================================
# Inference Species Distribution Model
#===================================================

# Importing data-sets
Data.present <- read.table(file="Data.present.txt",sep="\t",header=TRUE)
Data.present$geol <- as.factor(Data.present$geol)

#=================
# GLM/GAM Loops on
#=================
# (i) Species
# (ii) Statistic model: GLM, GAM
# (iii) Covariates: Complete, Only climatic variables

library(BIOMOD)
names(Data.present)
Initial.State(Response=Data.present[,c(3:5)],Explanatory=Data.present[,c(25,27:32)],
              IndependentResponse=NULL,IndependentExplanatory=NULL)
## ls()
## head(DataBIOMOD)
## str(DataBIOMOD)

#= Models

Models(# Model type
       GLM = TRUE, TypeGLM = "poly", Test = "AIC",
       GBM = FALSE, No.trees = 5000,
       GAM = TRUE, Spline = 3,
       CTA = FALSE, CV.tree = 50,
       ANN = FALSE, CV.ann = 5,
       SRE = FALSE, quant=0.025,
       FDA = FALSE,
       MARS = FALSE,
       RF = FALSE,
       # Cross-validation
       NbRunEval = 5,
       DataSplit = 70,
       # Pseudo-absence
       NbRepPA=1, strategy="random", coor=NULL, distance=0, nb.absences=10000, Yweights = NULL,
       # Diagnostic
       VarImport = TRUE,
       Roc = TRUE,
       Optimized.Threshold.Roc = TRUE,
       Kappa = TRUE,
       TSS = TRUE,
       KeepPredIndependent = FALSE)

#= Diagnostic
ls()
Evaluation.results.Roc
Evaluation.results.TSS
VarImportance
save.image()

#========================
# Maxent Cross-Validation
#========================

species <- c("Ag","Ap","As")

for (sp in 1:length(species)) {
  eval(parse(text=paste("
    
    #####################
    # Preparing the data:

    Col <- c(1,3,2,5,6:11)

    # sample
    sample <- ",species[sp],".data[",species[sp],".data$",species[sp],"==1,Col]
    sample$ID <- \"",species[sp],"\"
    names(sample)[1:3] <- c(\"species\",\"longitude\",\"latitude\")
    write.table(sample,file=\"./maxent_data/swd/sample.csv\",sep=\",\",quote=FALSE,row.names=FALSE)

    # background
    background <- ",species[sp],".data[",species[sp],".data$",species[sp],"==0,Col]
    background$ID <- \"background\"
    names(background)[1:3] <- c(\"species\",\"longitude\",\"latitude\")
    write.table(background,file=\"./maxent_data/swd/background.csv\",sep=\",\",quote=FALSE,row.names=FALSE)

    # projection
    projection <- ",species[sp],".data[,Col]
    projection$ID <- \"background\"
    names(projection)[1:3] <- c(\"species\",\"longitude\",\"latitude\")
    write.table(projection,file=\"./maxent_data/layers/projection.csv\",sep=\",\",quote=FALSE,row.names=FALSE)

    ######################
    # Call to maxent program
    system(\"java -mx512m -jar ./maxent_prog/maxent.jar -r -a nowarnings noprefixes randomseed=true replicates=5 -X 30 replicatetype=subsample outputdirectory=./maxent_data/outputs1/ projectionlayers=./maxent_data/layers/projection.csv samplesfile=./maxent_data/swd/sample.csv environmentallayers=./maxent_data/swd/background.csv -t geol\")",sep="")))
}

#========================
# Maxent full data-set
#========================

species <- c("Ag","Ap","As")

for (sp in 1:length(species)) {
  eval(parse(text=paste("
    
    #####################
    # Preparing the data:

    Col <- c(1,3,2,5,6:11)

    # sample
    sample <- ",species[sp],".data[",species[sp],".data$",species[sp],"==1,Col]
    sample$ID <- \"",species[sp],"\"
    names(sample)[1:3] <- c(\"species\",\"longitude\",\"latitude\")
    write.table(sample,file=\"./maxent_data/swd/sample.csv\",sep=\",\",quote=FALSE,row.names=FALSE)

    # background
    background <- ",species[sp],".data[",species[sp],".data$",species[sp],"==0,Col]
    background$ID <- \"background\"
    names(background)[1:3] <- c(\"species\",\"longitude\",\"latitude\")
    write.table(background,file=\"./maxent_data/swd/background.csv\",sep=\",\",quote=FALSE,row.names=FALSE)

    # projection
    projection <- ",species[sp],".data[,Col]
    projection$ID <- \"background\"
    names(projection)[1:3] <- c(\"species\",\"longitude\",\"latitude\")
    write.table(projection,file=\"./maxent_data/layers/projection.csv\",sep=\",\",quote=FALSE,row.names=FALSE)

    ######################
    # Call to maxent program
    system(\"java -mx512m -jar ./maxent_prog/maxent.jar -r -a nowarnings noprefixes outputdirectory=./maxent_data/outputs2/ projectionlayers=./maxent_data/layers/projection.csv samplesfile=./maxent_data/swd/sample.csv environmentallayers=./maxent_data/swd/background.csv -t geol\")",sep="")))
}

#==================
# Fitted values
#==================

#= Variables
species <- c("Ag","Ap","As")
glmgam <- c("glm","gam")

#= Data
Data.fitted <- Data.present[,c(1:5,25,27:32)]

#==========================================
# Fitted values with GLM and GAM
#==========================================

#= Importing models
load("./models/Ag_GAM_PA1")
gam.Ag <- Ag_GAM_PA1
load("./models/Ag_GLM_PA1")
glm.Ag <- Ag_GLM_PA1
load("./models/Ap_GAM_PA1")
gam.Ap <- Ap_GAM_PA1
load("./models/Ap_GLM_PA1")
glm.Ap <- Ap_GLM_PA1
load("./models/As_GAM_PA1")
gam.As <- As_GAM_PA1
load("./models/As_GLM_PA1")
glm.As <- As_GLM_PA1
rm(list=ls(pattern="^A[gps]_G*"))

#= Function
fit.glmgam <- function() {
  cat("GLM and GAM fit...\n")
  # Explicative variables
  X <- Data.present[,c(25,27:32)]
  for (sp in 1:length(species)) {
    cat(paste(species[sp]," ",sep=""))
    for (g in 1:length(glmgam)) {
      cat(paste(glmgam[g]," ",sep=""))
      # Computing fitted values
      eval(parse(text=paste("
           Response <- predict.",glmgam[g],"(object=",glmgam[g],".",species[sp],",newdata=X,type=\"response\")
           Data.fitted$Fit.",species[sp],".",glmgam[g]," <- Response",sep="")))
    }
    cat("\n")
  }
  return(Data.fitted)
}

#==========================================
# Fitted values with Maxent
#==========================================

#= Function
fit.maxent <- function() {
  cat("Maxent fit...\n")
  for (sp in 1:length(species)) {
    cat(paste(species[sp]," maxent ",sep=""))
    # projection
    projection <- Data.present[,c(26,1:2,25,27:32)]
    projection$ID <- "background"
    names(projection) <- c("species","longitude","latitude","geol",paste("Axis",c(1:6),sep=""))
    write.table(projection,file="./maxent_data/layers/projection.csv",sep=",",quote=FALSE,row.names=FALSE)
    
    # Call to maxent program
    system(paste("java -cp ./maxent_prog/maxent.jar density.Project ./maxent_data/outputs2/",species[sp],".lambdas ./maxent_data/layers/projection.csv ./maxent_data/fitted/",species[sp],"_fitted asc -r",sep=""))

    # Importing Maxent predictions
    eval(parse(text=paste("
         Data.fitted$Fit.",species[sp],".maxent <- scan(file=\"./maxent_data/fitted/",species[sp],"_fitted.asc\",skip=6)",sep="")))
  }
  return(Data.fitted)
}

#===============================
# Simulation RUNS
#===============================

Data.fitted <- fit.glmgam()
Data.fitted <- fit.maxent()

#==================
# Model performance
#==================

#= Var
species <- c("Ag","Ap","As")
n.sp <- length(species)
model.stat <- c("glm","gam","maxent")
n.ms <- length(model.stat)
n.th <- 100
th <- seq(from=0,to=1,length.out=n.th)

#= Matrix to hold the results
Mat.Index <- as.data.frame(matrix(NA,nrow=n.sp*n.ms,ncol=7))
names(Mat.Index) <- c("Sp","Mod","AUC","TSS","MST","th.MST","LPT")
Mat.Index$Sp <- rep(species,each=n.ms)
Mat.Index$Mod <- rep(model.stat,n.sp)

#= AUC
# Def: Probability that the model will rank a randomly
# chosen species presence site higher than a
# randomly chosen absence site (Pearce and
# Ferrier, 2000)
library(ROCR)
for (sp in 1:n.sp) {
  for (ms in 1:n.ms) {
    eval(parse(text=paste("
    Mat.Index$AUC[",(sp-1)*n.ms+ms,"] <- round(performance(prediction(Data.fitted$Fit.",species[sp],".",model.stat[ms],",Data.fitted$",species[sp],"),\"auc\")@y.values[[1]],digits=4)",sep="")))
  }
}

## #= Kappa
## library(vcd)
## Kappa.calc <- function(obs,pred,th) {
##   # Threshold
##   pred01 <- rep(0,length(pred))
##   pred01[pred>=th] <- 1
  
##   # Confusion matrix
##   n11 <- sum(obs==1&pred01==1)
##   n01 <- sum(obs==0&pred01==1)
##   n00 <- sum(obs==0&pred01==0)
##   n10 <- sum(obs==1&pred01==0)
##   Mat <- matrix(c(n00,n01,n10,n11),ncol=2,byrow=TRUE)
##   # Kappa of Cohen
##   return(Kappa(Mat))
## }
## for (sp in 1:length(species)) {
##   for (ms in 1:length(model.stat)) {
##     for (c in 1:length(cov)) {
##       eval(parse(text=paste("
##          K.",species[sp],".",model.stat[ms],".",cov[c],".vect <- rep(0,n.th)
##          for (i in 1:n.th) {
##            K.",species[sp],".",model.stat[ms],".",cov[c],".vect[i] <- Kappa.calc(obs=Data.fitted$",species[sp],",pred=Data.fitted$Fit.",species[sp],".",model.stat[ms],".",cov[c],",th=th[i])$Unweighted[1]
##          }
##          Mat.Index$Kappa[",(sp-1)*(n.ms*n.cov)+(ms-1)*n.cov+c,"] <- max(K.",species[sp],".",model.stat[ms],".",cov[c],".vect) 
##          Mat.Index$th.K[",(sp-1)*(n.ms*n.cov)+(ms-1)*n.cov+c,"] <- th[which.max(K.",species[sp],".",model.stat[ms],".",cov[c],".vect)]",sep="")))
##     }
##   }
## }

#= Maximized Sum Threshold
MST.calc <- function(obs,pred,th) {
  # Threshold
  pred01 <- rep(0,length(pred))
  pred01[pred>=th] <- 1
  
  # Confusion matrix
  n11 <- sum(obs==1&pred01==1)
  n01 <- sum(obs==0&pred01==1)
  n00 <- sum(obs==0&pred01==0)
  n10 <- sum(obs==1&pred01==0)
  sensitivity <- n11/(n11+n10)
  specificity <- n00/(n00+n01)
  return(sensitivity+specificity)
}
#== When the MST is maximal for different values of threshold, we take the mean of these values
th.ghvi <- function(x) {
  X <- round(x,2)
  max.X <- max(X)
  th <- mean(th[which(X==max.X)])
  return(th)
}
#=
for (sp in 1:length(species)) {
  for (ms in 1:length(model.stat)) {
    eval(parse(text=paste("
         MST.",species[sp],".",model.stat[ms],".vect <- rep(0,n.th)
         for (i in 1:n.th) {
           MST.",species[sp],".",model.stat[ms],".vect[i] <- MST.calc(obs=Data.fitted$",species[sp],",pred=Data.fitted$Fit.",species[sp],".",model.stat[ms],",th=th[i])
         }
         Mat.Index$MST[",(sp-1)*n.ms+ms,"] <- max(MST.",species[sp],".",model.stat[ms],".vect) 
         Mat.Index$th.MST[",(sp-1)*n.ms+ms,"] <- th.ghvi(MST.",species[sp],".",model.stat[ms],".vect)",sep="")))
  }
}

#= TSS: True Skill Statistics (MST-1)
Mat.Index$TSS <- Mat.Index$MST-1

#= LPT: Lowest presence threshold
for (sp in 1:length(species)) {
  for (ms in 1:length(model.stat)) {
    eval(parse(text=paste("
         Mat.Index$LPT[",(sp-1)*n.ms+ms,"] <- min(Data.fitted$Fit.",species[sp],".",model.stat[ms],"[Data.fitted$",species[sp],"==1])",sep="")))
  }
}

#= Backup
Mat.Index.rd <- Mat.Index
Mat.Index.rd[,3:7] <- round(Mat.Index.rd[,3:7],digit=3)
write.table(Mat.Index.rd,file="Mat.Index.txt",row.names=FALSE,quote=FALSE,sep="\t")
write.table(Data.fitted,file="Data.fitted.txt",row.names=FALSE,quote=FALSE,sep="\t")

#====================================================================
# 1) Model averaging on probabilities
# 2) Converting probabilities in presence-absence for selected models
# 3) Model averaging on presence/absence
#====================================================================
# we used th.MST as threshold

#= Ag
ncol.Ag <- grep(pattern="^Fit.Ag.*",x=names(Data.fitted))
nrow.Ag <- c(1:3)
fit.Ag <- as.data.frame(Data.fitted[,ncol.Ag])
for (i in 1:length(ncol.Ag)) {
  fit.Ag[,i] <- ifelse(fit.Ag[,i]>=Mat.Index$th.MST[nrow.Ag[i]],1,0)
}
Data.fitted$Proba.Ag <- apply(as.data.frame(Data.fitted[,ncol.Ag]),1,mean)
Data.fitted$SDA.Ag <- apply(fit.Ag,1,mean)
Proba_Ag <- Data.fitted$Proba.Ag
SDA_Ag <- Data.fitted$SDA.Ag
#= Ap
ncol.Ap <- grep(pattern="^Fit.Ap.*",x=names(Data.fitted))
nrow.Ap <- c(1:3)
fit.Ap <- as.data.frame(Data.fitted[,ncol.Ap])
for (i in 1:length(ncol.Ap)) {
  fit.Ap[,i] <- ifelse(fit.Ap[,i]>=Mat.Index$th.MST[nrow.Ap[i]],1,0)
}
Data.fitted$Proba.Ap <- apply(as.data.frame(Data.fitted[,ncol.Ap]),1,mean)
Data.fitted$SDA.Ap <- apply(fit.Ap,1,mean)
Proba_Ap <- Data.fitted$Proba.Ap
SDA_Ap <- Data.fitted$SDA.Ap
#= As
ncol.As <- grep(pattern="^Fit.As.*",x=names(Data.fitted))
nrow.As <- c(1:3)
fit.As <- as.data.frame(Data.fitted[,ncol.As])
for (i in 1:length(ncol.As)) {
  fit.As[,i] <- ifelse(fit.As[,i]>=Mat.Index$th.MST[nrow.As[i]],1,0)
}
Data.fitted$Proba.As <- apply(as.data.frame(Data.fitted[,ncol.As]),1,mean)
Data.fitted$SDA.As <- apply(fit.As,1,mean)
Proba_As <- Data.fitted$Proba.As
SDA_As <- Data.fitted$SDA.As

#====================================================================
# Exporting results to SpatialGridDataFrame and to GRASS raster
#====================================================================

#= Point coordinates
Coords <- Data.fitted[,c(1,2)]
names(Coords)=c("X","Y")
species <- c("Ag","Ap","As")

#= Converting Proba
for (sp in 1:length(species)) {

  eval(parse(text=paste("
  # We transform the simulations in a SpatialPointsDataFrame
  Proba.",species[sp],".data <- cbind(Proba_",species[sp],",Coords)
  Proba.",species[sp],".vect <- SpatialPointsDataFrame(coords=Coords,
                                       data=Proba.",species[sp],".data,
                                       proj4string = CRS(\"+proj=longlat +ellps=WGS84 +datum=WGS84 +n_defs +towgs84=0,0,0\"), # Coordinate Reference System (CRS)
                                       match.ID = TRUE,
                                       bbox = NULL)

  # Exporting SpatialPointsDataFrame object to GRASS as a vector
  writeVECT6(Proba.",species[sp],".vect,
             vname=\"Proba_",species[sp],"_vect\",
             v.in.ogr_flags=\"overwrite\")
  
  # Converting vector to raster for fitted probability of presence
  system(\"v.to.rast input=Proba_",species[sp],"_vect output=Proba_",species[sp],"_rast col=Proba_",species[sp]," --o\")",sep="")))

}

#= Converting SDA
for (sp in 1:length(species)) {

  eval(parse(text=paste("
  # We transform the simulations in a SpatialPointsDataFrame
  SDA.",species[sp],".data <- cbind(SDA_",species[sp],",Coords)
  SDA.",species[sp],".vect <- SpatialPointsDataFrame(coords=Coords,
                                       data=SDA.",species[sp],".data,
                                       proj4string = CRS(\"+proj=longlat +ellps=WGS84 +datum=WGS84 +n_defs +towgs84=0,0,0\"), # Coordinate Reference System (CRS)
                                       match.ID = TRUE,
                                       bbox = NULL)

  # Exporting SpatialPointsDataFrame object to GRASS as a vector
  writeVECT6(SDA.",species[sp],".vect,
             vname=\"SDA_",species[sp],"_vect\",
             v.in.ogr_flags=\"overwrite\")
  
  # Converting vector to raster for fitted probability of presence
  system(\"v.to.rast input=SDA_",species[sp],"_vect output=SDA_",species[sp],"_rast col=SDA_",species[sp]," --o\")",sep="")))

}

#===============================================================================
# Simulations
#===============================================================================

# Loop for predictions regarding
#1: Species (3: Ag, As, Ap)
#2: Statistical model (3: glm, gam, Maxent)
#3: Covariates (1: suffix 1 for glm, gam and Maxent)
#4: Climatic model (3: m1, m2, m3)
#5: Scenario (2: s1, s2)
#6: Year (2: y1, y2)
# TOTAL = 36 predictions by species = 108 simulations  

#= Library
library(ade4)

#= Variables
Model <- c("CCCMA","CSIRO","HADCM3")
model <- c("1","2","3")
Scenario <- c("A2a","B2a")
scenario <- c("1","2")
Year <- c("2050","2080")
year <- c("1","2")
species <- c("Ag","Ap","As")
glmgam <- c("glm","gam")

#= PCA info to project supplementary individuals
load("PCA")

#= Set region
system("g.region res=00:01")

#= Initial data frame with original variables
DataPCA <- read.table(file="DataPCA.txt",header=FALSE,sep="|")
names.species <- c("Ag","Ap","As")
names.bio <- paste("bio",c(1:19),sep="")
names(DataPCA) <- c("Long","Lat","Maskout",names.species,names.bio,"geol")
Data.sim.init <- DataPCA[,c(1:2,4:6,26)]
write.table(Data.sim.init,file="Data.sim.init.txt",sep="|",quote=FALSE,row.names=FALSE,col.names=FALSE) # to have "|" separator

#====================================================
# Function for data process to include future climate
#====================================================

sim.data <- function(m,s,y) {
  cat("Data process...")
  
  #= Make a copy of Data.sim.init
  system("cp Data.sim.init.txt Data.sim.txt")

  #= We fill Data.sim with the new bioclimatic variables
  eval(parse(text=paste("
  Layer <- paste(\"bio\",c(1:19),\"_",Model[m],"_",Scenario[s],"_",Year[y],"\",sep=\"\")",sep=""))) #
  n.Layer <- length(Layer)
  for (i in 1:n.Layer) {
    system(paste("r.out.ascii ",Layer[i]," null=-9999 > Maskout.txt",sep=""))
    system("awk 'NR>=7{for (i=1;i<=(NF);i++) {print $i}}' Maskout.txt > Temp.txt")
    system("awk '{OFS=\"|\"; getline x < \"Temp.txt\"; print $0,x}' Data.sim.txt > Results.txt")
    system("mv Results.txt Data.sim.txt")
  }
  Data.sim <- read.table(file="Data.sim.txt",header=FALSE,sep="|")
  colnames(Data.sim) <- c("Long","Lat","Ag","Ap","As","geol",paste("bio",c(1:19),sep=""))

  #= We compute the coordinates of each pixel (supplementary individuals) on the 6 axis of the original PCA
  SupInd <- suprow(PCA,Data.sim[,c(7:25)])$lisup

  #= New data-set and cleaning
  Data.sim <- cbind(Data.sim[,c(1:6)],SupInd)
  colnames(Data.sim) <- c("Long","Lat","Ag","Ap","As","geol",paste("Axis",c(1:6),sep=""))
  Data.sim[Data.sim==-9999] <- NA
  Vect.Row.notisna <- !is.na(apply(Data.sim,1,mean)) # Identifying row with the presence of NA
  Data.sim <- Data.sim[Vect.Row.notisna,]
  Data.sim$geol <- as.factor(Data.sim$geol)
  system("rm Data.sim.txt")
  cat("done\n")
  return(Data.sim)
}

#==========================================
# Function for predictions with GLM and GAM
#==========================================

sim.glmgam <- function(m,s,y) {
  cat("GLM and GAM predictions...\n")
  # Explicative variables
  X <- Data.sim[,c(6:12)]
  for (sp in 1:length(species)) {
    cat(paste(species[sp]," ",sep=""))
    for (g in 1:length(glmgam)) {
      cat(paste(glmgam[g]," ",sep=""))
      # Computing fitted values
      eval(parse(text=paste("
           Response <- predict.",glmgam[g],"(object=",glmgam[g],".",species[sp],",newdata=X,type=\"response\")
           Data.sim$Sim.",species[sp],".",glmgam[g],".",model[m],"",scenario[s],"",year[y]," <- Response",sep="")))
    }
    cat("\n")
  }
  return(Data.sim)
}

#=====================================
# Function for predictions with Maxent
#=====================================

sim.maxent <- function(m,s,y) {
  cat("Maxent predictions...\n")
  for (sp in 1:length(species)) {
    cat(paste(species[sp]," maxent ",sep=""))
    # projection
    projection <- Data.sim[,c(3,1:2,6:12)]
    projection$Ag <- "background"
    names(projection) <- c("species","longitude","latitude","geol",paste("Axis",c(1:6),sep=""))
    write.table(projection,file="./maxent_data/layers/projection.csv",sep=",",quote=FALSE,row.names=FALSE)
    
    # Call to maxent program
    system(paste("java -cp ./maxent_prog/maxent.jar density.Project ./maxent_data/outputs2/",species[sp],".lambdas ./maxent_data/layers/projection.csv ./maxent_data/sim/",species[sp],"_projection asc -r",sep=""))

    # Importing Maxent predictions
    eval(parse(text=paste("
         Data.sim$Sim.",species[sp],".maxent.",model[m],"",scenario[s],"",year[y]," <- scan(file=\"./maxent_data/sim/",species[sp],"_projection.asc\",skip=6)",sep="")))
  }
  return(Data.sim)
}

#===============================
# Backup function
#===============================

sim.backup <- function(m,s,y) {
  cat("Backup...")
  eval(parse(text=paste("
  write.table(Data.sim[,c(1:2,13:ncol(Data.sim))],file=\"Data.sim_",Model[m],"_",Scenario[s],"_",Year[y],".txt\",row.names=FALSE,quote=FALSE,sep=\"\t\")",sep="")))
  cat("done\n")
}

#===============================
# Simulation RUNS
#===============================

for (m in 1:length(Model)) {
  for (s in 1:length(Scenario)) {
    for (y in 1:length(Year)) {
      cat("#======================================#\n")
      cat("# Model ",Model[m],", Scenario ",Scenario[s],", Year ",Year[y]," #\n",sep="")
      cat("#======================================#\n")
      Data.sim <- sim.data(m,y,s)
      Data.sim <- sim.glmgam(m,y,s)
      Data.sim <- sim.maxent(m,y,s)
      sim.backup(m,y,s)
      cat("\n")
    }
  }
}

#==============================
# Model averaging
#==============================

Data.av <- Data.present[,c(1:2)]
Mat.Index <- read.table(file="Mat.Index.txt",header=TRUE,sep="\t")

for (y in 1:length(Year)) {
  cat(paste("#=== Year: ",Year[y],"\n",sep=""))
  for (s in 1:length(Scenario)) {
    cat(paste("* Scenario: ",Scenario[s],"\n",sep=""))
    for (m in 1:length(Model)) {
      eval(parse(text=paste("
      Data.sim_",Model[m],"_",Scenario[s],"_",Year[y]," <- read.table(file=\"Data.sim_",Model[m],"_",Scenario[s],"_",Year[y],".txt\",header=TRUE,sep=\"\t\")",sep="")))
    }
  
    names.obj <- paste(ls(pattern="^Data.sim_[CH][CSA]*"),sep="",collapse=",")
    eval(parse(text=paste("Data.sim_",Scenario[s],"_",Year[y]," <- cbind(",names.obj,")",sep="")))
    eval(parse(text=paste("rm(",names.obj,")",sep="")))
    cat("** Species: ")
    for (sp in 1:length(species)) {
      cat(species[sp]," ",sep="")
      eval(parse(text=paste("
           Mat.Index.sp <- Mat.Index[Mat.Index$Sp==\"",species[sp],"\",]
           colsp <- grep(pattern=\"Sim.",species[sp],".*\",x=names(Data.sim_",Scenario[s],"_",Year[y],"),perl=FALSE)
           Data.calc <- Data.sim_",Scenario[s],"_",Year[y],"[,c(1:3,colsp)]",sep="")))
      # Mean for probability of presence
      Data.calc$Prob <- apply(Data.calc[,-c(1:3)],1,mean)
      # Cut-off for Species Distribution Area
      th.MST <- rep(Mat.Index.sp$th.MST,3)
      Data.calc.SDA <- Data.calc[-c(ncol(Data.calc))]
      for (i in 1:9) {
        Data.calc.SDA[Data.calc[,3+i]>=th.MST[i],3+i] <- 1
        Data.calc.SDA[Data.calc[,3+i]<th.MST[i],3+i] <- 0
      }
      names(Data.calc.SDA) <- gsub("Sim","SDA",names(Data.calc.SDA))
      Data.calc.SDA$SDA <- apply(Data.calc.SDA[,-c(1:3)],1,mean)
      # Result
      eval(parse(text=paste("
           Data.av$Prob.",species[sp],".",Scenario[s],".",Year[y]," <- Data.calc$Prob
           Data.av$SDA.",species[sp],".",Scenario[s],".",Year[y]," <- Data.calc.SDA$SDA",sep="")))
    }
  cat("\n")
  }
}

#==============================
# Backup
#==============================

write.table(Data.av,file="Data.av.txt",row.names=FALSE,quote=FALSE,sep="\t")

#===============================
# Export of simulations to GRASS
#===============================

Coords <- Data.av[,c(1:2)]
names(Coords)=c("X","Y")

# Prob
for (sp in 1:length(species)) {
  for (y in 1:length(Year)) {
    for (s in 1:length(Scenario)) {

      # We transform the Simulations in a SpatialPointsDataFrame
      eval(parse(text=paste("Prob <- Data.av$Prob.",species[sp],".",Scenario[s],".",Year[y],sep="")))
      Prob.data <- cbind(Prob,Coords)
      names(Prob.data) <- c("Prob","X","Y")
      Prob.vect <- SpatialPointsDataFrame(coords=Coords,
                                          data=Prob.data,
                                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +n_defs +towgs84=0,0,0"), # Coordinate Reference System (CRS)
                                          match.ID = TRUE,
                                          bbox = NULL)

      # Exporting SpatialPointsDataFrame object to GRASS as a vector
      writeVECT6(Prob.vect,
                 vname="Prob_vect",
                 v.in.ogr_flags="overwrite")
  
      # Converting vector to raster for fitted probability of presence
      system(paste("v.to.rast input=Prob_vect output=Prob_",species[sp],"_",Scenario[s],"_",Year[y],"_rast col=Prob --o",sep=""))
    }
  }
}

# SDA
for (sp in 1:length(species)) {
  for (y in 1:length(Year)) {
    for (s in 1:length(Scenario)) {

      # We transform the Simulations in a SpatialPointsDataFrame
      eval(parse(text=paste("SDA <- Data.av$SDA.",species[sp],".",Scenario[s],".",Year[y],sep="")))
      SDA.data <- cbind(SDA,Coords)
      names(SDA.data) <- c("SDA","X","Y")
      SDA.vect <- SpatialPointsDataFrame(coords=Coords,
                                          data=SDA.data,
                                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +n_defs +towgs84=0,0,0"), # Coordinate Reference System (CRS)
                                          match.ID = TRUE,
                                          bbox = NULL)

      # Exporting SpatialPointsDataFrame object to GRASS as a vector
      writeVECT6(SDA.vect,
                 vname="SDA_vect",
                 v.in.ogr_flags="overwrite")
  
      # Converting vector to raster for fitted probability of presence
      system(paste("v.to.rast input=SDA_vect output=SDA_",species[sp],"_",Scenario[s],"_",Year[y],"_rast col=SDA --o",sep=""))
    }
  }
}

#================================
# Under zero dispersal hypothesis
#================================

for (sp in 1:length(species)) {
  for (y in 1:length(Year)) {
    for (s in 1:length(Scenario)) {
      
      ##############
      ## Raster of the Species Distribution Area (SDA)
      system(paste("r.mapcalc 'SDAzD_",species[sp],"_",Scenario[s],"_",Year[y],"_rast=if(SDA_",species[sp],"_rast>0 && SDA_",species[sp],"_",Scenario[s],"_",Year[y],"_rast>0,SDA_",species[sp],"_",Scenario[s],"_",Year[y],"_rast,0)'",sep=""))
      
    }
  }
}

#================================
# Graphics scenario A2a
#================================

#== Replace zeros by null()

# Species distribution area 2010
system("r.mapcalc 'SDA_Ag_rast=if(SDA_Ag_rast==0,null(),SDA_Ag_rast)'")
system("r.mapcalc 'SDA_Ap_rast=if(SDA_Ap_rast==0,null(),SDA_Ap_rast)'")
system("r.mapcalc 'SDA_As_rast=if(SDA_As_rast==0,null(),SDA_As_rast)'")
# Species distribution area A2a_2050_zDisp
system("r.mapcalc 'SDAzD_Ag_A2a_2050_rast=if(SDAzD_Ag_A2a_2050_rast==0,null(),SDAzD_Ag_A2a_2050_rast)'")
system("r.mapcalc 'SDAzD_Ap_A2a_2050_rast=if(SDAzD_Ap_A2a_2050_rast==0,null(),SDAzD_Ap_A2a_2050_rast)'")
system("r.mapcalc 'SDAzD_As_A2a_2050_rast=if(SDAzD_As_A2a_2050_rast==0,null(),SDAzD_As_A2a_2050_rast)'")
# Species distribution area A2a_2080_zDisp
system("r.mapcalc 'SDAzD_Ag_A2a_2080_rast=if(SDAzD_Ag_A2a_2080_rast==0,null(),SDAzD_Ag_A2a_2080_rast)'")
system("r.mapcalc 'SDAzD_Ap_A2a_2080_rast=if(SDAzD_Ap_A2a_2080_rast==0,null(),SDAzD_Ap_A2a_2080_rast)'")
system("r.mapcalc 'SDAzD_As_A2a_2080_rast=if(SDAzD_As_A2a_2080_rast==0,null(),SDAzD_As_A2a_2080_rast)'")

#== Cutting layers

#= As
system("g.region vect=MadaOutline_cut")
system("r.mapcalc As_rast_cut=As_rast")
system("r.mapcalc Proba_As_rast_cut=Proba_As_rast")
system("r.mapcalc Prob_As_A2a_2050_rast_cut=Prob_As_A2a_2050_rast")
system("r.mapcalc Prob_As_A2a_2080_rast_cut=Prob_As_A2a_2080_rast")
system("r.mapcalc SDA_As_rast_cut=SDA_As_rast")
system("r.mapcalc SDAzD_As_A2a_2050_rast_cut=SDAzD_As_A2a_2050_rast")
system("r.mapcalc SDAzD_As_A2a_2080_rast_cut=SDAzD_As_A2a_2080_rast")
#= Ap
system("g.region vect=MadaOutline_cut2")
system("r.mapcalc Ap_rast_cut2=Ap_rast")
system("r.mapcalc Proba_Ap_rast_cut2=Proba_Ap_rast")
system("r.mapcalc Prob_Ap_A2a_2050_rast_cut2=Prob_Ap_A2a_2050_rast")
system("r.mapcalc Prob_Ap_A2a_2080_rast_cut2=Prob_Ap_A2a_2080_rast")
system("r.mapcalc SDA_Ap_rast_cut2=SDA_Ap_rast")
system("r.mapcalc SDAzD_Ap_A2a_2050_rast_cut2=SDAzD_Ap_A2a_2050_rast")
system("r.mapcalc SDAzD_Ap_A2a_2080_rast_cut2=SDAzD_Ap_A2a_2080_rast")
#= Original region
system("g.region vect=Geol")

#== Importing GRASS layers

# Madagascar frontiers
MadaOutline <- readVECT6(vname="MadaOutline",plugin=FALSE,mapset="PERMANENT")
MadaOutline.cut <- readVECT6(vname="MadaOutline_cut",plugin=FALSE,mapset="PERMANENT")
MadaOutline.cut2 <- readVECT6(vname="MadaOutline_cut2",plugin=FALSE,mapset="PERMANENT")
# Protected area
ProArea.vect <- readVECT6(vname="ProArea",plugin=FALSE,mapset="PERMANENT")
ProArea.cut.vect <- readVECT6(vname="ProArea_cut",plugin=FALSE,mapset="PERMANENT")
ProArea.cut2.vect <- readVECT6(vname="ProArea_cut2",plugin=FALSE,mapset="PERMANENT")
# Presence-absence data points
Ag.rast <- readRAST6(vname="Ag_rast",plugin=FALSE,mapset="PERMANENT")
Ap.rast <- readRAST6(vname="Ap_rast",plugin=FALSE,mapset="PERMANENT")
As.rast <- readRAST6(vname="As_rast",plugin=FALSE,mapset="PERMANENT")
# If problem with As.rast colors, use the raster library
library(raster)
As.rast <- raster(As.rast)
#Ag.vect <- readVECT6(vname="Ag",plugin=FALSE,mapset="PERMANENT")
Ap.vect <- readVECT6(vname="Ap",plugin=FALSE,mapset="PERMANENT")
#As.vect <- readVECT6(vname="As",plugin=FALSE,mapset="PERMANENT")

# Fitted probability of presence
Proba.Ag.rast <- readRAST6(vname="Proba_Ag_rast",plugin=FALSE,mapset="ghvi")
Proba.Ap.rast <- readRAST6(vname="Proba_Ap_rast_cut2",plugin=FALSE,mapset="ghvi")
Proba.As.rast <- readRAST6(vname="Proba_As_rast_cut",plugin=FALSE,mapset="ghvi")
# Simulated probability of presence A2a_2050
Prob.Ag.A2a.2050.rast <- readRAST6(vname="Prob_Ag_A2a_2050_rast",plugin=FALSE,mapset="ghvi")
Prob.Ap.A2a.2050.rast <- readRAST6(vname="Prob_Ap_A2a_2050_rast_cut2",plugin=FALSE,mapset="ghvi")
Prob.As.A2a.2050.rast <- readRAST6(vname="Prob_As_A2a_2050_rast_cut",plugin=FALSE,mapset="ghvi")
# Simulated probability of presence A2a_2080
Prob.Ag.A2a.2080.rast <- readRAST6(vname="Prob_Ag_A2a_2080_rast",plugin=FALSE,mapset="ghvi")
Prob.Ap.A2a.2080.rast <- readRAST6(vname="Prob_Ap_A2a_2080_rast_cut2",plugin=FALSE,mapset="ghvi")
Prob.As.A2a.2080.rast <- readRAST6(vname="Prob_As_A2a_2080_rast_cut",plugin=FALSE,mapset="ghvi")

# Species distribution area 2010
SDA.Ag.rast <- readRAST6(vname="SDA_Ag_rast",plugin=FALSE,mapset="ghvi")
SDA.Ap.rast <- readRAST6(vname="SDA_Ap_rast_cut2",plugin=FALSE,mapset="ghvi")
SDA.As.rast <- readRAST6(vname="SDA_As_rast_cut",plugin=FALSE,mapset="ghvi")
# Species distribution area A2a_2050_zDisp
SDAzD.Ag.A2a.2050 <- readRAST6(vname="SDAzD_Ag_A2a_2050_rast",plugin=FALSE,mapset="ghvi")
SDAzD.Ap.A2a.2050 <- readRAST6(vname="SDAzD_Ap_A2a_2050_rast_cut2",plugin=FALSE,mapset="ghvi")
SDAzD.As.A2a.2050 <- readRAST6(vname="SDAzD_As_A2a_2050_rast_cut",plugin=FALSE,mapset="ghvi")
# Species distribution area A2a_2080_zDisp
SDAzD.Ag.A2a.2080 <- readRAST6(vname="SDAzD_Ag_A2a_2080_rast",plugin=FALSE,mapset="ghvi")
SDAzD.Ap.A2a.2080 <- readRAST6(vname="SDAzD_Ap_A2a_2080_rast_cut2",plugin=FALSE,mapset="ghvi")
SDAzD.As.A2a.2080 <- readRAST6(vname="SDAzD_As_A2a_2080_rast_cut",plugin=FALSE,mapset="ghvi")

#================================================#
# Converting rasters with probabilities into lists
#================================================#

## # Fitted probability of presence
## List.Ag <- as.image.SpatialGridDataFrame(Proba.Ag.rast)
## List.Ap <- as.image.SpatialGridDataFrame(Proba.Ap.rast)
## List.As <- as.image.SpatialGridDataFrame(Proba.As.rast)
## # Simulated probability of presence A2a_2050
## List.Ag.A2a.2050 <- as.image.SpatialGridDataFrame(Prob.Ag.A2a.2050.rast)
## List.Ap.A2a.2050 <- as.image.SpatialGridDataFrame(Prob.Ap.A2a.2050.rast)
## List.As.A2a.2050 <- as.image.SpatialGridDataFrame(Prob.As.A2a.2050.rast)

#========================================================
# Maximal probability of presence in fitted or projection
#========================================================

pmax.Ag <- 1 #round(max(List.Ag$z,List.Ag.A2a.2050$z,na.rm=TRUE),1)
pmax.Ap <- 1 #round(max(List.Ap$z,List.Ap.A2a.2050$z,na.rm=TRUE),1)
pmax.As <- 1 #round(max(List.As$z,List.As.A2a.2050$z,na.rm=TRUE),1)

#=======
# Colors
#=======

library(fields)
gcolors.Ag <- colorRampPalette(c("transparent","orange","red","black"))
gcolors.Ap <- colorRampPalette(c("transparent","purple","blue","black"))
gcolors.As <- colorRampPalette(c("transparent","green","dark green","black"))
gcolors.Legend <- colorRampPalette(c("transparent","black"))

#===========================
# Plot 1
#===========================

pdf(file="Plot1.pdf",width=9,height=7)

Mat.plot <- matrix(c(1,2,3,4,0,1,2,3,4,13,5,6,7,8,13,9,10,11,12,0),ncol=5,nrow=4,byrow=TRUE)
layout(Mat.plot,widths=c(rep(680/4,4),80),heights=rep(400,4))
par(mar=c(0,0,0,0),oma=c(0,2,2,0),cex=1.2)

#= Ag
plot(MadaOutline,border=grey(0.5))
image(Ag.rast,col="red",add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
# Zoom
rect(xleft=48.42,ybottom=-13.62,xright=50.05,ytop=-11.89,col="transparent",border="black")
rect(xleft=47.0,ybottom=-16.2,xright=50.6,ytop=-11.89,col="transparent",border="black")
# Tropic of Capricorn
segments(x0=42,y0=-23.5,x1=49,y1=-23.5,lty=2)
text(x=49.2,y=-23.5,label="Tropic of\n Capricorn",cex=0.7)
# Scale
## arrows(x0=50,y0=-22.5,x1=51,y1=-22.5,length=0.05,angle=30,code=2)
## text(x=51.5,y=-22.5,label="1° E",cex=0.8)
## arrows(x0=50,y0=-22.5,x1=50,y1=-21.5,length=0.05,angle=30,code=2)
## text(x=50,y=-21.2,label="1° N",cex=0.8)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Proba.Ag.rast)
image(xyz,col=gcolors.Ag(64),zlim=c(0,pmax.Ag),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ag.A2a.2050.rast)
image(xyz,col=gcolors.Ag(64),zlim=c(0,pmax.Ag),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ag.A2a.2080.rast)
image(xyz,col=gcolors.Ag(64),zlim=c(0,pmax.Ag),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)

#= Ap
plot(MadaOutline.cut2,border=grey(0.5))
points(Ap.vect,col="blue",pch="+",cex=0.8)
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Proba.Ap.rast)
image(xyz,col=gcolors.Ap(64),zlim=c(0,pmax.Ap),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ap.A2a.2050.rast)
image(xyz,col=gcolors.Ap(64),zlim=c(0,pmax.Ap),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ap.A2a.2080.rast)
image(xyz,col=gcolors.Ap(64),zlim=c(0,pmax.Ap),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
       
#= As
plot(MadaOutline.cut,border=grey(0.5))
image(As.rast,col="dark green",add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Proba.As.rast)
image(xyz,col=gcolors.As(64),zlim=c(0,pmax.As),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.As.A2a.2050.rast)
image(xyz,col=gcolors.As(64),zlim=c(0,pmax.As),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.As.A2a.2080.rast)
image(xyz,col=gcolors.As(64),zlim=c(0,pmax.As),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)

#= Legend
par(mar=c(0,1,0,3.5),cex=0.8)
plot.new()
plot.window(xlim=c(0,3),ylim=c(0, 1))
rect(0,seq(0,1,length=65)[-65],
      1,seq(0,1,length=65)[-1],
      col=gcolors.Ag(64),border=NA)
rect(1,seq(0,1,length=65)[-65],
      2,seq(0,1,length=65)[-1],
      col=gcolors.Ap(64),border=NA)
rect(2,seq(0,1,length=65)[-65],
      3,seq(0,1,length=65)[-1],
      col=gcolors.As(64),border=NA)
rect(0,0,3,1)
axis(4,at=c(0,1),labels=c(0,1),las=3,line=0.5)

#=
mtext(text=c(expression(bold("Presence data 2010")),expression(bold("Prob. 2010")),
        expression(bold("Prob. 2050")), expression(bold("Prob. 2080"))),
      outer=TRUE,at=c(0.13,0.35,0.57,0.80),side=3,cex=1.2)
mtext(text=c(expression(italic("A. suarezensis")),expression(italic("A. perrieri")),
        expression(italic("A. grandidieri"))),outer=TRUE,at=c(0.125,0.375,0.75),side=2,line=0.5,cex=1.2)

dev.off()

#= Convert from pdf to png
system("convert -density 300 Plot1.pdf Plot1.png")

#===========================
# Plot 2
#===========================

pdf(file="Plot2.pdf")

Mat.plot <- matrix(c(1,2,3,0,1,2,3,10,4,5,6,10,7,8,9,0),ncol=4,nrow=4,byrow=TRUE)
layout(Mat.plot,widths=c(rep(500/3,3),80),heights=rep(500,4))
par(mar=c(0,0,0,0),oma=c(0,2,1.5,0),cex=1.2)

#= Ag
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDA.Ag.rast)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"red"),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
plot(ProArea.vect,add=TRUE,col="#8080804C")
# Zoom
rect(xleft=48.42,ybottom=-13.62,xright=50.05,ytop=-11.89,col="transparent",border="black")
rect(xleft=47.0,ybottom=-16.2,xright=50.6,ytop=-11.89,col="transparent",border="black")
# Tropic of Capricorn
segments(x0=42,y0=-23.5,x1=49,y1=-23.5,lty=2)
text(x=49.2,y=-23.5,label="Tropic of\n Capricorn",cex=0.7)
# Scale
## arrows(x0=50,y0=-22.5,x1=51,y1=-22.5,length=0.05,angle=30,code=2)
## text(x=51.5,y=-22.5,label="1° E",cex=0.8)
## arrows(x0=50,y0=-22.5,x1=50,y1=-21.5,length=0.05,angle=30,code=2)
## text(x=50,y=-21.2,label="1° N",cex=0.8)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ag.A2a.2050)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"red"),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
plot(ProArea.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ag.A2a.2080)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"red"),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
plot(ProArea.vect,add=TRUE,col="#8080804C")

#= Ap
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDA.Ap.rast)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"blue"),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
plot(ProArea.cut2.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ap.A2a.2050)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"blue"),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
plot(ProArea.cut2.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ap.A2a.2080)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"blue"),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
plot(ProArea.cut2.vect,add=TRUE,col="#8080804C")

#= As
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDA.As.rast)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
plot(ProArea.cut.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.As.A2a.2050)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
plot(ProArea.cut.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.As.A2a.2080)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
plot(ProArea.cut.vect,add=TRUE,col="#8080804C")

#= Legend
par(mar=c(0,1,0,3.5),cex=0.8)
plot.new()
plot.window(xlim=c(0,3),ylim=c(0, 1))
rect(0,seq(0,0.5,length=6),
      1,c(seq(0,0.5,length=6)[-1],1),
      col=c(grey(seq(0.9,0.5,-0.1)),"red"),border="black")
rect(1,seq(0,0.5,length=6),
      2,c(seq(0,0.5,length=6)[-1],1),
      col=c(grey(seq(0.9,0.5,-0.1)),"blue"),border="black")
rect(2,seq(0,0.5,length=6),
      3,c(seq(0,0.5,length=6)[-1],1),
      col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),border="black")
rect(0,0,3,1)
axis(4,at=c(seq(0,0.5,0.1),1),labels=c(seq(0,0.5,0.1),1),las=3,line=0.5)

#=
mtext(text=c(expression(bold("SDA 2010")),expression(bold("SDA 2050")),
        expression(bold("SDA 2080"))),outer=TRUE,at=c(0.16,0.46,0.76),side=3,cex=1.2)
mtext(text=c(expression(italic("A. suarezensis")),expression(italic("A. perrieri")),
        expression(italic("A. grandidieri"))),outer=TRUE,at=c(0.125,0.375,0.75),side=2,line=0,cex=1.2)

dev.off()

#= Convert from pdf to png
system("convert -density 300 Plot2.pdf Plot2.png")

#==============================================================
# Overlap between species distribution area and protected areas
#==============================================================

Mat.Pro.SDA <- as.data.frame(matrix(NA,nrow=9,ncol=6))
names(Mat.Pro.SDA) <- c("Sp","Year","SDA","Pro.SDA","Perc","UZ")
Mat.Pro.SDA$Sp <- rep(c("Ag","Ap","As"),each=3)
Mat.Pro.SDA$Year <- rep(c(2010,2050,2080),3)

species <- c("Ag","Ap","As")

#= SDA
for (sp in 1:length(species)) {
  system(paste("r.mapcalc 'Rast_temp=if(SDA_",species[sp],"_rast>=0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=SDAcalc.txt")
  Temp <- read.table(file="SDAcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2010] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_A2a_2050_rast>=0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=SDAcalc.txt")
  Temp <- read.table(file="SDAcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2050] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  if (sp!=3) {
    system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_A2a_2080_rast>=0.5,1,null())'",sep=""))
    system("r.stats -an input=Rast_temp fs=tab output=SDAcalc.txt")
    Temp <- read.table(file="SDAcalc.txt",header=FALSE,sep="\t")
    eval(parse(text=paste("Mat.Pro.SDA$SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2080] <- round(as.numeric(Temp[2])/1000000)",sep="")))
  }
}

#= Pro.SDA
for (sp in 1:length(species)) {
  system(paste("r.mapcalc 'Rast_temp=if(SDA_",species[sp],"_rast>=0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp,ProArea_rast fs=tab output=ProSDAcalc.txt")
  Temp <- read.table(file="ProSDAcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$Pro.SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2010] <- round(as.numeric(Temp[3])/1000000)",sep=""))) # in km^2

  if (sp!=3) {
    system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_A2a_2050_rast>=0.5,1,null())'",sep=""))
    system("r.stats -an input=Rast_temp,ProArea_rast fs=tab output=ProSDAcalc.txt")
    Temp <- read.table(file="ProSDAcalc.txt",header=FALSE,sep="\t")
    eval(parse(text=paste("Mat.Pro.SDA$Pro.SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2050] <- round(as.numeric(Temp[3])/1000000)",sep=""))) # in km^2

    system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_A2a_2080_rast>=0.5,1,null())'",sep=""))
    system("r.stats -an input=Rast_temp,ProArea_rast fs=tab output=ProSDAcalc.txt")
    Temp <- read.table(file="ProSDAcalc.txt",header=FALSE,sep="\t")
    eval(parse(text=paste("Mat.Pro.SDA$Pro.SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2080] <- round(as.numeric(Temp[3])/1000000)",sep="")))
  }
}

#= Uncertainty.zone
for (sp in 1:length(species)) {
  system(paste("r.mapcalc 'Rast_temp=if(SDA_",species[sp],"_rast<0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=UZcalc.txt")
  Temp <- read.table(file="UZcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$UZ[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2010] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_A2a_2050_rast<0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=UZcalc.txt")
  Temp <- read.table(file="UZcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$UZ[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2050] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_A2a_2080_rast<0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=UZcalc.txt")
  Temp <- read.table(file="UZcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$UZ[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2080] <- round(as.numeric(Temp[2])/1000000)",sep="")))
}

#= Perc
Mat.Pro.SDA$Perc[seq(1,7,3)] <- round(100*Mat.Pro.SDA$Pro.SDA[seq(1,7,3)]/Mat.Pro.SDA$SDA[seq(1,7,3)],1)
Mat.Pro.SDA$Perc[seq(2,8,3)] <- round(100*Mat.Pro.SDA$Pro.SDA[seq(2,8,3)]/Mat.Pro.SDA$SDA[seq(1,7,3)],1)
Mat.Pro.SDA$Perc[seq(3,9,3)] <- round(100*Mat.Pro.SDA$Pro.SDA[seq(3,9,3)]/Mat.Pro.SDA$SDA[seq(1,7,3)],1)

#= Replacing NA with 0
Mat.Pro.SDA[is.na(Mat.Pro.SDA)] <- 0

#= Backup
write.table(Mat.Pro.SDA,file="Mat.Pro.SDA.txt",row.names=FALSE,quote=FALSE,sep="\t")

#====================================================
# A posteriori environmental niche of species in 2010
#====================================================

# Precip (bio12)
#= Mada
system("r.mapcalc 'niche_Precip_Mada=if(!isnull(Geol_rast),bio12@climpres,null())'")
niche.Precip.Mada <- readRAST6(vname="niche_Precip_Mada",plugin=FALSE,mapset="ghvi")
niche.Precip.Mada <- as.image.SpatialGridDataFrame(niche.Precip.Mada)
min.Precip <- min(niche.Precip.Mada$z,na.rm=TRUE)
max.Precip <- max(niche.Precip.Mada$z,na.rm=TRUE)
#= Ag
system("r.mapcalc 'niche_Precip_Ag=if(SDA_Ag_rast>=0.5,bio12@climpres,null())'")
niche.Precip.Ag <- readRAST6(vname="niche_Precip_Ag",plugin=FALSE,mapset="ghvi")
niche.Precip.Ag <- as.image.SpatialGridDataFrame(niche.Precip.Ag)
#= Ap
system("r.mapcalc 'niche_Precip_Ap=if(SDA_Ap_rast>=0.5,bio12@climpres,null())'")
niche.Precip.Ap <- readRAST6(vname="niche_Precip_Ap",plugin=FALSE,mapset="ghvi")
niche.Precip.Ap <- as.image.SpatialGridDataFrame(niche.Precip.Ap)
#= As
system("r.mapcalc 'niche_Precip_As=if(SDA_As_rast>=0.5,bio12@climpres,null())'")
niche.Precip.As <- readRAST6(vname="niche_Precip_As",plugin=FALSE,mapset="ghvi")
niche.Precip.As <- as.image.SpatialGridDataFrame(niche.Precip.As)


# Temp (bio1)
#= Mada
system("r.mapcalc 'niche_Temp_Mada=if(!isnull(Geol_rast),bio1@climpres,null())'")
niche.Temp.Mada <- readRAST6(vname="niche_Temp_Mada",plugin=FALSE,mapset="ghvi")
niche.Temp.Mada <- as.image.SpatialGridDataFrame(niche.Temp.Mada)
min.Temp <- min(niche.Temp.Mada$z,na.rm=TRUE)
max.Temp <- max(niche.Temp.Mada$z,na.rm=TRUE)
#= Ag
system("r.mapcalc 'niche_Temp_Ag=if(SDA_Ag_rast>=0.5,bio1@climpres,null())'")
niche.Temp.Ag <- readRAST6(vname="niche_Temp_Ag",plugin=FALSE,mapset="ghvi")
niche.Temp.Ag <- as.image.SpatialGridDataFrame(niche.Temp.Ag)
#= Ap
system("r.mapcalc 'niche_Temp_Ap=if(SDA_Ap_rast>=0.5,bio1@climpres,null())'")
niche.Temp.Ap <- readRAST6(vname="niche_Temp_Ap",plugin=FALSE,mapset="ghvi")
niche.Temp.Ap <- as.image.SpatialGridDataFrame(niche.Temp.Ap)
#= As
system("r.mapcalc 'niche_Temp_As=if(SDA_As_rast>=0.5,bio1@climpres,null())'")
niche.Temp.As <- readRAST6(vname="niche_Temp_As",plugin=FALSE,mapset="ghvi")
niche.Temp.As <- as.image.SpatialGridDataFrame(niche.Temp.As)


# PreSeas (bio15)
#= Mada
system("r.mapcalc 'niche_PreSeas_Mada=if(!isnull(Geol_rast),bio15@climpres,null())'")
niche.PreSeas.Mada <- readRAST6(vname="niche_PreSeas_Mada",plugin=FALSE,mapset="ghvi")
niche.PreSeas.Mada <- as.image.SpatialGridDataFrame(niche.PreSeas.Mada)
min.PreSeas <- min(niche.PreSeas.Mada$z,na.rm=TRUE)
max.PreSeas <- max(niche.PreSeas.Mada$z,na.rm=TRUE)
#= Ag
system("r.mapcalc 'niche_PreSeas_Ag=if(SDA_Ag_rast>=0.5,bio15@climpres,null())'")
niche.PreSeas.Ag <- readRAST6(vname="niche_PreSeas_Ag",plugin=FALSE,mapset="ghvi")
niche.PreSeas.Ag <- as.image.SpatialGridDataFrame(niche.PreSeas.Ag)
#= Ap
system("r.mapcalc 'niche_PreSeas_Ap=if(SDA_Ap_rast>=0.5,bio15@climpres,null())'")
niche.PreSeas.Ap <- readRAST6(vname="niche_PreSeas_Ap",plugin=FALSE,mapset="ghvi")
niche.PreSeas.Ap <- as.image.SpatialGridDataFrame(niche.PreSeas.Ap)
#= As
system("r.mapcalc 'niche_PreSeas_As=if(SDA_As_rast>=0.5,bio15@climpres,null())'")
niche.PreSeas.As <- readRAST6(vname="niche_PreSeas_As",plugin=FALSE,mapset="ghvi")
niche.PreSeas.As <- as.image.SpatialGridDataFrame(niche.PreSeas.As)


# TempSeas (bio4)
#= Mada
system("r.mapcalc 'niche_TempSeas_Mada=if(!isnull(Geol_rast),bio4@climpres,null())'")
niche.TempSeas.Mada <- readRAST6(vname="niche_TempSeas_Mada",plugin=FALSE,mapset="ghvi")
niche.TempSeas.Mada <- as.image.SpatialGridDataFrame(niche.TempSeas.Mada)
min.TempSeas <- min(niche.TempSeas.Mada$z,na.rm=TRUE)
max.TempSeas <- max(niche.TempSeas.Mada$z,na.rm=TRUE)
#= Ag
system("r.mapcalc 'niche_TempSeas_Ag=if(SDA_Ag_rast>=0.5,bio4@climpres,null())'")
niche.TempSeas.Ag <- readRAST6(vname="niche_TempSeas_Ag",plugin=FALSE,mapset="ghvi")
niche.TempSeas.Ag <- as.image.SpatialGridDataFrame(niche.TempSeas.Ag)
#= Ap
system("r.mapcalc 'niche_TempSeas_Ap=if(SDA_Ap_rast>=0.5,bio4@climpres,null())'")
niche.TempSeas.Ap <- readRAST6(vname="niche_TempSeas_Ap",plugin=FALSE,mapset="ghvi")
niche.TempSeas.Ap <- as.image.SpatialGridDataFrame(niche.TempSeas.Ap)
#= As
system("r.mapcalc 'niche_TempSeas_As=if(SDA_As_rast>=0.5,bio4@climpres,null())'")
niche.TempSeas.As <- readRAST6(vname="niche_TempSeas_As",plugin=FALSE,mapset="ghvi")
niche.TempSeas.As <- as.image.SpatialGridDataFrame(niche.TempSeas.As)

#=====
# Ranges
#=====
sink("RangeNiche.txt")
cat(paste("Ag, precip: ",range(niche.Precip.Ag$z,na.rm=TRUE),"\n",sep=""))
cat(paste("Ag, temp: ",range(niche.Temp.Ag$z,na.rm=TRUE),"\n",sep=""))
cat(paste("Ag, preseas: ",range(niche.PreSeas.Ag$z,na.rm=TRUE),"\n",sep=""))
cat(paste("Ag, tempseas: ",range(niche.TempSeas.Ag$z,na.rm=TRUE),"\n",sep=""))
cat("===========\n")
cat(paste("Ap, precip: ",range(niche.Precip.Ap$z,na.rm=TRUE),"\n",sep=""))
cat(paste("Ap, temp: ",range(niche.Temp.Ap$z,na.rm=TRUE),"\n",sep=""))
cat(paste("Ap, preseas: ",range(niche.PreSeas.Ap$z,na.rm=TRUE),"\n",sep=""))
cat(paste("Ap, tempseas: ",range(niche.TempSeas.Ap$z,na.rm=TRUE),"\n",sep=""))
cat("===========\n")
cat(paste("As, precip: ",range(niche.Precip.As$z,na.rm=TRUE),"\n",sep=""))
cat(paste("As, temp: ",range(niche.Temp.As$z,na.rm=TRUE),"\n",sep=""))
cat(paste("As, preseas: ",range(niche.PreSeas.As$z,na.rm=TRUE),"\n",sep=""))
cat(paste("As, tempseas: ",range(niche.TempSeas.As$z,na.rm=TRUE),"\n",sep=""))
sink()

#=====
# Plot
#=====

pdf(file="Plot3.pdf")

Mat.plot <- matrix(c(5,5,1,2,3,4),ncol=2,nrow=3,byrow=TRUE)
layout(Mat.plot,widths=rep(400,2),heights=c(50,rep(400,2)))
par(mar=c(5,4,0,1),oma=c(0,0,0,0),cex=1.2)

#= Precip
hist(niche.Precip.Mada$z,freq=FALSE,
     xlim=c(0,3500),
     ylim=c(0,0.005),
     xlab="Annual precipitations\n(mm.yr-1)",
     main="",
     axes=FALSE)
axis(side=1,at=seq(0,3500,500),labels=seq(0,3500,500))
axis(side=2,at=seq(0,0.005,0.001),labels=seq(0,0.005,0.001))
hist(niche.Precip.Ag$z,add=TRUE,freq=FALSE,border="red")
hist(niche.Precip.Ap$z,add=TRUE,freq=FALSE,border="blue")
hist(niche.Precip.As$z,add=TRUE,freq=FALSE,border="dark green")

#= PreSeas
hist(niche.PreSeas.Mada$z,freq=FALSE,
     xlim=c(30,150),
     ylim=c(0,0.05),
     xlab="Precipitation seasonality\n(mean/sd)",
     main="",
     axes=FALSE)
axis(side=1,at=seq(30,150,20),labels=seq(30,150,20))
axis(side=2,at=seq(0,0.05,0.01),labels=seq(0,0.05,0.01))
hist(niche.PreSeas.Ag$z,add=TRUE,freq=FALSE,border="red")
hist(niche.PreSeas.Ap$z,add=TRUE,freq=FALSE,border="blue")
hist(niche.PreSeas.As$z,add=TRUE,freq=FALSE,border="dark green")

#= Temp
hist(niche.Temp.Mada$z,freq=FALSE,
     xlim=c(100,300),
     ylim=c(0,0.05),
     xlab="Annual mean temperature\n(°C)",
     main="",
     axes=FALSE)
axis(side=1,at=seq(100,300,50),labels=seq(10,30,5))
axis(side=2,at=seq(0,0.05,0.01),labels=seq(0,0.05,0.01))
hist(niche.Temp.Ag$z,add=TRUE,freq=FALSE,border="red")
hist(niche.Temp.Ap$z,add=TRUE,freq=FALSE,border="blue")
hist(niche.Temp.As$z,add=TRUE,freq=FALSE,border="dark green")

#= TempSeas
hist(niche.TempSeas.Mada$z,freq=FALSE,
     xlim=c(700,3300),
     ylim=c(0,0.005),
     xlab="Temperature seasonality\n(sd in °C)",
     main="",
     axes=FALSE)
axis(side=1,at=seq(700,3300,400),labels=seq(0.7,3.3,0.4))
axis(side=2,at=seq(0,0.005,0.0025),labels=seq(0,0.005,0.0025))
hist(niche.TempSeas.Ag$z,add=TRUE,freq=FALSE,border="red")
hist(niche.TempSeas.Ap$z,add=TRUE,freq=FALSE,border="blue")
hist(niche.TempSeas.As$z,add=TRUE,freq=FALSE,border="dark green")

#= Legend
par(mar=c(0.5,0,0,0))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
legend(x=0.5,y=0.5,horiz=TRUE,
       legend=c(expression(italic("A. grandidieri")),expression(italic("A. perrieri")),expression(italic("A. suarezensis"))),
lty=rep(1,3),lwd=c(2,2,2),col=c("red","blue","dark green"),bty="n",cex=1,xjust=0.5,yjust=0.5)

dev.off()

#= Convert from pdf to png
system("convert -density 300 Plot3.pdf Plot3.png")

#======================================================================================
#======================================================================================


#================================
# Graphics scenario B2a
#================================

#== Replace zeros by null()

# Species distribution area 2010
system("r.mapcalc 'SDA_Ag_rast=if(SDA_Ag_rast==0,null(),SDA_Ag_rast)'")
system("r.mapcalc 'SDA_Ap_rast=if(SDA_Ap_rast==0,null(),SDA_Ap_rast)'")
system("r.mapcalc 'SDA_As_rast=if(SDA_As_rast==0,null(),SDA_As_rast)'")
# Species distribution area B2a_2050_zDisp
system("r.mapcalc 'SDAzD_Ag_B2a_2050_rast=if(SDAzD_Ag_B2a_2050_rast==0,null(),SDAzD_Ag_B2a_2050_rast)'")
system("r.mapcalc 'SDAzD_Ap_B2a_2050_rast=if(SDAzD_Ap_B2a_2050_rast==0,null(),SDAzD_Ap_B2a_2050_rast)'")
system("r.mapcalc 'SDAzD_As_B2a_2050_rast=if(SDAzD_As_B2a_2050_rast==0,null(),SDAzD_As_B2a_2050_rast)'")
# Species distribution area B2a_2080_zDisp
system("r.mapcalc 'SDAzD_Ag_B2a_2080_rast=if(SDAzD_Ag_B2a_2080_rast==0,null(),SDAzD_Ag_B2a_2080_rast)'")
system("r.mapcalc 'SDAzD_Ap_B2a_2080_rast=if(SDAzD_Ap_B2a_2080_rast==0,null(),SDAzD_Ap_B2a_2080_rast)'")
system("r.mapcalc 'SDAzD_As_B2a_2080_rast=if(SDAzD_As_B2a_2080_rast==0,null(),SDAzD_As_B2a_2080_rast)'")

#== Cutting layers

#= As
system("g.region vect=MadaOutline_cut")
system("r.mapcalc As_rast_cut=As_rast")
system("r.mapcalc Proba_As_rast_cut=Proba_As_rast")
system("r.mapcalc Prob_As_B2a_2050_rast_cut=Prob_As_B2a_2050_rast")
system("r.mapcalc Prob_As_B2a_2080_rast_cut=Prob_As_B2a_2080_rast")
system("r.mapcalc SDA_As_rast_cut=SDA_As_rast")
system("r.mapcalc SDAzD_As_B2a_2050_rast_cut=SDAzD_As_B2a_2050_rast")
system("r.mapcalc SDAzD_As_B2a_2080_rast_cut=SDAzD_As_B2a_2080_rast")
#= Ap
system("g.region vect=MadaOutline_cut2")
system("r.mapcalc Ap_rast_cut2=Ap_rast")
system("r.mapcalc Proba_Ap_rast_cut2=Proba_Ap_rast")
system("r.mapcalc Prob_Ap_B2a_2050_rast_cut2=Prob_Ap_B2a_2050_rast")
system("r.mapcalc Prob_Ap_B2a_2080_rast_cut2=Prob_Ap_B2a_2080_rast")
system("r.mapcalc SDA_Ap_rast_cut2=SDA_Ap_rast")
system("r.mapcalc SDAzD_Ap_B2a_2050_rast_cut2=SDAzD_Ap_B2a_2050_rast")
system("r.mapcalc SDAzD_Ap_B2a_2080_rast_cut2=SDAzD_Ap_B2a_2080_rast")
#= Original region
system("g.region vect=Geol")

#== Importing GRASS layers

# Madagascar frontiers
MadaOutline <- readVECT6(vname="MadaOutline",plugin=FALSE,mapset="PERMANENT")
MadaOutline.cut <- readVECT6(vname="MadaOutline_cut",plugin=FALSE,mapset="PERMANENT")
MadaOutline.cut2 <- readVECT6(vname="MadaOutline_cut2",plugin=FALSE,mapset="PERMANENT")
# Protected area
ProArea.vect <- readVECT6(vname="ProArea",plugin=FALSE,mapset="PERMANENT")
ProArea.cut.vect <- readVECT6(vname="ProArea_cut",plugin=FALSE,mapset="PERMANENT")
ProArea.cut2.vect <- readVECT6(vname="ProArea_cut2",plugin=FALSE,mapset="PERMANENT")
# Presence-absence data points
Ag.rast <- readRAST6(vname="Ag_rast",plugin=FALSE,mapset="PERMANENT")
Ap.rast <- readRAST6(vname="Ap_rast",plugin=FALSE,mapset="PERMANENT")
As.rast <- readRAST6(vname="As_rast",plugin=FALSE,mapset="PERMANENT")
#Ag.vect <- readVECT6(vname="Ag",plugin=FALSE,mapset="PERMANENT")
Ap.vect <- readVECT6(vname="Ap",plugin=FALSE,mapset="PERMANENT")
#As.vect <- readVECT6(vname="As",plugin=FALSE,mapset="PERMANENT")

# Fitted probability of presence
Proba.Ag.rast <- readRAST6(vname="Proba_Ag_rast",plugin=FALSE,mapset="ghvi")
Proba.Ap.rast <- readRAST6(vname="Proba_Ap_rast_cut2",plugin=FALSE,mapset="ghvi")
Proba.As.rast <- readRAST6(vname="Proba_As_rast_cut",plugin=FALSE,mapset="ghvi")
# Simulated probability of presence B2a_2050
Prob.Ag.B2a.2050.rast <- readRAST6(vname="Prob_Ag_B2a_2050_rast",plugin=FALSE,mapset="ghvi")
Prob.Ap.B2a.2050.rast <- readRAST6(vname="Prob_Ap_B2a_2050_rast_cut2",plugin=FALSE,mapset="ghvi")
Prob.As.B2a.2050.rast <- readRAST6(vname="Prob_As_B2a_2050_rast_cut",plugin=FALSE,mapset="ghvi")
# Simulated probability of presence B2a_2080
Prob.Ag.B2a.2080.rast <- readRAST6(vname="Prob_Ag_B2a_2080_rast",plugin=FALSE,mapset="ghvi")
Prob.Ap.B2a.2080.rast <- readRAST6(vname="Prob_Ap_B2a_2080_rast_cut2",plugin=FALSE,mapset="ghvi")
Prob.As.B2a.2080.rast <- readRAST6(vname="Prob_As_B2a_2080_rast_cut",plugin=FALSE,mapset="ghvi")

# Species distribution area 2010
SDA.Ag.rast <- readRAST6(vname="SDA_Ag_rast",plugin=FALSE,mapset="ghvi")
SDA.Ap.rast <- readRAST6(vname="SDA_Ap_rast_cut2",plugin=FALSE,mapset="ghvi")
SDA.As.rast <- readRAST6(vname="SDA_As_rast_cut",plugin=FALSE,mapset="ghvi")
# Species distribution area B2a_2050_zDisp
SDAzD.Ag.B2a.2050 <- readRAST6(vname="SDAzD_Ag_B2a_2050_rast",plugin=FALSE,mapset="ghvi")
SDAzD.Ap.B2a.2050 <- readRAST6(vname="SDAzD_Ap_B2a_2050_rast_cut2",plugin=FALSE,mapset="ghvi")
SDAzD.As.B2a.2050 <- readRAST6(vname="SDAzD_As_B2a_2050_rast_cut",plugin=FALSE,mapset="ghvi")
# Species distribution area B2a_2080_zDisp
SDAzD.Ag.B2a.2080 <- readRAST6(vname="SDAzD_Ag_B2a_2080_rast",plugin=FALSE,mapset="ghvi")
SDAzD.Ap.B2a.2080 <- readRAST6(vname="SDAzD_Ap_B2a_2080_rast_cut2",plugin=FALSE,mapset="ghvi")
SDAzD.As.B2a.2080 <- readRAST6(vname="SDAzD_As_B2a_2080_rast_cut",plugin=FALSE,mapset="ghvi")

#================================================#
# Converting rasters with probabilities into lists
#================================================#

## # Fitted probability of presence
## List.Ag <- as.image.SpatialGridDataFrame(Proba.Ag.rast)
## List.Ap <- as.image.SpatialGridDataFrame(Proba.Ap.rast)
## List.As <- as.image.SpatialGridDataFrame(Proba.As.rast)
## # Simulated probability of presence B2a_2050
## List.Ag.B2a.2050 <- as.image.SpatialGridDataFrame(Prob.Ag.B2a.2050.rast)
## List.Ap.B2a.2050 <- as.image.SpatialGridDataFrame(Prob.Ap.B2a.2050.rast)
## List.As.B2a.2050 <- as.image.SpatialGridDataFrame(Prob.As.B2a.2050.rast)

#========================================================
# Maximal probability of presence in fitted or projection
#========================================================

pmax.Ag <- 1 #round(max(List.Ag$z,List.Ag.B2a.2050$z,na.rm=TRUE),1)
pmax.Ap <- 1 #round(max(List.Ap$z,List.Ap.B2a.2050$z,na.rm=TRUE),1)
pmax.As <- 1 #round(max(List.As$z,List.As.B2a.2050$z,na.rm=TRUE),1)

#=======
# Colors
#=======

library(fields)
gcolors.Ag <- colorRampPalette(c("transparent","orange","red","black"))
gcolors.Ap <- colorRampPalette(c("transparent","purple","blue","black"))
gcolors.As <- colorRampPalette(c("transparent","green","dark green","black"))
gcolors.Legend <- colorRampPalette(c("transparent","black"))

#===========================
# Plot 1
#===========================

pdf(file="Plot1-B2a.pdf",width=9,height=7)

Mat.plot <- matrix(c(1,2,3,4,0,1,2,3,4,13,5,6,7,8,13,9,10,11,12,0),ncol=5,nrow=4,byrow=TRUE)
layout(Mat.plot,widths=c(rep(680/4,4),80),heights=rep(400,4))
par(mar=c(0,0,0,0),oma=c(0,2,2,0),cex=1.2)

#= Ag
plot(MadaOutline,border=grey(0.5))
image(Ag.rast,col="red",add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
# Zoom
rect(xleft=48.42,ybottom=-13.62,xright=50.05,ytop=-11.89,col="transparent",border="black")
rect(xleft=47.0,ybottom=-16.2,xright=50.6,ytop=-11.89,col="transparent",border="black")
# Tropic of Capricorn
segments(x0=42,y0=-23.5,x1=49,y1=-23.5,lty=2)
text(x=49.2,y=-23.5,label="Tropic of\n Capricorn",cex=0.7)
# Scale
## arrows(x0=50,y0=-22.5,x1=51,y1=-22.5,length=0.05,angle=30,code=2)
## text(x=51.5,y=-22.5,label="1° E",cex=0.8)
## arrows(x0=50,y0=-22.5,x1=50,y1=-21.5,length=0.05,angle=30,code=2)
## text(x=50,y=-21.2,label="1° N",cex=0.8)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Proba.Ag.rast)
image(xyz,col=gcolors.Ag(64),zlim=c(0,pmax.Ag),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ag.B2a.2050.rast)
image(xyz,col=gcolors.Ag(64),zlim=c(0,pmax.Ag),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ag.B2a.2080.rast)
image(xyz,col=gcolors.Ag(64),zlim=c(0,pmax.Ag),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)

#= Ap
plot(MadaOutline.cut2,border=grey(0.5))
points(Ap.vect,col="blue",pch="+",cex=0.8)
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Proba.Ap.rast)
image(xyz,col=gcolors.Ap(64),zlim=c(0,pmax.Ap),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ap.B2a.2050.rast)
image(xyz,col=gcolors.Ap(64),zlim=c(0,pmax.Ap),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.Ap.B2a.2080.rast)
image(xyz,col=gcolors.Ap(64),zlim=c(0,pmax.Ap),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
       
#= As
plot(MadaOutline.cut,border=grey(0.5))
image(As.rast,col="dark green",add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Proba.As.rast)
image(xyz,col=gcolors.As(64),zlim=c(0,pmax.As),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.As.B2a.2050.rast)
image(xyz,col=gcolors.As(64),zlim=c(0,pmax.As),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(Prob.As.B2a.2080.rast)
image(xyz,col=gcolors.As(64),zlim=c(0,pmax.As),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)

#= Legend
par(mar=c(0,1,0,3.5),cex=0.8)
plot.new()
plot.window(xlim=c(0,3),ylim=c(0, 1))
rect(0,seq(0,1,length=65)[-65],
      1,seq(0,1,length=65)[-1],
      col=gcolors.Ag(64),border=NA)
rect(1,seq(0,1,length=65)[-65],
      2,seq(0,1,length=65)[-1],
      col=gcolors.Ap(64),border=NA)
rect(2,seq(0,1,length=65)[-65],
      3,seq(0,1,length=65)[-1],
      col=gcolors.As(64),border=NA)
rect(0,0,3,1)
axis(4,at=c(0,1),labels=c(0,1),las=3,line=0.5)

#=
mtext(text=c(expression(bold("Presence data 2010")),expression(bold("PSH 2010")),
        expression(bold("PSH 2050")), expression(bold("PSH 2080"))),
      outer=TRUE,at=c(0.13,0.35,0.57,0.80),side=3,cex=1.2)
mtext(text=c(expression(italic("A. suarezensis")),expression(italic("A. perrieri")),
        expression(italic("A. grandidieri"))),outer=TRUE,at=c(0.125,0.375,0.75),side=2,line=0.5,cex=1.2)

dev.off()

#= Convert from pdf to png
system("convert -density 300 Plot1-B2a.pdf Plot1-B2a.png")

#===========================
# Plot 2
#===========================

pdf(file="Plot2-B2a.pdf")

Mat.plot <- matrix(c(1,2,3,0,1,2,3,10,4,5,6,10,7,8,9,0),ncol=4,nrow=4,byrow=TRUE)
layout(Mat.plot,widths=c(rep(500/3,3),80),heights=rep(500,4))
par(mar=c(0,0,0,0),oma=c(0,2,1.5,0),cex=1.2)

#= Ag
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDA.Ag.rast)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"red"),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
plot(ProArea.vect,add=TRUE,col="#8080804C")
# Zoom
rect(xleft=48.42,ybottom=-13.62,xright=50.05,ytop=-11.89,col="transparent",border="black")
rect(xleft=47.0,ybottom=-16.2,xright=50.6,ytop=-11.89,col="transparent",border="black")
# Tropic of Capricorn
segments(x0=42,y0=-23.5,x1=49,y1=-23.5,lty=2)
text(x=49.2,y=-23.5,label="Tropic of\n Capricorn",cex=0.7)
# Scale
## arrows(x0=50,y0=-22.5,x1=51,y1=-22.5,length=0.05,angle=30,code=2)
## text(x=51.5,y=-22.5,label="1° E",cex=0.8)
## arrows(x0=50,y0=-22.5,x1=50,y1=-21.5,length=0.05,angle=30,code=2)
## text(x=50,y=-21.2,label="1° N",cex=0.8)
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ag.B2a.2050)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"red"),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
plot(ProArea.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ag.B2a.2080)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"red"),add=TRUE)
plot(MadaOutline,border=grey(0.5),add=TRUE)
plot(ProArea.vect,add=TRUE,col="#8080804C")

#= Ap
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDA.Ap.rast)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"blue"),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
plot(ProArea.cut2.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ap.B2a.2050)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"blue"),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
plot(ProArea.cut2.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut2,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.Ap.B2a.2080)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"blue"),add=TRUE)
plot(MadaOutline.cut2,border=grey(0.5),add=TRUE)
plot(ProArea.cut2.vect,add=TRUE,col="#8080804C")

#= As
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDA.As.rast)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
plot(ProArea.cut.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.As.B2a.2050)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
plot(ProArea.cut.vect,add=TRUE,col="#8080804C")
#=
plot(MadaOutline.cut,border=grey(0.5))
xyz <- as.image.SpatialGridDataFrame(SDAzD.As.B2a.2080)
image(xyz,breaks=c(seq(0,0.5,0.1),1),col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),add=TRUE)
plot(MadaOutline.cut,border=grey(0.5),add=TRUE)
plot(ProArea.cut.vect,add=TRUE,col="#8080804C")

#= Legend
par(mar=c(0,1,0,3.5),cex=0.8)
plot.new()
plot.window(xlim=c(0,3),ylim=c(0, 1))
rect(0,seq(0,0.5,length=6),
      1,c(seq(0,0.5,length=6)[-1],1),
      col=c(grey(seq(0.9,0.5,-0.1)),"red"),border="black")
rect(1,seq(0,0.5,length=6),
      2,c(seq(0,0.5,length=6)[-1],1),
      col=c(grey(seq(0.9,0.5,-0.1)),"blue"),border="black")
rect(2,seq(0,0.5,length=6),
      3,c(seq(0,0.5,length=6)[-1],1),
      col=c(grey(seq(0.9,0.5,-0.1)),"dark green"),border="black")
rect(0,0,3,1)
axis(4,at=c(seq(0,0.5,0.1),1),labels=c(seq(0,0.5,0.1),1),las=3,line=0.5)

#=
mtext(text=c(expression(bold("SDA 2010")),expression(bold("SDA 2050")),
        expression(bold("SDA 2080"))),outer=TRUE,at=c(0.16,0.46,0.76),side=3,cex=1.2)
mtext(text=c(expression(italic("A. suarezensis")),expression(italic("A. perrieri")),
        expression(italic("A. grandidieri"))),outer=TRUE,at=c(0.125,0.375,0.75),side=2,line=0,cex=1.2)

dev.off()

#= Convert from pdf to png
system("convert -density 300 Plot2-B2a.pdf Plot2-B2a.png")

#==============================================================
# Overlap between species distribution area and protected areas
#==============================================================

Mat.Pro.SDA <- as.data.frame(matrix(NA,nrow=9,ncol=6))
names(Mat.Pro.SDA) <- c("Sp","Year","SDA","Pro.SDA","Perc","UZ")
Mat.Pro.SDA$Sp <- rep(c("Ag","Ap","As"),each=3)
Mat.Pro.SDA$Year <- rep(c(2010,2050,2080),3)

species <- c("Ag","Ap","As")

#= SDA
for (sp in 1:length(species)) {
  system(paste("r.mapcalc 'Rast_temp=if(SDA_",species[sp],"_rast>=0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=SDAcalc.txt")
  Temp <- read.table(file="SDAcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2010] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_B2a_2050_rast>=0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=SDAcalc.txt")
  Temp <- read.table(file="SDAcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2050] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  if (sp!=3) {
    system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_B2a_2080_rast>=0.5,1,null())'",sep=""))
    system("r.stats -an input=Rast_temp fs=tab output=SDAcalc.txt")
    Temp <- read.table(file="SDAcalc.txt",header=FALSE,sep="\t")
    eval(parse(text=paste("Mat.Pro.SDA$SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2080] <- round(as.numeric(Temp[2])/1000000)",sep="")))
  }
}

#= Pro.SDA
for (sp in 1:length(species)) {
  system(paste("r.mapcalc 'Rast_temp=if(SDA_",species[sp],"_rast>=0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp,ProArea_rast fs=tab output=ProSDAcalc.txt")
  Temp <- read.table(file="ProSDAcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$Pro.SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2010] <- round(as.numeric(Temp[3])/1000000)",sep=""))) # in km^2

  if (sp!=3) {
    system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_B2a_2050_rast>=0.5,1,null())'",sep=""))
    system("r.stats -an input=Rast_temp,ProArea_rast fs=tab output=ProSDAcalc.txt")
    Temp <- read.table(file="ProSDAcalc.txt",header=FALSE,sep="\t")
    eval(parse(text=paste("Mat.Pro.SDA$Pro.SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2050] <- round(as.numeric(Temp[3])/1000000)",sep=""))) # in km^2

    system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_B2a_2080_rast>=0.5,1,null())'",sep=""))
    system("r.stats -an input=Rast_temp,ProArea_rast fs=tab output=ProSDAcalc.txt")
    Temp <- read.table(file="ProSDAcalc.txt",header=FALSE,sep="\t")
    eval(parse(text=paste("Mat.Pro.SDA$Pro.SDA[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2080] <- round(as.numeric(Temp[3])/1000000)",sep="")))
  }
}

#= Uncertainty.zone
for (sp in 1:length(species)) {
  system(paste("r.mapcalc 'Rast_temp=if(SDA_",species[sp],"_rast<0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=UZcalc.txt")
  Temp <- read.table(file="UZcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$UZ[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2010] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_B2a_2050_rast<0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=UZcalc.txt")
  Temp <- read.table(file="UZcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$UZ[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2050] <- round(as.numeric(Temp[2])/1000000)",sep=""))) # in km^2

  system(paste("r.mapcalc 'Rast_temp=if(SDAzD_",species[sp],"_B2a_2080_rast<0.5,1,null())'",sep=""))
  system("r.stats -an input=Rast_temp fs=tab output=UZcalc.txt")
  Temp <- read.table(file="UZcalc.txt",header=FALSE,sep="\t")
  eval(parse(text=paste("Mat.Pro.SDA$UZ[Mat.Pro.SDA$Sp==\"",species[sp],"\"&Mat.Pro.SDA$Year==2080] <- round(as.numeric(Temp[2])/1000000)",sep="")))
}

#= Perc
Mat.Pro.SDA$Perc[seq(1,7,3)] <- round(100*Mat.Pro.SDA$Pro.SDA[seq(1,7,3)]/Mat.Pro.SDA$SDA[seq(1,7,3)],1)
Mat.Pro.SDA$Perc[seq(2,8,3)] <- round(100*Mat.Pro.SDA$Pro.SDA[seq(2,8,3)]/Mat.Pro.SDA$SDA[seq(1,7,3)],1)
Mat.Pro.SDA$Perc[seq(3,9,3)] <- round(100*Mat.Pro.SDA$Pro.SDA[seq(3,9,3)]/Mat.Pro.SDA$SDA[seq(1,7,3)],1)

#= Replacing NA with 0
Mat.Pro.SDA[is.na(Mat.Pro.SDA)] <- 0

#= Backup
write.table(Mat.Pro.SDA,file="Mat.Pro.SDA.B2a.txt",row.names=FALSE,quote=FALSE,sep="\t")

#=============================================================
# Connectivity of the protected area network from graph theory
#=============================================================

library(maptools) # to convert sp objects to spatstat object
library(spatstat) # to convert sp objects to spatstat object
library(spatgraphs) # to draw graphs for spatial points
#library(igraph)

#== Defining colors
library(fields)
gcolors.Red <- colorRampPalette(c("transparent","red"))
gcolors.Blue <- colorRampPalette(c("transparent","lightblue","blue"))
gcolors.Green <- colorRampPalette(c("transparent","dark green"))

#== Importing GRASS layers
MadaOutline <- readVECT6(vname="MadaOutline",plugin=FALSE,mapset="PERMANENT")
MadaOutline.cut <- readVECT6(vname="MadaOutline_cut",plugin=FALSE,mapset="PERMANENT")
MadaOutline.cut2 <- readVECT6(vname="MadaOutline_cut2",plugin=FALSE,mapset="PERMANENT")

# loop index
type <- c("SDA_Ag_rast","SDAzD_Ag_A2a_2050_rast","SDAzD_Ag_A2a_2080_rast",
          "SDA_Ap_rast_cut2","SDAzD_Ap_A2a_2050_rast_cut2","SDAzD_Ap_A2a_2080_rast_cut2",
          "SDA_As_rast_cut","SDAzD_As_A2a_2050_rast_cut","SDAzD_As_A2a_2080_rast_cut")

Connectivity <- rep(0,9)

pdf(file="Plot-Connectivity.pdf")

Mat.plot <- matrix(c(1,2,3,0,1,2,3,10,4,5,6,10,7,8,9,0),ncol=4,nrow=4,byrow=TRUE)
layout(Mat.plot,widths=c(rep(500/3,3),80),heights=rep(500,4))
par(mar=c(0,0,0,0),oma=c(0,1,2,0),cex=1.2)

for (t in 1:7) {

  #== Protected.Ag polygons
  system("g.region res=00:01")
  system(paste("r.mapcalc 'Protected_rast=if(",type[t],">=0.5 && !(isnull(ProArea_rast)),1,null())'",sep=""))   
  system("r.to.vect -s input=Protected_rast output=Protected_vect_temp feature=area --overwrite")
  system("v.clean input=Protected_vect_temp output=Protected_vect_clean tool=rmarea thresh=5e+07 --overwrite") # Remove small areas < 10 km2
  system("v.extract -d in=Protected_vect_clean out=Protected_vect --overwrite") # Remove islands
  system("v.to.rast input=Protected_vect output=Protected_rast column=cat --overwrite") # Remove islands

  #== Import
  Protected.vect <- readVECT6(vname="Protected_vect",plugin=FALSE,mapset="ghvi")
  Protected.rast <- readRAST6(vname="Protected_rast",plugin=FALSE,mapset="ghvi")

  #== Creating nodes
  system("g.region res=00:05")
  system("r.to.vect input=Protected_rast output=Nodes feature=point --overwrite")
  #== Import into R
  Nodes <- readVECT6(vname="Nodes",plugin=FALSE,type="point",mapset="ghvi")
  #== Converting to spatstat objects
  Nodes.ss <- as(Nodes["cat"],"ppp")
  #== Computing edge with maximal distance  
  edge <- spatgraph(Nodes.ss,"geometric",par=0.15)
  Connectivity[t] <- length(unlist(edge$edges))
  #== Coordinates for text
  bb.x <- function(bb) {bb[1,1]+0.08*(bb[1,2]-bb[1,1])}
  bb.y <- function(bb) {bb[2,1]+0.90*(bb[2,2]-bb[2,1])}
    
  #== Plot
  #= Ag
  if (t %in% c(1:3)) {
    plot(Protected.vect,col=gcolors.Red(20)[15])
    # plot(Protected.vect,col="#8080804C")
    plot.sg(edge,Nodes.ss,add.points=TRUE,
            points.col="black", points.pch=19, points.cex=0.4, lines.col="black",add=TRUE)
    plot(MadaOutline,border=grey(0.5),add=TRUE)
    bb <- bbox(Protected.vect)
    text(x=bb.x(bb),y=bb.y(bb),labels=paste("C = ",Connectivity[t],sep=""),cex=0.8)
  }
  #= Ap
  if (t %in% c(4:6)) {
    plot(MadaOutline.cut,border=grey(0.5))
    plot(Protected.vect,col=gcolors.Blue(20)[15],add=TRUE)
    # plot(Protected.vect,col="#8080804C",add=TRUE) 
    plot.sg(edge,Nodes.ss,add.points=TRUE,
            points.col="black", points.pch=19, points.cex=0.4, lines.col="black")
    bb <- bbox(MadaOutline.cut)
    text(x=bb.x(bb),y=bb.y(bb),labels=paste("C = ",Connectivity[t],sep=""),cex=0.8)
  }
  #= As
  if (t==7) {
    plot(MadaOutline.cut,border=grey(0.5))
    plot(Protected.vect,col=gcolors.Green(20)[15],add=TRUE)
    # plot(Protected.vect,col="#8080804C",add=TRUE) 
    plot.sg(edge,Nodes.ss,add.points=TRUE,
            points.col="black", points.pch=19, points.cex=0.4, lines.col="black")
    bb <- bbox(MadaOutline.cut)
    text(x=bb.x(bb),y=bb.y(bb),labels=paste("C = ",Connectivity[t],sep=""),cex=0.8)
  }
  
}
#= two last empty plots for As
bb <- bbox(MadaOutline.cut)
plot(MadaOutline.cut,border=grey(0.5))
text(x=bb.x(bb),y=bb.y(bb),labels=paste("C = ",Connectivity[8],sep=""),cex=0.8)
plot(MadaOutline.cut,border=grey(0.5))
text(x=bb.x(bb),y=bb.y(bb),labels=paste("C = ",Connectivity[9],sep=""),cex=0.8)

#= Madagascar output and zooms
plot(MadaOutline,border=grey(0.5))
rect(xleft=48.42,ybottom=-13.62,xright=50.05,ytop=-11.89,col="transparent",border="black")
rect(xleft=42.70,ybottom=-23.80,xright=45.10,ytop=-19.70,col="transparent",border="black")

#=
mtext(text=c(expression(bold("PPs connectivity\n2010")),expression(bold("PPs connectivity\n2050")),
        expression(bold("PPs connectivity\n2080"))),outer=TRUE,at=c(0.16,0.46,0.76),side=3,line=0,cex=1.2,adj=0.5)
mtext(text=c(expression(italic("A. suarezensis")),expression(italic("A. perrieri")),
        expression(italic("A. grandidieri"))),outer=TRUE,at=c(0.125,0.375,0.75),side=2,line=0,cex=1.2)

dev.off()

#= Convert from pdf to png
system("convert -density 300 Plot-Connectivity.pdf Plot-Connectivity.png")

#= Backup connectivity
write.table(Connectivity,file="Connectivity.txt",row.names=FALSE,quote=FALSE,sep="\t")

#==================================================================================================
#= END OF SCRIPT
#==================================================================================================




