library(devtools)
install_github("jingshuw/grapple")
plink_refdat <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"
# set the following to wherever datasets are stored
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/UNIVERSITY/part iii/
      Essay /datasets and similar for analysis")
library(GRAPPLE)
library(TwoSampleMR)
library(ggplot2)

# this function gets the data into the right format for input into univariable 
# analysis. it returns two elements, the first can be used for TwoSampleMR 
# methods (e.g. mregger, ivw), the second is used for GRAPPLE methods such as 
# grappleRobustEst, which runs MR RAPS
preprocessforunivariable <- function(sel.filehere,exp.filehere,out.filehere,
                                     exposure,outcome){
  data.test <- getInput(sel.filehere, exp.filehere, out.filehere, 
                        plink_refdat, p.thres = 1,
                        plink_exe = '~/Downloads/plink_mac_20241022/plink')
  dataforgrap<-data.test
  data.test$data$id.exposure <- NA
  data.test$data$id.outcome <- NA
  data.test$data$exposure <- exposure
  data.test$data$outcome <- outcome
  data.test$data$mr_keep <- TRUE
  
  # the issue of why MR GRAPPLE doesn't run is the column names of input, 
  # so these are changed here
  colnames(data.test$data)[colnames(data.test$data) == "gamma_exp1"] <- 
    "beta.exposure"
  colnames(data.test$data)[colnames(data.test$data) == "se_exp1"] <- 
    "se.exposure"
  colnames(data.test$data)[colnames(data.test$data) == "gamma_out1"] <- 
    "beta.outcome"
  colnames(data.test$data)[colnames(data.test$data) == "se_out1"] <- 
    "se.outcome"
  
  # View the updated dataframe
  print(head(data.test$data))
  return(list(data.test,dataforgrap) )
}


# the function below runs univariable MR for a selection of methods.
# the data are inputted in two formats (one for grapple, one for the rest)
# for grapple the output has to be reformatted to fit with the rest 
# it allows input of a vector of p value selection criteria and outputs
# a dataframe of results and prints estimates with their uncertainties
rununivariableMRs <- function(data.test,dataforgrap,pvalsels,exposure,outcome){
  
  # analyse with selection of univariable MR methods
  methods <- c("mr_egger_regression", "mr_weighted_median", "mr_ivw")
  parameters <- default_parameters()
  parameters$shrinkage <- TRUE
  out.all <- data.frame()
  
  for (pval.sel in pvalsels) {
    out <- mr(subset(data.test$data, selection_pvals< pval.sel),
              parameters = parameters,method_list = methods)
    grapp <- grappleRobustEst(dataforgrap$data,p.thres = pval.sel,
                              diagnosis = FALSE)
    grapprow <- c(NA,NA,outcome,exposure,"RAPS (via GRAPPLE)",out$nsnp[1],
                  grapp$beta.hat,sqrt(grapp$beta.var),
                  grapp$beta.p.value)
    out <- rbind(out,grapprow)
    out$pval.sel <- pval.sel
    out.all <- rbind(out.all, out)
  }
  
  out.all$pval.adjusted <- p.adjust(out.all$pval, "bonferroni")
  levels(out.all$method) <- c("IVW", "Egger", "RAPS", "W. Median")
  out.all$se <- as.numeric(out.all$se)
  out.all$b <- as.numeric(out.all$b)
  out.all$pval<-as.numeric(out.all$pval)
  out.all$nsnp<- as.numeric(out.all$nsnp)
  library(ggplot2)
  plot.new
  print(ggplot(out.all) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                              ymax = b + 1.96 * se,
                              col = pval.adjusted < 0.05) +
          geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                      linetype = "dashed") +
          facet_grid(exposure ~ pval.sel) + coord_flip() +
          theme_light(base_size = 18) + theme(legend.position = "bottom")) 
  print(out.all)
  return(out.all)
}

# check that the function above produces the same results as the GRAPPLE paper, 
# for CRP on CAD. All match except causal effects of CRP on CAD with
# RAPS, at 0.001 level the paper suggests 0.014, and this analysis gets 0.032. 
# They would match if 0.0135 was misread for 0.0315
sel.file <- "datasets_grapple/CRP-Prins17.csv"
exp.file <- "datasets_grapple/CRP-Dehghan11.csv"
out.file <- "datasets_grapple/CAD-Nelson17.csv"
crpcaddata <- preprocessforunivariable(sel.file,exp.file,out.file,"CRP","CAD")
rununivariableMRs(crpcaddata[1][[1]],crpcaddata[2][[1]],c(1e-8,1e-5,1e-3),
                  "CRP","CAD")


# load in data for new analyses and check they have the right column contents 
# and column names, renaming where necessary. The following datasets were 
# downloaded from the GWAS summary statistics
# as outlined in the manuscript of the part III essay. 
# Datasets downloaded from the GRAPPLE GitHub do not need editing.


APOBdavyson23 <- read.delim("APOB-davyson23.tsv")
APOBkarjalainen24 <- read.delim("APOB-karjalainen24.tsv")
APOBgudjonsson22 <- read.delim("APOB-Gudjonsson22.tsv")
APOBkettunen16 <- read.delim("APOB-Kettunen16.tsv")
HDLgraham21 <- read.delim("HDL-graham21.tsv")
CADkoskeridis24 <- read.delim("Koskeridis2024.tsv")
LDLgraham21 <- read.delim("LDL-graham21.tsv")
TGgraham21 <- read.delim("TG-graham21.tsv")
ISmishra22 <- read.delim("IS-mishra22.tsv")
ISmalik18<-read.delim("IS-malik18.tsv")
APOA1davyson23 <- read.delim("APOA1-davyson23.tsv")
APOA1kettunen16 <- read.delim("APOA1-Kettunen16.tsv")
smallHDLkettunen16 <- read.delim("smallHDL-kettunen16.tsv")
smallHDLdavyson23 <- read.delim("smallHDL-davyson23.tsv")
mediumHDLkettunen16 <- read.delim("mediumHDL-kettunen16.tsv")
mediumHDLdavyson23 <- read.delim("mediumHDL-davyson23.tsv")
largeHDLkettunen16 <- read.delim("largeHDL-kettunen16.tsv")
largeHDLdavyson23 <- read.delim("largeHDL-davyson23.tsv")

# check dataframes
head(APOBdavyson23)
head(APOBkarjalainen24)
# karjalainen does not have the right format of SNP IDs required
head(APOBgudjonsson22)
head(APOBkettunen16)
head(HDLgraham21)
head(CADkoskeridis24)
head(LDLgraham21)
head(TGgraham21)
head(ISmishra22)
# Mishra doesn't have the right format of SNP IDs required
head(ISmalik18)
head(APOA1davyson23)
head(APOA1kettunen16)
head(smallHDLkettunen16)
head(smallHDLdavyson23)
head(mediumHDLkettunen16)
head(mediumHDLdavyson23)
head(largeHDLkettunen16)
head(largeHDLdavyson23)

#change column names
colnames(APOBgudjonsson22)[colnames(APOBgudjonsson22)=="standard_error"]<- "se"
colnames(APOBgudjonsson22)[colnames(APOBgudjonsson22)=="p_value"]<- "pval"
colnames(APOBgudjonsson22)[colnames(APOBgudjonsson22)=="variant_id"]<-"SNP"
colnames(APOBdavyson23)[colnames(APOBdavyson23)=="rs_id"]<-"SNP"
colnames(APOBdavyson23)[colnames(APOBdavyson23)=="p_value"]<-"pval"
colnames(APOBdavyson23)[colnames(APOBdavyson23)=="standard_error"]<-"se"
colnames(APOBkettunen16)[colnames(APOBkettunen16)=="rs_id"]<-"SNP"
colnames(APOBkettunen16)[colnames(APOBkettunen16)=="p_value"]<-"pval"
colnames(APOBkettunen16)[colnames(APOBkettunen16)=="standard_error"]<-"se"
colnames(HDLgraham21)[colnames(HDLgraham21)=="standard_error"]<- "se"
colnames(HDLgraham21)[colnames(HDLgraham21)=="p_value"]<- "pval"
colnames(HDLgraham21)[colnames(HDLgraham21)=="variant_id"]<-"SNP"
colnames(CADkoskeridis24)[colnames(CADkoskeridis24)=="standard_error"]<- "se"
colnames(CADkoskeridis24)[colnames(CADkoskeridis24)=="p_value"]<- "pval"
colnames(CADkoskeridis24)[colnames(CADkoskeridis24)=="rs_id"]<-"SNP"
colnames(LDLgraham21)[colnames(LDLgraham21)=="standard_error"]<- "se"
colnames(LDLgraham21)[colnames(LDLgraham21)=="p_value"]<- "pval"
colnames(LDLgraham21)[colnames(LDLgraham21)=="variant_id"]<-"SNP"
colnames(TGgraham21)[colnames(TGgraham21)=="standard_error"]<- "se"
colnames(TGgraham21)[colnames(TGgraham21)=="p_value"]<- "pval"
colnames(TGgraham21)[colnames(TGgraham21)=="variant_id"]<-"SNP"
colnames(ISmalik18)[colnames(ISmalik18)=="standard_error"]<- "se"
colnames(ISmalik18)[colnames(ISmalik18)=="p_value"]<- "pval"
colnames(ISmalik18)[colnames(ISmalik18)=="variant_id"]<-"SNP"
colnames(APOA1davyson23)[colnames(APOA1davyson23)=="rs_id"]<-"SNP"
colnames(APOA1davyson23)[colnames(APOA1davyson23)=="p_value"]<-"pval"
colnames(APOA1davyson23)[colnames(APOA1davyson23)=="standard_error"]<-"se"
colnames(APOA1kettunen16)[colnames(APOA1kettunen16)=="rs_id"]<-"SNP"
colnames(APOA1kettunen16)[colnames(APOA1kettunen16)=="p_value"]<-"pval"
colnames(APOA1kettunen16)[colnames(APOA1kettunen16)=="standard_error"]<-"se"
colnames(smallHDLdavyson23)[colnames(smallHDLdavyson23)=="rs_id"]<-"SNP"
colnames(smallHDLdavyson23)[colnames(smallHDLdavyson23)=="p_value"]<-"pval"
colnames(smallHDLdavyson23)[colnames(smallHDLdavyson23)=="standard_error"]<-"se"
colnames(smallHDLkettunen16)[colnames(smallHDLkettunen16)=="rs_id"]<-"SNP"
colnames(smallHDLkettunen16)[colnames(smallHDLkettunen16)=="p_value"]<-"pval"
colnames(smallHDLkettunen16)[colnames(smallHDLkettunen16)=="standard_error"]<-
  "se"
colnames(mediumHDLdavyson23)[colnames(mediumHDLdavyson23)=="rs_id"]<-"SNP"
colnames(mediumHDLdavyson23)[colnames(mediumHDLdavyson23)=="p_value"]<-"pval"
colnames(mediumHDLdavyson23)[colnames(mediumHDLdavyson23)=="standard_error"]<-
  "se"
colnames(mediumHDLkettunen16)[colnames(mediumHDLkettunen16)=="rs_id"]<-"SNP"
colnames(mediumHDLkettunen16)[colnames(mediumHDLkettunen16)=="p_value"]<-"pval"
colnames(mediumHDLkettunen16)[colnames(mediumHDLkettunen16)=="standard_error"]<-
  "se"
colnames(largeHDLdavyson23)[colnames(largeHDLdavyson23)=="rs_id"]<-"SNP"
colnames(largeHDLdavyson23)[colnames(largeHDLdavyson23)=="p_value"]<-"pval"
colnames(largeHDLdavyson23)[colnames(largeHDLdavyson23)=="standard_error"]<-"se"
colnames(largeHDLkettunen16)[colnames(largeHDLkettunen16)=="rs_id"]<-"SNP"
colnames(largeHDLkettunen16)[colnames(largeHDLkettunen16)=="p_value"]<-"pval"
colnames(largeHDLkettunen16)[colnames(largeHDLkettunen16)=="standard_error"]<-"se"


#save back to file
write.csv(APOBgudjonsson22,"APOB-gudjonsson22new.csv", row.names = FALSE)
write.csv(APOBdavyson23,"APOB-davyson23new.csv", row.names = FALSE)
write.csv(APOBkettunen16,"APOB-kettunen16new.csv",row.names = FALSE)
write.csv(HDLgraham21,"HDL-graham21new.csv", row.names = FALSE)
write.csv(CADkoskeridis24,"CAD-koskeridis24new.csv",row.names = FALSE)
write.csv(LDLgraham21,"LDL-graham21new.csv", row.names = FALSE)
write.csv(TGgraham21,"TG-graham21new.csv", row.names = FALSE)
write.csv(ISmalik18,"IS-malik18new.csv", row.names = FALSE)
write.csv(APOA1davyson23,"APOA1-davyson23new.csv", row.names = FALSE)
write.csv(APOA1kettunen16,"APOA1-kettunen16new.csv",row.names = FALSE)
write.csv(smallHDLdavyson23,"smallHDL-davyson23new.csv", row.names = FALSE)
write.csv(smallHDLkettunen16,"smallHDL-kettunen16new.csv",row.names = FALSE)
write.csv(mediumHDLdavyson23,"mediumHDL-davyson23new.csv", row.names = FALSE)
write.csv(mediumHDLkettunen16,"mediumHDL-kettunen16new.csv",row.names = FALSE)
write.csv(largeHDLdavyson23,"largeHDL-davyson23new.csv", row.names = FALSE)
write.csv(largeHDLkettunen16,"largeHDL-kettunen16new.csv",row.names = FALSE)

##################################
# Analysis of lipids on CAD

# univariable analyses of HDL on CAD
sel.file.hdl <- "datasets_grapple/HDL-gera18.csv"
exp.file.hdl <- "HDL-graham21new.csv"
out.file.hdl <- "CAD-koskeridis24new.csv"
hdlcaddata <- preprocessforunivariable(sel.file.hdl,exp.file.hdl,out.file.hdl,
                                       "HDL","CAD")
univrunHDLCAD <- rununivariableMRs(hdlcaddata[1][[1]],hdlcaddata[2][[1]],
                                   c(1e-8,1e-5,1e-3),"HDL","CAD")
#seek to find diagnostic plots for RAPS
hdlcadforplots <- grappleRobustEst(hdlcaddata[2][[1]]$data,p.thres = 1e-05,
                                   diagnosis = TRUE)
# grapple methodology: seek modes in profile likelihood. The p-value threshold 
# can be adjusted to look for further modes at different thresholds
diagnosis <- findModes(hdlcaddata[2][[1]]$data, p.thres = 1e-8, 
                       marker.data = hdlcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers

# try an analysis of APOB on CAD with gudjonsson, davyson, koskeridis. 
sel.file.apob <- "APOB-gudjonsson22new.csv"
exp.file.apob <- "APOB-davyson23new.csv"
out.file.apob <- "CAD-koskeridis24new.csv"
apobcaddata <- preprocessforunivariable(sel.file.apob,exp.file.apob,
                                        out.file.apob,"APOB","CAD")
rununivariableMRs(apobcaddata[1][[1]],apobcaddata[2][[1]],c(1e-8,1e-5,1e-3),
                  "APOB","CAD")
# run grapple and look for modes: significant multimodality is found
diagnosis <- findModes(apobcaddata[2][[1]]$data, p.thres = 1e-5, 
                       marker.data = apobcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("apob-cad-profile-10-5.pdf")

# try an analysis of APOB on CAD with kettunen, davyson, koskeridis
sel.file.apob2 <- "APOB-kettunen16new.csv"
exp.file.apob <- "APOB-davyson23new.csv"
out.file.apob <- "CAD-koskeridis24new.csv"
apobcaddata2 <- preprocessforunivariable(sel.file.apob2,exp.file.apob,
                                         out.file.apob,"APOB","CAD")
univrunAPOBCAD <-rununivariableMRs(apobcaddata2[1][[1]],apobcaddata2[2][[1]],
                                   c(1e-8,1e-5,1e-3),"APOB","CAD")
#seek to find diagnostic plots for RAPS
apobcadforplots <- grappleRobustEst(apobcaddata[2][[1]]$data,p.thres = 1e-05,
                                    diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(apobcaddata2[2][[1]]$data, p.thres = 1e-3, 
                       marker.data = apobcaddata2[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers

# try an analysis of APOA1 on CAD with kettunen, davyson, koskeridis
sel.file.apoa1 <- "APOA1-kettunen16new.csv"
exp.file.apoa1 <- "APOA1-davyson23new.csv"
out.file.apoa1 <- "CAD-koskeridis24new.csv"
apoa1caddata <- preprocessforunivariable(sel.file.apoa1,exp.file.apoa1,
                                         out.file.apoa1,"APOA1","CAD")
univrunAPOA1CAD <-rununivariableMRs(apoa1caddata[1][[1]],apoa1caddata[2][[1]],
                                    c(1e-8,1e-5,1e-3),"APOA1","CAD")
#seek to find diagnostic plots for RAPS
apoa1cadforplots <- grappleRobustEst(apoa1caddata[2][[1]]$data,p.thres = 1e-05,
                                     diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(apoa1caddata[2][[1]]$data, p.thres = 1e-8, 
                       marker.data = apoa1caddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("apoa1_profile_10-3.pdf")

# analysis of small HDL concentration on CAD with kettunen, davyson, koskeridis
sel.file.hdl.small <- "smallHDL-kettunen16new.csv"
exp.file.hdl.small <- "smallHDL-davyson23new.csv"
out.file.hdl.small <- "CAD-koskeridis24new.csv"
smallhdlcaddata <- preprocessforunivariable(sel.file.hdl.small,
                                            exp.file.hdl.small,
                                            out.file.hdl.small,
                                            "small_HDL","CAD")
univrunsmallHDLCAD <-rununivariableMRs(smallhdlcaddata[1][[1]],
                                       smallhdlcaddata[2][[1]],
                                       c(1e-8,1e-5,1e-3),"small_HDL","CAD")
ggsave("smallhdl-cad-results.pdf")
#seek to find diagnostic plots for RAPS
smallHDLcadforplots <- grappleRobustEst(smallhdlcaddata[2][[1]]$data,
                                        p.thres = 1e-05,diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(smallhdlcaddata[2][[1]]$data, p.thres = 1e-5, 
                       marker.data = smallhdlcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("smallhdl_cad_profile_10-5.pdf")

# analysis of medium HDL concentration on CAD with kettunen, davyson, koskeridis
sel.file.hdl.medium <- "mediumHDL-kettunen16new.csv"
exp.file.hdl.medium <- "mediumHDL-davyson23new.csv"
out.file.hdl.medium <- "CAD-koskeridis24new.csv"
mediumhdlcaddata <- preprocessforunivariable(sel.file.hdl.medium,
                                             exp.file.hdl.medium,
                                             out.file.hdl.medium,
                                             "medium_HDL","CAD")
univrunmediumHDLCAD <-rununivariableMRs(mediumhdlcaddata[1][[1]],
                                        mediumhdlcaddata[2][[1]],
                                        c(1e-8,1e-5,1e-3),"medium_HDL","CAD")
ggsave("mediumhdl_cad_results.pdf")
#seek to find diagnostic plots for RAPS
mediumHDLcadforplots <- grappleRobustEst(mediumhdlcaddata[2][[1]]$data,
                                         p.thres = 1e-05,diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(mediumhdlcaddata[2][[1]]$data, p.thres = 1e-3, 
                       marker.data = mediumhdlcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("mediumhdl_cad_profile_10-3.pdf")

# analysis of large HDL concentration on CAD with kettunen, davyson, koskeridis
sel.file.hdl.large <- "largeHDL-kettunen16new.csv"
exp.file.hdl.large <- "largeHDL-davyson23new.csv"
out.file.hdl.large <- "CAD-koskeridis24new.csv"
largehdlcaddata <- preprocessforunivariable(sel.file.hdl.large,
                                            exp.file.hdl.large,
                                            out.file.hdl.large,
                                            "large_HDL","CAD")
univrunlargeHDLCAD <-rununivariableMRs(largehdlcaddata[1][[1]],
                                       largehdlcaddata[2][[1]],
                                       c(1e-8,1e-5,1e-3),"large_HDL","CAD")
ggsave("largehdl_cad_results.pdf")
#seek to find diagnostic plots for RAPS
largeHDLcadforplots <- grappleRobustEst(largehdlcaddata[2][[1]]$data,
                                        p.thres = 1e-05,diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(largehdlcaddata[2][[1]]$data, p.thres = 1e-5, 
                       marker.data = largehdlcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("largehdl_cad_profile_10-5.pdf")

# univariable analyses of HDL on CAD
sel.file.hdl <- "datasets_grapple/HDL-gera18.csv"
exp.file.hdl <- "HDL-graham21new.csv"
out.file.hdl <- "CAD-koskeridis24new.csv"
hdlcaddata <- preprocessforunivariable(sel.file.hdl,exp.file.hdl,out.file.hdl,
                                       "HDL","CAD")
univrunHDLCAD <- rununivariableMRs(hdlcaddata[1][[1]],hdlcaddata[2][[1]],
                                   c(1e-8,1e-5,1e-3),"HDL","CAD")
#seek to find diagnostic plots for RAPS
hdlcadforplots <- grappleRobustEst(hdlcaddata[2][[1]]$data,p.thres = 1e-05,
                                   diagnosis = TRUE)
# grapple methodology: seek modes in profile likelihood. The p-value threshold 
# can be adjusted to look for further modes at different thresholds
diagnosis <- findModes(hdlcaddata[2][[1]]$data, p.thres = 1e-8, 
                       marker.data = hdlcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers


# univariable analyses of LDL on CAD
sel.file.ldl <- "datasets_grapple/LDL-gera18.csv"
exp.file.ldl <- "LDL-graham21new.csv"
out.file.ldl <- "CAD-koskeridis24new.csv"
ldlcaddata <- preprocessforunivariable(sel.file.ldl,exp.file.ldl,out.file.ldl,
                                       "LDL","CAD")
univrunLDLCAD <- rununivariableMRs(ldlcaddata[1][[1]],ldlcaddata[2][[1]],
                                   c(1e-8,1e-5,1e-3),"LDL","CAD")
#seek to find diagnostic plots for RAPS
ldlcadforplots <- grappleRobustEst(ldlcaddata[2][[1]]$data,p.thres = 1e-05,
                                   diagnosis = TRUE)
#grapple methodology: seek modes in profile likelihood
diagnosis <- findModes(ldlcaddata[2][[1]]$data, p.thres = 1e-5, 
                       marker.data = ldlcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers

# univariable analyses of TG on CAD
sel.file.tg <- "datasets_grapple/TG-gera18.csv"
exp.file.tg <- "TG-graham21new.csv"
out.file.tg <- "CAD-koskeridis24new.csv"
tgcaddata <- preprocessforunivariable(sel.file.tg,exp.file.tg,out.file.tg,
                                      "TG","CAD")
univrunTGCAD <- rununivariableMRs(tgcaddata[1][[1]],tgcaddata[2][[1]],
                                  c(1e-8,1e-5,1e-3),"TG","CAD")
#seek to find diagnostic plots for RAPS
tgcadforplots <- grappleRobustEst(tgcaddata[2][[1]]$data,p.thres = 1e-05,
                                  diagnosis = TRUE)
# grapple methodology: seek modes in profile likelihood. The p-value threshold
# can be adjusted to look for further modes at different thresholds
diagnosis <- findModes(tgcaddata[2][[1]]$data, p.thres = 1e-8, 
                       marker.data = tgcaddata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers

# multivariable analysis of LDL, HDL and TG on CAD
sel.file.many <- c("datasets_grapple/LDL-gera18.csv",
                   "datasets_grapple/HDL-gera18.csv",
                   "datasets_grapple/TG-gera18.csv")
exp.file.many <- c("LDL-graham21new.csv","HDL-graham21new.csv",
                   "TG-graham21new.csv")
out.file.many <- "CAD-koskeridis24new.csv"
multicaddata <- getInput(sel.file.many, exp.file.many, out.file.many, 
                         plink_refdat, max.p.thres = 0.01,
                         plink_exe = '~/Downloads/plink_mac_20241022/plink')
multiruncad <- grappleRobustEst(multicaddata$data, p.thres = 1e-4)
multiruncad2 <-grappleRobustEst(multicaddata$data, p.thres = 1e-2)
multiruncad3 <- grappleRobustEst(multicaddata$data, p.thres = 1e-6)

# plot results of multivariate LDL HDL TG on CAD
dattoplot <- data.frame(b=c(multiruncad$beta.hat,multiruncad2$beta.hat,
                            multiruncad3$beta.hat),
                        se=c(sqrt(diag(multiruncad$beta.var)),
                             sqrt(diag(multiruncad2$beta.var)),
                             sqrt(diag(multiruncad3$beta.var))),
                        pval=c(multiruncad$beta.p.value,
                               multiruncad2$beta.p.value,
                               multiruncad3$beta.p.value),
                        exposure=rep(c("LDL","HDL","TG"),3),
                        method = rep("Multivariable GRAPPLE",9),
                        pval.sel = c(rep(10^(-4),3),rep(10^(-2),3),
                                     rep(10^(-6),3)),
                        pval.adjusted= rep(0,9))
dattoplot$pval.adjusted <- p.adjust(dattoplot$pval,"bonferroni")
print(ggplot(dattoplot) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                              ymax = b + 1.96 * se,
                              col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot,"ldl-hdl-tg-cad-multivar-results.csv")

# run multivariate APOB and LDL on CAD
sel.file.apob.ldl <-c("APOB-kettunen16new.csv", 
                      "datasets_grapple/LDL-gera18.csv")
exp.file.apob.ldl <- c("APOB-davyson23new.csv", "LDL-graham21new.csv")
out.file.apob.ldl <- "CAD-koskeridis24new.csv"
ldlapobcaddata <- getInput(sel.file.apob.ldl, exp.file.apob.ldl, 
                           out.file.apob.ldl, plink_refdat, max.p.thres = 0.01,
                           plink_exe = '~/Downloads/plink_mac_20241022/plink')
multirunapoldlcad <- grappleRobustEst(ldlapobcaddata$data, p.thres = 1e-4)
multirunapoldlcad2 <- grappleRobustEst(ldlapobcaddata$data,p.thres = 1e-6)
multirunapoldlcad3 <- grappleRobustEst(ldlapobcaddata$data,p.thres = 1e-3)
# plot multivariate results of APOB and LDL on CAD
dattoplotall <- data.frame(b=c(multirunapoldlcad$beta.hat,
                               multirunapoldlcad2$beta.hat,
                               multirunapoldlcad3$beta.hat),
                           se=c(sqrt(diag(multirunapoldlcad$beta.var)),
                                sqrt(diag(multirunapoldlcad2$beta.var)),
                                sqrt(diag(multirunapoldlcad3$beta.var))),
                           pval=c(multirunapoldlcad$beta.p.value,
                                  multirunapoldlcad2$beta.p.value,
                                  multirunapoldlcad3$beta.p.value),
                           exposure=rep(c("ApoB","LDL"),3),
                           method = rep("Multivariable GRAPPLE",6),
                           pval.sel = c(rep(10^(-4),2),rep(10^(-6),2),
                                        rep(10^(-3),2)),
                           pval.adjusted= rep(0,6))
dattoplotall$pval.adjusted <- p.adjust(dattoplotall$pval,"bonferroni")
library(ggplot2)
print(ggplot(dattoplotall) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                                 ymax = b + 1.96 * se,
                                 col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplotall,"apob-ldl-cad-multivar-results.csv")
ggsave("apob-ldl-cad-multivar-plot.pdf")

# HDL and ApoA1 on "CAD-koskeridis24new.csv"
sel.file.apoa1.hdl <-c("APOA1-kettunen16new.csv", 
                       "datasets_grapple/HDL-gera18.csv")
exp.file.apoa1.hdl <- c("APOA1-davyson23new.csv", "HDL-graham21new.csv")
out.file.cad <- "CAD-koskeridis24new.csv"
hdlapoa1caddata <- getInput(sel.file.apoa1.hdl, exp.file.apoa1.hdl, 
                            out.file.cad, plink_refdat, max.p.thres = 0.01,
                            plink_exe = '~/Downloads/plink_mac_20241022/plink')
multirunapohdlcad <- grappleRobustEst(hdlapoa1caddata$data, p.thres = 1e-4)
multirunapohdlcad2 <- grappleRobustEst(hdlapoa1caddata$data,p.thres = 1e-6)
multirunapohdlcad3 <- grappleRobustEst(hdlapoa1caddata$data,p.thres = 1e-2)
multirunapohdlcad4 <- mutlirunapohdlcad3
dattoplot2 <- data.frame(b=c(multirunapohdlcad$beta.hat,
                             multirunapohdlcad2$beta.hat,
                             multirunapohdlcad3$beta.hat),
                         se=c(sqrt(diag(multirunapohdlcad$beta.var)),
                              sqrt(diag(multirunapohdlcad2$beta.var)),
                              sqrt(diag(multirunapohdlcad3$beta.var))),
                         pval=c(multirunapohdlcad$beta.p.value,
                                multirunapohdlcad2$beta.p.value,
                                multirunapohdlcad3$beta.p.value),
                         exposure=rep(c("ApoA1","HDL"),3),
                         method = rep("Multivariable GRAPPLE",6),
                         pval.sel = c(rep(10^(-4),2),rep(10^(-6),2),
                                      rep(10^(-2),2)),
                         pval.adjusted= rep(0,6))
dattoplot2$pval.adjusted <- p.adjust(dattoplot2$pval,"bonferroni")
library(ggplot2)
print(ggplot(dattoplot2) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot2, "apoa1-hdl-cad-results.csv", row.names=FALSE)
ggsave("apoa1-hdl-cad-multivar-plot.pdf")


#multivariate small HDL with LDL, HDL and TG on CAD
sel.file.smallhdl.hdl.ldl.tg <- c("smallHDL-kettunen16new.csv",
                                  "datasets_grapple/LDL-gera18.csv",
                                  "datasets_grapple/HDL-gera18.csv",
                   "datasets_grapple/TG-gera18.csv")
exp.file.smallhdl.hdl.ldl.tg <- c("smallHDL-davyson23new.csv",
                                  "LDL-graham21new.csv","HDL-graham21new.csv",
                                  "TG-graham21new.csv")
out.file.cad <- "CAD-koskeridis24new.csv"
smallhdl_multi_caddata <- getInput(sel.file.smallhdl.hdl.ldl.tg, 
                                   exp.file.smallhdl.hdl.ldl.tg, 
                                   out.file.cad, plink_refdat, 
                                   max.p.thres = 0.01,
                                   plink_exe = 
                                     '~/Downloads/plink_mac_20241022/plink')
smallhdl_multi_cad <- grappleRobustEst(smallhdl_multi_caddata$data, 
                                       p.thres = 1e-4)
smallhdl_multi_cad2 <-grappleRobustEst(smallhdl_multi_caddata$data, 
                                       p.thres = 1e-2)
smallhdl_multi_cad3 <- grappleRobustEst(smallhdl_multi_caddata$data, 
                                        p.thres = 1e-6)

# plot results of multivariate small HDL LDL HDL TG on CAD
dattoplot4 <- data.frame(b=c(smallhdl_multi_cad$beta.hat,
                             smallhdl_multi_cad2$beta.hat,
                             smallhdl_multi_cad3$beta.hat),
                        se=c(sqrt(diag(smallhdl_multi_cad$beta.var)),
                             sqrt(diag(smallhdl_multi_cad2$beta.var)),
                             sqrt(diag(smallhdl_multi_cad3$beta.var))),
                        pval=c(smallhdl_multi_cad$beta.p.value,
                               smallhdl_multi_cad2$beta.p.value,
                               smallhdl_multi_cad3$beta.p.value),
                        exposure=rep(c("smallHDL","LDL","HDL","TG"),3),
                        method = rep("Multivariable GRAPPLE",12),
                        pval.sel = c(rep(10^(-4),4),rep(10^(-2),4),
                                     rep(10^(-6),4)),
                        pval.adjusted= rep(0,12))
dattoplot4$pval.adjusted <- p.adjust(dattoplot4$pval,"bonferroni")
print(ggplot(dattoplot4) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                              col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot4,"smallhdl-ldl-hdl-tg-cad-multivar-results.csv")
ggsave("smallhdl-ldl-hdl-tg-cad-multivar-plot.pdf")

#multivariate medium HDL with LDL, HDL and TG on CAD
sel.file.mediumhdl.hdl.ldl.tg <- c("mediumHDL-kettunen16new.csv",
                                   "datasets_grapple/LDL-gera18.csv",
                                   "datasets_grapple/HDL-gera18.csv",
                                  "datasets_grapple/TG-gera18.csv")
exp.file.mediumhdl.hdl.ldl.tg <- c("mediumHDL-davyson23new.csv",
                                   "LDL-graham21new.csv","HDL-graham21new.csv",
                                   "TG-graham21new.csv")
out.file.cad <- "CAD-koskeridis24new.csv"
mediumhdl_multi_caddata <- getInput(sel.file.mediumhdl.hdl.ldl.tg, 
                                    exp.file.mediumhdl.hdl.ldl.tg, 
                                    out.file.cad, plink_refdat, 
                                    max.p.thres = 0.01,plink_exe =
                                      '~/Downloads/plink_mac_20241022/plink')
mediumhdl_multi_cad <- grappleRobustEst(mediumhdl_multi_caddata$data,
                                        p.thres = 1e-4)
mediumhdl_multi_cad2 <-grappleRobustEst(mediumhdl_multi_caddata$data, 
                                        p.thres = 1e-2)
mediumhdl_multi_cad3 <- grappleRobustEst(mediumhdl_multi_caddata$data,
                                         p.thres = 1e-6)

# plot results of multivariate medium HDL LDL HDL TG on CAD
dattoplot5 <- data.frame(b=c(mediumhdl_multi_cad$beta.hat,
                             mediumhdl_multi_cad2$beta.hat,
                             mediumhdl_multi_cad3$beta.hat),
                         se=c(sqrt(diag(mediumhdl_multi_cad$beta.var)),
                              sqrt(diag(mediumhdl_multi_cad2$beta.var)),
                              sqrt(diag(mediumhdl_multi_cad3$beta.var))),
                         pval=c(mediumhdl_multi_cad$beta.p.value,
                                mediumhdl_multi_cad2$beta.p.value,
                                mediumhdl_multi_cad3$beta.p.value),
                         exposure=rep(c("mediumHDL","LDL","HDL","TG"),3),
                         method = rep("Multivariable GRAPPLE",12),
                         pval.sel = c(rep(10^(-4),4),rep(10^(-2),4),
                                      rep(10^(-6),4)),
                         pval.adjusted= rep(0,12))
dattoplot5$pval.adjusted <- p.adjust(dattoplot5$pval,"bonferroni")
print(ggplot(dattoplot5) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot5,"mediumhdl-ldl-hdl-tg-cad-multivar-results.csv")
ggsave("mediumhdl-ldl-hdl-tg-cad-multivar-plot.pdf")


#multivariate large HDL with LDL, HDL and TG
sel.file.largehdl.hdl.ldl.tg <- c("largeHDL-kettunen16new.csv",
                                  "datasets_grapple/LDL-gera18.csv",
                                  "datasets_grapple/HDL-gera18.csv",
                                   "datasets_grapple/TG-gera18.csv")
exp.file.largehdl.hdl.ldl.tg <- c("largeHDL-davyson23new.csv",
                                  "LDL-graham21new.csv","HDL-graham21new.csv",
                                  "TG-graham21new.csv")
out.file.cad <- "CAD-koskeridis24new.csv"
largehdl_multi_caddata <- getInput(sel.file.largehdl.hdl.ldl.tg, 
                                   exp.file.largehdl.hdl.ldl.tg, 
                                   out.file.cad, plink_refdat, 
                                   max.p.thres = 0.01,plink_exe =
                                     '~/Downloads/plink_mac_20241022/plink')
largehdl_multi_cad <- grappleRobustEst(largehdl_multi_caddata$data, 
                                       p.thres = 1e-4)
largehdl_multi_cad2 <-grappleRobustEst(largehdl_multi_caddata$data, 
                                       p.thres = 1e-2)
largehdl_multi_cad3 <- grappleRobustEst(largehdl_multi_caddata$data, 
                                        p.thres = 1e-6)

# plot results of multivariate large HDL LDL HDL TG on CAD
dattoplot5 <- data.frame(b=c(largehdl_multi_cad$beta.hat,
                             largehdl_multi_cad2$beta.hat,
                             largehdl_multi_cad3$beta.hat),
                         se=c(sqrt(diag(largehdl_multi_cad$beta.var)),
                              sqrt(diag(largehdl_multi_cad2$beta.var)),
                              sqrt(diag(largehdl_multi_cad3$beta.var))),
                         pval=c(largehdl_multi_cad$beta.p.value,
                                largehdl_multi_cad2$beta.p.value,
                                largehdl_multi_cad3$beta.p.value),
                         exposure=rep(c("largeHDL","LDL","HDL","TG"),3),
                         method = rep("Multivariable GRAPPLE",12),
                         pval.sel = c(rep(10^(-4),4),rep(10^(-2),4),
                                      rep(10^(-6),4)),
                         pval.adjusted= rep(0,12))
dattoplot6$pval.adjusted <- p.adjust(dattoplot6$pval,"bonferroni")
print(ggplot(dattoplot6) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot6,"largehdl-ldl-hdl-tg-cad-multivar-results.csv")
ggsave("largehdl-ldl-hdl-tg-cad-multivar-plot.pdf")



###################################
# we now analyse ischaemic stroke (which might display similar behaviour to CAD)
#HDL on IS
sel.file.hdl <- "datasets_grapple/HDL-gera18.csv"
exp.file.hdl <- "HDL-graham21new.csv"
out.file.is <- "IS-malik18new.csv"
hdlisdata <- preprocessforunivariable(sel.file.hdl,exp.file.hdl,out.file.is,
                                      "HDL","IS")
univrunHDLIS <- rununivariableMRs(hdlisdata[1][[1]],hdlisdata[2][[1]],
                                  c(1e-8,1e-5,1e-3),"HDL","IS")
#seek to find diagnostic plots for RAPS
hdlisforplots <- grappleRobustEst(hdlisdata[2][[1]]$data,p.thres = 1e-05,
                                  diagnosis = TRUE)
#grapple methodology: seek modes in profile likelihood
diagnosis <- findModes(hdlisdata[2][[1]]$data, p.thres = 1e-6, 
                       marker.data = hdlisdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers

# analysis of small HDL concentration on IS with kettunen, davyson, malik
sel.file.hdl.small <- "smallHDL-kettunen16new.csv"
exp.file.hdl.small <- "smallHDL-davyson23new.csv"
out.file.is <- "IS-malik18new.csv"
smallhdlisdata <- preprocessforunivariable(sel.file.hdl.small,
                                           exp.file.hdl.small,
                                           out.file.is,"small_HDL","IS")
univrunsmallHDLIS <-rununivariableMRs(smallhdlisdata[1][[1]],
                                      smallhdlisdata[2][[1]],c(1e-8,1e-5,1e-3),
                                      "small_HDL","IS")
ggsave("smallhdl_results_is.pdf")
#seek to find diagnostic plots for RAPS
smallHDLisforplots <- grappleRobustEst(smallhdlisdata[2][[1]]$data,
                                       p.thres = 1e-05,diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(smallhdlisdata[2][[1]]$data, p.thres = 1e-5, 
                       marker.data = smallhdlisdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("smallhdl_is_profile_10-5.pdf")

# analysis of medium HDL concentration on IS with kettunen, davyson, malik
sel.file.hdl.medium <- "mediumHDL-kettunen16new.csv"
exp.file.hdl.medium <- "mediumHDL-davyson23new.csv"
out.file.is <- "IS-malik18new.csv"
mediumhdlisdata <- preprocessforunivariable(sel.file.hdl.medium,
                                            exp.file.hdl.medium,
                                            out.file.is,"medium_HDL","IS")
univrunmediumHDLIS <-rununivariableMRs(mediumhdlisdata[1][[1]],
                                       mediumhdlisdata[2][[1]],
                                       c(1e-8,1e-5,1e-3),"medium_HDL","IS")
ggsave("mediumhdl_is_results.pdf")
#seek to find diagnostic plots for RAPS
mediumHDLisforplots <- grappleRobustEst(mediumhdlisdata[2][[1]]$data,
                                        p.thres = 1e-05,diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(mediumhdlisdata[2][[1]]$data, p.thres = 1e-8, 
                       marker.data = mediumhdlisdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("mediumhdl_is_profile_10-8.pdf")

# analysis of large HDL concentration on IS with kettunen, davyson, malik
sel.file.hdl.large <- "largeHDL-kettunen16new.csv"
exp.file.hdl.large <- "largeHDL-davyson23new.csv"
out.file.is <- "IS-malik18new.csv"
largehdlisdata <- preprocessforunivariable(sel.file.hdl.large,
                                           exp.file.hdl.large,out.file.is,
                                           "large_HDL","IS")
univrunlargeHDLIS <-rununivariableMRs(largehdlisdata[1][[1]],
                                      largehdlisdata[2][[1]],c(1e-8,1e-5,1e-3),
                                      "large_HDL","IS")
ggsave("largeHDL_is_results.pdf")
#seek to find diagnostic plots for RAPS
largeHDLisforplots <- grappleRobustEst(largehdlisdata[2][[1]]$data,
                                       p.thres = 1e-05,diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(largehdlisdata[2][[1]]$data, p.thres = 1e-8, 
                       marker.data = largehdlisdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("largehdl_is_profile_10-8.pdf")

#APOB on IS
sel.file.apob2 <- "APOB-kettunen16new.csv"
exp.file.apob <- "APOB-davyson23new.csv"
out.file.is <- "IS-malik18new.csv"
apoisdata <- preprocessforunivariable(sel.file.apob2,exp.file.apob,out.file.is,
                                      "APOB","IS")
univrunAPOBIS <- rununivariableMRs(apoisdata[1][[1]],apoisdata[2][[1]],
                                   c(1e-8,1e-5,1e-3),"APOB","IS")
#seek to find diagnostic plots for RAPS
apoisforplots <- grappleRobustEst(apoisdata[2][[1]]$data,p.thres = 1e-05,
                                  diagnosis = TRUE)
#grapple methodology: seek modes in profile likelihood
diagnosis <- findModes(apoisdata[2][[1]]$data, p.thres = 1e-5, 
                       marker.data = apoisdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers

# APOA1 on IS 
sel.file.apoa1 <- "APOA1-kettunen16new.csv"
exp.file.apoa1 <- "APOA1-davyson23new.csv"
out.file.is <- "IS-malik18new.csv"
apoa1isdata <- preprocessforunivariable(sel.file.apoa1,exp.file.apoa1,
                                        out.file.is,"APOA1","IS")
univrunAPOA1IS <-rununivariableMRs(apoa1isdata[1][[1]],apoa1isdata[2][[1]],
                                   c(1e-8,1e-5,1e-3),"APOA1","IS")
ggsave("apoa1-is-univariateresultschart.pdf")
#seek to find diagnostic plots for RAPS
apoa1isforplots <- grappleRobustEst(apoa1isdata[2][[1]]$data,p.thres = 1e-05,
                                    diagnosis = TRUE)
# run grapple and look for modes
diagnosis <- findModes(apoa1isdata[2][[1]]$data, p.thres = 1e-8, 
                       marker.data = apoa1isdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers
ggsave("apoa1_profile_10-8.pdf")

#LDL on IS
sel.file.ldl <- "datasets_grapple/LDL-gera18.csv"
exp.file.ldl <- "LDL-graham21new.csv"
out.file.is <- "IS-malik18new.csv"
ldlisdata <- preprocessforunivariable(sel.file.ldl,exp.file.ldl,out.file.is,
                                      "LDL","IS")
univrunLDLIS <- rununivariableMRs(ldlisdata[1][[1]],ldlisdata[2][[1]],
                                  c(1e-8,1e-5,1e-3),"LDL","IS")
#seek to find diagnostic plots for RAPS
ldlisforplots <- grappleRobustEst(ldlisdata[2][[1]]$data,p.thres = 1e-05,
                                  diagnosis = TRUE)
#grapple methodology: seek modes in profile likelihood
diagnosis <- findModes(ldlisdata[2][[1]]$data, p.thres = 1e-6, 
                       marker.data = ldlisdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers


#TG on IS
sel.file.tg <- "datasets_grapple/TG-gera18.csv"
exp.file.tg <- "TG-graham21new.csv"
out.file.is <- "IS-malik18new.csv"
tgisdata <- preprocessforunivariable(sel.file.tg,exp.file.tg,out.file.is,
                                     "TG","IS")
univrunTGIS <- rununivariableMRs(tgisdata[1][[1]],tgisdata[2][[1]],
                                 c(1e-8,1e-5,1e-3),"TG","IS")
#seek to find diagnostic plots for RAPS
tgisforplots <- grappleRobustEst(tgisdata[2][[1]]$data,p.thres = 1e-05,
                                 diagnosis = TRUE)
#grapple methodology: seek modes in profile likelihood
diagnosis <- findModes(tgisdata[2][[1]]$data, p.thres = 1e-3, 
                       marker.data = tgisdata[2][[1]]$marker.data)
diagnosis$p
diagnosis$markers

#save the univariable results
write.csv(univrunAPOBIS,"univariable_apob_IS.csv",row.names = FALSE)
write.csv(univrunAPOBCAD,"univariable_apob_CAD.csv",row.names = FALSE)
write.csv(univrunAPOA1IS,"univariable_apob_IS.csv",row.names = FALSE)
write.csv(univrunAPOA1CAD,"univariable_apob_CAD.csv",row.names = FALSE)
write.csv(univrunHDLIS,"univariable_HDL_IS.csv",row.names = FALSE)
write.csv(univrunHDLCAD,"univariable_HDL_CAD.csv",row.names = FALSE)
write.csv(univrunLDLCAD,"univariable_LDL_CAD.csv",row.names = FALSE)
write.csv(univrunLDLIS,"univariable_LDL_IS.csv",row.names = FALSE)
write.csv(univrunTGCAD,"univariable_TG_CAD.csv",row.names = FALSE)
write.csv(univrunTGIS,"univariable_TG_IS.csv",row.names = FALSE)
write.csv(univrunsmallHDLCAD,"univariable_smallHDL_CAD.csv",row.names = FALSE)
write.csv(univrunsmallHDLIS,"univariable_smallHDL_IS.csv",row.names = FALSE)
write.csv(univrunmediumHDLCAD,"univariable_mediumHDL_CAD.csv",row.names = FALSE)
write.csv(univrunmediumHDLIS,"univariable_mediumHDL_IS.csv",row.names = FALSE)
write.csv(univrunlargeHDLCAD,"univariable_largeHDL_CAD.csv",row.names = FALSE)
write.csv(univrunlargeHDLIS,"univariable_largeHDL_IS.csv",row.names = FALSE)



# multivariable IS
# LDL and ApoB on IS
sel.file.apob.ldl <-c("APOB-kettunen16new.csv", 
                      "datasets_grapple/LDL-gera18.csv")
exp.file.apob.ldl <- c("APOB-davyson23new.csv", "LDL-graham21new.csv")
out.file.is <- "IS-malik18new.csv"
ldlapobisdata <- getInput(sel.file.apob.ldl, exp.file.apob.ldl, out.file.is, 
                          plink_refdat, max.p.thres = 0.01,
                          plink_exe = '~/Downloads/plink_mac_20241022/plink')
multirunapoldlis <- grappleRobustEst(ldlapobisdata$data, p.thres = 1e-4)
multirunapoldlis2 <- grappleRobustEst(ldlapobisdata$data,p.thres = 1e-6)
multirunapoldlis3 <- grappleRobustEst(ldlapobisdata$data,p.thres = 1e-3)
dattoplot2 <- data.frame(b=c(multirunapoldlis$beta.hat,
                             multirunapoldlis2$beta.hat,
                             multirunapoldlis3$beta.hat),
                         se=c(sqrt(diag(multirunapoldlis$beta.var)),
                              sqrt(diag(multirunapoldlis2$beta.var)),
                              sqrt(diag(multirunapoldlis3$beta.var))),
                         pval=c(multirunapoldlis$beta.p.value,
                                multirunapoldlis2$beta.p.value,
                                multirunapoldlis3$beta.p.value),
                         exposure=rep(c("ApoB","LDL"),3),
                         method = rep("Multivariable GRAPPLE",6),
                         pval.sel = c(rep(10^(-4),2),rep(10^(-6),2),
                                      rep(10^(-3),2)),
                         pval.adjusted= rep(0,6))
dattoplot2$pval.adjusted <- p.adjust(dattoplot2$pval,"bonferroni")
library(ggplot2)
print(ggplot(dattoplot2) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot2, "apob-ldl-is-results.csv", row.names=FALSE)
ggsave("apob-ldl-is-graph.pdf")

# LDL,HDL and TG on IS
sel.file.many <- c("datasets_grapple/LDL-gera18.csv",
                   "datasets_grapple/HDL-gera18.csv",
                   "datasets_grapple/TG-gera18.csv")
exp.file.many <- c("LDL-graham21new.csv","HDL-graham21new.csv",
                   "TG-graham21new.csv")
out.file.is <- "IS-malik18new.csv"
multiisdata <- getInput(sel.file.many, exp.file.many, out.file.is, 
                        plink_refdat, max.p.thres = 0.01,
                        plink_exe = '~/Downloads/plink_mac_20241022/plink')
multirunis <- grappleRobustEst(multiisdata$data, p.thres = 1e-4)
multirunis2 <-grappleRobustEst(multiisdata$data, p.thres = 1e-2)
multirunis3 <- grappleRobustEst(multiisdata$data, p.thres = 1e-6)

dattoplot3 <- data.frame(b=c(multirunis$beta.hat,multirunis2$beta.hat,
                             multirunis3$beta.hat),
                         se=c(sqrt(diag(multirunis$beta.var)),
                              sqrt(diag(multirunis2$beta.var)),
                              sqrt(diag(multirunis3$beta.var))),
                         pval=c(multirunis$beta.p.value,
                                multirunis2$beta.p.value,
                                multirunis3$beta.p.value),
                         exposure=rep(c("LDL","HDL","TG"),3),
                         method = rep("Multivariable GRAPPLE",9),
                         pval.sel = c(rep(10^(-4),3),rep(10^(-2),3),
                                      rep(10^(-6),3)),
                         pval.adjusted= rep(0,9))
dattoplot3$pval.adjusted <- p.adjust(dattoplot3$pval,"bonferroni")
print(ggplot(dattoplot3) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot3,"ldl-hdl-tg-is-results.csv", row.names=FALSE)


#multivariate small HDL with LDL, HDL and TG on IS
sel.file.smallhdl.hdl.ldl.tg <- c("smallHDL-kettunen16new.csv",
                                  "datasets_grapple/LDL-gera18.csv",
                                  "datasets_grapple/HDL-gera18.csv",
                                  "datasets_grapple/TG-gera18.csv")
exp.file.smallhdl.hdl.ldl.tg <- c("smallHDL-davyson23new.csv",
                                  "LDL-graham21new.csv","HDL-graham21new.csv",
                                  "TG-graham21new.csv")
out.file.is <- "IS-malik18new.csv"
smallhdl_multi_isdata <- getInput(sel.file.smallhdl.hdl.ldl.tg, 
                                  exp.file.smallhdl.hdl.ldl.tg, out.file.is, 
                                  plink_refdat, max.p.thres = 0.01,
                                  plink_exe =
                                    '~/Downloads/plink_mac_20241022/plink')
smallhdl_multi_is <- grappleRobustEst(smallhdl_multi_isdata$data, 
                                      p.thres = 1e-4)
smallhdl_multi_is2 <-grappleRobustEst(smallhdl_multi_isdata$data, 
                                      p.thres = 1e-2)
smallhdl_multi_is3 <- grappleRobustEst(smallhdl_multi_isdata$data, 
                                       p.thres = 1e-6)

# plot results of multivariate small HDL LDL HDL TG on IS
dattoplot4 <- data.frame(b=c(smallhdl_multi_is$beta.hat,
                             smallhdl_multi_is2$beta.hat,
                             smallhdl_multi_is3$beta.hat),
                         se=c(sqrt(diag(smallhdl_multi_is$beta.var)),
                              sqrt(diag(smallhdl_multi_is2$beta.var)),
                              sqrt(diag(smallhdl_multi_is3$beta.var))),
                         pval=c(smallhdl_multi_is$beta.p.value,
                                smallhdl_multi_is2$beta.p.value,
                                smallhdl_multi_is3$beta.p.value),
                         exposure=rep(c("smallHDL","LDL","HDL","TG"),3),
                         method = rep("Multivariable GRAPPLE",12),
                         pval.sel = c(rep(10^(-4),4),rep(10^(-2),4),
                                      rep(10^(-6),4)),
                         pval.adjusted= rep(0,12))
dattoplot4$pval.adjusted <- p.adjust(dattoplot4$pval,"bonferroni")
print(ggplot(dattoplot4) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot4,"smallhdl-ldl-hdl-tg-is-multivar-results.csv")
ggsave("smallhdl-ldl-hdl-tg-is-multivar-plot.pdf")

#multivariate medium HDL with LDL, HDL and TG
sel.file.mediumhdl.hdl.ldl.tg <- c("mediumHDL-kettunen16new.csv",
                                   "datasets_grapple/LDL-gera18.csv",
                                   "datasets_grapple/HDL-gera18.csv",
                                   "datasets_grapple/TG-gera18.csv")
exp.file.mediumhdl.hdl.ldl.tg <- c("mediumHDL-davyson23new.csv",
                                   "LDL-graham21new.csv","HDL-graham21new.csv",
                                   "TG-graham21new.csv")
out.file.is <- "IS-malik18new.csv"
mediumhdl_multi_isdata <- getInput(sel.file.mediumhdl.hdl.ldl.tg, 
                                   exp.file.mediumhdl.hdl.ldl.tg, out.file.is, 
                                   plink_refdat, max.p.thres = 0.01,
                                   plink_exe = 
                                     '~/Downloads/plink_mac_20241022/plink')
mediumhdl_multi_is <- grappleRobustEst(mediumhdl_multi_isdata$data, 
                                       p.thres = 1e-4)
mediumhdl_multi_is2 <-grappleRobustEst(mediumhdl_multi_isdata$data, 
                                       p.thres = 1e-2)
mediumhdl_multi_is3 <- grappleRobustEst(mediumhdl_multi_isdata$data, 
                                        p.thres = 1e-6)

# plot results of multivariate medium HDL LDL HDL TG on IS
dattoplot5 <- data.frame(b=c(mediumhdl_multi_is$beta.hat,
                             mediumhdl_multi_is2$beta.hat,
                             mediumhdl_multi_is3$beta.hat),
                         se=c(sqrt(diag(mediumhdl_multi_is$beta.var)),
                              sqrt(diag(mediumhdl_multi_is2$beta.var)),
                              sqrt(diag(mediumhdl_multi_is3$beta.var))),
                         pval=c(mediumhdl_multi_is$beta.p.value,
                                mediumhdl_multi_is2$beta.p.value,
                                mediumhdl_multi_is3$beta.p.value),
                         exposure=rep(c("mediumHDL","LDL","HDL","TG"),3),
                         method = rep("Multivariable GRAPPLE",12),
                         pval.sel = c(rep(10^(-4),4),rep(10^(-2),4),
                                      rep(10^(-6),4)),
                         pval.adjusted= rep(0,12))
dattoplot5$pval.adjusted <- p.adjust(dattoplot5$pval,"bonferroni")
print(ggplot(dattoplot5) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot5,"mediumhdl-ldl-hdl-tg-is-multivar-results.csv")
ggsave("mediumhdl-ldl-hdl-tg-is-multivar-plot.pdf")


# HDL and ApoA1 on IS
sel.file.apoa1.hdl <-c("APOA1-kettunen16new.csv", 
                       "datasets_grapple/HDL-gera18.csv")
exp.file.apoa1.hdl <- c("APOA1-davyson23new.csv", "HDL-graham21new.csv")
out.file.is <- "IS-malik18new.csv"
hdlapoa1isdata <- getInput(sel.file.apoa1.hdl, exp.file.apoa1.hdl, out.file.is, 
                           plink_refdat, max.p.thres = 0.01,
                           plink_exe = '~/Downloads/plink_mac_20241022/plink')
multirunapohdlis <- grappleRobustEst(hdlapoa1isdata$data, p.thres = 1e-4)
multirunapohdlis2 <- grappleRobustEst(hdlapoa1isdata$data,p.thres = 1e-6)
multirunapohdlis3 <- grappleRobustEst(hdlapoa1isdata$data,p.thres = 1e-2)
multirunapohdlis4 <- mutlirunapohdlis3
dattoplot2 <- data.frame(b=c(multirunapohdlis$beta.hat,
                             multirunapohdlis2$beta.hat,
                             multirunapohdlis3$beta.hat),
                         se=c(sqrt(diag(multirunapohdlis$beta.var)),
                              sqrt(diag(multirunapohdlis2$beta.var)),
                              sqrt(diag(multirunapohdlis3$beta.var))),
                         pval=c(multirunapohdlis$beta.p.value,
                                multirunapohdlis2$beta.p.value,
                                multirunapohdlis3$beta.p.value),
                         exposure=rep(c("ApoA1","HDL"),3),
                         method = rep("Multivariable GRAPPLE",6),
                         pval.sel = c(rep(10^(-4),2),rep(10^(-6),2),
                                      rep(10^(-2),2)),
                         pval.adjusted= rep(0,6))
dattoplot2$pval.adjusted <- p.adjust(dattoplot2$pval,"bonferroni")
library(ggplot2)
print(ggplot(dattoplot2) + aes(x = method, y = b, ymin = b - 1.96 * se, 
                               ymax = b + 1.96 * se,
                               col = pval.adjusted < 0.05) +
        geom_point() + geom_errorbar() + geom_hline(yintercept = 0, 
                                                    linetype = "dashed") +
        facet_grid(exposure ~ pval.sel) + coord_flip() +
        theme_light(base_size = 18) + theme(legend.position = "bottom"))
write.csv(dattoplot2, "apoa1-hdl-is-results.csv", row.names=FALSE)
ggsave("apoa1-hdl-is-multivar-plot.pdf")
