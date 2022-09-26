setwd("C:/Users/ac361/desktop/QBIO_490/qbio_490_avinash/analysis_data") #set directory at beginning to premade folder
clinical <- read.csv("brca_clinical_data.csv")
library(BiocManager)
library(TCGAbiolinks)

clin_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")

clinical.drug <- GDCprepare_clinic(query = clin_query,clinical.info = "drug")

clinical.rad <- GDCprepare_clinic(query = clin_query,clinical.info = "radiation") #prepare the drug and radiation data which has already been queried

surgery <- clinical[,c(1,34)]
surgery <- surgery[surgery$breast_carcinoma_surgical_procedure_name != "",]
rad.region <- clinical.rad[,c(1,11)]
rad.region <- rad.region[rad.region$anatomic_treatment_site != "",] #helps to remove whitespace, as the null command does not pick these entries up

region.v.surgery <- rad.region
region.v.surgery$breast_carcinoma_surgical_procedure_name <- NA

for(i in 1:length(region.v.surgery$anatomic_treatment_site)) {
  surgery.entry <- surgery[which(surgery$bcr_patient_barcode == region.v.surgery$bcr_patient_barcode[i]),]
  if (nrow(surgery.entry) != 0) {
    region.v.surgery[i,3] <- surgery.entry[1,2] #combines the region and surgery info thru this for loop
  }
}

region.v.surgery <- region.v.surgery[!is.na(region.v.surgery$breast_carcinoma_surgical_procedure_name),]
region.v.surgery$anatomic_treatment_site <- as.character(region.v.surgery$anatomic_treatment_site) #wanted to convert from integer to character to get the name of the surgery type

rad.site <- c(rep("Primary Tumor Field", 4),rep("Regional site",4), rep("Local Recurrence", 4), rep("Distant Recurrence", 4), rep("Distant site", 4))
surgery.type <- rep(c("Lumpectomy", "Simple Mastectomy", "Modified Radical Mastectomy", "Other"), 5)

test <- data.frame(rad.site, surgery.type, count= NA)


for(i in 1:length(test$rad.site)){
  all.entries <- region.v.surgery[which((region.v.surgery$anatomic_treatment_site == test[i,1]) & (region.v.surgery$breast_carcinoma_surgical_procedure_name == test[i,2])),]
  test[i,3] <- nrow(all.entries) #get the total counts for each of the pairs/entry for surgical/radiation
}

library(ggplot2)

jpeg("region.v.surgery.jpg")
ggplot(test, aes(fill=surgery.type, y=count, x=rad.site)) + geom_bar(position= "dodge", stat="identity")
dev.off() #save the plot to computer


#KM PLOTS
library(survival)
library(survminer)

clinical$survival_time <- ifelse(!is.na(clinical$days_to_death), clinical$days_to_death, clinical$days_to_last_followup)
clinical$death_event <- ifelse(as.character(clinical$vital_status) == "Alive", FALSE, TRUE)

#initialize survival object by incorporating all the surgical information
surv_object_surgery <- Surv(time = clinical$survival_time, event = clinical$death_event)

#this creates the fit object
surgery_fit <- surv_fit(surv_object_surgery ~ breast_carcinoma_surgical_procedure_name, data=clinical)

#the ggtheme and legend arguments affect formatting.
?ggsurvplot
survplot_surgery = ggsurvplot(surgery_fit, pval = TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

#formatting minutiae
KM_plot_surgery = survplot_surgery$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title=element_text(size=5), legend.text = element_text(size=5))





clinical.rad$survival_time <- NA
clinical.rad$death_event <- NA

for (i in 1:length(clinical.rad$bcr_patient_barcode)) {
  exist <- clinical[which(clinical$bcr_patient_barcode == clinical.rad$bcr_patient_barcode[i]),]
  if (nrow(exist) != 0) {
    clinical.rad[i, 21] <- exist[1, 116]
    clinical.rad[i, 22] <- exist[1, 117]
  } #loop through element via the for loop, as only some patients are in the radiation database; so, we have to get their info from the clinical database and insert into radiation database
}


  surv_object_region <- Surv(time = clinical.rad$survival_time, event = clinical.rad$death_event)

#Create a fit object
region_fit <- surv_fit(surv_object_region ~ anatomic_treatment_site, data=clinical.rad)

#the ggtheme and legend arguments are for formatting.
?ggsurvplot
survplot_region = ggsurvplot(region_fit, pval = TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")


KM_plot_region = survplot_region$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title=element_text(size=5), legend.text = element_text(size=5))


#PRINT SURGERY GRAPH
jpeg("surgery_km.jpg")
KM_plot_surgery
dev.off()

#PRINT REGION GRAPH
jpeg("region_km.jpg")
KM_plot_region
dev.off()


write.csv(clinical, "brac_clinical_data.csv")
write.csv(clinical.drug, "brac_clinical_drug.csv")
write.csv(clinical.rad, "brac_clinical_rad.csv")






