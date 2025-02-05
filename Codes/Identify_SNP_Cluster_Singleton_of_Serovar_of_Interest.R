##-----------------------------------------------------------
## Project: USDA - Salmonella virulence
## Serovar: S. Typhimurium
## Script purpose: identification of candidate S. Typhimurium SNP clusters and singletons for downstream analysis.

## Start Date:  August 22, 2023
##-----------------------------------------------------------
## Notes: 
##
##-----------------------------------------------------------
## Load packages.
##-----------------------------------------------------------
rm(list=ls())
# Install essential packages if haven't already.
required_packages <- c('plyr', 'tidyverse', 'readr', "openxlsx", "readxl")
for (p in required_packages) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p, dep = TRUE)
  }
}
# Load the packages.
library(plyr);library(tidyverse);library(readr);library(openxlsx);library(readxl)
##-----------------------------------------------------------
## Load the merged metadata.
## Metadata below can be find at https://github.com/by967/USDA_Salmonella_Phylogeny_Project/tree/e46ebb5adae4c7243419cc62f52c811685f96f94/Metadata
## Change your file name accordingly if you are using a different input
Metadata_Merged <- read_csv("Metadata_Merged_BY080823.csv")
##-----------------------------------------------------------
## Select all isolates with computed_types indicating "Typhimurium"
## or the corresponding antigenic formula
## I 1(underscore),4,[5],12:i:1,2
ComputedType_Unique <- unique(Metadata_Merged$computed_types)
grep("Typhimurium", ComputedType_Unique, value = TRUE, ignore.case = TRUE)
grep("4:i:1,2", ComputedType_Unique, value = TRUE, ignore.case = TRUE)
# 4:i:1,2

# Typhimurium_Terms <- c("Typhimurium", "4:i:1,2")
# Extract all putative Typhimurium isolates.
Metadata_Merged_Typhimurium <- Metadata_Merged %>% filter(grepl("Typhimurium", computed_types) | grepl("4:i:1,2", computed_types))


##-----------------------------------------------------------
## Select all isolates in SNP clusters that have one or more isolates with computed_types "Typhimurium" and/or "4:i:1,2".
CandidateSNPClusters <- unique(Metadata_Merged_Typhimurium$PDS_acc)[!is.null(unique(Metadata_Merged_Typhimurium$PDS_acc))]
CandidateSNPCluster_Isolate <- Metadata_Merged %>% filter(PDS_acc %in% CandidateSNPClusters)
# 139795 rows.

## Retrieve all singleton isolates.
CandidateSingleton_Isolate <- Metadata_Merged_Typhimurium %>% filter(is.na(PDS_acc))
# 5107 rows.

## Merge isolates in SNP clusters and singletons to get all candidate Typhimurium isolates.
CandidateTyphimurium_Isolate <- rbind(CandidateSNPCluster_Isolate, CandidateSingleton_Isolate)
CandidateTyphimurium_Isolate$PDS_acc <- CandidateTyphimurium_Isolate$PDS_acc %>% replace_na("Not_Assigned")



##-----------------------------------------------------------
## Check the distribution of isolates across "computed_types" for each candidate SNP cluster.
CandidateTyphimurium_Isolate_Summary <- ddply(CandidateTyphimurium_Isolate, c("PDS_acc", "computed_types"), summarise, N = sum(!is.na(target_acc)))

## For the isolates with an empty cell for computed_types, fill in "Unknown" for now.
CandidateTyphimurium_Isolate_Summary$computed_types <- CandidateTyphimurium_Isolate_Summary$computed_types %>% replace_na("Unknown")
##-----------------------------------------------------------
## Get high-confident and questionable SNP clusters.
HC_SNPCluster <- data.frame()
Questionable_SNPCluster <- data.frame()

for (cluster in unique(CandidateTyphimurium_Isolate_Summary$PDS_acc)) {
  target_content = CandidateTyphimurium_Isolate_Summary %>% filter(PDS_acc == cluster)
  
  target_content_Typhimurium = sum(target_content$N[which(grepl("Typhimurium", target_content$computed_types, ignore.case = TRUE) | grepl("4:i:1,2", target_content$computed_types, ignore.case = TRUE))])
  target_content_Typhimurium = ifelse(length(target_content_Typhimurium)==0, 0, target_content_Typhimurium)
  target_content_Others = sum(target_content$N) - target_content_Typhimurium
  
  if (target_content_Typhimurium > target_content_Others) {
    HC_SNPCluster <- rbind(HC_SNPCluster, target_content)
  } else {
    Questionable_SNPCluster <- rbind(Questionable_SNPCluster, target_content)
  }
}

# Write out two files, one for high-confident and one for questionable SNP clusters.
write.xlsx(Questionable_SNPCluster, file = "Typhimurium_SNPCluster_Questionable_BY08222023.xlsx")
write.xlsx(HC_SNPCluster, file = "Typhimurium_SNPCluster_HighConfident_BY08222023.xlsx")
# Manually edit the spreadsheets.
# All high-confident SNP clusters should be classified as "Maybe".
# Questionable SNP clusters can be classified as "No" or "Maybe", 
# depending on if their proportion of Typhimurium isolates is less than 
# or more than a threshold (e.g., 10%).
##-----------------------------------------------------------
## Get all candidate SNP cluster IDs.

# Read in SNP cluster files.
Typhimurium_SNPCluster_HighConfident <- read_excel("Typhimurium_SNPCluster_HighConfident_BY08222023.xlsx")
Typhimurium_SNPCluster_Questionable <- read_excel("Typhimurium_SNPCluster_Questionable_BY08222023.xlsx")

# Combine the files while removing SNP clusters classified as "No".
Typhimurium_CandidateSNPCluster <- rbind(Typhimurium_SNPCluster_HighConfident, subset(Typhimurium_SNPCluster_Questionable, Status == "Maybe"))
Typhimurium_CandidateSNPCluster_Name <- unique(Typhimurium_CandidateSNPCluster$PDS_acc)

##-----------------------------------------------------------
## Select one representative isolate for each candidate SNP cluster.
RepIso_Df <- expand.grid("SNP_Cluster" = Typhimurium_CandidateSNPCluster_Name[1:length(Typhimurium_CandidateSNPCluster_Name)])
RepIso_Df$Representative <- NA
RepIso_Df$Contig_Number <- NA
for (cluster in RepIso_Df$SNP_Cluster) {
  target_content = Metadata_Merged %>% filter(PDS_acc == cluster)
  target_content_ordered = target_content[order(target_content$asm_stats_n_contig),]
  if (!is.na(target_content_ordered$asm_acc[1])) {
    RepIso_Df$Representative[which(RepIso_Df$SNP_Cluster==cluster)] = target_content_ordered$asm_acc[1]
    RepIso_Df$Contig_Number[which(RepIso_Df$SNP_Cluster==cluster)] = target_content_ordered$asm_stats_n_contig[1]
  } else {
    RepIso_Df$Representative[which(RepIso_Df$SNP_Cluster==cluster)] = target_content_ordered$Run[1]
    RepIso_Df$Contig_Number[which(RepIso_Df$SNP_Cluster==cluster)] = target_content_ordered$asm_stats_n_contig[1]
  }
  
}

write.xlsx(RepIso_Df, file = "Typhimurium_CandidateSNPCluster_RepresentativeIsolate_BY08282023.xlsx", overwrite = TRUE)
# A total of 3027 candidate SNP clusters have been identified in this case. 
# Here I don't need to remove any SNP clusters as the highest 
# number of contigs for one representative isolate is 450 which is lower than our cutoff (= or > 500 contigs)
# Remember to save your file!

##-----------------------------------------------------------
## Get all Typhimurium singletons.
Singleton_Df <- expand.grid("Singleton" = CandidateSingleton_Isolate$target_acc)
Singleton_Df$GB_Acc_Run <- NA
Singleton_Df$Contig_Number <- NA
for (singleton in Singleton_Df$Singleton) {
  target_content = CandidateSingleton_Isolate %>% filter(target_acc == singleton)
  if (!is.na(target_content$asm_acc[1])) {
    Singleton_Df$Contig_Number[which(Singleton_Df$Singleton==singleton)] = target_content$asm_stats_n_contig[1]
    Singleton_Df$GB_Acc_Run[which(Singleton_Df$Singleton==singleton)] = ifelse(target_content$asm_stats_n_contig[1]<500, target_content$asm_acc[1], target_content$Run[1])
  } else {
    Singleton_Df$Contig_Number[which(Singleton_Df$Singleton==singleton)] = target_content$asm_stats_n_contig[1]
    Singleton_Df$GB_Acc_Run[which(Singleton_Df$Singleton==singleton)] = target_content$Run[1]
  }
  
}

write.xlsx(Singleton_Df, file = "Typhimurium_CandidateSingleton_BY08282023.xlsx", overwrite = TRUE)

## Manually inspect the table to fill in empty cells.
# Remove any isolates if they have more than 500 contigs (including 500). 
# In this case, 61 isolates have been removed due to low sequence quality.  
# A total of 5046 singletons are left
# Remember to save your file! 

##-----------------------------------------------------------
## Here in excel, I changed the column name "GB_Acc_Run" to "Representative".
## I also changed the column name "SNP_Cluster" and "Singleton" to "SNPCluster_Singleton".
## This is to facilitate the combination of the two files.
PutativeTyphimurium_SNPCluster <- read_excel("Typhimurium_CandidateSNPCluster_RepresentativeIsolate_BY08282023.xlsx")
PutativeTyphimurium_Singleton <- read_excel("Typhimurium_CandidateSingleton_BY08282023.xlsx")
PutativeTyphimurium_SNPClusterSingleton <- rbind(PutativeTyphimurium_SNPCluster, PutativeTyphimurium_Singleton)
write.xlsx(PutativeTyphimurium_SNPClusterSingleton, file = "Typhimurium_CandidateSNPClusterSingleton_BY08282023.xlsx", overwrite = TRUE)
##-----------------------------------------------------------
## In the spreadsheet, retrieve all SRR/ERR numbers and save them into a .txt file;
## and retrieve all GCA numbers and save them into another .txt file.
## Use the GCA numbers to download assemblies from NCBI PD, and use the SRR/ERR
## numbers to run SKESA for de novo assembling.
##-----------------------------------------------------------




