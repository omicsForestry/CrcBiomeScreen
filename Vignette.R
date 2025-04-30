setwd("/home/CRCscreening/CRCscreening-Workflow/")
setwd("/mnt/scratch/ngzh5554/CRCscreening-Workflow")
# Load required libraries
source("Environment.R")
source("Workflow_function.R")
# Load the data_relative abundance
# toydata <- readRDS("toy_data.rds")
toydata <- curatedMetagenomicData(
            "ZellerG_2014.relative_abundance"
            , dryrun = FALSE, rownames = "short")  
# otu.data <- toydata[[1]]@assays@data@listData$relative_abundance
# rownames(otu.data) <- toydata[[1]]@rowLinks$nodeLab

CrcBiomeScreenObject <- CreateCrcBiomeScreenObject(RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
                                                  TaxaData = toydata[[1]]@rowLinks$nodeLab,
                                                  SampleData = toydata[[1]]@colData)
CrcBiomeScreenObject <- SplitTaxas(CrcBiomeScreenObject)
CrcBiomeScreenObject <- KeepGenusLevel(CrcBiomeScreenObject)
CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject, method = "GMPR")
CrcBiomeScreenObject <- SplitDataSet(CrcBiomeScreenObject, label = c("control","CRC"), partition = 0.7)







taxa.data <- tibble(variable = toydata[[1]]@rowLinks$nodeLab)
taxa.data %>%
  separate(variable, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep = '\\|', fill = "right") %>%
  mutate(across(everything(), ~ ifelse(. == "", NA, sub("^[a-z]__", "", .)))) -> taxa.temp
taxa.temp <- as.matrix(taxa.temp)
rownames(taxa.temp) <- toydata[[1]]@rowLinks$nodeLab

otu.sc <- otu_table(otu.data, taxa_are_rows = TRUE)
taxa.sc <- tax_table(taxa.temp)
sam.sc <- sample_data(as.data.frame(toydata[[1]]@colData))

toydata_phyloseq <- phyloseq(otu.sc, taxa.sc, sam.sc)

# Save the data in one object
screening_results$OriginalData <- toydata_phyloseq
# rownames(screening_results$OriginalData@otu_table) <- screening_results$OriginalData@tax_table@.Data[, 6]
screening_results$originalAbundance <- as.data.frame(t(screening_results$OriginalData@otu_table))

# Preprocess the data
## Select the taxas at Genus level
## Change to the absolute abundance
screening_results$preprocessing <- screening_results$originalAbundance
Genus_name <- screening_results$OriginalData@tax_table@.Data[, 6]
duplicated_colnames <- Genus_name[duplicated(Genus_name)]

for (colname in unique(duplicated_colnames)) {

  cols_to_merge <- screening_results$preprocessing[, which(Genus_name == colname)]  

  cols_to_merge <- as.data.frame(lapply(cols_to_merge, function(x) as.numeric(as.character(x))))
  merged_col <- rowSums(cols_to_merge, na.rm = TRUE)

  screening_results$preprocessing <- screening_results$preprocessing[, -which(Genus_name == colname)]

  screening_results$preprocessing[[colname]] <- merged_col
  
  Genus_name <- Genus_name[-which(Genus_name == colname)]
  Genus_name <- c(Genus_name, colname)
  
}
colnames(screening_results$preprocessing) <- Genus_name


screening_results$preprocessing <- screening_results$preprocessing %>%
  mutate(across(seq_len(dim(screening_results$preprocessing)[2]), ~ (. * screening_results$OriginalData@sam_data$number_reads / 100)))

# Preprocess the data
## Normalize the data
screening_results$NormalizedData <- ChooseNormalize(screening_results$preprocessing, method = "GMPR",taxa_col = seq_len(dim(screening_results$preprocessing)[2]),TaskName = "ToyData")
# Perform the analysis
screening_results$ModelData <- SplitDataSet(screening_results$NormalizedData,label = "StudyCondition", partition = 0.7)

## RF model
best.params <- RF_Modeling(ModelData = toy_data_ModelData,
                        label = "StudyCondition",
                        k.rf = 10,
                        taxa_col = c(3:150),
                        TaskName = "ToyData_RF",
                        TrueLabel = "Cancer")

TestModel(ModelData = toy_data_ModelData,
                      label = "StudyCondition",
                      taxa_col = c(3:150),
                      best.params = best.params,
                      TaskName = "ToyData_RF_Test",
    ModelData <- list(train,test)
                      TrueLabel = "Cancer")

## XGBoost model
XGBoost_Modeling(ModelData = toy_data_ModelData,
                label = "StudyCondition",
                taxa_col = c(3:150),
                TaskName = "ToyData_XGBoost",
                TrueLabel = "Cancer")








