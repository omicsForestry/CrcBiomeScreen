NormalizeData <- function(CrcBiomeScreenObject = NULL,
                            method = NULL,
                            TaskName = NULL){
    Data <- CrcBiomeScreenObject$GenusLevelData
    if(method == "TSS"){
                # Transforming into relative abundance
                Data <-  t(normalize(t(Data), method = "TSS"))
        }else if(method == "GMPR"){
            size.factor <- GMPR(t(Data))
            size.factor[is.na(size.factor)] <- mean(size.factor,na.rm = TRUE)
            Data <- Data / size.factor 
        }

    # Add the normalization method to the Data
    attr(Data, "NormalizationMethod") <- method
    attr(Data, "TaskName") <- TaskName
    attr(Data, "Timestamp") <- Sys.time()
    CrcBiomeScreenObject$NormalizedData <- Data
    # Save the Data
    saveRDS(CrcBiomeScreenObject,paste0("CrcBiomeScreenObject_",TaskName,".rds"))
    return(CrcBiomeScreenObject)
}