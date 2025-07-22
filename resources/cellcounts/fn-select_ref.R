validate_and_select_reference <- function(tissue, methylation_array, age) {
  # Define valid combinations
  blood_tissues <- c("blood", "pbmc", "pbl", "cord")
  saliva_tissues <- c("saliva", "buccal")
  arrays <- c("epic", "epic2", "450k")
  age_arrange <- c(">18", "<18", "all")
  
  # Validate inputs
  if (!tissue %in% c(blood_tissues, saliva_tissues)) {
    stop("Unknown tissue type. Please specify it again in config file.")
  }
  if (!methylation_array %in% arrays) {
    stop("Unknown methylation array type. Please specify it again in config file.")
  }
  if (!age %in% age_arrange) {
    stop("Unknown age group. Please specify it again in config file.")
  }
  
  # Select appropriate reference matrix based on tissue and array
  if (tissue %in% blood_tissues) {
    if (methylation_array %in% c("epic", "epic2")) {
      if (age == ">18") {
        return(list(salas = cent12CT.m, unilife = centUniLIFE.m))
      } else if (age == "<18") {
        return(list(unilife = centUniLIFE.m))
      } else if (age == "all") {
        return(list(unilife = centUniLIFE.m))
      }
    } else if (methylation_array == "450k") {
      if (age == ">18") {
        return(list(salas = cent12CT450k.m, unilife = centUniLIFE.m))
      } else if (age == "<18") {
        return(list(unilife = centUniLIFE.m))
      } else if (age == "all") {
        return(list(unilife = centUniLIFE.m))
      }
    }
  }
  
  if (tissue == "saliva") {
    if (methylation_array %in% c("epic", "epic2", "450k")) {
      if (age == ">18") {
        return(list(zheng = centEpiFibIC.m))
      } else if (age == "<18") {
        return(list(meffil = "saliva gse147318"))
      } else if (age == "all") {
        stop("No reference available for saliva tissue in all age range with EPIC array.")
      }
    }
  }
  
  if (tissue == "buccal") {
    if (methylation_array %in% c("epic", "epic2", "450k")) {
      if (age == ">18") {
        return(list(zheng = centEpiFibIC.m))
      }
    }
  }
}