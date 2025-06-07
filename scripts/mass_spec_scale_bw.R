# R function to generate norm factors
ms_norm <- function(bw, MS,species, exported_file_name,outdir) {
  bw <- import.bw(bw)
  if (species == "mouse"){
    blacklist <- import.bed("mm10_blacklist.bed")
  } else if (species == "human"){
    blacklist <- import.bed("hg38_blacklist.bed")
  }
  bw <- bw[!overlapsAny(bw,blacklist)]
  bw <- bw[!is.na(bw$score)]
  N <- as.numeric(length(bw))
  sum_bw <- as.numeric(sum(bw$score))
  MS <- as.numeric(paste0(MS))
  norm_factor <- MS * N / (sum_bw)
  bw$score <- bw$score * norm_factor
  filt_bw <- bw
  filt_bw$score[filt_bw$score > 100] <- as.numeric(100)
  
  print(paste0("Max value = ", max(filt_bw$score)))
  print(paste0("Post-processing quantiles:", quantile(filt_bw$score)))
  print(paste0("Norm factor = ", norm_factor))
  print(paste0("Signal mean = ", mean(filt_bw$score)))
  export.bw(filt_bw, sprintf("%s/%s",outdir,exported_file_name))
}

# example 
ms_norm(
  bw = "samp.bw",
  MS = MS_value,
  species = "mouse",
  outdir = outdir,
  exported_file_name = "samp"
)