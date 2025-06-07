raw_counts <- fread("featureCounts_output.counts") %>%
  as.data.frame() 

gene_length <- raw_counts$Length

mat <- raw_counts %>%
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()

metadata <- data.frame(kind = colnames(mat))

metadata$condition <- metadata$kind
# run DESeq2
metadata <- tibble::column_to_rownames(metadata, "kind")
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = metadata,
  design = ~condition
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

dds <- DESeq(dds)

resNames <- resultsNames(dds) %>% as.data.frame()

resLFC <- lfcShrink(dds=dds,coef=2,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("Geneid") %>% na.omit()