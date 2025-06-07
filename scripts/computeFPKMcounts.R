hg38_protein_coding_genes <- fread("hg38_protein_coding_genes.bed")
mm10_protein_coding_genes <- fread("mm10_protein_coding_genes.bed")

species=c("human","mouse")

if (species="human") {
raw_counts <- fread("featureCounts_output.counts") %>%
  as.data.frame() %>%
  dplyr::filter(Geneid %in% hg38_protein_coding_genes$Geneid)
} else if (species="mouse") {
  raw_counts <- fread("featureCounts_output.counts") %>%
  as.data.frame() %>%
  dplyr::filter(Geneid %in% mm10_protein_coding_genes$Geneid)
}

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

# FPKM counts
mcols(dds)$basepairs <- gene_length

FPKM_counts <- DESeq2::fpkm(dds,robust = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("Geneid")
