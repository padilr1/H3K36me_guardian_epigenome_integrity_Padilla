# plot DNA methylation within genes
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
computeMatrix scale-regions -R ${regions} \
-S PA.bw TKO.bw QKO.bw QuiKO.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
-p 6 --samplesLabel PA TKO QKO QuiKO --verbose \
-o out.mat.gz
Height=3.5
Width=4.5
plotProfile -m out.mat.gz -o aggregate_profile.pdf --dpi 600 --colors blue darkkhaki orange maroon --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "% CpG methylation" --startLabel "TSS" --endLabel "TES" --legendLocation best --labelRotation 90 --plotFileFormat pdf
```

# plot H3K4me1 within genes
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
computeMatrix scale-regions -R ${regions} \
-S PA_H3K4me1_merged.bw TKO_H3K4me1_merged.bw QKO_H3K4me1_merged.bw QuiKO_H3K4me1_merged.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
-p 6 --samplesLabel PA TKO QKO QuiKO --verbose \
-o out.mat.gz
Height=3.5
Width=4.5
plotProfile -m out.mat.gz -o aggregate_profile.pdf --dpi 600 --colors blue darkkhaki orange maroon --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "MS-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation best --labelRotation 90 --plotFileFormat pdf
```

# plot H3K9me3 within genes
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
computeMatrix scale-regions -R ${regions} \
-S PA_H3K9me3_merged.bw DKO_H3K9me3_merged.bw TKO_H3K9me3_merged.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
-p 6 --samplesLabel PA DKO TKO --verbose \
-o out.mat.gz
Height=3.5
Width=4.5
plotProfile -m out.mat.gz -o aggregate_profile.pdf --dpi 600 --colors blue darkkhaki orange maroon --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "MS-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation best --labelRotation 90 --plotFileFormat pdf
```

# plot H3K27ac within genes
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
computeMatrix scale-regions -R ${regions} \
-S PA_H3K27ac_merged.bw TKO_H3K27ac_merged.bw QKO_H3K27ac_merged.bw QuiKO_H3K27ac_merged.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
-p 6 --samplesLabel PA TKO QKO QuiKO --verbose \
-o out.mat.gz
Height=3.5
Width=4.5
plotProfile -m out.mat.gz -o aggregate_profile.pdf --dpi 600 --colors blue darkkhaki orange maroon --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "MS-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation best --labelRotation 90 --plotFileFormat pdf
```

# plot H3K4me3 within genes
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
computeMatrix scale-regions -R ${regions} \
-S PA_H3K4me3_merged.bw TKO_H3K4me3_merged.bw QKO_H3K4me3_merged.bw QuiKO_H3K4me3_merged.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
-p 6 --samplesLabel PA TKO QKO QuiKO --verbose \
-o out.mat.gz
Height=3.5
Width=4.5
plotProfile -m out.mat.gz -o aggregate_profile.pdf --dpi 600 --colors blue darkkhaki orange maroon --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "MS-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation best --labelRotation 90 --plotFileFormat pdf
```

# peakiness scores for H3K4me1
```r
# sample script
bl <- import.bed("mm10_blacklist.bed")
rx <- fread("H3K4me1_rx.csv")
b <- fread("10T_PA_rep1_H3K4me1.1kb.bed",col.names = c('chr', 'start', 'end', 'score'))
samp <- "PA_1"
librarySize <- rx$chip[rx$samp == samp]
gr <- makeGRangesFromDataFrame(b,keep.extra.columns = TRUE)
gr$score <- (gr$score/librarySize)*1000000
top1 <- gr[gr$score > (quantile(gr$score,c(0.99)))]
top1_filt <- top1[!overlapsAny(top1,bl)]
p <- data.frame(samp=samp,
                score=mean(top1_filt$score))
# for subsequent samples
b <- fread("10T_PA_rep2_H3K4me1.1kb.bed",col.names = c('chr', 'start', 'end', 'score'))
samp <- "PA_2"
librarySize <- rx$chip[rx$samp == samp]
gr <- makeGRangesFromDataFrame(b,keep.extra.columns = TRUE)
gr$score <- (gr$score/librarySize)*1000000
top1 <- gr[gr$score > (quantile(gr$score,c(0.99)))]
top1_filt <- top1[!overlapsAny(top1,bl)]
s <- data.frame(samp=samp,
                score=mean(top1_filt$score))
p <- p %>% rows_insert(s)

# after running the previous scripts for all samples, then refactor the dataframe based on desired order of samples
p$condition <- factor(p$condition,levels=c("QuiKO","QKO","TKO","PA"))

stats <- aggregate(score ~ condition, p, function(x) c(mean = mean(x), sd = sd(x)))

s <- stats$score %>% as.data.frame() %>% mutate(condition = stats$condition)

ggplot() + stat_summary(mapping=aes(x=condition,y=score,fill=condition),data = p,geom="col",fun=mean,show.legend = FALSE) +
  geom_jitter(mapping=aes(x=condition,y=score,fill=condition),data = p,show.legend = FALSE,size=0.7) +
  geom_errorbar(mapping = aes(y=mean,x=condition,ymin = mean - sd, ymax = mean + sd),data = s, width = 0.2) +
  coord_flip() +
  scale_fill_manual(values=c("maroon3","orange","khaki","blue")) +
  labs(x="",y="") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.y = element_text(family="Helvetica",size=9,colour = "black"),
    axis.text.x = element_text(family="Helvetica",size=9,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=9,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 9),
    plot.margin=unit(c(0,0,0,0), "cm"))

ggsave(filename = "H3K4me1_peakiness_scores.pdf",path="figs",device = "pdf",dpi = 600,bg="white")
```

# plot PA, QKO and QuiKO H3K9me3 within genes
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
computeMatrix scale-regions -R ${regions} \
-S PA_H3K9me3_merged.bw QKO_H3K9me3_merged.bw QuiKO_H3K9me3_merged.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
-p 6 --samplesLabel PA QKO QuiKO --verbose \
-o out.mat.gz
Height=3.5
Width=4.5
plotProfile -m out.mat.gz -o aggregate_profile.pdf --dpi 600 --colors blue orange maroon --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "MS-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation best --labelRotation 90 --plotFileFormat pdf
```