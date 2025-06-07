volc <- function(r, x, y, ylab,xlab, ttl,cutoff) {
  d <- as.data.frame(r) %>%
    dplyr::rename(x = !!x, y = !!y) %>% dplyr::filter(!is.na(x) & !is.na(y)) %>%
    mutate(kind = case_when((x > cutoff & y < .05) ~ paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = ""),
                            (x < -cutoff & y < .05) ~ paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = ""),(cutoff > abs(x) & y < .05) ~ ('FDR<0.05'),
                            y > .05 ~ 'NS'),
           y = -log10(y))
  d$kind <- factor(d$kind, levels = c(paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = ""),paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = ""),'FDR<0.05','NS'))
  ct_up <- d %>% 
    dplyr::filter(kind == paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = "")) %>%
    mutate(up = x > 0) %>%
    dplyr::count(up) %>%
    mutate(x = ifelse(up, Inf, -Inf),
           y = Inf,
           h = as.numeric(up))
  ct_down <- d %>% 
    dplyr::filter(kind == paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = "")) %>%
    mutate(up = x > 0) %>%
    dplyr::count(up) %>%
    mutate(x = ifelse(up, Inf, -Inf),
           y = Inf,
           h = as.numeric(up))
  ggplot(d, aes(x, y, color = kind)) +
    geom_vline(xintercept = c(-cutoff, cutoff), linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05),linetype = 'dashed') +
    geom_point(size=0.2) +
    geom_label(aes(x = x, y = y, label = n, hjust = h),
              vjust = 1, data = ct_up, inherit.aes = F,size=2.5) + scale_y_continuous() + geom_label(aes(x = x, y = y, label = n, hjust = h),
              vjust = 1, data = ct_down, inherit.aes = F,size=2.5) +
    # geom_text_repel(aes(label = TE_subfamily), data = d %>% dplyr::filter((kind %in% c(paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = ""),paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = "")))),max.overlaps = 5,show.legend = F, min.segment.length = 0) +
    scale_color_manual(values = c('purple','darkorange4','forestgreen','darkgrey'),labels=c("FDR<0.05\nlog2FC>2","FDR<0.05\nlog2FC<-2","FDR<0.05","NS")) +
    labs(x = xlab, y = ylab, title = ttl) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=6,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=6,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=6,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 8, color = "black",family = "Helvetica"),
    # legend
    legend.text=element_text(size=6,family = "Helvetica",color = "black"),
    legend.title=element_blank(),
    legend.background = element_rect(fill="white",size = 5),
    legend.key=element_rect(fill="white",size = 5),
    legend.margin = unit(c(0,0,0,0),"mm"),
    legend.spacing = unit(c(0,0,0,0),"mm"),
    legend.position = "top",
    legend.justification = "right",
    legend.box.margin = margin(-5,-5,-5,-5),
    legend.key.size = unit(0,"mm"),
    plot.margin = unit(c(2.4,2.4,2.4,2.4), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'))
}

volc(r=resLFC,x = 'log2FoldChange',y = 'padj',ylab = '-log10(padj)',xlab = "Log2 (FC) gene expression",ttl = "",cutoff = 2)
