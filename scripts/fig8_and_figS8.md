# compute saddleplots
```r
saddlePlot <- function(samp,cline,ttl,fontsize=6,isWT=TRUE){
  ll <- 2
  mi <- 2
  ma <- 50
  bsz <- '100000'
  
  cis <- sprintf('%s_%s_%s.%s.saddledump.npz',cline,samp, bsz, "cis") %>%
    np$load(allow_pickle = T)
  saddledata <- cis['saddledata'] %>% as.data.frame() %>% mutate(y = n():1) %>% pivot_longer(-y, names_to = 'x', values_to = 'z') %>%
    mutate(x = sub('^V', '', x))
  saddledata$x <- as.integer(saddledata$x)
  sdl <- saddledata
  
  bb <- sdl %>% dplyr::filter(x <10 & y >42) %>% dplyr::filter(z >= as.numeric(quantile(.$z,0.80)))
  aa <- sdl %>% dplyr::filter(x > 42 & y < 10) %>% dplyr::filter(z >= as.numeric(quantile(.$z,0.80)))
  ab_total <- sdl %>% dplyr::filter(x <10 & y < 10 | x > 42 & y > 42) %>% dplyr::filter(z <= as.numeric(quantile(.$z,0.20)))

  sum(bb$z)/sum(ab_total$z)
  sum(aa$z)/sum(ab_total$z)
  bb_str <- round(sum(bb$z)/sum(ab_total$z),digits = 2)
  aa_str <- round(sum(aa$z)/sum(ab_total$z),digits=2)
  
  if(isWT==TRUE){
    p <- sdl %>%
      ggplot(aes(x, y, fill = z)) +
      geom_raster() +
      scale_fill_gradientn('O/E', colors = pals::coolwarm(25),
                           oob = scales::squish,
                           guide = guide_colorbar(barwidth = 0.3,barheight = 4), limits = c(1/ll,ll),
                           breaks = scales::pretty_breaks(4), trans = 'log2') +
      annotate("text", x = 5, y = 38, label = "B-B",size=2) +
      annotate("text", x = 48, y = 15, label = "A-A",size=2) +
      annotate("text", x = 48, y = 38, label = "B-A",size=2) +
      annotate("text", x = 5, y = 15, label = "A-B",size=2) +
      annotate("text",label=paste0(bb_str),x=5,y=47,size=2,colour="white")+
      annotate("text",label=paste0(aa_str),x=47,y=5,size=2,colour="black")+
      geom_rect(xmin=0.5,xmax=10,ymin=0.5,ymax=10,fill=NA,colour="white",linetype="dashed",linewidth=0.1) +
      geom_rect(xmin=0.5,xmax=10,ymin=42,ymax=52.5,fill=NA,colour="white",linetype="dashed",linewidth=0.1) +
      geom_rect(xmin=42,xmax=52.5,ymin=42,ymax=52.5,fill=NA,colour="white",linetype="dashed",linewidth=0.1) +
      geom_rect(xmin=42,xmax=52.5,ymin=0.5,ymax=10,fill=NA,colour="white",linetype="dashed",linewidth=0.1) +
      labs(title=ttl)+
      theme(plot.background = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            plot.title=element_text(size=fontsize+0.5,family = "Helvetica",color = "black",hjust = 0.5,vjust=-5),
            # legend
            legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
            legend.title=element_text(size=fontsize,family = "Helvetica",color = "black"),
            legend.background = element_rect(fill="white",size = 4),
            legend.key=element_rect(fill="white",size = 4),
            legend.margin = unit(c(0,0,0,0),"mm"),
            legend.spacing = unit(c(0,0,0,0),"mm"),
            legend.position = "right",
            legend.justification = "right",
            plot.margin = unit(c(2.4,2.4,2.4,2.4), "mm"),
            panel.spacing = unit(0.001,'cm'),
            panel.spacing.y = unit(0.001,'cm'),
            panel.spacing.x = unit(0.001,'cm'))
  }
  
  if(isWT==FALSE){
  p <-  sdl %>%
    ggplot(aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_gradientn('O/E', colors = pals::coolwarm(25),
                         oob = scales::squish,
                         guide = guide_colorbar(barwidth = 0.3,barheight = 4), limits = c(1/ll,ll),
                         breaks = scales::pretty_breaks(4), trans = 'log2') +
    annotate("text",label=paste0(bb_str),x=5,y=47,size=2,colour="white")+
    annotate("text",label=paste0(aa_str),x=47,y=5,size=2,colour="black")+
    labs(title=ttl)+
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title=element_text(size=fontsize+0.5,family = "Helvetica",color = "black",hjust = 0.5,vjust=-5),
          # legend
          legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
          legend.title=element_text(size=fontsize,family = "Helvetica",color = "black"),
          legend.background = element_rect(fill="white",size = 4),
          legend.key=element_rect(fill="white",size = 4),
          legend.margin = unit(c(0,0,0,0),"mm"),
          legend.spacing = unit(c(0,0,0,0),"mm"),
          legend.position = "right",
          legend.justification = "right",
          plot.margin = unit(c(2.4,2.4,2.4,2.4), "mm"),
          panel.spacing = unit(0.001,'cm'),
          panel.spacing.y = unit(0.001,'cm'),
          panel.spacing.x = unit(0.001,'cm'))
  }
  
 plot(p)
}

saddlePlot(samp = "PA",cline = "mMSC",ttl = "mMSC PA",isWT = TRUE)
ggsave(filename = "mMSC_PA.saddleplot.100kb.pdf",path = "outdir",width = 4.8,height=3.8,units="cm",dpi = 600,device = "pdf",bg = "white")

saddlePlot(samp = "QuiKO",cline = "mMSC",ttl = "mMSC QuiKO",isWT = FALSE)
ggsave(filename = "mMSC_QuiKO.saddleplot.100kb.pdf",path = "outdir",width = 4.8,height=3.8,units="cm",dpi = 600,device = "pdf",bg = "white")

saddlePlot(samp = "PAmerged",cline = "Cal27",ttl = "Cal27 PA",isWT = TRUE)
ggsave(filename = "Cal27_PA_merged.saddleplot.100kb.pdf",path = "outdir",width = 4.8,height=3.8,units="cm",dpi = 600,device = "pdf",bg = "white")

saddlePlot(samp = "H3K36MOE",cline = "Cal27",ttl = "Cal27 H3K36M-OE",isWT = FALSE)
ggsave(filename = "Cal27_H3K36MOE.saddleplot.100kb.pdf",path = "outdir",width = 4.8,height=3.8,units="cm",dpi = 600,device = "pdf",bg = "white")
```

# compute compartment scores
```r
plotCompartmentScores <- function(dir,fontsize=7,factorOrder,colors=c("blue","maroon"),newSampLabels){
  bsz <- '100000'
  pd <- data.frame()
  m <- list.files(dir, pattern = paste0('^', bsz, '.cis.vecs.tsv'), 
                  full.names = T, recursive = T) %>%
    setNames(., basename(dirname(.))) %>%
    lapply(function(x) fread(x)$E1) %>%
    bind_cols() %>%
    na.omit()
  
  m_plot <- m %>%
    mutate(idx = 1:n()) %>%
    pivot_longer(-idx, names_to = 'samp', values_to = 'v') 
  
  m_plot$samp <- factor(x=m_plot$samp,levels=factorOrder)
  
  m_plot %>% ggplot(aes(x = samp,y  = v, color = samp, fill = samp)) +
    geom_hline(yintercept = 0) +
    geom_violin(color = NA, alpha = .5) +
    stat_pointinterval() +
    scale_x_discrete(name = "",
                     labels=c(newSampLabels)) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    labs(y = 'Compartment score') +
    coord_cartesian(ylim = c(-2, 2)) +
    scale_y_continuous(breaks = -2:2) +
    theme(plot.background = element_blank(),
          plot.title = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none',
          axis.line.x = element_line(color = 'black'),
          panel.grid.major.y = element_line(color = 'grey75', linetype = 'dashed'),
          strip.background = element_blank(),
          strip.text = element_text(color = 'white'),
          axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black"),
          axis.text.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
          axis.title.y = element_text(size=fontsize,family="Helvetica",colour = "black"),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(color = 'black'))
  
}

plotCompartmentScores(dir = "input_dir",fontsize = 6,factorOrder = c("PA","QuiKO"),colors = c("blue","maroon"),newSampLabels = c("PA","QuiKO"))
ggsave(filename = "mMSC_PA_QuiKO_compartmentScores.pdf",device = "pdf",path = "outdir",width = 3.5,height = 3.5,units = "cm",dpi = 600,bg = "white")

plotCompartmentScores(dir = "input_dir",fontsize = 6,factorOrder = c("PA","H3K36M_OE"),colors = c("blue","chocolate4"),newSampLabels = c("PA","H3K36M-OE"))
ggsave(filename = "Cal27_PA_merged_K36MOE_compartmentScores.pdf",device = "pdf",path = "outdir",width = 3.5,height = 3.5,units = "cm",dpi = 600,bg = "white")
```

# compartment shifts from A to B
```r
diffCompartmentAnalysis <- function(samp1_bw,samp2_bw,samp1Lab,samp2Lab,fontsize=6){
  samp1 <- import.bw(samp1_bw)
  samp2 <- import.bw(samp2_bw)
  # subset
  samp1_sub <- subsetByOverlaps(samp1,samp2) %>% as.data.frame() %>% na.omit() %>% mutate(coordinates=sprintf("%s:%d-%d",.$seqnames,.$start,.$end)) %>% dplyr::select(coordinates,score) %>% `names<-`(c('coordinates', 'samp1'))
  samp2_sub <- subsetByOverlaps(samp2,samp1) %>% as.data.frame() %>% na.omit() %>% mutate(coordinates=sprintf("%s:%d-%d",.$seqnames,.$start,.$end)) %>% dplyr::select(coordinates,score) %>% `names<-`(c('coordinates', 'samp2'))
  combined <- left_join(samp1_sub,samp2_sub,by="coordinates") %>% na.omit() %>% mutate(change = case_when(
    samp1 > 0 & samp2 < 0 ~ "A->B",
    samp1 < 0 & samp2 > 0 ~ "B->A",
    TRUE ~ "No change"
  ))
  
  samp1_counts <- samp1 %>% as.data.frame() %>% na.omit() %>% mutate(compartment=case_when(
    .$score < 0 ~ "B",
    .$score > 0 ~ "A"
  ))
  samp1_frequency <- samp1_counts %>% dplyr::count(compartment) %>% mutate(cond=samp1Lab)
  
  samp2_counts <- samp2 %>% as.data.frame() %>% na.omit() %>% mutate(compartment=case_when(
    .$score < 0 ~ "B",
    .$score > 0 ~ "A"
  ))
  samp2_frequency <- samp2_counts %>% dplyr::count(compartment) %>% mutate(cond=samp2Lab)
  
  combined <- rbind(samp1_frequency,samp2_frequency)
  combined$cond <- factor(levels = c(samp1Lab,samp2Lab),x = combined$cond)
  
  p <- ggplot(data=combined,aes(x=cond,y=n,fill=compartment)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values=c("darkorange3","darkviolet")) +
    labs(x="",y="Count",title="Compartments at\n 100kb windows") +
    theme(
      # panel
      panel.background = element_blank(),
      axis.line.x.bottom  = element_line(colour="black"),
      axis.line.x.top  = element_line(colour="black"),
      axis.line.y.left = element_line(colour="black"),
      plot.title = element_text(hjust = 0.5,color = "black",size=6,family="Helvetica"),
      # axis
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      # labels
      strip.background =element_blank(),
      strip.text.x = element_blank(),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
      # legend
      legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
      legend.title=element_text(size=fontsize,family="Helvetica",color="black"),
      legend.background = element_rect(fill="white"),
      legend.key=element_rect(fill="white",size=0.1),
      # margins
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
      panel.spacing = unit(0.5,'cm'),
      panel.spacing.y = unit(0.5,'cm'),
      panel.spacing.x = unit(0.5,'cm'),
      legend.position="right",
      legend.justification="right",
      legend.box.spacing = unit(-0.001, "cm"))
  plot(p)
}

diffCompartmentAnalysis(samp1_bw = "PA_100000.cis.bw",samp2_bw = "TKO_100000.cis.bw",samp1Lab = "PA",samp2Lab = "TKO")
ggsave(filename = "mMSC_PA_TKO.compartment_totals_100kb_res.pdf",path="outdir",device = "pdf",units = "cm",width = 3.3,height=4,dpi = 600,bg="white")

diffCompartmentAnalysis(samp1_bw = "PA_100000.cis.bw",samp2_bw = "H3K36MOE_100000.cis.bw",samp1Lab = "PA",samp2Lab = "H3K36M-OE")
ggsave(filename = "Cal27_PA_H3K36MOE.compartment_totals_100kb_res.pdf",path="outdir",device = "pdf",units = "cm",width = 3.3,height=4,dpi = 600,bg="white")
```

# Total TADs & TAD strengths
```r
getTADboundaries <- function(samp1Lab,samp2Lab,samp1_insMat,samp2_insMat,colors=c("blue","maroon"),fontsize){
  print("Using 100kb insulation boundaries")
  # 100kb
  samp1 <- fread(samp1_insMat) %>% dplyr::filter(is_boundary_100000 == "TRUE") %>% mutate(coordinates = sprintf("%s:%d-%d",.$chrom,.$start,.$end)) %>% dplyr::select(c("boundary_strength_100000","coordinates")) %>% `names<-`(c(samp1Lab, 'coordinates')) %>% na.omit()
  
  samp2 <- fread(samp2_insMat) %>% dplyr::filter(is_boundary_100000 == "TRUE") %>% mutate(coordinates = sprintf("%s:%d-%d",.$chrom,.$start,.$end)) %>% dplyr::select(c("boundary_strength_100000","coordinates")) %>% `names<-`(c(samp2Lab, 'coordinates')) %>% na.omit()
  
  combined <- left_join(samp1,samp2,by="coordinates") %>% na.omit() 
  lost_TAD_boundaries <- dplyr::anti_join(samp1,samp2,by="coordinates") %>% na.omit()
  gained_TAD_boundaries <- dplyr::anti_join(samp2,samp1,by="coordinates") %>% na.omit()
  
  # total number of TADs
  total_TADs <- tibble(
    cond=c(samp1Lab,samp2Lab),
    total_TADs= c(nrow(samp1),nrow(samp2))
  )
  
  total_TADs$cond <- factor(total_TADs$cond,levels = c(samp1Lab,samp2Lab))
  
  # TAD strength
  TAD_strength <- cbind(samp1 %>% dplyr::select(1) %>% slice_sample(n=nrow(samp1),replace=FALSE),samp2 %>% dplyr::select(1)%>% slice_sample(n=nrow(samp1),replace=FALSE)) %>% pivot_longer(c(1:2),names_to="cond",values_to="TAD_strength")
  
  TAD_strength$cond <- factor(TAD_strength$cond,levels = c(samp1Lab,samp2Lab))
  
  TADstrength_plot <- ggplot(data=TAD_strength,aes(x=cond,y=log2(TAD_strength),fill=cond)) +
    geom_boxplot(show.legend = FALSE,outlier.size = 0.1,linewidth=0.2) +
    scale_fill_manual(values=c(colors),name="") +
    labs(x="",y="log2(TAD strength)\n at 100kb windows",title="") +
    theme(
      # panel
      panel.background = element_blank(),
      axis.line.x.bottom  = element_line(colour="black"),
      axis.line.x.top  = element_line(colour="black"),
      axis.line.y.left = element_line(colour="black"),
      plot.title = element_text(hjust = 0.5,color = "black",size=fontsize,family="Helvetica"),
      # axis
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      # labels
      strip.background =element_blank(),
      strip.text.x = element_blank(),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
      # legend
      legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
      legend.title=element_text(size=fontsize,family="Helvetica",color="black"),
      legend.background = element_rect(fill="white"),
      legend.key=element_rect(fill="white"),
      legend.key.size = unit(2,"mm"),
      # margins
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
      panel.spacing = unit(0.5,'cm'),
      panel.spacing.y = unit(0.5,'cm'),
      panel.spacing.x = unit(0.5,'cm'),
      legend.position="top",
      legend.justification="right",
      legend.box.spacing = unit(-0.001, "cm"))
  
  totalTAD_plot <- ggplot(data=total_TADs,aes(x=cond,y=total_TADs,fill=cond)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values=c(colors),name="") +
    labs(x="",y="Total number of TADs\n at 100kb windows",title="") +
    theme(
      # panel
      panel.background = element_blank(),
      axis.line.x.bottom  = element_line(colour="black"),
      axis.line.x.top  = element_line(colour="black"),
      axis.line.y.left = element_line(colour="black"),
      plot.title = element_text(hjust = 0.5,color = "black",size=fontsize,family="Helvetica"),
      # axis
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      # labels
      strip.background =element_blank(),
      strip.text.x = element_blank(),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
      # legend
      legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
      legend.title=element_text(size=fontsize,family="Helvetica",color="black"),
      legend.background = element_rect(fill="white"),
      legend.key=element_rect(fill="white"),
      legend.key.size = unit(2,"mm"),
      # margins
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
      panel.spacing = unit(0.5,'cm'),
      panel.spacing.y = unit(0.5,'cm'),
      panel.spacing.x = unit(0.5,'cm'),
      legend.position="top",
      legend.justification="right",
      legend.box.spacing = unit(-0.001, "cm"))
  
  return(list(TADstrength_plot=TADstrength_plot,totalTAD_plot=totalTAD_plot))
}

TADboundaries <- getTADboundaries(samp1Lab = "samp1",samp2Lab = "samp2",samp1_insMat = "samp1_10000.100000.tsv",samp2_insMat = "samp2_10000.100000.tsv",fontsize = 5)

mMSC_TADboundaries$totalTAD_plot
ggsave(filename = "100kb_windows_total_number_TADs.pdf",path="outdir",device = "pdf",units = "cm",width = 3.5,height=3.5,dpi = 600,bg="white")

mMSC_TADboundaries$TADstrength_plot
ggsave(filename = "100KB_TAD_strength.pdf",path="outdir",device = "pdf",units = "cm",width = 3.5,height=3.5,dpi = 600,bg="white")
```

# loop distribution
```r
loopDistribution <- function(samp1Lab,samp2Lab,samp1Dots,samp2Dots,fontsize=6,ttl="Loops at 25kb resolution"){
  print("Loops at 25kb resolution")
  
  samp1 <- fread(samp1Dots) %>% mutate(coordinates=sprintf("%s:%d-%d_%s:%d-%d",.$chrom1,.$start1,.$end1,.$chrom2,.$start2,.$end2)) %>%
    as.data.frame() %>% dplyr::select(1:6) %>% mutate(center1=(rowMeans(.[grep("start1|end1", names(.), value = TRUE)]))) %>% mutate(center2=(rowMeans(.[grep("start2|end2", names(.), value = TRUE)]))) %>% mutate(distance=as.integer(center2-center1)) %>% group_by(distance) %>% summarize(freq=n()) %>% `names<-`(c('distance', 'samp1'))
  samp1_total=sum(samp1$samp1)
  
  samp2 <- fread("/Users/padilr1/Documents/10T_downstream/data/HiC_data/dot/QuiKO/25000") %>% mutate(coordinates=sprintf("%s:%d-%d_%s:%d-%d",.$chrom1,.$start1,.$end1,.$chrom2,.$start2,.$end2)) %>%
    as.data.frame()%>% dplyr::select(1:6) %>% mutate(center1=(rowMeans(.[grep("start1|end1", names(.), value = TRUE)]))) %>% mutate(center2=(rowMeans(.[grep("start2|end2", names(.), value = TRUE)]))) %>% mutate(distance=as.integer(center2-center1)) %>% group_by(distance) %>% summarize(freq=n()) %>% `names<-`(c('distance', 'samp2'))
  samp2_total=sum(samp2$samp2)
  
  
  samp1_pct <- samp1 %>% mutate(range = ifelse(distance <= 500000, "short_range", "long_range")) %>% group_by(range) %>%
    summarize(counts=sum(samp1)) %>% mutate(total=samp1_total) %>% mutate(pct=(.$counts/.$total)*100) %>% mutate(cond=samp1Lab)
  samp2_pct <- samp2 %>% mutate(range = ifelse(distance <= 500000, "short_range", "long_range")) %>% group_by(range) %>%
    summarize(counts=sum(samp2)) %>% mutate(total=samp2_total) %>% mutate(pct=(.$counts/.$total)*100)  %>% mutate(cond=samp2Lab)

  # aggregate
  agg <- rbind(samp1_pct,samp2_pct)
  agg$cond <- factor(x=agg$cond,levels=c(samp1Lab,samp2Lab))
  agg$range <- factor(x=agg$range,levels=c("long_range","short_range"))
  
  
  
  ggplot(data=agg,aes(x=cond,y=pct,fill=range)) +
    geom_col() +
    scale_fill_manual(values=c("aquamarine4","darkgoldenrod"),name="",labels=c("long-range\nloops > 500kb","short-range\nloops =< 500kb")) +
    labs(x="",y="Percentage of loops",title=ttl)+
    theme(
      # panel
      panel.background = element_blank(),
      axis.line.x.bottom  = element_line(colour="black"),
      axis.line.x.top  = element_line(colour="black"),
      axis.line.y.left = element_line(colour="black"),
      plot.title = element_text(hjust = 0.5,color = "black",size=fontsize,family="Helvetica",vjust=-5),
      # axis
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      # labels
      strip.background =element_blank(),
      strip.text.x = element_blank(),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
      # legend
      legend.text=element_text(size=fontsize-2,family = "Helvetica",color = "black"),
      legend.title=element_text(size=fontsize,family="Helvetica",color="black"),
      legend.background = element_rect(fill="white"),
      legend.key=element_rect(fill="white"),
      legend.key.size = unit(2,"mm"),
      # margins
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
      panel.spacing = unit(0.5,'cm'),
      panel.spacing.y = unit(0.5,'cm'),
      panel.spacing.x = unit(0.5,'cm'),
      legend.position="top",
      legend.justification="right",
      legend.box.spacing = unit(-0.001, "cm"))
}

loopDistribution(samp1Lab = "samp1",samp2Lab = "samp2",samp1Dots = "samp1_25000",samp2Dots = "samp2_25000",fontsize = 6,ttl = "")

ggsave(filename = "25kb_loopDistribution.pdf",path="outdir",device = "pdf",units = "cm",width = 3,height=4.5,dpi = 600,bg="white")
```

# Gene expression changes associated with loss of long-range loops
```r
samp1Dots = "samp1_25000"
samp2Dots = "samp2_25000"

# replace with mm10 fantom and ensembl enhancers for mMSC
hg38_fantom_enhancers <- import.bed("hg38FantomEnhancers.bed")
hg38_ensembl_enhancers <- import.bed("hg38EnsemblEnhancers.bed")

hg38_enhancers <- c(hg38_fantom_enhancers,hg38_ensembl_enhancers)%>% as.data.frame() %>% dplyr::mutate(start = .$start - 1) %>% dplyr::mutate(enhancerID = paste0("enhancer_",1:nrow(.))) %>% makeGRangesFromDataFrame()

#samp1
samp1 <- fread(samp1Dots) %>% mutate(coordinates=sprintf("%s:%d-%d_%s:%d-%d",.$chrom1,.$start1,.$end1,.$chrom2,.$start2,.$end2)) %>%
    as.data.frame() %>% dplyr::select(1:6) %>% mutate(center1=(rowMeans(.[grep("start1|end1", names(.), value = TRUE)]))) %>% mutate(center2=(rowMeans(.[grep("start2|end2", names(.), value = TRUE)]))) %>% mutate(distance=as.integer(center2-center1)) %>% dplyr::filter(distance > 500000)

samp1_anchor1 <- samp1 %>% dplyr::select(c("chrom1","start1","end1")) %>% makeGRangesFromDataFrame(seqnames.field = "chrom1",start.field = "start1",end.field = "end1")
samp1_anchor2 <- samp1 %>% dplyr::select(c("chrom2","start2","end2")) %>% makeGRangesFromDataFrame(seqnames.field = "chrom2",start.field = "start2",end.field = "end2")

samp1_agg_anchors <- c(samp1_anchor1,samp1_anchor2)

# samp2
samp2 <- fread(samp2Dots) %>% mutate(coordinates=sprintf("%s:%d-%d_%s:%d-%d",.$chrom1,.$start1,.$end1,.$chrom2,.$start2,.$end2)) %>%
    as.data.frame() %>% dplyr::select(1:6) %>% mutate(center1=(rowMeans(.[grep("start1|end1", names(.), value = TRUE)]))) %>% mutate(center2=(rowMeans(.[grep("start2|end2", names(.), value = TRUE)]))) %>% mutate(distance=as.integer(center2-center1)) %>% dplyr::filter(distance > 500000)

samp2_anchor1 <- samp2 %>% dplyr::select(c("chrom1","start1","end1")) %>% makeGRangesFromDataFrame(seqnames.field = "chrom1",start.field = "start1",end.field = "end1")
samp2_anchor2 <- samp2 %>% dplyr::select(c("chrom2","start2","end2")) %>% makeGRangesFromDataFrame(seqnames.field = "chrom2",start.field = "start2",end.field = "end2")

samp2_agg_anchors <- c(samp2_anchor1,samp2_anchor2)

# lost loops
exclusive_samp1LoopAnchors <- samp1_agg_anchors[!overlapsAny(samp1_agg_anchors,samp2_agg_anchors)]
exclusiveEnhancers_at_exlusiveSamp1Anchors <- subsetByOverlaps(hg38_enhancers,exclusive_samp1LoopAnchors,type="within")

# replace with mm10 protein coding genes for mMSC
hg38_protein_coding_genes <- import.bed("hg38_protein_coding_genes.bed")

nearest_hits <- distanceToNearest(exclusiveEnhancers_at_exlusiveSamp1Anchors,hg38_protein_coding_genes )

nearest_gene_indices <- subjectHits(nearest_hits)


filt_hg38_protein_coding_genes <- hg38_protein_coding_genes[nearest_gene_indices] %>% as.data.frame() %>% .[!duplicated(.[,c("name")]),] %>% mutate(ensembl_id = .$name)

load("RNASeq_dds.RData")

resultsNames <- resultsNames(dds) %>% as.data.frame()

resLFC <- lfcShrink(dds=dds,coef=2,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("ensembl_id") %>% na.omit()

filt_resLFC <- resLFC %>% dplyr::filter(ensembl_id %in% filt_hg38_protein_coding_genes$ensembl_id)

# plot
agg <- filt_resLFC %>% mutate(type="H3K36M-OE/PA") %>% dplyr::select(c("type","log2FoldChange"))
fontsize=4.5
sigFontSize=1.8
first_group <- agg
ggplot(data = agg,aes(x=type,y=log2FoldChange,fill=type)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),linewidth=0.3) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.05,alpha=0.7) +
  geom_boxplot(width=0.2, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("chocolate4")) +
  geom_hline(yintercept =0,linetype="dashed",linewidth=0.1) +
  labs(x="",y="Gene expression\nlog2 (FC)",title="") +
  annotate("text",y=1.5, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange > 0,])),size=sigFontSize) +
  annotate("text",y=-1.5, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange < 0,])),size=sigFontSize) +
  geom_hline(yintercept = 2,colour="white")+
  coord_flip() +
    theme(
      # panel
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.5),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black",angle=90,hjust=0.5),
      plot.title = element_text(hjust = 0.5,color = "black",size=fontsize,family="Helvetica"),
      # background
      strip.background =element_rect(fill="white"),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
      # legend
      legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
      legend.title=element_text(size=fontsize,family="Helvetica",color="black"),
      legend.background = element_rect(fill="white"),
      legend.key=element_rect(fill="white"),
      legend.position="right",
      legend.justification="right",
      legend.box.spacing = unit(-0.001, "cm"),
      # margins
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      panel.spacing = unit(0.1,'cm'),
      panel.spacing.y = unit(0,'cm'),
      panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "plot.pdf",path="outdir",device = "pdf",units = "cm",width = 3,height=3,dpi = 600,bg="white")
```