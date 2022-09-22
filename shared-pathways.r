library(ggpubr)

allm <- read.table(file = "allm.tsv", sep = "\t", header = T)
allc <- read.table(file = "allc.tsv", sep = "\t", header = T)

#Adding a column for shared/unique pathways:
allm$shared <- allm$progress
allm[allm$pathway %in% intersect(allm[allm$progress == "progressing",]$pathway, allm[allm$progress == "recovering",]$pathway),]$shared <- "shared"

allc$shared <- allc$progress
allc[allc$pathway %in% intersect(allc[allc$progress == "progressing",]$pathway, allc[allc$progress == "recovering",]$pathway),]$shared <- "shared"

allm$shared <- factor(allm$shared, levels = c("shared","progressing", "recovering"))
allc$shared <- factor(allc$shared, levels = c("shared","progressing", "recovering"))

#Setting facet labels
labels <- c("shared" = "Common Pathways", "progressing" = "Pathways in Progressing Patients", "recovering" = "Pathways in Recovering Patients")

#Plotting bubble plots, divided by progress, including the shared pathways:
shared <- ggplot(transform(allm[allm$shared == "shared",]),
                 aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +
  facet_grid(vars(shared), labeller = labeller(shared = labels), drop = TRUE, space = "free") +
  ylab("Pathways") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_size(range = c(3, 8), name="genes") +
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4), limits = c(-4,4)) +
  theme_bw() +
  theme(plot.subtitle=element_text(size=12, face="italic", color="black")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

progressing <- ggplot(transform(allm[allm$shared == "progressing",]),
                      aes(x = NES, y = reorder(pathway,NES), size = size)) +
  geom_point(alpha=0.7, color='#F8766D') +  
  facet_grid(vars(shared), labeller = labeller(shared = labels), drop = TRUE, space = "free") +
  theme(strip.text.x = element_blank()) +
  ylab("Pathways") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_size(range = c(3, 8), name="genes") +
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4), limits = c(-4,4)) +
  theme_bw() + 
  theme(axis.text=element_text(size=12), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.title=element_text(size=14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


recovering <- ggplot(transform(allm[allm$shared == "recovering",]),
                     aes(x = NES, y = reorder(pathway,NES), size = size)) +
  geom_point(alpha=0.7, color='#00BFC4') +  
  facet_grid(vars(shared), labeller = labeller(shared = labels), drop = TRUE, space = "free") +
  theme(strip.text.x = element_blank()) +
  ylab("Pathways") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_size(range = c(3, 8), name="genes") +
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4), limits = c(-4,4)) +
  theme_bw() +
  theme(axis.text=element_text(size=12), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.title=element_text(size=14))

bubble <- ggarrange(shared + rremove("xlab") + rremove("ylab"), progressing + rremove("xlab") + rremove("ylab"), recovering + rremove("ylab"),
                    ncol=1, nrow=3, common.legend = TRUE, legend="right",
                    labels = NULL, heights = c(0.5, 0.25, 0.25),
                    align = "v", 
                    font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))

bubble <- annotate_figure(bubble, left = textGrob("Pathways", rot = 90, vjust = 1, gp = gpar(fontsize = 14)))

tiff(filename = "graphs/pathways-bubble-moderate.tiff",
     width = 15, height = 9, units = "in", res = 300)
bubble
dev.off()
