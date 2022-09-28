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

allm$category <- "immune"
allm$category[allm$pathway == "MYC TARGETS V1" | allm$pathway == "MYC TARGETS V2" | allm$pathway == "P53 PATHWAY" | allm$pathway == "G2M CHECKPOINT" | allm$pathway == "E2F TARGETS"] <- "proliferation"
allm$category[allm$pathway == "KRAS SIGNALING DN" | allm$pathway == "TNFA SIGNALING VIA NFKB" | allm$pathway == "MTORC1 SIGNALING" | allm$pathway == "ANDROGEN RESPONSE"] <- "signaling"
allm$category[allm$pathway == "FATTY ACID METABOLISM" | allm$pathway == "OXIDATIVE PHOSPHORYLATION" | allm$pathway == "GLYCOLYSIS"] <- "metabolism"
allm$category[allm$pathway == "HYPOXIA" | allm$pathway == "UNFOLDED PROTEIN RESPONSE"] <- "pathway"
allm$category[allm$pathway == "EPITHELIAL MESENCHYMAL TRANSITION"] <- "development"
allm$category <- factor(allm$category, levels = c("immune", "signaling", "metabolism", "proliferation", "pathway",  "development"))

allc$category <- "immune"
allc$category[allc$pathway == "MYC TARGETS V1" | allc$pathway == "MYC TARGETS V2" | allc$pathway == "P53 PATHWAY" | allc$pathway == "G2M CHECKPOINT" | allc$pathway == "E2F TARGETS"] <- "proliferation"
allc$category[allc$pathway == "KRAS SIGNALING DN" | allc$pathway == "TNFA SIGNALING VIA NFKB" | allc$pathway == "MTORC1 SIGNALING" | allc$pathway == "ANDROGEN RESPONSE"] <- "signaling"
allc$category[allc$pathway == "FATTY ACID METABOLISM" | allc$pathway == "OXIDATIVE PHOSPHORYLATION" | allc$pathway == "GLYCOLYSIS"] <- "metabolism"
allc$category[allc$pathway == "HYPOXIA" | allc$pathway == "UNFOLDED PROTEIN RESPONSE"] <- "pathway"
allc$category[allc$pathway == "EPITHELIAL MESENCHYMAL TRANSITION"] <- "development"
allc$category <- factor(allc$category, levels = c("immune", "signaling", "metabolism", "proliferation", "pathway",  "development"))


#Setting facet labels
labels <- c("immune" = "Immune", "proliferation" = "Proliferation", "signaling" = "Signaling",
            "metabolism" = "Metabolism", "pathway" = "Pathway", "development" = "Development")

#Plotting bubble plots, divided by progress, including the shared pathways:
immune <- ggplot(transform(allm[allm$category == "immune",]),
                 aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +
  facet_grid(vars(category), labeller = labeller(category = labels), drop = TRUE, space = "free") +
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

signaling <- ggplot(transform(allm[allm$category == "signaling",]),
                      aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +  
  facet_grid(vars(category), labeller = labeller(category = labels), drop = TRUE, space = "free") +
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


proliferation <- ggplot(transform(allm[allm$category == "proliferation",]),
                     aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +  
  facet_grid(vars(category), labeller = labeller(category = labels), drop = TRUE, space = "free") +
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


metabolism <- ggplot(transform(allm[allm$category == "metabolism",]),
                     aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +  
  facet_grid(vars(category), labeller = labeller(category = labels), drop = TRUE, space = "free") +
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


pathway <- ggplot(transform(allm[allm$category == "pathway",]),
                     aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +  
  facet_grid(vars(category), labeller = labeller(category = labels), drop = TRUE, space = "free") +
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


development <- ggplot(transform(allm[allm$category == "development",]),
                     aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +  
  facet_grid(vars(category), labeller = labeller(category = labels), drop = TRUE, space = "free") +
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



bubble <- ggarrange(immune + rremove("xlab") + rremove("ylab"), signaling + rremove("xlab") + rremove("ylab"), metabolism + rremove("xlab") + rremove("ylab"), 
                    pathway + rremove("xlab") + rremove("ylab"), proliferation + rremove("xlab") + rremove("ylab"), development + rremove("ylab"),
                    ncol=1, nrow=6, common.legend = TRUE, legend="right",
                    labels = NULL, 
                    align = "v", 
                    font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))

bubble <- annotate_figure(bubble, left = textGrob("Pathways", rot = 90, vjust = 1, gp = gpar(fontsize = 14)))

tiff(filename = "graphs/pathways-bubble-moderate.tiff",
     width = 15, height = 9, units = "in", res = 300)
bubble
dev.off()
