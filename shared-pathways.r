library(ggpubr)

allm <- read.table(file = "allm.tsv", sep = "\t", header = T)

#Adding a column for shared/unique pathways:
allm <- allm %>% mutate(group = progress, progress = NULL)
allm$shared <- allm$group
allm[allm$pathway %in% intersect(allm[allm$group == "progressing",]$pathway, allm[allm$group == "recovering",]$pathway),]$shared <- "shared"

allm$shared <- factor(allm$shared, levels = c("shared","progressing", "recovering"))

allm$category <- "immune"
allm$category[allm$pathway == "MYC TARGETS V1" | allm$pathway == "MYC TARGETS V2" | allm$pathway == "P53 PATHWAY" | allm$pathway == "G2M CHECKPOINT" | allm$pathway == "E2F TARGETS"] <- "proliferation"
allm$category[allm$pathway == "KRAS SIGNALING DN" | allm$pathway == "TNFA SIGNALING VIA NFKB" | allm$pathway == "MTORC1 SIGNALING" | allm$pathway == "ANDROGEN RESPONSE"] <- "signaling"
allm$category[allm$pathway == "FATTY ACID METABOLISM" | allm$pathway == "OXIDATIVE PHOSPHORYLATION" | allm$pathway == "GLYCOLYSIS"] <- "metabolism"
allm$category[allm$pathway == "HYPOXIA" | allm$pathway == "UNFOLDED PROTEIN RESPONSE" | allm$pathway == "EPITHELIAL MESENCHYMAL TRANSITION" | allm$pathway == "APICAL JUNCTION"] <- "other"
allm$category <- factor(allm$category, levels = c("immune", "signaling", "metabolism", "proliferation", "other"))


#Setting facet labels
labels <- c("immune" = "Immune", "proliferation" = "Proliferation", "signaling" = "Signaling",
            "metabolism" = "Metabolic", "other" = "Other")

#Plotting bubble plots, divided by progress, including the shared pathways:
immune <- ggplot(transform(allm[allm$category == "immune",]),
                 aes(x = NES, y = reorder(pathway,NES), color = group, size = size)) +
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
        axis.ticks.x=element_blank(),
        strip.text.y = element_text(size = 14)) 

proliferation <- ggplot(transform(allm[allm$category == "proliferation",]),
                        aes(x = NES, y = reorder(pathway,NES), color = group, size = size)) +
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
        axis.ticks.x=element_blank(),
        strip.text.y = element_text(size = 14))

signaling <- ggplot(transform(allm[allm$category == "signaling",]),
                      aes(x = NES, y = reorder(pathway,NES), color = group, size = size)) +
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
        axis.ticks.x=element_blank(),
        strip.text.y = element_text(size = 14))

metabolism <- ggplot(transform(allm[allm$category == "metabolism",]),
                     aes(x = NES, y = reorder(pathway,NES), color = group, size = size)) +
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
        axis.ticks.x=element_blank(),
        strip.text.y = element_text(size = 14))

other <- ggplot(transform(allm[allm$category == "other",]),
                     aes(x = NES, y = reorder(pathway,NES), color = group, size = size)) +
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
        strip.text.y = element_text(size = 14))



bubble <- ggarrange(immune + rremove("xlab") + rremove("ylab"), proliferation + rremove("xlab") + rremove("ylab"),
                    signaling + rremove("xlab") + rremove("ylab"), metabolism + rremove("xlab") + rremove("ylab"), 
                    other + rremove("ylab"),
                    ncol=1, nrow=5, common.legend = TRUE, legend="right",
                    labels = NULL, 
                    align = "v", 
                    font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))

bubble <- annotate_figure(bubble, left = text_grob("Pathways", rot = 90, vjust = 1, size = 14))

tiff(filename = "graphs/pathways-category-bubble-moderate.tiff",
     width = 15, height = 10, units = "in", res = 300)
bubble
dev.off()
