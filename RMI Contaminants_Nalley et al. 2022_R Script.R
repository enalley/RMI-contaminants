### An assessment of the ecological and human health risks 
### associated with measured contaminant concentrations 
### in reef and pelagic fishes from the Republic of the Marshall Islands  
### Nalley et al. 2022

rm(list= ls())
## loading libraries
theme_set(theme_bw())
library("htmltools"); library("webshot"); library("sf"); library("ggplot2")
library("tidyverse"); library("rnaturalearth"); library("rnaturalearthdata")
library("ggspatial"); library("googleway"); library("lwgeom"); library("rgeos")
library("maps"); library("mapdata"); library("ggmap"); library("spData")
library("ggrepel"); library("ggsn"); library("dplyr"); library("lubridate")
library("ggmap"); library("MASS"); library("boot"); library("pscl"); library("forcats")
library("reshape2"); library("grid"); library("formattable"); library("vegan"); library("FSA")

## load data
locations <- read.csv("./Locations.csv", header = T)
trophic <- read.csv("./Trophic.csv", header = T)
trophic_colors <- c("Pelagic" = "red",
                    "Carnivore" = "darkorchid4",
                    "Benthic Invertivore" = "dodgerblue4",
                    "Detritivore" = "darkorange3",
                    "Herbivore" = "seagreen")     
flag_shapes <- c("J" = 1,
                 " " = 16)

### METALS #### 
metals <- read.csv("./RMI_Metals_2.7.21.csv", header = T)
t.metals <- as.data.frame(t(metals))
names(t.metals) <- lapply(t.metals[1,], as.character)
t.metals <- t.metals[-1,]
t.metals$Sample <- rownames(t.metals)
t.metals2 <- t.metals %>% 
  mutate_at(vars(c(2,4,6,8,10,12,14,16,18,20,22,24,26)), as.character) %>% 
  mutate_at(vars(c(2,4,6,8,10,12,14,16,18,20,22,24,26)), as.numeric)
t.metals2 <- arrange(t.metals2, Sample)
t.metals2 <- cbind(locations, t.metals2)
t.metals2 <- inner_join(t.metals2, trophic)
t.metals2$Trophic.Group <- ordered(t.metals2$Trophic.Group, 
                                   levels = c("Detritivore", "Herbivore",
                                              "Benthic Invertivore", "Carnivore", "Pelagic"))
t.metals2$Location <- as.character(t.metals2$Location)
t.metals2[t.metals2 == "MC Oceanside - Kwaj"] <- "MC Kwaj."
t.metals2$Location <- ordered(t.metals2$Location,
                              levels = c("Kwajalein", "MC Kwaj.", "Ebeye",
                                         "Jaluit", "Majuro", "Rongelap", 
                                         "Utirik", "Wotje"))
t.metals2$Species.Ab <- as.character(t.metals2$Species.Ab)
unique(t.metals2$Species.Ab)
t.metals2$Species.Ab <- ordered(t.metals2$Species.Ab,
                                levels = c("C. microrhinos", "C. sordidus" , 
                                           "H. longiceps", "A. triostegus",
                                           "A. lineatus",
                                           "N. lituratus", "L. gibbus",
                                           "S. spiniferum", "L. bohar",
                                           "P. laevis", "V. louti",
                                           "E. polyphekadion",
                                           "T. albacares", "K. pelamis",
                                           "G. unicolor", "E. bipinnulata"))

## nmds
nmds.sites <- as.data.frame(t.metals2[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27)])
nmds.species <- as.data.frame(t.metals2[,c(31,3,5,7,9,11,13,15,17,19,21,23,25,27)])
nmds.trophic <- as.data.frame(t.metals2[,c(30,3,5,7,9,11,13,15,17,19,21,23,25,27)])

metal.nmds.data <- as.data.frame(t.metals2[,c(29,1,31,30,3,5,7,9,11,13,15,17,19,21,23,25,27)])
rownames(metal.nmds.data) <- metal.nmds.data$Sample
metal.nmds.metadata <- metal.nmds.data[,c(1:4)]
metal.nmds.data2 <- metal.nmds.data[,-c(1:4)]

set.seed(1)
metal.nmds <- metaMDS(metal.nmds.data2, "bray", k = 3)
metal.site.scrs <- as.data.frame(scores(metal.nmds, display = "sites"))
metal.site.scrs <- cbind(metal.site.scrs, Location = metal.nmds.metadata$Location)
metal.site.scrs <- cbind(metal.site.scrs, Species = metal.nmds.metadata$Species.Ab)
metal.site.scrs <- cbind(metal.site.scrs, Trophic = metal.nmds.metadata$Trophic.Group)

## plotting
metal.nmds.plot_trophic <- ggplot(metal.site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Trophic)), size = 5, alpha = 0.6) + 
  stat_ellipse(aes(NMDS1, NMDS2, colour = factor(Trophic)), show.legend = FALSE) +
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", 
                                        size = 1, linetype = "solid"),
        legend.position = "top",
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  labs(colour = "") +
  scale_color_manual(values=c("darkorange3", "seagreen", "dodgerblue4",
                              "darkorchid4", "red"))

##
metal.nmds.plot_site <- ggplot(metal.site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Location)), size = 4) + 
  stat_ellipse(aes(NMDS1, NMDS2, colour = factor(Location)), show.legend = FALSE) +
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", 
                                        size = 1, linetype = "solid")) +
  labs(colour = "Location") 
##

## PERMANOVA
## sites 
sites.perm <- adonis(metal.nmds.data2 ~ metal.nmds.metadata$Location,
                     method = "bray", perm = 999)
sites.perm$aov.tab
## looking at BC diff
sites.bray <- vegdist(metal.nmds.data2, method = "bray")
sites.mod <- betadisper(sites.bray, metal.nmds.metadata$Location)
sites.mod

## trophic
trophic.perm <- adonis(metal.nmds.data2 ~ metal.nmds.metadata$Trophic.Group,
                       method = "bray", perm = 999)
trophic.perm$aov.tab
## looking at BC diff
trophic.bray <- vegdist(metal.nmds.data2, method = "bray")
trophic.mod <- betadisper(sites.bray, metal.nmds.metadata$Trophic.Group)
trophic.mod


#### Looking at individual metals ####
### aluminum
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Al, 
           color = Trophic.Group, shape = Al_Flag)) + 
  geom_hline(yintercept = 3294.12, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 941.17, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 235.29, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6, width = 0.25, height = 0) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Aluminum (ug/g wet weight)",
       color = "Trophic Group:") + 
  scale_color_manual(values = trophic_colors) +
  guides(shape = FALSE) + scale_y_log10()

al.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Al),
            Max = max(Al),
            Min = min (Al),
            SD = sd(Al),
            Med = median(Al))

### arsenic
mean(t.metals2$As)
as.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(As)/10,
            Max = max(As)/10,
            Min = min (As)/10)

Tot.as.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(As),
            Max = max(As),
            Min = min (As),
            SD = sd(As),
            Med = median(As))

ggplot(t.metals2, 
       aes(x = Species.Ab,  y = As*.1, 
           color = Trophic.Group, shape = As_Flag)) + 
  geom_hline(yintercept = .07, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .28, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .99, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6, width = 0.25, height = 0) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "top") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Inorganic Arsenic as 10% of Total Arsenic (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  guides(shape = FALSE) 
##
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = As*.01, 
           color = Trophic.Group, shape = As_Flag)) + 
  geom_hline(yintercept = .07, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .28, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .99, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6, width = 0.25, height = 0) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Inorganic Arsenic as 1% Total Arsenic (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  guides(shape = FALSE) 

nrow(t.metals2[t.metals2$As > (.99/.01), ])
nrow(t.metals2[t.metals2$As > (.99/.1), ])
nrow(t.metals2)

as.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(As)/10,
            Max = max(As)/10,
            Min = min (As)/10)

## anova
aov <- aov(As ~ Trophic.Group + Location, data = t.metals2)
summary(aov)
plot(aov)
#
leveneTest(As ~ Trophic.Group, data = t.metals2)
kruskal.test(As ~ Trophic.Group, data = t.metals2)
dunnTest(As ~ Trophic.Group, data = t.metals2, method="bonferroni")
#
leveneTest(As ~ Location, data = t.metals2)
aov2 <- aov(As ~ Location, data = t.metals2)
summary(aov2)


### cadmium
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Cd, 
           color = Trophic.Group, shape = Cd_Flag)) + 
  geom_hline(yintercept = .12, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .47, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 1.65, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) +
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Cadmium (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

cd.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Cd),
            Max = max(Cd),
            Min = min (Cd),
            SD = sd(Cd),
            Med = median(Cd))

## anova
aov <- aov(Cd ~ Trophic.Group + Location, data = t.metals2)
summary(aov)
plot(aov)
#
leveneTest(Cd ~ Trophic.Group, data = t.metals2)
kruskal.test(Cd ~ Trophic.Group, data = t.metals2)
dunnTest(Cd ~ Trophic.Group, data = t.metals2, method="bonferroni")
plot(Cd ~ Trophic.Group, data = t.metals2)
#
leveneTest(Cd ~ Location, data = t.metals2)
aov2 <- aov(Cd ~ Location, data = t.metals2)
summary(aov2)

### chromium
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Cr, 
           color = Trophic.Group, shape = Cr_Flag)) +
  geom_hline(yintercept = 5.5, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic")) +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Total Chromium (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors)  +
  guides(shape = FALSE) 

cr.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Cr),
            Max = max(Cr),
            Min = min (Cr),
            SD = sd(Cr),
            Med = median(Cr))

## anova
aov <- aov(Cr ~ Trophic.Group + Location, data = t.metals2)
summary(aov)
plot(aov)
#
leveneTest(Cr ~ Trophic.Group, data = t.metals2)
kruskal.test(Cr ~ Trophic.Group, data = t.metals2)
dunnTest(Cr ~ Trophic.Group, data = t.metals2, method="bonferroni")
plot(Cr ~ Trophic.Group, data = t.metals2)
#
leveneTest(Cr ~ Location, data = t.metals2)
aov <- aov(Cr ~ Location, data = t.metals2)
summary(aov)


### copper 
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Cu, 
           color = Trophic.Group, shape = Cu_Flag)) + 
  geom_hline(yintercept = 2.35, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 9.41, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 32.94, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Copper (ug/g wet weight)",
       color = "Trophic Group:") + 
  scale_color_manual(values = trophic_colors) +
  guides(shape = FALSE) + scale_y_log10() 

cu.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Cu),
            Max = max(Cu),
            Min = min (Cu),
            SD = sd(Cu),
            Med = median(Cu))

## anova
aov <- aov(Cu ~ Trophic.Group + Location, data = t.metals2)
summary(aov)
plot(aov)
#
leveneTest(Cu ~ Trophic.Group, data = t.metals2) 
kruskal.test(Cu ~ Trophic.Group, data = t.metals2) 
dunnTest(Cu ~ Trophic.Group, data = t.metals2, method="bonferroni")
plot(Cu ~ Trophic.Group, data = t.metals2)
#
leveneTest(Cu ~ Location, data = t.metals2) 
kruskal.test(Cu ~ Location, data = t.metals2) 
dunnTest(Cu ~ Location, data = t.metals2, method="bonferroni")
plot(Cu ~ Location, data = t.metals2)


### manganese
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Mn, 
           color = Trophic.Group, shape = Mn_Flag)) + 
  geom_hline(yintercept = 32.94, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 131.76, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 461.18, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Manganese (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

mn.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Mn),
            Max = max(Mn),
            Min = min (Mn),
            SD = sd(Mn),
            Med = median(Mn))


### nickel
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Ni, 
           color = Trophic.Group, shape = Ni_Flag)) + 
  geom_hline(yintercept = 4.71, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 18.82, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 65.88, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Nickel (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

ni.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Ni),
            Max = max(Ni),
            Min = min (Ni),
            SD = sd(Ni),
            Med = median(Ni))


### lead
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Pb, 
           color = Trophic.Group, shape = Pb_Flag)) + 
  geom_hline(yintercept = .336, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Lead (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

pb.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Pb),
            Max = max(Pb),
            Min = min (Pb),
            SD = sd(Pb),
            Med = median(Pb))

## anova
aov <- aov(Pb ~ Trophic.Group + Location, data = t.metals2)
summary(aov)
plot(aov)
#
leveneTest(Pb ~ Trophic.Group, data = t.metals2) 
kruskal.test(Pb ~ Trophic.Group, data = t.metals2) 
dunnTest(Pb ~ Trophic.Group, data = t.metals2, method="bonferroni")
plot(Pb ~ Trophic.Group, data = t.metals2)
#
leveneTest(Pb ~ Location, data = t.metals2) 
aov <- aov(Pb ~ Location, data = t.metals2)
summary(aov) 
tukey <- TukeyHSD(aov)
tukey


### selenium
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Se, 
           color = Trophic.Group, shape = Se_Flag)) + 
  geom_hline(yintercept = 1.18, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 4.71, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 16.47, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Selenium (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

se.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Se),
            Max = max(Se),
            Min = min (Se),
            SD = sd(Se),
            Med = median(Se))


### silver
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Ag, 
           color = Trophic.Group, shape = Ag_Flag)) + 
  geom_hline(yintercept = 1.18, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 4.71, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 16.47, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Silver (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

ag.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Ag),
            Max = max(Ag),
            Min = min (Ag),
            SD = sd(Ag),
            Med = median(Ag))


### tin
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Sn, 
           color = Trophic.Group, shape = Sn_Flag)) + 
  geom_hline(yintercept = 70.59, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 282.35, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 988.24, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Tin (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

sn.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Sn),
            Max = max(Sn),
            Min = min (Sn),
            SD = sd(Sn),
            Med = median(Sn))


### zinc
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Zn, 
           color = Trophic.Group, shape = Zn_Flag)) + 
  geom_hline(yintercept = 70.59, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 282.35, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 988.24, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Zinc (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

zn.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Zn),
            Max = max(Zn),
            Min = min (Zn),
            SD = sd(Zn),
            Med = median(Zn))


### iron
ggplot(t.metals2, 
       aes(x = Species.Ab,  y = Fe, 
           color = Trophic.Group, shape = Fe_Flag)) + 
  geom_hline(yintercept = 1852.94, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 529.41, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 132.35, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Iron (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1)) +
  guides(shape = FALSE) 

fe.summary <- t.metals2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique(Sample)),
            Mean = mean(Fe),
            Max = max(Fe),
            Min = min (Fe),
            SD = sd(Fe),
            Med = median(Fe))


### HG ####
Hg <- read.csv("./RMI_Hg.csv", header = T)
str(Hg)
t.Hg <- as.data.frame(t(Hg))
names(t.Hg) <- lapply(t.Hg[1,], as.character)
t.Hg <- t.Hg[-1,]
t.Hg2 <- t.Hg[,c(4,5,8)]
t.Hg2$conc <- as.numeric(as.character(t.Hg2$conc))
t.Hg2 <- arrange(t.Hg2, `Sample ID`)
t.Hg2 <- cbind(locations, t.Hg2)
t.Hg2 <- inner_join(t.Hg2, trophic)
t.Hg2$Trophic.Group <- ordered(t.Hg2$Trophic.Group, 
                               levels = c("Detritivore", "Herbivore",
                                          "Benthic Invertivore", "Carnivore", "Pelagic"))
t.Hg2$Location <- as.character(t.Hg2$Location)
t.Hg2[t.Hg2 == "MC Oceanside - Kwaj"] <- "MC Kwaj."
t.Hg2$Location <- ordered(t.Hg2$Location,
                          levels = c("Kwajalein", "MC Kwaj.", "Ebeye",
                                     "Jaluit", "Majuro", "Rongelap", 
                                     "Utirik", "Wotje"))
t.Hg2$Species.Ab <- ordered(t.Hg2$Species.Ab,
                            levels = c("C. microrhinos", "C. sordidus" , 
                                       "H. longiceps", "A. triostegus",
                                       "A. lineatus",
                                       "N. lituratus", "L. gibbus",
                                       "S. spiniferum", "L. bohar",
                                       "P. laevis", "V. louti",
                                       "E. polyphekadion",
                                       "T. albacares", "K. pelamis",
                                       "G. unicolor", "E. bipinnulata"))


ggplot(t.Hg2, 
       aes(x = Species.Ab,  y = conc, 
           color = Trophic.Group)) + 
  geom_hline(yintercept = .5, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 1, col = "darkgrey", linetype='solid', size = .75) +  ##
  geom_jitter(size = 5, alpha = .6) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic")) +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Total Mercury (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) 

hg.summary <- t.Hg2 %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique('Sample ID')),
            Mean = mean(conc),
            Max = max(conc),
            Min = min (conc),
            SD = sd(conc),
            Med = median(conc))

## anova
aov <- aov(conc ~ Trophic.Group + Location, data = t.Hg2)
summary(aov)
plot(aov)
#
leveneTest(conc ~ Trophic.Group, data = t.Hg2) 
kruskal.test(conc ~ Trophic.Group, data = t.Hg2) 
dunnTest(conc ~ Trophic.Group, data = t.Hg2, method="bonferroni")
plot(conc ~ Trophic.Group, data = t.Hg2)
#
leveneTest(conc ~ Location, data = t.Hg2) 
aov <- aov(conc ~ Location, data = t.Hg2)
summary(aov) 


#### PCBs ####
PCB <- read.csv("./RMI_PCB.csv", header=T)
str(PCB)
t.PCB <- as.data.frame(t(PCB))
names(t.PCB) <- lapply(t.PCB[1,], as.character)
t.PCB <- t.PCB[-1,]
t.PCB2 <- t.PCB[t.PCB$Matrix == "Tissue", ]
t.PCB2 <- t.PCB2[,-179]
t.PCB3 <- as.data.frame(t.PCB2[,c(1,3,22:179)])
t.PCB3 <- t.PCB3 %>% 
  mutate_at(vars(c(3:ncol(t.PCB3))), as.character) %>% 
  mutate_at(vars(c(3:ncol(t.PCB3))), as.numeric)

t.PCB3.melt <- melt(t.PCB3, id.vars = c(1:2),
                    measure.vars = c(3:ncol(t.PCB3)),
                    variable.name = "Compound")

PCB.ug.g <- t.PCB3[,c(3:ncol(t.PCB3))]
PCB.ug.g[,] <- PCB.ug.g[,]/1000
PCB.ug.g <- cbind(t.PCB3[,1:2], PCB.ug.g)
PCB.ug.g <- arrange(PCB.ug.g, `Sample ID`)
PCB.ug.g <- cbind(locations, PCB.ug.g)
PCB.ug.g <- inner_join(PCB.ug.g, trophic)
PCB.ug.g$Trophic.Group <- ordered(PCB.ug.g$Trophic.Group,
                                  levels = c("Detritivore", "Herbivore",
                                  "Benthic Invertivore", "Carnivore", "Pelagic"))
PCB.ug.g$PCB20 <- PCB.ug.g$`PCB 8/5` + PCB.ug.g$`PCB 18` + 
  PCB.ug.g$`PCB 28/31` + PCB.ug.g$`PCB 44` +
  PCB.ug.g$`PCB 52` + PCB.ug.g$`PCB 66/80` + 
  PCB.ug.g$`PCB 77` + PCB.ug.g$`PCB 101/84/90` +
  PCB.ug.g$`PCB 105/127` + PCB.ug.g$`PCB 118/108` +
  PCB.ug.g$`PCB 126` + PCB.ug.g$`PCB 128/167` + 
  PCB.ug.g$`PCB 138/164/163` + PCB.ug.g$`PCB 153/168` + 
  PCB.ug.g$`PCB 170/190` + PCB.ug.g$`PCB 180/193` + 
  PCB.ug.g$`PCB 187/182` + PCB.ug.g$`PCB 195` + 
  PCB.ug.g$`PCB 206` + PCB.ug.g$`PCB 209`
PCB.ng.g <- PCB.ug.g[,4:161]*1000
PCB.ng.g$PCB20 <- PCB.ug.g$PCB20*1000
PCB.ng.g <- cbind(PCB.ug.g[,c(1:3,162),], PCB.ng.g)

PCB.ug.g$Trophic.Group <- ordered(PCB.ug.g$Trophic.Group, 
                                  levels = c("Detritivore", "Herbivore",
                                             "Benthic Invertivore", "Carnivore", "Pelagic"))
PCB.ug.g$Location <- as.character(PCB.ug.g$Location)
PCB.ug.g[PCB.ug.g == "MC Oceanside - Kwaj"] <- "MC Kwaj."
PCB.ug.g$Location <- ordered(PCB.ug.g$Location,
                             levels = c("Kwajalein", "MC Kwaj.", "Ebeye",
                                        "Jaluit", "Majuro", "Rongelap", 
                                        "Utirik", "Wotje"))

pcb.summary <- PCB.ug.g %>% 
  group_by(Location, Trophic.Group, Species.Ab) %>% 
  summarise(N = length(unique('Sample ID')),
            Mean = mean(`Total PCB`),
            Max = max(`Total PCB`),
            Min = min (`Total PCB`),
            SD = sd(`Total PCB`),
            Med = median(`Total PCB`))

ggplot(PCB.ug.g, 
       aes(x = Species.Ab,  y = `Total PCB`, 
           color = Trophic.Group)) + 
  geom_hline(yintercept = .0025, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_jitter(size = 4, alpha = .7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Total PCBs (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors)

## anova
aov <- aov(`Total PCB` ~ Trophic.Group + Location, data = PCB.ug.g)
summary(aov)
plot(aov)
#
leveneTest(`Total PCB` ~ Trophic.Group, data = PCB.ug.g) 
kruskal.test(`Total PCB` ~ Trophic.Group, data = PCB.ug.g) 
dunnTest(`Total PCB` ~ Trophic.Group, data = PCB.ug.g, method="bonferroni")
plot(`Total PCB` ~ Trophic.Group, data = PCB.ug.g)
#
leveneTest(`Total PCB` ~ Location, data = PCB.ug.g) 
kruskal.test(`Total PCB` ~ Location, data = PCB.ug.g) 
dunnTest(`Total PCB` ~ Location, data = PCB.ug.g, method="bonferroni")
plot(`Total PCB` ~ Location, data = PCB.ug.g)


#### PAHS ####
PAH <- read.csv("./RMI_PAH.csv", header = T)
t.PAH <- as.data.frame(t(PAH))
names(t.PAH) <- lapply(t.PAH[1,], as.character)
t.PAH <- t.PAH[-1,]
t.PAH2 <- t.PAH[t.PAH$Matrix == "Tissue", ]
t.PAH2 <- t.PAH2[,-86]
t.PAH3 <- as.data.frame(t.PAH2[,c(1,3,22:86,93:119)])
t.PAH3 <- t.PAH3 %>% 
  mutate_at(vars(c(3:ncol(t.PAH3))), as.character) %>% 
  mutate_at(vars(c(3:ncol(t.PAH3))), as.numeric)
t.PAH3.melt <- melt(t.PAH3, id.vars = c(1:2),
                    measure.vars = c(3:ncol(t.PAH3)),
                    variable.name = "Compound")
PAH.ug.g <- t.PAH3[,c(3:ncol(t.PAH3))]
PAH.ug.g <- PAH.ug.g[,]/1000
PAH.ug.g <- cbind(t.PAH3[,1:2], PAH.ug.g)
PAH.ug.g$PAH4 <- PAH.ug.g$`Benz(a)anthracene` + PAH.ug.g$`Chrysene/Triphenylene` +
  PAH.ug.g$`C1-Chrysenes` + PAH.ug.g$`C2-Chrysenes` + PAH.ug.g$`C3-Chrysenes`+
  PAH.ug.g$`C4-Chrysenes` +
  PAH.ug.g$`Benzo(a)fluoranthene` + PAH.ug.g$`Benzo(a)pyrene`
PAH.ug.g <- arrange(PAH.ug.g, `Sample ID`)
PAH.ug.g <- cbind(locations, PAH.ug.g)
PAH.ug.g <- inner_join(PAH.ug.g, trophic)

PAH.summ <- PAH.ug.g %>% 
  group_by(Location, Species) %>% 
  summarise(N = length(unique(`Sample ID`)),
            Mean = mean(`Total PAHs`),
            Max = max(`Total PAHs`),
            Min = min(`Total PAHs`),
            SD = sd(`Total PAHs`),
            Median = median(`Total PAHs`))
PAH.summ <- inner_join(PAH.summ, trophic)

## anova
aov <- aov(`Total PAHs` ~ Trophic.Group + Location, data = PAH.ug.g)
summary(aov)
plot(aov)
#
leveneTest(`Total PAHs` ~ Trophic.Group, data = PAH.ug.g) 
#
leveneTest(`Total PAHs` ~ Location, data = PAH.ug.g) 


PAH.ug.g$Trophic.Group <- ordered(PAH.ug.g$Trophic.Group, levels = c("Detritivore", "Herbivore",
                                                                     "Benthic Invertivore", "Carnivore", "Pelagic"))

perc.napthalene <- PAH.ug.g$Naphthalene/PAH.ug.g$`Total PAHs`*100
perc.napthalene <- cbind(PAH.ug.g[,c(1,2,8,67)], perc.napthalene)
summary(perc.napthalene$perc.napthalene)
length(perc.napthalene[which(perc.napthalene$perc.napthalene > 0),])

PAH.ug.g$Location <- as.character(PAH.ug.g$Location)
PAH.ug.g[PAH.ug.g == "MC Oceanside - Kwaj"] <- "MC Kwaj."
PAH.ug.g$Location <- ordered(PAH.ug.g$Location,
                             levels = c("Kwajalein", "MC Kwaj.", "Ebeye",
                                        "Jaluit", "Majuro", "Rongelap", 
                                        "Utirik", "Wotje"))
ggplot(PAH.ug.g, 
       aes(x = Species.Ab,  y = `Total PAHs`, 
           color = Trophic.Group)) + 
  geom_jitter(size = 4, alpha = .7) + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Total PAHs (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) 

ggplot(PAH.ug.g, 
       aes(x = Species.Ab,  y = Naphthalene, 
           color = Trophic.Group)) + 
  geom_jitter(size = 4, alpha = .7) + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Naphthalene (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) 


#### OC pesticides #### 
OC <- read.csv("./RMI_OC.csv", header = T)
t.OC <- as.data.frame(t(OC))
names(t.OC) <- lapply(t.OC[1,], as.character)
t.OC <- t.OC[-1,]
t.OC2 <- t.OC[seq(1, 271, 2),]
t.OC_flag <- t.OC[seq(2, 272, 2), ]
colnames(t.OC_flag) <- paste(colnames(t.OC_flag), "Flag", sep = "_")
OC_all <- cbind(t.OC2, t.OC_flag)
OC_all2 <- OC_all[,c(1:58)] %>% 
  dplyr::mutate_at(vars(c(22:58)), as.character) %>% 
  dplyr::mutate_at(vars(c(22:58)), as.numeric)
OC_all2 <- cbind(OC_all2, OC_all[,c(251:287)])
OC.ug.g <- OC_all2[,c(22:58)]
OC.ug.g <- OC.ug.g[,]/1000
OC.ug.g <- cbind(OC_all2[,c(1:21)], OC.ug.g, OC_all[,c(251:287)])

OC.ug.g <- arrange(OC.ug.g, `Sample ID`)
OC.ug.g <- cbind(locations, OC.ug.g)
OC.ug.g <- inner_join(OC.ug.g, trophic)
OC.ug.g$Trophic.Group <- as.character(OC.ug.g$Trophic.Group)
OC.ug.g$Trophic.Group <- ordered(OC.ug.g$Trophic.Group, 
                                 levels = c("Detritivore", "Herbivore",
                                            "Benthic Invertivore", "Carnivore", "Pelagic"))
OC.ug.g$Location <- as.character(OC.ug.g$Location)
OC.ug.g[OC.ug.g == "MC Oceanside - Kwaj"] <- "MC Kwaj."
OC.ug.g$Location <- ordered(OC.ug.g$Location,
                            levels = c("Kwajalein", "MC Kwaj.", "Ebeye",
                                       "Jaluit", "Majuro", "Rongelap", 
                                       "Utirik", "Wotje"))
OC.ug.g$Species.Ab <- as.character(OC.ug.g$Species.Ab)
unique(OC.ug.g$Species.Ab)
OC.ug.g$Species.Ab <- ordered(OC.ug.g$Species.Ab,
                              levels = c("C. microrhinos", "C. sordidus" , 
                                         "H. longiceps", "A. triostegus",
                                         "A. lineatus",
                                         "N. lituratus", "L. gibbus",
                                         "S. spiniferum", "L. bohar",
                                         "P. laevis", "V. louti",
                                         "E. polyphekadion",
                                         "T. albacares", "K. pelamis",
                                         "G. unicolor", "E. bipinnulata"))

OC.ug.g$DDT.D.E <- (OC.ug.g$`4,4'-DDE` + OC.ug.g$`2,4'-DDE` + 
                      OC.ug.g$`4,4'-DDD` + OC.ug.g$`2,4'-DDD` + 
                      OC.ug.g$`4,4'-DDT` + OC.ug.g$`2,4'-DDT`)

ggplot(OC.ug.g, 
       aes(x = Species.Ab,  y = DDT.D.E, 
           color = Trophic.Group)) + 
  geom_hline(yintercept = .02, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .09, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .33, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_point(size = 4, alpha = .7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Sum of DDT, DDD, & DDE (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) 

ggplot(OC.ug.g, 
       aes(x = Species.Ab,  y = Mirex, 
           color = Trophic.Group, shape = Mirex_Flag)) + 
  geom_hline(yintercept = .07, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .28, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .98, col = "darkgrey", linetype='solid', size = .75) +  ##  
  geom_point(size = 4, alpha = .7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Mirex (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1,4)) +
  guides(shape = FALSE)

ggplot(OC.ug.g, 
       aes(x = Species.Ab,  y = Chlordecone, 
           color = Trophic.Group, shape = Chlordecone_Flag)) + 
  geom_hline(yintercept = .2, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .85, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = 2.96, col = "darkgrey", linetype='solid', size = .75) +  ##   
  geom_point(size = 4, alpha = .7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Chlordecone (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1,4)) +
  guides(shape = FALSE)

ggplot(OC.ug.g, 
       aes(x = Species.Ab,  y = Hexachlorobenzene, 
           color = Trophic.Group, shape = Hexachlorobenzene_Flag)) + 
  geom_hline(yintercept = .23, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .07, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .02, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_point(size = 4, alpha = .7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Hexachlorobenzene (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,1,4)) +
  guides(shape = FALSE)

summary(OC.ug.g$Chlorpyrifos)
ggplot(OC.ug.g, 
       aes(x = Species.Ab,  y = Chlorpyrifos, 
           color = Trophic.Group, shape = Chlorpyrifos_Flag)) + 
  geom_hline(yintercept = 2.82, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .94, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .23, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_point(size = 4, alpha = .7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Chlorpyrifos (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,4)) +
  guides(shape = FALSE)

ggplot(OC.ug.g, 
       aes(x = Species.Ab,  y = `Beta-HCH`, 
           color = Trophic.Group, shape = `Beta-HCH_Flag`)) + 
  geom_hline(yintercept = 1.98, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .56, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_hline(yintercept = .14, col = "darkgrey", linetype='solid', size = .75) +  ## 
  geom_point(size = 4, alpha = .7) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, 
                                   size = 7, face = "italic"),
        legend.justification = "bottom") +
  facet_grid(~ Location, scales = "free") +
  labs(x = "" , y = "Beta-HCH (ug/g wet weight)",
       color = "Trophic Group:") + scale_y_log10() +
  scale_color_manual(values = trophic_colors) +
  scale_shape_manual(values = c(16,4)) +
  guides(shape = FALSE)





