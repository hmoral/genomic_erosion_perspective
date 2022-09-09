####### load in the packages####
library(plink)
library(plyr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)
library(magrittr)
library(pivottabler)
library(openxlsx)
library(ggpattern)
library(ggVennDiagram)
library(arsenal)
library(ggpubr)
library(cowplot)
library(see)
library(reshape)
library(GGally)
library(multcomp)
library(cowplot)
####### import the UCE CADD scores from the different species ######

setwd("~/Documents/PhD-Conservation-genomic-endangered-birds/Genomes/CADD-multi-species-corrolation/CADD_inputs/")
Pig_CADD <- read.delim("Pig_all_chr_UCE_shared_CADD_head.bed")
Chicken_CADD <- read.delim("Chicken_CADD_scores_multi-sp-share_head.bed")
Human_CADD <- read.delim(("Human_all_chr_UCE_shared_CADD_head.bed"))

Pig_CADD_dup_rmv <- Pig_CADD %>% distinct()
Chicken_CADD_dup_rmv <- Chicken_CADD %>% distinct()
Human_CADD_dup_rmv <- Human_CADD %>% distinct()

###### select the cols that are needed ###
Pig_scores <- Pig_CADD_dup_rmv %>% dplyr::select(UCE_ID ,X.Chrom, Pos, PHRED)
Chicken_scores <- dplyr::select(Chicken_CADD_dup_rmv,UCE_ID ,X.Chrom, Pos, PHRED)
Human_scores <- dplyr::select(Human_CADD_dup_rmv,UCE_ID ,X.Chrom, Pos, PHRED)

Human_scores["Species"] <- c("Human")
Pig_scores["Species"] <- c("Pig")
Chicken_scores["Species"] <- c("Chicken")

max(Chicken_CADD$PHRED)

multi_sp_CADD <- rbind(Pig_scores,Chicken_scores,Human_scores)

UCE_agg <- aggregate(x = multi_sp_CADD$PHRED ,          
                     by = list(multi_sp_CADD$Species, multi_sp_CADD$UCE_ID),
                     FUN = sum) 

colnames(UCE_agg) <-  c("Species","UCE_ID","SUM")

UCE_agg_ln <- join(UCE_lengths,UCE_agg)

UCE_agg_ln["Average"] <- UCE_agg_ln$SUM/ UCE_agg_ln$n

UCE_avg <- select(UCE_agg_ln, UCE_ID, Species, Average)

UCE_spread <- tidyr::spread(UCE_avg,Species,Average)

###########################
######## corrr plot #######
###########################
##### functions############
###########################

library(GGally)
install.packages("GGally")

color_cor <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme(text = element_text(size = 60),panel.background = element_rect(fill=fill, colour=NA),
          panel.grid.major = element_blank()) 
}

color_smooth <- function (data, mapping, ..., method = "lm", formula = y ~ x, 
                          se = TRUE, shrink = TRUE) {
  p <- ggplot(data = data, mapping)
  p <- p + geom_point()  
  if (!is.null(mapping$color) || !is.null(mapping$colour)) {
    p <- p + geom_smooth(method = method, se = se, formula = formula)
  }
  else {
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    corr <- cor(x, y, method="p", use="pairwise")
    p <- p + geom_smooth(method = 'lm', se = FALSE, formula = formula, 
                         colour = I(cbPalette[4])) 
  }
  if (isTRUE(shrink)) {
    p <- p + coord_cartesian(ylim = range(eval_data_col(data, 
                                                        mapping$y), na.rm = TRUE))
  }
  p
}

####plot #####
ms_coorr <- ggpairs(data = UCE_spread, columns = 2:4, , 
                    upper = list(continuous = wrap("cor", size = 8)),
                    lower = list(continuous = color_smooth),
                    diag = list(continuous = "densityDiag" ))+ 
  theme(text = element_text(size = 50),panel.grid = element_blank(), panel.grid.major = element_blank())

ms_coorr 

##### make a small test df with a few UCE ####

multi_sp_CADD$UCE_ID %>% unique() %>% head(100)

multi_sp_CADD %>% distinct() %>% count()

multi_sp_CADD %>% count()

subset_of_UCE <-  multi_sp_CADD %>% filter(UCE_ID == c("uce_1004"))

######## colour blind pallet ######

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot(cbPalette)

######### Pic CADD scores ######
min(filter(subset_of_UCE, Species == c("Pig"))$Pos)  
max(filter(subset_of_UCE, Species == c("Pig"))$Pos)  

pig <- data_frame(x1=c(-1*min(filter(subset_of_UCE, Species == c("Pig"))$Pos) ,-67032568,-67032688), x2=c(-67032568,-67032688,-1*max(filter(subset_of_UCE, Species == c("Pig"))$Pos)),y1=c(0,0,0), y2=c(50,50,50), type=c('flank','UCE','flank'))

Pig <- ggplot(filter(subset_of_UCE, Species == c("Pig"))) + 
  geom_rect(pig,mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill= type), color="black", alpha=0.5) +
  geom_point(aes(x= -1*Pos, y = PHRED, colour = Species),colour = cbPalette[3]) + 
  scale_y_continuous(limits = c(0,50),name = "CADD-score")  + 
  scale_x_continuous(name = "Position",breaks = seq(-67033500, -67031500, by = 500), labels = c("670335","670330","670325","670320","670315")) + 
  scale_fill_manual(values = alpha(c("grey20", "grey"), .5)) +
  theme_classic2(base_size = 20) +
  theme(legend.position="none")
Pig

###### human CADD ######
min(filter(subset_of_UCE, Species == c("Human"))$Pos)  
max(filter(subset_of_UCE, Species == c("Human"))$Pos)  

human <- data_frame(x1=c(min(filter(subset_of_UCE, Species == c("Human"))$Pos) ,32695128,32695248), x2=c(32695128,32695248,max(filter(subset_of_UCE, Species == c("Human"))$Pos)),y1=c(0,0,0), y2=c(50,50,50), type=c('flank','UCE','flank'))

Human <- ggplot(filter(subset_of_UCE, Species == c("Human"))) + 
  geom_rect(human,mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill= type), color="black", alpha=0.5) +
  geom_point(aes(x= Pos, y = PHRED, colour = Species),colour = cbPalette[4]) + 
  scale_x_continuous(name = "Position", limits = c(32694128,32696248),breaks = seq(32694500, 32696000, by = 500), labels = c("326945","326950","326955","329660")) + 
  scale_y_continuous(limits = c(0,50),name = "CADD-score") + 
  scale_fill_manual(values = alpha(c("grey20", "grey"), .5)) +
  theme_classic2(base_size = 20)+
  theme(legend.position="none") 
Human

#### chicken UCE-1004 35164016-35164136 CADD ######
min(filter(subset_of_UCE, Species == c("Chicken"))$Pos)  
max(filter(subset_of_UCE, Species == c("Chicken"))$Pos)  

chicken <- data_frame(x1=c(35163000,35164016,35164136), x2=c(35164016,35164136,35165136),y1=c(0,0,0), y2=c(50,50,50), type=c('flank','UCE','flank'))

Chicken <- ggplot(filter(subset_of_UCE, Species == c("Chicken"))) + 
  geom_rect(chicken,mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill= type), color="black", alpha=0.5) +
  geom_point(aes(x= Pos, y = PHRED, colour = Species, line = "black"),colour = cbPalette[2]) + 
  scale_y_continuous(limits = c(0,50),name = "CADD-score") + 
  scale_x_continuous(name = "Position", labels = c("351630","351635","351640","351645","351650"))+
  scale_fill_manual(values = alpha(c("grey20", "grey"), .5)) +
  theme_classic2(base_size = 20) +
  theme(legend.position="none")

######combi plot #####

per_base_CADD_no_pic <- ggarrange(Chicken,Human,Pig, nrow = 3, ncol = 1, labels = c("Chicken","Human","Pig"), label.y = 1.1, label.x = c(-.01,0,0.04)) + theme(legend.position="none") 

per_base_CADD

############################################
#### plotting DNA seq for uce-1004 #########
############################################

###### import the sequences
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("YuLab-SMU/ggmsa")

install.packages("Biostrings")

library(ggmsa)
library(Biostrings)
library(ggplot2)

nt_sequences <- readDNAMultipleAlignment(c("~/Documents/PhD-Conservation-genomic-endangered-birds/Genomes/CADD-multi-species-corrolation/UCE_1004_flangs_align.msa"))

setwd("~/Documents/PhD-Conservation-genomic-endangered-birds/Genomes/CADD-multi-species-corrolation/")
setwd("~/Downloads/")
nt_sequences <- readDNAMultipleAlignment("prank-S20220411-111850-0641-43206619-p2m.fas")

nt_sequences <- readDNAMultipleAlignment(c("~/Documents/PhD-Conservation-genomic-endangered-birds/Genomes/CADD-multi-species-corrolation/UCE-1004_all3_sp.fa"))

complete_msa <- tidy_msa(nt_sequences, 0,2578)

############################
########plotting a consensus ############
###### how to plot just the geom_ggmsaBar aspect of this plot:

msa_plot <- ggplot() + geom_msa(complete_msa) + geom_msaBar()

consensus_plot <- msa_plot$plotlist[[1]]

######needs the top part as well i presume 
###### so can I crop it out ?

cropped_consensus_plot <- consensus_plot + theme(plot.margin = unit(c (+0.1,+0.1,-0.4,-0.3), "cm"))

######## convert the corr plot for ggarrange
grid_no_pic <- plot_grid(
  ggmatrix_gtable(ms_coorr),
  nrow = 1
)

########### add white space to the cropped_consensus

cropped_consensus_plot_white_space <- ggarrange("", cropped_consensus_plot ,"", ncol = 3, nrow = 1, widths = c(0.15,1,0.05),
                                                common.legend = TRUE)
########## put the consensus above the per base CADD

per_base_CADD_no_pic_full <- ggarrange(cropped_consensus_plot_white_space,per_base_CADD_no_pic,  ncol = 1, nrow = 2,heights = c(0.25,2),
                                       common.legend = TRUE)
##### plot 2 with no icons 
##### this is the one used in the paper
ggarrange(per_base_CADD_no_pic_full,grid_no_pic,  ncol = 2, nrow = 1, widths = c(1,1),labels = c("A.","B."),
          common.legend = TRUE)

############