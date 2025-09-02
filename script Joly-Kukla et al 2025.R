R.version
sessionInfo()

# Install the packages if needed

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("lefser")

BiocManager::install("phyloseq")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
BiocManager::install("ggtree")
BiocManager::install("readr")
BiocManager::install("phytools")

install.packages("writexl")
install.packages("readxl")
install.packages("moonBook")
install.packages("webr")
install.packages("paletteer")
install.packages("emmeans")
install.packages("pscl")
install.packages("tidyverse")
install.packages("devtools")
install.packages("usethis")


# read the libraries o the packages twice
library(devtools)
library(usethis)
library(dplyr)
library(phyloseq) 
library(readr)
library(ggplot2) 
library(writexl)
library(readxl)
library(ape)
library(stringr)
library(paletteer)
library(emmeans)
library(patchwork)
library(webr)
library(tidyr)
library(ggtree) 
library(vegan)
library(phytools)
library(stringr)
library(pscl)
library(car)
library(lefser)
library(tidyverse)
library(microbiome)
library(rstatix)

# Import data and format ####
## biom ####
seq3biom<-import_biom(BIOMfilename="biom/sq3biomfile.biom1")
seq3biom

## Metadata ####
metadataseq3 <- read_delim("biom/metadataseq3.csv", ";", escape_double = FALSE,
                           trim_ws = TRUE, show_col_types = FALSE)
# change into a dataframe
metadataseq3 <- as.data.frame(metadataseq3)
# Set the rownames with the "Samples" column
row.names(metadataseq3) = metadataseq3$echantillons
# complete your phyloseq object with the sample_data() function
seq3biom@sam_data=sample_data(metadataseq3)
#	seq3biom  is now a complete phyloseq object ! let's explore it.
seq3biom
# use the colnames and c() functions to change the tax colnames
colnames(seq3biom@tax_table)=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

## Subsets ####
# Subsetting the phyloseq object
# We created a subset of the BIOM file to include only the samples relevant to the present article.
cut_subset_seq3=subset_samples(seq3biom, lab_process =="decoupe")
cut_subset_seq3

temporel_subset_seq3=subset_samples(seq3biom, analyse_type =="temporel2022")
temporel_subset_seq3

## Update of genera names ####
tax_temp <- temporel_subset_seq3@tax_table
tax_temp <- as.data.frame(tax_temp)
write_xlsx(tax_temp, path = "dataframes/tax_temp.xlsx")

tax_temp_ok <- read_xlsx("dataframes/tax_temp_ok.xlsx")

tax_temp_ok <- tax_temp_ok %>% 
  tibble::column_to_rownames("otu")

tax_mat <- as.matrix(tax_temp_ok)

temporel_subset_seq3@tax_table=tax_table(tax_mat)
temporel_subset_seq3@tax_table

## Tree (Fig S1) #### 
tree <- ape::read.tree("biom/phyl_tree.nhx")
plot(tree)

temporel_subset_seq3@phy_tree=phy_tree(tree)
temporel_subset_seq3

tax_table <- temporel_subset_seq3@tax_table
tax_table <-as.data.frame(tax_table)

write.table(tax_table, "dataframes/tax_table1.txt", row.names = T)
tax_table2 <- read.delim("dataframes/tax_table2.txt")
tax_table2
# we updated the names of the multi-affiliated ASVs for more clarity. For example, we renamed
# the multi-affiliated ASV from the Moraxellaceae family : "M.affASVMoraxellaceaefamily" 

# combined name (ASV and genera)
for (i in 1:length(tree$tip.label)) {
  cluster_name <- tree$tip.label[i]
  genus_name <- tax_table2$Genus[which(tax_table2$Cluster == cluster_name)]
  
  if (length(genus_name) > 0) {
    tree$tip.label[i] <- paste(cluster_name, genus_name, sep = "_")
  }
}
write.tree(tree, file = "biom/updated_phyl_tree.nhx")
tree2 <- ape::read.tree("biom/updated_phyl_tree.nhx")

p <- ggtree(tree2, layout = "circular")
tree_data <- p$data

# extraction of the genera
tree_data <- tree_data %>%
  mutate(genus = str_extract(label, "[^_]+$"))

# create a dataframe from tax table with genus and phylum
genus_phylum <- tax_table2 %>%
  select(Genus, Phylum) %>%
  distinct()

# add the phylum column to the tree data
tree_data <- left_join(tree_data, genus_phylum, by = c("genus" = "Genus"))

tree2$tip.label <- sub("^Cluster_", "ASV_", tree2$tip.label)

png("figures/Figure S1 .png", units = 'in', width = 5.37, height = 4.79, res = 300) 
ggtree(tree2, layout = "circular", branch.length = "none") %<+% tree_data +
  geom_tree() +
  geom_tiplab2(aes(angle = angle, label = label, color = Phylum), size = 1.5, offset = 2, hjust = 0.1) +
  theme_tree() +
  theme(
    plot.margin = margin(40, 10, 90, 20),
    legend.position = "none"
  )
dev.off()

## Rel abund transofmation####
temp_rel = transform_sample_counts(temporel_subset_seq3, function(x) x/sum(x)*100)

# . ####
# Run the functions that will be needed ####

# summarySE

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# PieDonutCustom function
PieDonutCustom <- function (data, mapping, start = getOption("PieDonut.start", 
                                                             0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, 
                            showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 
                                                                                                     0.02), labelposition = getOption("PieDonut.labelposition", 
                                                                                                                                      2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0", 
                                                                                                                                                                                       0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2", 
                                                                                                                                                                                                                                              1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                            color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                            showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                            pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                            explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
                            family = getOption("PieDonut.family", ""), palette_name="Dark2")
{
  (cols = colnames(data))
  if (use.labels) 
    data = moonBook::addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- moonBook::getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = moonBook::getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = moonBook::getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = moonBook::getMapping(mapping, "x"))
  (donuts = moonBook::getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = moonBook::getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = moonBook::getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
                                                             "\n(", scales::percent(df$ratio), ")"), 
                      as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  mainCol = RColorBrewer::brewer.pal(nrow(df), name=palette_name)
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, "(", df3$ratio, 
                           ")")
      else df3$label = paste0(df3$label, "\n(", df3$ratio, 
                              ")")
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * 
                                                                     pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + ggforce::theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                             r0 = as.character(r0), r = as.character(r1), start = "start1", 
                                             end = "end1", fill = pies), alpha = pieAlpha, color = color, 
                                  data = df) + transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                       y = "segy", xend = "segxend", yend = "segyend"), 
                            data = df) + geom_text(aes_string(x = "segxend", 
                                                              y = "segyend", label = "label", hjust = "hjust", 
                                                              vjust = "vjust"), size = pieLabelSize, data = df, 
                                                   family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                       y = "segy", xend = "segxend", yend = "segyend"), 
                            data = df[df$ratio < labelpositionThreshold, ]) + 
      geom_text(aes_string(x = "segxend", y = "segyend", 
                           label = "label", hjust = "hjust", 
                           vjust = "vjust"), size = pieLabelSize, 
                data = df[df$ratio < labelpositionThreshold, 
                ], family = family) + geom_text(aes_string(x = "labelx", 
                                                           y = "labely", label = "label"), size = pieLabelSize, 
                                                data = df[df$ratio >= labelpositionThreshold, ], 
                                                family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                    label = "label"), size = pieLabelSize, data = df, 
                         family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                 start = "start1", end = "end1", fill = "no", 
                                                 explode = "focus"), alpha = donutAlpha, 
                                      color = color, data = df3)
    }
    else {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                 start = "start1", end = "end1", fill = "no"), 
                                      alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                         y = "segy", xend = "segxend", yend = "segyend"), 
                              data = df3) + geom_text(aes_string(x = "segxend", 
                                                                 y = "segyend", label = "label", hjust = "hjust", 
                                                                 vjust = "vjust"), size = donutLabelSize, 
                                                      data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", 
                                      y = "labely", label = "label"), size = donutLabelSize, 
                           data = df3, family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                         y = "segy", xend = "segxend", yend = "segyend"), 
                              data = df3[df3$ratio1 < labelpositionThreshold, 
                              ]) + geom_text(aes_string(x = "segxend", 
                                                        y = "segyend", label = "label", hjust = "hjust", 
                                                        vjust = "vjust"), size = donutLabelSize, 
                                             data = df3[df3$ratio1 < labelpositionThreshold, 
                                             ], family = family) + geom_text(aes_string(x = "labelx", 
                                                                                        y = "labely", label = "label"), size = donutLabelSize, 
                                                                             data = df3[df3$ratio1 >= labelpositionThreshold, 
                                                                             ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, 
                          label = title, size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, 
                          y = r3, label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    grid::grid.newpage()
    print(p1, vp = grid::viewport(height = 1, width = 1))
    print(p3, vp = grid::viewport(height = 1, width = 1))
  }
  else {
    p1
  }
}

# . ####
# Tick sex analysis ####
#subset of male and female ticks
temporel_females=subset_samples(temporel_subset_seq3, tick_sex =="F")

temporel_males=subset_samples(temporel_subset_seq3, tick_sex =="M")

# total number of sequences in males and females 
tp_nseq_f <- sum(taxa_sums(temporel_females))
tp_nseq_f # 227373

tp_nseq_m <- sum(taxa_sums(temporel_males))
tp_nseq_m # 229204
# approximately same number of sequences between males and females

# number of endosymbiotic sequences
tax_f <- tax_table(temporel_females)
midi_seqf <- taxa_names(temporel_females)[tax_f[, "Genus"] == "Midichloria"]
number_f <- sum(taxa_sums(temporel_females)[midi_seqf])
number_f
# 117988

tax_m <- tax_table(temporel_males)
midi_seqm <- taxa_names(temporel_males)[tax_m[, "Genus"] == "Midichloria"]
number_m <- sum(taxa_sums(temporel_males)[midi_seqm])
number_m
# 64012

# 117988 / 64012 = 1. 84. There are about 1.9 x more Midichloria sequences in females than in males

## Pie charts (Fig 1 AB)#### 

# Key question: What are the nine most common genera that comprise the majority of sequences in females and males ?

### Females (F)####

# dataframe for pie chart
t_data_pie_f <- tax_glom(temporel_females, taxrank = 'Genus', NArm = FALSE)
t_data_pie_f2 <- psmelt(t_data_pie_f)
as.data.frame(t_data_pie_f2)
write_xlsx(t_data_pie_f2, path = "dataframes/t_data_pie_f2.xlsx")

# load of pietemp_tot_f
# This dataframe is the result of manual analysis on the t_data_pie_f2 excel using the excel function "pivot table"
# This analysis was performed to calculate the sequences percentages of genera and phyla in order to plot the pie chart
pietemp_tot_f <- read_excel("dataframes/pietemp_tot_f.xlsx")
pie_temp_tot_f <- as.data.frame(pietemp_tot_f)
pie_temp_tot_f = pietemp_tot_f %>% group_by(Phylum, Genus) %>% summarise(Percentage = sum(Percentage_genus))
#pie_temp_tot_f$Phylum <- as.factor(pie_temp_tot_f$Phylum)
#pie_temp_tot_f$Genus <- as.factor(pie_temp_tot_f$Genus)
pietemp_tot_f$Percentage <- as.numeric (pietemp_tot_f$Percentage_genus)

"Phylum" %in% names(pie_temp_tot_f)
"Genus" %in% names(pie_temp_tot_f)
"Percentage" %in% names(pie_temp_tot_f)

png("figures/Figure 1A raw.png", units = 'in', width = 6.11, height = 7.87, res = 300) 
PieDonutCustom(pie_temp_tot_f, aes(Phylum, Genus, count=Percentage), 
               ratioByGroup = FALSE,
               showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0),
               maxx = 1.7,
               pieLabelSize = 0.001,
               donutLabelSize = 0.001,
               r0 = 0, r1 = 1, r2 = 1.6,
               labelpositionThreshold = 1,
               showPieName=FALSE
)
dev.off()

### Males (M)####

t_data_pie_m <- tax_glom(temporel_males, taxrank = 'Genus', NArm = FALSE)
t_data_pie_m2 <- psmelt(t_data_pie_m)
as.data.frame(t_data_pie_m2)
write_xlsx(t_data_pie_m2, path = "dataframes/t_data_pie_m2.xlsx")

# load of pietemp_tot_m
# This dataframe is the result of manual analysis on the t_data_pie_m2 excel using the excel function "pivot table"
# This analysis was performed to calculate the sequences percentages of genera and phyla in order to plot the pie chart
pietemp_tot_m <- read_excel("dataframes/pietemp_tot_m.xlsx")
pie_temp_tot_m <- as.data.frame(pietemp_tot_m)
pie_temp_tot_m = pietemp_tot_m %>% group_by(Phylum, Genus) %>% summarise(Percentage = sum(Percentage_genus))
pie_temp_tot_m$Phylum <- as.factor(pie_temp_tot_m$Phylum)
pie_temp_tot_m$Genus <- as.factor(pie_temp_tot_m$Genus)
pietemp_tot_m$Percentage <- as.numeric (pietemp_tot_m$Percentage_genus)

png("figures/Figure 1B raw.png", units = 'in', width = 6.11, height = 7.87, res = 300) 
PieDonutCustom(pie_temp_tot_m, aes(Phylum, Genus, count=Percentage), 
               ratioByGroup = FALSE,
               showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0),
               maxx = 1.7,
               pieLabelSize = 0.001,
               donutLabelSize = 0.001,
               r0 = 0, r1 = 1, r2 = 1.6,
               labelpositionThreshold = 1,
               showPieName=FALSE
)
dev.off()

## Heat map (Fig 1 CD)####

### F ####

# In the female pie chart, there are nine genera for which the number of sequence is above 0,5% of the 
# total sequence number. Therefore, we selected these nine genera to analysis their ubiquity in the samples

# Key question : Are these genera present in all of the individual female ticks ?

temporel_females@sam_data
t_gen_f <- tax_glom(temporel_females, taxrank = 'Genus', NArm = FALSE)
t_heat_f2 <- psmelt(t_gen_f)
t_heat_f2 

t_bin_f <- t_heat_f2 %>%
  mutate(Presence = ifelse(Abundance == 0, "0", "1"))
t_bin_f 

as.data.frame(t_bin_f)

t_bin_f$Genus <- ifelse(
  t_bin_f$Genus == "M.-aff. ASV (Comamonadaceae family)",
  "M.-aff. ASV\n(Comamonadaceae\nfamily)",
  t_bin_f$Genus
)

t_bin_f$Genus <- ifelse(
  t_bin_f$Genus == "M.-aff. ASV (Moraxellaceae family)",
  "M.-aff. ASV\n(Moraxellaceae\nfamily)",
  t_bin_f$Genus
)

t_bin_f$Genus <- ifelse(
  t_bin_f$Genus == "M.-aff. ASV (Bacilli class)",
  "M.-aff. ASV\n(Bacilli class)",
  t_bin_f$Genus
)

unique(t_bin_f$Genus)

order_gen3 <- c('Midichloria',"Francisella","Corynebacterium","Rickettsia","M.-aff. ASV\n(Bacilli class)","Staphylococcus",
                "Mycobacterium","M.-aff. ASV\n(Moraxellaceae\nfamily)","Williamsia")

t_bin_f$Genus <- factor(t_bin_f$Genus, levels = order_gen3)

png("figures/Figure 1C.png", units = 'in', width = 4.93, height = 4.87, res = 300) 
t_bin_f %>%
  dplyr::filter(Genus %in% c("Corynebacterium", "Mycobacterium", "Williamsia", "M.-aff. ASV\n(Bacilli class)", 
                      "Staphylococcus", "Midichloria", "Francisella", "M.-aff. ASV\n(Moraxellaceae\nfamily)", "Rickettsia")) %>%
  dplyr::mutate(Genus = factor(Genus, levels = order_gen3)) %>% 
  ggplot() +
  aes(x = echantillons, y = Genus, fill = Presence) +
  geom_tile() +
  scale_fill_manual(
    values = c(`0` = "#EAEAEA",
               `1` = "#8B8B8B"))+
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black", size = 13,face = "italic"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_blank()
  )+
  labs(y="", x="")
dev.off()

### M ####

t_gen_m <- tax_glom(temporel_males, taxrank = 'Genus', NArm = FALSE)
t_heat_m2 <- psmelt(t_gen_m)

t_bin_m <- t_heat_m2 %>%
  mutate(Presence = ifelse(Abundance == 0, "0", "1"))
t_bin_m 

as.data.frame(t_bin_m)

t_bin_m$Genus <- ifelse(
  t_bin_m$Genus == "M.-aff. ASV (Comamonadaceae family)",
  "M.-aff. ASV\n(Comamonadaceae\nfamily)",
  t_bin_m$Genus
)

t_bin_m$Genus <- ifelse(
  t_bin_m$Genus == "M.-aff. ASV (Moraxellaceae family)",
  "M.-aff. ASV\n(Moraxellaceae\nfamily)",
  t_bin_m$Genus
)

t_bin_m$Genus <- ifelse(
  t_bin_m$Genus == "M.-aff. ASV (Bacilli class)",
  "M.-aff. ASV\n(Bacilli class)",
  t_bin_m$Genus
)

order_gen3 <- c('Midichloria',"Francisella","Corynebacterium","Rickettsia","M.-aff. ASV\n(Bacilli class)","Staphylococcus",
                "Mycobacterium","M.-aff. ASV\n(Moraxellaceae\nfamily)","Williamsia")

t_bin_m$Genus <- factor(t_bin_m$Genus, levels = order_gen3)

png("figures/Figure 1D.png", units = 'in', width = 4.93, height = 4.87, res = 300) 
t_bin_m %>%
  filter(Genus %in% c("Corynebacterium", "Mycobacterium", "Williamsia", "M.-aff. ASV\n(Bacilli class)", 
                      "Staphylococcus", "Midichloria", "Francisella", "M.-aff. ASV\n(Moraxellaceae\nfamily)", "Rickettsia")) %>%
  ggplot() +
  aes(x = echantillons, y = Genus, fill = Presence) +
  geom_tile() +
  scale_fill_manual(
    values = c(`0` = "#EAEAEA",
               `1` = "#8B8B8B"))+
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black", size = 13,face = "italic"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_blank()
  )+
  labs(y="", x="")
dev.off()

## Alpha diversity ####

# Key question: Are the number of observed genera and their abundances influenced by the tick sex ?
# Do females have different number of genera compared to males ? What about their abundances ? 

temporel_gen <- tax_glom(temporel_subset_seq3, taxrank="Genus")
temporel_gen
ntaxa(temporel_subset_seq3); ntaxa(temporel_gen)
# there are 96 ASV and 30 genera

# Compute alpha-diversity indices using the estimate_richness() function
Richness<-estimate_richness(temporel_subset_seq3 ,measures = "Observed")
Shannon<-estimate_richness(temporel_subset_seq3 ,measures = "Shannon")
InvSimpson<-estimate_richness(temporel_subset_seq3 ,measures = "InvSimpson")

# Create a dataframe to summarize all alpha-diversity indices
# using cbind to merge multiple vectors/dataframes and add columns
adiv_temp <- cbind(temporel_subset_seq3@sam_data, 
                   Richness=Richness$Observed,
                   Shannon=Shannon$Shannon,
                   InvSimpson =InvSimpson)
              
adiv_temp

### Values (Table 1) ####

# Key question : What are the values of Richness, Shannon and Inverse Simpson indexes in males and females ?

SUM_rich_tickS <-summarySE(adiv_temp, measurevar = "Richness", groupvars = "tick_sex")
SUM_rich_tickS

SUM_shan_tickS <-summarySE(adiv_temp, measurevar = "Shannon", groupvars = "tick_sex")
SUM_shan_tickS

SUM_invsimp_tickS <-summarySE(adiv_temp, measurevar = "InvSimpson", groupvars = "tick_sex")
SUM_invsimp_tickS

SE_alpha_div <- list(SUM_rich_tickS,SUM_shan_tickS,SUM_invsimp_tickS)

# This file was used for the Table 1

### Statistical analysis (Table 1)####

# Key question : Are the alpha diversity indexes values statistically influenced by the tick sex ?

# Richness

# normality test

shapiro.test(adiv_temp$Richness)
# p-value non-significant so richess follows a normal distribution
hist(adiv_temp$Richness)
adiv_temp

# we evaluate the influence of the tick sex on the richness but we
# have to take into account the month in the model, since this maximal model cannot be reduced
# because the interaction between the variables "month" and "tick-sex" is significant

model_rich <- lm(Richness ~ tick_sex*month, data = adiv_temp)
car::Anova(model_rich, type = "III")
summary(model_rich)

# normality verification for the model
residuals <- residuals(model_rich, type = "pearson")

# residuals
ggplot(data.frame(residuals), aes(x = seq_along(residuals), y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Observations", y = "Residuals (Pearson)")

qqnorm(residuals)
qqline(residuals, col = 2)
shapiro.test(residuals)

# Shannon

#normality test
shapiro.test(adiv_temp$Shannon)
# p-value significant so the Shannon index does not follow a normal distribution
hist(adiv_temp$Shannon)

summary(adiv_temp$Shannon)

model_shannon <- glm(Shannon ~ tick_sex*month, family=Gamma(link = "identity"), data = adiv_temp)
car::Anova(model_shannon, type = "III")
summary(model_shannon)

residuals <- residuals(model_shannon, type = "pearson")
ggplot(data.frame(residuals), aes(x = seq_along(residuals), y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Observations", y = "Residuals (Pearson)")

qqnorm(residuals)
qqline(residuals, col = 2)
shapiro.test(residuals)

# InvSimpson

# normality test

shapiro.test(adiv_temp$InvSimpson)
# p-value significant so InvSimpson does not follow a normal distribution

summary(adiv_temp$InvSimpson)
hist(adiv_temp$InvSimpson)

model_invs <- glm(InvSimpson ~ tick_sex*month, family = Gamma(link = "log"), data = adiv_temp)
car::Anova(model_invs, type = "III")
summary(model_invs)

residuals <- residuals(model_invs, type = "pearson")
ggplot(data.frame(residuals), aes(x = seq_along(residuals), y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Observations", y = "Residuals (Pearson)")

qqnorm(residuals)
qqline(residuals, col = 2)
shapiro.test(residuals)

## Beta diversity####

# Key question : Is the bacterial structure (presence and abundance of taxa) influenced by the tick sex ?
# Different microbiota accordiing to the sex of the ticks ?

# Bray-Curtis distance

### NMDS (Fig 1E) ####

temp.BC2<-ordinate(temp_rel,method = "NMDS",distance = "bray", k=4)
temp.BC2

DF<-plot_ordination(temp_rel, temp.BC2, color = "tick_sex", justDF = TRUE)

palette_sex <- c("F" = "deeppink4", "M" = "cyan4")

png("figures/Figure 1E.png", units = 'in', width = 5.25, height = 5.57, res = 300) 
ggplot(DF) +
  aes(x = NMDS1, y = NMDS2, colour = tick_sex) +
  geom_point(shape = "circle", size = 3, alpha = 0.6) +
  scale_color_manual(
    values = palette_sex,
    labels = c("F" = "Female", "M" = "Male")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14, color = "black"),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = .8),
    title = element_text(size = 16),
  )+
  labs(title="Bray-Curtis")
dev.off()

### Statistical analysis ####

# bray curtis distance 

braydist <- phyloseq::distance(temp_rel, method = "bray")
sampledf <- data.frame(temp_rel@sam_data)
adonis2(braydist ~ tick_sex, data = sampledf,
        permutations = 999,padj = "bonferroni")

bd<-betadisper(braydist, sampledf$tick_sex)
anova(bd)
permutest(bd)

## SIMPER analysis (Table S1) ####

# Key question : What are the ASV and genera that contribute the most to the microbiota composition difference between 
# males and females ?

s <- vegan::simper(t(otu_table(temp_rel)), sample_data(temp_rel)$tick_sex)
summary(s)
summary_s3 <- summary(s)
average_contrib <- summary_s3[[1]]$average # extract values from column 'average'
total_contrib <- sum(average_contrib)      # Calculate total sum of mean contribution
percentage_contrib <- (average_contrib / total_contrib) * 100   # Calculate percentage contribution for each taxon
summary_s3[[1]]$percent_contrib <- percentage_contrib     # add the percentages values to the result dataframe
print(summary_s3[[1]])

# this file was exported into an excel file nammes "Simper_analysis" 
# and was formatted to create the Table S1

#  .   ####
# Temporal analysis ####
## Heat map (Fig S3) ####

# Key questions: What is the monthly variation in the detection of the genera ?
# Are there any genera that are not detected in all of the ticks collected during at least one month ?

### F ####

# dataframe preparation 

t_bin_f2 <- t_heat_f2 %>%
  mutate(Presence = ifelse(Abundance == 0, "0", "1"))
t_bin_f2 

as.data.frame(t_bin_f2)

t_bin_f2$Genus <- ifelse(
  t_bin_f2$Genus == "M.-aff. ASV (Comamonadaceae family)",
  "M.-aff. ASV\n(Comamonadaceae\nfamily)",
  t_bin_f2$Genus
)

t_bin_f2$Genus <- ifelse(
  t_bin_f2$Genus == "M.-aff. ASV (Moraxellaceae family)",
  "M.-aff. ASV\n(Moraxellaceae\nfamily)",
  t_bin_f2$Genus
)

t_bin_f2$Genus <- ifelse(
  t_bin_f2$Genus == "M.-aff. ASV (Bacilli class)",
  "M.-aff. ASV\n(Bacilli class)",
  t_bin_f2$Genus
)

# We identify the the genera that are not detected in all the ticks of at least one month in females

absent_taxa_f <- t_bin_f2 %>%
  group_by(month, Genus) %>%
  summarise(max_presence = max(Presence), .groups = "keep") %>%
  filter(max_presence == 0)
print(absent_taxa_f, n = 22)

# we plot he the genera that are not detected in all the  ticks of at least one month in males

t_bin_f2$month = factor(t_bin_f2$month,
                       levels = c("february","march", "april","may", "june", "july",
                                  "august", "september"))


png("figures/Figure S3A.png", units = 'in', width = 9.13, height = 4.19, res = 300)
t_bin_f2 %>%
  dplyr::filter(Genus %in% c("Brevundimonas","Caenimonas","Caviibacter","Chryseobacterium","Comamonas","Halomonas",
                      "Helcococcus","Mannheimia","Nocardioides", "Peptoniphilus","Porphyromonas","Pseudomonas","Rhodococcus",
                      "Sphingomonas","Streptococcus","Trueperella")) %>%
  ggplot() +
  aes(x = echantillons, y = Genus, fill = Presence) +
  geom_tile() +
  scale_fill_manual(
    values = c(`0` = "#EAEAEA",
               `1` = "#8B8B8B"))+
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black", face = "italic"),
    strip.background = element_rect(colour="black", fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.ticks.x = element_blank()
  )+
  labs(y="", x="") +
  facet_grid(vars(), vars(month), scales = "free_x")+
  theme(panel.spacing.x = unit(0.8, "cm"))
dev.off()

### M ####

# We identify the the genera that are not detected in all the  ticks of at least one month in males

genres_absents <- t_bin_m %>%
  group_by(month, Genus) %>%
  summarise(max_presence = max(Presence), .groups = "drop") %>%
  filter(max_presence == 0)
genres_absents

# we plot he the genera that are not detected in all the  ticks of at least one month in males

t_bin_m$month = factor(t_bin_m$month,
                       levels = c("february","march", "april","may", "june", "july",
                                  "august", "september"))

png("figures/Figure S3B.png", units = 'in', width = 9.13, height = 4.19, res = 300)
t_bin_m %>%
  filter(Genus %in% c("Caenimonas", "Caviibacter", "Helcococcus", "Mannheimia", 
                      'M.-aff. ASV\n(Comamonadaceae\nfamily)', "Nocardioides", "Peptoniphilus", "Porphyromonas", 
                      "Pseudomonas","Trueperella","Williamsia")) %>%
  ggplot() +
  aes(x = echantillons, y = Genus, fill = Presence) +
  geom_tile() +
  scale_fill_manual(
    values = c(`0` = "#EAEAEA",
               `1` = "#8B8B8B"))+
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black",face = "italic"),
    strip.background = element_rect(colour="black", fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.ticks.x = element_blank()
  )+
  labs(y="", x="") +
  facet_grid(vars(), vars(month), scales = "free_x")+
  theme(panel.spacing.x = unit(0.8, "cm"))
dev.off()

## Alpha diversity (Fig S2) ####

# Key question: Are the number of observed genera and their abundances influenced by the month ?

adiv_temp$month = factor(adiv_temp$month,
                         levels = c("february","march", "april","may", "june", "july",
                                    "august", "september"))

richness_t <- ggplot(adiv_temp) +
  aes(x = month, y = Richness, fill = tick_sex, color = tick_sex) +
  geom_boxplot(size = .75, width = .6, alpha = 0.4, colour = "#4B4E4D",outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3), size = 2, alpha = 0.6)+
  scale_fill_manual(values = palette_sex,labels = c("F" = "Female", "M" = "Male"))+
  scale_color_manual(values = palette_sex,labels = c("F" = "Female", "M" = "Male")) + 
  expand_limits(y = max(adiv_temp$t_chao) * 0.065) +
  theme_bw()+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust= 1, size = 14,color= "black"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = .8))+
  labs(title = "Richness",y = "")
richness_t

shann_t <- ggplot(adiv_temp) +
  aes(x = month, y = Shannon, fill = tick_sex, color = tick_sex) +
  geom_boxplot(size = .75, width = .6, alpha = 0.4, colour = "#4B4E4D",outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3), size = 2, alpha = 0.6)+
  scale_fill_manual(values = palette_sex,labels = c("F" = "Female", "M" = "Male"))+
  scale_color_manual(values = palette_sex,labels = c("F" = "Female", "M" = "Male")) + 
  expand_limits(y = max(adiv_temp$t_chao) * 0.018) +
  theme_bw()+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust= 1, size = 14,color= "black"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = .8))+
  labs(title = "Shannon",y = "")
shann_t 

invsimp  <- ggplot(adiv_temp) +
  aes(x = month, y = InvSimpson, fill = tick_sex, color = tick_sex) +
  geom_boxplot(size = .75, width = .6, alpha = 0.4, colour = "#4B4E4D",outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3), size = 2, alpha = 0.6)+
  scale_fill_manual(values = palette_sex,labels = c("F" = "Female", "M" = "Male"))+
  scale_color_manual(values = palette_sex,labels = c("F" = "Female", "M" = "Male")) + 
  expand_limits(y = max(adiv_temp$t_chao) * 0.065) +
  theme_bw()+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust= 1, size = 14,color= "black"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = .8))+
  labs(title = "Inverse Simpson index",y = "")
invsimp 

png("figures/Figure S2.png", units = 'in', width = 12.01, height = 4.17, res = 300)
richness_t+shann_t+invsimp+
  plot_layout(guides = 'collect', ncol = 3)+                    
  theme(plot.tag = element_text(face = 'bold', size = 25))  
dev.off()

### Statistical analysis month (Table 2) ####

#### F ####

# Graphic representations

tRichness_f<-estimate_richness(temporel_females ,measures = "Observed")
tShannon_f<-estimate_richness(temporel_females ,measures = "Shannon")
tInvSimpson_f<-estimate_richness(temporel_females ,measures = "InvSimpson")

adiv_temp_f <- cbind(temporel_females@sam_data, 
                     tRichness_f=tRichness_f$Observed,
                     tShannon_f=tShannon_f$Shannon,
                     tInvSimpson_f=tInvSimpson_f$InvSimpson)

adiv_temp_f

# Mean observed values 

SUM_rich_mf <-summarySE(adiv_temp_f, measurevar = "tRichness_f", groupvars = "month")
SUM_rich_mf

SUM_shan_mf <-summarySE(adiv_temp_f, measurevar = "tShannon_f", groupvars = "month")
SUM_shan_mf

SUM_invsimp_mf <-summarySE(adiv_temp_f, measurevar = "tInvSimpson_f", groupvars = "month")
SUM_invsimp_mf

SE_alpha_div_mf <- list(SUM_rich_mf,SUM_shan_mf,SUM_invsimp_mf)
#This file was used for the Table 2 in the manuscript

# richness
shapiro.test(adiv_temp_f$tRichness_f)
# richness follows a normal law

model_rich_f <- lm(tRichness_f ~ month, data = adiv_temp_f)
car::Anova(model_rich_f, type = "III")

# shannon
shapiro.test(adiv_temp_f$tShannon_f)
# shannon does not follow a normal law

model_shannon_f <- glm(tShannon_f ~ month, family=Gamma(link = "identity"), data = adiv_temp_f)
car::Anova(model_shannon_f, type = "III")

# inverse simpson
shapiro.test(adiv_temp_f$tInvSimpson_f)
# inverse simpson does not follow a normal law
model_invs_f <- glm(tInvSimpson_f ~ month, family = Gamma(link = "log"), data = adiv_temp_f)
car::Anova(model_invs_f, type = "III")

#### M ####

# Graphic representations

tRichness_m<-estimate_richness(temporel_males ,measures = "Observed")
tShannon_m<-estimate_richness(temporel_males ,measures = "Shannon")
tInvSimpson_m<-estimate_richness(temporel_males ,measures = "InvSimpson")

adiv_temp_m <- cbind(temporel_males@sam_data, 
                     tRichness_m=tRichness_m$Observed,
                     tShannon_m=tShannon_m$Shannon,
                     tInvSimpson_m=tInvSimpson_m$InvSimpson)

adiv_temp_m

# mean observed values
SUM_rich_mm <-summarySE(adiv_temp_m, measurevar = "tRichness_m", groupvars = "month")
SUM_rich_mm

SUM_shan_mm <-summarySE(adiv_temp_m, measurevar = "tShannon_m", groupvars = "month")
SUM_shan_mm

SUM_invsimp_mm <-summarySE(adiv_temp_m, measurevar = "tInvSimpson_m", groupvars = "month")
SUM_invsimp_mm

SE_alpha_div_mm <- list(SUM_rich_mm,SUM_shan_mm,SUM_invsimp_mm)

#This file was used for the Table 2 in the manuscript

# richness
shapiro.test(adiv_temp_m$tRichness_m)
# richness follows a normal law
hist(adiv_temp_m$tRichness_m)

model_rich_m <- lm(tRichness_m ~ month, data = adiv_temp_m)
car::Anova(model_rich_m, type = "III")

# shannon
shapiro.test(adiv_temp_m$tShannon_m)
# inverse simpson does not follow a normal law
hist(adiv_temp_m$tShannon_m)

model_shannon_m<- lm(tShannon_m ~ month, data = adiv_temp_m)
car::Anova(model_shannon_m, type = "III")

# inverse simpson
shapiro.test(adiv_temp_m$tInvSimpson_m)
# inverse simpson does not follow a normal law
hist(adiv_temp_m$tInvSimpson_m)

model_invs_m <- glm(tInvSimpson_m ~ month, family = Gamma(link = "log"), data = adiv_temp_m)
car::Anova(model_invs_m, type = "III")

### Statistical analysis season ####

# Is the microbiota of ticks more diverse in winter, spring or summer ?

#### F ####
#richness
model_rich_sf <- lm(tRichness_f ~ season, data = adiv_temp_f)
car::Anova(model_rich_sf, type = "III")
summary(model_rich_sf)

emmeans_rich_sf<-emmeans(model_rich_sf,pairwise~season,type="response",padj = "bonferroni")  
emmeans_rich_sf

# Shannon
model_shannon_sf <- glm(tShannon_f ~ season, family=Gamma(link = "identity"), data = adiv_temp_f)
car::Anova(model_shannon_sf, type = "III")

# InvSimpson
model_invs_sf <- glm(tInvSimpson_f ~ season, family = Gamma(link = "log"), data = adiv_temp_f)
car::Anova(model_invs_sf, type = "III")

#### M ####
#richness
model_rich_sm <- lm(tRichness_m ~ season, data = adiv_temp_m)
car::Anova(model_rich_sm, type = "III")
summary(model_rich_sm)

emmeans_rich_sm<-emmeans(model_rich_sm,pairwise~season,type="response",padj = "bonferroni")  
emmeans_rich_sm

# Shannon
model_shannon_sm <- lm(tShannon_m ~ season, data = adiv_temp_m)
car::Anova(model_shannon_sm, type = "III")

emmeans_rich_msh<-emmeans(model_shannon_sm,pairwise~season,type="response",padj = "bonferroni")  
emmeans_rich_msh

# InvSimpson
model_invs_sm <- glm(tInvSimpson_m ~ season, family = Gamma(link = "log"), data = adiv_temp_m)
car::Anova(model_invs_sm, type = "III")

emmeans_rich_msim<-emmeans(model_invs_sm,pairwise~season,type="response", padj = "bonferroni")  
emmeans_rich_msim

## Beta diversity (Fig 3AC)####

# Key question : Is the bacterial structure (presence and abundance of taxa) influenced by the month ?

### F ####

temporel_females=subset_samples(temporel_subset_seq3, tick_sex =="F")
temp_rel_f = transform_sample_counts(temporel_females, function(x) x/sum(x)*100)

temp.BC3<-ordinate(temp_rel_f,method = "NMDS",distance = "bray", k=4)
temp.BC3

palette_month <- c("february" = "#461DEC", "march" = "#99CDF5", "april" = "#3ED1A7", "may"= "#2BD552", "june"= "#BBE33E", 
                   "july"= "#E3613E", "august"="#D05A78", "september"="#E497C7")

DF3<-plot_ordination(temp_rel_f, temp.BC3, color = "month", justDF = TRUE)
DF3$month = factor(DF3$month,
                   levels = c("february","march", "april","may", "june", "july",
                              "august", "september"))


png("figures/Figure 3A.png", units = 'in', width = 6.41, height = 4.70, res = 300) 
NMDS_month_f<-
  ggplot(DF3) +
  aes(x = NMDS1, y = NMDS2, colour = month) +
  geom_point(shape = "circle", size = 5, alpha = 0.7) +
  scale_color_paletteer_d("awtools::a_palette")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position.inside = c(0.8, 0.2),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.x = element_text(size=16, color = "black"),
    axis.title.y = element_text(size=16, color = "black")
  )
NMDS_month_f
dev.off()


### M ####

temporel_males=subset_samples(temporel_subset_seq3, tick_sex =="M")
temp_rel_m = transform_sample_counts(temporel_males, function(x) x/sum(x)*100)

# Compute the Jaccard distance matrix
temp.jaccard_m<-ordinate(temp_rel_m,method = "PCoA",distance = "jaccard",binary=T)

temp.BC4<-ordinate(temp_rel_m,method = "NMDS",distance = "bray", k=4)
temp.BC4

palette_month <- c("february" = "#461DEC", "march" = "#99CDF5", "april" = "#3ED1A7", "may"= "#2BD552", "june"= "#BBE33E", 
                   "july"= "#E3613E", "august"="#D05A78", "september"="#E497C7")

DF4<-plot_ordination(temp_rel_m, temp.BC4, color = "month", justDF = TRUE)
DF4$month = factor(DF4$month,
                   levels = c("february","march", "april","may", "june", "july",
                              "august", "september"))


png("figures/Figure 3B.png", units = 'in', width = 6.41, height = 4.70, res = 300)
NMDS_month_m <-
  ggplot(DF4) +
  aes(x = NMDS1, y = NMDS2, colour = month) +
  geom_point(shape = "circle", size = 5, alpha = 0.7) +
  scale_color_paletteer_d("awtools::a_palette")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position.inside = c(0.8, 0.2),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16, color = "black"),
    axis.title.y = element_text(size=16, color = "black"))
NMDS_month_m
dev.off()

### Statistical analysis ####

# Key question : Is the bacterial structure (presence and abundance of taxa) statistically influenced by the month?

#### F ####
set.seed(123)
braydist5 <- phyloseq::distance(temp_rel_f, method = "bray")
sampledf5 <- data.frame(temp_rel_f@sam_data)
adonis2(braydist5 ~ month, data = sampledf5,
                         permutations = 999)

#### M ####
braydist1 <- phyloseq::distance(temp_rel_m, method = "bray")
sampledf1 <- data.frame(temp_rel_m@sam_data)
adonis2(braydist1 ~ month, data = sampledf1,
        permutations = 999)

## Taxonomic composition ####
### Barplots ####
#### F (Fig 3B)####

# Key question : What is the composition of the microbiota at the level of the genus across the months ? 
# What proportion of relative abundance represent each genus ?

t6_compo <- psmelt(temporel_subset_seq3)
t6_compo
t6_compo_f <-  subset(t6_compo, tick_sex == "F")
t6_agg_f_month <- aggregate(Abundance ~ month + Phylum + Class + Order + Family + Genus, data = t6_compo_f, FUN = sum)
t6_agg_f_month$rel_abund <- with(t6_agg_f_month, Abundance / ave(Abundance, month,FUN = sum))
t6_agg_f_month$rel_abund <- t6_agg_f_month$rel_abund * 100

# Relative abundance of genera across month
t6_top_phy_f <- t6_agg_f_month %>%
  group_by(month,Phylum, Class,Order,Family,Genus) %>%
  summarize(rel_abund_med = median(rel_abund))
t6_top_phy_f

# if the above function does not work please try to shut down R studio, reopen it and recharge the libraries.

# Identify phyla for which relative abundance median is above 2.5% of total abundance
t6_keep_f <- unique(t6_top_phy_f$Genus[t6_top_phy_f$rel_abund_med > 2.5])
t6_top_phy_f$Genus[!(t6_top_phy_f$Genus %in% t6_keep_f)] <- "< 2.5%"

# aggregrate data per month
t6_agg_sum_f <- t6_top_phy_f %>%
  group_by(month, Genus) %>%
  summarize(rel_abund_sum = sum(rel_abund_med))

t6_agg_sum_f$month = factor(t6_agg_sum_f$month,
                            levels = c("february","march", "april","may", "june", "july",
                                       "august", "september"))

genus_order_f <- c(
  "< 2.5%",  
  "Corynebacterium",
  "M.-aff. ASV (Bacilli class)",
  'Midichloria',
  "Francisella",
  "Rickettsia"
)

t6_agg_sum_f$Genus <- factor(t6_agg_sum_f$Genus, levels = genus_order_f)

png("figures/Figure 3B.png", units = 'in', width = 8.08, height = 6.84, res = 300) 
barplot_f <-
  ggplot(t6_agg_sum_f, aes(x = month, y = rel_abund_sum, fill = Genus)) +
  geom_col(width = 0.8) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(
    values = c(`< 2.5%` = "#000000",
               Corynebacterium = "#12DDCC",
               `M.-aff. ASV (Bacilli class)` = "#98F273",
               `Midichloria` = "#E3BCF6",
               Francisella = "#FFD8F1",
               Rickettsia = "#DA0062"),
    guide = guide_legend(
      title = NULL,
      title.position = "top",
      label.position = "right",
      keywidth = 1, 
      keyheight = 1.1, 
      default.unit = "cm")
  )+
  labs(x = "", y = "Relative abundance (%)") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(size = 18, angle = 50, hjust = 1, vjust = 1,colour = "black"),
    axis.text.y = element_text(size = 18,colour = "black"),
    axis.title.y = element_text(size = 20,colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(size = 18),
    legend.text= element_text(size = 15, face = "italic")
  )
barplot_f
dev.off()

#### M (Fig 3D) ####
t6_compo_m <-  subset(t6_compo, tick_sex == "M")
t6_agg_m_month <- aggregate(Abundance ~ month + Phylum + Class + Order + Family + Genus, data = t6_compo_m, FUN = sum)
t6_agg_m_month$rel_abund <- with(t6_agg_m_month, Abundance / ave(Abundance, month,FUN = sum))
t6_agg_m_month$rel_abund <- t6_agg_m_month$rel_abund * 100

t6_top_phy_m <- t6_agg_m_month %>%
  group_by(month,Phylum, Class,Order,Family,Genus) %>%
  summarize(rel_abund_med = median(rel_abund))

t6_keep_m <- unique(t6_top_phy_m$Genus[t6_top_phy_m$rel_abund_med > 2.5])
t6_top_phy_m$Genus[!(t6_top_phy_m$Genus %in% t6_keep_m)] <- "< 2.5%"

t6_agg_sum_m <- t6_top_phy_m %>%
  group_by(month, Genus) %>%
  summarize(rel_abund_sum = sum(rel_abund_med))

t6_agg_sum_m$month = factor(t6_agg_sum_m$month,
                            levels = c("february","march", "april","may", "june", "july",
                                       "august", "september"))


t6_agg_sum_m$Genus <- ifelse(t6_agg_sum_m$Genus == "M.-aff. ASV (Moraxellaceae family)", 
                             "M.-aff. ASV\n(Moraxellaceae family)", 
                             t6_agg_sum_m$Genus)

genus_order <- c(
  "< 2.5%",  
  "Mycobacterium",  
  "Corynebacterium",
  "M.-aff. ASV (Bacilli class)",
  "Acinetobacter",
  'Midichloria',
  "Francisella",
  "M.-aff. ASV\n(Moraxellaceae family)",
  "Rickettsia"
)

t6_agg_sum_m$Genus <- factor(t6_agg_sum_m$Genus, levels = genus_order)

png("figures/Figure 3D.png", units = 'in', width = 8.08, height = 6.84, res = 300) 
barplot_m <-
  ggplot(t6_agg_sum_m, aes(x = month, y = rel_abund_sum, fill = Genus)) +
  geom_col(width = 0.8) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(
    values = c(`< 2.5%` = "#000000",
               Corynebacterium = "#12DDCC",
               Mycobacterium = "#0EADA2",
               `M.-aff. ASV (Bacilli class)` = "#98F273",
               Acinetobacter = "#C277DA",
               `Midichloria` = "#E3BCF6",
               Francisella = "#FFD8F1",
               `M.-aff. ASV\n(Moraxellaceae family)` = "#FFA3E7",
               Rickettsia = "#DA0062"),
    guide = guide_legend(
      title = NULL,
      title.position = "top",
      label.position = "right",
      keywidth = 1, # Ajuster la largeur des clés de légende
      keyheight = 1.1, # Ajuster la hauteur des clés de légende
      default.unit = "cm")
  )+
  labs(x = "", y = "Relative abundance (%)") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(size = 18, angle = 50, hjust = 1, vjust = 1,colour = "black"),
    axis.text.y = element_text(size = 18,colour = "black"),
    axis.title.y = element_text(size = 20,colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(size = 18),
    legend.text= element_text(size = 15, face = "italic")
  )
barplot_m
dev.off()

### Boxplots ####
#### Statistics (Table S2) ####
##### F ####

# Key question : Are the relative abundance of the most abundant taxa influenced by the month ?

# Kruskal Wallis tests

df_KWfemales <- psmelt(temporel_females) %>%
  mutate(month = factor(month)) %>%
  group_by(Sample, month, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(relative_abundance = (Abundance / sum(Abundance)) * 100) %>%
  ungroup()
df_KWfemales

genera_keep <- c("Francisella", "Midichloria", "Rickettsia","Corynebacterium", "Acinetobacter","Mycobacterium", "Williamsia","M.-aff. ASV (Moraxellaceae family)", "M.-aff. ASV (Bacilli class)" )
df_KWfemales_filt <- df_KWfemales %>%
  filter(Genus %in% genera_keep)

results_KWF <- df_KWfemales_filt %>%
  group_by(Genus) %>%
  kruskal_test(relative_abundance ~ month)

results_KWF <- results_KWF %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

results_KWF

##### M ####
df_KWmales <- psmelt(temporel_males) %>%
  mutate(month = factor(month)) %>%
  group_by(Sample, month, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(relative_abundance = (Abundance / sum(Abundance)) * 100) %>%
  ungroup()
df_KWmales

genera_keep <- c("Francisella", "Midichloria", "Rickettsia","Corynebacterium", "Acinetobacter","Mycobacterium", "Williamsia", "M.-aff. ASV (Moraxellaceae family)", "M.-aff. ASV (Bacilli class)" )
df_KWmales_filt <- df_KWmales %>%
  filter(Genus %in% genera_keep)

results_KWM <- df_KWmales_filt %>%
  group_by(Genus) %>%
  kruskal_test(relative_abundance ~ month)

results_KWM <- results_KWM %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

results_KWM


#### Graphic representation ####

##### F (Fig S4A) #####
df_KWfemales_filt$month = factor(df_KWfemales_filt$month,
                                 levels = c("february","march", "april","may", "june", "july",
                                            "august", "september"))

df_KWfemales_filt$Genus <- ifelse(df_KWfemales_filt$Genus == "M.-aff. ASV (Moraxellaceae family)", 
                                  "M.-aff. ASV\n(Moraxellaceae family)",
                                  df_KWfemales_filt$Genus)

df_KWfemales_filt$Genus <- ifelse(df_KWfemales_filt$Genus == "M.-aff. ASV (Bacilli class)", 
                                  "M.-aff. ASV\n(Bacilli class)", 
                                  df_KWfemales_filt$Genus)

ordered_genera <- c(
  "Francisella", 
  "Midichloria", 
  "Rickettsia", 
  "Corynebacterium", 
  "Acinetobacter", 
  "Mycobacterium", 
  "Williamsia",
  "M.-aff. ASV\n(Bacilli class)", 
  "M.-aff. ASV\n(Moraxellaceae family)"
)

# Convertir Genus en facteur avec l’ordre imposé
df_KWfemales_filt$Genus <- factor(df_KWfemales_filt$Genus, levels = ordered_genera)

# Générer le graphique avec facettes dans l'ordre voulu
png("figures/Figure S4A.png", units = 'in', width = 8.97, height = 7.40, res = 300) 
boxplots_F <- df_KWfemales_filt %>%
  ggplot(aes(x = month, y = relative_abundance)) +
  geom_boxplot(aes(fill = month), size = 0.6, outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(color = month), alpha = 0.8, size = 2, width = 0.2) +
  scale_color_paletteer_d("awtools::a_palette") +
  scale_fill_paletteer_d("awtools::a_palette") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 70, hjust = 1, vjust = 1, colour = "black"),
    strip.text = element_text(face = "italic", size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  ) +
  labs(y = "Relative Abundance (%)") +
  facet_wrap(vars(Genus), scales = "free_y")
boxplots_F
dev.off()

##### M (Fig S4B) ####
df_KWmales_filt$month = factor(df_KWmales_filt$month,
                               levels = c("february","march", "april","may", "june", "july",
                                          "august", "september"))

df_KWmales_filt$Genus <- ifelse(df_KWmales_filt$Genus == "M.-aff. ASV (Moraxellaceae family)", 
                                "M.-aff. ASV\n(Moraxellaceae family)",
                                df_KWmales_filt$Genus)

df_KWmales_filt$Genus <- ifelse(df_KWmales_filt$Genus == "M.-aff. ASV (Bacilli class)", 
                                "M.-aff. ASV\n(Bacilli class)", 
                                df_KWmales_filt$Genus)

df_KWmales_filt$Genus <- factor(df_KWmales_filt$Genus, levels = ordered_genera)

png("figures/Figure S4B.png", units = 'in', width = 8.97, height = 7.40, res = 300) 
boxplots_M <- df_KWmales_filt %>%
  ggplot(aes(x = month, y = relative_abundance)) +
  geom_boxplot(aes(fill = month), size = 0.6, outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(color = month), alpha = 0.8, size = 2, width = 0.2) +
  scale_color_paletteer_d("awtools::a_palette") +
  scale_fill_paletteer_d("awtools::a_palette") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 70, hjust = 1, vjust = 1, colour = "black"),
    strip.text = element_text(face = "italic", size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  ) +
  labs(y = "Relative Abundance (%)") +
  facet_wrap(vars(Genus), scales = "free_y")
boxplots_M
dev.off()

