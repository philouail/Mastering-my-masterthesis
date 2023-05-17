library(ggpubr)
library(data.table)
library(tidyverse)

setwd("C:/Users/phili/OneDrive - Uppsala universitet/Master_thesis/QuPath_output")

###Panel analysis####
#Malignant - detection
setwd("malignant")

resultsTotal <- data.frame()
listtxt <- dir(pattern = "*.txt")

## Extracting results from QuPath output
for (s in 1:length(listtxt)) {
  print(s)
  annotation <- read.delim(listtxt[s], header=TRUE)
  rownames(annotation) <- annotation$Name
  numberofcell<- annotation["PathAnnotationObject","Num.Detections"]
  CD8Tcell_percent <- annotation["PathAnnotationObject", "Num.CD8A"]/numberofcell
  Tregcell_percent <- annotation["PathAnnotationObject", "Num.FOXP3"]/numberofcell
  Bcell_percent <- annotation["Stroma", "Num.CD79A"]/numberofcell #Only in stroma bc weird signal in tumor cells
  Macro_percent <- annotation["PathAnnotationObject", "Num.CD163"]/numberofcell
  
  #PDL1 info 
  PDL1_percent <- annotation["PDL1", "Area.µm.2"]/annotation["PathAnnotationObject", "Area.µm.2"]
  CD163_PDL1 <- annotation["PDL1", "Num.CD163"]/annotation["PDL1","Num.Detections"]
  CD8_PDL1 <- annotation["PDL1", "Num.CD8A"]/annotation["PDL1","Num.Detections"]
  other_PDL1 <-  1 - (CD163_PDL1 + CD8_PDL1)

  #infiltrated cells in tumor area
  CD8_infiltrated <- annotation["Tumor", "Num.CD8A"]/annotation["Tumor", "Num.Detections"]
  CD163_infiltrated <- annotation["Tumor", "Num.CD163"]/annotation["Tumor", "Num.Detections"]
  Tumor_area <- annotation["Tumor", "Area.µm.2"]/annotation["PathAnnotationObject", "Area.µm.2"]
  Tumor_cell <- 1 - (CD163_infiltrated + CD8_infiltrated)
 
  resultsslide <-  data.frame(numberofcell, 
                      CD8Tcell_percent = ifelse(length(CD8Tcell_percent) == 0, NA, CD8Tcell_percent), 
                      Tregcell_percent = ifelse(length(Tregcell_percent) == 0, NA, Tregcell_percent) , 
                      Bcell_percent = ifelse(length(Bcell_percent) == 0, NA, Bcell_percent), 
                      Macro_percent = ifelse(length(Macro_percent) == 0, NA, Macro_percent), 
                      PDL1_percent = ifelse(length(PDL1_percent) == 0, NA, PDL1_percent),
                      CD163_PDL1 = ifelse(length(CD163_PDL1) == 0, NA, CD163_PDL1), 
                      CD8_PDL1 = ifelse(length(CD8_PDL1) == 0, NA, CD8_PDL1), 
                      other_PDL1 = ifelse(length(other_PDL1) == 0, NA, other_PDL1),
                      Tumor_area = ifelse(length(Tumor_area) == 0, NA, Tumor_area),
                      Tumor_cell = ifelse(length(Tumor_cell) == 0, NA, Tumor_cell),
                      CD163_infiltrated= ifelse(length(CD163_infiltrated) == 0, NA, CD163_infiltrated) , 
                      CD8_infiltrated= ifelse(length(CD8_infiltrated) == 0, NA, CD8_infiltrated),
                      row.names = annotation[1,1])
  if ((nrow(resultsTotal) == 0)) {
    resultsTotal <- resultsslide
  } else {
    resultsTotal <- bind_rows(resultsTotal, resultsslide)
  }
}


write.csv(resultsTotal, "Annotation_malignant.csv")
##here manually open the excel file and just add a patient column. 
##Ideally you organise your slide better to begin with so you don't have to do that. 


setwd("C:/Users/phili/OneDrive - Uppsala universitet/Master_thesis/QuPath_output")
#Non- Malignant
setwd("Non_malignant")
##same code but not tumor stroma info 
resultsTotal <- data.frame()
listtxt <- dir(pattern = "*.txt")

## Extracting results from QuPath output
for (s in 1:length(listtxt)) {
  print(s)
  annotation <- read.delim(listtxt[s], header=TRUE)
  rownames(annotation) <- annotation$Name
  numberofcell<- annotation["PathAnnotationObject","Num.Detections"] 
  CD8Tcell_percent <- annotation["PathAnnotationObject", "Num.CD8A"]/numberofcell
  Tregcell_percent <- annotation["PathAnnotationObject", "Num.FOXP3"]/numberofcell
  Bcell_percent <- annotation["PathAnnotationObject", "Num.CD79A"]/numberofcell #Only in stroma bc weird signal in tumor cells
  Macro_percent <- annotation["PathAnnotationObject", "Num.CD163"]/numberofcell
  
  #PDL1 info 
  PDL1_percent <- annotation["PDL1", "Area.µm.2"]/annotation["PathAnnotationObject", "Area.µm.2"]
  CD163_PDL1 <- annotation["PDL1", "Num.CD163"]/annotation["PDL1","Num.Detections"]
  CD8_PDL1 <- annotation["PDL1", "Num.CD8A"]/annotation["PDL1","Num.Detections"]
  other_PDL1 <-  1 - (CD163_PDL1 + CD8_PDL1)
  
  resultsslide <- data.frame(numberofcell, 
                             CD8Tcell_percent = ifelse(length(CD8Tcell_percent) == 0, NA, CD8Tcell_percent), 
                             Tregcell_percent = ifelse(length(Tregcell_percent) == 0, NA, Tregcell_percent) , 
                             Bcell_percent = ifelse(length(Bcell_percent) == 0, NA, Bcell_percent), 
                             Macro_percent = ifelse(length(Macro_percent) == 0, NA, Macro_percent), 
                             PDL1_percent = ifelse(length(PDL1_percent) == 0, NA, PDL1_percent),
                             CD163_PDL1 = ifelse(length(CD163_PDL1) == 0, NA, CD163_PDL1), 
                             CD8_PDL1 = ifelse(length(CD8_PDL1) == 0, NA, CD8_PDL1), 
                             other_PDL1 = ifelse(length(other_PDL1) == 0, NA, other_PDL1),
                             row.names = annotation[1,1])
  if ((nrow(resultsTotal) == 0)) {
    resultsTotal <- resultsslide
  } else {
    resultsTotal <- bind_rows(resultsTotal, resultsslide)
  }
}

write.csv(resultsTotal, "Annotation_non_malignant.csv")

setwd("C:/Users/phili/OneDrive - Uppsala universitet/Master_thesis/QuPath_output")
#load both file
Mal_results <- fread("Analysed_data/Annotation_malignant.csv")
NonMal_results <- fread("Analysed_data/Annotation_non_malignant.csv")


#Combine
All_results <- bind_rows(Mal_results, NonMal_results)
##will use both combined and non-combined results downstream


#graph now 

##Immune cell composition for each patient 

graphlist <- c("CD8Tcell_percent","Tregcell_percent", "Bcell_percent", "Macro_percent")
titles <- c("CD8 T cells", "Regulatory Tcells", "B cells", "Macrophages")
plot_list <- list()
for (g in 1:length(graphlist)){

  formula <- as.formula(paste(graphlist[g], "~ Patient_ID"))
  
  results <- compare_means( formula , data= All_results)
  significant_pairs <- results %>% 
    filter(p < 0.05) %>%
    select(group1,group2)
  
  significant_pairs_list <- as.list(as.data.frame(t(significant_pairs)))
  #make a list of all the significant results - can be looped
  
  p <- ggboxplot(All_results,
            x="Patient_ID",
            y= graphlist[g],
            add = "jitter", 
            title = titles[g],
            order= c("77.3","734.5", "1039.4","749.1"),
            color = "Patient_ID",
            ylab = "Percentage of all Cells",
            xlab= "Patient ID",
            add.params = list(fill = "white", size = 1.5, alpha = .2),
            palette =  c("77.3" = "#8302a2","734.5" = "#31788e","1039.4" = "#d60f23","749.1" = "grey"),
            legend = "none", 
            repel= TRUE) + 
    scale_y_continuous(labels = scales::percent) + 
    stat_compare_means(comparisons = significant_pairs_list)
  
    plot_list[[g]] <- p 
  
  pdf(paste0(titles[g], ".pdf"))
  print(p)
  dev.off()
}

tp <- ggarrange(plotlist = plot_list, labels = c("A", "B", "C", "D" ))

tp


pdf(paste0("Immune_cell_composition.pdf"), width = 8 , height = 11)
print(tp)
dev.off()


#PDL1_cells

graphlist <- c( "PDL1_percent", "CD163_PDL1", "CD8_PDL1")
titles <- c( "PDL1", "Macrophages", "CD8" )
plot_list <- list()
for (g in 1:length(graphlist)){
  
  formula <- as.formula(paste(graphlist[g], "~ Patient_ID"))
  
  results <- compare_means( formula , data= All_results)
  significant_pairs <- results %>% 
    filter(p < 0.05) %>%
    select(group1,group2)
  
  significant_pairs_list <- as.list(as.data.frame(t(significant_pairs)))
  #make a list of all the significant results - can be looped
  
  p <- ggboxplot(All_results,
                 x="Patient_ID",
                 y= graphlist[g],
                 add = "jitter", 
                 title = titles[g],
                 order= c("77.3","734.5", "1039.4","749.1"),
                 color = "Patient_ID",
                 ylab = paste("Proportion of", titles[g], "positive cells"),
                 xlab= "Patient ID",
                 add.params = list(fill = "white", size = 1.5, alpha = .2),
                 palette =  c("77.3" = "#8302a2","734.5" = "#31788e","1039.4" = "#d60f23","749.1" = "grey"),
                 legend = "none") + 
    scale_y_continuous(labels = scales::percent) + 
    stat_compare_means(comparisons = significant_pairs_list)
  
  pdf(paste0(graphlist[g], ".pdf"))
  print(p)
  dev.off()
  
  plot_list[[g]] <- p 
  
}

tp <- ggarrange(plotlist = plot_list, labels = c("A", "B", "C"), common.legend = TRUE, legend = "bottom", nrow = 1 )

tp
pdf(paste0("PDL1_cell.pdf"), width = 11 , height = 8)
print(tp)
dev.off()

#Infiltrated cells
graphlist <- c("Tumor_area","CD163_infiltrated" , "CD8_infiltrated")
titles <- c("Tumor","Macrophages", "CD8")
plot_list <- list()

for (g in 1:length(graphlist)){
  
  formula <- as.formula(paste(graphlist[g], "~ Patient_ID"))
  
  results <- compare_means( formula , data= Mal_results)
  significant_pairs <- results %>% 
    filter(p < 0.05) %>%
    select(group1,group2)
  
  significant_pairs_list <- as.list(as.data.frame(t(significant_pairs)))
  #make a list of all the significant results - can be looped
  
  p <- ggboxplot(Mal_results,
                 x="Patient_ID",
                 y= graphlist[g],
                 add = "jitter", 
                 title = titles[g],
                 order= c("77.3","734.5", "1039.4","749.1"),
                 color = "Patient_ID",
                 ylab = paste("Proportion of", titles[g], "positive cells"),
                 xlab= "Patient ID",
                 add.params = list(fill = "white", size = 1.5, alpha = .2),
                 palette =  c("77.3" = "#8302a2","734.5" = "#31788e","1039.4" = "#d60f23","749.1" = "grey"),
                 legend = "none") + 
    scale_y_continuous(labels = scales::percent) + 
    stat_compare_means(comparisons = significant_pairs_list)
  
  pdf(paste0(graphlist[g], ".pdf"))
  print(p)
  dev.off()
  
  plot_list[[g]] <- p 
  
}

tp <- ggarrange(plotlist = plot_list, labels = c("A", "B", "C"), common.legend = TRUE, legend = "bottom", nrow=1 )

tp 
pdf(paste0("Immune_cell_infiltration.pdf"), width = 11 , height = 8)
print(tp)
dev.off()



###Unknown marker analysis###
#This part of the analysis will need to be adapted to the unknown marker, and how you want to analyize it
#SLAMF7 analysis 
## Extracting results from QuPath output
setwd("C:/Users/phili/OneDrive - Uppsala universitet/Master_thesis/QuPath_output/SLAMF7")

Percent_of_total_area <- data.frame()
composition <- data.frame()
listtxt <- dir(pattern = "*.txt")

for (s in 1:length(listtxt)) {
  print(s)
  annotation <- read.delim(listtxt[s], header=TRUE)
  rownames(annotation) <- annotation$Name
  total_area<- annotation["PathAnnotationObject","Area.µm.2"]
  SLAMF7_percent <- annotation["SLAMF7", "Area.µm.2"]/total_area
  
  
  resultsslide <-  data.frame(total_area, 
                              SLAMF7_percent = ifelse(is.na(SLAMF7_percent) == TRUE, 0, SLAMF7_percent),
                              row.names = annotation[1,1])
  if ((nrow(Percent_of_total_area) == 0)) {
    Percent_of_total_area <- resultsslide
  } else {
    Percent_of_total_area <- bind_rows(Percent_of_total_area, resultsslide)
  }
  
  
  Composition_line <- annotation["SLAMF7", c("Image", "Num.CD163", "Num.CD79A", "Num.CD8A", "Num.FOXP3", "Num.Non_immune.cells")]
  if (is.na(Composition_line[1,1]) == TRUE) {
    next
  }
  if (nrow(composition) == 0) {
    composition <- Composition_line
  } else {
    composition <- bind_rows(composition, Composition_line)
  }
}

write.csv(composition, "Composition.csv")
write.csv(Percent_of_total_area, "Percent_of_Total_area.csv")

#here again i have manually added a column "Patient ID" to both dataframe 

Composition <- fread("Composition.csv")
Percent_of_total_area <- fread("Percent_of_Total_area.csv")

results <- compare_means( SLAMF7_percent ~ Patient_ID , data= Percent_of_total_area)
significant_pairs <- results %>% 
  filter(p < 0.05) %>%
  select(group1,group2)

significant_pairs_list <- as.list(as.data.frame(t(significant_pairs)))
#nothing is significant
#make a list of all the significant results - can be looped

a <- ggboxplot(Percent_of_total_area,
               x="Patient_ID",
               y= "SLAMF7_percent",
               add = "jitter", 
               title = "Expression of SLAMF7 in different patients",
               order= c("77.3","734.5", "1039.4","749.1"),
               color = "Patient_ID",
               ylab = "Percentage of total Area",
               xlab= "Patient ID",
               add.params = list(fill = "white", size = 1.5, alpha = .2),
               palette =  c("77.3" = "#8302a2","734.5" = "#31788e","1039.4" = "#d60f23","749.1" = "grey"),
               legend = "none") + 
  scale_y_continuous(labels = scales::percent) + 
  stat_compare_means(comparisons = significant_pairs_list)

pdf(paste0("SLAMF7", ".pdf"))
print(p)
dev.off()

Comp_long <- Composition %>% pivot_longer(cols = c("Num.CD163", "Num.CD79A", "Num.CD8A", "Num.FOXP3", "Num.Non_immune.cells"), names_to = "Cell_types", values_to = "Cell_count")
Comp_long <- Comp_long[,3:5]

Comp_long$Cell_types[Comp_long$Cell_types == "Num.FOXP3"] <- "FOXP3"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.CD163"] <- "CD163"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.CD79A"] <- "CD79A"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.Non_immune.cells"] <- "Others"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.CD8A"] <- "CD8A"

fill_colors <- c("FOXP3" = "#8302a2", "CD163" = "pink", "CD79A" = "#d60f23", "Others" = "grey", "CD8A" = "black")

b <- ggbarplot(Comp_long, x = "Patient_ID" , y = "Cell_count", fill= "Cell_types", color = "Cell_types", 
          order= c("77.3","734.5", "1039.4","749.1"), add = c("mean"),
          ylim= c(0,750), ylab = "Cell count",
          xlab= "Patient ID",
          palette = fill_colors,legend= "bottom", legend.title = "" , title = "Cells expressing SLAMF7")

pdf(paste0("Comp_SLAMF7", ".pdf"))
print(g)
dev.off()


#GZMK analysis 
## Extracting results from QuPath output
setwd("C:/Users/phili/OneDrive - Uppsala universitet/Master_thesis/QuPath_output/GZMK")

Percent_of_total_area <- data.frame()
composition <- data.frame()
listtxt <- dir(pattern = "*.txt")

for (s in 1:length(listtxt)) {
  print(s)
  annotation <- read.delim(listtxt[s], header=TRUE)
  rownames(annotation) <- annotation$Name
  total_area<- annotation["PathAnnotationObject","Area.µm.2"]
  GZMK_percent <- annotation["GZMK", "Area.µm.2"]/total_area
  
  
  resultsslide <-  data.frame(total_area, 
                              GZMK_percent = ifelse(is.na(GZMK_percent) == TRUE, 0, GZMK_percent),
                              row.names = annotation[1,1])
  if ((nrow(Percent_of_total_area) == 0)) {
    Percent_of_total_area <- resultsslide
  } else {
    Percent_of_total_area <- bind_rows(Percent_of_total_area, resultsslide)
  }
  
  
  Composition_line <- annotation["GZMK", c("Image", "Num.CD163", "Num.CD79A", "Num.CD8A", "Num.FOXP3", "Num.Non_immune.cells")]
  if (is.na(Composition_line[1,1]) == TRUE) {
    next
  }
  if (nrow(composition) == 0) {
    composition <- Composition_line
  } else {
    composition <- bind_rows(composition, Composition_line)
  }
}

write.csv(composition, "Composition.csv")
write.csv(Percent_of_total_area, "Percent_of_Total_area.csv")

#here again i have manually added a column "Patient ID" to both dataframe 

Composition <- fread("Composition.csv")
Percent_of_total_area <- fread("Percent_of_Total_area.csv")

results <- compare_means( GZMK_percent ~ Patient_ID , data= Percent_of_total_area)
significant_pairs <- results %>% 
  filter(p < 0.05) %>%
  select(group1,group2)

significant_pairs_list <- as.list(as.data.frame(t(significant_pairs)))
#nothing is significant
#make a list of all the significant results - can be looped

p <- ggboxplot(Percent_of_total_area,
               x="Patient_ID",
               y= "GZMK_percent",
               add = "jitter", 
               title = "Expression of GZMK in different patients",
               order= c("77.3","734.5", "1039.4","749.1"),
               color = "Patient_ID",
               ylab = "Percentage of total Area",
               xlab= "Patient ID",
               add.params = list(fill = "white", size = 1.5, alpha = .2),
               palette =  c("77.3" = "#8302a2","734.5" = "#31788e","1039.4" = "#d60f23","749.1" = "grey"),
               legend = "none") + 
  scale_y_continuous(labels = scales::percent) + 
  stat_compare_means(comparisons = significant_pairs_list)
p
pdf(paste0("GZMK", ".pdf"))
print(p)
dev.off()

Comp_long <- Composition %>% pivot_longer(cols = c("Num.CD163", "Num.CD79A", "Num.CD8A", "Num.FOXP3", "Num.Non_immune.cells"), names_to = "Cell_types", values_to = "Cell_count")
Comp_long <- Comp_long[,3:5]

Comp_long$Cell_types[Comp_long$Cell_types == "Num.FOXP3"] <- "FOXP3"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.CD163"] <- "CD163"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.CD79A"] <- "CD79A"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.Non_immune.cells"] <- "Others"
Comp_long$Cell_types[Comp_long$Cell_types == "Num.CD8A"] <- "CD8A"

fill_colors <- c("FOXP3" = "#8302a2", "CD163" = "pink", "CD79A" = "#d60f23", "Others" = "grey", "CD8A" = "black")

g <- ggbarplot(Comp_long, x = "Patient_ID" , y = "Cell_count", fill= "Cell_types", color = "Cell_types", 
               order= c("77.3","734.5", "1039.4","749.1"), add = c("mean"),
               ylim= c(0,750), ylab = "Cell count",
               xlab= "Patient ID",
               palette = fill_colors,legend= "bottom", legend.title = "" , title = "Cells expressing GZMK")
g
pdf(paste0("Comp_GZMK", ".pdf"))
print(g)
dev.off()




pg <- ggarrange(a, b, p, g  , labels = c("A", "B", "C", "D") )

pdf(paste0("Test_marker.pdf"), width = 10 , height = 11)
print(pg)
dev.off()

 
