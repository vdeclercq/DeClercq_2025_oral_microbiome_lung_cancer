---
  title: "LungCancer_Study"
author: "Vanessa DeClercq"
date: "2025-12-01"
output: html_document
---
  
## R Markdown for Oral Microbiome Lung Cancer Study
  
#Read in files
#read metadata files
metadata <- read_csv("metadata.csv")
rownames(metadata) <- metadata$Saliva_sampleID
View(metadata) 

#check that variables are ready as numeric or string variables
print(sapply(metadata, class))


```

## Metadata Analysis
-normality test 
-Descriptive analysis (counts, mean, median, SD, IQR) of covariates - sex, AGE, BMI, ETC.

Examples:
  mean(metadata$age_at_surgery)
sd(metadata$age_at_surgery)
median(metadata$age_at_surgery)
IQR(metadata$age_at_surgery)
range(metadata$age_at_surgery)

```{r}
#check for normality
shapiro.test(metadata$age_at_surgery)

#counts of each gp
table(metadata$stage_cat)

#summary table of participants age by sex
library(PMCMRplus)
aggregate(age_at_surgery ~ Sex, data = metadata, summary)
kruskal.test(metadata$age_at_surgery, 
             metadata$Sex)
kwAllPairsDunnTest(metadata$age_at_surgery, 
                   as.factor(metadata$Sex), "BH")
t.test(metadata$age_at_surgery~metadata$Sex)

mean(metadata$age_at_surgery)
range(metadata$age_at_surgery)

#TABLES FOR COUNTS, PERCENTAGES - CATEGORICAL VARIABLES
#install.packages("gmodels")
library(gmodels)
library(reporttools)

Smoking_Status <- CrossTable(metadata$Smoking_Status, 
                             metadata$Sex, 
                             prop.r=TRUE, chisq = TRUE)
pairwise.fisher.test(metadata$SampleType, metadata$stage_group,
                     p.adjust.method = "BH")

Stage <- CrossTable(metadata$stage_cat, 
                    metadata$Sex, 
                    prop.r=TRUE, chisq = TRUE)

Substype <- CrossTable(metadata$histological_type_cat, 
                       metadata$Sex, 
                       prop.r=TRUE, chisq = TRUE)

PostTreamtment <- CrossTable(metadata$therapy_after, 
                             metadata$Sex, 
                             prop.r=TRUE, chisq = TRUE)

PDL1 <- CrossTable(metadata$PDL1, 
                   metadata$Sex, 
                   prop.r=TRUE, chisq = TRUE)

mutations <- CrossTable(metadata$`Mutations list`,
                        metadata$Sex, 
                        prop.r = TRUE, chisq = TRUE)

mutations2 <- CrossTable(metadata$Mutations,
                         metadata$Sex, 
                         prop.r = TRUE, chisq = TRUE)

EFS3mon <- CrossTable(metadata$EFS_3mon, 
                      metadata$Sex, 
                      prop.r=TRUE, chisq = TRUE)

EFS6mon <- CrossTable(metadata$EFS_6mon, 
                      metadata$Sex, 
                      prop.r=TRUE, chisq = TRUE)

EFS1yr <- CrossTable(metadata$EFS_1yr, 
                     metadata$Sex, 
                     prop.r=TRUE, chisq = TRUE)

abx_6mon <- CrossTable(metadata$abx_6mon_prior, 
                       metadata$Sex, 
                       prop.r=TRUE, chisq = TRUE)

abx_3mon <- CrossTable(metadata$abx_3mon_prior, 
                       metadata$Sex, 
                       prop.r=TRUE, chisq = TRUE)



```

##Diversity Analysis

#Alpha Diversity
Use QIIME2 generated alpha diversity metric. 
First, test for normality. Then perform parametric or non-parametric test for differences in alpha diversity among sample types. 
Four different alpha diversity metrics will be used; richness, evenness, Faiths PD, Shannon diversity.
```{r}
#Reading Alpha Diversity file
alpha <- read_csv("alpha_diversity.csv")
rownames(alpha) <- alpha$...1
View(alpha)
#rename column
colnames(alpha) [colnames(alpha) == '...1'] <- 'SampleID'

#match alpha diversity file to metadata file
intersect(rownames(metadata),rownames(alpha))
samples_to_keep <- intersect(rownames(metadata), rownames(alpha))
filtered_alpha <- data.frame(alpha[samples_to_keep,])
rownames(filtered_alpha) <- filtered_alpha$SampleID
View(filtered_alpha)

intersect(rownames(filtered_alpha),rownames(metadata))
samples_to_keep2 <- intersect(rownames(filtered_alpha),rownames(metadata))
filtered_metadata <- data.frame(metadata[samples_to_keep2,])
rownames(filtered_metadata) <- filtered_metadata$Saliva_sampleID
View(filtered_metadata)


#test for normality.
for (i in 2:5){print(shapiro.test(filtered_alpha[,i]))}
#results show none are normal so I will have to use non-parametric tests moving forward.

#I don't want the 'unknown' response from the 'stage_cat'column. Replace the unknowns with blanks.
filtered_metadata$abx_6mon_prior <-replace(filtered_metadata$abx_6mon_prior, filtered_metadata$abx_6mon_prior=="n/a", NA)
filtered_metadata$abx_3mon_prior <-replace(filtered_metadata$abx_3mon_prior, filtered_metadata$abx_3mon_prior=="n/a", NA)

#use Kruska wallis test to test for significant differences bewteen sample types.
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$Sex))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$Smoking._tatus))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$stage_cat))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$histological_type_cat))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$therapy_after))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$adjuvant))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$PDL1))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$Mutations))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$EFS_3mon))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$EFS_6mon))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$EFS_1yr))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$abx_6mon_prior))}
for (i in 2:5){print(kruskal.test(filtered_alpha[,i], filtered_metadata$abx_3mon_prior))}


#post-hoc test for between group differences.
library(PMCMRplus)
for (i in 2:5){print(kwAllPairsDunnTest(filtered_alpha[,i], as.factor(filtered_metadata$stage_cat), "BH"))}


#create box plot for alpha diversity measures
#boxplot with ggplot
#merge the 'filtered_alpha' and 'filtered_metadata' dataframes.
alpha_metadata <- merge(filtered_metadata, filtered_alpha, by = 'row.names', all=TRUE)
alpha_metadata <- data.frame(alpha_metadata, row.names = 1)
View(alpha_metadata)

#ggplot
library(ggplot2)

ggplot(data=subset(alpha_metadata, !is.na(stage_cat)), aes(x=stage_cat, y=observed_features, fill = stage_cat)) +geom_point() + 
  ylim(c(0,370)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
  labs(x = "Cancer Stage", y = "Observed Features") +
  labs(fill = "Cancer Stage") +
  scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) +
  theme_bw() +
  ggtitle("Observed Features") +
  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
  theme(axis.title.x = element_text(size =16, face = "bold"))+
  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
  theme(axis.text = element_text(size =14)) +
  theme(legend.title = element_text(size = 14, face = "bold"))+
  theme(legend.text = element_text(size = 14))
ggsave("observed_features_STAGE.pdf",  width = 6, height = 4)
ggsave("observed_features.png")



#correlation with alpha diversity
plot(filtered_metadata$age_at_surgery, filtered_alpha$observed_features, 
     xlab="Age (years)", ylab="Observed Features", 
     xlim = c(50, 90), ylim = c(0, 400), 
     pch=19, cex.axis = 1, cex.lab = 1.5, las=1)
abline(lm(filtered_alpha$observed_features ~ filtered_metadata$age_at_surgery),
       col = "blue", lwd=7)
test <- cor.test(filtered_alpha$observed_features, filtered_metadata$age_at_surgery, 
                 method=c("pearson"))    
test
test2 <- cor.test(filtered_alpha$observed_features, filtered_metadata$age_at_surgery,
                  method=c("spearman"))    
test2


```

#Beta Diversity
#use of QIIME2 generated Bray Curtis dissimilarity matrix or Weighted UniFrac distance matrix.
# Bray Curtis
```{r}
#read in beta diversity files
bray_curtis <- read_csv("bray_curtis.csv")
bray_curtis <- as.data.frame(bray_curtis)
rownames(bray_curtis) <- bray_curtis[,1]
bray_curtis <- as.data.frame(bray_curtis[,-c(1:1)])
View(bray_curtis)

weighted_unifrac <- read_csv("weighted_unifrac.csv")
weighted_unifrac <- as.data.frame(weighted_unifrac)
rownames(weighted_unifrac) <- weighted_unifrac[,1]
weighted_unifrac <- as.data.frame(weighted_unifrac[,-c(1:1)])
View(weighted_unifrac)

jaccard <- read_csv("jaccard.csv")
jaccard <- as.data.frame(jaccard)
rownames(jaccard) <- jaccard[,1]
jaccard <- as.data.frame(jaccard[,-c(1:1)])
View(jaccard)

unweighted_unifrac <- read_csv("unweighted_unifrac.csv")
unweighted_unifrac <- as.data.frame(unweighted_unifrac)
rownames(unweighted_unifrac) <- unweighted_unifrac[,1]
unweighted_unifrac <- as.data.frame(unweighted_unifrac[,-c(1:1)])
View(unweighted_unifrac)


#match beta diversity data to metadata
intersect(rownames(filtered_metadata),rownames(bray_curtis))
samples_to_keep <- intersect(rownames(filtered_metadata),rownames(bray_curtis))
filtered_bray_curtis <- as.data.frame(bray_curtis[samples_to_keep,])
rownames(filtered_bray_curtis) <- filtered_bray_curtis[,1]
View(filtered_bray_curtis)

filtered_bray_curtis <- t(filtered_bray_curtis)
intersect(rownames(filtered_metadata),rownames(filtered_bray_curtis))
samples_to_keep <- intersect(rownames(filtered_metadata),rownames(filtered_bray_curtis))
filtered_bray_curtis <- as.data.frame(filtered_bray_curtis[samples_to_keep,])
rownames(filtered_bray_curtis) <- filtered_bray_curtis[,1]
View(filtered_bray_curtis)

intersect(rownames(bray_curtis),rownames(filtered_metadata))
samples_to_keep2 <- intersect(rownames(bray_curtis),rownames(filtered_metadata))
filtered_metadata_bray <- as.data.frame(filtered_metadata[samples_to_keep2,])
rownames(filtered_metadata_bray) <- filtered_metadata_bray[,1]
View(filtered_metadata_bray)


#adonis test - permutational ANOVA of dissimilarities
library(vegan)

set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$age_at_surgery,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$age_at_surgery,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$age_at_surgery,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$age_at_surgery,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$Sex,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$Sex,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$Sex,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$Sex,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$Smoking_Status,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$Smoking_Status,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$Smoking_Status,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$Smoking_Status,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$stage_cat,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$stage_cat,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$stage_cat,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$stage_cat,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$histological_type_cat,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$histological_type_cat,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$histological_type_cat,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$histological_type_cat,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$PDL1,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$PDL1,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$PDL1,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$PDL1,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$adjuvant,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$adjuvant,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$adjuvant,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$adjuvant,
        permutations = 1000, by="margin")



set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$Mutations,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$Mutations,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$Mutations,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$Mutations,na.action=na.exclude,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$EFS_3mon,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$EFS_3mon,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$EFS_3mon,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$EFS_3mon,na.action=na.exclude,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$EFS_6mon,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$EFS_6mon,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$EFS_6mon,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$EFS_6mon,na.action=na.exclude,
        permutations = 1000, by="margin")

set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$EFS_1yr,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$XEFS_1yr,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$EFS_1yr,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$EFS_1yr,na.action=na.exclude,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$abx_6mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$abx_6mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$abx_6mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$abx_6mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")


set.seed(23)
adonis2(filtered_bray_curtis ~ filtered_metadata$abx_3mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_weighted_unifrac ~ filtered_metadata$abx_3mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_jaccard ~ filtered_metadata$abx_3mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")
set.seed(23)
adonis2(filtered_unweighted_unifrac ~ filtered_metadata$abx_3mon_prior,na.action=na.exclude,
        permutations = 1000, by="margin")


#plot beta diversity
library(devtools)
library(ggordiplots)
library(ggplot2)
bray_curtis_pcoa <- cmdscale(filtered_bray_curtis, k=2, eig = TRUE)
barplot(bray_curtis_pcoa$eig[1:10])
component1 <- bray_curtis_pcoa$eig[1]/sum(bray_curtis_pcoa$eig)
component2 <- bray_curtis_pcoa$eig[2]/sum(bray_curtis_pcoa$eig)
component1*100
component2*100

plot_data <- data.frame(pc1=bray_curtis_pcoa$points[ ,1],
                        pc2=bray_curtis_pcoa$points[ ,2],
                        Sample_type=filtered_metadata_bray$stage_cat)

bray_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +  geom_point(size=2) +
  stat_ellipse(linewidth=2) + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), 
                                                 name="Stage", labels=c("I", "II", "III"))+
  ggtitle("Bray-Curtis") +
  theme_bw() + xlab("PC1(4.6%)") + ylab("PC2(3.3%)") + theme(legend.position = c(0.88,0.84)) +
  theme(text = element_text(size = 26)) + ylim(-0.5, 0.5)+ xlim(-0.5, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
bray_plot
ggsave("bray_curtis_stage.pdf",  width = 9, height = 7)


plot_data <- data.frame(pc1=bray_curtis_pcoa$points[ ,1],
                        pc2=bray_curtis_pcoa$points[ ,2],
                        Sample_type=filtered_metadata_bray$EFS_6mon)
bray_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +  geom_point(size=2) +
  stat_ellipse(linewidth=2) + scale_color_manual(values=c("#0072B2", "#D55E00"), 
                                                 name="EFS 6mon", labels=c("No", "Yes"))+
  ggtitle("Bray-Curtis") +
  theme_bw() + xlab("PC1(4.6%)") + ylab("PC2(3.3%)") + theme(legend.position = c(0.85,0.88)) +
  theme(text = element_text(size = 26)) + ylim(-0.5, 0.5)+ xlim(-0.5, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
bray_plot
ggsave("bray_curtis_EFS6mon.pdf",  width = 9, height = 7)


plot_data <- data.frame(pc1=bray_curtis_pcoa$points[ ,1],
                        pc2=bray_curtis_pcoa$points[ ,2],
                        Sample_type=filtered_metadata_bray$Smoking_Status)
bray_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +  geom_point(size=2) +
  stat_ellipse(linewidth=2) + scale_color_manual(values=c("#332288", "#882255", "#6699CC"), 
                                                 name="Smoking", labels=c("Current", "Former", "Never"))+
  ggtitle("Bray-Curtis") +
  theme_bw() + xlab("PC1(4.6%)") + ylab("PC2(3.3%)") + theme(legend.position = c(0.85,0.85)) +
  theme(text = element_text(size = 26)) + ylim(-0.5, 0.5)+ xlim(-0.5, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
bray_plot
ggsave("bray_curtis_smoking.pdf",  width = 9, height = 7)


plot_data <- data.frame(pc1=bray_curtis_pcoa$points[ ,1],
                        pc2=bray_curtis_pcoa$points[ ,2],
                        Sample_type=filtered_metadata_bray$abx_3mon_prior)
bray_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +  geom_point(size=2) +
  stat_ellipse(linewidth=2) + scale_color_manual(values=c("darkgrey","lightpink"), 
                                                 name="Abx 3mon", labels=c("No", "Yes"))+
  ggtitle("Bray-Curtis") +
  theme_bw() + xlab("PC1(4.6%)") + ylab("PC2(3.3%)") + theme(legend.position = c(0.85,0.85)) +
  theme(text = element_text(size = 26)) + ylim(-0.5, 0.5)+ xlim(-0.5, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
bray_plot
ggsave("bray_curtis_abx3mon.pdf",  width = 9, height = 7)




#weighted UniFrac
weighted_unifrac_pcoa <- cmdscale(filtered_weighted_unifrac, k=2, eig = TRUE)
barplot(weighted_unifrac_pcoa$eig[1:10])
component1 <- weighted_unifrac_pcoa$eig[1]/sum(weighted_unifrac_pcoa$eig)
component2 <- weighted_unifrac_pcoa$eig[2]/sum(weighted_unifrac_pcoa$eig)
component1*100
component2*100

plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$stage_cat)
weighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), 
                                            name="Stage", labels=c("I", "II", "III"))+
  ggtitle("Weighted_uniFrac") +
  theme_bw() + xlab("PC1(45.6%)") + ylab("PC2(22.1%)") + theme(legend.position = c(0.89,0.86)) +
  theme(text = element_text(size = 26)) + ylim(-0.5, 0.5)+ xlim(-0.5, 0.6)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
weighted_unifrac_plot
ggsave("Stage_weighted_unifrac_stage.pdf",  width = 9, height = 7)
ggsave("Stage_weighted_unifrac_plot.png")

weighted_unifrac_pcoa <- cmdscale(filtered_weighted_unifrac, k=2, eig = TRUE)
barplot(weighted_unifrac_pcoa$eig[1:10])
component1 <- weighted_unifrac_pcoa$eig[1]/sum(weighted_unifrac_pcoa$eig)
component2 <- weighted_unifrac_pcoa$eig[2]/sum(weighted_unifrac_pcoa$eig)
component1*100
component2*100



plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$X6mon_EFS)
weighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#0072B2", "#D55E00"), 
                                            name="EFS 6mon", labels=c("No", "Yes"))+
  ggtitle("Weighted_uniFrac") +
  theme_bw() + xlab("PC1(45.6%)") + ylab("PC2(22.1%)") + theme(legend.position = c(0.85,0.88)) +
  theme(text = element_text(size = 26)) + ylim(-0.3, 0.3)+ xlim(-0.7, 0.7)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
weighted_unifrac_plot
ggsave("6mon_EFS_weighted_unifrac.pdf",  width = 9, height = 7)
ggsave("ESF_weighted_unifrac_plot.png")


plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$Sex)
weighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#66A61E", "#E6AB02"), 
                                            name="EFS", labels=c("No", "Yes"))+
  ggtitle("Weighted_uniFrac") +
  theme_bw() + xlab("PC1(45.6%)") + ylab("PC2(22.1%)") + theme(legend.position = c(0.85,0.88)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.4)+ xlim(-0.5, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
weighted_unifrac_plot
ggsave("Sex_weighted_unifrac.pdf",  width = 9, height = 7)
ggsave("Sex_weighted_unifrac_plot.png")


plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$Smoking.Status)
weighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#332288", "#882255", "#6699CC"), 
                                            name="Smoking Status", labels=c("Current", "Former", "Never"))+
  ggtitle("Weighted_uniFrac") +
  theme_bw() + xlab("PC1(45.6%)") + ylab("PC2(22.1%)") + theme(legend.position = c(0.82,0.86)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.4)+ xlim(-0.5, 0.7)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
weighted_unifrac_plot
ggsave("smoking_weighted_unifrac.pdf",  width = 9, height = 7)
ggsave("smoking_weighted_unifrac_plot.png")


plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$abx_3mon_prior)
weighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("darkgreen", "darkblue","darkorange"), 
                                            name="Antibiotic Use", labels=c("No", "Yes"))+
  ggtitle("Weighted_uniFrac") +
  theme_bw() + xlab("PC1(45.6%)") + ylab("PC2(22.1%)") + theme(legend.position = c(0.82,0.86)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.4)+ xlim(-0.5, 0.7)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
weighted_unifrac_plot
ggsave("abx3mon_weighted_unifrac.pdf",  width = 9, height = 7)

plot_data <- data.frame(pc1=weighted_unifrac_pcoa$points[ ,1],
                        pc2=weighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$age_at_surgery)
weighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=4) +
  scale_color_gradient(low = "red" , high= "blue", name = "Age") +
  ggtitle("Weighted_uniFrac") +
  theme_bw() + xlab("PC1(45.6%)") + ylab("PC2(22.1%)") +
  theme(legend.position = c(0.9,0.84)) +
  theme(text = element_text(size = 26)) + ylim(-0.32, 0.3)+ xlim(-0.4, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
weighted_unifrac_plot
ggsave("age_weighted_unifrac.pdf",  width = 9, height = 7)
ggsave("age_weighted_unifrac_plot.png")


#Jaccard plots
jaccard_pcoa <- cmdscale(filtered_jaccard, k=2, eig = TRUE)
barplot(jaccard_pcoa$eig[1:10])
component1 <- jaccard_pcoa$eig[1]/sum(jaccard_pcoa$eig)
component2 <- jaccard_pcoa$eig[2]/sum(jaccard_pcoa$eig)
component1*100
component2*100

plot_data <- data.frame(pc1=jaccard_pcoa$points[ ,1],
                        pc2=jaccard_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$stage_cat)

jaccard_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), 
                                            name="Stage", labels=c("I", "II", "III"))+
  ggtitle("Jaccard Distance") +
  theme_bw() + xlab("PC1(2.7%)") + ylab("PC2(2.5%)") + theme(legend.position = c(0.88,0.84)) +
  theme(text = element_text(size = 26)) + ylim(-0.5, 0.5)+ xlim(-0.5, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
jaccard_plot
ggsave("stage_jaccard.pdf",  width = 9, height = 7)
ggsave("stage_jaccard.png")


#unweighted UniFrac
unweighted_unifrac_pcoa <- cmdscale(filtered_unweighted_unifrac, k=2, eig = TRUE)
barplot(unweighted_unifrac_pcoa$eig[1:10])
component1 <- unweighted_unifrac_pcoa$eig[1]/sum(unweighted_unifrac_pcoa$eig)
component2 <- unweighted_unifrac_pcoa$eig[2]/sum(unweighted_unifrac_pcoa$eig)
component1*100
component2*100


plot_data <- data.frame(pc1=unweighted_unifrac_pcoa$points[ ,1],
                        pc2=unweighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$stage_cat)
unweighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), 
                                            name="Stage", labels=c("I", "II", "III"))+
  ggtitle("Unweighted_UniFrac") +
  theme_bw() + xlab("PC1(26.3%)") + ylab("PC2(10.5%)") + theme(legend.position = c(0.89,0.86)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.4)+ xlim(-0.6, 0.6)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
unweighted_unifrac_plot
ggsave("Stage_unweighted_unifrac.pdf",  width = 9, height = 7)
ggsave("Stage_unweighted_unifrac_plot.png")


plot_data <- data.frame(pc1=unweighted_unifrac_pcoa$points[ ,1],
                        pc2=unweighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$X6mon_EFS)
unweighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#0072B2", "#D55E00"), 
                                            name="EFS", labels=c("No", "Yes"))+
  ggtitle("Unweighted_UniFrac") +
  theme_bw() + xlab("PC1(26.3%)") + ylab("PC2(10.5%)") + theme(legend.position = c(0.85,0.88)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.4)+ xlim(-0.6, 0.6)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
unweighted_unifrac_plot
ggsave("EFS6mon_unweighted_unifrac.pdf",  width = 9, height = 7)
ggsave("ESF_unweighted_unifrac_plot.png")


plot_data <- data.frame(pc1=unweighted_unifrac_pcoa$points[ ,1],
                        pc2=unweighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$Sex)
unweighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#66A61E", "#E6AB02"), 
                                            name="Sex", labels=c("Female", "Male"))+
  ggtitle("Unweighted_UniFrac") +
  theme_bw() + xlab("PC1(26.3%)") + ylab("PC2(10.5%)") + theme(legend.position = c(0.85,0.88)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.4)+ xlim(-0.5, 0.5)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
unweighted_unifrac_plot
ggsave("Sex_unweighted_unifrac.pdf",  width = 9, height = 7)
ggsave("Sex_unweighted_unifrac_plot.png")


plot_data <- data.frame(pc1=unweighted_unifrac_pcoa$points[ ,1],
                        pc2=unweighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$Smoking.Status)
unweighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=2) +
  stat_ellipse(size=2) + scale_color_manual(values=c("#332288", "#882255", "#6699CC"), 
                                            name="Smoking Status", labels=c("Current", "Former", "Never"))+
  ggtitle("Unweighted_UniFrac") +
  theme_bw() + xlab("PC1(26.3%)") + ylab("PC2(10.5%)") + theme(legend.position = c(0.82,0.86)) +
  theme(text = element_text(size = 26)) + ylim(-0.4, 0.5)+ xlim(-0.6, 0.8)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
unweighted_unifrac_plot
ggsave("smoking_unweighted_unifrac.pdf",  width = 9, height = 7)
ggsave("smoking_unweighted_unifrac_plot.png")


plot_data <- data.frame(pc1=unweighted_unifrac_pcoa$points[ ,1],
                        pc2=unweighted_unifrac_pcoa$points[ ,2],
                        Sample_type=filtered_metadata$age_at_surgery)
unweighted_unifrac_plot <- ggplot(data=subset(plot_data, !is.na(Sample_type)), aes(x=pc1, y=pc2, color=Sample_type)) +
  geom_point(size=4) +
  scale_color_gradient(low = "red" , high= "blue", name = "Age") +
  ggtitle("Unweighted_UniFrac") +
  theme_bw() + xlab("PC1(26.3%)") + ylab("PC2(10.5%)") +
  theme(legend.position = c(0.9,0.84)) +
  theme(text = element_text(size = 26)) + ylim(-0.32, 0.3)+ xlim(-0.3, 0.4)+
  theme(axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)))
ggsave("age_unweighted_unifrac.pdf",  width = 9, height = 7)
ggsave("age_unweighted_unifrac_plot.png")


```


##TAXA ANALYSIS


```{r}
#preprep data for analysis.
#create taxa table (remove any extra columns & set row names)
library(readr)
taxa_table <- read_csv("feature-table_w_tax_FLIP.csv")
taxa_table <- as.data.frame(taxa_table)
rownames(taxa_table) <- taxa_table[,1]
taxa_table <- as.data.frame(taxa_table[,-c(1:1)])
View(taxa_table)

#create taxa table with only samples of interest, then flip the taxa table
#(filter taxa to match samples in the 'metadata' file).
intersect(rownames(taxa_table), rownames(metadata))
samples_to_keep <- intersect(rownames(taxa_table), rownames(metadata))
filtered_taxa <- taxa_table[samples_to_keep,]
View(filtered_taxa)

filtered_metadata <- metadata[rownames(filtered_taxa),]
filtered_metadata <- as.data.frame(filtered_metadata)
rownames(filtered_metadata) <- filtered_metadata[,1]
View(filtered_metadata)

#flip table so taxa are rows, samples columns
rare_filter_table <- t(filtered_taxa)
rare_filter_table <- rare_filter_table[-c(109:110),]
View(rare_filter_table)

#(filter taxa to match samples in the 'metadata_abx' file).
intersect(rownames(taxa_table), rownames(metadata_abx))
samples_to_keep_abx <- intersect(rownames(taxa_table), rownames(metadata_abx))
filtered_taxa_abx <- taxa_table[samples_to_keep_abx,]
View(filtered_taxa_abx)

filtered_metadata_abx <- metadata[rownames(filtered_taxa_abx),]
filtered_metadata_abx <- as.data.frame(filtered_metadata_abx)
rownames(filtered_metadata_abx) <- filtered_metadata_abx[,1]
View(filtered_metadata_abx)

rare_filter_table_abx <- t(filtered_taxa_abx)
rare_filter_table_abx <- rare_filter_table_abx[-c(109:110),]
View(rare_filter_table_abx)


#(filter taxa to match samples in the 'metadata_EFS6mon' file).
intersect(rownames(taxa_table), rownames(metadata_EFS6mon))
samples_to_keep_EFS6mon <- intersect(rownames(taxa_table), rownames(metadata_EFS6mon))
filtered_taxa_EFS6mon <- taxa_table[samples_to_keep_EFS6mon,]
View(filtered_taxa_EFS6mon)

filtered_metadata_EFS6mon <- metadata[rownames(filtered_taxa_EFS6mon),]
filtered_metadata_EFS6mon <- as.data.frame(filtered_metadata_EFS6mon)
rownames(filtered_metadata_EFS6mon) <- filtered_metadata_EFS6mon[,1]
View(filtered_metadata_EFS6mon)

rare_filter_table_EFS6mon <- t(filtered_taxa_EFS6mon)
rare_filter_table_EFS6mon <- rare_filter_table_EFS6mon[-c(109:110),]
View(rare_filter_table_EFS6mon)

#(filter taxa to match samples in the 'metadata_EFS1yr' file).
intersect(rownames(taxa_table), rownames(metadata_EFS1yr))
samples_to_keep_EFS1yr <- intersect(rownames(taxa_table), rownames(metadata_EFS1yr))
filtered_taxa_EFS1yr <- taxa_table[samples_to_keep_EFS1yr,]
View(filtered_taxa_EFS1yr)

filtered_metadata_EFS1yr <- metadata[rownames(filtered_taxa_EFS1yr),]
filtered_metadata_EFS1yr <- as.data.frame(filtered_metadata_EFS1yr)
rownames(filtered_metadata_EFS1yr) <- filtered_metadata_EFS1yr[,1]
View(filtered_metadata_EFS1yr)

rare_filter_table_EFS1yr <- t(filtered_taxa_EFS1yr)
rare_filter_table_EFS1yr <- rare_filter_table_EFS1yr[-c(109:110),]
View(rare_filter_table_EFS1yr)

#MaAsLin2 DA analysis
library("devtools")
install_github("biobakery/maaslin3")
library(maaslin3)
install.packages("gtable")

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_Sex',
                    formula = '~ Sex',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_age',
                    formula = '~ age_at_surgery',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_smoking',
                    formula = '~ Smoking_Status',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('Smoking_Status,Current'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_stage',
                    formula = '~ stage_cat',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('stage_cat,I'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_stage_sex',
                    formula = '~ stage_cat + Sex',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('stage_cat,I; Sex,Female'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_subtype',
                    formula = '~ histological_type_cat',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('histological_type_cat,adenocarcinoma'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_PDL1',
                    formula = '~ PDL1',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('PDL1,<1%'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa,
                    input_metadata = filtered_metadata,
                    output = 'Maalin3_output_EFS3mon',
                    formula = '~ EFS_3mon',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('EFS_3mon,NO'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa_EFS6mon,
                    input_metadata = filtered_metadata_EFS6mon,
                    output = 'Maalin3_output_EFS6mon',
                    formula = '~ EFS_6mon',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('EFS_6mon,NO'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa_EFS1yr,
                    input_metadata = filtered_metadata_EFS1yr,
                    output = 'Maalin3_output_EFS1yr',
                    formula = '~ EFS_1yr',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('EFS_1yr,NO'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa_abx,
                    input_metadata = filtered_metadata_abx,
                    output = 'Maalin3_output_abx3mon',
                    formula = '~ abx_3mon_prior',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('abx_3mon_prior,No'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)

set.seed(1)
fit_out <- maaslin3(input_data = filtered_taxa_abx,
                    input_metadata = filtered_metadata_abx,
                    output = 'Maalin3_output_abx6mon',
                    formula = '~ abx_6mon_prior',
                    normalization = 'TSS',
                    transform = 'LOG',
                    reference = c('abx_6mon_prior,No'),
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 250,
                    cores = 1)



#ANCOM2 DA analysis
library(compositions)
library(tidyr)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(exactRankTests)
deps = c("exactRankTests", "nlme", "dplyr", "ggplot2", "compositions")
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(dep)
  }
  library(dep, character.only = TRUE)
}
source("~/Desktop/ANCOM2.R")



#preprocessing step - to deal with different types of zeros before performing differential abundance analysis 
preprocess <- feature_table_pre_process(feature_table = rare_filter_table, 
                                        meta_data = filtered_metadata, sample_var = "Saliva_sampleID", 
                                        out_cut = 0.05,zero_cut = 0.9, lib_cut = 1000,)

#run main ANCOM function with preprocessed data,adjusting p-values for multiple comparisons
# without covariates, use NULL
rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"age_at_surgery", "BH", 0.1, NULL)
View(rez[[1]])

rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"Sex", "BH", 0.1, NULL)
View(rez[[1]])

rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"Smoking_Status", "BH", 0.1, NULL)
View(rez[[1]])

rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"stage_cat", "BH", 0.1, NULL)
View(rez[[1]])

rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"stage_cat", "BH", 0.1, "stage_cat + Sex")
View(rez[[1]])

rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"histological_type_cat", "BH", 0.1, NULL)
View(rez[[1]])

rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"PDL1", "BH", 0.1, NULL)
View(rez[[1]])


rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"EFS_3mon", "BH", 0.1, NULL)
View(rez[[1]])


#for specific datasets with missing data
preprocess <- feature_table_pre_process(feature_table = rare_filter_table_EFS6mon, 
                                        meta_data = filtered_metadata_EFS6mon, sample_var = "Saliva_sampleID", 
                                        out_cut = 0.05,zero_cut = 0.9, lib_cut = 1000,)
rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"EFS_6mon", "BH", 0.1, NULL)
View(rez[[1]])

preprocess <- feature_table_pre_process(feature_table = rare_filter_table_EFS1yr, 
                                        meta_data = filtered_metadata_EFS1yr, sample_var = "Saliva_sampleID", 
                                        out_cut = 0.05,zero_cut = 0.9, lib_cut = 1000,)
rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"EFS_1yr", "BH", 0.1, NULL)
View(rez[[1]])


preprocess <- feature_table_pre_process(feature_table = rare_filter_table_abx, 
                                        meta_data = filtered_metadata_abx, sample_var = "Saliva_sampleID", 
                                        out_cut = 0.05,zero_cut = 0.9, lib_cut = 1000,)
rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"abx_3mon_prior", "BH", 0.1, NULL)
View(rez[[1]])

rez <- ANCOM(preprocess$feature_table, preprocess$meta_data, 
             preprocess$structure_zeros,"abx_6mon_prior", "BH", 0.1, NULL)
View(rez[[1]])




#Corncob differential abundance (DA) analysis
library(corncob)
library(phyloseq)
library(magrittr)
#create the object containing the metadata and taxa to be tested.
otu_tab <- otu_table(rare_filter_table,taxa_are_rows = TRUE)
sam_data <- sample_data(filtered_metadata)
phylo <- phyloseq(otu_tab, sam_data)

#runs corncob DA analysis, returns plot and, DA results
results <- differentialTest(formula = ~ Sex,
                            phi.formula = ~ Sex,
                            formula_null = ~ 1,
                            phi.formula_null = ~ Sex,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)
plot(results)
results$p
results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)

results <- differentialTest(formula = ~ age_at_surgery,
                            phi.formula = ~ age_at_surgery,
                            formula_null = ~ 1,
                            phi.formula_null = ~ age_at_surgery,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1,
                            full_output = TRUE)
plot(results)
results$p
results$p_fdr
results$significant_taxa
results$significant_models
results$all_models
which(results$p <0.05)
results$all_models[7]
results$all_models[36]
results$all_models[37]
results$all_models[53]
results$all_models[61]
results$all_models[69]
results$all_models[87]
results$all_models[90]
results$p[90]


results <- differentialTest(formula = ~ Smoking_Status,
                            phi.formula = ~ Smoking_Status,
                            formula_null = ~ 1,
                            phi.formula_null = ~ Smoking_Status,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)
results$p[6]
results$p[13]
results$p[55]
results$p[85]
results$all_models[6]
results$all_models[13]
results$all_models[55]
results$all_models[85]


results <- differentialTest(formula = ~ stage_cat,
                            phi.formula = ~ stage_cat,
                            formula_null = ~ 1,
                            phi.formula_null = ~ stage_cat,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)

#wtih covariates
results <- differentialTest(formula = ~ stage_cat + Sex,
                            phi.formula = ~ stage_cat + Sex,
                            formula_null = ~ 1,
                            phi.formula_null = ~ stage_cat + Sex,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)

results <- differentialTest(formula = ~ histological_type_cat,
                            phi.formula = ~ histological_type_cat,
                            formula_null = ~ 1,
                            phi.formula_null = ~ histological_type_cat,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)
results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)


results <- differentialTest(formula = ~ PDL1,
                            phi.formula = ~ PDL1,
                            formula_null = ~ 1,
                            phi.formula_null = ~ PDL1,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)


results <- differentialTest(formula = ~ EFS_3mon,
                            phi.formula = ~ EFS_3mon,
                            formula_null = ~ 1,
                            phi.formula_null = ~ EFS_3mon,
                            data = phylo,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)


results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)


### read in specific files for variables with missing/NA values

#(filter taxa to match samples in the metadata file).
intersect(rownames(filtered_metadata_EFS6mon), rownames(taxa_table))
samples_to_keep_EFS6mon <- intersect(rownames(filtered_metadata_EFS6mon), rownames(taxa_table))
filtered_taxa_EFS6mon <- taxa_table[samples_to_keep_EFS6mon,]
filtered_taxa_EFS6mon <- as.data.frame(filtered_taxa_EFS6mon[,-c(1:1)])
View(filtered_taxa_EFS6mon)

intersect(rownames(filtered_metadata_EFS1yr), rownames(taxa_table))
samples_to_keep_EFS1yr <- intersect(rownames(filtered_metadata_EFS1yr), rownames(taxa_table))
filtered_taxa_EFS1yr <- taxa_table[samples_to_keep_EFS1yr,]
View(filtered_taxa_EFS1yr)

intersect(rownames(filtered_metadata_abx), rownames(taxa_table))
samples_to_keep_abx <- intersect(rownames(filtered_metadata_abx), rownames(taxa_table))
filtered_taxa_abx <- taxa_table[samples_to_keep_abx,]
View(filtered_taxa_abx)



#flip table so taxa are rows, samples columns
rare_filter_table_EFS6mon <- t(filtered_taxa_EFS6mon)
rare_filter_table_EFS1yr <- t(filtered_taxa_EFS1yr)
rare_filter_table_abx <- t(filtered_taxa_abx)


#make phyloseq object
otu_tab_EFS6mon <- otu_table(rare_filter_table_EFS6mon,taxa_are_rows = TRUE)
sam_data_EFS6mon <- sample_data(filtered_metadata_EFS6mon)
phylo_EFS6mon <- phyloseq(otu_tab_EFS6mon, sam_data_EFS6mon)

otu_tab_EFS1yr <- otu_table(rare_filter_table_EFS1yr,taxa_are_rows = TRUE)
sam_data_EFS1yr <- sample_data(filtered_metadata_EFS1yr)
phylo_EFS1yr <- phyloseq(otu_tab_EFS1yr, sam_data_EFS1yr)

otu_tab_abx <- otu_table(rare_filter_table_abx,taxa_are_rows = TRUE)
sam_data_abx <- sample_data(filtered_metadata_abx)
phylo_abx <- phyloseq(otu_tab_abx, sam_data_abx)


#DA
results <- differentialTest(formula = ~ EFS_6mon,
                            phi.formula = ~ EFS_6mon,
                            formula_null = ~ 1,
                            phi.formula_null = ~ EFS_6mon,
                            data = phylo_EFS6mon,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)
results$p[19]
results$all_models[19]



results <- differentialTest(formula = ~ EFS_1yr,
                            phi.formula = ~ EFS_1yr,
                            formula_null = ~ 1,
                            phi.formula_null = ~ EFS_1yr,
                            data = phylo_EFS1yr,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

results$p[28]
results$all_models[28]
results$p[51]
results$all_models[51]


results <- differentialTest(formula = ~ abx_3mon_prior,
                            phi.formula = ~ abx_3mon_prior,
                            formula_null = ~ 1,
                            phi.formula_null = ~ abx_3mon_prior,
                            data = phylo_abx,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

results$p_fdr
results$significant_taxa
which(results$p <0.05)

results <- differentialTest(formula = ~ abx_6mon_prior,
                            phi.formula = ~ abx_6mon_prior,
                            formula_null = ~ 1,
                            phi.formula_null = ~ abx_6mon_prior,
                            data = phylo_abx,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.1)

plot(results)
results$p
results$p_fdr
results$significant_taxa
results$significant_models
which(results$p <0.05)


#ALDEx2 DA analysis
#install.packages("BiocManager")
#BiocManager::install("ALDEx2")
library(ALDEx2)

#create model containing the variables to be tested.
matrixmodel <- model.matrix(~Sex, filtered_metadata) 
View(matrixmodel)
#for uncorrected, remove covariates.

#generate Monte Carlo samples of the Dirichlet distribution, perform centred log-ratio transformation.
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
#calculates the expected values for each coefficient of a glm model on the data returned by aldex.clr
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

#
matrixmodel <- model.matrix(~age_at_surgery, filtered_metadata) 
View(matrixmodel)
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

matrixmodel <- model.matrix(~Smoking_Status, filtered_metadata) 
View(matrixmodel)
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

matrixmodel <- model.matrix(~stage_cat, filtered_metadata) 
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

matrixmodel <- model.matrix(~histological_type_cat, filtered_metadata) 
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

matrixmodel <- model.matrix(~PDL1, filtered_metadata) 
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)


matrixmodel <- model.matrix(~EFS_3mon, filtered_metadata) 
CLR <- aldex.clr(rare_filter_table, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

#
matrixmodel <- model.matrix(~EFS_6mon, filtered_metadata_EFS6mon) 
CLR <- aldex.clr(rare_filter_table_EFS6mon, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

matrixmodel <- model.matrix(~EFS_1yr, filtered_metadata_EFS1yr) 
CLR <- aldex.clr(rare_filter_table_EFS1yr, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

matrixmodel <- model.matrix(~abx_3mon_prior, filtered_metadata_abx) 
CLR <- aldex.clr(rare_filter_table_abx, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)

matrixmodel <- model.matrix(~abx_6mon_prior, filtered_metadata_abx) 
CLR <- aldex.clr(rare_filter_table_abx, matrixmodel, mc.samples = 128)
GLM <- aldex.glm(CLR, matrixmodel)
View(GLM)



```
#RA of specific taxa 
RA <- sweep(rare_filter_table, 2, colSums(rare_filter_table), "/")
colSums(RA)

flipRA <-data.frame(t(RA))
colnames(flipRA)
View(flipRA)
identical(rownames(flipRA), rownames(filtered_metadata))

#Box plot - RA of specific taxa against by cancer variables
library(ggplot2)
library(ggpubr) 

Smoking_Status. scale_color_manual(values=c("#332288", "#882255", "#6699CC"), name="Smoking Status", labels=c("Current", "Former", "Never"))+
  Stage.  scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3"))
Sex. ("#66A61E", "#E6AB02")
Substype. ("red", "#075AFF", "grey")"
 PDL1.  scale_fill_manual(values=c("#66A61E", "#E6AB02", "deeppink3")) +
EFS scale_color_manual(values=c("#0072B2", "#D55E00"), name="EFS", labels=c("No", "Yes"))
abx scale_color_manual(values=c("("darkgrey","lightpink")) ), labels=c("No", "Yes"))
 
#Haemophilus.
DF <- data.frame(flipRA[ ,19],filtered_metadata$stage_cat)
RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.stage_cat)), aes(x=filtered_metadata.stage_cat, y=flipRA...19., fill = filtered_metadata.stage_cat)) +geom_point() + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
  labs(x = "Cancer Stage", y = "Relative Abundance") +
  labs(fill = "Cancer Stage") +
  scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) +
                                theme_bw() +
                                  ggtitle("Haemophilus") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_stage_Haemophilus.pdf",  width = 6, height = 4)
                                
                                DF <- data.frame(flipRA[ ,19],filtered_metadata$stage_cat)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.stage_cat)), aes(x=filtered_metadata.stage_cat, y=flipRA...19., fill = filtered_metadata.stage_cat)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Cancer Stage", y = "Relative Abundance") +
                                  labs(fill = "Cancer Stage") +
                                  scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) +
                                  theme_bw() +
                                  ylim(0,0.3)
                                ggtitle("Haemophilus") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_stage_Haemophilus2.pdf",  width = 6, height = 4)
                                
                                
                                DF <- data.frame(flipRA[ ,19],filtered_metadata$EFS_3mon)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.EFS_3mon)), aes(x=filtered_metadata.EFS_3mon, y=flipRA...19., fill = filtered_metadata.EFS_3mon)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "EFS at 3 months", y = "Relative Abundance") +
                                  labs(fill = "EFS 3mon") +
                                  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
                                  theme_bw() +
                                  ggtitle("Haemophilus") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_EFS3mon_Haemophilus.pdf",  width = 6, height = 4)
                                
                                
                                DF <- data.frame(flipRA[ ,19],filtered_metadata$EFS_6mon)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.EFS_6mon)), aes(x=filtered_metadata.EFS_6mon, y=flipRA...19., fill = filtered_metadata.EFS_6mon)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "EFS at 6 months", y = "Relative Abundance") +
                                  labs(fill = "EFS 3mon") +
                                  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
                                  theme_bw() +
                                  ggtitle("Haemophilus") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_EFS6mon_Haemophilus.pdf",  width = 6, height = 4)
                                
                                #Cardiobacterium
                                DF <- data.frame(flipRA[ ,42],filtered_metadata$abx_3mon_prior)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.abx_3mon_prior)), aes(x=filtered_metadata.abx_3mon_prior, y=flipRA...42., fill = filtered_metadata.abx_3mon_prior)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Antibiotic Use in Past 3 Months", y = "Relative Abundance") +
                                  labs(fill = "Abx 3mon") +
                                  scale_fill_manual(values=c("darkgrey","lightpink")) +
                                  theme_bw() +
                                  ggtitle("Cardiobacterium") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_abx3mon_Cardiobacterium.pdf",  width = 6, height = 4)
                                
                                #Megasphera
                                DF <- data.frame(flipRA[ ,1],filtered_metadata$abx_3mon_prior)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.abx_3mon_prior)), aes(x=filtered_metadata.abx_3mon_prior, y=flipRA...1., fill = filtered_metadata.abx_3mon_prior)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Antibiotic Use in Past 3 Months", y = "Relative Abundance") +
                                  labs(fill = "Abx 3mon") +
                                  scale_fill_manual(values=c("darkgrey","lightpink")) +
                                  theme_bw() +
                                  ggtitle("Megasphera") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_abx3mon_Megasphera.pdf",  width = 6, height = 4)
                                
                                #Streptococcus
                                DF <- data.frame(flipRA[ ,6],filtered_metadata$abx_3mon_prior)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.abx_3mon_prior)), aes(x=filtered_metadata.abx_3mon_prior, y=flipRA...6., fill = filtered_metadata.abx_3mon_prior)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Antibiotic Use in Past 3 Months", y = "Relative Abundance") +
                                  labs(fill = "Abx 3mon") +
                                  scale_fill_manual(values=c("darkgrey","lightpink")) +
                                  theme_bw() +
                                  ggtitle("Streptococcus") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_abx3mon_Streptococcus.pdf",  width = 6, height = 4)
                                
                                #Catonella
                                DF <- data.frame(flipRA[ ,51],filtered_metadata$EFS_1yr)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.EFS_1yr)), aes(x=filtered_metadata.EFS_1yr, y=flipRA...51., fill = filtered_metadata.EFS_1yr)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "EFS at 1 year", y = "Relative Abundance") +
                                  labs(fill = "EFS 1yr") +
                                  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
                                  theme_bw() +
                                  ggtitle("Catonella") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_EFS1yr_Catonella.pdf",  width = 6, height = 4)
                                
                                #Parvimonas
                                DF <- data.frame(flipRA[ ,12],filtered_metadata$EFS_1yr)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.EFS_1yr)), aes(x=filtered_metadata.EFS_1yr, y=flipRA...12., fill = filtered_metadata.EFS_1yr)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "EFS at 1 year", y = "Relative Abundance") +
                                  labs(fill = "EFS 1yr") +
                                  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
                                  theme_bw() +
                                  ggtitle("Parvimonas") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_EFS1yr_Parvimonas.pdf",  width = 6, height = 4)
                                
                                #Gemella
                                DF <- data.frame(flipRA[ ,28],filtered_metadata$EFS_1yr)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.EFS_1yr)), aes(x=filtered_metadata.EFS_1yr, y=flipRA...28., fill = filtered_metadata.EFS_1yr)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "EFS at 1 year", y = "Relative Abundance") +
                                  labs(fill = "EFS 1yr") +
                                  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
                                  theme_bw() +
                                  ggtitle("Gemella") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_EFS1yr_Gemella.pdf",  width = 6, height = 4)
                                
                                DF <- data.frame(flipRA[ ,28],filtered_metadata$EFS_3mon)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.EFS_3mon)), aes(x=filtered_metadata.EFS_3mon, y=flipRA...28., fill = filtered_metadata.EFS_3mon)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "EFS at 1 year", y = "Relative Abundance") +
                                  labs(fill = "EFS 1yr") +
                                  scale_fill_manual(values=c("#0072B2", "#D55E00")) +
                                  theme_bw() +
                                  ggtitle("Gemella") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_EFS3mon_Gemella.pdf",  width = 6, height = 4)
                                
                                
                                #Solobacterium
                                DF <- data.frame(flipRA[ ,34],filtered_metadata$stage_cat)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.stage_cat)), aes(x=filtered_metadata.stage_cat, y=flipRA...34., fill = filtered_metadata.stage_cat)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Cancer Stage", y = "Relative Abundance") +
                                  labs(fill = "Cancer Stage") +
                                  scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) +
                                  theme_bw() +
                                  ggtitle("Solobacterium") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_stage_Solobacterium.pdf",  width = 6, height = 4)
                                
                                DF <- data.frame(flipRA[ ,34],filtered_metadata$stage_cat)
                                RA_plot2 <- ggplot(data=subset(DF, !is.na(filtered_metadata.stage_cat)), aes(x=filtered_metadata.stage_cat, y=flipRA...34., fill = filtered_metadata.stage_cat)) +geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Cancer Stage", y = "Relative Abundance") +
                                  labs(fill = "Cancer Stage") +
                                  scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) +
                                  theme_bw() +
                                  ylim(0,0.01)
                                ggtitle("Solobacterium") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_stage_Solobacterium2.pdf",  width = 6, height = 4)
                                
                                
                                #Peptostreptococcaceae_G1
                                DF <- data.frame(flipRA[ ,34],filtered_metadata$histological_type_cat)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.histological_type_cat)), aes(x=filtered_metadata.histological_type_cat, y=flipRA...34., fill = filtered_metadata.histological_type_cat)) +
                                  geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Histological Subtype", y = "Relative Abundance") +
                                  labs(fill = "Subtype") +
                                  scale_fill_manual(values=c("red", "#075AFF", "grey")) +
                                  theme_bw() +
                                  ggtitle("Peptostreptococcaceae G1") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))+
                                  theme(legend.position = "none")
                                ggsave("RA_stage_PeptostreptococcaceaeG1.pdf",  width = 6, height = 4)
                                
                                
                                #Bifidobacterium
                                DF <- data.frame(flipRA[ ,53],filtered_metadata$Smoking_Status)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.Smoking_Status)), aes(x=filtered_metadata.Smoking_Status, y=flipRA...53., fill = filtered_metadata.Smoking_Status)) +
                                  geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Smoking Status", y = "Relative Abundance") +
                                  labs(fill = "Smoking") +
                                  scale_fill_manual(values=c("#332288", "#882255", "#6699CC")) +
                                  theme_bw() +
                                  ggtitle("Bifidobacterium") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_smoking_Bifidobacterium.pdf",  width = 6, height = 4)
                                
                                #Rothia
                                DF <- data.frame(flipRA[ ,22],filtered_metadata$Smoking_Status)
                                RA_plot <- ggplot(data=subset(DF, !is.na(filtered_metadata.Smoking_Status)), aes(x=filtered_metadata.Smoking_Status, y=flipRA...22., fill = filtered_metadata.Smoking_Status)) +
                                  geom_point() + 
                                  stat_boxplot(geom ='errorbar', width = 0.5) +
                                  geom_boxplot(alpha=1, color=1, outlier.shape = NA) +
                                  labs(x = "Smoking Status", y = "Relative Abundance") +
                                  labs(fill = "Smoking") +
                                  scale_fill_manual(values=c("#332288", "#882255", "#6699CC")) +
                                  theme_bw() +
                                  ggtitle("Rothia") +
                                  theme(plot.title = element_text(size = rel(1.5), face = "bold")) +
                                  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
                                  theme(axis.title.x = element_text(size =16, face = "bold"))+
                                  theme(axis.title.y = element_text(size = 16, face = "bold", angle = 90))+
                                  theme(axis.text = element_text(size =14)) +
                                  theme(legend.title = element_text(size = 14, face = "bold"))+
                                  theme(legend.text = element_text(size = 14))
                                ggsave("RA_smoking_Rothia.pdf",  width = 6, height = 4)
                                
