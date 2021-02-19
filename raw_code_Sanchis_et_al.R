####

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GEOquery", "canvasXpress", "ggplot2",
                       "clinfun", "GGally", "factoextra"))

####
{
   library(DESeq2)
   library(GEOquery)
   library(canvasXpress)
   library(ggplot2)
   library(clinfun)
   library(GGally)
   library(factoextra)
}


####
data <- getGEO(GEO = "GSE152075")   #replace the text between [] with the GSE of your choice

head(data)
head(data$GSE152075_series_matrix.txt.gz@phenoData@data[,c(1,2,8,40,39,42)])

####

clindata <- data[["GSE152075_series_matrix.txt.gz"]]@phenoData@data
head(clindata[,c(1,2,8,40,39,42)])



####

raw_counts <- read.delim("C:/Users/File/Location/GSE152075_raw_counts_GEO.txt.gz", stringsAsFactors=FALSE, sep = " ")
head(raw_counts[,c(1:10)])

####

raw_counts <- as.matrix(raw_counts) 
rownames(clindata) <- clindata$title
all(rownames(clindata) %in% colnames(raw_counts))
all(colnames(raw_counts) %in% rownames(clindata))


colnames(clindata)[colnames(clindata) == "sequencing_batch:ch1"] <- "batch"
clindata$batch = as.factor(clindata$batch)
colnames(clindata)[colnames(clindata) == "n1_ct:ch1"] <- "ct"
colnames(clindata)[colnames(clindata) == "sars-cov-2 positivity:ch1"] <- "positivity"
clindata$positivity[clindata$positivity == "pos"] <- "COVID19"
clindata$positivity[clindata$positivity == "neg"] <- "HEALTHY"
clindata$positivity = as.factor(clindata$positivity)

head(clindata$positivity)


####

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = clindata,
                              design = ~ positivity + batch)

head(dds)


####
dds <- estimateSizeFactors(dds)
head(dds)


####

norm_counts <- counts(dds, normalized=TRUE)
head(norm_counts[,c(1:10)])

####


{
   clindata$viral_load <- clindata$ct
   clindata$viral_load[clindata$viral_load == "N/A"] <- "Negative"
   clindata$viral_load[clindata$viral_load > 24 & clindata$viral_load != "Unknown" & clindata$viral_load != "Negative"] <- "LOW"
   clindata$viral_load[clindata$viral_load <= 24 & clindata$viral_load >= 19] <- "MEDIUM"
   clindata$viral_load[clindata$viral_load < 19] <- "HIGH"
   clindata$viral_load <- as.factor(clindata$viral_load)
   clindata$viral_load <- factor(clindata$viral_load, levels = c("Negative", "LOW", "MEDIUM", "HIGH", "Unknown"))
   
   clindata$positivity <- factor(clindata$positivity, levels = c("HEALTHY", "COVID19"))
}


#stratify age

{
   
   clindata$age_cat <- clindata$`age:ch1`
   
   clindata$age_cat[clindata$`age:ch1` < 30] = "< 30"
   
   clindata$age_cat[clindata$`age:ch1` >= 30 & clindata$`age:ch1` < 40] ="30s"
   
   clindata$age_cat[clindata$`age:ch1` >= 40 & clindata$`age:ch1`< 50] ="40s"
   
   clindata$age_cat[clindata$`age:ch1` >= 50 & clindata$`age:ch1` < 60] ="50s"
   
   clindata$age_cat[clindata$`age:ch1` >= 60 & clindata$`age:ch1` < 70] ="60s"
   
   clindata$age_cat[clindata$`age:ch1` >= 70] ="70+"
   
   clindata$age_cat[clindata$`age:ch1` == "Unknown"] = NA
   
}


table(clindata$viral_load)
table(clindata$age_cat)


####

{
   MX1 <- ggplot(NULL, aes(x=clindata$positivity, y=log2(t(norm_counts["MX1",]+1)))) +
      geom_jitter(aes(shape=clindata$positivity, color=clindata$positivity), size=3)+
      xlab(NULL) +
      ylab("MX1 expression \n log2 (norm counts +1)") +
      theme(legend.position = "bottom") +
      theme_bw() +
      theme(axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            plot.title =element_text(size = 25),
            legend.position = 'none') +
      stat_summary(fun=mean,
                   geom="point",
                   shape= '_',
                   size=14,
                   colour= c('#b53d35', '#066e70'))
   
   MX1 
}


MX1stat <- wilcox.test(norm_counts["MX1",] ~ clindata$positivity, paired = FALSE)
MX1stat

#####

{
   p_trend_age <- jonckheere.test(x= log2(t(norm_counts["MX1",]+1))[clindata$positivity == "COVID19"],
                                  g= factor(clindata$age_cat[clindata$positivity == "COVID19"], ordered=TRUE),
                                  alternative = "decreasing",
                                  nperm = 500)
   p_trend_age
} 


#####

{
   pairwise_corr <- ggpairs(as.data.frame(log2(t(norm_counts+1))), columns = c("MX1", "MX2", "ACE2", "TMPRSS2"),
             upper = list(continuous = wrap('cor', method = "spearman", size = 3),
                          combo = "box_no_facet", discrete = "count",
                          na ="na"),
             ggplot2::aes(colour=clindata$positivity, shape=clindata$positivity, alpha = 0.01))
   
   pairwise_corr
}


#####

{
  MX1_MX2 <- ggplot(NULL, aes(x = log2(t(norm_counts["MX1",]+1)[which(clindata$positivity=="COVID19" & clindata$ct != "Unknown")]),                                 
                              y = log2(t(norm_counts["MX2",]+1))[which(clindata$positivity=="COVID19"& clindata$ct != "Unknown")],
                              color = as.integer(clindata$ct[(which(clindata$positivity=="COVID19" & clindata$ct != "Unknown"))]))) +
    geom_point(size = 4, na.rm = TRUE) +
    scale_color_gradientn(colours=c("red","white","blue"), name = "Viral
    load (ct)") +
    ylab("MX2 expression RNA-seq \n log2 (norm counts +1)") +
    xlab("MX1 expression RNA-seq \n log2 (norm counts +1)") +
    theme(legend.position = "bottom") +
    theme_bw() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          plot.title =element_text(size = 25))
  
  MX1_MX2 
  
}

#####

MX1_MX2stat <-cor.test(norm_counts["MX1",]  
                       [which(clindata$positivity=="COVID19")],
                       norm_counts["MX2",]
                       [which(clindata$positivity=="COVID19")],
                       method = "spearman")

MX1_MX2stat

#####

res.pca <- prcomp(t(log2(norm_counts[c("MX1","MX2","ACE2","TMPRSS2","HMOX1"),]+1)), 
                  scale = TRUE)



#####
{
   p<- fviz_pca_biplot(res.pca, col.ind = clindata$positivity,
                    geom = "point",
                    addEllipses = TRUE, 
                    palette = c('#F8766D', '#00BFC4'), 
                    title='Principal Component Analysis')
   p
}


#####

table_3Dplot <- cbind(t(log2(norm_counts[c("MX1","ACE2","BSG"),]+1)), as.data.frame(clindata$positivity))
colnames(table_3Dplot)[4] <- "positivity"

canvasXpress(
  data=t(log2(norm_counts[c("MX1","ACE2","BSG"),]+1)),
  varAnnot=as.data.frame(clindata$positivity, row.names=rownames(clindata)),
  axisTickScaleFontFactor=0.6,
  axisTitleScaleFontFactor=0.6,
  ellipseBy="clindata$positivity",
  colorBy="clindata$positivity",
  colorKey=list("clindata$positivity"=list("pos"="#F8766D", "neg"="#00BFC4")),
  graphType="Scatter3D",
  title="3D scatter plot",
  xAxis=list("ACE2"),
  yAxis=list("BSG"),
  zAxis=list("MX1"),
  showLoessFit = FALSE)




###########  troubleshootuing
common_names= intersect(rownames(clindata), colnames(raw_counts))
clindata = clindata[rownames(clindata) %in% common_names,]
raw_counts = raw_counts[,colnames(raw_counts) %in% common_names]


######## cdownload reads

url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152075&format=file&file=GSE152075%5Fraw%5Fcounts%5FGEO%2Etxt%2Egz"
download.file(url, "raw_reads.gz")

raw_counts <- read.delim("raw_reads.gz", stringsAsFactors=FALSE, sep = " ")



#see how many genes per sample have 0 reads

zero <- colSums(norm_counts == 0)/nrow(norm_counts)

hist(zero, breaks = 50)

abline(v = 0.7, col="red")

summary(zero)


#select only those samples with less than 70% of genes with zero

good_samples <- zero <0.7

norm_counts <- norm_counts[,good_samples]

