######## Script for Marasco et al. 2021 unders revisoin in Spectrum microbiology    ##########

#Find the paper here: xxxxxxx - this information will be updated upon the final acceptance of the paper

# Packages needed to run this code
library(ggplot2)
library(PerformanceAnalytics)
library(ggtern)
library(igraph)
library(NetSwan)
library(SpiecEasi)
library(RCurl)
library(vegan)
library(gclus)
library(RColorBrewer)
library(ade4)
library(adespatial)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(osfr)
library(vegan)



# download the data for the analysis: see paper for data availability

#Ternary plot analysis
#Adult
tax_tern_adult<-read.table("BD_Adult_Ready.txt",sep='\t',header = T,row.names = 1)

summary(tax_tern_adult)

ad<-ggtern(data=tax_tern_adult, aes(y=replacement, z=similairty, x=Δ.richness,colour=Diet)) + #define data sources
  geom_point(size=3) +          #define data geometry
  theme_bw()+
  scale_colour_manual(values = c("black", "orange", "forestgreen"))+
  ggtitle("Adult") +      #add title
  xlab("Δ.richness") +                       #replace default axis labels
  ylab("replacement") +
  zlab("similairty")
ad

#larvae
tax_tern_larvae<-read.table("BD_larvae_Ready.txt",sep='\t',header = T,row.names = 1)

summary(tax_tern_larvae)

la<-ggtern(data=tax_tern_larvae, aes(y=replacement, z=similairty, x=Δ.richness,colour=Diet)) + #define data sources
  geom_point(size=3) +          #define data geometry
  theme_bw()+
  scale_colour_manual(values = c("black", "orange", "forestgreen"))+
  ggtitle("Larvae") +      #add title
  xlab("Δ.richness") +                       #replace default axis labels
  ylab("replacement") +
  zlab("similairty")
la

#pupae
tax_tern_pupae<-read.table("BD_pupae_Ready.txt",sep='\t',header = T,row.names = 1)

summary(tax_tern_pupae)

pu<-ggtern(data=tax_tern_pupae, aes(y=replacement, z=similairty, x=Δ.richness,colour=Diet)) + #define data sources
  geom_point(size=3) +          #define data geometry
  theme_bw()+
  scale_colour_manual(values = c("black", "orange", "forestgreen"))+
  ggtitle("Pupae") +      #add title
  xlab("Δ.richness") +                       #replace default axis labels
  ylab("replacement") +
  zlab("similairty")
pu


###Netwrok Robustness Analysis
hp5edges <- read.table("{File edge choose}", sep=",", header=TRUE)
head(hp5edges)

hp5edgesmat <- as.matrix(hp5edges[,c(1,2)])

hp5net <- graph_from_edgelist(hp5edgesmat, directed=FALSE)
plot(hp5net, vertex.shape="none", vertex.label.cex=0.6, edge.arrow.size=0.4, layout=layout.kamada.kawai)

hp5netsimple <- simplify(hp5net)

f4<-swan_combinatory(hp5netsimple,10)
plot(f4[,1],f4[,5], type='o', col='yellow',xlab="Fraction of nodes removed",
     ylab="Connectivity loss")
lines(f4[,1],f4[,3], type='o', col='red')
lines(f4[,1],f4[,4], type='o', col='orange')
lines(f4[,1],f4[,2], type='o', col='blue')
legend('bottomright',c("Random", "Betweenness", "Degree", "Cascading"),
       lty=c(1,1,1,1), pch=c(1,1,1,1),
       col=c("yellow","blue","red", "orange"))

fruit<-data.frame(f4[,c(1,2,3,4)])
colnames(fruit)<-c("Nodes_removed","Betweenness","Degree","Cascading")
write.table(fruit,"~/Dropbox/Articoli/HERMETIA/202105_npj Biofilm and Microbiome/01_New network/pupa_sd_robustness.txt",sep='\t',row.names = FALSE)


###Variance partitioning

#####Import OTU TABLE (COUNT DATA)
otu_hermetia<-read.table("otu_transposed.txt",header=T,sep='\t')
attach(otu_hermetia)
head(otu_hermetia)

###IMPORT MAP FILE (be aware that sample id match!!)
stage_hermetia<-read.table("map_file.txt",header=T)

##Filter per Stage if needed subselect the stage to investigate the diet within the stage:  select the Larvae/Pupae/Adult - inthis example we filer per adults.

otu_hermetia_stage <- filter(otu_hermetia, stage_hermetia$Stage == "Adult")
stage_hermetia_stage<- filter(stage_hermetia, stage_hermetia$Stage == "Adult")

##check same number of rows
nrow(otu_hermetia)
nrow(stage_hermetia)


###CHECK NUMERIC
is.numeric(otu_hermetia)
###Transform in numeric
otu_hermetia_stage <-mapply(otu_hermetia_stage, FUN=as.numeric)
otu_hermetia <-mapply(otu_hermetia, FUN=as.numeric)

###Distances
otu_hermetia_dist<-vegdist(log1p(otu_hermetia), method='bray') #Bray-Curtis
otu_hermetia_dist_UPGMA <- hclust(otu_hermetia_dist, method="average")
otu_hermetia_dist_o <- reorder.hclust(otu_hermetia_dist_UPGMA, otu_hermetia_dist)


dend <- as.dendrogram(otu_hermetia_dist_o)
com.tx <- apply(t(otu_hermetia),1,sum) #common taxa for heat map plot
spe.com <- otu_hermetia[,which(com.tx > 20)]
or <- vegemite(spe.com, otu_hermetia_dist_o,scale='log')
heatmap(t(spe.com), Rowv=NA, Colv=dend, cexRow=0.5,cexCol=0.5, col=c('white', brewer.pal(5,"Greens")), scale="column", revC = T, margin=c(4.5,8))


summary(stage_hermetia)
div <- stage_hermetia[,3:4]
div$abund <- apply(otu_hermetia,MARGIN=1,FUN=sum)
div$richness <- specnumber(otu_hermetia)
div$r.richness <- rarefy(otu_hermetia,min(div$abund))
div$shannon <- diversity(otu_hermetia, index='shannon')

plot(specaccum(otu_hermetia), xlab = "streams sampled", ylab = "taxa sampled", ci.type='polygon', ci.col='light gray', ci.lty=0, lwd=2)
plot(specaccum(otu_hermetia[which(div$Stage  == 'Larvae'),]), add=TRUE, ci=0, lwd=2, col='red')
plot(specaccum(otu_hermetia[which(div$Stage == 'Adult'),]), add=TRUE, ci=0, lwd=2, col='blue')
plot(specaccum(otu_hermetia[which(div$Stage == 'Pupae'),]), add=TRUE, ci=0, lwd=2, col='orange')
points(x=c(34,11,23),y=c(57,44,48),pch=16,col=c(1,2,4),cex=1.2)
legend('bottom', legend = c('Larvae','Adult','Pupae','All'), bty='n', pch=16, col=c('red','blue','orange','black'))

par(mfrow=c(1,1))
#beta diversity and dispersion from centroid
beta.BC <- betadisper(otu_hermetia_dist, stage_hermetia$Stage, type = 'centroid')
plot(beta.BC) #ugly

mybetadivplot(beta.BC)

boxplot(beta.BC$distances ~ stage_hermetia$Stage, ylab="Distance to centroid")
plot(beta.BC)
permutest(beta.BC)

#Betadiversity component using beta.div.comp
Larvae.BC <- beta.div.comp(otu_hermetia[stage_hermetia$Stage == 'Larvae',], coef='J', quant=F)
Pupae.BC <- beta.div.comp(otu_hermetia[stage_hermetia$Stage == 'Pupae',], coef='J', quant=F)
Adult.BC <- beta.div.comp(otu_hermetia[stage_hermetia$Stage == 'Adult',], coef='J', quant=F)

BD.mat <- cbind(Larvae.BC$part[2:3],Pupae.BC$part[2:3],Adult.BC$part[2:3],Larvae.BC$part[4:5],Pupae.BC$part[4:5],Adult.BC$part[4:5])
barplot(BD.mat, names.arg=rep(c('Larvae','Pupae','Adult'),2), col=c(4,2), width=0.5, ylab='contribution to total dissimilarity')
axis(1, lwd=0, lwd.ticks=0, at=c(0.66,1.85), labels = c('absolute','relative'), line = 2, cex = 1.5)
legend('topleft', inset=c(0.1,0), legend=c(expression(paste(Delta,' richness',sep='')),'replacement'), pch=22, pt.bg=c(2,4), bty='n', y.intersp = 2)

par(mfrow=c(1,1))


triang1 = cbind(Larvae.BC$rich,(1-Larvae.BC$D), Larvae.BC$repl)
colnames(triang1)<-c("Δ richness","similairty","replacement")
triangle.plot(as.data.frame(triang1), show=F, addmean=T, scale=F,labeltriangle=F)
triang2 = cbind(Pupae.BC$rich,(1-Pupae.BC$D), Pupae.BC$repl)
colnames(triang2)<-c("Δ richness","similairty","replacement")
triangle.plot(as.data.frame(triang2), show=F, labeltriangle=F, addmean = T, scale=F)
triang3 = cbind(Adult.BC$rich,(1-Adult.BC$D), Adult.BC$repl)
colnames(triang3)<-c("Δ richness","similairty","replacement")
triangle.plot(as.data.frame(triang3), show=F, labeltriangle=F, addmean = T, scale=F)

write.table(triang1,"BD_Larvae.txt",sep='\t')
write.table(triang2,"BD_Pupae.txt",sep='\t')
write.table(triang3,"BD_Adult.txt",sep='\t')


#Variation partitioning analysis using Varpart following the previos chink of code
varia<-varpart(otu_hermetia, ~ Stage + Diet, ~ Stage * Diet,data=stage_hermetia)
showvarparts(3)

rda.Stage <- rda(otu_hermetia_stage ~ Stage, data=stage_hermetia_stage)
rda.Diet <- rda (otu_hermetia_stage ~ Diet, data=stage_hermetia_stage)
rda.interaction <- rda (otu_hermetia_stage ~ Stage*Diet, data=stage_hermetia_stage)

RsquareAdj(rda.Stage)
RsquareAdj(rda.Diet)
RsquareAdj(rda.interaction)


varp<-varpart (bc, ~ PlantSpecies*Fraction, ~ Fraction,data=map)
plot (varp, digits = 2)

plot(varia)


###Correlazion analysis among diet growth and taxonomy.

correlo2<-read.table("Larvae_W vs Taxa.txt",header=T,sep="\t")
chart.Correlation(correlo2, histogram = TRUE, method = c("pearson", "kendall", "spearman"))

correlo2<-read.table("Prepupae_W vs Taxa.txt",header=TRUE,sep="\t")
chart.Correlation(correlo2, histogram = TRUE, method = c("pearson", "kendall", "spearman"))

correlo2<-read.table("Adult_W vs Taxa.txt",header=TRUE,sep="\t")
chart.Correlation(correlo2, histogram = TRUE, method = c("pearson", "kendall", "spearman"))

correlo2<-read.table("Survival_larvae.txt",header=TRUE,row.names = 1,sep="\t")
chart.Correlation(correlo2, histogram = TRUE, method = c("pearson", "kendall", "spearman"))

correlo2<-read.table("prepupae_replicates.txt",header=TRUE,row.names = 1,sep="\t")
chart.Correlation(correlo2, histogram = TRUE, method = c("pearson", "kendall", "spearman"))

correlo2<-read.table("adult_replicates.txt",header=TRUE,row.names = 1,sep="\t")
chart.Correlation(correlo2, histogram = TRUE, method = c("pearson", "kendall", "spearman"))

correlo2<-read.table("Diet.txt",header=TRUE,row.names = 1,sep="\t")
chart.Correlation(correlo2, histogram = TRUE, method = c("pearson", "kendall", "spearman"))




##Centrality degree analysis  for co-occurrence netwrok analysis for each diet

centrality_all<-read.table("2019_Network_Paramter_All.txt",header=T,sep="\t")


names(centrality_all)

###Enter the centarlity measure to analyse replacing the response variable in the graph and in the test below
sp_degree<-ggplot(data=centrality_all, aes(x=Diet, y=degree)) +
  geom_boxplot() +
  theme_classic() +
  labs(y="xxx",x="Fractions")  +  theme(legend.position="none") + theme(plot.title = element_text(size=10))
sp_degree
sp_degree_model<-aov(centrality$degree~centrality$Fraction*centrality$kingdom)
par(mfrow=c(2,2))
plot(sp_degree_model)
summary(sp_degree_model)

