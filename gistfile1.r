#Sort packages and workspace
setwd("/home/me/sexyplots") 
install.pacakges(c("ggplot2", "stringr", "reshape2", "grid")
library(c("ggplot2", "stringr", "reshape2", "grid")

#Read in the OTU table and improve naming
otu = read.table("otu_table.txt", sep = "\t", row.names=1, header=T, check.names=F, comment.char="", skip = 1)
OTUID = str_c(str_extract(otu$ConsensusLineage, "[^;]+$"),"_",row.names(otu))
rownames(otu) = OTUID
otu = otu[,-c(length(names(otu)))]

#Subselect OTUs and rearrange data
top.otus = as.data.frame(head(sort(rowSums(otu), decreasing = T),20))
otu.20 = otu[rownames(top.otus),]
otu.20.melt = melt(cbind(rownames(otu.20),otu.20))
colnames(otu.20.melt)[1] = "OTU"
colnames(otu.20.melt)[2] = "Sample"
colnames(otu.20.melt)[3] = "Abundance"

#Subselect and percent abundance instead
s_num<-dim(otu)[2]
otu.perc<-matrix(rep("NA",times=(dim(otu)[1]*s_num)),nrow=dim(otu)[1],ncol=s_num)
rownames(otu.perc)<-rownames(otu)
colnames(otu.perc)=colnames(otu) totals<-colSums(otu)

for(s in c(1:s_num)){
vec<-(otu[,s]/totals[s])*100
otu.perc[,s]<-vec
rm(vec)
}

write.table(otu.perc, "otu_perc.txt", col.names = NA, row.names = T, sep = "\t")
otu.perc = read.table("otu_perc.txt", sep = "\t", row.names=1, header=T, check.names=F)

otu.20.perc = otu.perc[rownames(top.otus),]
otu.20.melt.perc = melt(cbind(rownames(otu.20.perc),otu.20.perc))

colnames(otu.20.melt.perc)[1] = "OTU"
colnames(otu.20.melt.perc)[2] = "Sample"
colnames(otu.20.melt.perc)[3] = "Abundance"

#Plot it
pdf("my_plot.pdf" height = 8, width = 8)
otu.plot.data.perc = ggplot(otu.20.melt.perc, aes(x=Sample, y = Abundance, fill = OTU))
otu.barplot.perc = otu.plot.data.perc + geom_bar(stat = "identity", position = "stack")
otu.barplot.perc + facet_grid(OTU ~.) + ggtitle("My Plot") + theme(legend.position = "null", strip.text.y = element_text(size = 10, angle = 0), axis.text.x = element_text(size = 8 , angle = 90), axis.text.y = element_text(size = 6), panel.grid.major = element_line(size = 0.1)) + scale_y_continuous(breaks = c(0,50))

dev.off()