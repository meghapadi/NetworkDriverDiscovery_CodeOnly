## Megha Padi
## July 9, 2015
## Code to accompany manuscript 
## "Integrating transcriptional and protein interaction networks to discover condition-specific master regulators"

# -------------- Read in yeast rapamycin data ---------------

## Downloaded from http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-412/ on May 5, 2015
## expression matrix is zipped in a package called E-MTAB-412.processed.1.zip

rap.time <- read.table("RapTC_20100113_noAnnotation.txt",header=T,sep="\t")

rap.time.ids <- as.character(rap.time[2:nrow(rap.time),1])
rap.time.ids <- sapply(strsplit(rap.time.ids,split=":"),function(x){x[4]})

rap.time <- rap.time[2:nrow(rap.time),2:ncol(rap.time)]
rownames(rap.time) <- rap.time.ids

## Download from http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-412/files/ on May 5, 2015
## array annotation file 

rap.time.annot <- read.table("A-AFFY-47.adf.txt",skip=18,header=T,sep="\t",quote="\"")

rap.annot.ids <- as.character(rap.time.annot[214:5929,1])
rap.annot.orf <- as.character(rap.time.annot[214:5929,2])
names(rap.annot.orf) <- rap.annot.ids

## all rap.time.ids are present within rap.annot.ids

rap.time.orf <- rap.annot.orf[rap.time.ids]
rap.time.orf <- sapply(strsplit(rap.time.orf,split=".S1"),function(x){x[1]})
#rownames(rap.time) <- rap.time.orf

rap.time1 <- t(apply(rap.time,1,as.numeric))

## some probesets map to the same ORF, so here we create a single row corresponding to each ORF by taking the median expression level across probesets

rap.time.unique <- unique(rap.time.orf)

rtr.dat <- NULL
for (i in 1:length(rap.time.unique))
{
	if (i %% 100 ==0) print(i)
	this.id <- rap.time.ids[rap.time.orf %in% rap.time.unique[i]]
	if (length(this.id)>1)
		rtr.dat <- rbind(rtr.dat,apply(rap.time1[this.id,],2,median)) else
		rtr.dat <- rbind(rtr.dat,rap.time1[this.id,])
}

colnames(rtr.dat) <- colnames(rap.time)
rownames(rtr.dat) <- rap.time.unique

write.table(rtr.dat,"RapamycinGeneExpression_collapsed.txt",col.names=T,row.names=T,sep="\t", quote=F)

write.table(rownames(rtr.dat),"YeastInput/Rapamycin_gene_universe.txt",col.names=F,row.names=F,quote=F,sep="\t")

rtr.dat <- read.table("RapamycinGeneExpression_collapsed.txt")

# ------------- Read in and process phenotype data --------------

## Drug phenotypes downloaded from "Genetic basis of individual differences in the response to small-molecule drugs in yeast" Nature Genetics, 2007 
## Supplementary Table 2, called ng1991-S3.xls

drug.data <- read.table("DrugResponseData/drug_response.txt")
drug.rows <- t(read.table("DrugResponseData/drug_response_rows.txt",sep="\t"))[1,]
drug.cols <- t(read.table("DrugResponseData/drug_response_cols.txt",sep="\t"))[1,]

rap.rows <- grep("rapamycin",drug.rows)

dat.rap <- drug.data[rap.rows,]
dat.rap.mean <- apply(dat.rap,2,mean)
names(dat.rap.mean) <- drug.cols

## get the strains corresponding to the time course data
## Download from http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-412/files/ on May 5, 2015

rap.time.samples <- read.table("E-MTAB-412.sdrf.txt",header=T,sep="\t")
rap.samp.cel <- as.character(rap.time.samples[,24])
rap.samp.strain <- as.character(rap.time.samples[,1])
names(rap.samp.strain) <- rap.samp.cel

rap.time.cel <- readLines("RapTC_20100113_noAnnotation.txt")[1]
rap.time.cel <- strsplit(rap.time.cel,split="\t")[[1]][2:583]

colnames(rtr.dat) <- rap.time.cel

rap.time.strain <- rap.samp.strain[rap.time.cel]
rap.unique.strains <- unique(rap.time.strain)

## make a vector of phenotypes corresponding to all the strains in the gene expression data

rtr.pheno <- NULL
for (i in 1:length(rap.time.strain))
{
	if (rap.time.strain[i] %in% drug.cols)
		rtr.pheno <- c(rtr.pheno,dat.rap.mean[rap.time.strain[i]]) else
		rtr.pheno <- c(rtr.pheno,NA) 
}

write.table(rtr.pheno,"RapamycinGeneExp_Pheno.txt",col.names=F,row.names=F,quote=F,sep="\t")

# -------------- Subset the data for the 30 most sensitive strains -----------------

pheno.dat <- read.table("RapamycinGeneExp_Pheno.txt")
pheno.notna <- pheno.dat[!(is.na(pheno.dat[,1])),]
dat.notna <- rtr.dat[,!(is.na(pheno.dat[,1]))]

rtr.sens <- dat.notna[,pheno.notna < -0.1]

# -------------- Read in rapamycin driver genes ----------------

## Downloaded Supplementary Table 2 from Xie et al (PNAS 2005) which is named 00297Table2.xls
## took all 396 ORF names (associated with any phenotype change under rapamycin, positive or negative)

rap.drivers <- as.character(read.table("RapamycinDrivers.txt")[,1])

# -------------- Reading in yeast annotation -----------------

## File was downloaded from URL http://downloads.yeastgenome.org/curation/chromosomal_feature/dbxref.tab  at 10:30am on May 5, 2015

dbxref <- read.table("SGD_May_5_2015_dbxref.tab.txt",sep="\t",fill=T,quote="\"") 
dbxref.sym <- dbxref[dbxref[,6]!="",c(4,6)]
dbxref.tab <- unique(apply(dbxref.sym,1,function(x){paste(x,collapse="\t")}))

dbxref.final <- NULL
for (i in 1:length(dbxref.tab))
dbxref.final <- rbind(dbxref.final,strsplit(dbxref.tab[i],split="\t")[[1]])

rownames(dbxref.final) <- dbxref.final[,2]

## Note that this final set of ORFs is not complete - this only represents ORFs which have gene symbols associated with them. This table should be used only to convert stray symbols into systematic ORF names, not as a complete universe of ORFs


# -------------- Reading in TF-target prior -----------------

## Downloaded from http://fraenkel.mit.edu/improved_map/orfs_by_factor.tar.gz on May 5, 2015
## file called "orfs_by_factor_p0.005_cons0.txt"

mac.tab <- readLines("Fraenkel_orfs_by_factor_p0.005_cons0.txt")

mac.list <- NULL
mac.tfs <- NULL
for (i in 1:length(mac.tab))
{
	this.split <- strsplit(mac.tab[i],split="\t")[[1]]
	mac.tfs <- c(mac.tfs,this.split[1])
	this.targ <- this.split[2:length(this.split)]
	targ.sym <- which(this.targ %in% dbxref.final[,2])
	if (length(targ.sym)>0)
	{
		this.targ[targ.sym] <- dbxref.final[this.targ[targ.sym],1]
	}
	mac.list <- append(mac.list,list(this.targ))
}

mac.tfs.orf <- NULL
for (i in 1:length(mac.tfs))
{
	if (mac.tfs[i] %in% dbxref.final[,2]) 
	mac.tfs.orf <- c(mac.tfs.orf,dbxref.final[mac.tfs[i],1]) else mac.tfs.orf <- c(mac.tfs.orf,mac.tfs[i])
}
	
names(mac.list) <- mac.tfs.orf

# --------------- Differential expression of rapamycin data and clustering ---------------
## Figure 2A

library(limma)

dat.limma <- rtr.sens

time.pts <- c("t0m","t10m","t20m","t30m","t40m","t50m")
lev <- time.pts
f <- factor(rep(time.pts,30),levels=lev)

design <- model.matrix(~0+f)
colnames(design) <- lev
fit <- lmFit(dat.limma,design)

cont.diff <- makeContrasts(t50m-t0m, levels=design)
fit2 <- contrasts.fit(fit,cont.diff)
fit2 <- eBayes(fit2)

rtr.sens.diff <- topTable(fit2,adjust="BH",number=nrow(dat.limma))
rownames(rtr.sens.diff) -> rtr.allgenes

write.table(rtr.allgenes,"YeastInput/Rap_time_allgenes.txt",col.names=F,row.names=F,quote=F,sep="\t")

rtr.allgenes <- t(read.table("YeastInput/Rap_time_allgenes.txt"))[1,]

## Figure 2A

library(mclust)

top20perc <- rtr.allgenes[1:1131]
rtr.sens.top <- rtr.sens[top20perc,]

rst.ctr <- rtr.sens.top-apply(rtr.sens.top,1,mean)

rst.ctr.sd <- NULL
for (i in 1:nrow(rst.ctr))
	rst.ctr.sd <- rbind(rst.ctr.sd,rst.ctr[i,]/sd(rst.ctr[i,]))

library(mclust)
rst.clust <- Mclust(rst.ctr.sd,G=2:30)

write.table(cbind(rownames(rtr.sens.top),rst.clust$classification),"YeastOutput/Top20percent_Mclust.txt",col.names=F,row.names=F,quote=F,sep="\t")

rst.ctr.sd1 <- NULL
for (i in 1:6)
rst.ctr.sd1 <- cbind(rst.ctr.sd1,as.matrix(rst.ctr.sd[,c(1:30)*6-6+i]))
rownames(rst.ctr.sd1) <- rownames(rtr.sens.top)

mclust.memb <- rst.clust$classification
names(mclust.memb) <- rownames(rtr.sens.top)

hm.dat <- NULL
for (i in 1:max(mclust.memb))
{
	this.clust <- names(mclust.memb)[mclust.memb==i]
	this.dat <- apply(rst.ctr.sd1[this.clust,],2,mean)
	hm.dat <- rbind(hm.dat,this.dat)
}

library(gplots)

heatmap.2(hm.dat,density.info=NULL,scale="none",symbreaks=TRUE,col=bluered(75),dendrogram="row", Rowv=T,Colv=F, trace="none", cexCol = 0.5, labCol=F, cexRow=0.5,labRow=c(1:25))

## Functional enrichment

library(GOstats)
library(org.Sc.sgd.db)

this.comm <- top20perc[mclust.memb==1]

params<- new("GOHyperGParams",geneIds = this.comm,
	universeGeneIds = rtr.allgenes,
	ontology= "BP",
	conditional = TRUE,
	testDirection="over",
	annotation="org.Sc.sgd.db")
	
hm.tab.true <- NULL
for (i in 1:max(mclust.memb))
{
	print(i)
	this.comm <- top20perc[mclust.memb==i]
	if (length(this.comm)>5)
	{
		geneIds(params) <- this.comm
		conditional(params) <- TRUE
		go.res <- hyperGTest(params)
		this.df <- summary(go.res)
		this.df.p <- p.adjust(this.df[,2],method="BH")
		if (sum(this.df.p<0.05)>0)
		{
			this.tab <- this.df[this.df.p<0.05,c(1:3,7)]
			hm.tab.true <- rbind(hm.tab.true,cbind(rep(i,nrow(this.tab)),this.tab))
		}
	}
	write.table(hm.tab.true,"YeastOutput/Fig2A_GOterms.txt",col.names=T,row.names=F,quote=F,sep="\t")
}

## Enrichment in TF binding sites
## Figure 2C

mclust.memb1 <- read.table("YeastOutput/Top20percent_Mclust.txt")
rap.univ <- t(read.table("YeastInput/Rapamycin_gene_universe.txt"))[1,]

tf.enrich <- NULL
for (i in 1:max(mclust.memb1[,2]))
{
	this.clus <- as.character(mclust.memb1[mclust.memb1[,2]==i,1])
	clus.pvals <- NULL
	for (j in 1:length(mac.list))
	{
		this.targ <- mac.list[[j]]
		
		a <- length(intersect(this.targ,this.clus))
		b <- length(intersect(rap.univ,this.targ))-a
		c <- length(this.clus) -a
		d <- length(rap.univ)-a-b-c
		
		this.pval <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
		clus.pvals <- c(clus.pvals,this.pval)
	}
	clus.pvals <- p.adjust(clus.pvals,method="BH")
	if (sum(clus.pvals<0.05)>0)
	tf.enrich[[i]] <- mac.tfs.orf[clus.pvals<0.05]	else tf.enrich <- append(tf.enrich,list(NULL))
}

tf.enr.tab <- sapply(tf.enrich,function(x){paste(x,collapse=" ")})
write.table(tf.enr.tab,"YeastOutput/Rapamycin_enrichedTFs.txt",col.names=F,row.names=T,quote=F,sep="\t")

all.tf.enr <- unique(unlist(tf.enrich))

a <- length(intersect(all.tf.enr,rap.drivers))
b <- length(all.tf.enr)-a
c <- length(intersect(mac.tfs.orf,rap.drivers))-a
d <- length(mac.tfs.orf)-a-b-c

fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")

venn.diagram(list(B=all.tf.enr,A=intersect(mac.tfs.orf,rap.drivers)),fill=c("purple","red"),alpha=c(0.5,0.5),cex=0,filename = "YeastOutput/Fig2C_Venn_EnrichTFs_Drivers.tiff",ext.text=F,cat.cex=0,rotation.degree=180)


# ---------- Enrichment of drivers in differentially expressed genes -----------
##  Figure 2B

rtr.topgenes <- rownames(rtr.sens.diff)[1:565]
rtr.allgenes <- rownames(rtr.sens.diff)

a <- length(intersect(rap.drivers,rtr.topgenes))
b <- length(intersect(rap.drivers,rtr.allgenes)) - a
c <- length(rtr.topgenes) - a
d <- length(rtr.allgenes) - a-b-c

fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")
fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="less")

library(VennDiagram)

venn.diagram(list(B=intersect(rap.drivers,rtr.allgenes),A=rtr.topgenes),fill=c("red","forestgreen"),alpha=c(0.5,0.5),cex=0,filename = "YeastOutput/Fig2B_Venn_DiffExp_Drivers.tiff",ext.text=F,cat.cex=0,rotation.degree=180)


## plot the enrichment of driver genes in the top "n" DEGs 

driv.pval <- NULL
for (i in 1:length(rtr.allgenes))
{
	if (i %% 500 ==0) print(i)
	a <- length(intersect(rap.drivers,rtr.allgenes[1:i]))
	b <- length(intersect(rap.drivers,rtr.allgenes)) - a
	c <- i - a
	d <- length(rtr.allgenes) - a-b-c
	this.pval <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
	driv.pval <- c(driv.pval,this.pval)	
}


plot(-1*log2(driv.pval),col="white",ylim=c(0,5),xlab="Number of DEGs",ylab="-log2(P)")
lines(-1*log2(driv.pval),col="black",lwd=2)
lines(c(0,length(rtr.allgenes)+300),c(-1*log2(0.05),-1*log2(0.05)),col="red",lty=2,lwd=2)

## ROC curve for all DEGs

tpr.vec <- fpr.vec <- NULL 
this.ord <- rtr.allgenes
for (j in 1:length(this.ord))
{
	tp <- length(intersect(this.ord[1:j],rap.drivers))
	fp <- j-tp
	fn <- length(intersect(rap.drivers,rtr.allgenes))-tp
	tn <- length(rtr.allgenes) -tp-fp-fn
		
	tpr.vec <- c(tpr.vec,tp/(tp+fn))
	fpr.vec <- c(fpr.vec,fp/(fp+tn))
}


plot(fpr.vec,tpr.vec,col="white",xlab="False positive rate",ylab="True positive rate")
lines(fpr.vec,tpr.vec,col="forestgreen",lwd=2)
lines(c(0,1),c(0,1),lty=2,col="black",lwd=2)

## enrichment statistics for ranking by differential expression

diff.scores <- rtr.sens.diff[,6]
names(diff.scores) <- rtr.allgenes

diff.driv <- intersect(rtr.allgenes,rap.drivers)
diff.nondriv <- rtr.allgenes[!(rtr.allgenes %in% rap.drivers)]

ks.test(diff.scores[diff.driv],diff.scores[diff.nondriv],exact=FALSE,alternative="less")

wilcox.test(diff.scores[diff.driv],diff.scores[diff.nondriv],exact=FALSE,alternative="greater")


# ------------ DREM results -------------

## Outputing data in the format for DREM

eachcols <- c("0m","10m","20m","30m","40m","50m")

for (i in 1:30)
{
	this.dat <- rtr.sens[,(6*i-5):(6*i)]
	this.dat <- cbind(rownames(this.dat),this.dat)
	colnames(this.dat) <- c("Gene Symbol",eachcols)
	write.table(this.dat,paste(c("YeastInput/DREM/RapSens_",i,".txt"),collapse=""),col.names=T,row.names=F,sep="\t",quote=F)
}

## DREM results

drem.out <- c("GCN4","GLN3","SKN7","MSN2","RTG3","FHL1","SWI6","RAP1")
drem.out <- unique(c("GCN4","MSN2","YAP7","GAT1","GLN3","MSN4","SKN7","RTG3","SKO1","GLN3","GAT1","GCN4","HSF1","SFP1","FHL1","RAP1","MBP1","ROX1","GCN4"))

drem.orf <- dbxref.final[dbxref.final[,2] %in% drem.out,1]

a <- length(intersect(drem.orf,rap.drivers))
b <- length(drem.orf)-a
c <- length(intersect(mac.tfs.orf,rap.drivers))-a
d <- length(mac.tfs.orf)-a-b-c

fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")

venn.diagram(list(B=drem.orf,A=intersect(mac.tfs.orf,rap.drivers)),fill=c("yellow","red"),alpha=c(0.5,0.5),cex=0,filename = "YeastOutput/Fig2C_Venn_DREM_Drivers.tiff",ext.text=F,cat.cex=0,rotation.degree=180)


# ------------ Infer yeast rapamycin transcriptional networks ------------

source("NetworkAnalysisFunctions.R")

## use GMIT to infer network

yeast.prior.genes <- unique(c(mac.tfs.orf,unlist(mac.list)))
rap.net.genes <- intersect(yeast.prior.genes,rtr.allgenes)

rap.disc <- DatDiscOut(rap.net.genes,rtr.sens)

write.table(rap.disc,"YeastInput/Rapamycin_disc.txt",col.names=F,row.names=F,quote=F,sep="\t")

mac.links <- NULL
for (i in 1:length(mac.list))
{
	mac.links <- rbind(mac.links,cbind(rep(mac.tfs.orf[i],length(mac.list[[i]])),mac.list[[i]]))
}

rap.prior <- PriorMatrix(rap.net.genes,rtr.sens,mac.links)

write.table(rap.prior,"YeastInput/rapamycin_prior.txt",col.names=F,row.names=F,quote=F,sep="\t")

rap.gmit.net <- read.table("YeastOutput/Rapamycin_GMITnet.txt",skip=1)

rap.gmit.edges <- NetworkReadGMIT(rap.net.genes,rap.gmit.net,mac.tfs.orf)

write.table(rap.gmit.edges,"YeastOutput/Rapamycin_GMIT_edges.txt",col.names=F,row.names=F,quote=F,sep="\t")

## use PANDA to infer network

rtr.dat <- read.table("RapamycinGeneExpression_collapsed.txt")

pheno.dat <- t(read.table("RapamycinGeneExp_Pheno.txt"))[1,]

pheno.notna <- pheno.dat[!(is.na(pheno.dat))]
dat.notna <- rtr.dat[,!(is.na(pheno.dat))]

rtr.sens <- dat.notna[,pheno.notna < -0.1]

rtr.sens.pan <- cbind(rownames(rtr.sens),rtr.sens)

write.table(rtr.sens.pan,"YeastInput/Rapamycin_dat_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")

rap.pan.ran <- NULL
for (i in 1:nrow(rtr.sens.pan))
{
	if (i %% 500 ==0) print(i)
	this.ran <- sample(rtr.sens.pan[i,],ncol(rtr.sens.pan),replace=F)
	rap.pan.ran <- rbind(rap.pan.ran,this.ran)
}
rap.pan.ran <- rap.pan.ran[sample(1:nrow(rap.pan.ran),nrow(rap.pan.ran),replace=F),]

rap.pan.ran <- cbind(rownames(rtr.sens),rap.pan.ran)

write.table(rap.pan.ran,"YeastInput/Rapamycin_dat_PANDA_random.txt",col.names=F,row.names=F,quote=F,sep="\t")
	

rtr.genes <- rownames(rtr.sens)

rap.pan.prior <- NULL
for (i in 1:length(mac.list))
{
	if (mac.tfs.orf[i] %in% rtr.genes)
	{
		this.targ <- intersect(rtr.genes,mac.list[[i]])
		rap.pan.prior <- rbind(rap.pan.prior,cbind(rep(mac.tfs.orf[i],length(this.targ)),this.targ))
	}
}

rap.pan.prior <- cbind(rap.pan.prior,rep(1,nrow(rap.pan.prior)))

write.table(rap.pan.prior,"YeastInput/rapamycin_motif_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")


## run Panda

## read in result and get degree in network

rap.pan <- read.table("~/Documents/Yeast/Manuscript/DriverCentrality_code/YeastOutput/Rapamycin_PANDA_Output_FinalNetwork.pairs",header=F,sep="\t")
rap.pan.ran <- read.table("~/Documents/Yeast/Manuscript/DriverCentrality_code/YeastOutput/Rapamycin_PANDA_random_Output_FinalNetwork.pairs",header=F,sep="\t")

rap.pan.edges <- NetworkReadPanda(rap.pan,rap.pan.ran)



# ----------------- Analyze rapamycin transcriptional network ------------------


## Figure 3A

library(igraph)
NetworkAnalysis(rap.gmit.edges,rap.drivers,"YeastOutput/Rapamycin/Y2Honly/",T,mac.tfs.orf,y2h.deg)

rap.gmit.edges <- read.table("YeastOutput/Rapamycin_GMIT_edges.txt")
rap.graph <- graph.edgelist(as.matrix(rap.gmit.edges),directed=F)
rap.degree <- degree(rap.graph)
rap.degord <- V(rap.graph)$name[order(rap.degree,decreasing=T)]
rap.driv.trxord <- rap.degord[rap.degord %in% rap.drivers]
rap.driv.trxord.deg <- rap.degree[rap.driv.trxord]

orf.to.name <- function(x){
	if (x %in% dbxref.final[,1]) dbxref.final[dbxref.final[,1] %in% x,2] else x
}

rap.driv.trxord.names <- sapply(rap.driv.trxord,orf.to.name)

write.table(cbind(rap.driv.trxord.names[1:10],rap.driv.trxord.deg[1:10]),"YeastOutput/Top10rapamycinDrivers_GMITdeg.txt",col.names=F,row.names=F,quote=F,sep="\t")

rap.gmit.tfs <- intersect(mac.tfs.orf,names(rap.degree))

all.scores <- array(0,dim=c(2,length(mac.tfs.orf)))
colnames(all.scores) <- mac.tfs.orf

all.scores[1,rap.gmit.tfs] <- rap.degree[rap.gmit.tfs]

rap.diff.tfs <- intersect(rtr.allgenes,mac.tfs.orf)
rtr.sens.score <- rtr.sens.diff[,6]
names(rtr.sens.score) <- rownames(rtr.sens.diff)

all.scores[2,rap.diff.tfs] <- rtr.sens.score[rap.diff.tfs]


tpr.tab <- array(1,dim=dim(all.scores))
fpr.tab <- array(1,dim=dim(all.scores))
for (i in 1:nrow(all.scores))
{
	this.score <- all.scores[i,]
	this.ord <- mac.tfs.orf[order(this.score,decreasing=T)]
	for (j in 1:length(this.ord))
	{
		tp <- length(intersect(this.ord[1:j],rap.drivers))
		fp <- j-tp
		fn <- length(intersect(rap.drivers,mac.tfs.orf))-tp
		tn <- length(mac.tfs.orf) -tp-fp-fn
		
		tpr.tab[i,j] <- tp/(tp+fn)
		fpr.tab[i,j] <- fp/(fp+tn)
	}
}
	
plot.colors <- c("violetred3","forestgreen")
	
pdf("YeastOutput/Rapamycin/Fig3A_ROC_degree_diffexp_small.pdf",width=4,height=4.2)
plot(fpr.tab[1,],tpr.tab[1,],col="white",xlab="False positive rate",ylab="True positive rate",cex.lab=1.3)
for (i in 1:nrow(all.scores))
	lines(fpr.tab[i,],tpr.tab[i,],col=plot.colors[i],lwd=3)
lines(c(0,1),c(0,1),lty=2,col="black",lwd=3)
dev.off()

rap.tf.driv <- intersect(mac.tfs.orf,rap.drivers)
rap.tf.nondriv <- mac.tfs.orf[!(mac.tfs.orf %in% rap.drivers)]

ks.pvals <- NULL
	for (i in 1:nrow(all.scores))
	ks.pvals <- c(ks.pvals,ks.test(all.scores[i,rap.tf.driv],all.scores[i,rap.tf.nondriv],exact=FALSE,alternative="less")$p.value)
	
	w.pvals <- NULL
for (i in 1:nrow(all.scores))
w.pvals <- c(w.pvals,wilcox.test(all.scores[i,rap.tf.driv],all.scores[i,rap.tf.nondriv],exact=FALSE,alternative="greater")$p.value)


## effect size of overlap with top 10% TFs

or.percs <- c(10,20,30)

tf.ord.or.list <- NULL
tf.ord.p.list <- NULL
for (k in 1:3)
{
this.num <- floor(or.percs[k]*length(mac.tfs.orf)/100)	
tf.ord.or <- NULL
tf.ord.p <- NULL
for (i in 1:nrow(all.scores))
{
	this.ord <- mac.tfs.orf[order(all.scores[i,],decreasing=T)]
		
	a <- length(intersect(this.ord[1:this.num],rap.tf.driv))
	b <- length(rap.tf.driv)-a
	c <- this.num-a
	d <- length(mac.tfs.orf)-a-b-c
		
	this.p <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
	this.or <- (a*d)/(b*c)
	tf.ord.p <- c(tf.ord.p,this.p)
	tf.ord.or <- c(tf.ord.or,this.or)
}
tf.ord.or.list[[k]] <- tf.ord.or
tf.ord.p.list[[k]] <- tf.ord.p
}

par(mfrow=c(1,3))
for (i in 1:3)
barplot(tf.ord.or.list[[i]],col=plot.colors,cex.axis=2,ylim=c(0,5.5))


## Essential genes downloaded from Stanford yeast deletion project website 
## http://www-sequence.stanford.edu/group/yeast_deletion_project/downloads.html
## on May 8, 2015

essential.file <- read.table("Essential_ORFs.txt",header=T,sep="\t",quote="\"")
essential.file <- essential.file[2:nrow(essential.file),]

ess.orfs <- as.character(essential.file[,2])
ess.orfs.split <- sapply(ess.orfs,function(x){strsplit(x,split=" ")[[1]][1]})
ess.orfs.fin <- toupper(ess.orfs.split)

## Visualize network

rap.gmit.tf <- rap.gmit.edges[rap.gmit.edges[,1] %in% mac.tfs.orf,]
rap.gmit.tf <- rap.gmit.tf[rap.gmit.tf[,2] %in% mac.tfs.orf,]

rap.gmit.tf.name <- cbind(sapply(rap.gmit.tf[,1],orf.to.name),sapply(rap.gmit.tf[,2],orf.to.name))

write.table(rap.gmit.tf.name,"YeastOutput/Cytoscape/Rapamycin_GMIT_TFnet.txt",col.names=T,row.names=F,quote=F,sep="\t")

ess.genes <- sapply(ess.orfs.fin,orf.to.name)

write.table(cbind(ess.genes,rep("essential",length(ess.genes))),"YeastOutput/Cytoscape/EssentialGenes.txt",col.names=T,row.names=F,quote=F,sep="\t")

rap.driv.genes <- sapply(rap.drivers,orf.to.name)

write.table(cbind(rap.driv.genes,rep("RapDrivers",length(rap.driv.genes))),"YeastOutput/Cytoscape/RapamycinDrivers.txt",col.names=T,row.names=F,quote=F,sep="\t")

rap.deg.genes <- sapply(names(rap.degree),orf.to.name)

write.table(cbind(rap.deg.genes,rap.degree),"YeastOutput/Cytoscape/RapamycinDegree.txt",col.names=T,row.names=F,quote=F,sep="\t")

# ------------------- Analyze PANDA network for rapamycin ---------------------

library(igraph)

NetworkAnalysis(as.matrix(rap.pan.edges),rap.drivers,"YeastOutput/RapamycinPANDA/Y2Honly/",TRUE,mac.tfs.orf,y2h.deg)

rap.pan.graph <- graph.edgelist(as.matrix(rap.pan.edges),directed=F)
rap.pan.degree <- degree(rap.pan.graph)

rap.pan.tfs <- intersect(mac.tfs.orf,names(rap.pan.degree))

all.scores <- array(0,dim=c(2,length(mac.tfs.orf)))
colnames(all.scores) <- mac.tfs.orf

all.scores[1,rap.pan.tfs] <- rap.pan.degree[rap.pan.tfs]

rap.diff.tfs <- intersect(rtr.allgenes,mac.tfs.orf)
rtr.sens.score <- rtr.sens.diff[,6]
names(rtr.sens.score) <- rownames(rtr.sens.diff)

all.scores[2,rap.diff.tfs] <- rtr.sens.score[rap.diff.tfs]


tpr.tab <- array(1,dim=dim(all.scores))
fpr.tab <- array(1,dim=dim(all.scores))
for (i in 1:nrow(all.scores))
{
	this.score <- all.scores[i,]
	this.ord <- mac.tfs.orf[order(this.score,decreasing=T)]
	for (j in 1:length(this.ord))
	{
		tp <- length(intersect(this.ord[1:j],rap.drivers))
		fp <- j-tp
		fn <- length(intersect(rap.drivers,mac.tfs.orf))-tp
		tn <- length(mac.tfs.orf) -tp-fp-fn
		
		tpr.tab[i,j] <- tp/(tp+fn)
		fpr.tab[i,j] <- fp/(fp+tn)
	}
}
	
plot.colors <- c("violetred3","forestgreen")
	
pdf("YeastOutput/RapamycinPANDA/Fig3ASupp_ROC_degree_diffexp.pdf",width=4,height=4.2)
plot(fpr.tab[1,],tpr.tab[1,],col="white",xlab="False positive rate",ylab="True positive rate",cex.lab=1.3)
for (i in 1:nrow(all.scores))
	lines(fpr.tab[i,],tpr.tab[i,],col=plot.colors[i],lwd=3)
lines(c(0,1),c(0,1),lty=2,col="black",lwd=3)
dev.off()

rap.tf.driv <- intersect(mac.tfs.orf,rap.drivers)
rap.tf.nondriv <- mac.tfs.orf[!(mac.tfs.orf %in% rap.drivers)]

ks.pvals <- NULL
	for (i in 1:nrow(all.scores))
	ks.pvals <- c(ks.pvals,ks.test(all.scores[i,rap.tf.driv],all.scores[i,rap.tf.nondriv],exact=FALSE,alternative="less")$p.value)

w.pvals <- NULL
for (i in 1:nrow(all.scores))
w.pvals <- c(w.pvals,wilcox.test(all.scores[i,rap.tf.driv],all.scores[i,rap.tf.nondriv],exact=FALSE,alternative="greater")$p.value)


or.percs <- c(10,20,30)

tf.ord.or.list <- NULL
tf.ord.p.list <- NULL
for (k in 1:3)
{
this.num <- floor(or.percs[k]*length(mac.tfs.orf)/100)	
tf.ord.or <- NULL
tf.ord.p <- NULL
for (i in 1:nrow(all.scores))
{
	this.ord <- mac.tfs.orf[order(all.scores[i,],decreasing=T)]
		
	a <- length(intersect(this.ord[1:this.num],rap.tf.driv))
	b <- length(rap.tf.driv)-a
	c <- this.num-a
	d <- length(mac.tfs.orf)-a-b-c
		
	this.p <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
	this.or <- (a*d)/(b*c)
	tf.ord.p <- c(tf.ord.p,this.p)
	tf.ord.or <- c(tf.ord.or,this.or)
}
tf.ord.or.list[[k]] <- tf.ord.or
tf.ord.p.list[[k]] <- tf.ord.p
}

par(mfrow=c(1,3))
for (i in 1:3)
barplot(tf.ord.or.list[[i]],col=plot.colors,cex.axis=2,ylim=c(0,4.5))


## Visualize network

rap.pan.tf <- rap.pan.edges[rap.pan.edges[,1] %in% mac.tfs.orf,]
rap.pan.tf <- rap.pan.tf[rap.pan.tf[,2] %in% mac.tfs.orf,]

rap.pan.tf.name <- cbind(sapply(rap.pan.tf[,1],orf.to.name),sapply(rap.pan.tf[,2],orf.to.name))

write.table(rap.pan.tf.name,"YeastOutput/Cytoscape/Rapamycin_PANDA_TFnet.txt",col.names=T,row.names=F,quote=F,sep="\t")

rap.pan.deg.genes <- sapply(names(rap.pan.degree),orf.to.name)

write.table(cbind(rap.pan.deg.genes,rap.pan.degree),"YeastOutput/Cytoscape/Rapamycin_PANDA_Degree.txt",col.names=T,row.names=F,quote=F,sep="\t")



# ---------------------- Yeast protein-protein interaction ------------------------


## Y2H only

krogan <- read.table("PPIdata/Krogan_2006_TableS7.txt",quote="\"")
krogan.ppi <- krogan[,c(2,4)]

krogan1 <- read.table("PPIdata/Krogan_2006_TableS8.txt",fill=T,sep="\t",quote="\"")
krogan.ppi1 <- krogan1[,c(2,4)]

krogan.ppi2 <- rbind(krogan.ppi,krogan.ppi1)

vidal <- read.table("PPIdata/Vidal_Yeast_Y2H_union.txt",quote="\"")

krogan.ppi2 <- as.matrix(krogan.ppi2)
vidal <- as.matrix(vidal)

y2h.yeast <- rbind(krogan.ppi2,vidal)
yy.self <- apply(y2h.yeast,1,function(x){x[1]==x[2]})
y2h.yeast <- y2h.yeast[!yy.self,]

yy.graph <- graph.edgelist(y2h.yeast,directed=F)
yy.simple <- simplify(yy.graph)
yy.deg <- degree(yy.simple)

write.table(yy.deg,"PPIdata/Y2H_Yeast_degree.txt",col.names=F,row.names=T,quote=F,sep="\t")

y2h.deg <- read.table("PPIdata/Y2H_Yeast_degree.txt")

yy.degord <- V(yy.simple)$name[order(yy.deg,decreasing=T)]
yy.driv <- intersect(V(yy.simple)$name,rap.drivers)
yy.nondriv <- V(yy.simple)$name[!(V(yy.simple)$name %in% rap.drivers)]

ks.test(yy.deg[yy.driv],yy.deg[yy.nondriv],exact=FALSE,alternative="less")$p.value

yy.deg.graph <- induced.subgraph(yy.simple,vids=intersect(rtr.allgenes[1:1000],V(yy.simple)$name))
yy.deg.deg <- degree(yy.deg.graph)
yy.ddord <- V(yy.deg.graph)$name[order(yy.deg.deg,decreasing=T)]


## Y2H + BioGRID (BioGRID datasets for yeast and human were downloaded May 18, 2015)

biogrid <- read.table("PPIdata/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.3.124.tab2.txt",sep="\t",header=T,quote="\"",fill=T)

yb.table <- as.matrix(biogrid[,c(6,7)])

ya.graph <- graph.edgelist(rbind(y2h.yeast,yb.table),directed=F)
ya.simple <- simplify(ya.graph)
ya.deg <- degree(ya.simple)

write.table(ya.deg,"PPIdata/All_Yeast_degree.txt",col.names=F,row.names=T,quote=F,sep="\t")

ya.deg <- read.table("PPIdata/All_Yeast_degree.txt")

ya.degord <- V(ya.simple)$name[order(ya.deg,decreasing=T)]
ya.driv <- intersect(V(ya.simple)$name,rap.drivers)
ya.nondriv <- V(ya.simple)$name[!(V(ya.simple)$name %in% rap.drivers)]

ks.test(ya.deg[ya.driv],ya.deg[ya.nondriv],exact=FALSE,alternative="less")$p.value


ya.deg.graph <- induced.subgraph(ya.simple,vids=intersect(rtr.allgenes[1:2000],V(ya.simple)$name))
ya.deg.deg <- degree(ya.deg.graph)
ya.ddord <- V(ya.deg.graph)$name[order(ya.deg.deg,decreasing=T)]
ya.dd.driv <- intersect(V(ya.deg.graph)$name,rap.drivers)
ya.dd.nondriv <- V(ya.deg.graph)$name[!(V(ya.deg.graph)$name %in% rap.drivers)]

ks.test(ya.deg.deg[ya.dd.driv],ya.deg.deg[ya.dd.nondriv],exact=FALSE,alternative="less")$p.value

ya.deg.edges <- get.edgelist(ya.deg.graph)
ya.deg.edges.name <- cbind(sapply(ya.deg.edges[,1],orf.to.name),sapply(ya.deg.edges[,2],orf.to.name))

write.table(ya.deg.edges.name,"YeastOutput/Cytoscape/PPInet.txt",col.names=T,row.names=F,quote=F,sep="\t")

ya.deg.name <- sapply(names(ya.deg.deg),orf.to.name)
write.table(cbind(ya.deg.name,ya.deg.deg),"YeastOutput/Cytoscape/PPInet_degree.txt",col.names=T,row.names=F,quote=F,sep="\t")

rap.driv.ppiord <- ya.ddord[ya.ddord %in% rap.drivers]
rap.driv.ppiord.deg <- ya.deg.deg[rap.driv.ppiord]
rap.driv.ppiord.names <- sapply(rap.driv.ppiord,orf.to.name)

write.table(cbind(rap.driv.ppiord.names[1:10],rap.driv.ppiord.deg[1:10]),"YeastOutput/Top10rapamycinDrivers_PPIdeg.txt",col.names=F,row.names=F,quote=F,sep="\t")

# -------------------- Collect all sets of drivers for yeast -----------------------

## Menadione and Hydrogen Peroxide
## data downloaded from the website http://depts.washington.edu/sfields/deletion/index.html
## on May 19, 2015

tucker.file <- readLines("YeastDrivers/Tucker_CompFunctGenom_2004/supplemental_table 2.txt")

tucker.cols <- strsplit(tucker.file[1],split="\t")[[1]]
tucker.rows <- NULL
tucker.mat <- NULL
for (i in 2:length(tucker.file))
{
	this.split <- strsplit(tucker.file[i],split="\t")[[1]]
	tucker.rows <- c(tucker.rows,this.split[1])
	this.num <- as.numeric(this.split[4:7])
	tucker.mat <- rbind(tucker.mat,this.num)
}

rownames(tucker.mat) <- tucker.rows
colnames(tucker.mat) <- tucker.cols[4:7]

tucker.men <- tucker.mat[,2]
tucker.men <- tucker.men[!(is.na(tucker.men))]

men.drivers <- names(tucker.men)[tucker.men>1.5]

tucker.per <- tucker.mat[,1]
tucker.per <- tucker.per[!(is.na(tucker.per))]

per.drivers <- names(tucker.per)[tucker.per>1.5]


## compare with raw data

men.file <- read.table("YeastDrivers/Tucker_CompFunctGenom_2004/MEN_rawdata.txt",header=T,sep="\t",fill=T,quote="\"")

men.unique <- unique(as.character(men.file[,1]))

men.drivers1 <- NULL
for (i in 1:length(men.unique))
{
	this.ratio <- men.file[men.file[,1] %in% men.unique[i],6]
	if (sum(!(is.na(this.ratio)))>0)
	{
		this.ratio <- mean(this.ratio[!(is.na(this.ratio))])
		if (this.ratio>1.5) men.drivers1 <- c(men.drivers1,men.unique[i])
	}
}


per.file <- read.table("YeastDrivers/Tucker_CompFunctGenom_2004/PER_rawdata.txt",header=T,sep="\t",fill=T,quote="\"")

per.unique <- unique(as.character(per.file[,1]))

per.drivers1 <- NULL
for (i in 1:length(per.unique))
{
	this.ratio <- per.file[per.file[,1] %in% per.unique[i],5]
	if (sum(!(is.na(this.ratio)))>0)
	{
		this.ratio <- mean(this.ratio[!(is.na(this.ratio))])
		if (this.ratio>1.5) per.drivers1 <- c(per.drivers1,per.unique[i])
	}
}

## Heat shock
## data from Gibney PNAS 2013 Tables 2,3,4 (pooled all the genes)
## and Jarolim Genes Genomes Genetics (G3) 2013, Supplementary File S1 (took all the genes)

gibney.file <- read.table("YeastDrivers/Gibney_PNAS_2013/Table2_table3_table4.txt",sep="\t",header=F)
hs.driv1 <- as.character(gibney.file[,1])

hs.res.stat <- t(read.table("YeastDrivers/Jarolim_G3_2013/StatResistant_JarolimG3_2013.txt"))[1,]
hs.sens.stat <- t(read.table("YeastDrivers/Jarolim_G3_2013/StatSensitive_JarolimG3_2013.txt"))[1,]
hs.res.exp <- t(read.table("YeastDrivers/Jarolim_G3_2013/ExpResistant_JarolimG3_2013.txt"))[1,]
hs.sens.exp <- t(read.table("YeastDrivers/Jarolim_G3_2013/ExpSensitive_JarolimG3_2013.txt"))[1,]

hs.drivers <- unique(c(hs.driv1,hs.res.stat,hs.sens.stat,hs.res.exp,hs.sens.exp))

## DTT  
## taken from Rand, Mol. Biol. Cell 2006, Table 1

dtt.drivers <- t(read.table("YeastDrivers/Rand_MolBiolCell_2006/DTTdrivers.txt"))[1,]


## diamide
## Thorpe PNAS 2004, Table 2, worksheet titled "Total"


thorpe.file <- read.table("YeastDrivers/Thorpe_PNAS_2004/Thorpe_PNAS2004_diamide_menadione_H2O2.txt",header=F,sep="\t",fill=T,quote="\"")

thorpe.cols <- c("Gene","ORF","Gene function","CHP","X1","Diamide","X2","H2O2","X3","LoaOOH","X4","Menadione","X5","Glycerol","Slow Growth")
colnames(thorpe.file) <- thorpe.cols

diam.drivers <- unique(as.character(thorpe.file[as.character(thorpe.file[,"Diamide"])!="","ORF"]))


## Sorbitol
## Downloaded from http://genomics.lbl.gov/YeastFitnessData/websitefiles/cel_index.html on May 19, 2015
## we use the cutoffs for statistical significance suggested in the supp info


sorb.file1 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen5_8_1.5M-SorbitolA_res.txt",header=F,sep="\t",quote="\"")
sorb.driv1 <- as.character(sorb.file1[sorb.file1[,2]>20,1])

sorb.file2 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen5_8_1.5M-SorbitolA_sen.txt",header=F,sep="\t",quote="\"")
sorb.driv2 <- as.character(sorb.file2[sorb.file2[,2]>20,1])

sorb.file3 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen5_51_sorbitolB_res.txt",header=F,sep="\t",quote="\"")
sorb.driv3 <- as.character(sorb.file3[sorb.file3[,2]>20,1])

sorb.file4 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen5_51_sorbitolB_sen.txt",header=F,sep="\t",quote="\"")
sorb.driv4 <- as.character(sorb.file4[sorb.file4[,2]>20,1])

sorb.file5 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen15_52_sorbitolB_res.txt",header=F,sep="\t",quote="\"")
sorb.driv5 <- as.character(sorb.file5[sorb.file5[,2]>100,1])

sorb.file6 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen15_52_sorbitolB_sen.txt",header=F,sep="\t",quote="\"")
sorb.driv6 <- as.character(sorb.file6[sorb.file6[,2]>100,1])

sorb.file7 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen20_29_1.5M-SorbitolA_res.txt",header=F,sep="\t",quote="\"")
sorb.driv7 <- as.character(sorb.file7[sorb.file7[,2]>100,1])

sorb.file8 <- read.table("YeastDrivers/Giaever_Nature_2002/Gen20_29_1.5M-SorbitolA_sen.txt",header=F,sep="\t",quote="\"")
sorb.driv8 <- as.character(sorb.file8[sorb.file8[,2]>100,1])

osm.drivers <- unique(c(sorb.driv1,sorb.driv2,sorb.driv3,sorb.driv4,sorb.driv5,sorb.driv6,sorb.driv7,sorb.driv8))


# ----------------- Running GMIT for all yeast datasets -----------------

## Gene expression downloaded from http://genome-www.stanford.edu/yeast_stress/data/rawdata/complete_dataset.txt  May 24, 2015

gasch.lines <- readLines("YeastInput/Gasch_complete_dataset.txt")

gasch.data <- NULL
gasch.genes <- NULL
for (i in 2:length(gasch.lines))
{
	if (i %% 200 == 0) print(i)
	this.split <- strsplit(gasch.lines[i],split="\t")[[1]]
	this.gene <- toupper(this.split[1])
	gasch.genes <- c(gasch.genes,this.gene)
	this.data <- this.split[4:length(this.split)]
	if (length(this.data)<173) 
		this.data <- c(this.data,rep("",173-length(this.data)))
	gasch.data <- rbind(gasch.data,this.data)	
}
rownames(gasch.data) <- gasch.genes

gasch.cols <- strsplit(gasch.lines[1],split="\t")[[1]]
gasch.cols <- gasch.cols[4:length(gasch.cols)]

gasch.data.na <- array(NA,dim=dim(gasch.data))
for (j in 1:nrow(gasch.data))
{
	if (j %% 200 ==0) print(j)
	for (i in 1:ncol(gasch.data))
		if (gasch.data[j,i]!="") gasch.data.na[j,i] <- as.numeric(gasch.data[j,i])
}

gasch.colna <- which(apply(gasch.data.na,2,function(x){sum(is.na(x))>0.8*nrow(gasch.data)}))

gasch.cols1 <- gasch.cols[setdiff(1:length(gasch.cols),gasch.colna)]

gasch.data.na1 <- gasch.data.na[,setdiff(1:length(gasch.cols),gasch.colna)]


library(impute)
gasch.dat <- impute.knn(as.matrix(gasch.data.na1))$data

write.table(gasch.dat,"YeastInput/Gasch_impute_data.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(gasch.cols1,"YeastInput/Gasch_impute_colnames.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(gasch.genes,"YeastInput/Gasch_impute_rownames.txt",col.names=F,row.names=F,quote=F,sep="\t")


gasch.dat <- read.table("YeastInput/Gasch_impute_data.txt")
gasch.genes <- t(read.table("YeastInput/Gasch_impute_rownames.txt"))[1,]
gasch.cond <- readLines("YeastInput/Gasch_impute_colnames.txt")
colnames(gasch.dat) <- gasch.cond
rownames(gasch.dat) <- gasch.genes

gasch.net.genes <- intersect(yeast.prior.genes,gasch.genes)


## hydrogen peroxide

h2o2.cols <- gasch.cond[grep("H2O2",gasch.cond)]
h2o2.cols <- h2o2.cols[1:10]
h2o2.dat <- gasch.dat[,h2o2.cols]

h2o2.disc <- DatDiscOut(gasch.net.genes,h2o2.dat)

write.table(h2o2.disc,"YeastInput/Peroxide_disc.txt",col.names=F,row.names=F,quote=F,sep="\t")

per.prior <- PriorMatrix(gasch.net.genes,h2o2.dat,mac.links)

write.table(per.prior,"YeastInput/Peroxide_prior.txt",col.names=F,row.names=F,quote=F,sep="\t")


per.gmit.net <- read.table("YeastOutput/Peroxide_GMITnet.txt",skip=1)

per.gmit.edges <- NetworkReadGMIT(gasch.net.genes,per.gmit.net,mac.tfs.orf)

NetworkAnalysis(per.gmit.edges,per.drivers,"YeastOutput/Peroxide/",T,mac.tfs.orf,ppi.deg)


## heat shock

hs.cols <- gasch.cond[1:29]
hs.dat <- gasch.dat[,hs.cols]

hs.disc <- DatDiscOut(gasch.net.genes,hs.dat)

write.table(hs.disc,"YeastInput/HeatShock_disc.txt",col.names=F,row.names=F,quote=F,sep="\t")

hs.prior <- PriorMatrix(gasch.net.genes,hs.dat,mac.links)

write.table(hs.prior,"YeastInput/HeatShock_prior.txt",col.names=F,row.names=F,quote=F,sep="\t")

hs.gmit.net <- read.table("YeastOutput/HeatShock_GMITnet.txt",skip=1)
hs.gmit.edges <- NetworkReadGMIT(gasch.net.genes,hs.gmit.net,mac.tfs.orf)

NetworkAnalysis(hs.gmit.edges,hs.drivers,"YeastOutput/HeatShock/",T,mac.tfs.orf,ppi.deg)


## sorbitol

osm.cols <- gasch.cond[78:90]
osm.dat <- gasch.dat[,osm.cols]

osm.disc <- DatDiscOut(gasch.net.genes,osm.dat)

write.table(osm.disc,"YeastInput/Sorbitol_disc.txt",col.names=F,row.names=F,quote=F,sep="\t")

osm.prior <- PriorMatrix(gasch.net.genes,osm.dat,mac.links)

write.table(osm.prior,"YeastInput/Sorbitol_prior.txt",col.names=F,row.names=F,quote=F,sep="\t")

osm.gmit.net <- read.table("YeastOutput/Sorbitol_GMITnet.txt",skip=1)

osm.gmit.edges <- NetworkReadGMIT(gasch.net.genes,osm.gmit.net,mac.tfs.orf)

NetworkAnalysis(osm.gmit.edges,osm.drivers,"YeastOutput/Sorbitol/",T,mac.tfs.orf,ppi.deg)



## diamide

diam.cols <- gasch.cond[grep("diamide",gasch.cond)]
diam.dat <- gasch.dat[,diam.cols]

diam.disc <- DatDiscOut(gasch.net.genes,diam.dat)

write.table(diam.disc,"YeastInput/Diamide_disc.txt",col.names=F,row.names=F,quote=F,sep="\t")

diam.prior <- PriorMatrix(gasch.net.genes,diam.dat,mac.links)

write.table(diam.prior,"YeastInput/Diamide_prior.txt",col.names=F,row.names=F,quote=F,sep="\t")

diam.gmit.net <- read.table("YeastOutput/Diamide_GMITnet.txt",skip=1)

diam.gmit.edges <- NetworkReadGMIT(gasch.net.genes,diam.gmit.net,mac.tfs.orf)

NetworkAnalysis(diam.gmit.edges,diam.drivers,"YeastOutput/Diamide/",T,mac.tfs.orf,ppi.deg)


## menadione

men.cols <- gasch.cond[grep("Menadione",gasch.cond)]
men.dat <- gasch.dat[,men.cols]

men.disc <- DatDiscOut(gasch.net.genes,men.dat)

write.table(men.disc,"YeastInput/Menadione_disc.txt",col.names=F,row.names=F,quote=F,sep="\t")

men.prior <- PriorMatrix(gasch.net.genes,men.dat,mac.links)

write.table(men.prior,"YeastInput/Menadione_prior.txt",col.names=F,row.names=F,quote=F,sep="\t")


men.gmit.net <- read.table("YeastOutput/Menadione_GMITnet.txt",skip=1)

men.gmit.edges <- NetworkReadGMIT(gasch.net.genes,men.gmit.net,mac.tfs.orf)


NetworkAnalysis(men.gmit.edges,men.drivers,"YeastOutput/Menadione/",T,mac.tfs.orf,ppi.deg)


## DTT

dtt.cols <- unique(gasch.cond[c(grep("DTT",gasch.cond),grep("dtt",gasch.cond))])

dtt.dat <- gasch.dat[,dtt.cols]

dtt.disc <- DatDiscOut(gasch.net.genes,dtt.dat)

write.table(dtt.disc,"YeastInput/DTT_disc.txt",col.names=F,row.names=F,quote=F,sep="\t")

dtt.prior <- PriorMatrix(gasch.net.genes,dtt.dat,mac.links)

write.table(dtt.prior,"YeastInput/DTT_prior.txt",col.names=F,row.names=F,quote=F,sep="\t")

dtt.gmit.net <- read.table("YeastOutput/DTT_GMITnet.txt",skip=1)

dtt.gmit.edges <- NetworkReadGMIT(gasch.net.genes,dtt.gmit.net,mac.tfs.orf)

NetworkAnalysis(dtt.gmit.edges,dtt.drivers,"YeastOutput/DTT/",T,mac.tfs.orf,ppi.deg)

# --------------------- Running PANDA for all yeast datasets -------------------

mac.links.panda <- mac.links[mac.links[,1] %in% gasch.net.genes,]
mac.links.panda <- mac.links.panda[mac.links.panda[,2] %in% gasch.net.genes,]

pan.prior <- cbind(mac.links.panda,rep(1,nrow(mac.links.panda)))

write.table(pan.prior,"YeastInput/MotifPrior_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")


## heat shock

hs.dat.pan <- cbind(rownames(hs.dat),hs.dat)
write.table(hs.dat.pan,"YeastInput/HeatShock_dat_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")

hs.dat.ran <- RandomizeData(hs.dat)
hs.dat.ran1 <- cbind(rownames(hs.dat),hs.dat.ran)

write.table(hs.dat.ran1,"YeastInput/HeatShock_dat_PANDA_random.txt",col.names=F,row.names=F,quote=F,sep="\t")


hs.pan <- read.table("YeastOutput/HeatShock_PANDA_Output_FinalNetwork.pairs",header=F,sep="\t")
hs.pan.ran <- read.table("YeastOutput/HeatShock_PANDA_random_Output_FinalNetwork.pairs",header=F,sep="\t")

hs.pan.edges <- NetworkReadPanda(hs.pan,hs.pan.ran)

NetworkAnalysis(hs.pan.edges,hs.drivers,"YeastOutput/HeatShockPANDA/",T,mac.tfs.orf,ppi.deg)


## hydrogen peroxide

per.dat.pan <- cbind(rownames(h2o2.dat),h2o2.dat)
write.table(per.dat.pan,"YeastInput/Peroxide_dat_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")

per.dat.ran <- RandomizeData(h2o2.dat)
per.dat.ran1 <- cbind(rownames(h2o2.dat),per.dat.ran)

write.table(per.dat.ran1,"YeastInput/Peroxide_dat_PANDA_random.txt",col.names=F,row.names=F,quote=F,sep="\t")


per.pan <- read.table("YeastOutput/Peroxide_PANDA_Output_FinalNetwork.pairs",header=F,sep="\t")
per.pan.ran <- read.table("YeastOutput/Peroxide_PANDA_random_Output_FinalNetwork.pairs",header=F,sep="\t")

per.pan.edges <- NetworkReadPanda(per.pan,per.pan.ran)

NetworkAnalysis(per.pan.edges,per.drivers,"YeastOutput/PeroxidePANDA/",T,mac.tfs.orf,ppi.deg)


## menadione

men.dat.pan <- cbind(rownames(men.dat),men.dat)
write.table(men.dat.pan,"YeastInput/Menadione_dat_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")

men.dat.ran <- RandomizeData(men.dat)
men.dat.ran1 <- cbind(rownames(men.dat),men.dat.ran)

write.table(men.dat.ran1,"YeastInput/Menadione_dat_PANDA_random.txt",col.names=F,row.names=F,quote=F,sep="\t")


men.pan <- read.table("YeastOutput/Menadione_PANDA_Output_FinalNetwork.pairs",header=F,sep="\t")
men.pan.ran <- read.table("YeastOutput/Menadione_PANDA_random_Output_FinalNetwork.pairs",header=F,sep="\t")

men.pan.edges <- NetworkReadPanda(men.pan,men.pan.ran)

NetworkAnalysis(men.pan.edges,men.drivers,"YeastOutput/MenadionePANDA/",T,mac.tfs.orf,ppi.deg)

## DTT

dtt.dat.pan <- cbind(rownames(dtt.dat),dtt.dat)
write.table(dtt.dat.pan,"YeastInput/DTT_dat_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")

dtt.dat.ran <- RandomizeData(dtt.dat)
dtt.dat.ran1 <- cbind(rownames(dtt.dat),dtt.dat.ran)

write.table(dtt.dat.ran1,"YeastInput/DTT_dat_PANDA_random.txt",col.names=F,row.names=F,quote=F,sep="\t")

dtt.pan <- read.table("YeastOutput/DTT_PANDA_Output_FinalNetwork.pairs",header=F,sep="\t")
dtt.pan.ran <- read.table("YeastOutput/DTT_PANDA_random_Output_FinalNetwork.pairs",header=F,sep="\t")

dtt.pan.edges <- NetworkReadPanda(dtt.pan,dtt.pan.ran)

NetworkAnalysis(dtt.pan.edges,dtt.drivers,"YeastOutput/DTTPANDA/",T,mac.tfs.orf,ppi.deg)

## Sorbitol

osm.dat.pan <- cbind(rownames(osm.dat),osm.dat)
write.table(osm.dat.pan,"YeastInput/Sorbitol_dat_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")

osm.dat.ran <- RandomizeData(osm.dat)
osm.dat.ran1 <- cbind(rownames(osm.dat),osm.dat.ran)

write.table(osm.dat.ran1,"YeastInput/Sorbitol_dat_PANDA_random.txt",col.names=F,row.names=F,quote=F,sep="\t")

osm.pan <- read.table("YeastOutput/Sorbitol_PANDA_Output_FinalNetwork.pairs",header=F,sep="\t")
osm.pan.ran <- read.table("YeastOutput/Sorbitol_PANDA_random_Output_FinalNetwork.pairs",header=F,sep="\t")

osm.pan.edges <- NetworkReadPanda(osm.pan,osm.pan.ran)

NetworkAnalysis(osm.pan.edges,osm.drivers,"YeastOutput/SorbitolPANDA/",T,mac.tfs.orf,ppi.deg)


## Diamide

diam.dat.pan <- cbind(rownames(diam.dat),diam.dat)
write.table(diam.dat.pan,"YeastInput/Diamide_dat_PANDA.txt",col.names=F,row.names=F,quote=F,sep="\t")

diam.dat.ran <- RandomizeData(diam.dat)
diam.dat.ran1 <- cbind(rownames(diam.dat),diam.dat.ran)

write.table(diam.dat.ran1,"YeastInput/Diamide_dat_PANDA_random.txt",col.names=F,row.names=F,quote=F,sep="\t")

diam.pan <- read.table("YeastOutput/Diamide_PANDA_Output_FinalNetwork.pairs",header=F,sep="\t")
diam.pan.ran <- read.table("YeastOutput/Diamide_PANDA_random_Output_FinalNetwork.pairs",header=F,sep="\t")

diam.pan.edges <- NetworkReadPanda(diam.pan,diam.pan.ran)

NetworkAnalysis(diam.pan.edges,diam.drivers,"YeastOutput/DiamidePANDA/",T,mac.tfs.orf,ppi.deg)



# ------------- Read in viral oncogene data --------------

dat <- read.table("HumanInput/BatchAll_mich_ComBat.txt")

full.samples <- t(read.table("HumanInput/ListofSamples_geo_6-2-12.txt"))[1,]
sample.nums <- c(1:7,14:26,33:38,40,44,45,47:51,53:60,70,74,75,77:79,81,82,85:88,91:93,96,98,99,102,107,109:115,117:125,127,129:132,135,136)

transform <- full.samples[c(5:12,14,15,17:20,25:28,31:36,42,43,52,53,55,56,58,59,61,64,66:77)]
non.transform <- full.samples[!(full.samples %in% transform)]
controls <- full.samples[c(1:4,21:24,41,51,60,65)]


transf.nums <- sample.nums[c(5:12,14,15,17:20,25:28,31:36,42,43,52,53,55,56,58,59,61,64,66:77)]
nontr.nums <- sample.nums[!(sample.nums %in% transf.nums)]
ctrl.nums <- sample.nums[c(1:4,21:24,41,51,60,65)]

dat.cols <- colnames(dat)
dat.cols.nums <- sapply(dat.cols,function(x){
	this.split <- strsplit(x,split="_")[[1]][1]
	this.split <- strsplit(this.split,split="")[[1]]
	if (this.split[length(this.split)-1]==".")
	paste(this.split[2:(length(this.split)-3)],collapse="") else
	paste(this.split[2:(length(this.split)-1)],collapse="")
	})
dat.cols.nums <- as.numeric(dat.cols.nums)

mich.annot <- read.table("HumanInput/ViralOncogene_annotation.txt")
rownames(mich.annot) <- mich.annot[,1]
mich.annot1 <- mich.annot[!is.na(mich.annot[,2]),]


cegs.univ <- as.character(mich.annot[rownames(dat),3])
cegs.univ <- unique(cegs.univ[!(is.na(cegs.univ))])

## load in motif prior

tftarget1 <- read.table("HumanInput/cegs_tf_to_gene_11-Jan-2012_25000_relaxed.txt",sep="\t",quote="\"")
tftarget1 <- tftarget1[!(is.na(tftarget1[,1])),]

all.tfs.motif1 <- as.character(unique(tftarget1[,1]))
all.tfs.motif1 <- all.tfs.motif1[!(is.na(all.tfs.motif1))]

viral.prior.genes <- unique(c(all.tfs.motif1,as.character(tftarget1[,2])))

enrich.tab1 <- read.table("HumanInput/cegs_tf_cluster_enrich_11-Jan-2012_ppi_express_0.1_1.25_1unique.txt",sep="\t")
cl.tfs1 <- list(NULL)
for (i in 1:dim(enrich.tab1)[1])
{
	this.tfs <- strsplit(as.character(enrich.tab1[i,2]),split=" ")[[1]]
	if (length(this.tfs)>0)
	{
		this.tfs.exp <- NULL
		this.tfs.ppi <- NULL
		for (j in 1:length(this.tfs))
		{
			this.tf <- strsplit(this.tfs[j],split="")[[1]]
			if (this.tf[length(this.tf)]=="^") {this.tfs.exp <- c(this.tfs.exp, paste(this.tf[1:length(this.tf)-1],collapse=""))} else this.tfs.ppi <- c(this.tfs.ppi, paste(this.tf[1:length(this.tf)-1],collapse=""))
		}
		cl.tfs1[[i]] <- list(this.tfs.ppi,this.tfs.exp)	
	} else cl.tfs1[[i]] <- list(this.tfs,this.tfs)
}
path.tfs1 <- unique(unlist(cl.tfs1))

viral.enr.tfs <- unique(unlist(sapply(cl.tfs1,function(x){x[[2]]})))

## Sanger cancer gene census downloaded on May 25, 2015

sanger.genes <- t(read.table("HumanInput/Sanger_5-25-15.txt"))[1,]

tf.sanger <- intersect(sanger.genes,all.tfs.motif1)

a <- length(intersect(viral.enr.tfs,tf.sanger))
b <- length(viral.enr.tfs)-a
c <- length(tf.sanger) -a
d <- length(all.tfs.motif1)-a-b-c
fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")

venn.diagram(list(B=viral.enr.tfs,A=tf.sanger),fill=c("purple","red"),alpha=c(0.5,0.5),cex=0,filename = "HumanOutput/Fig5A_Venn_EnrichTFs_Drivers.tiff",ext.text=F,cat.cex=0,rotation.degree=180)


# --------------- Run PANDA ----------------

viral.net.genes <- intersect(cegs.univ,viral.prior.genes)

dat.pan <- NULL
for (i in 1:length(viral.net.genes))
{
	if (i %% 200 ==0) print(i)
	this.ids <- as.character(mich.annot1[as.character(mich.annot1[,3]) %in% viral.net.genes[i],1])
	this.dat <- dat[this.ids,]
	if (length(this.ids)>1)
		this.dat <- apply(this.dat,2,median)
	dat.pan <- rbind(dat.pan,as.numeric(this.dat))
}

dp.trans <- cbind(viral.net.genes,dat.pan[,dat.cols.nums %in% transf.nums])
dp.ctrl <- cbind(viral.net.genes,dat.pan[,dat.cols.nums %in% ctrl.nums])

write.table(dp.trans,"HumanInput/Viral_transf_dat.txt",col.names=F,row.names=F,quote=F,sep="\t")
write.table(dp.ctrl,"HumanInput/Viral_ctrl_dat.txt",col.names=F,row.names=F,quote=F,sep="\t")

prior.tfs <- intersect(all.tfs.motif1,viral.net.genes)

pan.prior <- tftarget1[tftarget1[,1] %in% prior.tfs,]
pan.prior <- pan.prior[pan.prior[,2] %in% viral.net.genes,1:2]

pan.prior <- cbind(pan.prior,rep(1,nrow(pan.prior)))

write.table(pan.prior,"HumanInput/Viral_motifprior.txt",col.names=F,row.names=F,quote=F,sep="\t")

viral.pan <- read.table("HumanOutput/Viral_transf_Out_FinalNetwork.pairs",header=F,sep="\t")
viral.pan.ctrl <- read.table("HumanOutput/Viral_ctrl_Out_FinalNetwork.pairs",header=F,sep="\t")

viral.pan.edges <- NetworkReadPanda(viral.pan,viral.pan.ctrl)

viral.graph <- graph.edgelist(as.matrix(viral.pan.edges),directed=F)
viral.deg <- degree(viral.graph)

viral.deg.ord <- V(viral.graph)$name[order(viral.deg,decreasing=T)]
viral.deg.tfs <- viral.deg.ord[viral.deg.ord %in% all.tfs.motif1]
length(intersect(viral.deg.tfs[1:50],sanger.genes))

viral.deg.score <- rep(0,length(all.tfs.motif1))
names(viral.deg.score) <- all.tfs.motif1
deg.tfs <- intersect(names(viral.deg),all.tfs.motif1)
viral.deg.score[deg.tfs] <- viral.deg[deg.tfs]

tf.sanger <- intersect(sanger.genes,all.tfs.motif1)
tf.nonsanger <- all.tfs.motif1[!(all.tfs.motif1 %in% sanger.genes)] 

viral.deg.stats <- NULL
viral.deg.stats[[1]] <- wilcox.test(viral.deg.score[tf.sanger],viral.deg.score[tf.nonsanger],exact=FALSE,alternative="greater")$p.value
viral.deg.stats[[2]] <- ks.test(viral.deg.score[tf.sanger],viral.deg.score[tf.nonsanger],exact=FALSE,alternative="less")$p.value


a <- length(intersect(viral.deg.tfs[1:29],tf.sanger))
b <- 29-a
c <- length(tf.sanger) -a
d <- length(all.tfs.motif1)-a-b-c
viral.deg.stats[[3]] <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value 



## Differential expression

dp.trans <- read.table("HumanInput/Viral_transf_dat.txt")
dp.ctrl <- read.table("HumanInput/Viral_ctrl_dat.txt")
viral.limma <- cbind(dp.trans[,2:ncol(dp.trans)],dp.ctrl[,2:ncol(dp.ctrl)])
rownames(viral.limma) <- dp.trans[,1]


library(limma)

dat.limma <- viral.limma

lev <- c("trnf","ctrl")
f <- factor(c(rep(lev[1],242),rep(lev[2],59)),levels=lev)

design <- model.matrix(~0+f)
colnames(design) <- lev
fit <- lmFit(dat.limma,design)

cont.diff <- makeContrasts(trnf - ctrl, levels=design)
fit2 <- contrasts.fit(fit,cont.diff)
fit2 <- eBayes(fit2)

viral.tab <- topTable(fit2,adjust="BH",number=nrow(dat.limma))

limma.tfs <- rownames(viral.tab)[rownames(viral.tab) %in% all.tfs.motif1]
limma.scores <- viral.tab[rownames(viral.tab) %in% all.tfs.motif1,6]

a <- length(intersect(limma.tfs[1:29],tf.sanger))
b <- 29-a
c <- length(tf.sanger) -a
d <- length(all.tfs.motif1)-a-b-c
fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value 


write.table(cbind(limma.tfs[1:30],viral.tab[limma.tfs[1:30],5]),"HumanOutput/DifferentiallyExpressedTFs.txt",col.names=F,row.names=F,quote=F,sep="\t")

viral.limma.score <- rep(0,length(all.tfs.motif1))
names(viral.limma.score) <- all.tfs.motif1
viral.limma.score[limma.tfs] <- limma.scores

viral.limma.stats <- NULL
viral.limma.stats[[1]] <- wilcox.test(viral.limma.score[tf.sanger],viral.limma.score[tf.nonsanger],exact=FALSE,alternative="greater")$p.value
viral.limma.stats[[2]] <- ks.test(viral.limma.score[tf.sanger],viral.limma.score[tf.nonsanger],exact=FALSE,alternative="less")$p.value



## Human PPI network downloaded from http://interactome.dfci.harvard.edu/H_sapiens/download/HI-II-14.tsv

vidal.ppi <- read.table("HumanInput/HI-II-14.txt",header=T,sep="\t")
vidal.net <- as.matrix(vidal.ppi[,c(2,4)])

vidal.graph <- graph.edgelist(vidal.net,directed=F)
vidal.deg <- degree(vidal.graph)


viral.ppi.score <- rep(0,length(all.tfs.motif1))
names(viral.ppi.score) <- all.tfs.motif1
vidal.tfs <- intersect(names(vidal.deg),all.tfs.motif1)

viral.ppi.score[vidal.tfs] <- vidal.deg[vidal.tfs]

viral.ppi.stats <- NULL
viral.ppi.stats[[1]] <- wilcox.test(viral.ppi.score[tf.sanger],viral.ppi.score[tf.nonsanger],exact=FALSE,alternative="greater")$p.value
viral.ppi.stats[[2]] <- ks.test(viral.ppi.score[tf.sanger],viral.ppi.score[tf.nonsanger],exact=FALSE,alternative="less")$p.value

viral.ppi.ord <- all.tfs.motif1[order(viral.ppi.score,decreasing=T)]

a <- length(intersect(viral.ppi.ord[1:29],tf.sanger))
b <- 29-a
c <- length(tf.sanger) -a
d <- length(all.tfs.motif1)-a-b-c
viral.ppi.stats[[3]] <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value 

## combining the PPI and TRX networks

viral.ppi.rank <- rank(viral.ppi.score)
viral.deg.rank <- rank(viral.deg.score)
viral.comb.rank <- apply(cbind(viral.ppi.rank,viral.deg.rank),1,mean)

viral.comb.ord <- all.tfs.motif1[order(viral.comb.rank,decreasing=T)]

viral.comb.stats  <- NULL
viral.comb.stats[[1]] <- wilcox.test(viral.comb.rank[tf.sanger],viral.comb.rank[tf.nonsanger],exact=FALSE,alternative="greater")$p.value
viral.comb.stats[[2]] <- ks.test(viral.comb.rank[tf.sanger],viral.comb.rank[tf.nonsanger],exact=FALSE,alternative="less")$p.value

a <- length(intersect(viral.comb.ord[1:29],tf.sanger))
b <- 29-a
c <- length(tf.sanger) -a
d <- length(all.tfs.motif1)-a-b-c
viral.comb.stats[[3]] <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value ## 0.066


par(mfrow=c(1,3),mar=c(3,3,3,3))
barplot(-1*log2(viral.deg.stats),ylim=c(0,8),col="violetred3",cex.axis=2)
lines(c(0,6),rep(-1*log2(0.05),2),lty=2,lwd=3,col="black")
barplot(-1*log2(viral.ppi.stats),ylim=c(0,8),col="deepskyblue2",cex.axis=2)
lines(c(0,6),rep(-1*log2(0.05),2),lty=2,lwd=3,col="black")
barplot(-1*log2(viral.comb.stats),ylim=c(0,8),col="salmon1",cex.axis=2)
lines(c(0,6),rep(-1*log2(0.05),2),lty=2,lwd=3,col="black")


NetworkAnalysis(viral.pan.edges,sanger.genes,"HumanOutput/Viral/",T,all.tfs.motif1,vidal.deg.tab)


## Final table of top-ranked TFs

viral.comb.ord[1:30]





