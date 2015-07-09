DatDiscOut <- function(net.genes,exp.dat){
	
	## Input: net.genes is a vector of ORFs that are present in both the prior network and the data
	## Input: exp.dat is a matrix of gene expression values which have been log-normalized and whose rownames are ORFs (including all net.genes); colnames are not important
	## Output: discretized data that can be written to a file and inputted directly into GlobalMIT
	
	exp.dat.net <- exp.dat[net.genes,]
	
	exp.disc <- array(2,dim=c(length(net.genes),ncol(exp.dat.net)))
	for (i in 1:length(net.genes))
	{
		if (i %% 1000 ==0) print(i)
		this.dat <- as.numeric(exp.dat.net[i,])
		exp.disc[i,this.dat>quantile(this.dat,0.66666)] <- 3
		exp.disc[i,this.dat<quantile(this.dat,0.33333)] <- 1
	}
	
	exp.disc1 <- rbind(t(exp.disc),t(exp.disc))
	exp.out <- apply(exp.disc1,1,function(x){paste(x,collapse="\t")})
	exp.out <- c(paste(c(ncol(exp.disc),nrow(exp.disc)),collapse="\t"),exp.out)
		
	exp.out
}


PriorMatrix <- function(net.genes,exp.dat,prior.links){
	
	mac.innet <- prior.links[prior.links[,1] %in% net.genes,]
	mac.innet <- mac.innet[mac.innet[,2] %in% net.genes,]
	
	prior.net <- array(0,dim=c(length(net.genes),length(net.genes)))
	rownames(prior.net) <- colnames(prior.net) <- net.genes
	for (i in 1:length(net.genes))
	{
		if (i %% 1000 ==0) print(i)
		this.parents <- mac.innet[mac.innet[,2] %in% net.genes[i],1]
		if (length(this.parents)>10)
		{
			par.dat <- exp.dat[this.parents,]
			targ.dat <- exp.dat[net.genes[i],]
			
			this.cor <- cor(t(par.dat),as.numeric(targ.dat))
			fin.par <- this.parents[order(abs(this.cor),decreasing=T)][1:10]
			prior.net[fin.par,i] <- 1
		} else prior.net[this.parents,i] <- 1
	}
		
	prior.net
}

RandomizeData <- function(dat.mat){
	dat.ran <- NULL
	for (i in 1:nrow(dat.mat))
	{
		if (i %% 500 ==0) print(i)
		this.ran <- sample(dat.mat[i,],ncol(dat.mat),replace=F)
		dat.ran <- rbind(dat.ran,this.ran)
	}
	dat.ran <- dat.ran[sample(1:nrow(dat.ran),nrow(dat.ran),replace=F),]
	dat.ran
}

NetworkReadGMIT <- function(net.genes,net.mat,all.tfs){
	
	net.tfs <- intersect(all.tfs,net.genes)
	
	net.links <- NULL
	for (i in 1:length(net.tfs))
	{
		if (i %% 50 ==0) print(i)
		this.targ <- net.genes[net.mat[which(net.genes %in% net.tfs[i]),]>0]
		if (length(this.targ)>0)
		net.links <- rbind(net.links,cbind(rep(net.tfs[i],length(this.targ)),this.targ))
	}
		
	net.links
}

NetworkReadPanda <- function(net.pan,net.ran.pan){
	s1.pan <- net.pan
	s2.pan <- net.ran.pan

	pan.weights <- cbind(s1.pan[,4],s2.pan[,4])
	s1.score <- apply(pan.weights,1,function(x){pnorm(x[1])*pnorm(x[1]-x[2])})
	s2.score <- apply(pan.weights,1,function(x){pnorm(x[2])*pnorm(x[2]-x[1])})
	s1.spec <- s1.pan[s1.score>0.8,1:2]
	s1.spec
}


NetworkAnalysis <- function(net.links,driver.genes,folder.name,plot.logical,all.tfs,ppi.deg){

	mac.tfs <- all.tfs
	
	net.graph <- graph.edgelist(as.matrix(net.links),directed=F)
	net.deg <- degree(net.graph)
		
	ppi.deg.names <- as.character(ppi.deg[,1])
	ppi.deg <- as.numeric(ppi.deg[,2])
	names(ppi.deg) <- ppi.deg.names
	
	tf.ppi.deg <- ppi.deg[names(ppi.deg) %in% mac.tfs]
	ppi.deg.norm <- tf.ppi.deg/max(tf.ppi.deg)
	ppi.deg.norm1 <- rep(0,length(mac.tfs))
	names(ppi.deg.norm1) <- mac.tfs
	ppi.deg.norm1[names(ppi.deg.norm)] <- ppi.deg.norm

	tf.trx.deg <- net.deg[names(net.deg) %in% mac.tfs]
	trx.deg.norm <- tf.trx.deg/max(tf.trx.deg)
	trx.deg.norm1 <- rep(0,length(mac.tfs))
	names(trx.deg.norm1) <- mac.tfs
	trx.deg.norm1[names(trx.deg.norm)] <- trx.deg.norm
	
	ppi.rank <- rank(ppi.deg.norm1)
	trx.rank <- rank(trx.deg.norm1)

	ppi.trx.deg <- apply(cbind(trx.rank,ppi.rank),1,mean)

	all.scores <- array(0,dim=c(3,length(mac.tfs)))
	colnames(all.scores) <- mac.tfs
	all.scores[1,names(net.deg)[names(net.deg) %in% mac.tfs]] <- net.deg[names(net.deg) %in% mac.tfs]
	all.scores[2,names(ppi.deg)[names(ppi.deg) %in% mac.tfs]] <- ppi.deg[names(ppi.deg) %in% mac.tfs]
	all.scores[3,] <- ppi.trx.deg
	
	score.file <- paste(c(folder.name,"AllScores.txt"),collapse="")
	write.table(all.scores,score.file,row.names=T,col.names=T,quote=F,sep="\t")
	
	tpr.tab <- array(1,dim=dim(all.scores))
	fpr.tab <- array(1,dim=dim(all.scores))
	for (i in 1:nrow(all.scores))
	{
		this.score <- all.scores[i,]
		this.ord <- mac.tfs[order(this.score,decreasing=T)]
		for (j in 1:length(this.ord))
		{
			tp <- length(intersect(this.ord[1:j],driver.genes))
			fp <- j-tp
			fn <- length(intersect(driver.genes,mac.tfs))-tp
			tn <- length(mac.tfs) -tp-fp-fn
		
			tpr.tab[i,j] <- tp/(tp+fn)
			fpr.tab[i,j] <- fp/(fp+tn)
		}
	}

	auc.vec <- apply(tpr.tab,1,sum)/ncol(tpr.tab)
	
	tf.driv <- intersect(driver.genes,mac.tfs)
	tf.nondriv <- mac.tfs[!(mac.tfs %in% driver.genes)]
	
	ks.pvals <- NULL
	for (i in 1:nrow(all.scores))
	ks.pvals <- c(ks.pvals,ks.test(all.scores[i,tf.driv],all.scores[i,tf.nondriv],exact=FALSE,alternative="less")$p.value)
	
	ks.round <- sapply(ks.pvals,function(x){round(x,digits=2)})
	auc.round <- sapply(auc.vec,function(x){round(x,digits=2)})
	
	plot.colors <- c("violetred3","deepskyblue2","salmon1")
	
	roc.file <- paste(c(folder.name,"ROC_combined.pdf"),collapse="")
	fig.file <- paste(c(folder.name,"ROC_combined_figure.pdf"),collapse="")
	
	if (plot.logical){
	pdf(roc.file,width=5,height=5)
	plot(fpr.tab[1,],tpr.tab[1,],col="white",xlab="False positive rate",ylab="True positive rate",cex.lab=1.5,cex.axis=1.8)
	for (i in 1:nrow(all.scores))
	{
		lines(fpr.tab[i,],tpr.tab[i,],col=plot.colors[i],lwd=4)
		text(0.45,0.15-(i-1)*0.06,paste(c("P = ",ks.round[i],"  AUC = ",auc.round[i]),collapse=""), col=plot.colors[i],pos=4,cex=1.2)
	}
	lines(c(0,1),c(0,1),lty=2,col="black",lwd=4)
#legend("bottomright",col=rainbow(nrow(all.scores)),pch=20,c("TRXdeg","PPIdeg","TRX+PPI"))
	dev.off()
	
	pdf(fig.file,width=5,height=5)
	plot(fpr.tab[1,],tpr.tab[1,],col="white",xlab="False positive rate",ylab="True positive rate",cex.lab=1.5,cex.axis=1.8)
	for (i in 1:nrow(all.scores))
	{
		lines(fpr.tab[i,],tpr.tab[i,],col=plot.colors[i],lwd=4)
		##text(0.5,0.15-(i-1)*0.06,paste(c("P = ",ks.round[i],"  AUC = ",auc.round[i]),collapse=""), col=plot.colors[i],pos=4,cex=1.2)
	}
	lines(c(0,1),c(0,1),lty=2,col="black",lwd=4)
#legend("bottomright",col=rainbow(nrow(all.scores)),pch=20,c("TRXdeg","PPIdeg","TRX+PPI"))
	dev.off()
	}

	## odds ratios for intersections with top 20% (23) TFs
	tf.ord.or <- NULL
	tf.ord.p <- NULL
	tf.ord.spec <- NULL
	for (i in 1:nrow(all.scores))
	{
		this.ord <- mac.tfs[order(all.scores[i,],decreasing=T)]
		
		a <- length(intersect(this.ord[1:23],tf.driv))
		b <- length(tf.driv)-a
		c <- 23-a
		d <- length(mac.tfs)-a-b-c
		
		this.p <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
		this.or <- (a*d)/(b*c)
		tf.ord.p <- c(tf.ord.p,this.p)
		tf.ord.or <- c(tf.ord.or,this.or)
		tf.ord.spec <- c(tf.ord.spec,a/23)
	}
	
	bar.file <- paste(c(folder.name,"Barplot_OR_pval.pdf"),collapse="")
	pdf(bar.file,width=2.5,height=3)
	barplot(tf.ord.or,col=plot.colors)
	dev.off()
	
	## odds ratios for intersections with top 50 TFs
	tf.ord.or1 <- NULL
	tf.ord.p1 <- NULL
	tf.ord.spec1 <- NULL
	for (i in 1:nrow(all.scores))
	{
		this.ord <- mac.tfs[order(all.scores[i,],decreasing=T)]
		
		a <- length(intersect(this.ord[1:50],tf.driv))
		b <- length(tf.driv)-a
		c <- 50-a
		d <- length(mac.tfs)-a-b-c
		
		this.p <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
		this.or <- (a*d)/(b*c)
		tf.ord.p1 <- c(tf.ord.p1,this.p)
		tf.ord.or1 <- c(tf.ord.or1,this.or)
		tf.ord.spec1 <- c(tf.ord.spec1,a/50)
	}

	
	## odds ratios for intersection of top 50 TFs
	trx.ord <- mac.tfs[order(all.scores[1,],decreasing=T)]
	ppi.ord <- mac.tfs[order(all.scores[2,],decreasing=T)]
	top.int <- intersect(trx.ord[1:50],ppi.ord[1:50])
	
	a <- length(intersect(top.int,driver.genes))
	b <- length(tf.driv)-a
	c <- length(top.int)-a
	d <- length(mac.tfs)-a-b-c
	
	int.p <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
	int.or <- (a*d)/(b*c)
	int.spec <- a/length(top.int)
	
	or.p.tab <- cbind(c(tf.ord.or,tf.ord.or1,int.or), c(tf.ord.p,tf.ord.p1,int.p),c(tf.ord.spec,tf.ord.spec1,int.spec))
	
	inter.file <- paste(c(folder.name,"All_OR_pval.txt"),collapse="")
	
	write.table(or.p.tab,inter.file,col.names=F,row.names=F,quote=F,sep="\t")
	
}



