loadBed <- function(file,fname) {
	bed<-read.table(file,sep="",stringsAsFactors=F)
	bed.col.n<-as.integer(dim(bed)[2])
	if (bed.col.n<3) stop("Bed file must contain at least 3 columns!")
	bed.cols<-c("chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
	colnames(bed)<-bed.cols[1:bed.col.n]
	bed$chromStart <- as.integer(bed$chromStart)
	bed$chromEnd <- as.integer(bed$chromEnd)
	comment(bed) <- paste(fname)
	bed$source<-paste(fname)
	bed$source<-factor(bed$source)
	bed$chrom <- factor(bed$chrom, levels = unique(bed$chrom))
	return(bed)
}

gplotBed <- function(bed) {
	theTheme<-theme(
		text = element_text(family="Arial"),
		plot.title = element_text(size=rel(2),hjust=0.5,face="bold"),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_rect(colour = "black", fill=NA, size=0.75),
		panel.background = element_blank(),
		axis.line = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(size = rel(1.8)),
		axis.text = element_text(size = rel(1.5)),
		axis.text.x = element_text(angle=45,hjust=1),
		legend.position='none',
		plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
	)
	nChr <- length(unique(bed$chrom))
	scoreStats <- summary(bed$score)
	maxScore <- scoreStats[6]
	minScore <- scoreStats[1]
	iqr <- scoreStats[5]-scoreStats[2]
	top3 <- scoreStats[5]+3*iqr
	low3 <- scoreStats[2]-3*iqr
    boxplot <- ggplot(bed,aes(chrom,score)) 
	if (maxScore > top3) {
		boxplot <- boxplot + geom_rect(
			xmin=1-nChr*0.03,
			ymin=top3,
			xmax=nChr+nChr*0.03,
			ymax=maxScore+(maxScore-minScore)*0.02,
			fill="gray95",
			size=0
		)
	}
	if (minScore < low3) {
		boxplot <- boxplot + geom_rect(
			xmin=1-nChr*0.03,
			ymin=minScore-(maxScore-minScore)*0.02,
			xmax=nChr+nChr*0.03,
			ymax=low3,
			fill="gray95",
			size=0
		)
	}
	boxplot <- boxplot + geom_boxplot()
	boxplot <- boxplot + coord_cartesian(xlim=c(1,nChr),ylim=c(minScore,maxScore)) 
	boxplot <- boxplot + ggtitle("Raw reads distribution")
	boxplot <- boxplot + labs(y="Score value")
	boxplot <- boxplot + theTheme
    return(boxplot)
}

rmMax <- function(bed) {
	bed<-bed[-which(bed$score==max(bed$score)),]
	return(bed)
}

rmChr <- function(bed,chr) {
	bed<-subset(bed,chrom!=chr)
	bed$chrom <- factor(bed$chrom, levels = unique(bed$chrom))
	return(bed)
}

rmOut <- function(bed) {
	scoreStats <- summary(bed$score)
	iqr <- scoreStats[5]-scoreStats[2]
	bed<-bed[!(bed$score<scoreStats[2]-3*iqr | bed$score>scoreStats[5]+3*iqr),]
	return(bed)
}

makeRatio <- function(bedRep,bedNonRep) {
	repSum <- sum(bedRep$score)
	nonRepSum <- sum(bedNonRep$score)
	corrFactor <- repSum/nonRepSum
	byVector <- c("chrom","chromStart","chromEnd")
	ratioBed <- merge(bedRep,bedNonRep,by=byVector,all=F,sort=F)
	mergedNames <- (colnames(ratioBed))
# rename columns
	mergedNames <- gsub(".x",".rep",mergedNames)
	mergedNames <- gsub(".y",".nonRep",mergedNames)
	names(ratioBed) <- mergedNames
# deduplicate repetitive columns
	byCommon <- intersect(colnames(bedRep),colnames(bedNonRep))
	byOther <- setdiff(byCommon,byVector)
	for (colName in byOther) {
		repCol <- paste0(colName,".rep")
		nonRepCol <- paste0(colName,".nonRep")
		if (all(as.character(ratioBed[[repCol]])==as.character(ratioBed[[nonRepCol]]))) {
			names(ratioBed)[names(ratioBed) == repCol] <- colName
			ratioBed[[nonRepCol]] <- NULL
		}
	}
#~ Calculate ratio
	ratioBed$ratio <- ratioBed$score.rep/ratioBed$score.nonRep/corrFactor
	ratioBed$chrom <- factor(ratioBed$chrom, levels = unique(ratioBed$chrom))
	bin <- ratioBed$chromEnd[1]-ratioBed$chromStart[1]
#~ Smoothing
	ratioBed$fit <- paste("NA")
	ratioBed$group <- paste("NA")
	j <- 1
	for (chr in levels(ratioBed$chrom)) {
		ratioBed[ratioBed$chrom==chr,"fit"] <- round(as.numeric(smooth.spline(subset(ratioBed,chrom==chr,select=c(chromStart,ratio)))$y),3)
#~ Grouping lines
		starts <- as.vector(unlist(subset(ratioBed,chrom==chr,select=chromStart)))
		group <- j
		i <- 1
		for (i in 1:(length(starts)-1)) {
			i <- i + 1
			if ( (starts[i] - starts[i-1])/bin > 9 ) { j <- j+1 }
			group <- append(group,j)
		}
		ratioBed[ratioBed$chrom==chr,"group"] <- as.integer(group)
		j <- j+1
	}
	ratioBed$fit<-as.numeric(ratioBed$fit)
#~ Removing small groups
	for (group in unique(ratioBed$group)) {
		tmp<-which(ratioBed$group == group)
		if (length(tmp) < 3) { ratioBed$group[tmp] <- NA }
	}
	ratioBed$group <- factor(ratioBed$group)
	return(ratioBed)
}

gplotRatio <- function(ratio) {
	ymax <- max(hist(ratio, breaks=100, plot=FALSE)$counts)
	theTheme<-theme(
		text = element_text(family="Arial"),
		plot.title = element_text(size=rel(2),hjust=0.5,face="bold"),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_rect(colour = "black", fill=NA, size=0.75),
		panel.background = element_blank(),
		axis.line = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(size = rel(1.8)),
		axis.text = element_text(size = rel(1.5)),
		legend.position='none',
		plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
	)
    hist <- ggplot() + aes(ratio)
    hist <- hist + annotate("rect", xmin=1, xmax=2, ymin=0, ymax=Inf, alpha=0.5, fill="lightgreen")
	hist <- hist + geom_histogram(bins=100,fill='gray90',color='gray30')
	hist <- hist + ggtitle("Ratio distribution")
	hist <- hist + labs(y="Counts")
	hist <- hist + theTheme
    return(hist)
}

normaliseRatio <- function(ratioBed) {
	f <- function (factor,values,ceiling=2) abs(sum(factor*values[factor*values>ceiling]-ceiling)-sum(1-factor*values[factor*values<1]))
	xmin <- optimise(f,c(1,2),tol = 0.001,values=ratioBed$fit,maximum=F)
	ratioFactor <- xmin[[1]]
	ratioBed$tmpFit <- ratioBed$fit*ratioFactor
	ratioBed$tmpRatio <- ratioBed$ratio*ratioFactor
	comment(ratioBed) <- as.character(round(ratioFactor,3))
	return(ratioBed)
}

plotGenome <- function(ratioDF,genome,circles=NULL,pointers=NULL,rectangles=NULL) {
	require(ggplot2)
	if (dim(ratioDF)[1] == 0) stop("There is no data to plot")
	window <- ratioDF$chromEnd[1]-ratioDF$chromStart[1]
	genomeName <- deparse(substitute(genome))
	if (!is.null(circles)) circlesName <- deparse(substitute(circles))
##  Reorder the data based on the genome dataframe
	if (genomeName == "sacCer3") {
		ratioDF$chrom<-factor(ratioDF$chrom,levels=unique(genome$chrom))
		ratioDF<-ratioDF[order(ratioDF$chrom,ratioDF$chromStart),]
		rownames(ratioDF) <- 1:nrow(ratioDF)
	}
##  Make nice x-axis labels
	tmp <- max(genome$chromEnd)
	n <- round(log10(tmp/10),digits=0)
	step <- as.integer(round(tmp/10,digits=-n))
	if (tmp/step < 10) { step <- as.integer(step/2) }
	ticks <- seq(from = 0, to = tmp, by = step)
	labels <- vector(mode="character",length=0)
	labels<-ifelse(log10(ticks)>=3 & log10(ticks)<6, paste0(ticks/1000,"kb"), ifelse(log10(ticks)>=6 & log10(ticks)<9,paste0(ticks/1000000,"Mb"), ifelse(log10(ticks)>=9 & log10(ticks)<12,paste0(ticks/1000000,"Gb"),ticks)))
	labels <- data.frame(ticks,labels)
##  initialising legendary vectors here
	color.names <- vector(mode="character",length=0)
	color.values <- vector(mode="character",length=0)
	color.lines <- vector(mode="character",length=0)
	color.shapes <- vector(mode="numeric",length=0)
	color.sizes <- vector(mode="numeric",length=0)
	color.fills<-  vector(mode="character",length=0)
	fill.names <- vector(mode="character",length=0)
	fill.values <- vector(mode="character",length=0)
	fill.shapes <- vector(mode="numeric",length=0)
	fill.sizes <- vector(mode="numeric",length=0)
	fill.colors <- vector(mode="character",length=0)
	fill.lines<- vector(mode="character",length=0)
	line.names <-vector(mode="character",length=0)
	line.values <-vector(mode="character",length=0)
	line.colors <- vector(mode="character",length=0)
	line.shapes <- vector(mode="numeric",length=0)
	line.sizes <- vector(mode="numeric",length=0)
##  plotting parameters
	plot.theme<-theme(
		text = element_text(family="Arial"),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		plot.title = element_text(hjust = 0.5,size=rel(1.8)),
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.x = element_line(colour="gray50",size=.5),
		axis.title = element_text(size = rel(1.6),face="bold"),
		axis.text.y = element_blank(),
		axis.text.x = element_text(size=rel(1.7)),
		legend.position="bottom",
		legend.title=element_blank(),
		legend.text = element_text(size=rel(1.2)),
		legend.key = element_blank(),
		strip.background = element_blank(),
		strip.text = element_blank(),
		plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
	)
	plot<-ggplot(ratioDF)
	plot<-plot+scale_y_continuous(name="Relative copy number",limits = c(0.8, 2.2))  ##  y scale
	plot<-plot+annotate(geom="segment",x=0, xend=0, y=1, yend=2, colour="gray50", lwd=0.5)  ##  y scale line
	plot<-plot+scale_x_continuous(breaks=labels$ticks,labels=labels$labels,name="Chromosome coordinates",limits=c(-window,NA),expand=c(0.01,0.01))  ##  x scale
	plot<-plot+facet_grid(chrom ~ .,scales = "free", space = "free_x")  ##  Facet
	plot<-plot+geom_text(aes(x=chromEnd+5*window,y=1.5,label=chrom),data=genome,angle=270,size=5,vjust=0)  ##  Chromosome labels
	plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=1,yend=1),data=genome,linetype="dashed",size=.25,color="gray40")  ##  Lower y line
	plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=2,yend=2),data=genome,linetype="dashed",size=.25,color="gray40")  ##  Upper y line
	plot<-plot+geom_text(aes(x=-window/2,y=1,label="1"),data=genome,size=4,hjust=1)  ##  Lower y break
	plot<-plot+geom_text(aes(x=-window/2,y=2,label="2"),data=genome,size=4,hjust=1)  ##  Upper y break
	plot<-plot+geom_segment(aes(x=0,xend=genome$chromEnd,y=0.9,yend=0.9),data=genome,size=.7)  ##  Chromosome length line
	if ("cen" %in% colnames(genome)) {
		plot<-plot+geom_segment(aes(x=genome$cen,xend=genome$cen,y=1,yend=2),data=genome,color="green3")
		plot<-plot+geom_vline(aes(xintercept=-step,linetype="Centromere"),data=genome)  ##  dummy
		line.names<-c(line.names,"Centromere")
		line.values<-c(line.values,"solid")
		line.colors<-c(line.colors,"green3")
		line.shapes<-c(line.shapes,NA)
		line.sizes<-c(line.sizes,1)
	}
	if (genomeName=="sacCer3") {  ##  If it's sacCer3, plot the origins
		plot <- plot+geom_point(aes(x=chromStart+window/2,y=0.9,fill=circlesName),data=circles,color="black",size=2,shape=21)
		fill.names<-c(fill.names,circlesName)  ##  Circles legend values
		fill.values<-c(fill.values,"white")
		fill.shapes<-c(fill.shapes,21)
		fill.sizes<-c(fill.sizes,4)
		fill.colors<-c(fill.colors,"black")
		fill.lines<-c(fill.lines,"blank")
	}
		my.fills<-fill.values
		names(my.fills)<-fill.names
		bed.src<-levels(ratioDF$source.rep)[1]  ##  Retrieve bed file name
		plot<-plot+ggtitle(paste0(bed.src," replication profile"))  ##  Plot title
		plot<-plot+geom_line(aes(x=chromStart+500,y=fit,group=group,color=paste0(bed.src," smooth")),size=1)  ##  Smoothed data
		color.names<-c(color.names,paste0(bed.src," smooth"))  ##  Smoothed data legend values
		color.values<-c(color.values,"gray25")
		color.lines<-c(color.lines,"solid")
		color.shapes<-c(color.shapes,NA)
		color.sizes<-c(color.sizes,1.5)
		color.fills<-c(color.fills,NA)
		plot<-plot+geom_point(aes(x=chromStart+500,y=ratio,color=paste0(bed.src," raw")), size=.5)  ##  Raw data
		color.names<-c(color.names,paste0(bed.src," raw"))  ##  Raw data legend values
		color.values<-c(color.values,"royalblue1")
		color.lines<-c(color.lines,"blank")
		color.shapes<-c(color.shapes,20)
		color.sizes<-c(color.sizes,4)
		color.fills<-c(color.fills,NA)
		my.colors<-color.values  ##  Preparing named vector for colors
		names(my.colors)<-color.names
		my.lines<-line.values
		names(my.lines)<-line.names
		##  Supplying the manual color, fill and linetype scales - this is required for the legend
		plot<-plot+scale_colour_manual(name="a",values=my.colors,guide = guide_legend(override.aes = list(
			linetype=color.lines,shape= color.shapes,size=color.sizes,fill=color.fills)),breaks=names(my.colors))
		plot<-plot+scale_fill_manual(name="b",values=my.fills,guide = guide_legend(override.aes = list(
			size=fill.sizes,shape=fill.shapes,linetype=fill.lines,color=fill.colors)),breaks = names(my.fills))
		plot<-plot+scale_linetype_manual(name="c",values=my.lines,guide = guide_legend(override.aes = list(
			color=line.colors,size=line.sizes,shape=line.shapes)),breaks = names(my.lines))
		##  Giving hope and saving the pdf
		plot<-plot+plot.theme
		return(plot)
}

plotGenomes <- function(mergedRatioDFs,genome,circles=NULL,pointers=NULL,rectangles=NULL) {
	require(ggplot2)
	window <- mergedRatioDFs$chromEnd[1]-mergedRatioDFs$chromStart[1]
	genomeName <- deparse(substitute(genome))
	if (!is.null(circles)) circlesName <- deparse(substitute(circles))
	firstName <- levels(mergedRatioDFs$source.rep.x)[1]
	secondName <- levels(mergedRatioDFs$source.rep.y)[1]
	if (!is.null(circles)) circlesName <- deparse(substitute(circles))
##  Reorder the data based on the genome dataframe
	if (genomeName == "sacCer3") {
		mergedRatioDFs$chrom<-factor(mergedRatioDFs$chrom,levels=unique(genome$chrom))
		mergedRatioDFs<-mergedRatioDFs[order(mergedRatioDFs$chrom,mergedRatioDFs$chromStart),]
		rownames(mergedRatioDFs) <- 1:nrow(mergedRatioDFs)
	}
##  Make nice x-axis labels
	tmp <- max(genome$chromEnd)
	n <- round(log10(tmp/10),digits=0)
	step <- as.integer(round(tmp/10,digits=-n))
	if (tmp/step < 10) { step <- as.integer(step/2) }
	ticks <- seq(from = 0, to = tmp, by = step)
	labels <- vector(mode="character",length=0)
	labels<-ifelse(log10(ticks)>=3 & log10(ticks)<6, paste0(ticks/1000,"kb"), ifelse(log10(ticks)>=6 & log10(ticks)<9,paste0(ticks/1000000,"Mb"), ifelse(log10(ticks)>=9 & log10(ticks)<12,paste0(ticks/1000000,"Gb"),ticks)))
	labels <- data.frame(ticks,labels)
	##  initialising legendary vectors here
	color.names <- vector(mode="character",length=0)
	color.values <- vector(mode="character",length=0)
	color.lines <- vector(mode="character",length=0)
	color.shapes <- vector(mode="numeric",length=0)
	color.sizes <- vector(mode="numeric",length=0)
	color.fills<-  vector(mode="character",length=0)
	fill.names <- vector(mode="character",length=0)
	fill.values <- vector(mode="character",length=0)
	fill.shapes <- vector(mode="numeric",length=0)
	fill.sizes <- vector(mode="numeric",length=0)
	fill.colors <- vector(mode="character",length=0)
	fill.lines<- vector(mode="character",length=0)
	line.names <-vector(mode="character",length=0)
	line.values <-vector(mode="character",length=0)
	line.colors <- vector(mode="character",length=0)
	line.shapes <- vector(mode="numeric",length=0)
	line.sizes <- vector(mode="numeric",length=0)
##  plotting parameters
	plot.theme<-theme(
		text = element_text(family="Arial"),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		plot.title = element_text(hjust = 0.5,size=rel(1.8)),
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.x = element_line(colour="gray50",size=.5),
		axis.title = element_text(size = rel(1.6),face="bold"),
		axis.text.y = element_blank(),
		axis.text.x = element_text(size=rel(1.7)),
		legend.position="bottom",
		legend.title=element_blank(),
		legend.text = element_text(size=rel(1.2)),
		legend.key = element_blank(),
		strip.background = element_blank(),
		strip.text = element_blank(),
		plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
	)
	plot<-ggplot(mergedRatioDFs)
	plot<-plot+scale_y_continuous(name="Relative copy number",limits = c(0.8, 2.2))  ##  y scale
	plot<-plot+annotate(geom="segment",x=0, xend=0, y=1, yend=2, colour="gray50", lwd=0.5)  ##  y scale line
	plot<-plot+scale_x_continuous(breaks=labels$ticks,labels=labels$labels,name="Chromosome coordinates",limits=c(-window,NA),expand=c(0.01,0.01))  ##  x scale
	plot<-plot+facet_grid(chrom ~ .,scales = "free", space = "free_x")  ##  Facet
	plot<-plot+geom_text(aes(x=chromEnd+5*window,y=1.5,label=chrom),data=genome,angle=270,size=5,vjust=0)  ##  Chromosome labels
	plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=1,yend=1),data=genome,linetype="dashed",size=.25,color="gray40")  ##  Lower y line
	plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=2,yend=2),data=genome,linetype="dashed",size=.25,color="gray40")  ##  Upper y line
	plot<-plot+geom_text(aes(x=-window/2,y=1,label="1"),data=genome,size=4,hjust=1)  ##  Lower y break
	plot<-plot+geom_text(aes(x=-window/2,y=2,label="2"),data=genome,size=4,hjust=1)  ##  Upper y break
	plot<-plot+geom_segment(aes(x=0,xend=genome$chromEnd,y=0.9,yend=0.9),data=genome,size=.7)  ##  Chromosome length line
	if ("cen" %in% colnames(genome)) {
		plot<-plot+geom_segment(aes(x=genome$cen,xend=genome$cen,y=1,yend=2),data=genome,color="green3")
		plot<-plot+geom_vline(aes(xintercept=-step,linetype="Centromere"),data=genome)  ##  dummy
		line.names<-c(line.names,"Centromere")
		line.values<-c(line.values,"solid")
		line.colors<-c(line.colors,"green3")
		line.shapes<-c(line.shapes,NA)
		line.sizes<-c(line.sizes,1)
	}
	if (!is.null(circles)) {
		plot <- plot+geom_point(aes(x=chromStart+window/2,y=0.9,fill=circlesName),data=circles,color="black",size=2,shape=21)
		fill.names<-c(fill.names,circlesName)  ##  Circles legend values
		fill.values<-c(fill.values,"white")
		fill.shapes<-c(fill.shapes,21)
		fill.sizes<-c(fill.sizes,4)
		fill.colors<-c(fill.colors,"black")
		fill.lines<-c(fill.lines,"blank")
	}
	my.fills<-fill.values
	names(my.fills)<-fill.names
	plot<-plot+ggtitle(paste0(firstName,"-",secondName," replication profiles"))  ##  Plot title

##  If there are p-values, get the significant ones
	temp<-which(mergedRatioDFs$p_value <= 0.001)  ##  p<=0.001
	signLevel3<-mergedRatioDFs[temp,]
	temp<-which(mergedRatioDFs$p_value <= 0.01 & mergedRatioDFs$p_value > 0.001)  ##  0.01 < p < 0.001
	signLevel2<-mergedRatioDFs[temp,]
##  Plot significant p-values
	if (dim(signLevel3)[1]!=0) {
		plot<-plot+geom_segment(aes(x=chromStart+window/2,xend=chromStart+window/2,y=2.1,yend=2.2),data=signLevel3,color="gray25")
		plot<-plot+geom_vline(aes(xintercept=-step,linetype="p-value (< 0.001)"),data=signLevel3)  ##  dummy
		line.names<-c(line.names,"p-value (< 0.001)")
		line.values<-c(line.values,"solid")
		line.colors<-c(line.colors,"gray25")
		line.shapes<-c(line.shapes,NA)
		line.sizes<-c(line.sizes,1)
	}
	if (dim(signLevel2)[1]!=0) {
		plot<-plot+geom_segment(aes(x=chromStart+window/2,xend=chromStart+window/2,y=2.1,yend=2.2),data=signLevel2,color="gray50")
		plot<-plot+geom_vline(aes(xintercept=-step,linetype="p-value (0.01-0.001)"),data=signLevel2)  ##  dummy
		line.names<-c(line.names,"p-value (0.01-0.001)")
		line.values<-c(line.values,"solid")
		line.colors<-c(line.colors,"gray50")
		line.shapes<-c(line.shapes,NA)
		line.sizes<-c(line.sizes,1)
	}

	plot<-plot+geom_point(aes(x=chromStart+window/2,y=ratio.x,color=paste0(firstName," raw")), size=.4)  ##  Raw data
	color.names<-c(color.names,paste0(firstName," raw"))  ##  Raw data legend values
	color.values<-c(color.values,"gray50")
	color.lines<-c(color.lines,"blank")
	color.shapes<-c(color.shapes,20)
	color.sizes<-c(color.sizes,4)
	color.fills<-c(color.fills,NA)

	plot<-plot+geom_point(aes(x=chromStart+window/2,y=ratio.y,color=paste0(secondName," raw")), size=.4)  ##  Raw data
	color.names<-c(color.names,paste0(secondName," raw"))  ##  Raw data legend values
	color.values<-c(color.values,"royalblue1")
	color.lines<-c(color.lines,"blank")
	color.shapes<-c(color.shapes,20)
	color.sizes<-c(color.sizes,4)
	color.fills<-c(color.fills,NA)

	plot<-plot+geom_line(aes(x=chromStart+window/2,y=fit.x,group=group.x,color=paste0(firstName," smooth")),size=.8)  ##  Smoothed data
	color.names<-c(color.names,paste0(firstName," smooth"))  ##  Smoothed data legend values
	color.values<-c(color.values,"gray25")
	color.lines<-c(color.lines,"solid")
	color.shapes<-c(color.shapes,NA)
	color.sizes<-c(color.sizes,1.5)
	color.fills<-c(color.fills,NA)

	plot<-plot+geom_line(aes(x=chromStart+window/2,y=fit.y,group=group.y,color=paste0(secondName," smooth")),size=.8)  ##  Smoothed data
	color.names<-c(color.names,paste0(secondName," smooth"))  ##  Smoothed data legend values
	color.values<-c(color.values,"blue")
	color.lines<-c(color.lines,"solid")
	color.shapes<-c(color.shapes,NA)
	color.sizes<-c(color.sizes,1.5)
	color.fills<-c(color.fills,NA)

	my.colors<-color.values  ##  Preparing named vector for colors
	names(my.colors)<-color.names
	my.lines<-line.values
	names(my.lines)<-line.names
	##  Supplying the manual color, fill and linetype scales - this is required for the legend
	plot<-plot+scale_colour_manual(name="a",values=my.colors,guide = guide_legend(override.aes = list(
		linetype=color.lines,shape= color.shapes,size=color.sizes,fill=color.fills)),breaks=names(my.colors))
	plot<-plot+scale_fill_manual(name="b",values=my.fills,guide = guide_legend(override.aes = list(
		size=fill.sizes,shape=fill.shapes,linetype=fill.lines,color=fill.colors)),breaks = names(my.fills))
	plot<-plot+scale_linetype_manual(name="c",values=my.lines,guide = guide_legend(override.aes = list(
		color=line.colors,size=line.sizes,shape=line.shapes)),breaks = names(my.lines))

	plot<-plot+plot.theme
	return(plot)
}

plotTimeCourse <- function(ratioDFs,tCourseName,tCourse,genome,circles=NULL,pointers=NULL,rectangles=NULL) {
#~ tCourse: character vector containing sample names to plot
#~ genome: dataframe containing character "chrom", integer "chromEnd" and optional integer "cen"
	require(ggplot2)
	genomeName <- deparse(substitute(genome))
	if (!is.null(circles)) circlesName <- deparse(substitute(circles))

		# reorder based on means
	tCourseMeans <- vector(mode="numeric",length=0)
	for (sample in tCourse) {
		sample <- ratioDFs[[sample]]
		tCourseMeans <- c(tCourseMeans,as.numeric(summary(sample$fit)["Mean"]))
			##  Reorder the data based on the genome dataframe
		if (genomeName == "sacCer3") {
			sample$chrom <- factor(sample$chrom,levels=unique(genome$chrom))
			sample<-sample[order(sample$chrom,sample$chromStart),]
			rownames(sample) <- 1:nrow(sample)
		}
	}
	sortIndex <- sort(tCourseMeans,decreasing=T,index.return=T)$ix
	tamedTcourse <- tCourse[sortIndex]
	window <- sample$chromEnd[1]-sample$chromStart[1]
##  Make nice x-axis labels
	tmp <- max(genome$chromEnd)
	n <- round(log10(tmp/10),digits=0)
	step <- as.integer(round(tmp/10,digits=-n))
	if (tmp/step < 10) { step <- as.integer(step/2) }
	ticks <- seq(from = 0, to = tmp, by = step)
	labels <- vector(mode="character",length=0)
	labels<-ifelse(log10(ticks)>=3 & log10(ticks)<6, paste0(ticks/1000,"kb"), ifelse(log10(ticks)>=6 & log10(ticks)<9,paste0(ticks/1000000,"Mb"), ifelse(log10(ticks)>=9 & log10(ticks)<12,paste0(ticks/1000000,"Gb"),ticks)))
	labels <- data.frame(ticks,labels)

##  initialising legendary vectors here
	color.names <- vector(mode="character",length=0)
	color.values <- vector(mode="character",length=0)
	color.lines <- vector(mode="character",length=0)
	color.shapes <- vector(mode="numeric",length=0)
	color.sizes <- vector(mode="numeric",length=0)
	color.fills<-  vector(mode="character",length=0)
	fill.names <- vector(mode="character",length=0)
	fill.values <- vector(mode="character",length=0)
	fill.shapes <- vector(mode="numeric",length=0)
	fill.sizes <- vector(mode="numeric",length=0)
	fill.colors <- vector(mode="character",length=0)
	fill.lines<- vector(mode="character",length=0)
	line.names <-vector(mode="character",length=0)
	line.values <-vector(mode="character",length=0)
	line.colors <- vector(mode="character",length=0)
	line.shapes <- vector(mode="numeric",length=0)
	line.sizes <- vector(mode="numeric",length=0)
	tCourseColors <- c("gray60","purple","darkgreen","yellow","blue","red","darkcyan","orange")
##  plotting parameters
	plot.theme<-theme(
		text = element_text(family="Arial"),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		plot.title = element_text(hjust = 0.5,size=rel(1.8)),
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.x = element_line(colour="gray50",size=.5),
		axis.title = element_text(size = rel(1.6),face="bold"),
		axis.text.y = element_blank(),
		axis.text.x = element_text(size=rel(1.7)),
		legend.position="bottom",
		legend.title=element_blank(),
		legend.text = element_text(size=rel(1.2)),
		legend.key = element_blank(),
		strip.background = element_blank(),
		strip.text = element_blank(),
		plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
	)
	plot<-ggplot(sample)  # using mapping of the last sample in the timecourse
	i <- 0
	for (sampleName in tamedTcourse) {
		i <- i+1
		sampleData <- ratioDFs[[sampleName]]
		geom = paste0("geom_ribbon(aes(x=chromStart+500,ymin=1,ymax=fit,fill='",sampleName,"'),data=sampleData)")
		plot<-plot+eval(parse(text=geom))
		fill.names<-c(fill.names,sampleName)
		fill.values<-c(fill.values,tCourseColors[i])
		fill.shapes<-c(fill.shapes,22)
		fill.sizes<-c(fill.sizes,4)
		fill.colors<-c(fill.colors,NA)
		fill.lines<-c(fill.lines,"blank")
	}
	plot<-plot+scale_y_continuous(name="Relative copy number",limits = c(0.8, 2.2))  ##  y scale
	plot<-plot+annotate(geom="segment",x=0, xend=0, y=1, yend=2, colour="gray50", lwd=0.5)  ##  y scale line
	plot<-plot+scale_x_continuous(breaks=labels$ticks,labels=labels$labels,name="Chromosome coordinates",limits=c(-window,NA),expand=c(0.01,0.01))  ##  x scale
	plot<-plot+facet_grid(chrom ~ .,scales = "free", space = "free_x")  ##  Facet
	plot<-plot+geom_text(aes(x=chromEnd+5*window,y=1.5,label=chrom),data=genome,angle=270,size=5,vjust=0)  ##  Chromosome labels
	plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=1,yend=1),data=genome,linetype="dashed",size=.25,color="gray40")  ##  Lower y line
	plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=2,yend=2),data=genome,linetype="dashed",size=.25,color="gray40")  ##  Upper y line
	plot<-plot+geom_text(aes(x=-window/2,y=1,label="1"),data=genome,size=4,hjust=1)  ##  Lower y break
	plot<-plot+geom_text(aes(x=-window/2,y=2,label="2"),data=genome,size=4,hjust=1)  ##  Upper y break
	plot<-plot+geom_segment(aes(x=0,xend=genome$chromEnd,y=0.9,yend=0.9),data=genome,size=.7)  ##  Chromosome length line
	if ("cen" %in% colnames(genome)) {
		plot<-plot+geom_segment(aes(x=genome$cen,xend=genome$cen,y=1,yend=2),data=genome,color="green3")
		plot<-plot+geom_vline(aes(xintercept=-step,linetype="Centromere"),data=genome)  ##  dummy
		line.names<-c(line.names,"Centromere")
		line.values<-c(line.values,"solid")
		line.colors<-c(line.colors,"green3")
		line.shapes<-c(line.shapes,NA)
		line.sizes<-c(line.sizes,1)
	}
	if (!is.null(circles))  {  ##  If provided, plot the origins
		circlesName <- deparse(substitute(circles))
		plot <- plot+geom_point(aes(x=chromStart+window/2,y=0.9,fill=circlesName),data=circles,color="black",size=2,shape=21)
		fill.names<-c(fill.names,circlesName)  ##  Circles legend values
		fill.values<-c(fill.values,"white")
		fill.shapes<-c(fill.shapes,21)
		fill.sizes<-c(fill.sizes,4)
		fill.colors<-c(fill.colors,"black")
		fill.lines<-c(fill.lines,"blank")
	}
		plot<-plot+ggtitle(paste0(tCourseName," time course"))  ##  Plot title
		my.fills<-fill.values
		names(my.fills)<-fill.names
		my.colors<-color.values  ##  Preparing named vector for colors
		names(my.colors)<-color.names
		my.lines<-line.values
		names(my.lines)<-line.names
		##  Supplying the manual color, fill and linetype scales - this is required for the legend
		plot<-plot+scale_colour_manual(name="a",values=my.colors,guide = guide_legend(override.aes = list(
			linetype=color.lines,shape= color.shapes,size=color.sizes,fill=color.fills)),breaks=names(my.colors))
		plot<-plot+scale_fill_manual(name="b",values=my.fills,guide = guide_legend(override.aes = list(
			size=fill.sizes,shape=fill.shapes,linetype=fill.lines,color=fill.colors)),breaks = names(my.fills))
		plot<-plot+scale_linetype_manual(name="c",values=my.lines,guide = guide_legend(override.aes = list(
			color=line.colors,size=line.sizes,shape=line.shapes)),breaks = names(my.lines))
		plot<-plot+plot.theme
		return(plot)
}
