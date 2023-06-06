library(ggplot2)
library(ggsignif)
library(PMCMRplus)
library(reshape2)
library(RColorBrewer)


argv <- commandArgs(TRUE)
var <- read.table(argv[2], header=T, sep="\t", check.names=F,row.names=1)#variable.txt
rown <- colnames(var)
map <- read.table(argv[3], header=F, sep="\t",check.names=F)#map.txt
df <- var
fac <- as.factor(map[, 4])
mer <- data.frame(map[, 4], df)
df$Group<-as.factor(map[match(rownames(df),map$V1),4])

extract <- function(mer){
    out <- {}
    for ( k in 2 : ncol(mer)){
        #df1 <- wilcox.test(as.numeric(format(mer[, k], scientific=TRUE))~mer[, 1],data=mer,exact=FALSE,correct=FALSE,subset=mer[, 1] %in% c(i, j))  
        if(argv[1] == "wilcox"){
	    out <- c(out,wilcox.test(as.numeric(format(mer[, k], scientific=TRUE))~mer[, 1],data=mer,exact=FALSE,correct=FALSE,subset=mer[, 1] %in% c(i, j))$p.value)
	}else if(argv[1] == "welch-t"){
	    out <- c(out,t.test(as.numeric(format(mer[, k], scientific=TRUE))~mer[, 1],data=mer,exact=FALSE,correct=FALSE,subset=mer[, 1] %in% c(i, j))$p.value)
	}else{
	    out <- NULL
	}
        out[is.na(out)]=1
    }
    return(out)
}

dir.create(argv[1])

if(argv[1] == "wilcox" | argv[1] == "welch-t"){
res <- NULL
for ( i in c(levels(fac))){
    for ( j in c(levels(fac))){
        if ( i > j ){
            mi2 <- c(NULL)
            mj2 <- c(NULL)
            p <- extract(mer)
            p_adj <- p.adjust(p,method="fdr",length(p))
            meri <- subset(mer, mer[1]==i)
            merj <- subset(mer, mer[1]==j)
            for ( k in 2 :ncol(mer)){
                mi1 <- mean(meri[, k])
                mi2 <- c(mi2, mi1)
                mj1 <- mean(merj[, k])
                mj2 <- c(mj2, mj1)
            }
            rn <- data.frame(lapply(as.data.frame(rown), as.character), stringsAsFactors=FALSE)
            wt <- cbind(rn, mi2, mj2, p, p_adj)
	    res2 <- wt[,c(1,4)]
	    res2$compare <- paste(i,"-",j,sep="")
	    res <- rbind(res,res2)
            label <- c("Variable", paste("mean", i, sep="-"),  paste("mean", j, sep="-"), "p", "p.adj")
            wto <- rbind(label, wt)
            write.table(wto, file=paste(argv[1],"/test_", i,"_vs_", j, ".txt", sep=""), sep="\t", row.names=F, col.names=F,quote=F)
        }
    }
}
}else if(argv[1] == "anova"){
dir.create("anova/paired_compairisons")
res <- NULL
for(i in (1:(ncol(df)-1))){
    sink(paste("anova/paired_compairisons/",colnames(df)[i],".txt",sep=""))
    print(summary(aov(df[,i] ~ Group,df)))
    writeLines("")
    print(TukeyHSD(aov(df[,i] ~ Group,df)))
    sink()
    res2 <- as.data.frame(TukeyHSD(aov(df[,i] ~ Group,df))$Group)
    res2$rown <- colnames(df)[i]
    res2$compare <- rownames(res2)
    colnames(res2)[4] <- "p"
    res <- rbind(res,res2)
}

wto <- c("ID", "F value", "Pr(>F)")
for(i in (1:(ncol(df)-1))){
    stat <- c(colnames(df)[i], summary(aov(df[,i] ~ Group,df))[[1]][[4]][1], summary(aov(df[,i] ~ Group,df))[[1]][[5]][1])
    wto <- rbind(wto, stat)
}

write.table(wto, "anova/anova_results.xls", sep = "\t", quote = F, col.names = F, row.names = F)

}else if(argv[1] == "kruskal"){

dir.create("kruskal/paired_compairisons")

res <- NULL

for(i in (1:(ncol(df)-1))){
    out <- kwAllPairsDunnTest(df[,i], df$Group)
    sink(paste("kruskal/paired_compairisons/",colnames(df)[i],".txt",sep=""))
    print(kruskal.test(df[,i] ~ Group,df))
    writeLines("")
    print(summary(out))
    sink()
    res2 <- NULL
    for(m in 1:nrow(out$p.value)){
        for(n in 1:ncol(out$p.value)){
            if(m >= n){
                res2 <- rbind(res2,data.frame(paste(rownames(out$p.value)[m],"-",colnames(out$p.value)[n],sep=""),out$p.value[m,n]))
            }
        }
    }
    res2$rown <- colnames(df)[i]
    colnames(res2)[1:2] <- c("compare","p")
    res <- rbind(res,res2)
}

wto <- c("ID", "chi-squared", "p.value")
for(i in (1:(ncol(df)-1))){
    stat <- c(colnames(df)[i], kruskal.test(df[,i] ~ Group,df)$statistic, kruskal.test(df[,i] ~ Group,df)$p.value)
    wto <- rbind(wto, stat)
}
write.table(wto, "kruskal/kruskal_results.xls", sep = "\t", quote = F, col.names = F, row.names = F)
}


#boxplot
###############################################################################################################
plab <- data.frame(x=(1+length(levels(df$Group)))/2, y=Inf, lab="", Group=NA,
             variable=unique(res$rown))
class(plab$lab) <- "character"
#if(argv[1] == "wilcox" | argv[1] == "welch-t"){nc <- 4}else{nc <- 3}

for (i in 1:nrow(plab)){
    if(argv[1] == "anova" | argv[1] == "kruskal"){
        plab$lab[i] <- paste("p = ", round(as.numeric(wto[,3][wto[,1]==plab$variable[i]]),3), sep="")
    }else{
        plab$lab[i] <- " "
    }
}


df_plot <- melt(df,id="Group")
p <- list()
for(i in levels(df$Group)){
    for(j in levels(df$Group)){
        if(i > j && (paste(i,j,sep="-") %in% res$compare | paste(j,i,sep="-") %in% res$compare)){
            p[[paste(i,j,sep="-")]] <- c(i,j)
        }
    } 
}


num_x_text=max(nchar(levels(df$Group)))
if (num_x_text > 5) {angle=45;hjust=1;vjust=1}else{angle=0;hjust=0.5;vjust=1}

col_source=rep(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"),brewer.pal(8,"Set2")),t=100)
col=c(col_source[1:length(levels(df$Group))])



myplot <- ggplot(df_plot, aes(x=Group,y = value, group=Group, color=Group, fill=Group)) + 
    geom_boxplot(alpha=0.5) + geom_jitter(aes(color=Group),shape=19,size=1.5) +
    facet_wrap(.~variable, scales = "free") + 
    geom_signif(comparisons = p, step_increase = 0.12, annotations = 0.01, size = .6, textsize = 5, 
        vjust = .4, tip_length = 0.02, color="black") + 
    geom_text(aes(x, y, label=lab), color="black" ,
        data=plab, vjust=1, size = 4) +
    scale_color_manual(values = col,guide = guide_legend(nrow =min(8,nrow(df)), title = "Groups"))+
    scale_fill_manual(values = col,guide = "none")+ 
    scale_y_continuous(expand =expand_scale(mult=c(0.05,0.2)))+
    theme_bw() + 
    theme(axis.text.y = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 12,angle=angle,hjust=hjust,vjust=vjust),
	axis.title = element_blank(),axis.ticks.x = element_blank(),
        legend.justification = c(0.5,0.5),legend.position="right",
        legend.text=element_text(face="plain",size=12),legend.title = element_text(size=12,face = "bold"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.key.width=unit(0.8,'cm'),legend.key.height=unit(0.8,'cm'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(face = "bold"))
  
 
myplot2 <- ggplot_build(myplot)


res$panel <- 0
res$compare2 <- "name"
for (j in 1:nrow(res)){
    res$panel[j] <- c(1:length(unique(res$rown)))[res$rown[j] == unique(res$rown)]
    res$compare2[j] <- paste(res$compare[j], c(1:length(names(p)))[res$compare[j] == names(p)], sep="-")
}

for (i in 1:nrow(myplot2$data[[3]])){
    myplot2$data[[3]]$annotation[i] <- res$p[myplot2$data[[3]]$group[i] == res$compare2 & myplot2$data[[3]]$PANEL[i] == res$panel]
}

myplot2$data[[3]] <- myplot2$data[[3]][myplot2$data[[3]]$annotation < 0.05, ]

f1 <- myplot2$data[[3]]$annotation < 0.001
f2 <- myplot2$data[[3]]$annotation < 0.01 & myplot2$data[[3]]$annotation >= 0.001
f3 <- myplot2$data[[3]]$annotation < 0.05 & myplot2$data[[3]]$annotation >= 0.01
myplot2$data[[3]]$annotation[f1] <- "***"
myplot2$data[[3]]$annotation[f2] <- "**"
myplot2$data[[3]]$annotation[f3] <- "*"

myplot3 <- ggplot_gtable(myplot2)

pdf(paste(argv[1],"/boxplot.pdf",sep=""), 
    height = wrap_dims(length(unique(res$rown)))[1]*3, 
    width = wrap_dims(length(unique(res$rown)))[2]*length(levels(df$Group))*0.5+2.5)

plot(myplot3)

dev.off()
