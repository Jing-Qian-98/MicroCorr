library(permute)
library(vegan)
library(ggplot2)
map_file <- "map.txt"
in_file <- "otu_table_L6.txt"
output_filepath <- ""
argv <- commandArgs(TRUE)
if(length(argv) > 0){
    if (length(argv) != 2) stop ('The script needs two files: map.txt otu_modified_L6.txt')
    map_file <- argv[2]
    in_file <- argv[1]
}
df <- read.table(in_file, sep = "\t", header = T, row.names = 1, check.names = F)
Proportion_Explained <- df[nrow(df),]
map <- read.table(map_file, sep = "\t", header = T, comment.char = "@", stringsAsFactors = F, colClasses = "character")
df <- df[1:nrow(map),]
#df <- df[1:nrow(map),]
de <- c(-2,-3,-(ncol(map)))
map <- map[,de]
#map <- map[order(map[,1]),]
group<-map[,1]
#map <- map[, -1, drop = F]
#for (i in 1:ncol(map)){
    #dft <- df[map[,i] != NA,]
    #result <- rda(dft)
    #pca <- summary(result)$sites
    #Proportion_Explained <- summary(result)$cont$importance[2, 1:(ncol(pca))]
    xl <- paste("PC1 (" ,as.character(round(Proportion_Explained[1], 2)), "%)", 
            sep = "")
    yl <- paste("PC2 (" ,as.character(round(Proportion_Explained[2], 2)), "%)", 
            sep = "")
    pca <- as.data.frame(df)
    pca$Group <- map[,2];
    pca$names<-group
    library(ellipse)
    pca1 <- pca[-nrow(pca),]
    centroids <- aggregate(cbind(pca[,1],pca[,2])~Group,pca,mean)
    conf.rgn  <- do.call(rbind,lapply(unique(pca$Group),function(t)
    data.frame(Group=as.character(t),
    ellipse(cov(pca[pca$Group==t,1:2]),
        centre=as.matrix(centroids[centroids$Group==t,2:3]),
        level=0.95),
        stringsAsFactors=FALSE)))
#    q=ggplot(pca, aes(pca[,1], pca[,2], col= Group,shape=Group))+theme(text=element_text(family="mono"))+geom_point(size=3)+theme(legend.title=element_blank())+ theme(legend.key = element_rect(linetype='blank',fill = 'white'))+labs(title = "PCoA - PC1 vs PC2", x = xl, y = yl)
    colnames(pca)[1:2] <- c("PC1","PC2")
    colnames(conf.rgn) <- c("Group","PC1","PC2")
    q <- ggplot(pca, aes(PC1, PC2,colour = Group, 
        shape = Group)) +
        geom_point(size = 5) +theme(legend.text=element_text(size=12))+ geom_polygon(data=conf.rgn, aes(fill=Group,colour=NA), alpha=0.2) +
        geom_hline(data = pca, yintercept = 0, linetype = "dashed", colour = "black", size = 0.25) +
        geom_vline(data = pca, xintercept = 0, linetype = "dashed", colour = "black", size = 0.25) +
        theme_bw() + 
        theme(legend.key = element_rect(linetype='blank',fill = 'white'))+
        labs(x = xl, y = yl)+ggtitle("PCoA")+xlab(xl)+ylab(yl)
    if(length(unique(pca$Group)) <25){
    para <- 1               
    }else{
    para <- 2                
    }
    q=q + guides(col = guide_legend(ncol = para))
    if (length(unique(pca$Group)) <11){q=q+scale_colour_manual(values = c("red","blue","darkgreen","blueviolet","brown","cadetblue","chocolate1","slateblue","cyan","green","black"))+
      scale_fill_manual(values = c("red","blue","darkgreen","blueviolet","brown","cadetblue","chocolate1","slateblue","cyan","green","black"))
    }
    if (length(unique(pca$Group)) > 6){
        q <- q + scale_shape_manual(values = c(0:(length(unique(pca$Group)) - 1)))
    }
    q=q+theme_bw()+theme(text=element_text(family="mono"))+theme(legend.title=element_blank())+ 
        geom_hline(data = pca, yintercept = 0, linetype = "dashed", colour = "black", size = 0.25) +
        geom_vline(data = pca, xintercept = 0, linetype = "dashed", colour = "black", size = 0.25) +
        theme_bw() +         
	theme(legend.key = element_rect(linetype='blank',fill = 'white'))
    ggsave(paste("pcoa_2d.",colnames(map)[2],".pdf", sep = ""),width = 8, height = 7)
    #pca <- rbind(pca, Proportion_Explained)
    #rownames(pca)[nrow(pca)] <- "Proportion Explained"
    #pca[nrow(pca), ncol(pca)] <- NA
    #sample <- rownames(pca)
    #pca <- cbind(sample, pca)
    #pca$Group <- NULL
    #write.table(pca, paste(output_filepath, "pca.", colnames(map)[i], ".xls", sep = ""), 
    #            sep = "\t", quote = F, row.names = F)
#}
