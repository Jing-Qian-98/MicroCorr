argv <- commandArgs(TRUE);
in_file <- "otu_modified_L6.txt"
width <- 7
if (length(argv) > 0){
    if (length(argv) >2){
        stop('The script needs 1 file:otu_modified_L6.txt ,also can have a width defaulted 7.')
    }
    in_file <- argv[1]
    if (length(argv) ==2){width <-argv[2]}
}
df <- read.table(in_file, sep = "\t", header = T, 
                 stringsAsFactors = F, check.names = F,
                 quote = "")
df$sum <- apply(df[,2:ncol(df)], 1, sum)
df <- df[order(df$sum, decreasing = T),]
df<-subset(df,sum>0)
df$sum <- NULL
write.table(df, "temp.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)
#system("perl /home/meta/flzx/Scripts/perl/heatmap.pl temp.txt")

library(pheatmap)
data <- read.table('temp.txt', header=TRUE, check.names = F, sep = "\t", 
                  stringsAsFactors = F, quote = "")
data <- data[match(unique(data[,1]), data[,1]),]
if(nrow(data) < 50){
	t <- nrow(data)
}else{
	t <- 50
}
data <- data[1:t,]
rownames(data) <- data[,1]
data <- data[,-1]
if((ncol(data)>30 && (ncol(data)/6)+2)>width){
    width <- (ncol(data)/6)+2
}
pdf('heatmap.pdf',onefile=F,width=width)
pheatmap(data,scale="none",border_color = NA,cluster_col=F,fontsize_row =6,color = colorRampPalette(c("green", "black", "red"))(200)) 
         #cluster_row=T,
         #cellwidth =5,
         #fontsize_col =5,
         #cellheight = 5, 
         #fontsize_row =5)
dev.off()
system("inkscape -e heatmap.png heatmap.pdf -d 300")
