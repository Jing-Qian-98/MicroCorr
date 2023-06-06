library(plotrix)
library(scales)
library(RColorBrewer)

flower_plot2 <- function(sample, value, start, a, b,  
                    ellipse_col = rgb(135, 206, 235, 150, max = 255), 
                    circle_col = rgb(0, 162, 214, max = 255),
                    circle_text_cex = 1, labels=labels) {
par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
plot(c(0,10),c(0,10),type="n")
n   <- length(sample)
deg <- 360 / n
res <- lapply(1:n, function(t){
    ellipse_col <- ellipse_col[t]
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                          y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                          col = ellipse_col,
                          border = ellipse_col,
                          a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t]
    )

    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) - start,
             adj = 1,
             cex = circle_text_cex
        )

    } else {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) + start,
             adj = 0,
             cex = circle_text_cex
        )
    }           
})
plotrix::draw.circle(x = 5, y = 5, r = 1, col = circle_col, border = circle_col )

# tune location by x and y.
text(x = 5, y = 5, labels=labels)
}

argv <- commandArgs(TRUE)
system(paste("sed -i '1d' ", argv[1], seq=""))
system(paste("sed -i 's/#//g' ",argv[1],seq=""))
df <- read.table(argv[1],header=T,row.names=1,check.names=F,sep="\t")
print(str(df))

#df <- read.table("otu_modified_group.xls",header=T,row.names=1,check.names=F,sep="\t")        #for test
if(any(colnames(df)=="taxonomy")){df$taxonomy=NULL}

df$sum <- apply(df, 1, sum)
df <- df[df$sum>0,]

ll <- list()
cnt <- c()

for( i in 1:(ncol(df)-1)){
    ll[[i]] <- rownames(df)[df[,i]==df$sum]
    cnt <- c(cnt, length(ll[[i]]))
    names(cnt)[i] <- colnames(df)[i]
    ll[[i]] <- c(colnames(df)[i], ll[[i]])
}

core <- length(rownames(df)[apply(df, 1, min)>0])
ll[[length(ll)+1]] <- c("core",rownames(df)[apply(df, 1, min)>0])

lapply(ll, write, "temp_summary.txt", append=TRUE, ncolumns=1000000, sep="\t")
system("mv temp_summary.txt flower_summary.txt")

col_source=rep(c(brewer.pal(8,"Accent"),brewer.pal(12,"Set3"),brewer.pal(8,"Set2")),t=100)
col=c(col_source[2:(length(cnt)+1)])
  
pdf("flower_plot.pdf",bg="white")
flower_plot2 (names(cnt), cnt, 90, 0.5, 2, labels=core,
        ellipse_col = alpha(col, alpha = 0.5), 
        circle_col = alpha(col_source[1], alpha = 0.8))
dev.off()
