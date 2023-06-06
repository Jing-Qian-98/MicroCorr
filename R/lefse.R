library("tidyverse")
library(microbiomeViz)
library(ggtree)
library(tidyverse)
library(phyloseq)
vegan_otu <- function(physeq){
 OTU <- otu_table(physeq)
 if(taxa_are_rows(OTU)){
 OTU <- t(OTU)
 }
 return(as(OTU,"matrix"))
 }
 ### 添加OTU注释信息
 vegan_tax <- function(physeq){
 tax <- tax_table(physeq)

 return(as(tax,"matrix"))
 }
 # 导入otu表格

otu = read.delim("../../data/otutab.txt",row.names = 1)
head(otu)
#导入注释文件
tax = read.delim("../..//data/taxonomy.txt",row.names = 1)
head(tax)
#导入分组文件
map = read.delim("../../data/metadata.tsv",row.names = 1)
head(map)

# head(otu)
otu = as.matrix(otu)
str(otu)

tax = as.matrix(tax)
# taxa_names(tax)

ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE),
 sample_data(map) ,
 tax_table(tax)
 # phy_tree(tree)
)
# design = as.data.frame(sample_data(ps))
#合并otu表格和tax表格#---------
otu_tax = merge(otu_table,tax_table,by = "row.names",all = F)
# head(otu_tax)

 #--------添加门类标记，未注释的结果去掉#----
 otu_tax$Kingdom = paste("k",otu_tax$Kingdom,sep = "__")
 otu_tax$Kingdom[otu_tax$Kingdom == "k__NA"] =""

 otu_tax$Phylum = paste("|p",otu_tax$Phylum,sep = "__")
 otu_tax$Phylum[otu_tax$Phylum == "|p__NA"] = ""

 otu_tax$Class = paste("|c",otu_tax$Class,sep = "__")
 otu_tax$Class[otu_tax$Class == "|c__NA"] = ""

 otu_tax$Order = paste("|o",otu_tax$Order,sep = "__")
 otu_tax$Order[otu_tax$Order == "|o__NA"] = ""

 otu_tax$Family = paste("|f",otu_tax$Family,sep = "__")
 otu_tax$Family [otu_tax$Family == "|f__NA"] = ""

 otu_tax$Genus = paste("|g",otu_tax$Genus,sep = "__")
 otu_tax$Genus[ otu_tax$Genus == "|g__NA"] = ""

 otu_tax$Species = paste("|s",otu_tax$Species,sep = "__")
 otu_tax$Species[otu_tax$Species == "|s__NA"] = ""

 otu_tax$Row.names = paste("|t",otu_tax$Row.names,sep = "__")
 #合并得到结合全部门类的OTU名称#----
 OTU_name = paste(otu_tax$Kingdom,otu_tax$Phylum,otu_tax$Class,otu_tax$Order,otu_tax$Family,
 otu_tax$Genus,sep = "")
  #----不同等级注释修改#----
 otu_tax$Phylum=paste(otu_tax$Kingdom,otu_tax$Phylum,sep = "")
 otu_tax$Class=paste(otu_tax$Phylum,otu_tax$Class,sep = "")
 otu_tax$Order=paste(otu_tax$Class,otu_tax$Order,sep = "")
 otu_tax$Family=paste(otu_tax$Order,otu_tax$Family,sep = "")
 otu_tax$Genus=paste(otu_tax$Family,otu_tax$Genus,sep = "")
 otu_tax$Species=paste(otu_tax$Genus,otu_tax$Species,sep = "")


 #替换两个括号等特殊符号#--------
 library("tidyverse")
 OTU_name = str_replace(OTU_name,"[(]","")
 OTU_name = str_replace(OTU_name,"[)]","")
 # as.character(OTU_name )
 # OTU_name = gsub("(","",OTU_name[311])
 # row.names(otu_tax) = OTU_name
 # #将otu表格和tax文件分离#-----
 # otu_table = otu_tax[2:(ncol(otu_table)+1)]
 # tax_table = otu_tax[(ncol(otu_table)+2):(ncol(otu_table)+2+6)]
 #
  #-------将不同分类登记的也添加上去
 HA = otu_table
 # 按Kingdom合并
 grp <- otu_tax[rownames(otu_tax), "Kingdom", drop=F]
 merge=cbind(HA, grp)
 HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
 colnames(HA_Kingdom)[1]="Class"

 # 按Phylum合并
 grp <- otu_tax[rownames(otu_tax), "Phylum", drop=F]
 merge=cbind(HA, grp)
 HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
 colnames(HA_Phylum)[1]="Class"

 # 按Class合并
 grp <- otu_tax[rownames(otu_tax), "Class", drop=F]
 merge=cbind(HA, grp)
 HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
 colnames(HA_Class)[1]="Class"

 # 按Order合并
 grp <- otu_tax[rownames(otu_tax), "Order", drop=F]
 merge=cbind(HA, grp)
 HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
 colnames(HA_Order)[1]="Class"

 # 按Family合并
 grp <- otu_tax[rownames(otu_tax), "Family", drop=F]
 merge=cbind(HA, grp)
 HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
 colnames(HA_Family)[1]="Class"

 # 按Genus合并
 grp <- otu_tax[rownames(otu_tax), "Genus", drop=F]
 merge=cbind(HA, grp)
 HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
 colnames(HA_Genus)[1]="Class"
 # colnames(otu_tax)
 # # 按Species合并
 # grp <- otu_tax[rownames(otu_tax), "Species", drop=F]
 # merge=cbind(HA, grp)
 # HA_Species = merge %>% group_by(Species) %>% summarise_all(sum)
 # colnames(HA_Species)[1]="Class"
 ## OTU水平
 # merge=cbind(HA, OTU_name)
 # head(otu_table)
 # HA_OTU = otu_table

 # HA_OTU = data.frame(Class = row.names(otu_table),otu_table)
 # colnames(HA_OTU )
 # head(HA_OTU)
 # 合并6个分类级
 all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus)
 dim(all)


# head(OTU_name)
# 去除重复
all = distinct(all, Class, .keep_all = TRUE)
dim(all)

#-------------lefse分析构建
all1 = all
row.names(all1) = all1$Class
all1$Class = NULL
all1 = as.matrix(all1)
# head(all1)

#-构建phylose对象
ps_G_graphlan = phyloseq(otu_table(as.matrix(all1),taxa_are_rows = TRUE),
 sample_data(ps))
ps_G_graphlan

#----提取OTU表格
otu_table = as.data.frame((vegan_otu(ps_G_graphlan)))
otu_table[otu_table==0] <- 1

# row.names(otu_table)

# head(design)

design = as.data.frame(sample_data(ps_G_graphlan))
taxtree = resTable[clapvalues<=p.lvl & ldamean$LDAscore>=lda.lvl,]

delimiter = "\\|"
tax_split <- strsplit(row.names(taxtree), delimiter)
row.names(taxtree)<- vapply(tax_split, tail, n = 1, "")

# head(taxtree)

#-提取所需要的颜色
colour = c('darkgreen','red',"blue")

selececol = colour[1:length(levels(as.factor(taxtree$class)))]

names(selececol) = levels(as.factor(taxtree$class))
A = rep("a",length(row.names(taxtree)))

i = 1
for (i in 1:length(row.names(taxtree))) {
 A[i] = selececol [taxtree$class[i]]
}

# A

lefse_lists = data.frame(node=row.names(taxtree),
 color=A,
 stringsAsFactors = FALSE
)

# str(all)

## 计算均值用于呈现结点大小
dat <- data.frame(V1=all[,1], V2=rowMeans(all[,-1]), stringsAsFactors = FALSE)

# head(dat)
# dim(dat)
# write.csv(dat,"./tree_tax.csv")
# dat$V2 = NULL
# 用物种和丰度生成树骨架
tr <- parseMetaphlanTSV(dat, node.size.offset=2, node.size.scale=0.8)

# tr

# 构造树#------

p =tree.backbone(tr, size=1,layout= 'circular')
p
# ?tree.backbone

# 注释树
p <- clade.anno(p, lefse_lists, alpha=0.3)
p

ggsave("./cs.pdf",p,width = 10,height = 10)
