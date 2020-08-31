####################################################################################
###    following code contained result visualization and parameters of methods   ###
###    of "Longitudinal multi-omics transition associated with the mortality and ###
###    survival of critically ill patients with COVID-19"                        ###
###                                                                              ###
###    code author : Shiyou Wang, Jiafeng Li                                     ###
###                                                                              ###
####################################################################################

#########################################
############    DEseq2      #############
#########################################

#### input data
mRNA_exp_filter <- read.table("/Users/wangshiyou/Documents/work/BGI/program/新冠/同济时间序列/data/Time_series_data/RNA_timeSeries/mRNA_counts_filter1.txt",header = T,sep = "\t",row.names = 1)
mRNA_exp_filter <- round(mRNA_exp_filter,0)
####  entrez id was transformed to gene symbol
symbol.change <- read.table("/Users/wangshiyou/Documents/work/BGI/program/新冠/同济时间序列/data/Time_series_data/gene2symbol.txt",header = F,sep = "\t",stringsAsFactors = F)
length(intersect(symbol.change$V1,rownames(mRNA_exp_filter)))
s1 <- symbol.change$V2
names(s1) <- symbol.change$V1
rownames(mRNA_exp_filter) <- s1[rownames(mRNA_exp_filter)]

library(DESeq2)
#### construction of coldata including cluster information
all_patient_info2 = read.table("/Users/wangshiyou/Documents/work/BGI/program/新冠/同济时间序列/0713 时序样本_编号对应_英文版.txt",header = T,sep = "\t", check.names = F,stringsAsFactors = F,comment.char = "")
condition <- all_patient_info2[which(all_patient_info2$Analysis_ID_RNA %in% colnames(mRNA_exp_filter))  ,c(3,8,9,12,31:37)]
cluster_id <- read.table("/Users/wangshiyou/Documents/work/BGI/program/新冠/同济时间序列/result/4omics_trans/mRNA_cluster3_-4.txt",row.names = 1)
cluster_id2 <- cluster_id[condition$Analysis_ID_RNA,1]
condition$cluster <- cluster_id2

mRNA_exp_filter <- mRNA_exp_filter[,condition$Analysis_ID_RNA]
#### DESeq2
dds <- DESeqDataSetFromMatrix(mRNA_exp_filter, coldata, design= ~gender+age+termimal_status2) # 构建dds矩阵
head(dds)
dds2 <- DESeq(dds) 

resultsNames(dds2) 
res4 <- results(dds2,contrast = c("termimal_status2","Critical","Severe"),pAdjustMethod = "BH",alpha = 0.05) # 使用results()函数获取结果，并赋值给res
summary(res4)

rr4 <- as.data.frame(res4)
rr4[is.na(rr4)] <- 1
rr4$symbol <- rownames(rr4)
rr4$status <- c(rep("Critical_Severe",dim(as.data.frame(res4))[1]))
rr4 <- rr4[which(rr4$padj < 0.05),]
write.table(rr4,"mRNA_diff_significant_bystatus_113.txt",quote = F,sep = "\t",row.names = F)

#########################################
############   barplot      #############
#########################################

p <- ggplot(data,aes(x = ID,y=Age,fill = Gender))+
  geom_bar(stat = "identity", alpha = 0.7) + 
  geom_hline(yintercept =c(20,40,60,75),linetype=2,size=.25)+
  geom_text(aes(x=ID,y=Age,label=Age),family="myfont",size=3.5,lineheight=1)+
  scale_x_continuous(limits=c(0,66),expand=c(0,0))+
  scale_y_continuous(limits=c(-30,80))+
  coord_polar(start=-14.17)+theme_void()

#########################################
#####   boxplot with testing      #######
#########################################

my_comparisons1 <- list(c("Severe", "Critical"))
my_comparisons2 <- list(c("Critical", "Death"))
my_comparisons3 <- list(c("Severe", "Death"))
pb1<-ggboxplot(m1,x="Terminal_status",y=markers[i],color="Terminal_status",fill=NULL,add = "jitter",
               bxp.errorbar.width = 0.6,width = 0.4,size=0.4,font.label = list(size=30),
               palette = c("#4DAF4A","#377EB8","#E41A1C"))+
  stat_compare_means(method="wilcox.test",hide.ns = F,comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3),label="p.signif")+
  theme(panel.background =element_blank(),legend.position = "NA")+labs(title = markers[i])

#########################################
############     PCA        #############
#########################################

library(stats)
data1 = read.table("mRNA.test.txt",sep = "\t",check.names = F,stringsAsFactors = F,
                   header = T,row.names = 1,comment.char = "") 
set.seed(123)
final <- kmeans(scale(forehead),3,nstart=100)

fviz_cluster(final,data = scale(forehead),labelsize = 8,shape = 1)

pca <- stats::prcomp(scale(forehead), scale = FALSE, center = FALSE)
cluster <- as.factor(final$cluster) ## save PCA and cluster result

ind <- facto_summarize(pca, element = "ind", result = "coord", 
                       axes = c(1,2))
eig <- get_eigenvalue(pca)[c(1,2), 2]
xlab = paste0("Dim", 1, " (", round(eig[1], 
                                    1), "%)")
ylab = paste0("Dim", 2, " (", round(eig[2], 
                                    1), "%)")

## use ggplot2 to beautify PCA plot
cluster1 = read.table("/Users/wangshiyou/Documents/work/BGI/program/新冠/同济时间序列/data/Time_series_data/PCA_mRNA/cluster.txt",sep = "\t", check.names = F,stringsAsFactors = F,comment.char = "")
yy = as.character(c(cluster1[,2]))
r1 <- which(yy == 3)
r2 <- which(yy == 1)
yy[r1] <- 1
yy[r2] <- 3

p <- ggplot(ind)+
  xlab(xlab) +
  ylab(ylab) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme_classic()+
  geom_point(aes(x=ind[,2], y=ind[,3], color=cluster), size=3) +
  geom_text_repel(aes(x=ind[,2], y=ind[,3],color = yy,label = rownames(ind)), size=2.5)+
  stat_ellipse(aes(x=ind[,2], y=ind[,3], fill=cluster), size=1, geom="polygon", level=0.8, alpha=0.3)# +

pdf("mRNA_PCA_final_0714.pdf",width = 14,height = 11)
print(p)
dev.off()


#########################################
############     Mfuzz      #############
#########################################

m1 <- mRNA.FPKM[diff.genes,s1]
m1 <- data.matrix(m1)
#m1 <- na.omit(m1)
m1 <- m1[-which(rowSums(m1) == 0),]
m1 <- log2(m1+0.000001)
m2 <- new("ExpressionSet",exprs = m1)
m2 <- standardise(m2)
m <- mestimate(m2)
#### elbow plot to determine clusters number 
Dmin(m2, m=m, crange=seq(2,20,1), repeats=3, visu=TRUE)
c = 13
cl <- mfuzz(m2, c = c, m = m)
## print Gene ID
write.table(cl$cluster,paste(i,"Mfuzz.txt",sep = "_"),quote=F,row.names=T,col.names=F,sep="\t")
## plot broken line 
pdf(paste(i,"Mfuzz.pdf",sep = "_"),width = 10,height = 16)
mfuzz.plot(m2,cl,mfrow=c(5,3),new.window= FALSE)
dev.off()

#########################################
###  linear regression fitted line    ###
#########################################
r2_HB <- c()
r2_RBC <- c()
coefficients_RBC <- c()
coefficients_HB <- c()
NRMSE_RBC <- c()
NRMSE_HB <- c()
for (i in seq(unique(RBC_HB_15$patient_id))) {
  erythrocyte.numbers <- RBC_HB_15[which(RBC_HB_15$patient_id == unique(RBC_HB_15$patient_id)[i]),"RBC"]
  Hb.numbers <- RBC_HB_15[which(RBC_HB_15$patient_id == unique(RBC_HB_15$patient_id)[i]),"Hb"]
  da <- data.frame(number = erythrocyte.numbers,Hb.number = Hb.numbers,time = c(1:length(erythrocyte.numbers)))
  
  model <- lm(number~time,data = da)
  model2 <- lm(Hb.number~time,data = da)
  
  a <- summary(model)
  b <- summary(model2)
  
  r2_RBC <- c(r2_RBC,a$adj.r.squared)
  r2_HB <- c(r2_HB,b$adj.r.squared)
  
  coefficients_RBC <- c(coefficients_RBC,a$coefficients[2,1])
  coefficients_HB <- c(coefficients_HB,b$coefficients[2,1])
  
  prediction <- predict(model,da)
  prediction2 <- predict(model2,da)
  
  r <- rmse(da$number, prediction)
  r_Hb <- rmse(da$Hb.number, prediction2)
  
  NRMSE <- round(r/(max(da$number)-min(da$number)),3)
  NRMSE_Hb <- round(r_Hb/(max(da$Hb.number)-min(da$Hb.number)),3)
  
  NRMSE_RBC <- c(NRMSE_RBC,NRMSE)
  NRMSE_HB <- c(NRMSE_HB,NRMSE_Hb)
  
  p1 <- ggplot(da,aes(x = time,y = number))+geom_line(colour = "brown")+
    theme(panel.background = element_blank(),axis.line = element_line(size = 1))+ggtitle(unique(RBC_HB_15$patient_id)[i])+theme(plot.title = element_text(hjust = 0.1))+
    geom_smooth(method = lm,colour="#764C29",fill="#E7E1D7")+
    stat_poly_eq(
      aes(label = paste("RBC_NRMSE:",NRMSE,"  ", ..adj.rr.label.., sep = '~')),
      formula = y ~ x,  parse = TRUE,
      size = 5, 
      label.x = 0.99,  
      label.y = 0.95)
  p2 <- ggplot(da,aes(x = time,y = Hb.number))+geom_line(colour = "blue")+
    theme(panel.background = element_blank(),axis.line = element_line(size = 1))+
    geom_smooth(method = lm,colour="#80B1D3",fill="#80B1D3")+
    stat_poly_eq(
      aes(label = paste("Hb_NRMSE:",NRMSE_Hb,"  ",..adj.rr.label.., sep = '~')),
      formula = y ~ x,  parse = TRUE,
      size = 5, #text size
      label.x = 0.99,  #text position
      label.y = 0.8)
  
  # we integrated two clinical parameters to one plot
  
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], 
                       pp$t,pp$l, pp$b, pp$l)
  #axis tweaks
  ia <- which(g2$layout$name == "axis-l")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths))
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  grid.newpage() 
  pdf(paste(unique(RBC_HB_15$patient_id)[i],".pdf",sep = ""))
  grid.draw(g)
  dev.off()
}  

#########################################
###   construction of DC network      ###
#########################################

library(MEGENA)
library(DGCA)
design_multi<-matrix(ncol = 2,nrow = ncol(Exp_multi))
row.names(design_multi)<-colnames(Exp_multi)
design_multi[which(colnames(Exp_multi)%in%Cluster2),1]=1
design_multi[which(colnames(Exp_multi)%in%Cluster1),2]=1
design_multi[which(is.na(design_multi))]=0
colnames(design_multi)<-c("Cluster2","Cluster1")
############################################################################################################################################
#################################################---------Calculate DC------------------####################################################
# input parameters
n.cores <- 4; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "Spearman" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
### 
############################################################################################################################################
design_trans_c1c2<-matrix(ncol = 2,nrow = ncol(Exp_trans_c1_c2))
row.names(design_trans_c1c2)<-colnames(Exp_trans_c1_c2)
design_trans_c1c2[which(colnames(Exp_trans_c1_c2)%in%Cluster2),1]=1
design_trans_c1c2[which(colnames(Exp_trans_c1_c2)%in%Cluster1),2]=1
design_trans_c1c2[which(is.na(design_trans_c1c2))]=0
colnames(design_trans_c1c2)<-c("cluster2","cluster1")
#----------------------------------------------------------#
design_trans_ICU<-matrix(ncol = 2,nrow = ncol(Exp_trans))
row.names(design_trans_ICU)<-colnames(Exp_trans)
design_trans_ICU[which(colnames(Exp_trans)%in%severe),1]=1
design_trans_ICU[which(colnames(Exp_trans)%in%death_critical),2]=1
design_trans_ICU[which(is.na(design_trans_ICU))]=0
colnames(design_trans_ICU)<-c("N_ICU","ICU")
#------------------------------------------------------------#
edge_multi_c1c2 = ddcorAll(inputMat = Exp_trans_c1_c2, 
                           design = design_trans_c1c2, 
                           compare = c("cluster2","cluster1"),
                           corrType = "spearman",
                           adjust = "fdr",
                           nPerms = 0)
edge_multi_ICU = ddcorAll(inputMat = Exp_trans, 
                          design = design_trans_ICU, 
                          compare = c("N_ICU","ICU"),
                          corrType = "spearman",
                          adjust = "fdr",
                          nPerms = 0)
#### Used if NA Exist
if (length(which(is.na(edge_multi_ICU[[3]])|is.na(edge_multi_ICU[[5]])))!=0) {
  edge_multi_ICU<-edge_multi_ICU[-which(is.na(edge_multi_ICU[[3]])|is.na(edge_multi_ICU[[5]])),]
}
#module_multi<-ddMEGENA(ddcor_res = edge_multi,adjusted = FALSE,evalCompactness = TRUE)

ddcor_res=edge_multi_ICU
ddcor_res_sig = ddcor_res[ddcor_res$pValDiff<0.05,] 
ddcor_res_sig = ddcor_res_sig[-which(ddcor_res_sig[,1]%in%row.names(Exp_mRNA)&ddcor_res_sig[,2]%in%row.names(Exp_mRNA)),]
ddcor_res_sig = ddcor_res_sig[which(ddcor_res_sig[,3]<0&ddcor_res_sig[,5]<0),]
ddcor_res_sig_c1c2<-ddcor_res_sig
ddcor_res_sig_ICU<-ddcor_res_sig
a<-merge(ddcor_res_sig_ICU,ddcor_res_sig_c1c2,by = c("Gene1","Gene2"))
ddcor_res_megena = ddcor_res_sig[,c("Gene1","Gene2","zScoreDiff")] 
ddcor_res_megena$zScoreDiff = abs(ddcor_res_megena$zScoreDiff) 
pfn_res = calculate.PFN(a, doPar =TRUE, num.cores = 2) 
pfn_res$weight = (pfn_res$weight/max(pfn_res$weight)) * 0.999999999
g <- graph.data.frame(pfn_res,directed = FALSE)

##### perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
                           mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 10,max.size = vcount(g)/2,
                           doPar = TRUE,num.cores = n.cores,n.perm = hub.perm,
                           save.output = FALSE)
summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                       min.size = 10,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)

pnet.obj <- plot_module(output.summary = summary.output,PFN = g,subset.module = "c1_3",
                        layout = "kamada.kawai",label.hubs.only = TRUE,
                        gene.set = NULL,color.code =  "grey",
                        output.plot = FALSE,out.dir = "modulePlot",col.names = c("magenta","green","cyan"),label.scaleFactor = 20,
                        hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)
print(pnet.obj[[1]])
colnames(pfn_res)[1:2]<-colnames(edge_multi)[1:2]
module_list<-as.data.frame(MEGENA.output[["node.summary"]][["module.membership"]])
colnames(pfn_res)[1:2]<-c("Gene1","Gene2")
final_Result_c1c2<-merge(pfn_res,ddcor_res_sig,by = colnames(pfn_res)[1:2])
save(Exp_multi,design_multi,edge_multi,pfn_res,MEGENA.output, file = "c1_c2.RData")
identity<-merge(final_Result,final_Result_ICU,by = colnames(pfn_res)[1:2])
write.csv(final_Result_c1c2,"final_result_trans_c1c2.csv")
write.csv(module_list,"module_list_multi_trans_c1c2.csv")

