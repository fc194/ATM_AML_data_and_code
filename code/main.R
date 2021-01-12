##### methylation

library("RColorBrewer")
library("gplots")
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rna<-read.table("73CN-AML-RNA-TCGA.csv",sep=',',header=T,row.name=1)
lbl = as.numeric(rna['ANP32A',] <= median(as.numeric(rna['ANP32A',])))
lbl[lbl==0] = 2
label = matrix(ncol=2,nrow=dim(rna)[2])
label[1:dim(rna)[2],1] = seq(1, dim(rna)[2])
label[1:dim(rna)[2],2] = lbl
label = data.frame(label)
write.table(label, "anp32a_label.txt")

run_methylation <- function() {
  
  label<-read.table("anp32a_label.txt",header=T)
  data<-read.csv("73CN-AML-cn_meth-TCGA.csv",row.names=1,header=T,sep=",")
  
  low<-label[which(label[,2]==1),1]
  high<-label[which(label[,2]==2),1]
  data.low<-data[,low]
  data.high<-data[,high]
  dnaM_cmb<-cbind(data.low,data.high)
  dnaM_cmb<-dnaM_cmb[-which(is.na(rowMeans(dnaM_cmb))),]
  
  rna_pval<-matrix(nrow=dim(dnaM_cmb)[1],ncol=2)
  for (i in 1:dim(dnaM_cmb)[1]){
    test_res<-t.test(dnaM_cmb[i,1:37],dnaM_cmb[i,38:73])
    log2fc<-log2(rowMeans(dnaM_cmb[i,38:73])/rowMeans(dnaM_cmb[i,1:37]))#high/low
    rna_pval[i,1]<-test_res$p.value
    rna_pval[i,2]<-log2fc
  }
  rownames(rna_pval)<-rownames(dnaM_cmb)
  colnames(rna_pval)<-c("pval","log2fc(high/low)")
  
  sig_label<-matrix(nrow=dim(dnaM_cmb)[1],ncol=1)
  sig_label[which(rna_pval[,1]>=0.05),1]<-"UnsigDEGs"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2] > 1)),1]<-"Up"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2] < -1)),1]<-"Down"
  sig_label[intersect(intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2] >= -1)),which(rna_pval[,2] <= 1)),1]<-"Marginal sigDEGs"
  
  rna_pval_new<-cbind(rna_pval,sig_label)
  colnames(rna_pval_new)<-c("pval","log2fc(high/low)","sig_label")
  
  sig<-rbind(rna_pval_new[which(rna_pval_new[,3] == "Up"),],rna_pval_new[which(rna_pval_new[,3] == "Down"),])
  sig_value<-rbind(dnaM_cmb[which(rna_pval_new[,3] == "Up"),],dnaM_cmb[which(rna_pval_new[,3] == "Down"),])
  
  write.table(rna_pval_new,"meth_deg.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig,"meth_sig_p_0.05_fc_1.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(dnaM_cmb,"meth_dnaM_cmb.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig_value,"meth_sig_value.txt",quote=F,row.names=T,col.names=T,sep="\t")
  
}
run_methylation()

plot_scatter <- function() {
  a<-read.table("meth_deg.txt",header=T,row.names=1,sep="\t")
  ggplot(a,aes(x=log2fc.high.low.,y=-log10(a$pval),colour=sig_label))+geom_point()+scale_color_manual(breaks=c("Up","Down","Marginal sigDEGs","UnsigDEGs"),values=c("#00A14B","#2B3990","black","#ED1C24"))+labs(x="log2(Fold Change)",y="-log10(p-value)")+theme(panel.grid.major = element_blank())+theme_bw()+ggsave("meth_deg.png",width = 11, height = 8)
}
plot_scatter()

annotation <- function(){
  a<-read.table("meth_sig_p_0.05_fc_1.txt",header=T,sep="\t")
  anno<-read.csv("/meth_sig_p_0.05_fc_1_anno.txt",header=T,sep=",")
  anno_select<-anno[,c(1,12,13,22,24,25,26,28,32)]
  merge_sig<-merge(a[,c(1,4)],anno_select,by.y="IlmnID",by.x="X",sort=F)
  merge_sig_up<-merge_sig[1:31,]
  merge_sig_down<-merge_sig[32:764,]
  table(merge_sig_up$Relation_to_UCSC_CpG_Island)
  table(merge_sig_down$Relation_to_UCSC_CpG_Island)
  
  island <-matrix(c(21, 346, 10, 387),nrow = 2)
  sshore <-matrix(c(2, 94, 29, 639),nrow = 2)
  nshore <-matrix(c(7, 114, 24, 619),nrow = 2)
  sshelf <-matrix(c(0, 9, 31, 724),nrow = 2)
  nshelf <-matrix(c(0, 18, 31, 715),nrow = 2)
  other <-matrix(c(1, 152, 30, 581),nrow = 2)
  
  fisher.test(island)
  fisher.test(sshore)
  fisher.test(nshore)
  fisher.test(sshelf)
  fisher.test(nshelf)
  fisher.test(other)
  
  tss200 <-matrix(c(6, 138, 48, 1058),nrow = 2)
  tss1500 <-matrix(c(11, 219, 43, 977),nrow = 2)
  body1 <-matrix(c(16, 382, 38, 814),nrow = 2)
  utr5 <-matrix(c(9, 148, 45, 1048),nrow = 2)
  utr3 <-matrix(c(0, 40, 54, 1156),nrow = 2)
  exon <-matrix(c(7, 73, 47, 1123),nrow = 2)
  
  fisher.test(tss200)
  fisher.test(tss1500)
  fisher.test(body1)
  fisher.test(utr5)
  fisher.test(utr3)
  fisher.test(exon)
  
  write.table(merge_sig,"meth_sig_p_0.05_fc_1_anno_select.txt",quote=F,row.names=F,col.names=T,sep="\t")
}
annotation()

heatmap <- function(){
  a<-read.table("meth_sig_p_0.05_fc_1.txt",header=T,sep="\t")
  anno<-read.csv("meth_sig_p_0.05_fc_1_anno.txt",header=T,sep=",")
  anno_select<-anno[,c(1,12,13,22,24,25,26,28,32)]
  merge_sig<-merge(a[,c(1,4)],anno_select,by.y="IlmnID",by.x="X",sort=F)
  sig_value<-read.table("meth_sig_value.txt",header=T,row.names=1)
  
  merge_sig_island_value<-sig_value[which(merge_sig[,8] == "Island"),]
  merge_sig_promo_value<-sig_value[which(merge_sig[,10] == "Promoter_Associated"),]
  
  g<-heatmap.2(as.matrix(merge_sig_island_value-rowMeans(merge_sig_island_value)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",labCol=F,labRow=F)
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"meth_island_heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"meth_island_heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
  
  g<-heatmap.2(as.matrix(merge_sig_promo_value-rowMeans(merge_sig_promo_value)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",labCol=F,labRow=F)
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"meth_promo_heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"meth_promo_heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
}

heatmap()

#### mRNA
#1 low; 2 high
process_mRNA <- function(){
  label<-read.table("anp32a_label.txt",header=T)
  rna<-read.table("73CN-AML-RNA-TCGA.csv",sep=',',header=T,row.name=1)
  rna<-rna[-which(rowMeans(rna)<0.1),]
  rna<-log2(rna+1)
  
  low<-label[which(label[,2]==1),1]
  high<-label[which(label[,2]==2),1]
  rna_low<-rna[,low]
  rna_high<-rna[,high]
  
  rna_cmb<-cbind(rna_low,rna_high)
  rna_pval<-matrix(nrow=dim(rna_cmb)[1],ncol=4)
  for (i in 1:dim(rna_cmb)[1]){
    test_res<-t.test(rna_cmb[i,1:37],rna_cmb[i,38:73])
    log2fc<-log2(rowMeans(rna_cmb[i,38:73])/rowMeans(rna_cmb[i,1:37]))#high/low
    rna_pval[i,3]<-test_res$p.value
    rna_pval[i,4]<-log2fc
    rna_pval[i,1]<-rowMeans(rna_cmb[i,38:73]) #high
    rna_pval[i,2]<-rowMeans(rna_cmb[i,1:37]) #low
  }
  rownames(rna_pval)<-rownames(rna_cmb)
  padj<-p.adjust(rna_pval[,3],method="BY")
  
  rna_pval_new<-cbind(rna_pval,padj)
  colnames(rna_pval_new)<-c("Mean expression of ANP32A high", "Mean expression of ANP32A low", "pval","log2fc(high/low)","padj")
  
  sig<-rna_pval_new[which(rna_pval_new[,5]<0.1),]
  sig_label<-matrix(nrow=dim(rna_cmb)[1],ncol=1)
  sig_label[which(rna_pval_new[,5]>=0.1),1]<-"UnsigDEGs"
  sig_label[intersect (which(rna_pval_new[,5]<0.1),which(rna_pval_new[,2]>0)),1]<-"Up"
  sig_label[intersect (which(rna_pval_new[,5]<0.1),which(rna_pval_new[,2]<0)),1]<-"Down"
  sig_value<-rna_cmb[which(rna_pval_new[,5]<0.1),]
  rna_pval_new<-cbind(rna_pval_new,sig_label)
  colnames(rna_pval_new)<-c("Mean expression of ANP32A high", "Mean expression of ANP32A low","pval","log2fc(high/low)","padj","sig_label")
  
  write.table(rna_pval_new,"deg.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig,"sig_q_0.1.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(rna_cmb,"rna_cmb.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig_value,"sig_value.txt",quote=F,row.names=T,col.names=T,sep="\t")
  
  a<-read.table("kegg_plot.txt",header=T,sep="\t")
  ggplot(data = a, aes(x = Term, y = -log10(a$PValue), fill = Label)) +  geom_bar(stat = "identity")+scale_fill_manual(values=c("#00A14B","#ED1C24")) +labs(x="",y="-log10(p-value)")+theme(axis.text.x=element_text(vjust = 1, hjust = 1,size=12),plot.margin=margin(l=30,r=5,t=5,b=5,unit="pt"))+ scale_x_discrete(limits=rev(a$Term))+coord_flip()
  
  a<-read.table("deg.txt",header=T,row.names=1,sep="\t")
  ggplot(a,aes(x=log2fc.high.low.,y=-log10(a$pval),colour=sig_label))+geom_point()+scale_color_manual(breaks=c("Up","Down","UnsigDEGs"),values=c("#00A14B","black","#ED1C24"))+labs(x="log2(Fold Change)",y="-log10(p-value)")+theme(panel.grid.major = element_blank())+theme_bw()+xlim(-4,4)
  
  a<-read.table("sig_value.txt",header=T,row.names=1)
  g<-heatmap.2(as.matrix(a-rowMeans(a)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",margins=c(2,10),labCol=F,cexRow=0.1)
  
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(sig[orderid,],"heatmap_geneid_sig.txt",sep="\t",quote=F)
  
  anp32a<-rna_cmb[which(rownames(rna_cmb) == "ANP32A|8125_calculated"),]
  plot(as.numeric(anp32a[1,]),type="l")
}
process_mRNA()

#### miRNA

process_miRNA <- function(){
  #1 low; 2 high
  label<-read.table("anp32a_label.txt",header=T)
  rna<-read.table("73CN-AML-miRNA-TCGA.csv",sep=',',header=T,row.name=1)
  rna<-rna[-which(rowMeans(rna)<0.1),]
  rna<-log2(rna+1)
  
  low<-label[which(label[,2]==1),1]
  high<-label[which(label[,2]==2),1]
  rna_low<-rna[,low]
  rna_high<-rna[,high]
  
  rna_cmb<-cbind(rna_low,rna_high)
  
  a<-read.table("rna_cmb.txt",header=T,row.names=1)
  anp32a<-a[which(rownames(a) == "ANP32A|8125_calculated"),]
  
  rna_pval<-matrix(nrow=dim(rna_cmb)[1],ncol=2)
  for (i in 1:dim(rna_cmb)[1]){
    test_res<-cor.test(as.numeric(anp32a),as.numeric(rna_cmb[i,]))
    rna_pval[i,1]<-test_res$p.value
    rna_pval[i,2]<-test_res$estimate
  }
  rownames(rna_pval)<-rownames(rna_cmb)
  colnames(rna_pval)<-c("pval","corr1")
  
  sig<-rna_pval[which(rna_pval[,1]<0.05),]
  sig_value<-rna_cmb[which(rna_pval[,1]<0.05),]
  
  sig_label<-matrix(nrow=dim(rna_cmb)[1],ncol=1)
  sig_label[which(rna_pval[,1]>=0.05),1]<-"UnsigDEGs"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2]>0)),1]<-"Up"
  sig_label[intersect (which(rna_pval[,1]<0.05),which(rna_pval[,2]<0)),1]<-"Down"
  rna_pval<-cbind(rna_pval,sig_label)
  colnames(rna_pval)<-c("pval","corr1","sig_label")
  
  write.table(rna_pval,"mirna_deg.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig,"mirna_sig_p_0.05.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(rna_cmb,"mirna_rna_cmb.txt",quote=F,row.names=T,col.names=T,sep="\t")
  write.table(sig_value,"mirna_sig_value.txt",quote=F,row.names=T,col.names=T,sep="\t")
  
  a<-read.table("mirna_deg.txt",header=T,row.names=1,sep="\t")
  ggplot(a,aes(x=pval,y=corr1,colour=sig_label))+geom_point()+scale_color_manual(breaks=c("Up","Down","UnsigDEGs"),values=c("#00A14B","black","#ED1C24"))+labs(x="Correlation test (P-value)",y="Correlation coeffients")+theme(panel.grid.major = element_blank())+theme_bw()
  
  a<-read.table("mirna_sig_value.txt",header=T,row.names=1)
  g<-heatmap.2(as.matrix(a-rowMeans(a)), col = greenred, trace="none", scale="none",Colv=NULL,dendrogram="row",margins=c(2,10),labCol=F,cexRow=0.8)
  
  orderid<-rev(rownames(a)[g$rowInd])
  ordername<-colnames(a)[g$colInd]
  write.table(orderid,"mirna_heatmap_geneid.txt",sep="\t",quote=F,row.names=F,col.names=F)
  write.table(ordername,"mirna_heatmap_samplename.txt",sep="\t",quote=F,row.names=F,col.names=F)
}
process_miRNA()


#### mRNA & miRNA

process_mRNA_miRNA(){
  rna_cmb<-read.table("rna_cmb.txt",header=T,row.names=1)
  mirna_cmb<-read.table("mirna_rna_cmb.txt",header=T,row.names=1)
  relations<-read.table("all_target_only_del_p_0.4.txt",header=F,sep="\t")
  
  corr1<-matrix(nrow=dim(relations)[1],ncol=2)
  for (i in 1:dim(relations)[1]){
    if(length(grep(pattern=paste("^",relations[i,1],"\\","|",sep=""),rownames(rna_cmb))) != 0){
      mrna<-rna_cmb[grep(pattern=paste("^",relations[i,1],"\\","|",sep=""),rownames(rna_cmb)),]
      if(dim(mrna)[1]>1){
        mrna <- colMeans(mrna)
      }
      mirna<-mirna_cmb[which(rownames(mirna_cmb) == relations[i,3]),]
      test_res<-cor.test(as.numeric(mrna),as.numeric(mirna))
      corr1[i,1]<-test_res$p.value
      corr1[i,2]<-test_res$estimate
    }
  }
  relations_new <- cbind(relations,corr1)
  colnames(relations_new)<-c("gene","id","mirna","p_val","corr1")
  
  sig<-relations_new[which(relations_new[,4]<0.05),]
  sig_label<-matrix(nrow=dim(sig)[1],ncol=1)
  sig_label[which(sig[,5]>0),1]<-"Positive"
  sig_label[which(sig[,5]<0),1]<-"Negtive"
  sig_new<-cbind(sig,sig_label)
  colnames(sig_new)[6]<-c("sig_label")
  
  write.table(sig_new,"mirna_tar_sig_p_0.05_0.4.txt",quote=F,row.names=F,col.names=T,sep="\t")
}
process_mRNA_miRNA()

