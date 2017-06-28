#find RNA editing istes from NT sample only

# to calculate kai square test pval and odd ratio
#using ref alt count in each sample
#put number of mapped reads in the other sample for the position not in the sample

library("VariantAnnotation")
library(data.table)
library(gplots)
library(VennDiagram)
library(parallel)
require(bit64)

args <- commandArgs(TRUE)
print(args)
rnacntfile=args[1]
dnacntfile=args[2]
dnavarfile=args[3]

tprnacntfile=args[4]
tpdnacntfile=args[5]
tpdnavarfile=args[6]
outfile=args[7]

nonrefcut=5.0

cat(sprintf("RNA count   NT files: %s\n", rnacntfile))
cat(sprintf("DNA count   NT files: %s\n", dnacntfile))
cat(sprintf("DNA variant NT files: %s\n", dnavarfile))
cat(sprintf("RNA count   TP files: %s\n", tprnacntfile))
cat(sprintf("DNA count   TP files: %s\n", tpdnacntfile))
cat(sprintf("DNA variant TP files: %s\n", tpdnavarfile))


## Confidence of a real zero non-ref allele, 
## accounting for coverage and errors 
## 
## conf0 = -10*log10(upper-confidence bound of prob) 
## based on observing x successes out of n trials 
## Example: 
#> conf0(0,5); conf0(0,10); conf0(0,30) 
#[1] 3.460869 [1] 5.869316 [1] 10.22103 
# 
fn = function(th,x,n, alpha){ 
  pleft= pbinom(x,n,th) 
  return(pleft-alpha) 
} 
conf0 = function(x,n,alpha=0.05){ 
  if (x==n) return(0) 
  run = uniroot(fn,c(x/n,1),x=x,n=n, alpha=alpha) 
  return(-10*log10(run$root)) 
} 


# rnacntfile  ="/data/tempa0d9/TCGA-A7-A0D9_NT_star_tophat_merged_grp_dedup_karyotic_bq20_dp5_Arefbase.vcf"
# dnacntfile  ="/data/tempa0d9/TCGA-A7-A0D9_NT_bq25dp4_dnaaref_readcnt_allaz.vcf"
# dnavarfile  ="/data/tempa0d9/TCGA-A7-A0D9_NT_bqsr_1basesubstitute_ac1_dp5_q50_atog.vcf"
# 
# tprnacntfile="/data/tempa0d9/TCGA-A7-A0D9_TP_star_tophat_merged_grp_dedup_karyotic_bq20_dp5_Arefbase.vcf"
# tpdnacntfile="/data/tempa0d9/TCGA-A7-A0D9_TP_bq25dp4_dnaaref_readcnt_allaz.vcf"
# tpdnavarfile="/data/tempa0d9/TCGA-A7-A0D9_TP_bqsr_1basesubstitute_ac1_dp5_q50_atog.vcf"
# 
# outfile="/data/tempa0d9/TCGA-A7-A0D9_NTTP_bq25dp4_sb0.01_az5.0_fisherpval0.04_fc2.0_editedloci.vcf"

#reading variant vcf from exome
cat(sprintf("Working NT sample\n"))
cat(sprintf("Reading DNA variant file: %s\n", dnavarfile))
ntexome=readVcf(file=dnavarfile,"b37_lite")
ntexome_names=names(rowRanges(ntexome))
ntexome_names    =vapply(1:length(ntexome_names), function(i) ntexome_names[i]=strsplit(ntexome_names[i],"_")[[1]][1],FUN.VALUE=character(1))
rm(ntexome)
gc()

#real zero nonref G filter


cat(sprintf("Reading DNA cnt file: %s\n", dnacntfile))
exnt_alldp=fread(file=dnacntfile,header=FALSE, sep="\t", skip=1)
#exnt_pos=paste(exnt_alldp$chrom,exnt_alldp$position,sep=":")
exnt_alldp=as.data.frame(exnt_alldp)
colnames(exnt_alldp)=c("chr","pos","depth","gcnt","conf0")

workfunc=function(n1,n2) { return(conf0(n1,n2)) }
gc()



#Reading readcnt for all pos in each sample, merged start and tophat alignment
stnt_alldp=NULL
cat(sprintf("reading RNA count file: %s\n",rnacntfile))
stnt_alldp=fread(file=rnacntfile, select=c("chrom","position","q20_depth","Altbase","Refcnt","Altcnt","Alt_SB(pos/(pos+neg))","Alt_ER(alt/q20depth)","ref_plus","ref_minus","alt_plus","alt_minus"),header=TRUE,sep="\t")
#stnt_alldp$position=paste(stnt_alldp$chrom,stnt_alldp$position,sep=":")
stnt_alldp=as.data.frame(stnt_alldp)
#stnt_alldp=stnt_alldp[stnt_alldp$Altcnt>0,]
#stnt_alldp=subset(stnt_alldp, select=-chrom)
chrs=unique(stnt_alldp$chrom)
chrs=chrs[!grepl("GL", chrs)]
stnt_alldp=stnt_alldp[stnt_alldp$chrom %in% chrs,]

cat(sprintf("Reading RNA cnt file done\n"))
gc()
cat(sprintf("Ordering RNA cnt\n"))
stUth_alldp=stnt_alldp
stUth_alldp$position=paste(stUth_alldp$chrom,stUth_alldp$position,sep=":")
stUth_alldp=stUth_alldp[order(stUth_alldp$chrom,stUth_alldp$position),]
rm(stnt_alldp)

gc()

#remove DNA AG variants
stUth_alldp=stUth_alldp[!(stUth_alldp$position %in% ntexome_names),]

#stUth_alldp1=stUth_alldp2
#rm(stUth_alldp2)

#stUth_alldp$position=paste(stUth_alldp$chrom,stUth_alldp$position,sep=":")
#stuthntrefpos=stUth_alldp$position[stUth_alldp$Altbase=="-"] #AA
stuthntaltpos=stUth_alldp$position[stUth_alldp$Altbase=="G"] #AG
stuthnt_alt=stUth_alldp[stUth_alldp$position %in% stuthntaltpos,] # only positions with A-to-G loci
gc()
#strand bias filter fisher exact test using the cnt form both strand
stcntmat=cbind(stuthnt_alt$ref_plus,stuthnt_alt$ref_minus,stuthnt_alt$alt_plus,stuthnt_alt$alt_minus)
pmt=proc.time()
#pvals=vapply(1:NROW(cntmat), function(i) chisq.test(matrix(cntmat[i,],nrow=2,ncol=2),correct=TRUE)$p.value, FUN.VALUE=double(1))
stpvals=vapply(1:NROW(stcntmat), function(i) fisher.test(t(matrix(stcntmat[i,],nrow=2,ncol=2)))$p.value, FUN.VALUE=double(1))
pmt1=proc.time()-pmt
print(pmt1)
stuthnt_alt$sbfisherpval=stpvals
stuthnt_alt=stuthnt_alt[stuthnt_alt$sbfisherpval>=0.01,]
gc()
#write.table(stuthnt_alt,file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20_dp5_aadna_dnavar_sbfisher0.01.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)


#need to make list that has eratio(0.1) and azcut(5,10)
#need to remove AG variant a: position(1:129874), exnt_alldp(chr pos)
templist=unlist(strsplit(stuthnt_alt$position,":"))
stuthnt_alt$position=templist[seq(2,length(templist),2)]
stuthnt_alt$position=as.integer(stuthnt_alt$position)
#stuthnt_alt1=merge(stuthnt_alt,exnt_alldp[,c("V1","V2","V5")],by.x=c("chrom","position"),by.y=c("V1","V2"))
stuthnt_alt1=merge(stuthnt_alt,exnt_alldp,by.x=c("chrom","position"),by.y=c("chr","pos"))
stuthnt_alt1=stuthnt_alt1[stuthnt_alt1$conf0>=nonrefcut,]
stuthnt_alt1=stuthnt_alt1[order(stuthnt_alt1$chrom,stuthnt_alt1$position),]
rm(stuthnt_alt)
gc()

stUth_sb_filterposnt1=stuthnt_alt1
# #stUth_sb_filterposnt1=stUth_sb_filterposnt1[stUth_sb_filterposnt1$position %in% zeroconfnames_nt,]
# #stUth_sb_filterpostp1=stUth_sb_filterpostp1[stUth_sb_filterpostp1$position %in% bothzeroconf,]
# write.table(stUth_sb_filterposnt1,file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf5.0.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
# write.table(stUth_sb_filterposnt1[stUth_sb_filterposnt1[,14]>=10.0,],file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf10.0.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
# stUth_sb_filterposnt2=stUth_sb_filterposnt1[stUth_sb_filterposnt1[,14]>=10.0,]
# 
 #stUth_sb_filterposnt1$eratio=stUth_sb_filterposnt1$Altcnt/(stUth_sb_filterposnt1$Refcnt+stUth_sb_filterposnt1$Altcnt)
# stUth_sb_filterposnt2$eratio=stUth_sb_filterposnt2$Altcnt/(stUth_sb_filterposnt2$Refcnt+stUth_sb_filterposnt2$Altcnt)
# 
# write.table(stUth_sb_filterposnt1[stUth_sb_filterposnt1$eratio>=0.1,],file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf5.0_eratio0.1.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
# write.table(stUth_sb_filterposnt2[stUth_sb_filterposnt2$eratio>=0.1,],file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf10.0_eratio0.1.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 colnames(stUth_sb_filterposnt1)[5]="RNA_refcnt"
 colnames(stUth_sb_filterposnt1)[6]="RNA_altcnt"
 colnames(stUth_sb_filterposnt1)[13]="RNA_sbpval"
 colnames(stUth_sb_filterposnt1)[14]="DNA_depth"
 colnames(stUth_sb_filterposnt1)[15]="DNA_gcnt"
 colnames(stUth_sb_filterposnt1)[16]="DNA_conf0"
 
 #write.table(stUth_sb_filterposnt1[,c(1,2,3,4,5,6,13,14,15,16)],file=outfile,row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 

 cat(sprintf("Working TP sample\n"))
 cat(sprintf("Reading DNA variant file: %s\n", tpdnavarfile))
 tpexome=readVcf(file=tpdnavarfile,"b37_lite")
 tpexome_names=names(rowRanges(tpexome))
 tpexome_names    =vapply(1:length(tpexome_names), function(i) tpexome_names[i]=strsplit(tpexome_names[i],"_")[[1]][1],FUN.VALUE=character(1))
 rm(tpexome)
 gc()
 
 #real zero nonref G filter
 
 cat(sprintf("Reading DNA cnt file: %s\n", tpdnacntfile))
 extp_alldp=fread(file=tpdnacntfile,header=FALSE, sep="\t", skip=1)
 #exnt_pos=paste(exnt_alldp$chrom,exnt_alldp$position,sep=":")
 extp_alldp=as.data.frame(extp_alldp)
 colnames(extp_alldp)=c("chr","pos","depth","gcnt","conf0")
 gc()
 
 
 
 #Reading readcnt for all pos in each sample, merged start and tophat alignment
 sttp_alldp=NULL
 cat(sprintf("reading RNA count file: %s\n",tprnacntfile))
 sttp_alldp=fread(file=tprnacntfile, select=c("chrom","position","q20_depth","Altbase","Refcnt","Altcnt","Alt_SB(pos/(pos+neg))","Alt_ER(alt/q20depth)","ref_plus","ref_minus","alt_plus","alt_minus"),header=TRUE,sep="\t")
 #stnt_alldp$position=paste(stnt_alldp$chrom,stnt_alldp$position,sep=":")
 sttp_alldp=as.data.frame(sttp_alldp)
 #stnt_alldp=stnt_alldp[stnt_alldp$Altcnt>0,]
 #stnt_alldp=subset(stnt_alldp, select=-chrom)
 chrs=unique(sttp_alldp$chrom)
 chrs=chrs[!grepl("GL", chrs)]
 sttp_alldp=sttp_alldp[sttp_alldp$chrom %in% chrs,]
 
 cat(sprintf("Reading RNA cnt file done\n"))
 gc()
 cat(sprintf("Ordering RNA cnt\n"))
 tpstUth_alldp=sttp_alldp
 tpstUth_alldp$position=paste(tpstUth_alldp$chrom,tpstUth_alldp$position,sep=":")
 tpstUth_alldp=tpstUth_alldp[order(tpstUth_alldp$chrom,tpstUth_alldp$position),]
 rm(sttp_alldp)
 
 gc()
 
 #remove DNA AG variants
 tpstUth_alldp=tpstUth_alldp[!(tpstUth_alldp$position %in% tpexome_names),]
 
 #stUth_alldp1=stUth_alldp2
 #rm(stUth_alldp2)
 
 #stUth_alldp$position=paste(stUth_alldp$chrom,stUth_alldp$position,sep=":")
 #stuthntrefpos=stUth_alldp$position[stUth_alldp$Altbase=="-"] #AA
 tpstuthntaltpos=tpstUth_alldp$position[tpstUth_alldp$Altbase=="G"] #AG
 tpstuthnt_alt=tpstUth_alldp[tpstUth_alldp$position %in% tpstuthntaltpos,] # only positions with A-to-G loci
 gc()
 #strand bias filter fisher exact test using the cnt form both strand
 tpstcntmat=cbind(tpstuthnt_alt$ref_plus,tpstuthnt_alt$ref_minus,tpstuthnt_alt$alt_plus,tpstuthnt_alt$alt_minus)
 pmt=proc.time()
 #pvals=vapply(1:NROW(cntmat), function(i) chisq.test(matrix(cntmat[i,],nrow=2,ncol=2),correct=TRUE)$p.value, FUN.VALUE=double(1))
 stpvals=vapply(1:NROW(tpstcntmat), function(i) fisher.test(t(matrix(tpstcntmat[i,],nrow=2,ncol=2)))$p.value, FUN.VALUE=double(1))
 pmt1=proc.time()-pmt
 print(pmt1)
 tpstuthnt_alt$sbfisherpval=stpvals
 tpstuthnt_alt=tpstuthnt_alt[tpstuthnt_alt$sbfisherpval>=0.01,]
 
 #write.table(stuthnt_alt,file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20_dp5_aadna_dnavar_sbfisher0.01.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 
 
 #need to make list that has eratio(0.1) and azcut(5,10)
 #need to remove AG variant a: position(1:129874), exnt_alldp(chr pos)
 templist=unlist(strsplit(tpstuthnt_alt$position,":"))
 tpstuthnt_alt$position=templist[seq(2,length(templist),2)]
 tpstuthnt_alt$position=as.integer(tpstuthnt_alt$position)
 #stuthnt_alt1=merge(stuthnt_alt,exnt_alldp[,c("V1","V2","V5")],by.x=c("chrom","position"),by.y=c("V1","V2"))
 tpstuthnt_alt1=merge(tpstuthnt_alt,extp_alldp,by.x=c("chrom","position"),by.y=c("chr","pos"))
 tpstuthnt_alt1=tpstuthnt_alt1[tpstuthnt_alt1$conf0>=nonrefcut,]
 tpstuthnt_alt1=tpstuthnt_alt1[order(tpstuthnt_alt1$chrom,tpstuthnt_alt1$position),]
 rm(tpstuthnt_alt)
 gc()
 
 tpstUth_sb_filterposnt1=tpstuthnt_alt1
 # #stUth_sb_filterposnt1=stUth_sb_filterposnt1[stUth_sb_filterposnt1$position %in% zeroconfnames_nt,]
 # #stUth_sb_filterpostp1=stUth_sb_filterpostp1[stUth_sb_filterpostp1$position %in% bothzeroconf,]
 # write.table(stUth_sb_filterposnt1,file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf5.0.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 # write.table(stUth_sb_filterposnt1[stUth_sb_filterposnt1[,14]>=10.0,],file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf10.0.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 # stUth_sb_filterposnt2=stUth_sb_filterposnt1[stUth_sb_filterposnt1[,14]>=10.0,]
 # 
 #stUth_sb_filterposnt1$eratio=stUth_sb_filterposnt1$Altcnt/(stUth_sb_filterposnt1$Refcnt+stUth_sb_filterposnt1$Altcnt)
 # stUth_sb_filterposnt2$eratio=stUth_sb_filterposnt2$Altcnt/(stUth_sb_filterposnt2$Refcnt+stUth_sb_filterposnt2$Altcnt)
 # 
 # write.table(stUth_sb_filterposnt1[stUth_sb_filterposnt1$eratio>=0.1,],file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf5.0_eratio0.1.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 # write.table(stUth_sb_filterposnt2[stUth_sb_filterposnt2$eratio>=0.1,],file=paste(ntoutdir,"/",barcode,"_nt_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf10.0_eratio0.1.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 colnames(tpstUth_sb_filterposnt1)[5]="RNA_refcnt"
 colnames(tpstUth_sb_filterposnt1)[6]="RNA_altcnt"
 colnames(tpstUth_sb_filterposnt1)[13]="RNA_sbpval"
 colnames(tpstUth_sb_filterposnt1)[14]="DNA_depth"
 colnames(tpstUth_sb_filterposnt1)[15]="DNA_gcnt"
 colnames(tpstUth_sb_filterposnt1)[16]="DNA_conf0"
 
 #finish reading rna cnt, dna cnt, dna var and apply sbfilter, remove dna ag vars, filter by conf0 vscore for each NT and TP sample independently
 #need to compute fisher pval and fc of 2x2 table from NT and TP 
 
 both_alt=merge(stUth_sb_filterposnt1,tpstUth_sb_filterposnt1[,c(1,2,3,4,5,6)],by.x=c(1,2),by.y=c(1,2))
 both_alt=both_alt[order(both_alt$chrom,both_alt$position),]
 
 cnames=colnames(both_alt)
 cnames=gsub("x","nt",cnames)
 cnames=gsub("y","tp",cnames)
 colnames(both_alt)=cnames
 rm(cnames)
 
 cat(sprintf("Computing fisher exact test of NT and TP counts\n"))
 stcntmat=cbind(both_alt$RNA_refcnt.nt,both_alt$RNA_altcnt.nt,both_alt$RNA_refcnt.tp,both_alt$RNA_altcnt.tp)
 pmt=proc.time()
 #pvals=vapply(1:NROW(cntmat), function(i) chisq.test(matrix(cntmat[i,],nrow=2,ncol=2),correct=TRUE)$p.value, FUN.VALUE=double(1))
 stpvals=vapply(1:NROW(stcntmat), function(i) fisher.test(matrix(stcntmat[i,],nrow=2,ncol=2))$p.value, FUN.VALUE=double(1))
 pmt1=proc.time()-pmt
 print(pmt1)
 
 cat(sprintf("Computing Odds ratio\n"))
 pmt=proc.time()
 stcntmat1    =t(vapply(1:NROW(stcntmat), function(i) if(0%in%stcntmat[i,]) {stcntmat[i,]+1} else{stcntmat[i,]},FUN.VALUE=double(4)))
 pmt1=proc.time()-pmt
 print(pmt1)
 rm(stcntmat)
 
 pmt=proc.time()
 stfcs=vapply(1:NROW(stcntmat1), function(i) (stcntmat1[i,1]*stcntmat1[i,4])/(stcntmat1[i,2]*stcntmat1[i,3]), FUN.VALUE=double(1))
 pmt1=proc.time()-pmt
 print(pmt1)
 rm(stcntmat1)
 
 both_alt$stpvals=stpvals  
 both_alt$stfcs=log2(stfcs)
 pcut=0.05
 fccut1=log2(2.0)
 fccut2=log2(0.5)
 
 #write.table(both_alt,file=paste(nttpoutdir,"/",barcode,"_stuth_bq20dp5_dnavar_sbfisher0.01_zeroconf5.0_all.txt",sep=""),row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 write.table(both_alt,file=outfile,row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 






