#find RNA editing istes from NT sample only
#setwd("/pico/home/userexternal/whwang00/tcga_rnaediting/Rfiles")
#source("/pico/home/userexternal/whwang00/tcga_rnaediting/Rfiles/list14abcd2txt.R")
#source("/pico/home/userexternal/whwang00/tcga_rnaediting/Rfiles/list15abcd2txt.R")
#source("/pico/home/userexternal/whwang00/tcga_rnaediting/Rfiles/list7abcd2txt.R")
#source("/pico/home/userexternal/whwang00/tcga_rnaediting/Rfiles/list7abcd2annovartxt.R")
# to calculate kai square test pval and odd ratio
#using ref alt count in each sample
#put number of mapped reads in the other sample for the position not in the sample

library("VariantAnnotation")
library(data.table)
library(gplots)
library(VennDiagram)
library(parallel)
require(bit64)
#Rscript rnaediting_star_tophat_unoin_pico_computeaz.R $barcode1 $numps
args <- commandArgs(TRUE)
print(args)
rnacntfile=args[1]
dnacntfile=args[2]
outfile=args[3]
#dnavarfile=args[3]
#barcode="TCGA-A7-A0D9"
#numprocs=4
#print(paste("args: ",barcode))
cat(sprintf("RNA count   files: %s\n", rnacntfile))
cat(sprintf("DNA count   files: %s\n", dnacntfile))


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

#rnacntfile="/data/tempa0d9/TCGA-A7-A0D9_NT_star_tophat_merged_grp_dedup_karyotic_bq20_dp5_Arefbase.vcf"
#dnacntfile="/data/tempa0d9/TCGA-A7-A0D9_NT_bq25dp4_dnaaref_readcnt.vcf"
#outfile="/data/tempa0d9/TCGA-A7-A0D9_NT_bq25dp4_dnaaref_allaz.vcf"
#reading readcnt file from nt exome for A/A bq25 dp4 and compute az score
#ntexfiles=list.files(exome_nt_subdir)
#exnt_alldp=fread(paste(exome_nt_subdir,"/",excntfile_nt,sep=""), select=c("chrom","position"),header=TRUE,sep="\t")
cat(sprintf("Reading DNA cnt file: %s\n", dnacntfile))
exnt_alldp=fread(file=dnacntfile,header=FALSE, sep="\t", skip=1, select=c(1,2,3,6))
cat(sprintf("Reading DNA cnt file done\n"))
exnt_alldp=as.data.frame(exnt_alldp)
#exnt_alldp=exnt_alldp[,-7]
colnames(exnt_alldp)=c("chr","pos",   "coverage", "gcnt")
#exnt_alldp$position=paste(exnt_alldp$chrom,exnt_alldp$position,sep=":")
exnt_alldp=exnt_alldp[exnt_alldp[,3]>6,]
#real zero nonref G filter
nonrefcut=5.0
#temp_alldp=exnt_alldp[1:100000,]
workfunc=function(n1,n2) { return(conf0(n1,n2)) }
#tempalldp_sub=temp_alldp[,3:4]
#numprocs=4
#numprocs=as.integer(numprocs)
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
stUth_alldp=stnt_alldp[order(stnt_alldp$chrom,stnt_alldp$position),]
rm(stnt_alldp)
#stUth_sb_filterpostp1=stUth_sb_filterpostp1[order(stUth_sb_filterpostp1$chrom,stUth_sb_filterpostp1$position),]

#exchrpos=paste(exnt_alldp$chr,exnt_alldp$pos,sep=":")
#stuthpos=paste(stUth_alldp[,1],stUth_alldp[,2],sep=":")
cat(sprintf("Merging DNA and RNA cnt by their position\n"))
exnt_alldp=merge(exnt_alldp,stUth_alldp[,c("chrom","position")],by.x=c("chr","pos"),by.y=c("chrom","position"))
exnt_alldp=exnt_alldp[order(exnt_alldp$chr,exnt_alldp$pos),]
rm(stUth_alldp)
gc()

cat(sprintf("Computing confidence zero score\n"))
pmt=proc.time()
#exnt_alldp$zerononref=mcmapply(workfunc,exnt_alldp[,4],exnt_alldp[,3],mc.cores=numprocs)
exnt_alldp$zerononref=vapply(1:NROW(exnt_alldp), function(i) conf0(exnt_alldp$gcnt[i],exnt_alldp$coverage[i]), FUN.VALUE=double(1))
#temp_alldp$zerononref=mcmapply(workfunc,temp_alldp[,4],temp_alldp[,3],mc.cores=numprocs)
pmt1=proc.time()-pmt
print(pmt1)
  
#exnt_zconf5=exnt_alldp[exnt_alldp$zerononref>=nonrefcut,]
#rm(exnt_alldp)
gc()


write.table(exnt_alldp,file=outfile,row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)








