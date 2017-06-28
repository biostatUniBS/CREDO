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
outfile=args[4]

cat(sprintf("RNA count   files: %s\n", rnacntfile))
cat(sprintf("DNA count   files: %s\n", dnacntfile))
cat(sprintf("DNA variant files: %s\n", dnavarfile))


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
#dnacntfile="/data/tempa0d9/TCGA-A7-A0D9_NT_bq25dp4_dnaaref_readcnt_allaz.vcf"
#dnavarfile="/data/tempa0d9/TCGA-A7-A0D9_NT_bqsr_1basesubstitute_ac1_dp5_q50_atog.vcf"
#outfile="/data/tempa0d9/TCGA-A7-A0D9_NT_bq25dp4_sb0.01_az5.0_editedloci.vcf"

#reading variant vcf from exome
cat(sprintf("Reading DNA variant file: %s\n", dnavarfile))
ntexome=readVcf(file=dnavarfile,"b37_lite")
ntexome_names=names(rowRanges(ntexome))
ntexome_names    =vapply(1:length(ntexome_names), function(i) ntexome_names[i]=strsplit(ntexome_names[i],"_")[[1]][1],FUN.VALUE=character(1))
rm(ntexome)
gc()

#real zero nonref G filter
nonrefcut=5.0

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
 
 write.table(stUth_sb_filterposnt1[,c(1,2,3,4,5,6,13,14,15,16)],file=outfile,row.names = FALSE,col.names = TRUE,sep="\t",quote=FALSE)
 









