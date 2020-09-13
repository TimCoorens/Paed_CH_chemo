#Filtering script for "Clonal hematopoiesis and therapy-related myeloid neoplasms following neuroblastoma treatment" (Coorens, Collord, et al.)
#Tim Coorens, September 2020

options(stringsAsFactors = F)
patient="PD31013"
opt="snp"
gender=NA

library(data.table)

source("/lustre/scratch117/casm/team268/tc16/Scripts/R_scripts/germline_exact_binom.R")
source("/lustre/scratch117/casm/team268/tc16/Scripts/R_scripts/beta_binom_flt.R")
source("/lustre/scratch117/casm/team268/tc16/Scripts/R_scripts/plot_spectrum.R")

#Read in output from cgpVAF (any allele counts will do)
data = fread(paste0(patient,".",opt,".tsv"),header=T,data.table=F)
Muts = paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_")
Genotype = data[,grepl("VAF",colnames(data))&colnames(data)!="PDv37is_VAF"]
NR = data[,grepl("DEP",colnames(data))&colnames(data)!="PDv37is_DEP"]
NV = data[,grepl("MTR",colnames(data))&colnames(data)!="PDv37is_MTR"]
rownames(Genotype)=rownames(NV)=rownames(NR)=Muts
samples = colnames(Genotype)=colnames(NR)=colnames(NV)=gsub("_VAF","",colnames(Genotype))

#Select samples from this patient
samples_patient=samples[grepl(patient,samples)]

#Select samples without copy number changes
if(patient=="PD31012")samples_select="PD31012b"
if(patient=="PD31013")samples_select=("PD31013a","PD31013c","PD31013g")
if(patient=="PD42747")samples_select=("PD42747a","PD42747c")

XY_chromosomal = grepl("X|Y",Muts)
autosomal = !XY_chromosomal
xy_depth=mean(NR[XY_chromosomal,samples_select])
autosomal_depth=mean(NR[autosomal,samples_select])

gender='male'
if(xy_depth>0.8*autosomal_depth) gender='female'

if(patient=="PD31012"){
  pdf("PD31012_depth.pdf")
  hist(NR[,samples_select],breaks=150,col='steelblue',xlim=c(0,150),main="PD31012 depth",xlab="Coverage")
  abline(col='red',lwd=2,v=c(50,150),lty='dashed')
  dev.off()
  depth_flt=NR[,samples_select]>50&NR[,samples_select]<150
  
}else{
  pdf(paste0(patient,"_depth.pdf"))
  hist(rowMeans(NR[,samples_select]),breaks=150,col='steelblue',xlim=c(0,150),main=paste0(patient," depth"),xlab="Coverage")
  abline(col='red',lwd=2,v=c(50,110),lty='dashed')
  dev.off()
  depth_flt=rowMeans(NR[,samples_select])>50&rowMeans(NR[,samples_select])<110
}

#Filter out germline variants via an exact binomial test
germline=exact.binomial(gender=gender,NV=NV[,samples_select],NR=NR[,samples_select],cutoff = -5) #determine which variants are germline
NR_flt=NR[!germline&depth_flt,]
NV_flt=NV[!germline&depth_flt,]

NR_flt_nonzero=NR_flt
NR_flt_nonzero[NR_flt_nonzero==0]=1
Genotype_flt=NV_flt/NR_flt_nonzero

#Filter out artefacts by calculating the beta-binomial overdispersion coefficient ("rho")
shared_muts=rowSums(NV_flt>0)>1
rho=beta.binom.filter(NR=NR_flt[shared_muts,samples_patient],NV=NV_flt[shared_muts,samples_patient])

flt_rho=log10(rho)<(-1)
rho_filtered_out = rownames(NR_flt[shared_muts,])[flt_rho]
write.table(rho_filtered_out,"muts_bbinom_filtered_out.txt")
write.table(rho,"rho_est.txt")

NR_flt_2 = NR_flt[!rownames(NR_flt)%in%rho_filtered_out,]
NV_flt_2 = NV_flt[!rownames(NV_flt)%in%rho_filtered_out,]
write.table(NR_flt_2,"NR_filtered_all.txt")
write.table(NV_flt_2,"NV_filtered_all.txt")

#Create bed file to check variants from normal samples in Jbrowse
variants_check=rownames(Genotype_flt2)[rowSums(NV_flt_2[,samples_select])>0]
con = textConnection(variants_check)
Muts_coord = read.table(con,sep="_")
intv=50
Muts_bed = data.frame(Chr=Muts_coord$V1,
                      Start=Muts_coord$V2-intv,
                      End=Muts_coord$V2+intv)
write.table(Muts_bed,"Jbrowse_normal.bed",row.names = F,col.names = F,quote=F,sep="\t")

