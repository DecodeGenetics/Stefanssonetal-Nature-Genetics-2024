a <- commandArgs(TRUE)

# usage $ Rscript /.../run.R <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7> <arg8> <arg9> <arg10> <arg11> <arg12> <arg13> <arg14> <arg15> <arg16> <arg17> <arg18>

# arg1: path to file containing index markers (in house variant names in 'Name' column).... tab-separated file format: Chrom	Pos	Name
# arg2: name of file containing index markers (Chrom, Pos, Name)
# arg3: path to file containing LD information (Chrom=chromosome name, Pos=index marker position, mk1=index marker name, mk2=correlated to index marker, pos2=position of mk2, maf1=frequency of mk1)
# arg4: folder for writing out temporary files
# arg5: number of re-samplings used to determine empirical p-value (min=200)
# arg6: name of sequence variant annotation
# arg7: run identifier (to make it easy to collect output files)
# arg8: path to folder where output files are to be written
# arg9: number for set.seed
# arg10: logical TRUE/FALSE statement: print out the resampling tables? Rarely used, please set this option to FALSE. 
# arg11: logical TRUE/FALSE statement: print out annotation table? Rarely used, please set this option to FALSE. 
# arg12: path to .gord file containing the annotations to use in enrichment analysis
# arg13: name of annotation source
# arg14: isSNP logical TRUE/FALSE statement: is annotation SNP-based? TRUE=SNP-based, FALSE=Region-based
# arg15: if arg14='TRUE' then set this argument as 0.01 (MAF threshold; >1%)
# arg16: set max number of lead association variants in input marker file (set in arg2); set to 'Inf' in most cases.  If marker number is very high (e.g. >1000) the re-sampling process will be very slow in which case this argument will help speed up the run time. 
# arg17: set type of association list; i.e. "QTL" or "GWAS"... if set to "QTL" then arg18 is required
# arg18: if arg17="QTL"... path to file containing "bag of random markers" defined using command 'make-random-marker-lists.R'
# ...... if arg17="GWAS".. set arg18 to "NA"

# Example run:
# $ Rscript /.../run.R /.../enrichment-of-annotation/ markers-lead_spaced-by-1Mb.gor /.../LD-to-markers-lead_spaced-by-1Mb.gor /.../tmp/ 1000 AlleleDB_CTCF asm-qtls /.../enrichment-of-annotation/ 999 FALSE FALSE /.../AlleleDB-daisy.gord AlleleDB-ASBs_in_daisy-Lymphoblasts TRUE 0.01 1000 QTL /.../LD-to-bag-of-random-markers.gor

# Testing arguments (QTL setting)
a <- c("/.../enrichment-of-annotation/", 
       "markers-lead_spaced-by-1Mb.gor", 
       "/.../LD-to-markers-lead_spaced-by-1Mb.gor",
       ".../tmp/",
       1000,
       "AlleleDB_CTCF",
       "asm-qtls",
       "/.../enrichment-of-annotation/",
       999,
       FALSE,
       FALSE,
       "/.../AlleleDB-daisy.gord",
       "AlleleDB-ASBs_in_daisy-Lymphoblasts",
       TRUE,
       0.01,
       2000,
       "QTL",
       "/.../LD-to-bag-of-random-markers.gor")

# Detailed documentation of arguments: 

# <arg1>: the folder name is specified in <arg1>

# <arg2>: list of lead markers for each GWAS or QTL association signal, the lead markers should be separated by at least 1Mb which can be performed using "make-sparse-marker-lists.R". 
# the data format of "markers-lead_spaced-by-1Mb.gor" under argument <arg2> is as follows: 
# Chrom      Pos             Name     nlog10p
# chr1 24065476 chr1:24065476:SG    74

# <arg3>: the data format of "/.../LD-to-markers-lead_spaced-by-1Mb.gor" under argument <arg3> is as follows: 
# fread("/.../LD-to-markers-lead_spaced-by-1Mb.gor")[1]
# Chrom   Pos     Name       pos2    mk1     mk2     maf1    maf2    rsq
# chr10 31660483 chr10:31660483:SG      3019822 chr10:31660483:SG  chr10:31623818:M:2  0.19   0.2  0.93

# <arg4>: specify the folder name where temporary files will be written

# <arg5>: specify the number of sampled sets to be performed, by default set to 1000

# <arg6>: specify the name of the sequence variant annotation to test in <arg6> corresponding to an ID in .gord file as specified under <arg12>

# <arg7>: specify a name that will be added to the output file in <arg7> (this can be any name specified by the user)

# <arg8>: the folder name where output files will be written is specified in <arg8>

# <arg9>: a random number choosen by the user to be used in set.seed() command is specified in <arg9>

# <arg10> is a rarely used option, by default it should be set to FALSE, but this option can be set to TRUE to perform large number of sampled sets (e.g., 50000). 

# <arg11> is a rarely used option, by default it should be set to FALSE, but it can be set to TRUE to see annotations for the input markers.  

# <arg12> is the .gord file containing the sequence variant annotation specified in <arg6> is specified in <arg12>
# the data format of the "/.../AlleleDB-daisy.gord" file is as follows:
# fread("/.../AlleleDB-daisy.gord", header=F)
# V1                    V2
# <char>                <char>
# 1: /.../ASB-to-daisy_hg38_CTCF.gor         AlleleDB_CTCF
# 2: /.../ASB-to-daisy_hg38_EBF.gor          AlleleDB_EBF
# 3: /.../ASB-to-daisy_hg38_POL2.gor         AlleleDB_POL2
# 4: /.../ASB-to-daisy_hg38_PU1.gor          AlleleDB_PU1
# 5: /.../ASB-to-daisy_hg38_RPB2.gor         AlleleDB_RPB2
# 6: /.../ASB-to-daisy_hg38_SA1.gor          AlleleDB_SA1
# 7: /.../ASB-to-daisy_hg38.gor              AlleleDB_full-dataset

# the data format of the .gor files, e.g., "/.../ASB-to-daisy_hg38_CTCF.gor" is as follows: 
# fread("/.../ASB-to-daisy_hg38_CTCF.gor")[1:3]
# Chrom     Pos              Name
# <char>   <int>            <char>
# 1:   chr1 3806639   chr1:3806639:SG
# 2:   chr1 3806663 chr1:3806663:SG:1
# 3:   chr1 4027648   chr1:4027648:SG

# <arg13>: name of annotation source in file "/.../list_of_annotations.tab"
# the data format of the "/.../list_of_annotations.tab" file is: 
# fread("/.../list_of_annotations.tab")[dataset_id %in% c("AlleleDB-ASBs_in_daisy-Lymphoblasts","ENCODE-cCRE_V3-CellTypeAgnostic")]
# dataset_id                                                                  path2gord     N  isSNP freeze   MAF
# <char>                                                                     <char> <int> <lgcl> <char> <num>
# 1: ENCODE-cCRE_V3-CellTypeAgnostic  /.../cCRE_V3-cell-type-agnostic-hg38.gord     9  FALSE    any    NA
# 2: AlleleDB-ASBs_in_daisy-Lymphoblasts /.../AlleleDB-daisy.gord     7   TRUE  daisy  0.01
# Description
# <char>
# 1: Encode classification of open chromatin regions into promoter-like, enhancer-like or insulator regulatory elements
# 2:                                MAF > 0.01, AlleleDB: Allele-specific binding sites (ASB) in B-lymphoblastoid cells

# <arg14> is a logical variable (TRUE/FALSE) which should be set accoding to the "isSNP" column in "/.../list_of_annotations.tab" file. 
# ... specifically, if the .gord annotation file referred to in <arg12> refers to .gor files that are based on sequence variant associations (Chrom,Pos,Name) then <arg14> should be set to TRUE, otherwise the .gord file is refering to .gor files that are based on genomic regions (Chrom,Start,End) in which case <arg14> should be set to FALSE. 

# <arg15> should be set to 0.01 if <arg14> is TRUE, i.e., the MAF column in "/.../list_of_annotations.tab" file, which thereby restricts sequence variants to MAF >1%, but <arg15> should otherwise be set to NA

# <arg16>: set the max number of lead association variants in input marker file (set in arg2); set to 'Inf' in most cases.  If marker number is very high (e.g. >1000) the re-sampling process will be very slow in which case this argument will help speed up the run time. 

# <arg17>: set type of association list; i.e. "QTL" or "GWAS"... if set to "QTL" then arg18 is required

# <arg18>: if <arg17> is set to "QTL", then <arg18> is required and should be set to the path to the file containing a "bag of random markers" (which can be defined using command 'make-random-marker-lists.R')
# ......   if <arg17> is set to "GWAS" then arg18 is not required and is set to "NA"

######################
# required libraries #
######################

suppressPackageStartupMessages(library(GenomicRanges)); suppressPackageStartupMessages(library(data.table)); suppressPackageStartupMessages(library(dplyr)); suppressPackageStartupMessages(library(tictoc))

################################################################################################################
# help function (used to ensure sequence variants in randomly re-sampled marker list will be spaced 1Mb apart) #
################################################################################################################

make.sparse.markerlist <- function(p, minDist){
  
  u <- unique(p$Chrom)
  
  pin.c <- NULL
  
  for(i in 1:length(u)){
    
    pp <- p[Chrom%in%u[i], ]
    
    pin <- pp[1, ]
    
    if(nrow(pp)>1){  
      
        for(z in 2:nrow(pp)){
          
          pc <- pp[z, ]
          
          len <- abs(pc$Pos-pin[,Pos]); if(all(len>minDist)) pin <- rbindlist(list(pin,pc))
          
        }
      
    }
    
    pin.c <- rbindlist(list(pin.c, pin))
    
  }
  
  return(pin.c)
  
}


#################
# set arguments #
#################

# path to association markers and filename
path2file <- as.character(a[1]); filename <- as.character(a[2])

# path to LD information for association markers
path2ld <- as.character(a[3])

# set number of the column containing "in house" sequence variant identifiers
marker_column <- 3

# set location of path for writing out temporary files: 
tmp_path <- as.character(a[4])

# set number of re-samplings
n_samples <- as.numeric(a[5])

# set tissue / cell type 
eName <- a[6]

# 'gord' file
gord <- as.character(a[12]); cstats <- fread(gord, header=F); colnames(cstats) <- c("path2gorfile","id")

# set output path
path2output <- as.character(a[8])

# path-to-results by name of annotation source ('dataset_id' column in file: "list_of_annotations.tab")
annotation_id <- as.character(a[13])

path2results <- paste0(path2output,annotation_id,"/"); if(file.exists(path2results)==FALSE) system(paste0("mkdir ", path2results))

if(file.exists(path2results)==FALSE) system(paste0("mkdir ",path2results))

# make random number for temporary file
randomN <- sample(1:1000000,1)

# set run identifier
run_id <- as.character(a[7])

# number for 'set.seed' command
setseed <- as.numeric(a[9])

# logical;TRUE=write out resmpling tables for polishing the P-value; set to FALSE for scanning
writeTable <- as.logical(as.character(a[10]))

# logical; TRUE=annotation table will be written out and the program will then exit without any computations carried out; set to FALSE for scanning
annotationTable <- as.logical(as.character(a[11]))

# logical; TRUE=annotation is SNP-based e.g. eQTLs,sQTLs...etc, FALSE=annotation is regions-based e.g. DHS regions in reticulocytes
isSNP <- as.logical(as.character(a[14]))

# if SNP-based annotation, then specify the MAF threshold as either 0.01 or 0.05
mafSNP <- as.numeric(as.character(a[15]))

# set max number of association markers to be included; set to 2000 in most cases (>2000 will severely impact the run time)
maxLeads <- as.numeric(as.character(a[16])); if(is.infinite(maxLeads) | maxLeads>2000) maxLeads <- 2000

# set to either "GWAS" or "QTL"
isBackground <- as.character(a[17])

# if 'isBackground'="QTL" then set 'bag_of_randomvariants':
bag_of_randomvariants <- as.character(a[18])

# if 'isBackground'="GWAS" then set classification of GWAS signals as 'gwas_classification' (NOTE: set to NA for including all GWAS signals): 
gwas_classification <- as.character(a[18]); if(gwas_classification=="NA") gwas_classification <- NA

########################
# load neccessary data #
########################

# GWAS catalog, lead variants and variants in strong LD (r2>0.80)
if(is.na(isBackground)) isBackground <- "GWAS"

if(isBackground=="GWAS"){

  gwas <- fread("/.../GWAS_catalog_EuropeanAncestry_P1E-9_1Mb-spacing-TraitClass_LD-tidy0.gor")
  
    if(is.na(gwas_classification)==FALSE){ print(paste0("Classification: ", gwas_classification)); gwas_classification <- unlist(tstrsplit(gwas_classification,",")); gwas <- gwas[Classification%in%gwas_classification] }
    
  }

# QTLs as the background
if(isBackground=="QTL"){ print(paste0("loading bag of random variants: ", bag_of_randomvariants)); gwas <- fread(bag_of_randomvariants); colnames(gwas)[c(1,4,3,6,7,8,9)] <- c("Chrom","inLD_pos","lead_name","inLD_name","ImpMAF_lead","inLD_maf","r2"); gwas <- left_join(gwas, gwas[,.(LD_class=.N),lead_name][,Phenotype:=paste0("make-believe-phenotype",1:length(unique(gwas$lead_name)))][,signal:=paste0(lead_name,"-",Phenotype)], by="lead_name"); gwas[ImpMAF_lead>0.5, ImpMAF_lead:=1-ImpMAF_lead]; gwas[inLD_maf>0.5, inLD_maf:=1-inLD_maf] }

# define 'gwas_index'
gwas_index <- gwas[lead_name==inLD_name, ]

# set LD class category
gwas_index[,LD_class_category:=cut(LD_class,breaks=c(0,1,5,10,20,50,100,200,Inf))]

# load index markers
print(paste0("loading marker file: ", paste0(path2file, filename)))
ind <- fread(paste0(path2file, filename)); colnames(ind)[1:2] <- c("Chrom","Pos")

# set "SNP" column
colnames(ind)[marker_column] <- "SNP"

# max 2000 lead variants
if(is.na(maxLeads)) maxLeads <- 2000; if(nrow(ind)<maxLeads) maxLeads <- nrow(ind); if(nrow(ind>maxLeads)){ ind <- ind[sample(1:nrow(ind), maxLeads)]  } 

# load LD file for association markers
print(paste0("loading LD file: ", path2ld))
ann <- fread(path2ld); ann <- ann[mk1%in%ind[,SNP], ]

# maf representation
ann[maf1>0.5, maf1:=1-maf1]; ann[maf2>0.5, maf2:=1-maf2]

# LD class
ann_ldclass <- ann[,.N,mk1]; colnames(ann_ldclass)[2] <- "LD_class"; ann <- left_join(ann,ann_ldclass, by="mk1")

lead <- ann[mk1==mk2 & mk1%in%ind[,SNP], ][duplicated(paste0(mk1,"-",mk2))==F]

# classify LD_class among "lead" variants
ldclass_find <- NULL; for(w in 1:nrow(lead)){ di <- abs(gwas_index$LD_class-lead$LD_class[w]); ldclass_find <- c(ldclass_find, as.character(gwas_index[which(di==min(di))[1], LD_class_category])) }

# set LD class to category
lead[, LD_class_category:=ldclass_find]

# load annotation files
elistnames <- as.character(cstats[,.N,id][,id])

# if SNP-based annotation, then use common variants based on threshold set in 'mafSNP' (either >0.05 or >0.01)
if(isSNP==TRUE){
  
  # filter on maf 
  if(is.na(mafSNP)) mafSNP <- 0.01; lead <- lead[maf1>mafSNP, ]

  ann <- ann[mk1 %in% lead$mk1, ]
  ann_ldclass <- ann_ldclass[mk1 %in% lead$mk1, ]
  ind <- ind[SNP %in% lead$mk1, ]
  
  # filter on maf 
  gwas <- gwas[ImpMAF_lead>mafSNP, ]
  gwas_index <- gwas_index[ImpMAF_lead>mafSNP, ]
  
}

########
# run  #
########

nOut <- NULL; fr <- NULL; N2save <- NULL
  
genome.x <- cstats[id%in%eName, ]; path2genome.x <- as.character(genome.x[,path2gorfile])

#  define " number of annotated loci " for associations

annoFile <- paste0(tmp_path, "anno_", eName,"_",randomN,".gor")

print("write to disk: intersection to annotation file (gorpipe command line)")

# GWAS
if(isSNP==FALSE){
sys.cmd <- paste0("gorpipe \"", path2ld, " | select Chrom,pos2,mk2,mk1 | sort genome | join -snpseg ", path2genome.x,"\" > ", annoFile)
writeLines(sys.cmd)
system(sys.cmd)
}

# QTL
if(isSNP==TRUE){
  sys.cmd <- paste0("gorpipe \"", path2ld, " | select Chrom,pos2,mk1,mk2,maf1 | sort genome | join -snpsnp ", path2genome.x,"\" > ", annoFile)
  writeLines(sys.cmd)
  system(sys.cmd)
}

# load annotation file
print("loading gor file containing the observed intersections")
ann_genome <- fread(annoFile); ann_genome <- ann_genome[mk1%in%lead[,mk1]]

# if annotation is SNP-based then ONLY allow "common" SNPs

if(isSNP){ 
  if(("ImpMAF"%in%colnames(ann_genome))==FALSE){ colnames(ann_genome)[grep("maf1", colnames(ann_genome))] <- "ImpMAF" }
  }

if(annotationTable){
  
  ann_genome[,pseudoCol:=rep(as.character(eName),nrow(ann_genome))]
    atable <- dcast(ann_genome, Col4~pseudoCol,value.var="Col5", fun.aggregate=function(x)paste0(x[duplicated(x)==FALSE],collapse=","))
      annTmp <- ann[Col4==Col5,.(Col4,Col6,LD_class)]; colnames(annTmp)[2] <- "MAF_lead"
  
  atable <- left_join(atable,annTmp,by="Col4"); colnames(atable)[1] <- "lead"; atable <- atable[order(lead), ]
  
  a.outpath <- paste0(path2output,"remapPeaks-annotation-tables/")
    
  if(file.exists(a.outpath)==FALSE) system(paste0("mkdir ", a.outpath))

  file.out <- paste0(a.outpath,run_id,"_chromatin_",eName,".tab"); write.table(atable, file=file.out, quote=FALSE, sep="\t", row.names=FALSE)
  
}

# stop running if only annotation table is needed
if(annotationTable) stop(paste0("annotation table written here: ", file.out))

print(paste0("loading full annotation file: ", path2genome.x))

# GWAS
if(isSNP==FALSE){
  v <- fread(path2genome.x); colnames(v)[1:3] <- c("Chrom","Start","End")
    v.gr <- reduce(makeGRangesFromDataFrame(v))
      annotationLength <- sum(width(reduce(v.gr)))
}

# QTL
if(isSNP==TRUE){
  v <- fread(path2genome.x); colnames(v)[1:2] <- c("Chrom","Pos"); v[, coords:=paste0(v$Chrom,":",v$Pos)]
    v <- v[duplicated(coords)==FALSE, ]
  
  v.gr <- GRanges(seqnames=v$Chrom, IRanges(start=v$Pos,end=v$Pos)); annotationLength <- length(unique(v$coords))
}

# store info
fr <- data.table(facet=eName, coverage_bp=annotationLength)

# define 'n'
n <- data.table(facet=eName, N=nrow(ann_genome[,.N,by=mk1])); colnames(n) <- c("V1","V2")

# define obsIntersection (fixed number)
obsIntersection <- n$V2/length(unique(lead$mk1)); obsIntersection_N <- n$V2

# simulation
set.seed(setseed)

# set max number of re-sampling (for loops) in "maxloops", can be greater than requested in "n_samples"... if p-value exactly 0, then for loop will continue until no longer exactly zero but for no more than "maxloops"
maxloops <- n_samples; if(n_samples>maxloops){ print(paste0("max re-sampling number is: ", maxloops)); n_samples <- maxloops }

# min 200 iterations
if(n_samples<200) { n_samples <- 200; print(paste0("min re-sampling number is: ", 200))}

if(n_samples>=200) p_break <- seq(from=100, to=maxloops, by=100); p_tmp <- 0; p_tmp_depl <- 0; tic.start <- c(1, p_break[-length(p_break)]+1) 

N <- NULL; is_p_non_zero <- FALSE; pDone <- FALSE

print(paste0("re-sampling: ", n_samples))
print(paste0("number of lead markers: ", nrow(lead)))

#####################
# start re-sampling #
#####################

for(a in 1:maxloops){

  if(a%in%tic.start){ tic() }
  
  gwas_tmp0 <- GRanges(seqnames=gwas_index$Chrom, IRanges(start=gwas_index$inLD_pos, end=gwas_index$inLD_pos)); mcols(gwas_tmp0) <- gwas_index[,-(1:2)]
  gwas_sim <- NULL; gwas_sim_pheno <- NULL
  
  ldBins <- lead[,.N,LD_class_category]

  # select "random" sequence variants into "gwas_sim" to estimate the expected overlap
  
  while(length(gwas_sim)<nrow(lead)){
    
      gwas_tmp <- as.data.table(gwas_tmp0); gwas_tmp[,end:=NULL]; colnames(gwas_tmp)[1:2] <- c("Chrom","inLD_pos")

      # ensure "random" sequence variants have the same LD class as the observed "lead" variants
        for(i in 1:nrow(ldBins)){
          
          tmp <- gwas_tmp[LD_class_category %in% ldBins[i, LD_class_category], ]
          
          Nl <- ldBins[i, N]
          
          # randomly select two-times more (Nl*2) than needed to reduce run speed (Note, "gwas_sim" will later be reduced to the same length as number of lead markers)
          resampling.ids <- sample(1:nrow(tmp), Nl*2, replace=TRUE)
          
          randomLead <- tmp[resampling.ids, lead_name]; randomPheno <- tmp[resampling.ids, Phenotype]
          
          gwas_sim <- c(gwas_sim, randomLead);  gwas_sim_pheno <- c(gwas_sim_pheno, randomPheno)
          
        }
     
      tmpDT <- data.table(Chrom=unlist(tstrsplit(gwas_sim,":",keep=1)), Pos=as.numeric(unlist(tstrsplit(gwas_sim,":",keep=2))), Name=gwas_sim, gwas_sim_pheno); tmpDT <- tmpDT[order(Chrom,Pos), ]
      
      # ensure "random" sequence variants are separated by at least 1Mb as is the requirement for the "lead" variants
      
      tmpDT.gr <- GRanges(seqnames=tmpDT$Chrom, IRanges(start=tmpDT$Pos, end=tmpDT$Pos)); tmpDT.gr$Name <- tmpDT$Name
      tmpDT.r <- reduce(tmpDT.gr, min.gap = 1e+6)
      tmpDT <- tmpDT[paste0(Chrom,":", Pos) %in% paste0(seqnames(tmpDT.r), ":", start(tmpDT.r)), ]
      
      # tmpDT <- make.sparse.markerlist(tmpDT, 1e+6)
      
      gwas_sim <- as.character(tmpDT$Name); gwas_sim_pheno <- as.character(tmpDT$gwas_sim_pheno)
      
      tmpDT.gr <- GRanges(seqnames=tmpDT$Chrom, IRanges(start=tmpDT$Pos-0.5e+6, end=tmpDT$Pos+0.5e+6)); mcols(tmpDT.gr) <- tmpDT[,-(1:2)]
       
      gwas_tmp0 <- gwas_tmp0[(gwas_tmp0%over%tmpDT.gr)==FALSE]
      
    }
  
  # ensure 'gwas_sim' has the same length as number of lead markers
  finalsel <- sample(1:length(gwas_sim), length(unique(lead$mk1))); gwas_sim <- gwas_sim[finalsel]; gwas_sim_pheno <- gwas_sim_pheno[finalsel]
  
  rsq <- gwas[r2>0.8 & lead_name%in%gwas_sim & Phenotype %in% gwas_sim_pheno, ]
 
  rsq[,uniq:=paste0(rsq[,inLD_name],"-",rsq[,signal])]; rsq <- rsq[duplicated(uniq)==FALSE, ]
  
    r.gr <- GRanges(seqnames=rsq[,Chrom], IRanges(start=rsq[,inLD_pos], end=rsq[, inLD_pos]))
  
  mcols(r.gr) <- rsq[,lead_name]
  
  lNum <- length(unique(r.gr[r.gr%over%v.gr]$X))

  tmp <- data.frame(facet=eName, count=lNum)
 
  N <- rbind(N,tmp)

  ind_resampled <- ind[sample(1:nrow(ind),nrow(ind), replace=T), ]
  
  if(n_samples>=100){
    
    if(a %in% p_break){
      
      # quick look at p-value after 'k' re-sampling has been carried out
      p_tmp <- nrow(N[N$count >= n$V2, ]) / nrow(N); p_tmp_depl <- nrow(N[N$count <= n$V2, ]) / nrow(N)
      
      # display some progress information!
      print( paste0("P-value after k=", a, " re-samplings for enrichment: ", p_tmp, " and for depletion: ", p_tmp_depl) )
      print( "timing reported each 100 iterations:"  ); toc()
      
    }
    
  }
  
  # see if further re-sampling is worthwile given status of p-value after 'k' re-samplings have been carried out, but only if writeTable='FALSE'
  if(writeTable==FALSE){

          if(p_tmp>=0.05 & p_tmp_depl>=0.05 & a >= 200) pDone <- TRUE
          
          if(p_tmp>=0.01 & p_tmp_depl>=0.01 & a >= 500) pDone <- TRUE
          
          if(p_tmp>=0.0025 & p_tmp_depl>=0.0025 & a >= 800) pDone <- TRUE

          # if there is one, or less than one observed intersection then stop running and print out the results after 200 samples
          if(obsIntersection_N <=1 & a>=200) pDone <- TRUE

  }

is_p_non_zero <- a>=n_samples & p_tmp!=0 & p_tmp_depl!=0

if(is_p_non_zero) break

if(pDone) break

}

print("processing re-sampling results")

# process results 

N.keep <- N; N <- as.data.table(N); setkey(N, facet)
setkey(n, V1); N <- n[N,]

# N
pval <- NULL; count_above <- NULL; pval_depl <- NULL; count_below <- NULL

  pval <- c(pval, nrow(N[V1==n$V1 & count >= V2, ]) / nrow(N[V1==n$V1, ]))
  count_above <- c(count_above, nrow(N[V1==n$V1 & count >= V2, ])) 

  pval_depl <- c(pval_depl, nrow(N[V1==n$V1 & count <= V2, ]) / nrow(N[V1==n$V1, ]))
  count_below <- c(count_below, nrow(N[V1==n$V1 & count <= V2, ])) 
  
n$enr_p_value <- signif(pval, digits=4)
n$count_above <- count_above

n$depl_p_value <- signif(pval_depl, digits=4)
n$count_below <- count_below

n$N_resampling <- a

n$N_markers <- nrow(ind)

# confidence intervals
ci <- NULL; m <- NULL

  exp_ratio <- N[V1==n$V1, count]
  
  frac.asis <-  exp_ratio / length(unique(lead$mk1))
  frac.asis <- frac.asis[is.infinite(frac.asis)==FALSE & is.na(frac.asis)==FALSE]
  
  ci <- rbind(ci, quantile(frac.asis, probs=c(0.025,0.975)))
  
  # "m" = expected intersection: the sum of how many found in intersection over all samples divided by the total number of loci sampled
  m <- sum(N$count)/(nrow(N)*length(unique(lead$mk1)))
  
colnames(ci) <- c("lower","upper")
n <- cbind(n,ci)
n$expected_intersection <- signif(m, digits=4)
n$expected_CI95 <- paste0(signif(n$"lower",digits=4), ",", signif(n$"upper", digits=4))
n$observed_intersection <- signif(n$V2/length(unique(lead$mk1)), digits=4)

n[,lower:=NULL]; n[,upper:=NULL]

colnames(n)[1:2] <- c("facet","N_loci_intersection")

final <- left_join(n, fr, by="facet")

# enrichment as computed in version 4
final$fold_enrichment <- signif(final$observed_intersection/final$expected_intersection, digits=4)

final$pvalue_dbinom_enr <- sum(dbinom(final$N_loci_intersection:final$N_markers, final$N_markers, prob = final$expected_intersection))

final$pvalue_dbinom_depl <- sum(dbinom(0:final$N_loci_intersection, final$N_markers, prob = final$expected_intersection))

#########################
# write results to disk #
#########################

print(paste0("writing results to disk: ", path2results))

# clean-up temporary file
system(paste0("rm ", annoFile))

# write out results: 

if(writeTable==FALSE) write.table(final, file=paste0(path2results,run_id,"_", eName, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

colNames <- paste0(path2results,"column_names.txt"); if(file.exists(colNames)==FALSE) write.table(colnames(final), file=colNames, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

if(writeTable==TRUE){

  tableFolder <- paste0(path2results,"resampling-tables/")

  if(file.exists(tableFolder)==FALSE) system(paste0("mkdir ", tableFolder))
  
  N.keep$N_loci_intersection <- final$N_loci_intersection
  N.keep$N_leads_used <- length(unique(lead$mk1))

write.table(N.keep, file=paste0(tableFolder,"Pval-", run_id,"_output_seed", setseed,"_", eName, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

obsTable <- data.table(resampl_observed_intersection=signif(obs_ci, digits=5), observed_intersection=signif(final$observed_intersection, digits=5))
write.table(obsTable, file=paste0(tableFolder,"ObsIntersection-", run_id,"_output_seed", setseed,"_", eName, ".txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

}

print("finished!")
