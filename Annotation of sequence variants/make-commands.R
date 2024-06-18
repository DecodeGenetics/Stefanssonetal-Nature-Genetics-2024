#######################################################################
# define working folder and filename of index markers (ImpMAF > 0.01) #
#######################################################################

output_folder <- "/.../annotation/"

fileName <- "markers.gor"

name_of_markerlist <- "asm-qtls"

markers_file <- paste0(output_folder, "markers-lead_spaced-by-1Mb.tab")

marker_filename <- gsub(".tab",".gor",markers_file)

path_to_ldfile <- gsub("markers-","LD-to-markers-", marker_filename)

full_set_filename <- "all-markers.gor"

# for bag of random markers, define the following: 

genomic_attributes_file <- "/.../MDSs.gor"

bagName <- "bag-of-random-markers.tab"

path2outputfile <- paste0(output_folder,bagName)

markerset2use <- "/.../markers_v3-varjoin.gor"

path_to_ldrandomfile <- paste0(output_folder, "LD-to-",gsub(".tab",".gor", bagName))

##################################
# define markers (ImpMAF > 0.01) #
##################################

f <- fread(paste0(output_folder,fileName))f[,.(Chrom,Pos,Name,t_abs,ImpMAF)]

if(file.exists(output_folder)==FALSE) system(paste0("mkdir ", output_folder))

options(scipen=999)

write.table(f, file=paste0(output_folder,full_set_filename), quote=FALSE, sep="\t", row.names=FALSE)

# restrict to MAF >1%. 
f <- f[ImpMAF > 0.01, ]; f[,ImpMAF:=NULL]

write.table(f[order(Chrom,Pos), ], file=paste0(output_folder, fileName), quote=FALSE, sep="\t", row.names=FALSE)

#######################
# make sparse markers #
#######################

sys.cmd <- paste0("Rscript /.../make-sparse-marker-lists.R ", output_folder, fileName, " 1e+6 TRUE ", markers_file)
writeLines(sys.cmd)
system(sys.cmd)

# convert output to .gor format

sys.cmd <- paste0("gorpipe \"", markers_file, " | sort genome\" > ", marker_filename)
writeLines(sys.cmd)
system(sys.cmd)

# clean-up
system(paste0("rm ", markers_file))

##############################
# make bag of random markers #
##############################

sys.cmd <- paste0("Rscript /.../make-random-marker-lists.R ", output_folder, fileName, " ", genomic_attributes_file, " 10000 200000 ", path2outputfile, " ", markerset2use)
writeLines(sys.cmd)
system(sys.cmd)

# Show data format of "markerset2use"
# Chrom    Pos             Name ImpMAF             REF                ALT
# <char>  <int>           <char>  <num>          <char>             <char>
# 1: chr1 586844   chr1:586844:SG  0.012               G                  A
# 2: chr1 587429   chr1:587429:SG  0.061               G                  A

#################
# make commands #
#################

library(data.table)

annList <- fread("/.../list_of_annotations.tab")[(MAF==0.01 | is.na(MAF)) & freeze %in% c("any","daisy"), ]

# Show data format of "annList"
# head(annList,2)
# dataset_id
# <char>
# 1: ENCODE-ChIPseq-ChromatinAttributes
# 2: ENCODE-cCRE_V3-CellTypeAgnostic
# path2gord     N  isSNP freeze
# <char> <int> <lgcl> <char>
# 1: /.../path2gorfiles_ids_full.gord  2696  FALSE    any
# 2: /.../cCRE_V3-cell-type-agnostic-hg38.gord     9  FALSE    any
# MAF
# <num>
# 1:    NA
# 2:    NA
# Description
# <char>
# 1:                                     DNA binding proteins, histone marked regions in various cell types and tissues
# 2: Encode classification of open chromatin regions into promoter-like, enhancer-like or insulator regulatory elements

# filter on dataset_id
annList <- annList[dataset_id %in% c("ENCODE-cCRE_V3-CellTypeAgnostic","AlleleDB-ASBs_in_daisy-Lymphoblasts")]

##########################################################
# list out settings for running the enrichment algorithm #
##########################################################

target_list <- NULL; gord <- NULL; so <- NULL; issnp <- NULL; maf_setting <- NULL; dupl <- NULL

for(i in 1:nrow(annList)){
  
tmp <- as.character(fread(as.character(annList$path2gord[i]), header=FALSE)$V2)
target_list <- c(target_list, tmp)

# detect duplicates in .gord files
dupl <- c(dupl, tmp[duplicated(tmp)])

gord <- c(gord, rep(annList$path2gord[i], length(tmp)))

so <- c(so, rep(annList$dataset_id[i], length(tmp)))

issnp <- c(issnp, rep(annList$isSNP[i], length(tmp)))

maf_setting <- c(maf_setting, rep(annList$MAF[i], length(tmp)))

}

# arg1: 
path_to_markers <- output_folder

# arg2: already defined above
marker_filename <- gsub(path_to_markers, "", marker_filename)

# Show data format of "marker_filename" file
# Chrom       Pos                Name     nlog10P
# <char>     <int>              <char> <num>
# 1: chr1    901812      chr1:901812:SG    35
# 2: chr1   1955018 chr1:1955018:IG.0:1   119

# arg3: already defined above
path_to_ldfile

# Show data format of "path_to_ldfile"
# Chrom     Pos                Name     P    pos2                 mk1               mk2  maf1
# <char>   <int>              <char> <num>   <int>              <char>            <char> <num>
# 1:   chr1  901812      chr1:901812:SG    35  901812      chr1:901812:SG    chr1:901812:SG 0.078
# 2:   chr1 1955018 chr1:1955018:IG.0:1   119 1950071 chr1:1955018:IG.0:1 chr1:1950071:SG:1 0.315
# maf2   rsq
# <num> <num>
# 1: 0.078  1.00
# 2: 0.285  0.85

# arg4: 
path_tmp_folder <- paste0(output_folder,"tmp/")

# arg5
Nresamplings <- 1000

# arg6: 
head(target_list)

# arg7:
name_of_markerlist

# arg8: already defined above
output_folder

# arg9: 
nseed <- 999

# arg10: print out resampling table?
logical1 <- FALSE
  
# arg11 print out annotation table?
logical2 <- FALSE
  
# arg12
head(gord)

# arg13
head(so)

# arg14
head(issnp)

# arg15
head(maf_setting)

# arg16
maxlead <- 2000

# arg17; what type of markers are in "arg2"?
markerType <- "QTL"

# arg18
background_snps_list <- path_to_ldrandomfile

# Show data format of "background_snps_list"
# Chrom    Pos           Name   pos2            mk1            mk2  maf1  maf2   rsq
# <char>  <int>         <char>  <int>         <char>         <char> <num> <num> <num>
# 1: chr1 899452 chr1:899452:SG 897843 chr1:899452:SG chr1:897843:SG  0.25  0.25     1
# 2: chr1 899452 chr1:899452:SG 897922 chr1:899452:SG chr1:897922:SG  0.25  0.25     1

# set arguments
arg_settings <- NULL; for(i in 1:length(gord)) { arg_settings <- c(arg_settings, paste0(path_to_markers, " ", marker_filename, " ", path_to_ldfile, " ", path_tmp_folder, " ", Nresamplings  ," ", target_list[i], " ", name_of_markerlist, " ", output_folder, " ", nseed, " ", logical1, " ", logical2, " ", gord[i], " ", so[i], " ", issnp[i], " ", maf_setting[i], " ", maxlead," ", markerType, " ", background_snps_list)) }

# must be "TRUE"
length(target_list)==length(arg_settings)

# set path to script:
script2use <- "/.../run.R"

# make command lines
sys.cmd <- paste("Rscript", script2use, arg_settings)

# show command lines
print(sys.cmd)
