a <- commandArgs(TRUE)

# usage $ Rscript /.../annotate-signals.R <arg1> <arg2> <arg3>

# arg1: folder containing file listing out the markers in .gor format (Chrom, Pos, Name)
# arg2: name of .gor file
# arg3: <optional> file (no header) containing a list of annotation dataset_id in "/.../list_of_annotations.tab"

# Example run: $ Rscript /.../annotate-signals.R /.../enrichment-of-annotation/ markers.txt /.../annotate-signals.list

# "/.../annotate-signals.list" is a list (without header) corresponding to dataset_id in "/.../list_of_annotations.tab"
# $ head -n3 /.../annotate-signals.list
# ENCODE-ChIPseq-ChromatinAttributes
# ENCODE-cCRE_V3-CellTypeAgnostic
# Remap2022-ChIPseq-ChromatinAttributes

# Show data format of "/.../list_of_annotations.tab"
# dataset_id
# <char>
# 1: ENCODE-ChIPseq-ChromatinAttributes
# 2: ENCODE-cCRE_V3-CellTypeAgnostic
# path2gord     N  isSNP
# <char> <int> <lgcl>
# 1: /.../path2gorfiles_ids_full.gord  2696  FALSE
# 2: /.../cCRE_V3-cell-type-agnostic-hg38.gord     9  FALSE
# freeze   MAF
# <char> <num>
# 1:  any    NA
# 2:  any    NA
# Description
# <char>
# 1: DNA binding proteins, histone marked regions in various cell types and tissues
# 2: Encode classification of open chromatin regions into promoter-like, enhancer-like or insulator regulatory elements

suppressPackageStartupMessages(library(data.table))

# specify input folder containing .gor file listing out the markers to be annotated
input_folder <- as.character(a[1])

# specify filename of .gor file containing markers to be annotated
markers_filename <- as.character(a[2])

# specify output folder
output_folder <- paste0(input_folder, "output/")

# specify rp folder
rp_folder <- paste0(input_folder, "rp_folder/")

# make folders if not already existing
if(file.exists(output_folder)==F) system(paste0("mkdir ", output_folder))

if(file.exists(rp_folder)==F) system(paste0("mkdir ", rp_folder))

# specify path to .gor file
markers_gor <- paste0(input_folder, markers_filename)

# annotation files to be used
l <- fread(a[3], header=F)$V1

# if there is no file is specified in <arg3> then annotate according to "/.../list_of_annotations.tab"
if(is.na(l[1])){ annList <- fread("/.../list_of_annotations.tab")[freeze %in% c("any","daisy"), ] } else { annList <- fread("/.../list_of_annotations.tab")[dataset_id %in% l] }

# exclude annotation files that have already been analysed according to files present in "output_folder"
d <- dir(output_folder); annList <- annList[ (dataset_id %in% gsub(".gor","", d))==FALSE ]

if(nrow(annList)==0){stop("Sequence variants already annotated in the output folder")}

# specify LD file to use
ld_file <- "/.../rsq_high.gord"

# Show data format of "rs_high.gord"
# fread("/.../rsq_high.gord", header=F)
# V1
# <char>
# 1:                                           v3/rsq_high_v0.gorz
# 2:                         v3/ld_ONT_large_tandem_repeat_v1.gorz
# 3:                            v3/rp_results_all_v2a.new/sym.gorz
# 4:                            v3/rp_results_all_v2b.new/sym.gorz
# 5:                            v3/rp_results_all_v2c.new/sym.gorz
# 6: v3/rp_results_phasedseq_all_b4_all_b5_loadbalance_v3/sym.gorz
# 7:                                     v3/sym_v3_singletons.gorz
# 8:                          v4/rp_results_all_b7_all_b5/sym.gorz
# 9:                                     v4/sym_v4_singletons.gorz

# Show data format of an example .gorz file from "/.../rsq_high.gord": 
# $ gorpipe " /.../rsq_high_v0.gorz | first 3"
# chrom	pos1	pos2	mk1	mk2	maf1	maf2	rsq
# chr1	55385	55385	chr1:55385:SG	chr1:55385:SG	0.000016	0.000016	1.000000
# chr1	55385	191818	chr1:55385:SG	chr1:191818:SG	0.000016	0.000016	1.000000
# chr1	55385	788807	chr1:55385:SG	chr1:788807:IG:1	0.000016	0.000017	1.000000

# specify filename for markers in LD to markers in .gor file
ld_markers <- paste0(input_folder, "ld-to-", markers_filename)

# get markers in LD to markers in .gor file
toWrite <- FALSE
if(file.exists(ld_markers)==FALSE){
sys.cmd <- paste0("gorpipe \"", markers_gor, " | join -snpsnp ", ld_file, "\" > ", ld_markers)
writeLines(sys.cmd)
system(sys.cmd)
toWrite <- TRUE
}

ld_tidy <- fread(ld_markers)

if(toWrite==TRUE){
  
ld_tidy <- ld_tidy[Name==mk1, ]

ld_markers_tab <- gsub(".gor",".tab",ld_markers)

write.table(ld_tidy[,.(Chrom, Pos=pos2,Name=mk2,Index=Name,rsq)][duplicated(Name)==F], 
            file=ld_markers_tab, 
            quote=FALSE, sep="\t", row.names=FALSE)

# convert to gor
sys.cmd <- paste0("gorpipe \"", ld_markers_tab, " | sort genome\" > ", ld_markers)
system(sys.cmd)

# cleanup
system(paste0("rm ", ld_markers_tab))

}

# make command lines
sys.cmd0 <- NULL

for(i in 1:nrow(annList)){
  
ann <- annList[i, ]

  if(ann$isSNP=="TRUE"){ join_parameter <- "-snpsnp" } else { join_parameter <- "-segsnp" }
  
    output_file <- paste0(output_folder, ann$dataset_id,".gor")
  
      sys.cmd <- paste0("gorpipe \" -sID ", ann[, as.character(path2gord)], " | join ", join_parameter, " ", ld_markers, "\" > ", output_file)

        sys.cmd0 <- c(sys.cmd0, sys.cmd)

  }

# print out the number of command to execute
print(paste0("Number of commands to execute: ", length(sys.cmd0)))

# write out command lines
write.table(sys.cmd0, file=paste0(rp_folder, "commands.txt"), col.names=FALSE, row.names=FALSE, sep="\t", quote=F)
