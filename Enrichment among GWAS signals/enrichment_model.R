#### ---- Package loading and argument parsing ####
library(argparser)
library(dplyr)
library(tidyr)
library(data.table)
#Set the working directory to the folder containing utils.R
source('utils.R')

#### ---- Input arguments ####
p <- arg_parser("Run enrichment model presented in Sveinbjornsson et.al. (2015)","enrichment_optim.R")
p <- add_argument(p,"--anno_path",nargs=1,type='character',
                  help="Path to a folder with annotation files. Folder must include the files maf-distrib-of_variant-annotations.tsv and signals-info.tab")
p <- add_argument(p,"--N",nargs=1,type='numeric',default=1000,
                  help="Number of bootstrap iterations to run")
p <- add_argument(p,"--out_path",nargs=1,type='character',
                  help="Path to which output is to be written")

args <- parse_args(p)

#### ---- Uploading and processing of input data ---- 
# signal-info.tab has the following columns
# - signal: ID of the association signal in the format {variant_id of lead variant of the signal}-{phenotype_id}
# - Phenotype: Phenotype ID
# - Name: Variant ID
# - inLD_pos: Position of variant in column Name
# - inLD_X2: X2 association statistic of the variant in Name column with the phenotype
# - Annotation: Annotation of the variant in name
# - lead_name: Variant ID of the lead variant in the signal
# - MAF_bin_lead: Minor allele frequency cateogory of the lead

signal_LD_dat <- fread(paste0(args$anno_path,'/signals-info.tab')) %>% 
                 distinct() %>%
                 group_by(signal) %>%
                 dplyr::mutate(Annotation_lead=Annotation[Name==lead_name]) %>%
                 ungroup()

signalinfo <- signal_LD_dat[signal_LD_dat$lead_name==signal_LD_dat$Name,1:2]

# maf-distrib-of-variant-annotations.tsv is located in a subfolder under args$anno_path called summary. 
# This table is a matrix containing the number of variants for each annotation in each MAF bin category (MAF_bin categories and annotations must match those in the signals-info.tab).
# The first column is called Max_Consequence and the remaining column names match each of the MAF_bin category names

MAF_freq_dat <- fread(paste0(args$anno_path,'/summary/maf-distrib-of-variant-annotations.tsv')) %>%
                as_tibble() %>%
                dplyr::rename(Annotation=Max_Consequence) %>% 
                mutate(tot_Freq=rowSums(across(matches('MAF')))) %>% 
                select(-tot_Freq) %>%
                pivot_longer(cols=matches('MAF'),names_to='MAF_bin',values_to='Freq') %>%
                mutate(Freq=if_else(is.na(Freq),0L,Freq))
#summarise the number of variants in each MAF bin category among lead variants in association signals              
signal_freq_dat <- distinct(signal_LD_dat,signal,MAF_bin_lead) %>%
                   count(MAF_bin_lead) %>%
                   dplyr::rename(Freq=n)

MAF_bins <- sort(unique(MAF_freq_dat$MAF_bin))
Annotation_sorted <- sort(unique(MAF_freq_dat$Annotation))
freq_summary_dat <- inner_join(MAF_freq_dat,signal_freq_dat,by=c('MAF_bin'='MAF_bin_lead'),suffix=c('_MAF','_signal')) %>%
                    mutate(Freq_MAF=as.double(Freq_MAF),
                           Freq_signal=as.double(Freq_signal)) %>%
                    group_by(MAF_bin) %>%
                    mutate(Freq_MAF_bin=sum(Freq_MAF),
                           q_MAF=Freq_MAF/Freq_MAF_bin) %>%
                    group_by(Annotation) %>%
                    mutate(q=sum(Freq_MAF)/sum(Freq_MAF_bin)) %>%
                    ungroup() %>%
                    mutate(MAF_bin=factor(MAF_bin,levels=MAF_bins),
                           Annotation=factor(Annotation,levels=Annotation_sorted)) %>%
                    arrange(MAF_bin,Annotation)

signal_summary_default <- crossing(signal=unique(signal_LD_dat$signal),Annotation=Annotation_sorted) %>%
                          inner_join(distinct(signal_LD_dat,signal,lead_name,MAF_bin_lead,Annotation_lead),by='signal')
#precompute the product of the q's using the logsumexp trick before fitting the model
signal_summary_dat <- inner_join(signal_LD_dat,select(freq_summary_dat,Annotation,MAF_bin,q_MAF,q),by=c('Annotation','MAF_bin_lead'='MAF_bin')) %>%
                      group_by(signal) %>%
                      mutate(log_q_MAF_prod=sapply(1:n(),function(i) sum(log(q_MAF[-i]))),
                             q_MAF_prod_adj=exp(log_q_MAF_prod-min(log_q_MAF_prod)),
                             weight_MAF=exp(inLD_X2-max(inLD_X2))*q_MAF_prod_adj) %>%
                      ungroup() %>%
                      group_by(signal,lead_name,MAF_bin_lead,Annotation_lead,Annotation) %>%
                      summarise(weight_MAF_tot=sum(weight_MAF),
                                topMarkers=paste(head(Name[order(-inLD_X2)]),collapse=',')) %>%
                      left_join(signal_summary_default,.,by=c('signal','lead_name','MAF_bin_lead','Annotation_lead','Annotation')) %>%
                      mutate(weight_MAF_tot=if_else(is.na(weight_MAF_tot),0,weight_MAF_tot)) %>%
                      pivot_wider(id_cols=c('signal','lead_name','MAF_bin_lead','Annotation_lead'),names_from='Annotation',values_from=c('weight_MAF_tot','topMarkers')) %>%
                      mutate(MAF_bin_lead=factor(MAF_bin_lead,levels=MAF_bins),
                             Annotation_lead=factor(Annotation_lead,levels=Annotation_sorted))

signal_MAF_mat <- as.matrix(select(signal_summary_dat,matches('weight_MAF_tot')))
q_mat <- matrix(freq_summary_dat$q_MAF,nrow=length(MAF_bins),byrow = T)
rownames(q_mat) <- MAF_bins

#### ---- Optimization and post-processing -----
if(!file.exists(args$out_path)){
  dir.create(args$out_path)
}
# Sveinbjornsson revised model
idx=lapply(MAF_bins,function(m) which(signal_summary_dat$MAF_bin_lead==m))
names(idx) <- MAF_bins
ref_anno_idx <- distinct(freq_summary_dat,Annotation,q) %>% arrange(Annotation) %>% summarise(ref_anno_idx=which.max(q)) %>% .$ref_anno_idx
point_estimate_rmod <- get_E_rmod(E_init=rep(1,length(Annotation_sorted)), idx=idx, anno_vec=Annotation_sorted, MAF_bins=MAF_bins, signal_mat=signal_MAF_mat, q_mat=q_mat, ref_anno_idx=ref_anno_idx)
if(args$N==1){
  bootstrap_rmod <- select(point_estimate_rmod,Annotation) %>% mutate(enr=as.numeric(NA),conv=0,nfuneval=0,iter=1)
}else{
  bootstrap_rmod <- get_bootstrap_rmod(signal_mat=signal_MAF_mat,q_mat=q_mat,anno_vec=Annotation_sorted,MAF_bins=MAF_bins,signal_summary_dat,N=args$N,ref_anno_idx=ref_anno_idx)
}
enr_summary_rmod <- get_enr_summary(bootstrap_rmod,point_estimate_rmod,freq_summary_dat)
write.table(bootstrap_rmod,file = paste0(args$out_path,'/bootstrap_mod_MAF_adjusted.tsv'),sep='\t',quote=F,row.names=F)
write.table(enr_summary_rmod,file = paste0(args$out_path,'/enr_summary_mod_MAF_adjusted.tsv'),sep='\t',quote=F,row.names=F)


#### ---- Plot results ------
enr_plot_rmod <- get_enr_plot(enr_summary_rmod)
ggsave(enr_plot_rmod,filename = paste0(args$out_path,'/enr_mod_MAF_adjusted.png'),width=21,height=10)
#Plot the figure with y axis on the log scale
lim_log_vec <- c(floor(log10(min(enr_summary_rmod$enr_lower))),ceiling(log10(max(enr_summary_rmod$enr_upper))))
ggsave(enr_plot_rmod + scale_y_continuous(trans = scales::log10_trans(),
                                          breaks = 10^(seq(lim_log_vec[1],lim_log_vec[2])),
                                          labels = function(x) parse(text=paste0('10^',log10(x)))),filename = paste0(args$out_path,'/enr_mod_MAF_adjusted_logy.png'),width=21,height=10)
