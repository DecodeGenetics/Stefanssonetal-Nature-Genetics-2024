library(Rsolnp)
library(ggplot2)
library(gridExtra)

#compute the likelihood for a given enrichment vector E
log_likelihood_MAF <- function(E,idx,signal_mat,q_mat,MAF_bins){
  loglik <- 0
  for(m in MAF_bins){
    loglik <- loglik + sum(log((signal_mat[idx[[m]],]) %*% (E*q_mat[m,]/sum(E*q_mat[m,]))))
  }
  return(loglik)
}

#Run maximum liklehiood estimation of the enrichment parameter E
get_E_rmod <- function(E_init,idx,anno_vec,MAF_bins,signal_mat,q_mat,ref_anno_idx){
  equality_constraint <- function(theta){theta[ref_anno_idx]}
  MAF_bins <- names(idx)
  #This is the optimization step
  optim_obj <- solnp(E_init,
                     fun=function(x){-log_likelihood_MAF(x,idx,signal_mat,q_mat,MAF_bins)},
                     eqfun=equality_constraint,
                     eqB=1,
                     LB=rep(0,length(E_init)),
                     control=list(trace=0))
  optim_dat <- tibble(Annotation=anno_vec,enr=optim_obj$pars,conv=optim_obj$convergence,nfuneval=optim_obj$nfuneval)
  return(optim_dat)
}

#Run the enrichment model iteratively using bootstrapping to compute enrichment confidence intervals
get_bootstrap_rmod <- function(signal_mat,q_mat,anno_vec,MAF_bins,signal_summary_dat,N,ref_anno_idx){
  E_init <- rep(1,length(anno_vec))
  signal_summary_unique <- distinct(signal_summary_dat,lead_name,MAF_bin_lead)
  signal_idx <- lapply(signal_summary_unique$lead_name,function(s){
    which(signal_summary_dat$lead_name==s)
  })
  opt_dat <- lapply(1:N,function(i){
    bootstrap_idx <- lapply(MAF_bins,function(m){
      uniq_idx <- sample(which(signal_summary_unique$MAF_bin_lead==m),size=sum(signal_summary_unique$MAF_bin_lead==m),replace=T)
      return(unlist(sapply(uniq_idx,function(j) signal_idx[[j]])))
    })
    names(bootstrap_idx) <- MAF_bins
    get_E_rmod(E_init,bootstrap_idx,anno_vec,MAF_bins,signal_mat,q_mat,ref_anno_idx) %>% 
      mutate(iter=i)
  }) %>% bind_rows()
  return(opt_dat)
}

#Summarize enrichment point estimates for each annotation and extract confidence intervals from bootstrapping
get_enr_summary <- function(bootstrap_dat,point_estimate_dat,freq_summary_dat){
  if('MAF_bin' %in% names(bootstrap_dat)){
    enr_summary <- inner_join(bootstrap_dat,freq_summary_dat,by=c('Annotation','MAF_bin')) %>%
      group_by(iter,Annotation) %>%
      summarise(enr=sum(p)/sum(q_MAF),.groups='drop')
  }else if('p' %in% names(bootstrap_dat)){
    enr_summary <- inner_join(bootstrap_dat,distinct(freq_summary_dat,Annotation,q),by='Annotation') %>%
      mutate(enr=p/q)
  }else{
    enr_summary <- bootstrap_dat
  }
  enr_summary <- inner_join(enr_summary,select(point_estimate_dat,Annotation,enr),by=c('Annotation'),suffix=c('','_hat')) %>%
    mutate(log_enr=log(enr)) %>%
    group_by(Annotation) %>%
    summarise(enr_lower=quantile(enr,0.025,na.rm=T),
              enr_upper=quantile(enr,0.975,na.rm=T),
              enr_hat=mean(enr_hat),
              log_enr_hat=log(enr_hat),
              boot_sd=sd(log_enr),
              pval_emp=if(log_enr_hat<=0) sum(log_enr>0)/n() else sum(log_enr<=0)/n(),
              z_stat=log_enr_hat/boot_sd,
              pval_z=2*pnorm(abs(z_stat),lower.tail=F),
              .groups='drop') %>%
    arrange(Annotation) %>%
    mutate(pval_threshold=0.05/n()) %>%
    select(Annotation,enr_hat,enr_lower,enr_upper,pval_emp,pval_z,pval_threshold)
  return(enr_summary)
}

### Enrichment model visualization
get_enr_plot <- function(summary_dat){
  ggplot(summary_dat) +
    geom_point(aes(Annotation,enr_hat)) +
    geom_errorbar(aes(Annotation,ymin=enr_lower,ymax=enr_upper),width=0.2) +
    geom_hline(aes(yintercept=1),col='red') +
    theme_bw() +
    theme(axis.text.x = element_text(size=12,angle = 45,hjust=1),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14)) +
    theme(axis.text.x = element_text(angle = 45,hjust=1),
          plot.margin = unit(c(0.5,0.5,0.5,0.5),'cm')) +
    xlab('') +
    ylab('Enrichment')  
}

