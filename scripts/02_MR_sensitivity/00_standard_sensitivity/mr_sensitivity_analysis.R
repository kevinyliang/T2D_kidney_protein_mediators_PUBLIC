#' Compute the sensitivity analyses
#' includes fstat for full MR, per snp, het, pleio
#' 
#' Assumes the same input file format -> project specific and requires parameter files
#' 
library(doMC)
source("get_harmonized_res.R")

doMC::registerDoMC(15)
compute_isq  <- function(Q,Q_df){
  I_sq  <- 100*(Q - Q_df)/Q
  return(I_sq)
}

compute_sensitivity  <- function(test_input,UKB=F){
  all_egger_res  <- foreach(each_row = 1:nrow(test_input),.combine='rbind')%dopar%{
    print(glue("{each_row}/{nrow(test_input)}"))
    if (UKB){
      hres  <- get_harmonized_res_UKB(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }else{
      hres  <- get_harmonized_res(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }
    mr_egger_res  <- mr_pleiotropy_test(hres)
    filtered_hres  <- hres %>% filter(mr_keep)
    if (nrow(mr_egger_res) == 0){
      mr_egger_res  <- data.frame(
          id.exposure = hres$id.exposure %>% unique(),
          id.outcome = hres$id.outcome %>% unique(),
          outcome = hres$outcome %>% unique(),
          exposure = hres$exposure %>% unique(),
          egger_intercept=NA,
          se = NA,
          pval=NA,
          disease_on_prot = test_input$disease_on_prot[each_row]
        )
    }else if (is.na(mr_egger_res$egger_intercept)){
      mr_egger_res  <- data.frame(
          id.exposure = hres$id.exposure %>% unique(),
          id.outcome = hres$id.outcome %>% unique(),
          outcome = hres$outcome %>% unique(),
          exposure = hres$exposure %>% unique(),
          egger_intercept=NA,
          se = NA,
          pval=NA,
          disease_on_prot = test_input$disease_on_prot[each_row]
        )
    }else{
      mr_egger_res  <- mr_egger_res %>% mutate(
          disease_on_prot = test_input$disease_on_prot[each_row]
        )
    }
    mr_egger_res
  }

  all_heterogeniety_res  <- foreach(each_row = 1:nrow(test_input),.combine='rbind')%dopar%{
    print(glue("{each_row}/{nrow(test_input)}"))
    if (UKB){
      hres  <- get_harmonized_res_UKB(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }else{
      hres  <- get_harmonized_res(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }
    het_res  <- mr_heterogeneity(hres)
    filtered_hres  <- hres %>% filter(mr_keep)
    #' more than 1 SNP gives IVW and Egger regressino ->
    #' we don't want that, just heterogeniety of wald ratios.
    mr_snp  <- mr_singlesnp(hres) %>% filter(!grepl("All",SNP))
    if (nrow(het_res)==0){
      het_res <- data.frame(
          id.exposure = hres$id.exposure %>% unique(),
          id.outcome = hres$id.outcome %>% unique(),
          outcome = hres$outcome %>% unique(),
          exposure = hres$exposure %>% unique(),
          method=NA,
          Q = NA,
          Q_df = NA,
          Q_pval = NA,
          disease_on_prot = test_input$disease_on_prot[each_row],
          I_sq_formula = NA,
          I_sq_isqfunc = NA
        )
    }else{
      het_res  <- het_res %>% mutate(
          disease_on_prot = test_input$disease_on_prot[each_row],
          I_sq_formula = max(compute_isq(Q = het_res$Q,Q_df = het_res$Q_df),0),
          I_sq_isqfunc = Isq(mr_snp$b,mr_snp$se)
        )
    }
    het_res
  }


  all_snp_fstat_results  <- foreach(each_row = 1:nrow(test_input),.combine='rbind')%dopar%{
    print(glue("{each_row}/{nrow(test_input)}"))
    if (UKB){
      hres  <- get_harmonized_res_UKB(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }else{
      hres  <- get_harmonized_res(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }
    filtered_hres  <- hres %>% filter(mr_keep)
    if (nrow(filtered_hres) != 0){
      each_snp_fstat_results  <- filtered_hres %>% 
        mutate(
          r2 = compute_var_explained(
            beta = beta.exposure,
            maf = eaf.exposure,
            se = se.exposure,
            N = samplesize.exposure
          ),
          fstat = compute_fstat(r2,samplesize.exposure,1),
          disease_on_prot = test_input$disease_on_prot[each_row]
        ) %>% select(outcome,exposure,SNP,r2,fstat,disease_on_prot)
    }else{
      each_snp_fstat_results  <- data.frame(
        outcome  = NA,
        exposure = NA,
        SNP= NA,
        r2 = NA,
        fstat = NA,
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }
    each_snp_fstat_results
  }

  all_fstat_results  <- foreach (each_row = 1:nrow(test_input),.combine='rbind')%dopar%{
    print(glue("{each_row}/{nrow(test_input)}"))
    if (UKB){
      hres  <- get_harmonized_res_UKB(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }else{
      hres  <- get_harmonized_res(
        exposure = test_input$exposure[each_row],
        outcome = test_input$outcome[each_row],
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
    }
    filtered_hres  <- hres %>% filter(mr_keep)
    if (nrow(filtered_hres) != 0){
      total_r2  <- compute_var_explained(
        beta = filtered_hres$beta.exposure,
        maf = filtered_hres$eaf.exposure,
        se = filtered_hres$se.exposure,
        N = filtered_hres$samplesize.exposure
      )
      total_fstat = compute_fstat(
        total_r2,
        max(filtered_hres$samplesize.exposure),
        nrow(filtered_hres)
      )
      each_fstat_results  <- data.frame(
          outcome = unique(filtered_hres$outcome),
          exposure = unique(filtered_hres$exposure),
          r2 = total_r2,
          fstat = total_fstat,
          disease_on_prot = test_input$disease_on_prot[each_row]
        )
    }else{
      each_fstat_results  <- data.frame(
        outcome = NA,
        exposure =NA,
        r2 = NA,
        fstat = NA,
        disease_on_prot = test_input$disease_on_prot[each_row]
      )
      
    }
    each_fstat_results
  }
  return(
    list(
      "all_fstat_results" = all_fstat_results,
      "all_snp_fstat_results" = all_snp_fstat_results,
      "all_heterogeniety_res" = all_heterogeniety_res,
      "all_egger_res" = all_egger_res
    )
  )
}
