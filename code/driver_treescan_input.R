# Top-level code for execution of data request

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splitTools))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(lubridate))


# Need to do this for assignInNamespace to work
suppressPackageStartupMessages(library(dbplyr))

# Required for execution using Rscript
suppressPackageStartupMessages(library(methods))

#' Set up the execution environment
#'
#' The .load() function sources the R files needed to execute the query
#' and sets up the execution environment.  In particular, all of the base
#' framework files, as well as files inthe code_dir with names matching
#' `cohort_*.R` or `analyze_*.R` will be sourced.
#'
#' This function is usually run automatically when the `run.R` file is sourced
#' to execute the request.  It may also be executed manually during an
#' interactive session to re-source changed code or to re-establish a connection
#' to the database.
#'
#' **N.B.** You will almost never have to edit this function.
#'
#' @param here The name of the top-level directory for the request.  The default
#'   is `config('base_dir')` if the config function has been set up, or the
#'   global variable `base_dir` if not.
#'
#' @return The value of `here`.
#' @md
.load <- function(here = ifelse(typeof(get('config')) == 'closure',
                                config('base_dir'), base_dir)) {
  source(file.path(here, 'code', 'config.R'))
  source(file.path(here, 'code', 'req_info.R'))
  source(config('site_info'))
  source(file.path(here, config('subdirs')$code_dir, 'setup.R'))
  source(file.path(here, config('subdirs')$code_dir, 'codesets.R'))
  for (fn in list.files(file.path(here, config('subdirs')$code_dir),
                        'util_.+\\.R', full.names = TRUE))
    source(fn)
  for (fn in list.files(file.path(here, config('subdirs')$code_dir),
                        'cohort_.+\\.R', full.names = TRUE))
    source(fn)
  for (fn in list.files(file.path(here, config('subdirs')$code_dir),
                        'analyze_.+\\.R', full.names = TRUE))
    source(fn)
  source(file.path(here, config('subdirs')$code_dir, 'cohorts.R'))
  
  .env_setup()
  
  for (def in c('retain_intermediates', 'results_schema')) {
    if (is.na(config(def)))
      config(def, config(paste0('default_', def)))
  }
  
  here
}

#' Execute the request
#'
#' This function presumes the environment has been set up, and executes the
#' steps of the request.
#'
#' In addition to performing queries and analyses, the execution path in this
#' function should include periodic progress messages to the user, and logging
#' of intermediate totals and timing data through [append_sum()].
#'
#' This function is also typically executed automatically, but is separated from
#' the setup done in [.load()] to facilitate direct invocation during
#' development and debugging.
#'
#' @param base_dir The name of the top-level directory for the request.  The default
#'   is `config('base_dir')`, which should always be valid after execution of
#'   [.load()].
#'
#' @return The return value is dependent on the content of the request, but is
#'   typically a structure pointing to some or all of the retrieved data or
#'   analysis results.  The value is not used by the framework itself.
#' @md
.run  <- function(base_dir = config('base_dir')) {
  
  message('Starting execution with framework version ',
          config('framework_version'))
  
  init_sum(cohort = 'Start', persons = 0)
  rslt <- list()
  
  ## parameters
  min_date="2020-11-01"
  max_date="2022-09-01" 
  min_window=-90
  max_window=-28
  n_visits=2
  
  ## Output rules for when to include B94.8 code.
  read.csv("./specs/post_infectious_disorder_rules.csv") %>%
    output_tbl("post_infectious_disorder_rules", index=c("condition_source_concept_id"))
  
  ## Output cluster master codeset
  load_codeset("cluster_master") %>% output_tbl("cluster_master", index="concept_id")
  
  ## Make cohort: COVID positive+PASC+MIS-C diagnosed patients, require at least 2 visits following index date
  cohort_pasc_or_misc_untested<-make_cohort(max_date=max_date) %>%
    filter_visits(n_visits=n_visits)
  
  cohort_pasc_or_misc_untested_demo<-cohort_pasc_or_misc_untested %>% impute_dates(min_window=min_window, max_window=max_window) %>%
  join_cohort_demo() %>% filter_dates(min_date, max_date) %>% get_cohort_entry_period()
  
  output_tbl(cohort_pasc_or_misc_untested_demo, "cohort_pasc_or_misc_untested", index="person_id")
  

  ## create 5-fold cross validation split and output cohort tables
  split_ids<-split_cohort(cohort_demo=results_tbl("cohort_pasc_or_misc_untested") %>% collect, seed=11015)
  train_id_list<-split_ids[[1]]
  test_id_list<-split_ids[[2]]
  for (i in 1:5){
    as.data.frame(train_id_list[i], col.names=c("person_id")) %>% output_tbl(paste0("cohort", i, "_train_ids", sep=""))
    as.data.frame(test_id_list[i], col.names=c("person_id")) %>% output_tbl(paste0("cohort", i, "_test_ids", sep=""))
    
    results_tbl("cohort_pasc_or_misc_untested") %>% 
      inner_join(results_tbl(paste0("cohort", i, "_train_ids", sep="")), by="person_id") %>% 
      output_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep=""), index="person_id")
    
    results_tbl("cohort_pasc_or_misc_untested") %>% 
      inner_join(results_tbl(paste0("cohort", i, "_test_ids", sep="")), by="person_id") %>% 
      output_tbl(paste("cohort_pasc_or_misc_untested_test_", i, sep=""), index="person_id")
    
    results_tbl(paste0("cohort_pasc_or_misc_untested_train_", i, sep="")) %>% collect %>% write.csv(paste0("./specs/cohort_features/cohort_untested_demo_train_", i, ".csv",sep=""), row.names=FALSE)
    results_tbl(paste0("cohort_pasc_or_misc_untested_test_", i, sep="")) %>% collect %>% write.csv(paste0("./specs/cohort_features/cohort_untested_demo_test_", i, ".csv",sep=""), row.names=FALSE)
    
  }
  
  ######### FEATURE SELECTION
  ######### Create input for TreeScan to construct features in the 
  ######## condition, drug, procedure, and lab domains
  
  
  ## TreeScan input parameters
  days_delta=-28
  follow_months=6
  
  ###### Conditions 
  snomed_tree<-make_tree(type="SNOMED", remove_postviral=TRUE)
  snomed_tree %>% output_tbl("snomed_tree", indices=c("ancestor_concept_id", "descendant_concept_id"))
  
  snomed_tree %>%
    collect %>% 
    distinct(ancestor_concept_name, ancestor_concept_id, descendant_concept_name, descendant_concept_id) %>%
    dplyr::union(data.frame(descendant_concept_id=c(441840L), descendant_concept_name="Clinical finding", ancestor_concept_id=NA, ancestor_concept_name="")) %>%
    dplyr::union(data.frame(descendant_concept_id=c(4194404L), descendant_concept_name="Disorder confirmed", ancestor_concept_id=NA, ancestor_concept_name="")) %>%
    write.csv("./specs/treescan_input/snomed_tree.csv", row.names=FALSE, na="")

  snomed_counts_all<-get_snomed_counts(tree=snomed_tree, 
                                            cohort=results_tbl("cohort_pasc_or_misc_untested"), 
                                            days_delta=days_delta,
                                            follow_months=follow_months,
                                            max_date=max_date) %>%
    ungroup %>%
    distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
    select(node, node_concept_id, pop_count, case_count, control_count) %>% 
    write.csv("./specs/treescan_input/snomed_counts_all.csv", row.names=FALSE)
  
  
  snomed_counts_list<-list()
  for (i in 1:5){
    snomed_counts_list[[i]]=get_snomed_counts(tree=snomed_tree,
                                            cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")) ,
                                            days_delta=days_delta,
                                            follow_months=follow_months, max_date=max_date)%>%
      ungroup %>%
      distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
      select(node, node_concept_id, pop_count, case_count, control_count)
    write.csv(snomed_counts_list[i], paste("./specs/treescan_input/snomed_counts_", i, ".csv", sep=""), row.names=FALSE)
  }
  
 
  ####### Labs
  loinc_tree<-make_tree(type="LOINC")
  loinc_tree %>% output_tbl("loinc_tree", indices=c("ancestor_concept_id", "descendant_concept_id"))
  
  
  
  ancestors_loinc<-loinc_tree %>% distinct(ancestor_concept_id) %>% pull(ancestor_concept_id)
  descendants_loinc<-loinc_tree %>% distinct(descendant_concept_id) %>% pull(descendant_concept_id)
  
  roots<-setdiff(ancestors_loinc, descendants_loinc)
  
  missing_nodes<-data.frame(ancestor_concept_id=NA, ancestor_concept_name="", descendant_concept_id=roots, descendant_concept_name="")
  
  loinc_tree %>%
    collect %>% 
    mutate(ancestor_concept_id=as.character(ancestor_concept_id)) %>%
    dplyr::union(missing_nodes) %>% 
    distinct(ancestor_concept_name, ancestor_concept_id, descendant_concept_name, descendant_concept_id) %>%
    write.csv("./specs/treescan_input/loinc_tree.csv", row.names=FALSE, na="")
  


  get_loinc_counts(tree=loinc_tree, 
                   cohort=results_tbl("cohort_pasc_or_misc_untested"), 
                   days_delta=days_delta,
                   follow_months=follow_months) %>%
    ungroup %>%
    distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
    write.csv("./specs/treescan_input/loinc_counts_all.csv", row.names=FALSE)
  
  for (i in 1:5){
    get_loinc_counts(tree=loinc_tree, 
                     cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")), 
                     days_delta=days_delta,
                     follow_months=follow_months) %>%
      ungroup %>%
      distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
      write.csv(paste("./specs/treescan_input/loinc_counts_", i, ".csv", sep=""), row.names=FALSE)
  }
  
  
  
  
  ##### Procedures
  
  proc_tree<-make_tree(type="proc")
  proc_tree %>% output_tbl("proc_tree", indices=c("ancestor_concept_id", "descendant_concept_id"))
  
  
  ancestors_proc<-proc_tree %>% distinct(ancestor_concept_id) %>% pull(ancestor_concept_id)
  descendants_proc<-proc_tree %>% distinct(descendant_concept_id) %>% pull(descendant_concept_id)
  
  roots<-setdiff(ancestors_proc, descendants_proc)
  
  missing_nodes<-data.frame(ancestor_concept_id=NA, ancestor_concept_name="", descendant_concept_id=roots, descendant_concept_name="")
  
  proc_tree %>%
    collect %>% 
    mutate(ancestor_concept_id=as.character(ancestor_concept_id)) %>% 
    dplyr::union(missing_nodes) %>% 
    distinct(ancestor_concept_name, ancestor_concept_id, descendant_concept_name, descendant_concept_id) %>%
    write.csv("./specs/treescan_input/proc_tree.csv", row.names=FALSE, na="")
  
  
  
  
  get_proc_counts(tree=proc_tree, 
                  cohort=results_tbl("cohort_pasc_or_misc_untested"), 
                  days_delta=days_delta,
                  follow_months=follow_months) %>%
    ungroup %>%
    distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
    write.csv("./specs/treescan_input/proc_counts_all.csv", row.names=FALSE)
  
  for (i in 1:5){
    get_proc_counts(tree=proc_tree, 
                    cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")), 
                    days_delta=days_delta,
                    follow_months=follow_months) %>%
      ungroup %>%
      distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
      write.csv(paste("./specs/treescan_input/proc_counts_", i, ".csv", sep=""), row.names=FALSE)
  }
  
  
  
  ##### Drugs
  
  drug_tree<-make_tree(type="drug")
  drug_tree %>% output_tbl("drug_tree", indices=c("ancestor_concept_id", "descendant_concept_id"))
  
  
  ancestors_drug<-drug_tree %>% distinct(ancestor_concept_id) %>% pull(ancestor_concept_id)
  descendants_drug<-drug_tree %>% distinct(descendant_concept_id) %>% pull(descendant_concept_id)
  
  roots<-setdiff(ancestors_drug, descendants_drug)
  
  missing_nodes<-data.frame(ancestor_concept_id=NA, ancestor_concept_name="", descendant_concept_id=roots, descendant_concept_name="")
  
  drug_tree %>%
    collect %>% 
    mutate(ancestor_concept_id=as.character(ancestor_concept_id)) %>% 
    dplyr::union(missing_nodes) %>% 
    distinct(ancestor_concept_name, ancestor_concept_id, descendant_concept_name, descendant_concept_id) %>%
    write.csv("./specs/treescan_input/drug_tree.csv", row.names=FALSE, na="")
  
  
  
  
  
  get_drug_counts(tree=drug_tree, 
                  cohort=results_tbl("cohort_pasc_or_misc_untested"), 
                  days_delta=days_delta,
                  follow_months=follow_months) %>%
    ungroup %>%
    distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
    write.csv("./specs/treescan_input/drug_counts_all.csv", row.names=FALSE)
  
  for (i in 1:5){
    get_drug_counts(tree=drug_tree, 
                    cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")), 
                    days_delta=days_delta,
                    follow_months=follow_months) %>%
      ungroup %>%
      distinct(node, node_concept_id, pop_count, case_count, control_count) %>%
      write.csv(paste("./specs/treescan_input/drug_counts_", i, ".csv", sep=""), row.names=FALSE)
  }
  
  
  
  
  #### Treescan analysis parameters
  print_split_probs(max_date=max_date, follow_months=follow_months)
  
  
  
  
  message('Done.')
  
  
}

#' Set up and execute a data request
#'
#' This function encapsulates a "production" run of the data request.  It sets
#' up the environment, executes the request, and cleans up the environment.
#'
#' Typically, the `run.R` file calls run_request() when in a production mode.
#'
#' @param base_dir Path to the top of the data request files.  This is
#'   typically specified in `run.R`.
#'
#' @return The result of [.run()].
#' @md
run_request <- function(base_dir) {
  base_dir <- .load(base_dir)
  on.exit(.env_cleanup())
  .run(base_dir)
}

  