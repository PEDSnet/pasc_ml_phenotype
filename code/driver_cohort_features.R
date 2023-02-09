# Top-level code for execution of data request

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splitTools))
suppressPackageStartupMessages(library(VennDiagram))


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
.run_features  <- function(base_dir = config('base_dir')) {
  
  message('Starting execution with framework version ',
          config('framework_version'))
  
  init_sum(cohort = 'Start', persons = 0)
  rslt <- list()
  
  ## parameters
  max_date="2022-09-01"
  days_delta=-28
  follow_months=6
  split_dates=TRUE
  days_min=27
  days_max=180
  
  if (split_dates){
    str_append=""
  }else{
    str_append="_no_split"
  }
  
  
  COND_TREE_LEVEL<-3
  DRUG_TREE_LEVEL<-2
  # LAB_TREE_LEVEL<-2
  LAB_TREE_LEVEL<-3
  # PROC_TREE_LEVEL<-1
  PROC_TREE_LEVEL<-3
  
  # COND_P_VAL<-.05
  COND_P_VAL<-.01
  # DRUG_P_VAL<-.05
  DRUG_P_VAL<-.01
  # LAB_P_VAL<-.05
  LAB_P_VAL<-.01
  # PROC_P_VAL<-.05
  PROC_P_VAL<-.01
  
  COND_N_OBS<-1600
  # DRUG_N_OBS<-200
  DRUG_N_OBS<-500
  # LAB_N_OBS<-1000
  LAB_N_OBS<-10000
  # PROC_N_OBS<-50
  PROC_N_OBS<-100
  
  ## Read in trees
  snomed_tree<-results_tbl("snomed_tree")
  loinc_tree<-results_tbl("loinc_tree")
  proc_tree<-results_tbl("proc_tree")
  drug_tree<-results_tbl("drug_tree")
  
  
  ######## Construct cohort features for each domain
  ########
  ########
  
  
  ### Util counts
  
  
  cohort_untested_util<- results_tbl("cohort_pasc_or_misc_untested")  %>%
    get_util(days_delta=days_delta)
  
  cohort_untested_util_wide<-cohort_untested_util %>%
    partition_dates_and_pivot_wide_util(split_dates=split_dates, days_min=days_min, days_max=days_max)
  
  cohort_untested_util_final<-results_tbl("cohort_pasc_or_misc_untested") %>%  
    distinct(person_id, index_date, pasc_flag) %>%  
    full_join(cohort_untested_util_wide %>% select(-observation_type), 
              by=c("person_id", "index_date", "pasc_flag"))
  
  myList <- setNames(lapply(vector("list", ncol(cohort_untested_util_final)-3), function(x) x <- 0), setdiff(colnames(cohort_untested_util_final), c("person_id", "index_date", "pasc_flag")))
  
  cohort_untested_util_final_b<-cohort_untested_util_final %>%
    replace_na(myList) %>%compute_new(index="person_id")
  
  cohort_untested_util_final_b %>%
    output_tbl("cohort_untested_util_features_all", index="person_id")
  
#  cohort_untested_util_final[,4:ncol(cohort_untested_util_final)]=cohort_untested_util_final[,4:ncol(cohort_untested_util_final)] %>% replace(is.na(.), 0)
#  cohort_untested_util_final %>% write.csv(paste("./specs/cohort_features/cohort_untested_util_all", str_append, ".csv", sep=""), row.names=FALSE)
#  cohort_util_temp<-cohort_untested_util_final %>% 
#    output_tbl(paste("cohort_util_temp", str_append, sep=""), temp=TRUE, index="person_id") 
#  rm(cohort_untested_util_final)
  
  for (i in 1:5){
    results_tbl("cohort_untested_util_features_all") %>% 
      inner_join(results_tbl(paste("cohort", i, "_train_ids", sep="")), by="person_id") %>%
      output_tbl(paste("cohort_util_features_train_", i, sep=""), index="person_id")
    results_tbl("cohort_untested_util_features_all") %>% 
      inner_join(results_tbl(paste("cohort", i, "_test_ids", sep="")), by="person_id") %>%
      output_tbl(paste("cohort_util_features_test_", i, sep=""), index="person_id")
  }

  
  
  #  Lab counts
  
  cohort_untested_lab_cts<- results_tbl("cohort_pasc_or_misc_untested") %>%
    get_lab_cts(days_delta=days_delta)
  
  cohort_untested_lab_cts_wide<-cohort_untested_lab_cts %>%
    partition_dates_and_pivot_wide_lab_cts(split_dates=split_dates, days_min=days_min, days_max=days_max)
  
  cohort_untested_lab_cts_final<-results_tbl("cohort_pasc_or_misc_untested") %>%  
    distinct(person_id, index_date, pasc_flag) %>%  
    full_join(cohort_untested_lab_cts_wide %>% select(-observation_type), 
              by=c("person_id", "index_date", "pasc_flag"))# %>% collect
  
  myList <- setNames(lapply(vector("list", ncol(cohort_untested_lab_cts_final)-3), function(x) x <- 0), setdiff(colnames(cohort_untested_lab_cts_final), c("person_id", "index_date", "pasc_flag")))
  
  cohort_untested_lab_cts_final_b<-cohort_untested_lab_cts_final %>%
    replace_na(myList) %>%compute_new(index="person_id")
  
  cohort_untested_lab_cts_final_b %>%
    output_tbl("cohort_untested_lab_cts_features_all", index="person_id")
  
  for (i in 1:5){
    results_tbl("cohort_untested_lab_cts_features_all") %>% 
      inner_join(results_tbl(paste("cohort", i, "_train_ids", sep="")), by="person_id") %>%
      output_tbl(paste("cohort_lab_cts_train_", i,  sep=""), index="person_id")
    results_tbl("cohort_untested_lab_cts_features_all") %>% 
      inner_join(results_tbl(paste("cohort", i, "_test_ids", sep="")), by="person_id") %>%
      output_tbl(paste("cohort_lab_cts_test_", i,  sep=""), index="person_id")
  }
  
  
  

  
#  cohort_untested_lab_cts<-results_tbl("cohort_pasc_or_misc_untested") %>%
#    get_lab_cts(days_delta=days_delta)
#  
#  cohort_untested_lab_cts_wide<-cohort_untested_lab_cts %>%
#    partition_dates_and_pivot_wide_lab_cts()
#  
#  
#  cohort_untested_lab_cts_final<-results_tbl("cohort_pasc_or_misc_untested") %>%  
#    distinct(person_id, index_date, pasc_flag) %>%  
#    full_join(cohort_untested_lab_cts_wide %>% select(-observation_type), 
#              by=c("person_id", "index_date", "pasc_flag")) %>% collect
#  
#  cohort_untested_lab_cts_final[,4:ncol(cohort_untested_lab_cts_final)]=cohort_untested_lab_cts_final[,4:ncol(cohort_untested_lab_cts_final)] %>% replace(is.na(.), 0)
#  cohort_untested_lab_cts_final %>% write.csv("./specs/cohort_features/cohort_untested_lab_cts_all.csv", row.names=FALSE)
#  cohort_lab_cts_temp<-cohort_untested_lab_cts_final %>% 
#    output_tbl("cohort_lab_cts_temp", temp=TRUE, index="person_id") 
#  rm(cohort_untested_lab_cts_final)
#  
#  for (i in 1:5){
#    cohort_lab_cts_temp %>% 
#      inner_join(results_tbl(paste("cohort", i, "_train_ids", sep="")), by="person_id") %>%
#      collect %>%
#      write.csv(paste("./specs/cohort_features/cohort_lab_cts_train_", i, ".csv", sep=""), row.names=FALSE)
#    cohort_lab_cts_temp %>% 
#      inner_join(results_tbl(paste("cohort", i, "_test_ids", sep="")), by="person_id") %>%
#      collect %>%
#      write.csv(paste("./specs/cohort_features/cohort_lab_cts_test_", i, ".csv", sep=""), row.names=FALSE)
#  }
  
  
  
  
  
  
  
  
  
  
  
  ##### Conditions
  
  all_descendants_snomed<-snomed_tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors_snomed<-snomed_tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  all_nodes_snomed<-dplyr::union(all_descendants_snomed, all_ancestors_snomed) %>% distinct() %>% output_tbl("all_nodes_snomed", index="node_concept_id", temp=TRUE) 
  
  cuts_all<-read.csv("./specs/treescan_output/treescan_conditions_all.csv") %>% output_tbl("treescan_condition_all", temp=TRUE) %>% 
    left_join(all_nodes_snomed, by=c("Node.Identifier"="node_concept_id")) %>%
    select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
    filter(!is.na(P.value), Tree.Level>COND_TREE_LEVEL, Observations>COND_N_OBS, P.value<COND_P_VAL) %>%
    arrange(desc(Log.Likelihood.Ratio)) #%>% head(500)
  
  #cuts_all %>% write.csv("./specs/treescan_output/condition_cuts_all.csv", row.names=FALSE)
  make_features(cuts_all, all_nodes_snomed, type="SNOMED") %>%
    output_tbl("condition_features_untested_all", indices=c("ancestor_concept_id", "descendant_concept_id))"))
  #results_tbl("condition_features_untested_all") %>% collect %>% write.csv("./specs/condition_features_untested_all.csv", row.names=FALSE)
  
  condition_cut_list<-list()
  for (i in 1:5){
    condition_cut_list[[i]]<-read.csv(paste("./specs/treescan_output/treescan_conditions_", i, ".csv", sep="")) %>% output_tbl(paste("treescan_condition_", i, sep=""), temp=TRUE) %>% 
      left_join(all_nodes_snomed, by=c("Node.Identifier"="node_concept_id")) %>%
      select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
      filter(!is.na(P.value), Tree.Level>COND_TREE_LEVEL, Observations>COND_N_OBS, P.value<COND_P_VAL) %>%
      arrange(desc(Log.Likelihood.Ratio)) #%>% head(500)
    
    make_features(condition_cut_list[[i]], all_nodes_snomed, type="SNOMED") %>%
      output_tbl(paste("condition_features_untested_", i, sep=""), indices=c("ancestor_concept_id", "descendant_concept_id"))
  }
  
  
  venn.diagram(
    x = list(condition_cut_list[[1]] %>% pull(Node.Identifier), 
             condition_cut_list[[2]] %>% pull(Node.Identifier), 
             condition_cut_list[[3]] %>% pull(Node.Identifier), 
             condition_cut_list[[4]] %>% pull(Node.Identifier),
             condition_cut_list[[5]] %>% pull(Node.Identifier)),
    category.names = c("Fold 1" , "Fold 2 " , "Fold 3", "Fold 4", "Fold 5"),
    filename = './results/conditions_venn_diagram.png',
    output=TRUE
  )
  
  #all
  
  condition_features_untested_all=results_tbl("condition_features_untested_all") %>% compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
  a=Sys.time()
  cohort_untested_cond_features_all<-make_feature_tbl(cohort=results_tbl("cohort_pasc_or_misc_untested"),
                                                      type="cond", features=condition_features_untested_all,
                                                      days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
  
  output_tbl(cohort_untested_cond_features_all, "cohort_cond_features_all", index="person_id")
  print(Sys.time()-a)
#  write.csv(cohort_untested_cond_features_all, "./specs/cohort_features", str_append, "/cohort_cond_features_all.csv", row.names=FALSE)
#  rm(cohort_untested_cond_features_all)
  
  ##by split
  
  #  condition_feature_list<-list()
  #  cohort_cond_feature_train_list<-list()
  #  cohort_cond_feature_test_list<-list()
  
  
  for (i in 1:5){
    print(i)
    condition_features=results_tbl(paste("condition_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_cond_feature_train<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")),
                                                type="cond", features=condition_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_cond_feature_train %>% output_tbl(paste("cohort_cond_features_train_", i, sep=""), index="person_id")
    #write.csv(cohort_cond_feature_train, paste("./specs/cohort_features",str_append,"/cohort_cond_features_train_", i,".csv", sep=""), row.names=FALSE)
    #rm(cohort_cond_feature_train)
    #rm(cohort_untested_cond_features_train_1)
    print(i)
    condition_features=results_tbl(paste("condition_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_cond_feature_test<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_test_", i, sep="")),
                                                type="cond", features=condition_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_cond_feature_test %>% output_tbl(paste("cohort_cond_features_test_", i, sep=""), index="person_id")
    #write.csv(cohort_cond_feature_test, paste("./specs/cohort_features",str_append,"/cohort_cond_features_test_", i,".csv", sep=""), row.names=FALSE)
    #rm(cohort_cond_feature_test)
    #rm(cohort_untested_cond_features_test_1)
  }
  
  
  
  
  
  
  #### Labs
  all_descendants_loinc<-loinc_tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors_loinc<-loinc_tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  all_nodes_loinc<-dplyr::union(all_descendants_loinc, all_ancestors_loinc) %>% distinct() 
  
  
  cuts_all<-read.csv("./specs/treescan_output/treescan_labs_all.csv") %>% output_tbl("treescan_loinc_all", temp=TRUE) %>% 
    left_join(all_nodes_loinc, by=c("Node.Identifier"="node_concept_id")) %>%
    select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
    filter(!is.na(P.value), Tree.Level>=LAB_TREE_LEVEL, Observations>LAB_N_OBS, P.value<LAB_P_VAL) %>%
    arrange(desc(Log.Likelihood.Ratio)) #%>% head(500)
  #cuts_1 %>% output_tbl("treescan_lab_cuts_1", index="Node.Identifier")
  
  
#  cuts_all %>% write.csv("./specs/treescan_output/lab_cuts_all.csv", row.names=FALSE)
  make_features(cuts_all, all_nodes_loinc, type="LOINC") %>%
    output_tbl("lab_features_untested_all", indices=c("ancestor_concept_id", "descendant_concept_id"))
  
  
  
  lab_cut_list<-list()
  for (i in 1:5){
    lab_cut_list[[i]]<-read.csv(paste("./results/treescan_output/treescan_loinc_", i, ".csv", sep="")) %>% output_tbl(paste("treescan_loinc_", i, sep=""), temp=TRUE) %>% 
      left_join(all_nodes_loinc, by=c("Node.Identifier"="node_concept_id")) %>%
      select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
      filter(!is.na(P.value), Tree.Level>=LAB_TREE_LEVEL, Observations>LAB_N_OBS, P.value<LAB_P_VAL) %>%
      arrange(desc(Log.Likelihood.Ratio))
    
    make_features(lab_cut_list[[i]], all_nodes_loinc, type="LOINC") %>%
      output_tbl(paste("lab_features_untested_", i, sep=""), indices=c("ancestor_concept_id", "descendant_concept_id))"))
    
  }
  
  venn.diagram(
    x = list(lab_cut_list[[1]] %>% pull(Node.Identifier), 
             lab_cut_list[[2]] %>% pull(Node.Identifier), 
             lab_cut_list[[3]] %>% pull(Node.Identifier), 
             lab_cut_list[[4]] %>% pull(Node.Identifier),
             lab_cut_list[[5]] %>% pull(Node.Identifier)),
    category.names = c("Fold 1" , "Fold 2 " , "Fold 3", "Fold 4", "Fold 5"),
    filename = './results/labs_venn_diagram.png',
    output=TRUE
  )
  
  
  
  lab_features_untested_all=results_tbl("lab_features_untested_all") %>% compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
  cohort_untested_lab_features_all<-make_feature_tbl(cohort=results_tbl("cohort_pasc_or_misc_untested"),
                                                     type="lab", features=lab_features_untested_all, days_delta=days_delta,
                                                     split_dates=split_dates, days_min=days_min, days_max=days_max)
  
  cohort_untested_lab_features_all %>%
    output_tbl("cohort_untested_lab_features_all", index="person_id")
  
#  write.csv(cohort_untested_lab_features_all, "./specs/cohort_features",str_append,"/cohort_lab_features_all.csv", row.names=FALSE)
#  rm(cohort_untested_lab_features_all)
  
  #  lab_feature_list<-list()
  #  cohort_lab_feature_train_list<-list()
  #  cohort_lab_feature_test_list<-list()
  for (i in 1:5){
    print(i)
    lab_features=results_tbl(paste("lab_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_lab_feature_train<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")),
                                                type="lab", features=lab_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_lab_feature_train %>% output_tbl(paste("cohort_lab_features_train_", i, sep=""), index="person_id")

    print(i)
    lab_features=results_tbl(paste("lab_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_lab_feature_test<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_test_", i, sep="")),
                                               type="lab", features=lab_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_lab_feature_test %>% output_tbl(paste("cohort_lab_features_test_", i, sep=""), index="person_id")
  }
  
  
  
  
  #### Procedures   
  
  all_descendants_proc<-proc_tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors_proc<-proc_tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  all_nodes_proc<-dplyr::union(all_descendants_proc, all_ancestors_proc) %>% distinct() 
  
  
  cuts_all<-read.csv("./specs/treescan_output/treescan_procs_all.csv") %>% output_tbl("treescan_proc_all", temp=TRUE) %>% 
    left_join(all_nodes_proc, by=c("Node.Identifier"="node_concept_id")) %>%
    select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
    filter(!is.na(P.value), Tree.Level>=PROC_TREE_LEVEL, Observations>PROC_N_OBS, P.value<PROC_P_VAL) %>%
    arrange(desc(Log.Likelihood.Ratio)) #%>% head(500)
  #cuts_1 %>% output_tbl("treescan_proc_cuts_1", index="Node.Identifier")
  
#  cuts_all %>% write.csv("./specs/treescan_output/proc_cuts_all.csv", row.names=FALSE)
  make_features(cuts_all, all_nodes_proc, type="proc") %>%
    output_tbl("proc_features_untested_all", indices=c("ancestor_concept_id", "descendant_concept_id"))
  
  
  
  proc_cut_list<-list()
  for (i in 1:5){
    proc_cut_list[[i]]<-read.csv(paste("./results/treescan_output/treescan_proc_", i, ".csv", sep="")) %>% output_tbl(paste("treescan_proc_", i, sep=""), temp=TRUE) %>% 
      left_join(all_nodes_proc, by=c("Node.Identifier"="node_concept_id")) %>%
      select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
      filter(!is.na(P.value), Tree.Level>=PROC_TREE_LEVEL, Observations>PROC_N_OBS, P.value<PROC_P_VAL) %>%
      arrange(desc(Log.Likelihood.Ratio))
    
    make_features(proc_cut_list[[i]], all_nodes_proc, type="proc") %>%
      output_tbl(paste("proc_features_untested_", i, sep=""), indices=c("ancestor_concept_id", "descendant_concept_id))"))
    
  }
  
  require(VennDiagram)
  venn.diagram(
    x = list(proc_cut_list[[1]] %>% pull(Node.Identifier), 
             proc_cut_list[[2]] %>% pull(Node.Identifier), 
             proc_cut_list[[3]] %>% pull(Node.Identifier), 
             proc_cut_list[[4]] %>% pull(Node.Identifier),
             proc_cut_list[[5]] %>% pull(Node.Identifier)),
    category.names = c("Fold 1" , "Fold 2 " , "Fold 3", "Fold 4", "Fold 5"),
    filename = './results/proc_venn_diagram.png',
    output=TRUE
  )
  
  
  
  proc_features_untested_all=results_tbl("proc_features_untested_all") %>% compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
  cohort_untested_proc_features_all<-make_feature_tbl(cohort=results_tbl("cohort_pasc_or_misc_untested"),
                                                      type="proc", features=proc_features_untested_all, days_delta=days_delta,
                                                      split_dates=split_dates, days_min=days_min, days_max=days_max)
  cohort_untested_proc_features_all %>% 
    output_tbl("cohort_untested_proc_features_all", index="person_id")
  
#  write.csv(cohort_untested_proc_features_all, "./specs/cohort_features",str_append,"/cohort_proc_features_all.csv", row.names=FALSE)
#  rm(cohort_untested_proc_features_all)
  
  #  proc_feature_list<-list()
  #  cohort_proc_feature_train_list<-list()
  #  cohort_proc_feature_test_list<-list()
  
  
  
  for (i in 1:5){
    print(i)
    proc_features=results_tbl(paste("proc_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_proc_feature_train<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")),
                                               type="proc", features=proc_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_proc_feature_train %>% output_tbl(paste("cohort_proc_features_train_", i, sep=""), index="person_id")
    
    print(i)
    proc_features=results_tbl(paste("proc_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_proc_feature_test<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_test_", i, sep="")),
                                              type="proc", features=proc_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_proc_feature_test %>% output_tbl(paste("cohort_proc_features_test_", i, sep=""), index="person_id")
  }
  
  
  
  
  #### Drugs
  
  
  all_descendants_drug<-drug_tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors_drug<-drug_tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  all_nodes_drug<-dplyr::union(all_descendants_drug, all_ancestors_drug) %>% distinct() 
  
  cuts_all<-read.csv("./specs/treescan_output/treescan_drugs_all.csv") %>% output_tbl("treescan_drug_all", temp=TRUE) %>% 
    left_join(all_nodes_drug, by=c("Node.Identifier"="node_concept_id")) %>%
    select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
    filter(!is.na(P.value), Tree.Level>=2, Observations>200, P.value<0.05) %>%
    arrange(desc(Log.Likelihood.Ratio)) #%>% head(500)
  #cuts_1 %>% output_tbl("treescan_drug_cuts_1", index="Node.Identifier")
  
#  cuts_all %>% write.csv("./specs/treescan_output/drug_cuts_all.csv", row.names=FALSE)
  make_features(cuts_all, all_nodes_drug, type="drug") %>%
    output_tbl("drug_features_untested_all", indices=c("ancestor_concept_id", "descendant_concept_id"))
  
  
  
  drug_cut_list<-list()
  for (i in 1:5){
    drug_cut_list[[i]]<-read.csv(paste("./results/treescan_output/treescan_drug_", i, ".csv", sep="")) %>% output_tbl(paste("treescan_drug_", i, sep=""), temp=TRUE) %>% 
      left_join(all_nodes_drug, by=c("Node.Identifier"="node_concept_id")) %>%
      select(Cut.No., Node.Identifier, node_concept, Tree.Level, Relative.Risk, Log.Likelihood.Ratio, Observations, Expected, P.value) %>%
      filter(!is.na(P.value), Tree.Level>=DRUG_TREE_LEVEL, Observations>DRUG_N_OBS, P.value<DRUG_P_VAL) %>%
      arrange(desc(Log.Likelihood.Ratio))
    
    make_features(drug_cut_list[[i]], all_nodes_drug, type="drug") %>%
      output_tbl(paste("drug_features_untested_", i, sep=""), indices=c("ancestor_concept_id", "descendant_concept_id))"))
    
  }
  
  require(VennDiagram)
  venn.diagram(
    x = list(drug_cut_list[[1]] %>% pull(Node.Identifier), 
             drug_cut_list[[2]] %>% pull(Node.Identifier), 
             drug_cut_list[[3]] %>% pull(Node.Identifier), 
             drug_cut_list[[4]] %>% pull(Node.Identifier),
             drug_cut_list[[5]] %>% pull(Node.Identifier)),
    category.names = c("Fold 1" , "Fold 2 " , "Fold 3", "Fold 4", "Fold 5"),
    filename = './results/drug_venn_diagram.png',
    output=TRUE
  )
  
  
  
  drug_features_untested_all=results_tbl("drug_features_untested_all") %>% compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
  cohort_untested_drug_features_all<-make_feature_tbl(cohort=results_tbl("cohort_pasc_or_misc_untested"),
                                                      type="drug", features=drug_features_untested_all, days_delta=days_delta,
                                                      split_dates=split_dates, days_min=days_min, days_max=days_max)
  cohort_untested_drug_features_all %>%
    output_tbl("cohort_untested_drug_features_all", index="person_id")
  
#  write.csv(cohort_untested_drug_features_all, "./specs/cohort_features",str_append,"/cohort_drug_features_all.csv", row.names=FALSE)
#  rm(cohort_untested_drug_features_all)
  
  #  drug_feature_list<-list()
  #  cohort_drug_feature_train_list<-list()
  #  cohort_drug_feature_test_list<-list()
  
  for (i in 1:5){
    print(i)
    drug_features=results_tbl(paste("drug_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_drug_feature_train<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")),
                                               type="drug", features=drug_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_drug_feature_train %>% output_tbl(paste("cohort_drug_features_train_", i, sep=""), index="person_id")
    
    print(i)
    drug_features=results_tbl(paste("drug_features_untested_", i, sep="")) 
    a=Sys.time()
    cohort_drug_feature_test<-make_feature_tbl(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_test_", i, sep="")),
                                              type="drug", features=drug_features, days_delta=days_delta, split_dates=split_dates, days_min=days_min, days_max=days_max)
    print(Sys.time()-a)
    cohort_drug_feature_test %>% output_tbl(paste("cohort_drug_features_test_", i, sep=""), index="person_id")
  }
  
  
  
#  
#  ## Cluster features (sensitivity analysis)
#  rules<-results_tbl("post_infectious_disorder_rules")
#  
#  
#  clust_visits_train_list<-list()
#  clust_visits_test_list<-list()
#  
#  for (i in 1:5){
#    clust_visits_train<-get_cluster_visits(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep=""))) %>% remove_pasc_obs_clust()
#    clust_visits_train_wide_temp <-partition_dates_and_pivot_wide_clust_features(cohort=clust_visits_train)
#    clust_visits_train_list[[i]]<-results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")) %>% distinct(person_id, index_date, observation_type, pasc_flag) %>% collect %>%
#      left_join(clust_visits_train_wide_temp, by=c("person_id", "index_date", "observation_type", "pasc_flag"))
#    clust_visits_train_list[[i]][,5:ncol(clust_visits_train_list[[i]])]<- clust_visits_train_list[[i]][,5:ncol(clust_visits_train_list[[i]])] %>% replace(is.na(.), 0)
#    write.csv(clust_visits_train_list[[i]], paste("./specs/cohort_features/cohort_clust_features_train_", i, ".csv", sep=""), row.names=FALSE)
#    
#    clust_visits_test<-get_cluster_visits(cohort=results_tbl(paste("cohort_pasc_or_misc_untested_test_", i, sep=""))) %>% remove_pasc_obs_clust()
#    clust_visits_test_wide_temp <-partition_dates_and_pivot_wide_clust_features(cohort=clust_visits_test)
#    clust_visits_test_list[[i]]<-results_tbl(paste("cohort_pasc_or_misc_untested_test_", i, sep="")) %>% distinct(person_id, index_date, observation_type, pasc_flag) %>% collect %>%
#      left_join(clust_visits_test_wide_temp, by=c("person_id", "index_date", "observation_type", "pasc_flag"))
#    clust_visits_test_list[[i]][,5:ncol(clust_visits_test_list[[i]])]<- clust_visits_test_list[[i]][,5:ncol(clust_visits_test_list[[i]])] %>% replace(is.na(.), 0)
#    write.csv(clust_visits_test_list[[i]], paste("./specs/cohort_features/cohort_clust_features_test_", i, ".csv", sep=""), row.names=FALSE)
#  }
  
  
  
  
  fix_demo_cols<-function(columns){

    columns[grepl("cchmc", columns)]="Site H"
    columns[grepl("chop", columns)]="Site C"
    columns[grepl("colorado", columns)]="Site B"
    columns[grepl("lurie", columns)]="Site D"
    columns[grepl("nationwide", columns)]="Site I"
    columns[grepl("nemours", columns)]="Site A"
    columns[grepl("seattle", columns)]="Site G"
    columns[grepl("stanford", columns)]="Site F"
    columns[grepl("Male", columns)]="Sex: Male"
    columns[grepl("Female", columns)]="Sex: Female"
    columns[grepl("Asian/PI", columns)]="Race/ethnicity: NH Asian/PI"
    columns[grepl("Black/AA", columns)]="Race/ethnicity: NH Black/AA"
    columns[grepl("Hispanic", columns)]="Race/ethnicity: Hispanic"
    columns[grepl("White", columns)]="Race/ethnicity: NH White"
    columns[grepl("Other/unknown", columns)]="Race/ethnicity: Other/unknown"
    columns[grepl("Multiple", columns)]="Race/ethnicity: Multiple"
    columns[grepl("age", columns)]="Age"
    columns[grepl("mar_jul_22", columns)]="Cohort entry: Mar-Jul 2022"
    columns[grepl("mar_jun_21", columns)]="Cohort entry: Mar-Jun 2021"
    columns[grepl("jul_oct_21", columns)]="Cohort entry: Jul-Oct 2021"
    columns[grepl("nov_feb_21", columns)]="Cohort entry: Nov 2020-Feb 2021"
    columns[grepl("nov_feb_22", columns)]="Cohort entry: Nov 2021-Feb 2022"
    
    return(columns)
  }
  
  fix_cond_cols<-function(string){
    date=substr(string, 1, 8)
    date[date=="1m_to_2m"]="(1 to 2 m)"
    date[date=="3m_to_6m"]="(3 to 6 m)"
    date[date=="-1m_to_0"]="(-1 to 0 m)"
    date[date=="0m_to_1m"]="(0 to 1 m)"
    date[date=="2m_to_3m"]="(2 to 3 m)"
    rest=substr(string, 10, nchar(string))
    rest[date=="(-1 to 0 m)"]=substr(string, 11, nchar(string))[date=="(-1 to 0 m)"]
    text_temp=str_extract(rest, "^(.+?),")
    text=substr(text_temp, 1, nchar(text_temp)-1)
    string_new=paste(text, date, sep=" ")
    return(string_new)
  }
  
  fix_drug_cols<-function(string){
    date=substr(string, 1, 8)
    date[date=="1m_to_2m"]="(1 to 2 m)"
    date[date=="3m_to_6m"]="(3 to 6 m)"
    date[date=="-1m_to_0"]="(-1 to 0 m)"
    date[date=="0m_to_1m"]="(0 to 1 m)"
    date[date=="2m_to_3m"]="(2 to 3 m)"
    rest=substr(string, 10, nchar(string))
    rest[date=="(-1 to 0 m)"]=substr(string, 11, nchar(string))[date=="(-1 to 0 m)"]
    text_temp=str_extract(rest, "^(.+?),")
    text_temp2=substr(text_temp, 1, nchar(text_temp)-1)
    text_temp3=str_extract(text_temp2, "_[^_]+$")
    text=str_to_sentence(substr(text_temp3, 2, nchar(text_temp3)))
    string_new=paste(text, date, sep=" ")
    return(string_new)
  }
  
  fix_util_cols<-function(string){
    
    date=substr(string, nchar(string)-7, nchar(string))
    date[date=="1m_to_2m"]="(1 to 2 m)"
    date[date=="3m_to_6m"]="(3 to 6 m)"
    date[date=="1m_to_0m"]="(-1 to 0 m)"
    date[date=="0m_to_1m"]="(0 to 1 m)"
    date[date=="2m_to_3m"]="(2 to 3 m)"
    rest=substr(string, 1, nchar(string)-9)
    rest[date=="(-1 to 0 m)"]=substr(string, 1, nchar(string)-10)[date=="(-1 to 0 m)"]
    rest[rest=="Inpatient_non_ICU_"]="Inpatient non-ICU"
    rest[rest=="Inpatient_ICU_"]="Inpatient ICU"
    rest[rest=="Inpatient_ICU"]="Inpatient ICU"
    rest[rest=="ED_ICU"]="ED ICU"
    rest[rest=="ED_ICU_"]="ED ICU"
    rest[rest=="ED_non_ICU_"]="ED non-ICU"
    rest[rest=="ED_non_ICU"]="ED non-ICU"

    string_new=paste(rest, "visits", date, sep=" ")
    
    return(string_new)
  }
  
  fix_lab_cts_cols<-function(columns){
    columns[grepl("0m_to_1m", columns)]="Lab count (0 to 1 m)"
    columns[grepl("-1m_to_0m", columns)]="Lab count (-1 to 0 m)"
    columns[grepl("1m_to_2m", columns)]="Lab count (1 to 2 m)"
    columns[grepl("2m_to_3m", columns)]="Lab count (2 to 3 m)"
    columns[grepl("3m_to_6m", columns)]="Lab count (3 to 6 m)"
    return(columns)
  }

  ###### Clean up column names
  all_cols<-read.csv("./specs/X_test_columns.csv") %>% pull(X0)
  demo_cols<-read.csv("./specs/demo_columns.csv") %>% pull(X0) %>% fix_demo_cols() %>% write.csv("./specs/demo_columns_fixed.csv", row.names=FALSE)
  cond_cols<-read.csv("./specs/cond_columns.csv") %>% pull(X0) %>% fix_cond_cols()%>% write.csv("./specs/cond_columns_fixed.csv", row.names=FALSE)
  drug_cols<-read.csv("./specs/drug_columns.csv") %>% pull(X0) %>% fix_cond_cols()%>% write.csv("./specs/drug_columns_fixed.csv", row.names=FALSE)
  lab_cols<-read.csv("./specs/lab_columns.csv") %>% pull(X0) %>% fix_cond_cols()%>% write.csv("./specs/lab_columns_fixed.csv", row.names=FALSE)
  proc_cols<-read.csv("./specs/proc_columns.csv") %>% pull(X0) %>% fix_cond_cols()%>% write.csv("./specs/proc_columns_fixed.csv", row.names=FALSE)
  util_cols<-read.csv("./specs/util_columns.csv") %>% pull(X0) %>% fix_util_cols()%>% write.csv("./specs/util_columns_fixed.csv", row.names=FALSE)
  lab_cts_cols<-read.csv("./specs/lab_cts_columns.csv") %>% pull(X0) %>% fix_lab_cts_cols %>% write.csv("./specs/lab_cts_columns_fixed.csv", row.names=FALSE)
  
  
  
  
  c("Lab count (0 to 1 m)", "Lab count (-1 to 0 m)", "Lab count (1 to 2 m)", "Lab count (2 to 3 m)")
  cond_new=as.character(sapply(cond, FUN=window_transform_cond))
  
  
  ## Chart review code    columns[grepl("mar_jul_22")]="Cohort entry: Mar-Jul 2022"

  #  chart_review_ids_a<-read.csv("./specs/chart_review_val.csv")$X0
  #  chart_review_ids_b<-read.csv("./specs/chart_review_val_b.csv")$X0
  #  
  #  chart_review_ids<-union(chart_review_ids_a, chart_review_ids_b)
  #  
  #  obs_der_chart_review<-cdm_tbl("observation_derivation_recover") %>% filter(person_id %in% chart_review_ids) %>% compute_new(index="person_id")
  #  
  #  chart_review_persons<-cdm_tbl('person') %>% filter(person_id %in% chart_review_ids) %>% compute_new(index="person_id")
  #  
  #  chart_review_sample<-chart_review_persons %>% 
  #    group_by(site) %>%
  #    slice_sample(n=2) %>%
  #    select(person_id) %>%
  #    ungroup %>%
  #    compute_new(index="person_id")
  #  
  #  
  #  chart_review_sample %>% inner_join(pasc_dx, by="person_id") %>% summarize(n_distinct(person_id))
  #  
  #  
  #  
  #  chart_review_sample %>% inner_join(cdm_tbl('person') %>% select(person_id, site_id), by="person_id") %>% output_tbl("chart_review_sample", index="person_id") 
  #  
  #  chart_review_sample %>% inner_join(cdm_tbl('person') %>% select(person_id, site_id), by="person_id") %>% write.csv("chart_review_sample.csv", row.names=FALSE)
  
  
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
  .run_features(base_dir)
}

