#'
#' This file contains functions to identify cohorts in a study.  Its
#' contents, as well as the contents of any R files whose name begins
#' with "cohort_", will be automatically sourced when the request is
#' executed.
#'
#' For simpler requests, it will often be possible to do all cohort
#' manipulation in functions within this file.  For more complex
#' requests, it may be more readable to separate different cohorts or
#' operations on cohorts into logical groups in different files.
#'


#' Constructs ML phenotype cohort
#' Inclusion criteria:
#' Patients who were COVID positive either by PCR or antigen or qualifying serology test
#' or had a MIS-C or PASC diagnosis.
#' 
#' For patients who had a positive test or COVID diagnosis, index date is date of first positive event.
#' For PASC patients with no prior positive COVID test or COVID diagnosis, index date
#' is currently imputed as 28 days prior to earliest PASC/MIS-C dx. 
#' 
#' PASC patients are included by a U09.9 diagnosis or a more general post infectious disorder diagnosis--
#' the latter are not included in observation_derivation_recover so we use a separate rules set (based on source concept id
#' and source value) to decided which should be included.
#' 
#' 
#' @param max_date The index date cutoff (update when running on new versions of the RECOVER data)
#' 
#' @return Cohort table with the following variables: `person_id`, `pasc_flag` (whether patient had a PASC or MIS-C diagnosis),
#' `pasc_not_misc_flag` (whether patient had a PASC diagnosis and no MIS-C diagnosis), `index_date` (as described above), 
#' `observation_type` (whether earliest positive observation was PCR, antigen, or serology test or PASC or MIS-C diagnosis),
#' `index_visit_occurrence_id` (visit_occurrence_id of index event)


make_cohort<-function(max_date){
  
 obs_der<-cdm_tbl("observation_derivation_recover") %>% filter(observation_date<=as.Date(max_date)) 
  pcr_positive<-obs_der %>%
    filter(observation_concept_id==2000001530L, value_as_concept_id %in% c(9191L, 2000001526L))
  antigen_positive<-obs_der %>%
    filter(observation_concept_id==2000001529L, value_as_concept_id %in% c(9191L, 2000001526L))
  serology_positive<-obs_der %>%
    filter(observation_concept_id==2000001528L, value_as_concept_id %in% c(9191L, 2000001526L))
#  covid_dx_specific<-obs_der %>%
#    filter(observation_concept_id==2000001527L, value_as_concept_id==2000001525L)
#  covid_dx_complication<-obs_der %>%
#    filter(observation_concept_id==2000001527L, value_as_concept_id==2000001523L)
#  covid_dx_history<-obs_der %>%
#    filter(observation_concept_id==2000001527L, value_as_concept_id==2000001522L)
#  covid_dx_any<-covid_dx_specific %>%
#    dplyr::union(covid_dx_complication) %>%
#    dplyr::union(covid_dx_history)
  pasc_dx<-obs_der %>%
    filter(observation_concept_id==2000001527L, value_as_concept_id==2000001520L)
  misc_dx<-obs_der %>%
    filter(observation_concept_id==2000001527L, value_as_concept_id==703578L)
  
  
  positive_any<-pcr_positive %>%
    dplyr::union(antigen_positive) %>%
    dplyr::union(serology_positive)
  
#  pasc_not_misc<-pasc_dx %>%
#    anti_join(misc_dx, by="person_id")
  
  pasc_non_misc<-pasc_dx %>% 
    anti_join(misc_dx, by="person_id") %>%
    mutate(pasc_not_misc_flag=1) %>%
    distinct(person_id, pasc_not_misc_flag)
  cases<-pasc_dx %>% 
    dplyr::union(misc_dx) %>%
    mutate(pasc_flag=1) %>%
    left_join(pasc_non_misc, by="person_id") %>%
    mutate(pasc_not_misc_flag=case_when(
      is.na(pasc_not_misc_flag)~0,
      TRUE~pasc_not_misc_flag
    )) %>% 
    group_by(person_id) %>%
    slice_min(observation_date, with_ties=FALSE) %>%
    ungroup %>%
    mutate(index_date=observation_date) %>%
    mutate(observation_type=case_when(
      observation_concept_id==2000001527L & value_as_concept_id==2000001520L~"pasc",
      observation_concept_id==2000001527L & value_as_concept_id==703578L~"misc"
    )) %>%
    compute_new(index="person_id")
  
  controls<-positive_any %>%
    anti_join(cases, by="person_id") %>%
    group_by(person_id) %>%
    slice_min(observation_date, with_ties=FALSE) %>%
    ungroup %>%
    mutate(pasc_flag=0, pasc_not_misc_flag=0) %>%
    mutate(observation_type=case_when(
      observation_concept_id==2000001530L~"pcr",
      observation_concept_id==2000001529L~"antigen",
      observation_concept_id==2000001528L~"serology"#,
      #      observation_concept_id==2000001527L~"covid_dx"
    ))  %>%
    mutate(index_date=observation_date) %>% 
    compute_new(index="person_id")
  
  cohort=cases %>% dplyr::union(controls) %>%
    rename(index_visit_occurrence_id=visit_occurrence_id) %>% 
    distinct(person_id, pasc_flag, pasc_not_misc_flag, index_date, observation_type, index_visit_occurrence_id) %>% 
#    filter(index_date>=as.Date("2021-01-01")) %>% 
    compute_new(index="person_id")
  

    
  
  message("adjusting to include post infectious disorder patients")
  cohort_final<-code_adjust(cohort=cohort, max_date=max_date) 
  
#  ## Inelegant but this bit of code sets the index date to be date of first positive test for
#  ## those PASC patients who had a prior positive test.
#  pasc_misc_fix<-cohort_final %>% 
#    filter(pasc_flag==1) %>% 
#    inner_join(positive_any %>% select(person_id, observation_concept_id, observation_date, visit_occurrence_id), by="person_id") %>%
#    filter(observation_date<index_date+days(28), observation_date>index_date-months(6)) %>%
#    group_by(person_id) %>%
#    slice_min(observation_date, with_ties=FALSE) %>%
#    ungroup %>%
#    mutate(index_date=observation_date) %>% 
#    mutate(index_visit_occurrence_id=visit_occurrence_id) %>%
#    mutate(observation_type=case_when(
#      observation_concept_id==2000001530L~"pcr",
#      observation_concept_id==2000001529L~"antigen",
#      observation_concept_id==2000001528L~"serology"#,
#      #      observation_concept_id==2000001527L~"covid_dx"
#    )) %>%
#    distinct(person_id, pasc_flag, pasc_not_misc_flag, index_date, observation_type, index_visit_occurrence_id) %>%
#    compute_new(index="person_id")
#  
#  cohort_fixed=(cohort_final %>%
#                  anti_join(pasc_misc_fix, by="person_id")) %>%
#    dplyr::union(pasc_misc_fix) %>%
#    filter(index_date>=as.Date("2021-01-01")) %>% 
#    compute_new(index="person_id")
  
      

  cohort_final %>% return()
    
  
}


#' Adjusts ML phenotype cohort to include post-infectious disorder diagnosed patients among cases
#' and exclude them from controls
#' 
#' @param cohort cohort table with `pasc_flag`, `pasc_not_misc_flag` `index_date`, and `observation_type` columns
#' @param rules for descendants of post infectious disorder code in the SNOMED hierarchy a
#' set of rules for each combination of condition_source_concept_id and condition_source_value
#' whether to categorize as pasc, ignore, or exclude.
#' 
#' @return Cohort table with the same columns as input table


code_adjust<-function(cohort, rules=results_tbl("post_infectious_disorder_rules"), max_date){
  cond_visit_tbl<-cdm_tbl('condition_occurrence') %>% 
    inner_join(cdm_tbl('visit_occurrence') %>% 
                 select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id") %>%
    filter(visit_start_date<as.Date(max_date))
  
  controls_switch<-cohort %>% filter(pasc_flag==0) %>% 
    inner_join(cond_visit_tbl, by="person_id") %>%
    filter(visit_start_date> index_date) %>% 
    inner_join(rules, by=c("condition_concept_id", "condition_source_concept_id", "condition_source_value")) %>%
    filter(action=="pasc") %>% 
    group_by(person_id) %>% 
    slice_min(visit_start_date, with_ties=FALSE) %>%
    ungroup %>%
    mutate(temp=TRUE) %>%
    rename(pseudo_pasc_date=visit_start_date) %>%
    distinct(person_id, temp, pseudo_pasc_date) %>% 
    compute_new(index="person_id")
  
  controls_drop<-cohort %>% filter(pasc_flag==0) %>% 
    inner_join(cond_visit_tbl, by="person_id") %>%
    filter(visit_start_date> index_date) %>% 
    inner_join(rules, by=c("condition_concept_id", "condition_source_concept_id", "condition_source_value")) %>%
    filter(action=="drop") %>% 
    distinct(person_id) %>%
    compute_new(index="person_id")
    
  message("fixing codes")
  cohort_fixed<-cohort %>% 
    anti_join(controls_drop, by="person_id") %>% 
    left_join(controls_switch, by="person_id") %>%
    mutate(pasc_flag=case_when(
      temp==TRUE~1,
      TRUE~pasc_flag
    ),
    pasc_not_misc_flag=case_when(
      temp==TRUE~1,
      TRUE~pasc_not_misc_flag
    ),
    index_date=case_when(
      temp==TRUE~pseudo_pasc_date,
      TRUE~index_date
    )) %>%
    select(-temp, -pseudo_pasc_date) %>%
    mutate(observation_type=case_when(
      (pasc_flag==1) & (pasc_not_misc_flag==1)~"non_misc_pasc",
      (pasc_flag==1)~"misc",
      pasc_flag==0~observation_type
      ))%>% 
    compute_new(index="person_id")
  return(cohort_fixed)
}
  


#' Filters to require specified number of post-index date visits in cohort
#' 
#' @param cohort cohort table with `index_date` column
#' @param n_visits number of visits to require post-index date
#' 
#' @return Cohort table with the same columns as input table

filter_visits<-function(cohort, n_visits){
  keep_ids<-cohort %>% 
    inner_join(
    cdm_tbl("visit_occurrence") %>% select(person_id, visit_start_date), by="person_id"
    ) %>%
    filter(visit_start_date>index_date) %>%
    group_by(person_id) %>%
    summarize(n_vis=n_distinct(visit_start_date)) %>%
    filter(n_vis>=n_visits) %>% distinct(person_id) %>% compute_new(index="person_id")
  
  cohort %>%
    inner_join(keep_ids, by="person_id") %>%
    compute_new(index="person_id")
  
}

#' Filters to include cohort entry dates in specified window
#' Because cases have index date imputed between 28 days and 90 days prior
#' to earliest PASC dx, we move the max date back by a month
#' (since otherwise there would be controls but not cases in that last month)
#' 
#' @param cohort cohort table with `index_date` column
#' @param min_date
#' @param max_date
#' @return Cohort table with the same columns as input table

filter_dates<-function(cohort, min_date, max_date){
  cohort %>%
    mutate(keep=case_when(
      (pasc_flag==1) & (index_date>=as.Date(min_date)) & (index_date<as.Date(max_date)-months(1))~1,
      (pasc_flag==0) & (index_date>=as.Date(min_date)) & (index_date<as.Date(max_date)-months(1))~1,
      TRUE~0
    )) %>%
    filter(keep==1) %>%
    select(-keep) %>%
    compute_new(index="person_id")
}



#' For patients with a PASC or MIS-C diagnosis, imputes index date by picking
#' a random date in the 28 days to 90 days prior to earliest dx
#' 
#' @param cohort cohort table with person_id, pasc_flag, and index_date columns
#' @param min_window min # of days for admissible range in which to impute index date
#' @param max_window max # of days for admissible range in which to impute index date
#' 
#' @return cohort with index dates imputed for patients with pasc_flag==1
#' 
impute_dates<-function(cohort, min_window=-90, max_window=-28, seed=1234){
  
  set.seed(seed)
  
  n=cohort %>% filter(pasc_flag==1) %>%
    summarize(n=n()) %>% pull(n)
  
  date_window=min_window:max_window
  
  shifts<-data.frame(date_shift=date_window) %>%
    output_tbl("shifts", temp=TRUE)
  
  pasc_adjust<-cohort %>%
    filter(pasc_flag==1) %>%
    merge(shifts) %>% 
    group_by(person_id) %>%
    slice_sample(n=1) %>%
    ungroup %>% 
    mutate(index_date=as.Date(index_date)+days(date_shift)) %>%
    select(-date_shift) %>% 
    output_tbl("pasc_adjust", temp=TRUE, index="person_id")
  
  cohort_final<-(cohort %>% 
                   filter(pasc_flag==0) %>%
                   anti_join(pasc_adjust, by="person_id")) %>%
    dplyr::union(pasc_adjust) %>%
    compute_new(index="person_id")
  
  return(cohort_final)
}



#' Creates 5-fold cross validation split of cohort person_ids
#' 
#' @param cohort_demo cohort table (with one row for each patient)
#' @param seed 
#' 
#' @return list of 1) list of train_id's, 2) list of test_id's
#' 

split_cohort<-function(cohort_demo, seed=0){
  set.seed(seed)
  folds <- create_folds(cohort_demo$pasc_flag, k = 5)
  fold1_train_ids<-cohort_demo[folds$Fold1,]$person_id
  fold1_test_ids<-cohort_demo[-folds$Fold1,]$person_id
  fold2_train_ids<-cohort_demo[folds$Fold2,]$person_id
  fold2_test_ids<-cohort_demo[-folds$Fold2,]$person_id
  fold3_train_ids<-cohort_demo[folds$Fold3,]$person_id
  fold3_test_ids<-cohort_demo[-folds$Fold3,]$person_id
  fold4_train_ids<-cohort_demo[folds$Fold4,]$person_id
  fold4_test_ids<-cohort_demo[-folds$Fold4,]$person_id
  fold5_train_ids<-cohort_demo[folds$Fold5,]$person_id
  fold5_test_ids<-cohort_demo[-folds$Fold5,]$person_id
  
  train_id_list=list(fold1_train_ids, fold2_train_ids, fold3_train_ids, fold4_train_ids, fold5_train_ids)
  test_id_list=list(fold1_test_ids, fold2_test_ids, fold3_test_ids, fold4_test_ids, fold5_test_ids)

  id_list=list(train_id_list, test_id_list)
  return(id_list)
}




#' Gets cohort visits for conditions during post-acute period (as specified by days_delta past index date) and categorizes
#' them by cluster using the cluster_master codeset.
#' 
#' @param cohort cohort table with `index_date` column
#' @param cluster_master cluster_master codeset which sorts ICD10CM codes into
#' clusters
#' @param cond_visit_tbl condition_occurrence joined to visit_occurrence along
#' visit_occurrence_id to get visit_start_dates
#' @param days_delta number of days post index date to include visits
#' 
#' @return Table with cohort visits and cluster designation for all ICD10CM codes.


get_cluster_visits<-function(cohort,
                             cluster_master=results_tbl('cluster_master') %>% filter(! cluster %in% c('pasc', 'multisystem inflammatory syndrome')),
                             cond_visit_tbl=cdm_tbl('condition_occurrence') %>% 
                               inner_join(cdm_tbl('visit_occurrence') %>% 
                                            select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id"),
                             days_delta=28){
  cohort %>% inner_join(cond_visit_tbl, by="person_id") %>% 
    inner_join(cluster_master, by=c("condition_source_concept_id"="concept_id")) %>%
    filter(visit_start_date>=index_date+days(days_delta)) %>%
    compute_new(index="person_id")
}




#' Same as above but includes visit_loc variable
get_cluster_visits_alt<-function(cohort,
                             cluster_master=results_tbl('cluster_master'),
                             cond_visit_tbl=cdm_tbl('condition_occurrence') %>% 
                               inner_join(cdm_tbl('visit_occurrence') %>% 
                                            select(visit_occurrence_id, visit_concept_id, visit_start_date), by="visit_occurrence_id") %>%
                               left_join(cdm_tbl('adt_occurrence') %>% 
                                           filter(service_concept_id %in% c(2000000079L,2000000080L,2000000078L)) %>%
                                           distinct(visit_occurrence_id) %>% mutate(icu=1), by="visit_occurrence_id"),
                             months_delta=0){
  
  cohort %>% inner_join(cond_visit_tbl, by="person_id") %>% 
    inner_join(cluster_master, by=c("condition_source_concept_id"="concept_id")) %>%
    filter(visit_start_date>=index_date+months(months_delta)) %>%
    mutate(visit_loc=case_when(
      (visit_concept_id %in% c(9201L,2000000088L,2000000048L) & icu==1)~"Inpatient_ICU",
      (visit_concept_id %in% c(9201L,2000000088L,2000000048L) & is.na(icu))~"Inpatient_non_ICU",
       (visit_concept_id %in% c(9203L) & icu==1)~"ED_ICU",
       (visit_concept_id %in% c(9203L) & is.na(icu))~"ED_non_ICU",
       visit_concept_id %in% c(9202L,581399L) ~ 'Outpatient Office',
       visit_concept_id %in% c(2000000469L,44814711L) ~ 'Outpatient: Test Only',
       icu==1~"ICU Unknown",
       TRUE ~ 'Other/Unknown'
      ))%>% 
    select(-visit_concept_id, -icu) %>%
    compute_new(index="person_id")
}


#' Gets cohort visits for drugs during post-acute period (as specified by days_delta past index date) at
#' and below the ingredient level in the drug hierarchy
#' @param cohort cohort table with `index_date` column
#' @param drug_tbl drug_exposure table
#' @param days_delta number of days post index date to include visits
#' 
#' @return Table with all drug exposures at the ingredient level or lower

get_drugs<-function(cohort,
                    drug_tbl=cdm_tbl("drug_exposure") %>% select(person_id, drug_concept_id, drug_visit_start_date=drug_exposure_start_date),
                    days_delta=28){
  ingredient_descendants<-vocabulary_tbl("concept") %>% 
    filter(domain_id=="Drug") %>% 
    filter(concept_class_id=="Ingredient") %>% 
    inner_join(vocabulary_tbl("concept_ancestor"), by=c("concept_id"="ancestor_concept_id")) %>% 
    distinct(concept_name, descendant_concept_id)
  
  cohort %>% inner_join(drug_tbl, by="person_id") %>%
    inner_join(ingredient_descendants, by=c("drug_concept_id"="descendant_concept_id")) %>%
    filter(drug_visit_start_date>=index_date+days(days_delta)) %>%
    compute_new(index="person_id")
}


#' Gets cohort visits by visit type during post-acute period (as specified by days_delta past index date) 
#' @param cohort cohort table with `index_date` column
#' @param adt_tbl adt_occurrence table to capture ICU admissions
#' @param visit_tbl visit_occurrence table
#' @param days_delta number of days post index date to include visits
#' 
#' @return Table with cohort visits by type and ICU/non-ICU designation

get_util<-function(cohort_tbl,
                   adt_tbl = cdm_tbl('adt_occurrence') %>% 
                     filter(service_concept_id %in% c(2000000079L,2000000080L,2000000078L)) %>%
                     distinct(visit_occurrence_id) %>% mutate(icu=1),
                     visit_tbl=cdm_tbl("visit_occurrence") %>% 
                       select(person_id, visit_occurrence_id, util_start_date=visit_start_date, visit_concept_id),
                     days_delta=28){
  
  visit_icu_tbl<-visit_tbl %>% 
    left_join(adt_tbl, by="visit_occurrence_id")%>%
    mutate(icu=case_when(
      is.na(icu)~0,
      TRUE~1
    )) %>% select(-visit_occurrence_id)

  cohort_tbl %>% left_join(visit_icu_tbl, by="person_id") %>%
    filter(util_start_date>=index_date+days(days_delta)) %>%
    mutate(util_loc=case_when(
      (visit_concept_id %in% c(9201L,2000000088L,2000000048L) & icu==1)~"Inpatient_ICU",
      (visit_concept_id %in% c(9201L,2000000088L,2000000048L) & icu==0)~"Inpatient_non_ICU",
      (visit_concept_id %in% c(9203L) & icu==1)~"ED_ICU",
      (visit_concept_id %in% c(9203L) & icu==0)~"ED_non_ICU",
      visit_concept_id %in% c(9202L,581399L) ~ 'Outpatient Office',
      visit_concept_id %in% c(2000000469L,44814711L) ~ 'Outpatient: Test Only',
      TRUE ~ 'Other/Unknown'
    ))%>% 
    select(-visit_concept_id, -icu) %>% 
    compute_new(index="person_id")
}


#' Partitions visits by start date and pivots to wide format (one column for each combination of visit_type and time period)
#' @param cohort table with a row for each visit and `util_start_date`, `util_loc`, `index_date`, `observation_type`, and `pasc_flag` columns.

#' @return Table with a row for each a patient and a column for each visit type and time window combination
#' with value the number of visits of that type during that time period.

partition_dates_and_pivot_wide_util<-function(cohort, split_dates=TRUE, days_min=-28, days_max=180){
  if (split_dates){
    cohort_wide<-cohort %>%
      mutate(util_window=case_when(
        #      util_start_date<index_date-days(30)~"-6m_to_-1m",
        util_start_date<index_date+days(0)~"-1m_to_0m",
        util_start_date<index_date+days(30)~"0m_to_1m",
        util_start_date<index_date+days(60)~"1m_to_2m",
        util_start_date<index_date+days(90)~"2m_to_3m",
        util_start_date<index_date+days(180)~"3m_to_6m",
        util_start_date>=index_date+days(180)~">6m"
      )) %>%
      filter(util_window!=">6m") %>%
      group_by(person_id, index_date, observation_type, pasc_flag, util_loc, util_window) %>%
      summarize(n_visits=n_distinct(util_start_date)) %>%
      ungroup %>%
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(util_loc, util_window), values_from=n_visits, values_fill=0) #%>%
      #compute_new(index="person_id")
    
    return(cohort_wide)
  }else{
    cohort_wide<-cohort %>%
      filter(util_start_date>index_date+days(days_min), util_start_date<index_date+days(days_max)) %>% 
      group_by(person_id, index_date, observation_type, pasc_flag, util_loc) %>%
      summarize(n_visits=n_distinct(util_start_date)) %>%
      ungroup %>%
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(util_loc), values_from=n_visits, values_fill=0) #%>%
      #compute_new(index="person_id")
  }
  
}






print_split_probs<-function(max_date, follow_months){
  for (i in 1:5){
    n<-results_tbl(paste("cohort_pasc_or_misc_untested_train_", i, sep="")) %>% 
      filter(index_date<as.Date(!!max_date)-months(!!follow_months)) %>%
      group_by(pasc_flag) %>% summarize(n=n()) %>% collect  %>%  arrange(desc(pasc_flag)) %>%
      pull(n)
    
    print(paste("Split i case probability=", n[1], "/", n[1]+n[2], sep=""))
  }
  n<-results_tbl("cohort_pasc_or_misc_untested") %>% 
    filter(index_date<as.Date(!!max_date)-months(!!follow_months)) %>%
    group_by(pasc_flag) %>% summarize(n=n()) %>% collect  %>%  arrange(desc(pasc_flag)) %>%
    pull(n)
  
  print(paste("Full cohort case probability=", n[1], "/", n[1]+n[2], sep=""))
  
}







partition_dates_and_pivot_wide_cond_features<-function(cohort, split_dates=TRUE, days_min=-28, days_max=180){
  if (split_dates){
    cohort_wide<-cohort %>%
      rename(cond_feature_id=ancestor_concept_id,
             cond_feature_name=ancestor_concept) %>% 
      #    filter(visit_start_date>=index_date+days(28)) %>% 
      mutate(clust_window=case_when(
        #      visit_start_date<index_date-days(30)~"-6m_to_-1m",
        visit_start_date<index_date+days(0)~"-1m_to_0m",
        visit_start_date<index_date+days(30)~"0m_to_1m",
        visit_start_date<index_date+days(60)~"1m_to_2m",
        visit_start_date<index_date+days(90)~"2m_to_3m",
        visit_start_date<index_date+days(180)~"3m_to_6m",
        visit_start_date>=index_date+days(180)~">6m"
      )) %>%
      filter(clust_window!=">6m") %>%
      group_by(person_id, index_date, observation_type, pasc_flag, cond_feature_id, cond_feature_name, clust_window) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup %>% #%>% collect
    #compute_new(index="person_id")
    #cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(cond_feature_name, cond_feature_id, clust_window), values_from=n_visits, values_fill=0, names_glue = "{clust_window}_{cond_feature_id}")# %>%
    #compute_new(index="person_id")
    return(cohort_wide)
  }else{
    cohort_wide<-cohort %>%
      rename(cond_feature_id=ancestor_concept_id,
             cond_feature_name=ancestor_concept) %>% 
      filter(visit_start_date>index_date+days(days_min), visit_start_date<index_date+days(days_max)) %>% 
      group_by(person_id, index_date, observation_type, pasc_flag, cond_feature_id, cond_feature_name) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(cond_feature_name, cond_feature_id), values_from=n_visits, values_fill=0, names_glue = "{cond_feature_id}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }
  
}

partition_dates_and_pivot_wide_clust_features<-function(cohort, split_dates=TRUE, days_min=-28, days_max=180){
  if (split_dates){
    cohort_wide<-cohort %>%
      mutate(clust_window=case_when(
        #      visit_start_date<index_date-days(30)~"-6m_to_-1m",
        visit_start_date<index_date+days(0)~"-1m_to_0m",
        visit_start_date<index_date+days(30)~"0m_to_1m",
        visit_start_date<index_date+days(60)~"1m_to_2m",
        visit_start_date<index_date+days(90)~"2m_to_3m",
        visit_start_date<index_date+days(180)~"3m_to_6m",
        visit_start_date>=index_date+days(180)~">6m"
      )) %>%
      filter(clust_window!=">6m") %>%
      group_by(person_id, index_date, observation_type, pasc_flag, cluster, clust_window) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(cluster, clust_window), values_from=n_visits, values_fill=0, names_glue = "{clust_window}_{cluster}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }else{
    cohort_wide<-cohort %>%
      filter(visit_start_date>index_date+days(days_min), visit_start_date<index_date+days(days_max)) %>% 
      group_by(person_id, index_date, observation_type, pasc_flag, cluster) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup# %>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(cluster), values_from=n_visits, values_fill=0, names_glue = "{cluster}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }
  
}


partition_dates_and_pivot_wide_lab_features<-function(cohort, split_dates=TRUE, days_min=-28, days_max=180){
  if (split_dates){
    cohort_wide<-cohort %>%
      rename(lab_feature_id=ancestor_concept_id,
             lab_feature_name=ancestor_concept) %>% 
      #    filter(visit_start_date>=index_date+days(28)) %>% 
      mutate(lab_window=case_when(
        #      visit_start_date<index_date-days(30)~"-6m_to_-1m",
        visit_start_date<index_date+days(0)~"-1m_to_0m",
        visit_start_date<index_date+days(30)~"0m_to_1m",
        visit_start_date<index_date+days(60)~"1m_to_2m",
        visit_start_date<index_date+days(90)~"2m_to_3m",
        visit_start_date<index_date+days(180)~"3m_to_6m",
        visit_start_date>=index_date+days(180)~">6m"
      )) %>%
      filter(lab_window!=">6m") %>%
      group_by(person_id, index_date, observation_type, pasc_flag, lab_feature_id, lab_feature_name, lab_window) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(lab_feature_name, lab_feature_id, lab_window), values_from=n_visits, values_fill=0, names_glue = "{lab_window}_{lab_feature_id}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }else{
    cohort_wide<-cohort %>%
      rename(lab_feature_id=ancestor_concept_id,
             lab_feature_name=ancestor_concept) %>% 
      filter(visit_start_date>index_date+days(days_min), visit_start_date<index_date+days(days_max)) %>% 
      group_by(person_id, index_date, observation_type, pasc_flag, lab_feature_id, lab_feature_name) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(lab_feature_name, lab_feature_id), values_from=n_visits, values_fill=0, names_glue = "{lab_feature_id}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }
  
  
}

partition_dates_and_pivot_wide_proc_features<-function(cohort, split_dates=TRUE, days_min=-28, days_max=180){
  if (split_dates){
    cohort_wide<-cohort %>%
      rename(proc_feature_id=ancestor_concept_id,
             proc_feature_name=ancestor_concept) %>% 
      #    filter(visit_start_date>=index_date+days(28)) %>% 
      mutate(proc_window=case_when(
        #      visit_start_date<index_date-days(30)~"-6m_to_-1m",
        visit_start_date<index_date+days(0)~"-1m_to_0m",
        visit_start_date<index_date+days(30)~"0m_to_1m",
        visit_start_date<index_date+days(60)~"1m_to_2m",
        visit_start_date<index_date+days(90)~"2m_to_3m",
        visit_start_date<index_date+days(180)~"3m_to_6m",
        visit_start_date>=index_date+days(180)~">6m"
      )) %>%
      filter(proc_window!=">6m") %>%
      group_by(person_id, index_date, observation_type, pasc_flag, proc_feature_id, proc_feature_name, proc_window) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(proc_feature_name, proc_feature_id, proc_window), values_from=n_visits, values_fill=0, names_glue = "{proc_window}_{proc_feature_id}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }else{
    cohort_wide<-cohort %>%
      rename(proc_feature_id=ancestor_concept_id,
             proc_feature_name=ancestor_concept) %>% 
      filter(visit_start_date>index_date+days(days_min), visit_start_date<index_date+days(days_max)) %>% 
      group_by(person_id, index_date, observation_type, pasc_flag, proc_feature_id, proc_feature_name) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(proc_feature_name, proc_feature_id), values_from=n_visits, values_fill=0, names_glue = "{proc_feature_id}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }
}


partition_dates_and_pivot_wide_drug_features<-function(cohort, split_dates=TRUE, days_min=-28, days_max=180){
  if (split_dates){
    cohort_wide<-cohort %>%
      rename(drug_feature_id=ancestor_concept_id,
             drug_feature_name=ancestor_concept) %>% 
      #    filter(visit_start_date>=index_date+days(28)) %>% 
      mutate(drug_window=case_when(
        #      visit_start_date<index_date-days(30)~"-6m_to_-1m",
        visit_start_date<index_date+days(0)~"-1m_to_0m",
        visit_start_date<index_date+days(30)~"0m_to_1m",
        visit_start_date<index_date+days(60)~"1m_to_2m",
        visit_start_date<index_date+days(90)~"2m_to_3m",
        visit_start_date<index_date+days(180)~"3m_to_6m",
        visit_start_date>=index_date+days(180)~">6m"
      )) %>%
      filter(drug_window!=">6m") %>%
      group_by(person_id, index_date, observation_type, pasc_flag, drug_feature_id, drug_feature_name, drug_window) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(drug_feature_name, drug_feature_id, drug_window), values_from=n_visits, values_fill=0, names_glue = "{drug_window}_{drug_feature_id}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }else{
    cohort_wide<-cohort %>%
      rename(drug_feature_id=ancestor_concept_id,
             drug_feature_name=ancestor_concept) %>% 
      filter(visit_start_date>index_date+days(days_min), visit_start_date<index_date+days(days_max)) %>% 
      group_by(person_id, index_date, observation_type, pasc_flag, drug_feature_id, drug_feature_name) %>%
      summarize(n_visits=n_distinct(visit_start_date)) %>%
      mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
      ungroup #%>% collect
    #compute_new(index="person_id")
    
    cohort_wide_final<-cohort_wide %>% 
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(drug_feature_name, drug_feature_id), values_from=n_visits, values_fill=0, names_glue = "{drug_feature_id}") #%>%
    #compute_new(index="person_id")
    return(cohort_wide_final)
  }
}




#partition_dates_and_pivot_wide<-function(cohort){
#  cohort_wide<-cohort %>%
#    mutate(clust_window=case_when(
##      visit_start_date<index_date-days(30)~"-6m_to_-1m",
#      visit_start_date<index_date+days(0)~"-1m_to_0m",
#      visit_start_date<index_date+days(30)~"0m_to_1m",
#      visit_start_date<index_date+days(60)~"1m_to_2m",
#      visit_start_date<index_date+days(90)~"2m_to_3m",
#      visit_start_date<index_date+days(180)~"3m_to_6m",
#      visit_start_date>=index_date+days(180)~">6m"
#    )) %>%
#    filter(clust_window!=">6m") %>%
#    group_by(person_id, index_date, observation_type, pasc_flag, cluster, clust_window) %>%
#    summarize(n_visits=n_distinct(visit_start_date)) %>%
#    mutate(n_visits=case_when(n_visits>0~1, TRUE~0)) %>%
#    ungroup %>%
#    pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
#                names_from=c(cluster, clust_window), values_from=n_visits, values_fill=0, names_glue = "{clust_window}_{cluster}") %>%
#    compute_new(index="person_id")
#}



join_cohort_demo<-function(cohort_tbl){
  cohort_demo<- cohort_tbl %>%
    inner_join(cdm_tbl("person") %>% dplyr::select(person_id, birth_date,site,
                                                   gender_concept_id, race_concept_id,
                                                   ethnicity_concept_id), by="person_id") %>%
    mutate(index_date=as.Date(index_date), birth_date=as.Date(birth_date)) %>% 
    mutate(entry_age=floor(as.numeric(index_date-birth_date)/365.25)) %>%
#    filter(entry_age<21) %>%
    mutate(sex_cat = case_when(gender_concept_id == 8507L ~ 'Male',
                               gender_concept_id == 8532L ~ 'Female',
                               TRUE ~ 'Other/unknown'),
           eth_cat = case_when(ethnicity_concept_id == 38003563L ~ 'Hispanic',
                               race_concept_id == 8516L ~ 'Black/AA',
                               race_concept_id %in% c(8515L, 8557L) ~
                                 'Asian/PI',
                               #                               race_concept_id == 8657L ~ 'Native American',
                               race_concept_id == 8527L ~ 'White',
                               race_concept_id == 44814659L ~ 'Multiple',
                               TRUE ~ 'Other/Unknown')) %>%
    dplyr::select(-gender_concept_id, -race_concept_id, -ethnicity_concept_id) %>%
#    filter(index_date<as.Date("2022-09-01")) %>% 
#    mutate(index_entry_month=paste(month(index_date), year(index_date), sep=",")) %>%
#    mutate(follow_months=ceiling((as.Date("2022-09-01")-index_date)/30)) %>% 
    filter(entry_age<21) %>%
    compute_new(index="person_id") 
}

get_cohort_entry_period<-function(cohort_tbl){
  cohort_demo<- cohort_tbl %>%
    mutate(index_date=as.Date(index_date)) %>% 
    mutate(index_month=month(index_date), index_year=year(index_date)) %>%
    mutate(cohort_entry_period=case_when(
      index_year==2020 | (index_year==2021 & index_month %in% c(1, 2))~"nov_feb_21",
      index_year==2021 & index_month %in% c(3, 4, 5, 6)~"mar_jun_21",
      index_year==2021 & index_month %in% c(7, 8, 9, 10)~"jul_oct_21",
      (index_year==2021 & index_month %in% c(11, 12))|(index_year==2022 & index_month %in% c(1, 2))~"nov_feb_22",
      (index_year==2022) & index_month %in% c(3, 4, 5, 6, 7)~"mar_jul_22"
    )) %>%
    filter(!is.na(cohort_entry_period)) %>% 
    compute_new(index="person_id") 
}





#' Makes table of ancestor-descendnat relationships representing the hierarchy in the 
#' desired vocabulary from the concept and concept_ancestor tables. In case of SNOMED diagnoses,
#' removes U09.9 and postviral disorder codes from the hierarchy from consideration
#' to avoid using these as features downstream.
#' 
#' @param type desired vocabulary for hierarchy: SNOMED, LOINC, proc (ICD10PCS, HCPCS, and CPT4),
#' or drug (RxNorm)
#' @param remove_postviral whether postviral disorder codes and their descendants should be removed from 
#' the tree

#' @return Table with `ancestor_concept_id`, `ancestor_concept_name`, `descendant_concept_id`, and `descendant_concept_name` columns
#' expressing the structure of the hierarchy.


make_tree<-function(type="SNOMED", remove_postviral=FALSE){
  if (type=="SNOMED"){
    snomed_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id=="SNOMED", domain_id=="Condition") %>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")
    
    snomed_tree<-snomed_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      filter(min_levels_of_separation==1) %>%
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
#      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    
    
    
    if(remove_postviral==TRUE){
      message('removing postviral codes')
      post_viral_descendants<-vocabulary_tbl("concept_ancestor") %>% 
        filter(ancestor_concept_id==444201L) %>% inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>%
        compute_new(index="descendant_concept_id")
      
      snomed_tree_temp<-(snomed_tree %>% 
        anti_join(post_viral_descendants, by=c("descendant_concept_id"))) %>%
        anti_join(post_viral_descendants, by=c("ancestor_concept_id"="descendant_concept_id")) %>%
        compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
      
      snomed_tree=snomed_tree_temp
    }
    
#    snomed_tree_all<-snomed_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
#      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
#      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
#      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
#      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
#      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))

    return(snomed_tree)
  }
  if (type=="LOINC"){
    loinc_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id=="LOINC", domain_id=="Measurement") %>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")
    loinc_tree<-loinc_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      filter(min_levels_of_separation==1) %>%
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
      #      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    return(loinc_tree)
  }
  if (type=="proc"){
    proc_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id %in% c("ICD10PCS", "HCPCS", "CPT4"), domain_id=="Procedure") %>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")
    proc_tree<-proc_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      filter(min_levels_of_separation==1) %>%
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
      #      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    return(proc_tree)
  }
  if(type=="drug"){
    drug_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id == 'RxNorm', domain_id=="Drug")%>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")
    drug_tree<-drug_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      filter(min_levels_of_separation==1) %>%
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
      #      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    return(drug_tree)
  }
}







get_snomed_counts<-function(tree, cohort, days_delta, method="bernoulli", follow_months=3, max_date){
  
  
  cond_visit_tbl=cdm_tbl('condition_occurrence') %>% 
    inner_join(cdm_tbl('visit_occurrence') %>% 
                 select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id")
  
  all_descendants<-tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors<-tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  
  all_nodes<-dplyr::union(all_descendants, all_ancestors)
  
  if (method=="bernoulli"){
    cohort_follow_up<-cohort %>% 
      filter(index_date<as.Date(!!max_date)-months(!!follow_months))
    
    cohort_visits<-cohort_follow_up %>% 
      inner_join(cond_visit_tbl, by="person_id") %>% 
      filter(visit_start_date>=index_date+days(!!days_delta), visit_start_date<index_date+months(!!follow_months)) %>% 
      inner_join(all_nodes, by=c("condition_concept_id"="node_concept_id")) %>%
      rename(node_concept_id=condition_concept_id) %>% 
      #group_by(person_id, pasc_flag, follow_months, node_concept_id, node_concept) %>%
      #summarize(n_person_visits=n_distinct(visit_start_date)) %>%
      #    mutate(follow_months=follow_months+1) %>% 
      #mutate(visit_count=floor(1000*n_person_visits/follow_months)) %>%
      group_by(pasc_flag, node_concept_id, node_concept) %>%
      summarize(n_persons=n_distinct(person_id)) %>%
      compute_new()
    
    counts<-cohort_visits %>% 
      collect %>%
      mutate(pasc_flag=as.factor(pasc_flag)) %>% 
      pivot_wider(id_cols=c("node_concept_id", "node_concept"), names_from=c("pasc_flag"), values_from=c("n_persons"), values_fill=0) %>%
      rename("positive"="1", "negative"="0") %>%
      mutate(pop_count=positive+negative, case_count=positive, control_count=negative,
             node=node_concept_id) %>%
      select(node, pop_count, case_count, control_count) 
  }
  
}


get_loinc_counts<-function(tree, cohort, days_delta, method="bernoulli", follow_months=3){
  
  
  lab_tbl=cdm_tbl('measurement_labs') %>% 
    inner_join(cdm_tbl('visit_occurrence') %>% 
                 select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id")
  
  all_descendants<-tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors<-tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  
  all_nodes<-dplyr::union(all_descendants, all_ancestors)
  
  if (method=="bernoulli"){
    cohort_follow_up<-cohort %>% 
      filter(index_date<as.Date("2022-09-01")-months(!!follow_months))
    
    cohort_visits<-cohort_follow_up %>% 
      inner_join(lab_tbl, by="person_id") %>% 
      filter(visit_start_date>=index_date+days(!!days_delta), visit_start_date<index_date+months(!!follow_months)) %>% 
      inner_join(all_nodes, by=c("measurement_concept_id"="node_concept_id")) %>%
      rename(node_concept_id=measurement_concept_id) %>% 
      #group_by(person_id, pasc_flag, follow_months, node_concept_id, node_concept) %>%
      #summarize(n_person_visits=n_distinct(visit_start_date)) %>%
      #    mutate(follow_months=follow_months+1) %>% 
      #mutate(visit_count=floor(1000*n_person_visits/follow_months)) %>%
      group_by(pasc_flag, node_concept_id, node_concept) %>%
      summarize(n_persons=n_distinct(person_id)) %>%
      compute_new()
    
    counts<-cohort_visits %>% 
#      mutate(pasc_flag=as.factor(pasc_flag)) %>% 
      pivot_wider(id_cols=c("node_concept_id", "node_concept"), names_from=c("pasc_flag"), values_from=c("n_persons"), values_fill=0) %>%
      rename("positive"="1", "negative"="0") %>%
      mutate(pop_count=positive+negative, case_count=positive, control_count=negative,
             node=node_concept_id) %>%
      select(node, pop_count, case_count, control_count) %>% collect
  }
  
}




get_proc_counts<-function(tree, cohort, days_delta, method="bernoulli", follow_months=3){
  
  
  proc_tbl=cdm_tbl('procedure_occurrence') %>% 
    inner_join(cdm_tbl('visit_occurrence') %>% 
                 select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id")
  
  all_descendants<-tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors<-tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  
  all_nodes<-dplyr::union(all_descendants, all_ancestors)
  
  if (method=="bernoulli"){
    cohort_follow_up<-cohort %>% 
      filter(index_date<as.Date("2022-09-01")-months(!!follow_months))
    
    cohort_visits<-cohort_follow_up %>% 
      inner_join(proc_tbl, by="person_id") %>% 
      filter(visit_start_date>=index_date+days(!!days_delta), visit_start_date<index_date+months(!!follow_months)) %>% 
      inner_join(all_nodes, by=c("procedure_concept_id"="node_concept_id")) %>%
      rename(node_concept_id=procedure_concept_id) %>% 
      #group_by(person_id, pasc_flag, follow_months, node_concept_id, node_concept) %>%
      #summarize(n_person_visits=n_distinct(visit_start_date)) %>%
      #    mutate(follow_months=follow_months+1) %>% 
      #mutate(visit_count=floor(1000*n_person_visits/follow_months)) %>%
      group_by(pasc_flag, node_concept_id, node_concept) %>%
      summarize(n_persons=n_distinct(person_id)) %>%
      compute_new()
    
    counts<-cohort_visits %>% 
      #      mutate(pasc_flag=as.factor(pasc_flag)) %>% 
      pivot_wider(id_cols=c("node_concept_id", "node_concept"), names_from=c("pasc_flag"), values_from=c("n_persons"), values_fill=0) %>%
      rename("positive"="1", "negative"="0") %>%
      mutate(pop_count=positive+negative, case_count=positive, control_count=negative,
             node=node_concept_id) %>%
      select(node, pop_count, case_count, control_count) %>% collect
  }
  
}



get_drug_counts<-function(tree, cohort, days_delta, method="bernoulli", follow_months=3){
  
  
  drug_tbl=cdm_tbl('drug_exposure') %>% 
    inner_join(cdm_tbl('visit_occurrence') %>% 
                 select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id")
  
  all_descendants<-tree %>% distinct(descendant_concept_name, descendant_concept_id) %>%
    mutate(node_concept=paste(descendant_concept_name, descendant_concept_id, sep=", "), node_concept_id=descendant_concept_id) %>%
    select(node_concept, node_concept_id)
  all_ancestors<-tree %>% distinct(ancestor_concept_name, ancestor_concept_id) %>%
    mutate(node_concept=paste(ancestor_concept_name, ancestor_concept_id, sep=", "), node_concept_id=ancestor_concept_id)%>%
    select(node_concept, node_concept_id)
  
  all_nodes<-dplyr::union(all_descendants, all_ancestors)
  
  if (method=="bernoulli"){
    cohort_follow_up<-cohort %>% 
      filter(index_date<as.Date("2022-09-01")-months(!!follow_months))
    
    cohort_visits<-cohort_follow_up %>% 
      inner_join(drug_tbl, by="person_id") %>% 
      filter(visit_start_date>=index_date+days(!!days_delta), visit_start_date<index_date+months(!!follow_months)) %>% 
      inner_join(all_nodes, by=c("drug_concept_id"="node_concept_id")) %>%
      rename(node_concept_id=drug_concept_id) %>% 
      #group_by(person_id, pasc_flag, follow_months, node_concept_id, node_concept) %>%
      #summarize(n_person_visits=n_distinct(visit_start_date)) %>%
      #    mutate(follow_months=follow_months+1) %>% 
      #mutate(visit_count=floor(1000*n_person_visits/follow_months)) %>%
      group_by(pasc_flag, node_concept_id, node_concept) %>%
      summarize(n_persons=n_distinct(person_id)) %>%
      compute_new()
    
    counts<-cohort_visits %>% 
      #      mutate(pasc_flag=as.factor(pasc_flag)) %>% 
      pivot_wider(id_cols=c("node_concept_id", "node_concept"), names_from=c("pasc_flag"), values_from=c("n_persons"), values_fill=0) %>%
      rename("positive"="1", "negative"="0") %>%
      mutate(pop_count=positive+negative, case_count=positive, control_count=negative,
             node=node_concept_id) %>%
      select(node, pop_count, case_count, control_count) %>% collect
  }
  
}


make_features<-function(cuts_b, all_nodes, type){
  sep_cuts(cuts=cuts_b, type=type, deduplicate=FALSE) %>%
    filter(ancestor_concept_id!=0, descendant_concept_id!=0) %>%
    left_join(all_nodes, by=c("ancestor_concept_id"="node_concept_id")) %>%
    rename(ancestor_concept=node_concept) %>%
    left_join(all_nodes, by=c("descendant_concept_id"="node_concept_id")) %>%
    rename(descendant_concept=node_concept) %>%
    distinct(ancestor_concept_id, ancestor_concept, descendant_concept_id, descendant_concept) %>%
    filter(!is.na(descendant_concept)) %>% 
    compute_new()
}

sep_cuts<-function(cuts, type="SNOMED", deduplicate=TRUE){
  if (type=="SNOMED"){
    snomed_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id=="SNOMED", domain_id=="Condition") %>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")
    snomed_tree_all<-snomed_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    #      collect
    #        
    if (deduplicate==FALSE){
      output<-cuts %>% distinct(Node.Identifier) %>% rename(ancestor_concept_id=Node.Identifier) %>% 
        inner_join(
          snomed_tree_all, by=c("ancestor_concept_id")
        ) %>% 
        select(ancestor_concept_id, descendant_concept_id) %>% 
        compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
      return(output)
    }else{
      cut_ids<-cuts %>% pull(Node.Identifier)
      
      codes_df<-data.frame(ancestor_concept_id=c(0L), descendant_concept_id=c(0L)) %>% output_tbl("codes_df", temp=TRUE)
      i=1
      for (cut_id in cut_ids){
        print(i)
        cut_desc<-snomed_tree_all %>% filter(ancestor_concept_id==cut_id) %>% compute_new(index="descendant_concept_id")
        cut_desc_new<-cut_desc %>% anti_join(codes_df, by=c("descendant_concept_id")) %>%
          distinct(ancestor_concept_id, descendant_concept_id) %>% compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
        codes_df<-codes_df %>%
          dplyr::union(cut_desc_new) %>% compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
        i=i+1
      }
      
      return(codes_df)
    }
    
    
    
  }
  if (type=="LOINC"){
    loinc_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id=="LOINC", domain_id=="Measurement") %>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")
    loinc_tree_all<-loinc_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    #         compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    if (deduplicate==FALSE){
      output<-cuts %>% distinct(Node.Identifier) %>% rename(ancestor_concept_id=Node.Identifier) %>% 
        inner_join(
          loinc_tree_all, by=c("ancestor_concept_id")
        ) %>% 
        select(ancestor_concept_id, descendant_concept_id) %>% 
        compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
      return(output)
    }else{
      cut_ids<-cuts %>% pull(Node.Identifier)
      
      codes_df<-data.frame(ancestor_concept_id=c(0L), descendant_concept_id=c(0L))
      i=1
      for (cut_id in cut_ids){
        print(i)
        cut_desc<-loinc_tree_all %>% filter(ancestor_concept_id==cut_id)
        cut_desc_new<-cut_desc %>% anti_join(codes_df, by=c("descendant_concept_id")) %>%
          distinct(ancestor_concept_id, descendant_concept_id)
        codes_df<-codes_df %>%
          rbind(cut_desc_new)
        i=i+1
      }
      
      return(codes_df)
    }
    
    
  }
  
  if (type=="proc"){
    proc_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id %in% c("ICD10PCS", "HCPCS", "CPT4"), domain_id=="Procedure") %>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")
    
    proc_tree_all<-proc_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
      #collect
      compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    if (deduplicate==FALSE){
      output<-cuts %>% distinct(Node.Identifier) %>% rename(ancestor_concept_id=Node.Identifier) %>% 
        inner_join(
          proc_tree_all, by=c("ancestor_concept_id")
        ) %>% 
        select(ancestor_concept_id, descendant_concept_id) %>% 
        compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
      return(output)
    }else{
      cut_ids<-cuts %>% pull(Node.Identifier)
      
      codes_df<-data.frame(ancestor_concept_id=c(0L), descendant_concept_id=c(0L))
      i=1
      for (cut_id in cut_ids){
        print(i)
        cut_desc<-proc_tree_all %>% filter(ancestor_concept_id==cut_id)
        cut_desc_new<-cut_desc %>% anti_join(codes_df, by=c("descendant_concept_id")) %>%
          distinct(ancestor_concept_id, descendant_concept_id)
        codes_df<-codes_df %>%
          rbind(cut_desc_new)
        i=i+1
      }
      
      return(codes_df)
    }
    
    
  }
  
  if (type=="drug"){
    drug_codes<-vocabulary_tbl("concept") %>% filter(vocabulary_id == 'RxNorm', domain_id=="Drug") %>% distinct(concept_id, concept_name) %>% compute_new(index="concept_id")    
    drug_tree_all<-drug_codes %>% rename(ancestor_concept_id=concept_id, ancestor_concept_name=concept_name) %>% 
      left_join(vocabulary_tbl("concept_ancestor"), by=c("ancestor_concept_id")) %>% 
      inner_join(vocabulary_tbl("concept"), by=c("descendant_concept_id"="concept_id")) %>% 
      select(ancestor_concept_id, ancestor_concept_name, descendant_concept_id, descendant_concept_name=concept_name) %>%
      mutate(ancestor_concept=ancestor_concept_id, descendant_concept=descendant_concept_id) %>%
#      collect
             compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
    if (deduplicate==FALSE){
      output<-cuts %>% distinct(Node.Identifier) %>% rename(ancestor_concept_id=Node.Identifier) %>% 
        inner_join(
          drug_tree_all, by=c("ancestor_concept_id")
        ) %>% 
        select(ancestor_concept_id, descendant_concept_id) %>% 
        compute_new(indices=c("ancestor_concept_id", "descendant_concept_id"))
      return(output)
    }else{
      cut_ids<-cuts %>% pull(Node.Identifier)
      
      codes_df<-data.frame(ancestor_concept_id=c(0L), descendant_concept_id=c(0L))
      i=1
      for (cut_id in cut_ids){
        print(i)
        cut_desc<-drug_tree_all %>% filter(ancestor_concept_id==cut_id)
        cut_desc_new<-cut_desc %>% anti_join(codes_df, by=c("descendant_concept_id")) %>%
          distinct(ancestor_concept_id, descendant_concept_id)
        codes_df<-codes_df %>%
          rbind(cut_desc_new)
        i=i+1
      }
      
      return(codes_df)
    }
    
  }
  
  
}

get_lab_cts<-function(cohort_tbl,
                      lab_tbl = cdm_tbl('measurement_labs'),
                      visit_tbl=cdm_tbl("visit_occurrence") %>% 
                        select(visit_occurrence_id, lab_start_date=visit_start_date),
                      days_delta=28){
  
  lab_visit_tbl<-lab_tbl %>% 
    inner_join(visit_tbl, by="visit_occurrence_id")%>%
    select(-visit_occurrence_id)
  
  cohort_tbl %>% left_join(lab_visit_tbl, by="person_id") %>%
    filter(lab_start_date>=index_date+days(days_delta)) %>%
    compute_new(index="person_id")
}



partition_dates_and_pivot_wide_lab_cts<-function(cohort, split_dates=TRUE, days_min=-28, days_max=180){
  
  if (split_dates){
    cohort_wide<-cohort %>%
      mutate(lab_window=case_when(
        lab_start_date<index_date+days(0)~"-1m_to_0m",
        lab_start_date<index_date+days(30)~"0m_to_1m",
        lab_start_date<index_date+days(60)~"1m_to_2m",
        lab_start_date<index_date+days(90)~"2m_to_3m",
        lab_start_date<index_date+days(180)~"3m_to_6m",
        lab_start_date>=index_date+days(180)~">6m"
      )) %>%
      filter(lab_window!=">6m") %>%
      group_by(person_id, index_date, observation_type, pasc_flag, lab_window, measurement_concept_id) %>%
      summarize(n_lab_visits=n_distinct(lab_start_date)) %>%
      ungroup %>%
      group_by(person_id, index_date, observation_type, pasc_flag, lab_window) %>%
      summarize(n_labs=sum(n_lab_visits)) %>% 
      ungroup %>%
      pivot_wider(id_cols=c(person_id, index_date, observation_type, pasc_flag), 
                  names_from=c(lab_window), values_from=n_labs, values_fill=0) %>%
      compute_new(index="person_id")
    
    return(cohort_wide)
  }else{
    cohort_wide<-cohort %>%
      filter(lab_start_date>index_date+days(days_min), lab_start_date<index_date+days(days_max)) %>% 
      group_by(person_id, index_date, observation_type, pasc_flag, measurement_concept_id) %>%
      summarize(n_lab_visits=n_distinct(lab_start_date)) %>%
      ungroup %>%
      group_by(person_id, index_date, observation_type, pasc_flag) %>%
      summarize(n_labs=sum(n_lab_visits)) %>% 
      ungroup %>%
      compute_new(index="person_id")
  }
}


get_feature_visits<-function(cohort,features,
                             cond_visit_tbl=cdm_tbl('condition_occurrence') %>% 
                               inner_join(cdm_tbl('visit_occurrence') %>% 
                                            select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id"),
                             days_delta=28){
  cond_feature_tbl<-cond_visit_tbl %>% 
    inner_join(features, by=c("condition_concept_id"="descendant_concept_id"))%>%
    rename(descendant_concept_id=condition_concept_id)
  
  cohort %>% inner_join(cond_feature_tbl, by="person_id") %>% 
    filter(visit_start_date>=index_date+days(days_delta)) %>%
    compute_new(index="person_id")
}


get_lab_feature_visits<-function(cohort,features,
                             lab_tbl=cdm_tbl('measurement_labs') %>% 
                               inner_join(cdm_tbl('visit_occurrence') %>% 
                                            select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id"),
                             days_delta=28){
  lab_tbl_temp<-lab_tbl %>% 
    inner_join(features, by=c("measurement_concept_id"="descendant_concept_id"))%>%
    rename(descendant_concept_id=measurement_concept_id)
  
  cohort %>% inner_join(lab_tbl_temp, by="person_id") %>% 
    filter(visit_start_date>=index_date+days(days_delta)) %>%
    compute_new(index="person_id")
}

remove_pasc_obs<-function(cohort){
  pasc_dx<-cdm_tbl("observation_derivation_recover") %>%
    filter(observation_concept_id==2000001527L, value_as_concept_id==2000001520L)
  misc_dx<-cdm_tbl("observation_derivation_recover") %>%
    filter(observation_concept_id==2000001527L, value_as_concept_id==703578L)
  pasc_or_misc_dx<-dplyr::union(pasc_dx, misc_dx)
  
  cohort %>% anti_join(pasc_or_misc_dx, by=c("condition_occurrence_id"="observation_source_concept_id")) %>%
    filter(ancestor_concept_id!=444201L) %>% 
                         compute_new(indices=c("person_id", "descendant_concept_id", "ancestor_concept_id")) %>%
                         return()
}

remove_pasc_obs_clust<-function(cohort){
  pasc_dx<-cdm_tbl("observation_derivation_recover") %>%
    filter(observation_concept_id==2000001527L, value_as_concept_id==2000001520L)
  misc_dx<-cdm_tbl("observation_derivation_recover") %>%
    filter(observation_concept_id==2000001527L, value_as_concept_id==703578L)
  pasc_or_misc_dx<-dplyr::union(pasc_dx, misc_dx)
  
  rules<-results_tbl("post_infectious_disorder_rules") %>% filter(action!="ignore")
  
  (cohort %>% anti_join(pasc_or_misc_dx, by=c("condition_occurrence_id"="observation_source_concept_id"))) %>%
    anti_join(rules, by=c("condition_concept_id", "condition_source_concept_id", "condition_source_value")) %>%
    filter(!grepl("Post-COVID", condition_source_value)) %>% 
    compute_new(indices=c("person_id")) %>%
    return()
}


get_proc_feature_visits<-function(cohort,features,
                                 proc_tbl=cdm_tbl('procedure_occurrence') %>% 
                                   inner_join(cdm_tbl('visit_occurrence') %>% 
                                                select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id"),
                                 days_delta=28){
  proc_tbl_temp<-proc_tbl %>% 
    inner_join(features, by=c("procedure_concept_id"="descendant_concept_id"))%>%
    rename(descendant_concept_id=procedure_concept_id)
  
  cohort %>% inner_join(proc_tbl_temp, by="person_id") %>% 
    filter(visit_start_date>=index_date+days(days_delta)) %>%
    compute_new(index="person_id")
}



get_drug_feature_visits<-function(cohort,features,
                                  drug_tbl=cdm_tbl('drug_exposure') %>% 
                                    inner_join(cdm_tbl('visit_occurrence') %>% 
                                                 select(visit_occurrence_id, visit_start_date), by="visit_occurrence_id"),
                                  days_delta=28){
  drug_tbl_temp<-drug_tbl %>% 
    inner_join(features, by=c("drug_concept_id"="descendant_concept_id"))%>%
    rename(descendant_concept_id=drug_concept_id)
  
  cohort %>% inner_join(drug_tbl_temp, by="person_id") %>% 
    filter(visit_start_date>=index_date+days(days_delta)) %>%
    compute_new(index="person_id")
}


make_feature_tbl<-function(cohort=results_tbl("cohort_pasc_or_misc_untested"), type, features, days_delta=-28, split_dates=TRUE, days_min=-28, days_max=180){
  if (type=="cond"){
    cohort_untested_cond_features<-get_feature_visits(cohort=cohort,
                                                      features=features, days_delta=days_delta) %>%
      remove_pasc_obs() %>% partition_dates_and_pivot_wide_cond_features(split_dates=split_dates, days_min=days_min, days_max=days_max) %>%
      compute_new(index="person_id")
    
    cohort_untested_cond_features_final<-cohort %>% #collect %>%
      distinct(person_id, index_date, pasc_flag) %>%
      full_join(cohort_untested_cond_features %>%  select(-observation_type), 
                by=c("person_id", "index_date", "pasc_flag")) %>% compute_new(index="person_id")#%>% collect
    myList <- setNames(lapply(vector("list", ncol(cohort_untested_cond_features_final)-3), function(x) x <- 0), setdiff(colnames(cohort_untested_cond_features_final), c("person_id", "index_date", "pasc_flag")))
    
    cohort_untested_cond_features_final_b<-cohort_untested_cond_features_final %>%
      replace_na(myList) %>% compute_new(index="person_id")
      
    #cohort_untested_cond_features_final[,4:ncol(cohort_untested_cond_features_final)]=cohort_untested_cond_features_final[,4:ncol(cohort_untested_cond_features_final)] %>% replace(is.na(.), 0)
    return(cohort_untested_cond_features_final_b)
    
  }
  if (type=="lab"){
    
    cohort_untested_lab_features<-get_lab_feature_visits(cohort=cohort,
                                                      features=features, days_delta=days_delta) %>%
      remove_pasc_obs() %>% partition_dates_and_pivot_wide_lab_features(split_dates=split_dates, days_min=days_min, days_max=days_max) 
    
    cohort_untested_lab_features_final<-cohort %>% #collect %>%
      distinct(person_id, index_date, pasc_flag) %>%
      full_join(cohort_untested_lab_features %>%  select(-observation_type), 
                by=c("person_id", "index_date", "pasc_flag")) #%>% compute_new(index="person_id")#%>% collect
    myList <- setNames(lapply(vector("list", ncol(cohort_untested_lab_features_final)-3), function(x) x <- 0), setdiff(colnames(cohort_untested_lab_features_final), c("person_id", "index_date", "pasc_flag")))
    
    cohort_untested_lab_features_final_b<-cohort_untested_lab_features_final %>%
      replace_na(myList) %>%compute_new(index="person_id")
    
    #cohort_untested_lab_features_final[,4:ncol(cohort_untested_lab_features_final)]=cohort_untested_lab_features_final[,4:ncol(cohort_untested_lab_features_final)] %>% replace(is.na(.), 0)
    return(cohort_untested_lab_features_final_b)
  }
  if (type=="proc"){
    cohort_untested_proc_features<-get_proc_feature_visits(cohort=cohort,
                                                     features=features, days_delta=days_delta) %>%
      remove_pasc_obs() %>% partition_dates_and_pivot_wide_proc_features(split_dates=split_dates, days_min=days_min, days_max=days_max) 
    
    cohort_untested_proc_features_final<-cohort %>% #collect %>%
      distinct(person_id, index_date, pasc_flag) %>%
      full_join(cohort_untested_proc_features %>%  select(-observation_type), 
                by=c("person_id", "index_date", "pasc_flag")) #%>% compute_new(index="person_id")#%>% collect
    myList <- setNames(lapply(vector("list", ncol(cohort_untested_proc_features_final)-3), function(x) x <- 0), setdiff(colnames(cohort_untested_proc_features_final), c("person_id", "index_date", "pasc_flag")))
    
    cohort_untested_proc_features_final_b<-cohort_untested_proc_features_final %>%
      replace_na(myList) %>%compute_new(index="person_id")
    
    #cohort_untested_proc_features_final[,4:ncol(cohort_untested_proc_features_final)]=cohort_untested_proc_features_final[,4:ncol(cohort_untested_proc_features_final)] %>% replace(is.na(.), 0)
    return(cohort_untested_proc_features_final_b)
  }
  if (type=="drug"){
    cohort_untested_drug_features<-get_drug_feature_visits(cohort=cohort,
                                                     features=features, days_delta=days_delta) %>%
      remove_pasc_obs() %>% partition_dates_and_pivot_wide_drug_features(split_dates=split_dates, days_min=days_min, days_max=days_max) 
    
    cohort_untested_drug_features_final<-cohort %>% #collect %>%
      distinct(person_id, index_date, pasc_flag) %>%
      full_join(cohort_untested_drug_features %>%  select(-observation_type), 
                by=c("person_id", "index_date", "pasc_flag")) #%>% compute_new(index="person_id")#%>% collect
    myList <- setNames(lapply(vector("list", ncol(cohort_untested_drug_features_final)-3), function(x) x <- 0), setdiff(colnames(cohort_untested_drug_features_final), c("person_id", "index_date", "pasc_flag")))
    
    cohort_untested_drug_features_final_b<-cohort_untested_drug_features_final %>%
      replace_na(myList) %>%compute_new(index="person_id")
    
    #cohort_untested_drug_features_final[,4:ncol(cohort_untested_drug_features_final)]=cohort_untested_drug_features_final[,4:ncol(cohort_untested_drug_features_final)] %>% replace(is.na(.), 0)
    return(cohort_untested_drug_features_final_b)
  }
}






colname_cleanup<-function(columns){
  
  window_transform_cond<-function(string){
    date=substr(string, 1, 8)
    if (date=="1m_to_2m"){date_new="(1 to 2 m)"}
    if (date=="2m_to_3m"){date_new="(2 to 3 m)"}
    if (date=="0m_to_1m"){date_new="(0 to 1 m)"}
    if (date=="-1m_to_0"){date_new="(-1 to 0 m)"}
    if (date_new=="(-1 to 0 m)"){rest=substr(string, 11, nchar(string))}
    else{rest=substr(string, 10, nchar(string))}
    text_temp=str_extract(rest, "^(.+?),")
    text=substr(text_temp, 1, nchar(text_temp)-1)
    string_new=paste(text, date_new, sep=" ")
    return(string_new)
  }
  window_transform_drug<-function(string){
    date=substr(string, 1, 8)
    if (date=="1m_to_2m"){date_new="(1 to 2 m)"}
    if (date=="2m_to_3m"){date_new="(2 to 3 m)"}
    if (date=="0m_to_1m"){date_new="(0 to 1 m)"}
    if (date=="-1m_to_0"){date_new="(-1 to 0 m)"}
    if (date_new=="(-1 to 0 m)"){rest=substr(string, 11, nchar(string))}
    else{rest=substr(string, 10, nchar(string))}
    text_temp=str_extract(rest, "^(.+?),")
    text_temp2=substr(text_temp, 1, nchar(text_temp)-1)
    text_temp3=str_extract(text_temp2, "_[^_]+$")
    text=str_to_sentence(substr(text_temp3, 2, nchar(text_temp3)))
    string_new=paste(text, date_new, sep=" ")
    return(string_new)
  }
  window_transform_util<-function(string){
    date=substr(string, nchar(string)-7, nchar(string))
    if (date=="1m_to_2m"){date_new="(1 to 2 m)"}
    if (date=="2m_to_3m"){date_new="(2 to 3 m)"}
    if (date=="0m_to_1m"){date_new="(0 to 1 m)"}
    if (date=="1m_to_0m"){date_new="(-1 to 0 m)"}
    if (date_new=="(-1 to 0 m)"){rest=substr(string, 1, nchar(string)-10)}
    else{rest=substr(string, 1, nchar(string)-9)}
    if (rest=="Inpatient_non_ICU_"){rest="Inpatient non-ICU"}
    if (rest=="Inpatient_ICU_"){rest="Inpatient ICU"}
    if (rest=="Inpatient_ICU"){rest="Inpatient ICU"}
    if (rest=="ED_ICU"){rest="ED ICU"}
    if (rest=="ED_non_ICU"){rest="ED non-ICU"}
    if (rest=="ED_ICU_"){rest="ED ICU"}
    if (rest=="ED_non_ICU_"){rest="ED non-ICU"}
    
    string_new=paste(rest, "visits", date_new, sep=" ")
    
    return(string_new)
  }
  
  
  sites=columns[grepl("site_", columns)]
  sex=columns[grepl("sex_", columns)]
  eth=columns[grepl("eth_", columns)]
  cohort_entry=columns[grepl("cohort_entry_", columns)]
  age=columns[columns=="entry_age"]
  cond=columns[36:1909]
  drug=columns[1910:5332]
  util=columns[5333:5361]
  lab=columns[5362:11171]
  proc=columns[11172:11931]
  lab_cts=columns[11932:11935]
  
  sites_new=c("Site H", "Site C", "Site B", "Site D", "Site I", "Site A", "Site G", "Site F")
  sex_new=c("Sex: Female", "Sex: Male")
  eth_new=c("Race/ethnicity: NH Asian/PI", "Race/ethnicity: NH Black/AA", "Race/ethnicity: Hispanic",
            "Race/ethnicity: Multiple", "Race/ethnicity: Other/unknown", "Race/ethnicity: NH White")
  cohort_entry_new=c("Cohort entry: Jan 2021",  "Cohort entry: Jan 2022",  "Cohort entry: Oct 2021", "Cohort entry: Nov 2021", "Cohort entry: Dec 2021", "Cohort entry: 2 2021", 
  "Cohort entry: Feb 2022",  "Cohort entry: Mar 2021",  "Cohort entry: Mar 2022",  "Cohort entry: Apr 2021",  "Cohort entry: Apr 2022",  "Cohort entry: May 2021", 
   "Cohort entry: May 2022",  "Cohort entry: Jun 2021",  "Cohort entry: Jun 2022",  "Cohort entry: Jul 2021",  "Cohort entry: Aug 2021",  "Cohort entry: Sep 2021" )
  age_new=c("Age")
  cond_new=as.character(sapply(cond, FUN=window_transform_cond))
  drug_new=as.character(sapply(drug, FUN=window_transform_drug))
  util_new=as.character(sapply(util, FUN=window_transform_util))
  lab_new=as.character(sapply(lab, FUN=window_transform_cond))
  proc_new=as.character(sapply(proc, FUN=window_transform_drug))
  lab_cts_new=c("Lab count (0 to 1 m)", "Lab count (-1 to 0 m)", "Lab count (1 to 2 m)", "Lab count (2 to 3 m)")
  
  columns_new=c(sites_new, sex_new, eth_new, cohort_entry_new, age_new, cond_new, drug_new, util_new, lab_new, proc_new, lab_cts_new)
  
  as.data.frame(columns_new) %>% write.csv("./specs/columns_new.csv", row.names=FALSE)
  return(columns_new)
}

