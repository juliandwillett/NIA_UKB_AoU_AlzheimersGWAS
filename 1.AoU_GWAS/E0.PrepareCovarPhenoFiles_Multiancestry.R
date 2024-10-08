# ICD codes for below. Most code obtained from the AoU example workbook for dementia. The breast cancer workbook was used as a reference for getting the family history

condition_codes_icd9_ad = c("331.0")
condition_codes_icd10_ad = c("G30.9")
condition_codes_icd9_subquery_ad=paste(sprintf('"%s"', condition_codes_icd9_ad),collapse = ", ")
condition_codes_icd10_subquery_ad=paste(sprintf('"%s"', condition_codes_icd10_ad),collapse = ", ")
query= str_glue("
SELECT
    DISTINCT person_id AS PERSON_ID,
    birth_datetime AS DATE_OF_BIRTH,
    c_race.concept_name AS RACE,
    c_sex.concept_name AS GENDER,
    c_ethn.concept_name AS ETHNICITY
FROM
    `{dataset}.person` p
    LEFT JOIN `{dataset}.concept` c_race
        ON p.race_concept_id = c_race.concept_id
    LEFT JOIN `{dataset}.concept` c_sex
        ON p.sex_at_birth_concept_id = c_sex.concept_id
    LEFT JOIN `{dataset}.concept` c_ethn
        ON p.ethnicity_concept_id = c_ethn.concept_id
")
cohort_all=dbGetQuery(con, query)

# Get early IDs
clean_cohort_all = cohort_all[,c('PERSON_ID','DATE_OF_BIRTH','RACE','GENDER','ETHNICITY')]
colnames(clean_cohort_all) <- c("Count", "Age","Race","Sex_at_Birth","Hispanic")
clean_cohort_all=clean_cohort_all %>%
mutate(Sex_at_Birth = replace(Sex_at_Birth, str_detect(Sex_at_Birth, "prefer not to answer|No matching"), "Unspecified")) %>%
mutate(Race = replace(Race, str_detect(Race, "PMI: Skip"), "Skip")) %>%
mutate(Race = replace(Race, str_detect(Race, "prefer not to answer|None of|No matching"), "Unspecified")) %>%
mutate(Hispanic = replace(Hispanic, str_detect(Hispanic, "PMI: Skip"), "Skip")) %>%
mutate(Hispanic = replace(Hispanic, str_detect(Hispanic, "Not To Answer|None Of|No matching"), "Unspecified"))

# Get age
library(eeptools)
clean_cohort_all$Age=floor(age_calc(as.Date(clean_cohort_all$Age),Sys.Date(),units = "years")) # default is age in months

df = data.frame(IID=clean_cohort_all$Count) %>% mutate(AD=ifelse(IID %in% cohort_PIDs_ad$person_id,1,0))
merged2 = merge(df,clean_cohort_all %>% rename(IID=Count),by='IID',all.x = T)

# Get the FHx of dementia
library(bigrquery)
library(questionr)
(CURRENT_DATASET <- Sys.getenv("WORKSPACE_CDR"))
my_bucket = Sys.getenv('WORKSPACE_BUCKET')

download_data <- function(query) {
    tb <- bq_project_query(Sys.getenv('GOOGLE_PROJECT'), query)
    bq_table_download(tb)
}
codes_ad <- paste("\"",c('NervousCondition_Dementia_yes'),"\"", collapse = ",", sep="")
query_ad <- sprintf("SELECT person_id, observation_source_value, observation_source_concept_id, value_source_concept_id,
                                  question.concept_class_id as qtype, question.concept_name as qname,
                                  answer.concept_class_id as atype, answer.concept_name as aname
                           FROM `%s.observation` 
                           LEFT JOIN `%s.concept` question on (question.concept_id=observation_source_concept_id)
                           LEFT JOIN `%s.concept` answer on (answer.concept_id=value_source_concept_id)
                           WHERE observation_source_value IN (%s)", 
                      CURRENT_DATASET, CURRENT_DATASET, CURRENT_DATASET, codes_ad)
fhx_dementia_all <- download_data(query_ad)
no_skips_fhx = fhx_dementia_all %>% filter(!str_detect(aname,'Skip'))
merged3 = merge(merged2,no_skips_fhx %>% rename(IID=person_id),by='IID',all.x = T)
ids_by_proxy = (merged3)$IID

# Start making final data frame
full_df = merged3 %>% select(-Race,-Hispanic,-observation_source_value,-observation_source_concept_id,
                  -value_source_concept_id,-qtype,-qname,-atype) %>%
rename(AD_by_proxy = aname,Sex = Sex_at_Birth) %>%
filter(!is.na(Age)) %>% 
mutate(Sex = ifelse(Sex == 'Female',0,
                   ifelse(Sex == 'Male',1,
                         ifelse(Sex == 'PMI: Skip',2,
                               ifelse(Sex == 'Unspecified',2,
                                     ifelse(Sex == 'None',2,3)))))) %>%
mutate(AD_by_proxy = ifelse(IID %in% ids_by_proxy,1,0)) %>% distinct()

# Bring in PC data
pcs = vroom('array_autosomes_postqc_pruned_onlysrwgs_pca_results_plink.eigenvec.eigenvec') %>% 
rename(IID = `#IID`)
full_df_pcs = merge(full_df,pcs,by="IID")

# Now make individual covar and pheno data frames
covar_df = data.frame(IID = full_df_pcs$IID,Age = full_df_pcs$Age,Sex = full_df_pcs$Sex,
                      PC1 = full_df_pcs$PC1,PC2 = full_df_pcs$PC2,PC3 = full_df_pcs$PC3,
                     PC4 = full_df_pcs$PC4, PC5 = full_df_pcs$PC5, PC6 = full_df_pcs$PC6,
                     PC7 = full_df_pcs$PC7, PC8 = full_df_pcs$PC8, PC9 = full_df_pcs$PC9,
                     PC10 = full_df_pcs$PC10, PC11 = full_df_pcs$PC11, PC12 = full_df_pcs$PC12,
                     PC13 = full_df_pcs$PC13, PC14 = full_df_pcs$PC14, PC15 = full_df_pcs$PC15,
                     PC16 = full_df_pcs$PC16, PC17 = full_df_pcs$PC17, PC18 = full_df_pcs$PC18,
                     PC19 = full_df_pcs$PC19, PC20 = full_df_pcs$PC20) %>%
mutate(FID = 0, .before = 'IID')

pheno_df_all = data.frame(IID = full_df_pcs$IID, AD = full_df_pcs$AD, 
                     Dementia_By_Proxy = full_df_pcs$AD_by_proxy) %>%
mutate(AD_any = ifelse(AD == 1 | Dementia_By_Proxy == 1,1,0)) %>%
mutate(FID = 0, .before = 'IID')
pheno_df = pheno_df_all %>% select(-Dementia,-Dementia_any,-Dementia_By_Proxy)

# Write the files
vroom_write(covar_df,'regenie_covar_revised.txt')
vroom_write(pheno_df,'regenie_pheno_revised.txt')

# Backup
system("gsutil -m cp -rn regenie_* $WORKSPACE_BUCKET/data/")
