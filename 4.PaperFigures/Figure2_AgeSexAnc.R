# To respect AoU data privacy, had to bring the UKB regenie covar/pheno file onto AoU. I de-identified subjects, as we do not need that information for this figure.

getUKBDataforFigure2 = function() {
  ukb_samples = vroom("ukb_psam_wgs.psam")
  ukb_covar = vroom("~/true_lab_storage/Data_Links/UKB_Covar_PhenoFiles/covar2") %>% 
    filter(IID %in% ukb_samples$IID)
  ukb_pheno = vroom("~/true_lab_storage/Data_Links/UKB_Covar_PhenoFiles/pheno2") %>% 
    filter(IID %in% ukb_samples$IID)
  ukb_anc = vroom("~/true_lab_storage/Data_Links/UKB_Covar_PhenoFiles/keep_white.tsv") %>% 
    filter(IID %in% ukb_samples$IID)
  
  table_for_aou = (merge(ukb_covar,ukb_pheno,by="IID") %>% mutate(White=(IID %in% ukb_anc$IID)) %>%
    select(-IID,-FID.x,-FID.y))
  vroom_write(table_for_aou[-c(3:23),],"ukb_covar_pheno_for_fig_2.txt")
}

#########
# Next pull this data onto AoU, then start running the following R code to organize data
# For personal reference, this code is in the following AoU notebook: duplicateofalzheimersgwastake5/notebooks/Part5_GetMetrics.ipynb
library(vroom) ; library(glue) ; library(ggplot2) ; library(tidyverse) ; library(viridis)

# Get AoU covar/pheno/ancestry data
system("gsutil cp -r -n bucket/data/regenie_pheno.txt .")
system("gsutil cp -r -n bucket/data/regenie_covar.txt .")
system("gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv .")
pheno_file = vroom('regenie_pheno.txt',show_col_types = FALSE)
covar_file = vroom('regenie_covar.txt',show_col_types = FALSE)
anc_file = vroom("ancestry_preds.tsv",show_col_types = FALSE) %>% rename(IID = research_id) %>% 
select(IID,ancestry_pred)

# Merge files to make it easier to work with
tmp_merged_df = merge(pheno_file,covar_file,by="IID")
merged_df = merge(tmp_merged_df,anc_file,by="IID")

# Get numbers of cases
print(glue("Number of Subjects with AD ICD codes: {length(which(merged_df$AD == 1))}"))
print(glue("Number of Subjects with AD-by-proxy: {length(which(merged_df$Dementia_By_Proxy == 1))}"))
print(glue("Number of Subjects with Either: {length(which(merged_df$Dementia_By_Proxy == 1 | merged_df$AD == 1))}"))

# Pull the ukb data
ukb_vars = vroom("ukb_covar_pheno_for_fig_2.txt")

##########
### Make the actual figures
## Age plot:
aou_ages = merged_df %>% mutate(AD_any = as.logical(AD_any),
                             AD_any = ifelse(AD_any,"Case","Control")) %>%
select(AD_any,Age,Sex) %>% mutate(Cohort = "AoU")

ukb_ages = ukb_vars %>% mutate(AD_any = ifelse(ad_proxy > 0,"Case","Control")) %>%
select(AD_any,Age_at_recruitment,Sex) %>% rename(Age=Age_at_recruitment) %>% 
mutate(Cohort = "UKB")

# Make plot
plt_df = aou_ages %>% add_row(ukb_ages)
ggplot(plt_df,aes(x=AD_any,y=Age,fill = Cohort)) + 
geom_violin(color = "blue", alpha = 0.7) +
  geom_boxplot(width=0.3,position=position_dodge(0.9)) +
  labs(x = "AD Status", y = "Age") +
  theme_minimal() + theme(axis.text = element_text(size = 16),   
    axis.title = element_text(size = 20),  
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16))

## Sex distribution
# Sex codes: Female 0, Male 1, Skip 2, Unspecified 2, None 2, Other 3
aou_perc_female_case = length(which(aou_ages$AD_any == "Case" & aou_ages$Sex == 0)) / length(which(aou_ages$AD_any == "Case"))
aou_perc_female_control = length(which(aou_ages$AD_any == "Control" & aou_ages$Sex == 0)) / length(which(aou_ages$AD_any == "Control"))
aou_new_df = data.frame(Group = c("Case","Control"), PercFemale = c(aou_perc_female_case,aou_perc_female_control)) %>%
mutate(Cohort = "AoU")

ukb_perc_female_case = length(which(ukb_ages$AD_any == "Case" & ukb_ages$Sex == 0)) / length(which(ukb_ages$AD_any == "Case"))
ukb_perc_female_control = length(which(ukb_ages$AD_any == "Control" & ukb_ages$Sex == 0)) / length(which(ukb_ages$AD_any == "Control"))
ukb_new_df = data.frame(Group = c("Case","Control"), PercFemale = c(ukb_perc_female_case,ukb_perc_female_control)) %>%
mutate(Cohort = "UKB")

# Make plot

plt_df = aou_new_df %>% add_row(ukb_new_df)
ggplot(plt_df,aes(x=Group,y=PercFemale,fill=Cohort)) + 
geom_bar(position = "dodge",stat = "identity", color = "blue",alpha = 0.7) +
  labs(x = "AD Status", y = "Percentage Female") +
  theme_minimal() + theme(axis.text = element_text(size = 16),   
    axis.title = element_text(size = 20),  
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16))

## Ancestry distribution
aou_perc_eur_case = length(which(merged_df$AD_any == 1 & merged_df$ancestry_pred == "eur")) / length(which(merged_df$AD_any == 1))
aou_perc_eur_control = length(which(merged_df$AD_any == 0 & merged_df$ancestry_pred == "eur")) / length(which(merged_df$AD_any == 0))
aou_new_df = data.frame(Group = c("Case","Control"), PercEUR = c(aou_perc_eur_case,aou_perc_eur_control)) %>%
mutate(Cohort = "AoU")

ukb_perc_eur_case = length(which(ukb_vars$ad_proxy == 1 & ukb_vars$White)) / length(which(ukb_vars$ad_proxy == 1))
ukb_perc_eur_control = length(which(ukb_vars$ad_proxy == 0 & ukb_vars$White)) / length(which(ukb_vars$ad_proxy == 0))
ukb_new_df = data.frame(Group = c("Case","Control"), PercEUR = c(ukb_perc_eur_case,ukb_perc_eur_control)) %>%
mutate(Cohort = "UKB")

# Make plot
plt_df = aou_new_df %>% add_row(ukb_new_df)
ggplot(plt_df,aes(x=Group,y=PercEUR,fill=Cohort)) + 
geom_bar(position = "dodge",stat = "identity", color = "blue",alpha = 0.7) +
  labs(x = "AD Status", y = "Percentage EUR") +
  theme_minimal() + theme(axis.text = element_text(size = 16),   
    axis.title = element_text(size = 20),  
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16))
