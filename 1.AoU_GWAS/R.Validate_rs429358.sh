# using regenie_input/regenie_covar_eur.txt as the covar file. Just running step 2
# detailed info on AD by proxy in: regenie_input/regenie_covar_pheno_no_pcs.txt

# AD ICD phenotype
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/rs429358 \
                    --phenoFile regenie_input/regenie_pheno_adonly_notproxy.txt \
                    --covarFile regenie_input/regenie_covar_eur.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1_singleanc/aou_step1_rg_common_eur_pred.list --ignore-pred \
                    --bsize 400 \
                    --out test/rs429358_icd_ad_remove_byproxy_otherwise \
                    --minMAC 20 \
                    --phenoCol AD

# Grandparents not cases
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/rs429358 \
                    --phenoFile regenie_input/regenie_pheno_adbyproxy_nograndparent.txt \
                    --covarFile regenie_input/regenie_covar_eur.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1_singleanc/aou_step1_rg_common_eur_pred.list --ignore-pred \
                    --bsize 400 \
                    --out test/rs429358_ad_by_proxy_no_grandparent \
                    --minMAC 20 \
                    --phenoCol AD_by_proxy

# Grandparents not cases, controls at least age 60
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/rs429358 \
                    --phenoFile regenie_input/regenie_pheno_adbyproxy_nograndparent_controlsover60.txt \
                    --covarFile regenie_input/regenie_covar_eur.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1_singleanc/aou_step1_rg_common_eur_pred.list --ignore-pred \
                    --bsize 400 \
                    --out test/rs429358_ad_by_proxy_no_grandparent_controlsover60 \
                    --minMAC 20 \
                    --phenoCol AD_by_proxy

# Grandparents by proxy cases only
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/rs429358 \
                    --phenoFile regenie_input/regenie_pheno_adbyproxy_onlygrandparent_proxy.txt \
                    --covarFile regenie_input/regenie_covar_eur.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1_singleanc/aou_step1_rg_common_eur_pred.list --ignore-pred \
                    --bsize 400 \
                    --out test/rs429358_ad_by_proxy_grandparentonly \
                    --minMAC 20 \
                    --phenoCol AD_by_proxy

# Grandparents by proxy cases only, controls over 60
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/rs429358 \
                    --phenoFile regenie_input/regenie_pheno_adbyproxy_onlygrandparent_proxy_controls_over60.txt \
                    --covarFile regenie_input/regenie_covar_eur.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1_singleanc/aou_step1_rg_common_eur_pred.list --ignore-pred \
                    --bsize 400 \
                    --out test/rs429358_ad_by_proxy_grandparentonly_controlsover60 \
                    --minMAC 20 \
                    --phenoCol AD_by_proxy
