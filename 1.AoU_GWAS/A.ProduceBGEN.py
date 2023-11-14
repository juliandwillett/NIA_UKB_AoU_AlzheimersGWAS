# set up environment
import hail as hl
import os
import sys
import subprocess
from datetime import datetime
my_bucket = os.getenv('WORKSPACE_BUCKET')
hl.init(default_reference='GRCh38', idempotent=True,tmp_dir = f"gs://fc-secure-f509fd3c-7cd4-49e0-b409-fc98ab38bb4f/tmp6/") # tmpdir to boost stability with preemptible workers

# get flagged samples for QC
flagged_samples_path = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv"
flagged_samples = hl.import_table(flagged_samples_path) # key already set to s
flagged_samples = flagged_samples.annotate(sample_id = flagged_samples.s).key_by('sample_id')

# isolate single chromosome, to enable parallelization and make failed jobs (due to insufficient RAM/HDD) less costly
curr_chr = "chr19" 
vds_path = os.getenv("WGS_VDS_PATH")
vds = hl.vds.read_vds(vds_path)
vds = hl.vds.filter_chromosomes(vds, keep = [curr_chr])
vds = hl.vds.filter_samples(vds, flagged_samples, keep = False, remove_dead_alleles = True) #remove flagged samples

# do additional processing steps suggested in AoU notebook on manipulating VDS files, then write out
mt = vds.variant_data.annotate_entries(AD = hl.vds.local_to_global(vds.variant_data.LAD,
                                                                   vds.variant_data.LA,
                                                                   n_alleles=hl.len(vds.variant_data.alleles),
                                                                   fill_value=0, number='R'))
mt = mt.annotate_entries(GT = hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
mt = mt.transmute_entries(FT = hl.if_else(mt.FT, "PASS", "FAIL"))
mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, mt))
mt = mt.annotate_rows(info = hl.agg.call_stats(mt.GT, mt.alleles))
mt = mt.drop('as_vqsr','LAD', 'LGT', 'LA',
            'tranche_data', 'truth_sensitivity_snp_threshold', 
             'truth_sensitivity_indel_threshold','snp_vqslod_threshold','indel_vqslod_threshold')

# split multiallelics
mt = hl.split_multi_hts(mt) # use split_multi_hts (not split_multi as in AoU example code) to properly update other columns (per Hail developers)
out_bgen = f'{my_bucket}/data/bgen_minimal_qc/plink_{curr_chr}_multi_split'
gp_values = hl.literal([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

# remove flagged variants
mt = mt.filter_rows(hl.len(mt.filters)==0)
mt = mt.transmute_entries(FT = hl.or_else(mt.FT, "PASS")) # put PASS to missing entries, per QC doc
mt = mt.filter_entries(mt.FT == "PASS") # should keep PASS and missing entries, cutting FAIL results, but be sure to verify

# remove monomorphic variants
mt = hl.methods.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.AF[1] > 0) & (mt.variant_qc.AF[1] < 1)) # remove monomorphic variants

# now write to bgen. Use parallel to address memory concerns. QC and bgen file production took: 20 minutes.
print(str(datetime.now()))
print('Writing bgen file')
hl.export_bgen(mt, out_bgen, gp=gp_values[mt.GT.n_alt_alleles()], parallel='header_per_shard')
print(str(datetime.now()))
