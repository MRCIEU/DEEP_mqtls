import onnx
import hail as hl
import random
import pandas as pd
import matplotlib.pyplot as plt
from gnomad.sample_qc.ancestry import (
    apply_onnx_classification_model,
    assign_population_pcs,
    pc_project,
)
from gnomad.utils.filtering import filter_to_adj
import sys

logfile = sys.argv[1]
vcf_file = sys.argv[2]

logfile = "/user/work/er20212/test/data/ancestry_infer_hail.log"

hl.init(backend = "spark", # use local
    local="local[*]",  # use as many cores as available
    min_block_size=128,  # minimum block size for Hail
    log=logfile) # log file for Hail

if "${genome_build}" == "37":
    reference_genome = 'GRCh37'
    dat = hl.import_vcf(vcf_file, reference_genome=reference_genome)
    chain_file = f"{scripts_directory}/resources/genetics/references_grch37_to_grch38.over.chain.gz"
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover(chain_file, rg38)
    dat = dat.annotate_rows(new_locus = hl.liftover(mt.locus, 'GRCh38'))
    dat = dat.filter_rows(hl.is_defined(dat.locus))

    dat = dat.key_rows_by(
        locus = dat.new_locus,
        alleles = dat.alleles)
    dat = dat.drop('new_locus')

elif "${genome_build}" == "38":
    reference_genome = 'GRCh38'
    dat = hl.import_vcf(vcf_file, reference_genome=reference_genome)

pca_loadings = hl.read_table(f'${scripts_directory}/resources/genetics/hgdp_tgp_pca_gbmi_snps_loadings.GRCh38.ht')

dat = dat.filter_rows(hl.is_defined(pca_loadings[mt38.locus, mt38.alleles]))

ht_projections = pc_project(dat, pca_loadings)
ht_projections = ht_projections.transmute(**{f"PC{i}": ht_projections.scores[i - 1] for i in range(1, 21)})

ht_projections = ht_projections.key_by()
ht_projections = ht_projections.select(
    **{"IID": hl.str(ht_projections.key)}, 
    **{f"PC{i}": ht_projections[f"PC{i}"] for i in range(1, 21)}
)

my_df = ht_projections.to_pandas()

ref = pd.read_csv(f'{scripts_directory}/resources/genetics/hgdp_tgp_unrel_pass_pca_gbmi_score.tsv', sep='\t')
print(ref.columns)

with gzip.open('/user/work/er20212/test/data/release_3.1.2_vcf_genomes_gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz', 'rt') as f:
    meta_df = pd.read_csv(f, sep='\t')

meta_df_filt = meta_df[meta_df['s'].isin(ref['s'])]
meta_subset = meta_df_filt[['s', 'hgdp_tgp_meta']]

def extract_genetic_region(hgdp_meta):
    """from hgdp_tgp_meta dictionary, extract genetic_region"""
    if pd.isna(hgdp_meta):
        return None
    try:
        if isinstance(hgdp_meta, str):
            meta_dict = json.loads(hgdp_meta)
        else:
            meta_dict = hgdp_meta
        return meta_dict.get('genetic_region', None)
    except (json.JSONDecodeError, AttributeError, TypeError):
        return None

meta_subset['genetic_region'] = meta_subset['hgdp_tgp_meta'].apply(extract_genetic_region)

print(meta_subset['genetic_region'].value_counts())

hgdp_tgp_ref = pd.merge(ref, meta_subset[['s', 'genetic_region']], on='s', how='left')

def extract_global_pca_scores(hgdp_meta):
    """from hgdp_tgp_meta dictionary, extract global_pca_scores"""
    if pd.isna(hgdp_meta):
        return None, None
    try:
        if isinstance(hgdp_meta, str):
            meta_dict = json.loads(hgdp_meta)
        else:
            meta_dict = hgdp_meta
        
        pca_scores = meta_dict.get('global_pca_scores', None)
        if pca_scores and len(pca_scores) >= 2:
            return pca_scores[0], pca_scores[1]  # PC1, PC2
        else:
            return None, None
    except (json.JSONDecodeError, AttributeError, TypeError):
        return None, None

pca_results = meta_df['hgdp_tgp_meta'].apply(extract_global_pca_scores)
meta_df['PC1'] = [result[0] for result in pca_results]
meta_df['PC2'] = [result[1] for result in pca_results]
meta_df['genetic_region'] = meta_df['hgdp_tgp_meta'].apply(extract_genetic_region)
meta_df = meta_df[1:]

plt.figure(figsize=(8, 6))
for pop in meta_df['genetic_region'].unique():
    sub = meta_df[meta_df['genetic_region'] == pop]
    plt.scatter(sub['PC1'], sub['PC2'], label=pop, s=10, alpha=0.5)

plt.scatter(my_df['PC1'], my_df['PC2'],
            color='black', marker='x', s=10, label='bib mother')

plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Projection of Samples onto gnomAD PCA space')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig(f${}'.png', dpi=300, bbox_inches='tight')