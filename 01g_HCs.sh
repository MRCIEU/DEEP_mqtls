# because PCs are constructed from total variance explained, many are associated with genetically localized structures rather than population structure.
# computing haplotype components (HCs) through PBWTpaint
source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_01g_logfile})
print_version

# vcf files needs to be cleaned with SNPs with samples and variants from 01b
# input vcf from input_data folder

# needs number of SNPs on each chromosome
# needs number of individuals
# input to combine file

# sample id order from vcf

for i in {1..22}; do

    echo "Computing Haplotype Components for chr${i}..."
    ${pbwt} -readVcfGT chr${i}_data.vcf.gz -paintSparse chr${i}_alnall 100 2 500

done