# module load vcftools/0.1.16-xtum
module load bcftools/1.19-openblas-2333
# module load samtools/1.19.2-l4co
# module load bwa/0.7.17-pf24

# VCF files are downloaded from TOPMed Freeze 8
# https://legacy.bravo.sph.umich.edu/freeze8/hg38/downloads

for chr in {1..22} X; do
    if [[ "$chr" == "X" ]]; then
        vcf_file="chrX.BRAVO_TOPMed_Freeze_8.vcf.gz"
        temp_file="temp_23.snplist"
    else
        vcf_file="chr${chr}.BRAVO_TOPMed_Freeze_8.vcf.gz"
        temp_file="temp_${chr}.snplist"
    fi

    # extract CHROM, POS, REF, ALT, AC, AF，change chrX to 23
    bcftools query -i 'FILTER="PASS"' -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AF\n' $vcf_file | sed 's/^chr//' | sed 's/^X/23/' > $temp_file

    snp_count_filtered=$(wc -l < $temp_file)
    snp_count_vcf=$(bcftools query -i 'FILTER="PASS"' -f '%CHROM:%POS\n' $vcf_file | wc -l)

    if [ "$snp_count_filtered" -eq "$snp_count_vcf" ]; then
        echo "chr${chr} SNP list is complete"
        # cat $temp_file >> $output_file  # merge SNP list
    else
        echo "chr${chr} SNP list is incomplete. Exiting..."
        exit 1  # if the SNP list is incomplete, exit the script
    fi
done

echo "All SNP lists have been merged into $output_file."

rm -f temp*

output_file="pass_unsorted.snplist"
> $output_file
echo -e "CHR\tPOS\tREF\tALT\tAC\tAF" > $output_file
cat temp_{1..23}.snplist >> "$output_file"

# filter maf >= 0.001 and maf <= 0.999
awk -F'\t' 'NR==1 || ($5 >= 5 && $6 >= 0.01 && $6 <= 0.999)' pass_unsorted.snplist > pass_filtered.snplist

Rscript ../script/topmed38snplist.R

gzip topmed.GRCh38.f8wgs.pass.mac5.maf001.tab.snplist

cp topmed.GRCh38.f8wgs.pass.mac5.maf001.tab.snplist.gz /projects/MRC-IEU/research/projects/ieu3/p3/022/working/data/

# liftOver 38 to 37 chain file
# wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
# liftOver 37 to 38 chain file
# wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

# liftOver chain files ieu space
# /mnt/storage/private/mrcieu/data/genomic_data/liftover/

# wget https://ftp.ncbi.nlm.nih.gov/snp/archive/b156/VCF/GCF_000001405.25.gz
# wget https://ftp.ncbi.nlm.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz

# after submit run_topmed38to37.bash
awk 'FNR==1 && NR!=1 {next} {print}' temp_*_liftoverto37.snplist > topmed.GRCh37.f8wgs.pass.mac5.maf001.tab.snplist
{ head -n 1 merged_file.snplist; tail -n +2 merged_file.snplist | sort -k1,1; } > merged_file_sorted.snplist
# then run R script merge_topmedhrc.R