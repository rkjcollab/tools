

# SDS modifying from SER's pipeline in ancestry_estimation.qmd.
# Step 1 of Admixture for ancestry estimation.
# NOTE: this script only needs to be run once.

# First we need to download the 1000 genomes genetic data and the demographic data.
# The genetic data is downloaded by chromosome.
# There is also a file of well established regions to subset on provided by 1000 Genomes
# where the quality of reading is up to par. We will first restrict to the well-genotyped
# regions before QC.

# get wd
DOWNLOAD_DIR="${RKJCOLLAB}/Immunogenetics_T1D/raw/genetic_maps/1000genomes-phase3"

# get the genetic files from 1000 genomes
echo "Get files from 1000 genomes..."
for chr in $(seq 1 22) 
do
	echo "Processing chr${chr}..."
	wget -P "${DOWNLOAD_DIR}" https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	wget -P "${DOWNLOAD_DIR}" https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
done

# get ancestry info
wget -P "${DOWNLOAD_DIR}" https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# strict mask .bed (extract only SNPs in the strict mask file)
wget -P "${DOWNLOAD_DIR}" https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed

echo "DONE..."