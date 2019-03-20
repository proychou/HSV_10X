#10x script
#Pavitra Roychoudhury
#Aug 2018
# sacct --format=JobID,Start,End,AllocCPUS,MaxRSS,Elapsed,State,Nodelist

cd /fh/fast/jerome_k/HSV_10x
module load cellranger/2.1.1

#Genomics core generated demultiplexed fastqs using cellranger
#These are at /fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/

#There are 4 libraries:
# Tgmixed_AAV
# SCGmixed_AAV
# NoAAV_TG
# NoAAV_SCG

#Run cellranger via grabnode -- this failed when the connection dropped. Boo!
#lib_dir='/fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/NoAAV_TG'
# cellranger count --id NoAAV_TG --fastqs $lib_dir --sample NoAAV_TG --transcriptome /fh/fast/jerome_k/HSV_10x/refdata-cellranger-hg19_and_mm10-1.2.0/

#lib_dir='/fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/Tgmixed_AAV'
# cellranger count --id Tgmixed_AAV --fastqs $lib_dir --sample Tgmixed_AAV --transcriptome /fh/fast/jerome_k/HSV_10x/refdata-cellranger-hg19_and_mm10-1.2.0/


#1. Using sbatch: NoAAV_TG: 28163820
lib_dir='/fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/NoAAV_TG'

sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_cellranger.sh -i NoAAV_TG -d $lib_dir -s NoAAV_TG

#2. Using sbatch: Tgmixed_AAV: 28164223
lib_dir='/fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/Tgmixed_AAV'

sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_cellranger.sh -i Tgmixed_AAV -d $lib_dir -s Tgmixed_AAV

#3. Using sbatch: NoAAV_SCG: 28164208
lib_dir='/fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/NoAAV_SCG'

sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_cellranger.sh -i NoAAV_SCG -d $lib_dir -s NoAAV_SCG

#4. Using sbatch: SCGmixed_AAV: 28164209
lib_dir='/fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/SCGmixed_AAV'

sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_cellranger.sh -i SCGmixed_AAV -d $lib_dir -s SCGmixed_AAV 




#########################################################################################
### Part II: Mapping reads to AAV and HSV
cd /fh/fast/jerome_k/HSV_10x

#SCGmixed_AAV
for f2 in `ls /fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/SCGmixed_AAV/*_R2_001.fastq.gz`; do

f1=`echo $f2| sed s/_R2/_R1/`; 
sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_bbduk_filter.sh -1 $f1 -2 $f2 -s 'SCGmixed_AAV'

done

#Tgmixed_AAV
for f2 in `ls /fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/Tgmixed_AAV/*_R2_001.fastq.gz`; do

f1=`echo $f2| sed s/_R2/_R1/`; 
sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_bbduk_filter.sh -1 $f1 -2 $f2 -s 'Tgmixed_AAV'

done

#NoAAV_TG
for f2 in `ls /fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/NoAAV_TG/*_R2_001.fastq.gz`; do

f1=`echo $f2| sed s/_R2/_R1/`; 
sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_bbduk_filter.sh -1 $f1 -2 $f2 -s 'NoAAV_TG'

done

#NoAAV_SCG
for f2 in `ls /fh/fast/jerome_k/SR/ngs/illumina/dstrongi/180824_D00300_0597_BCCFCEANXX/analysis/CCFCEANXX/outs/fastq_path/CCFCEANXX/NoAAV_SCG/*_R2_001.fastq.gz`; do

f1=`echo $f2| sed s/_R2/_R1/`; 
sbatch -t 2-0 -c 14 --mem=256000 -p largenode run_bbduk_filter.sh -1 $f1 -2 $f2 -s 'NoAAV_SCG'

done

#########################################################################################
### Part IIb: Mapping reads to AAV and HSV
## Method 2: Make alternate reference that includes our sequences 
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references

cd /fh/fast/jerome_k/HSV_10x
module load cellranger/2.1.1

cellranger mkref --genome=mm10 --fasta=/fh/fast/_CTR/BioLib/ref/cellranger/refdata-cellranger-mm10-1.2.0/fasta/genome.fa --genes=/fh/fast/_CTR/BioLib/ref/cellranger/refdata-cellranger-mm10-1.2.0/genes/genes.gtf --genome=HSV1 --fasta=/fh/fast/jerome_k/HSV_10x/refs/NC_001806.fa --genes=/fh/fast/jerome_k/HSV_10x/refs/NC_001806.gtf --genome=AAV8-NLS-mEGFP --fasta=/fh/fast/jerome_k/HSV_10x/refs/AAV8-NLS-mEGFP.fa --genes=/fh/fast/jerome_k/HSV_10x/refs/AAV8-NLS-mEGFP.gtf --genome=AAV1-NLS-mScarlet --fasta=/fh/fast/jerome_k/HSV_10x/refs/AAV1-NLS-mScarlet.fa --genes=/fh/fast/jerome_k/HSV_10x/refs/AAV1-NLS-mScarlet.gtf --genome=AAV_PHP_S_NLS-DsRed-Express2 --fasta=/fh/fast/jerome_k/HSV_10x/refs/AAV-PHP_S_NLS-DsRed-Express2.fa --genes=/fh/fast/jerome_k/HSV_10x/refs/AAV-PHP_S-NLS-DsRed-Express2.gtf --genome=AAV-Rh10-NLS-mTagBFP2 --fasta=/fh/fast/jerome_k/HSV_10x/refs/AAV-Rh10-NLS-mTagBFP2.fa --genes=/fh/fast/jerome_k/HSV_10x/refs/AAV-Rh10-NLS-mTagBFP2.gtf



#########################################################################################
### Part III. Run aggregator to merge libraries
cd /fh/fast/jerome_k/HSV_10x
module load cellranger/2.1.1

# cellranger aggr --id=hsv10x_agg \
#                   --csv='./hsv10x_libraries.csv' \
#                   --normalize=mapped
                  
                  
cellranger aggr --id=hsv10x_agg_nocontrols \
                  --csv='./hsv10x_libraries_removecontrols.csv' \
                  --normalize=mapped

##########################################################################################
### Part IV. Re-run tnse on subset of genes

cd /fh/fast/jerome_k/HSV_10x
module load cellranger/2.1.1

cellranger reanalyze --id=subset252genes_nocontrols --matrix=/fh/fast/jerome_k/HSV_10x/hsv10x_agg_nocontrols/subset252genes_nocontrols.h5






