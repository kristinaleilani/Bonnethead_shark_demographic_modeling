# Script is based on: https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh
# Reads were processed using the high-performance computing resources of the Texas Advanced Computing Center (TACC). 
# Necessary packages include: Cutadapt, Bowtie2, ANGSD, ngsRelate, GADMA, ngsF, and 2bRAD processing scripts from https://github.com/z0on/2bRAD_denovo/tree/master


#####----------- PREPROCESSING READS ---------###############

# find number of reads per sample
>reads
for file in *.fastq; do
echo "$file $(cat $file|wc -l)/4|bc" >> reads;
done

# for reference-based analysis, trim poor quality bases off ends:
>trimse
for file in *.fastq; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.fastq/}.trim $file > ${file}_trimlog.txt" >> trimse;
done



#####----------- IMPORT REFERENCE GENOME ---------###############

# import the great hammerhead shark reference genome (i.e., sphyrna mokarran) from NCBI
# Assembly: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_024679065.1/
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/679/065/GCA_024679065.1_ASM2467906v1/GCA_024679065.1_ASM2467906v1_genomic.fna.gz .
gunzip GCA_024679065.1_ASM2467906v1_genomic.fna.gz



#####------------ INDEX GENOME ---------###############

export GENOME_FASTA=GCA_024679065.1_ASM2467906v1_genomic.fna
echo "bowtie2-build $GENOME_FASTA $GENOME_FASTA" >btb2
samtools faidx $GENOME_FASTA




#####------------ MAPPING TO GENOME -----------------###############
# Mapping your reads to the genome, converting to bams, indexing
export GENOME_FASTA=GCA_024679065.1_ASM2467906v1_genomic.fna
>maps
for F in `ls *.trim`; do
echo "bowtie2 --local --no-unal -x $GENOME_FASTA -U ${F} -S ${F/.trim/}.sam && samtools sort -O bam -o ${F/.trim/}.bam ${F/.trim/}.sam && samtools index ${F/.trim/}.bam">>maps
done


# Get mapping stats
>alignmentRates
for F in `ls *trim`; do 
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.e* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "${F/.trim/}.bam $M">>alignmentRates;
done
sort -k 2 -n alignmentRates > alignmentRates_sorted





#####-------------- GENOTYPING SAMPLES --------------------###############
ls *.bam > bams
# quality assessment, removing bams with log(coverage)<3SD
# also imposing minimum number of individuals(MI) a locus must be seen in (genotyping rate cutoff - 50%)
# Use ANGSD installed locally (not conda version)
module load Rstats
export MinIndPerc=0.50
FILTERSQ="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1000 -minInd $MI"
echo 'export NIND=`cat bams | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo "source calc && angsd -b bams -GL 1 $FILTERSQ $TODOQ -P 1 -out dd && Rscript ~/bin/plotQC.R prefix=dd >qualRanks">a0
# look at quality of reads in dd.pdf, check how many bams retained in bams.qc
# Removed TB127.bam because it is most likely a different species


# Initial IBS production
export MinIndPerc=0.75
FILTERS1='-minInd $MI2 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.025 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO1='-doMajorMinor 1 -doMaf 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams.qc | wc -l`; export MI2=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc2
echo "source calc2 && angsd -b bams.qc -GL 1 $FILTERS1 $TODO1 -P 12 -out Bonnet1">g5
# Check technical replicates in hctree.pdf, check how many bams retained in bams.nr


# Final IBS production after removing replicates:
export MinIndPerc=0.75
FILTERS1='-minInd $MI2 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.025 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO1='-doMajorMinor 1 -doMaf 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams.nr | wc -l`; export MI2=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc3
echo "source calc3 && angsd -b bams.nr -GL 1 $FILTERS1 $TODO1 -P 12 -out Bonnet2">g5
# Download Bonnet2.ibsMat 






#####-------------- ADMIXTURE ---------###############
echo 'for K in `seq 4` ; do  NGSadmix -likes Bonnet2.beagle.gz -K $K -P 12 -o Bonnet2_${K}; done' >adm3
ls6_launcher_creator.py -j adm3 -n adm3 -a IBN21018 -e kblack@utexas.edu -t 0:30:00 -w 1 
sbatch adm3.slurm






#####-------------- RELATEDNESS ---------###############
# Rerun angsd one more time to produce a glf file for ngsRelate (must use -doGlf 3)
export MinIndPerc=0.75
FILTERS1='-minInd $MI2 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.025 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO1='-doMajorMinor 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 3'
echo 'export NIND=`cat bams.nr | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc3
echo "source calc3 && angsd -b bams.nr -GL 1 $FILTERS1 $TODO1 -P 12 -out Bonnet_rel">g4
# Relatedness with ngsRelate
echo 'export NIND2=`cat bams.qc | wc -l`; export NS=`zcat Bonnet_rel.mafs.gz | wc -l`' >calc
echo 'source calc && zcat Bonnet_rel.mafs.gz | cut -f5 |sed 1d >freq && ngsRelate  -g Bonnet_rel.glf.gz -n $NIND2 -f freq -O Bonnet.res >Bonnet.relatedness' >rel






#####----------- DEMOGRAPHIC MODELING ---------###############

# Rerun ANGSD on all well-genotyped sites for SFS production

# generating list of filtered SNP sites for SFS production (note: no filters that distort allele frequency!):
# sb_pval - strand bias filter; only use for 2bRAD, GBS or WGS (not for ddRAD or RADseq)
# hetbias - detects weird heterozygotes because they have unequal representation of alleles
# maxHetFreq - filters out lumped paralogs 

# ---- stage 1: finding well-genotytped sites (both variable and invariable) that are shared between two populations.
# we will run filtering command on each population, WITHOUT allele-frequency distorting filters (-minMaf and -snp_pval)
# make separate files listing bams for each population (without replicates)
grep '^TB*' bams.nr > TB_bams # 20
grep '^BB*' bams.nr > BB_bams # 69

# computing setting for the -minInd filter:
export GenRate=0.75 # desired genotyping rate
export N1=`wc -l TB_bams | cut -f 1 -d " "`
export N2=`wc -l BB_bams | cut -f 1 -d " "`
export MI1=`echo "($N1*$GenRate+0.5)/1" | bc`
export MI2=`echo "($N2*$GenRate+0.5)/1" | bc`

FILTERS='-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -maxHetFreq 0.5 -hetbias_pval 1e-3 -sb_pval 1e-3'
# add `-sb_pval 1e-3` (strand bias) to FILTERS if you have 2bRAD, GBS, or WGS data. Other types of RAD only sequence one strand so -sb_pval filter would remove everything.
export GENOME_REF="GCA_024679065.1_ASM2467906v1_genomic.fna" # reference to which the reads were mapped
TODO="-doHWE 1 -doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -dosnpstat 1 -doGeno 11 -doGlf 2 -anc $GENOME_REF -ref $GENOME_REF"
echo "angsd -b TB_bams -GL 1 -p 4 -minInd $MI1 $FILTERS $TODO -out Tampa
angsd -b BB_bams -GL 1 -P 4 -minInd $MI2 $FILTERS $TODO -out Biscayne" >afs

# collecting and indexing filter-passing sites
zcat Tampa.mafs.gz | cut -f 1,2 | tail -n +2 | sort | uniq >Tampa.sites
zcat Biscayne.mafs.gz | cut -f 1,2 | tail -n +2 | sort | uniq >Biscayne.sites
# collecting and indexing common sites:
comm -12 Tampa.sites Biscayne.sites | sort -k1 >TBSites # 4,971,503 sites
#get a list of chromosomes that have been genotyped in both
awk '{print $1}' TBSites | sort | uniq > chroms
#format this list for angsd by adding  colons
sed -i 's/$/:/' chroms
#index your new sites
angsd sites index TBSites

#now run angsd again to get your .saf files with this site and region file for each population - then can use them with gadma, realSFS, etc
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
echo "angsd -sites TBSites -rf chroms -b TB_bams -GL 1 -P 1 $TODO -out Tampa
angsd -sites TBSites -rf chroms -b BB_bams -GL 1 -P 1 $TODO -out Biscayne" > saf

# checking 2d sfs plot
# generating dadi-like posterior counts based on sfs priors
realSFS Tampa.saf.idx Biscayne.saf.idx >d2
SFSIZE="41 139" # 2N+1 for each population. 
echo $SFSIZE >d2.sfs;
cat d2 | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> d2.sfs;
# change numbers to 2 x 0.9 x number of samples for in each pop (can try different projections):
2dAFS.py d2.sfs TB BB 36 124 # 80%


# per site per generation mutation rates and generation times for great hammerhead, shortfin mako, white shark, brownbanded bamboo shark, cloudy catshark, and whale shark 
# species: mutation rate, generation length (Martin et al. 1999, Stanhope et al. 2023)
# great hammerhead: 3.92e-09, 24.5 years (estimated from mean synonymous site divergence of 50 random chosen single copy BUSCO genes and 55 MY separation from white shark)
# shortfin mako: 4.33e-09, 24.5 years
# white shark: 9.36e-09 (estimated as with mako), 53 years
# brownbanded bamboo shark: 1.02e-08, 9.5 years
# cloudy catshark: 7.43e-09, 9 years
# whale shark: 2.17e-08, 25 years
 
# pop0 is Tampa Bay, pop1 is Biscayne Bay
echo "#params_file
Input file: d2.sfs
Population labels: pop0, pop1
Projections: 36, 124 # change numbers to 2 x 0.9 x number of samples for in each pop 

Output directory: gadma_result_TampaBiscayne

Engine: moments

# If outgroup is False then SFS will be folded.
Outgroup: True

Theta0: Null
Relative parameters: False

Sequence length: 4971503
Mutation rate: 3.92e-09

Time for generation: 4.0


Initial structure: 1,1
Final structure: 2,1

Symmetric migrations: False
Split fractions: False

# We want output in stdout
Silence: False

# How many repeats to run and how many processes to use.
Number of repeats: 96
Number of processes: 48"  >params_file


# Make sure to remove gadma_result directory

echo "gadma -p params_file">gd
# run gd


# To look at the best model:
python3 best_logLL_model_moments_code.py






#####-------------- Check heterozygosity ------------###########
export GENOME_REF=GCA_024679065.1_ASM2467906v1_genomic.fna
FILTERS='-minInd 1 -setMinDepthInd 10 -uniqueOnly 1 -minMapQ 20 -minQ 20'
>hets
>mybams.het
for F in `cat bams.qc`; do
echo "angsd -i $F -anc $GENOME_REF $FILTERS -GL 1 -dosaf 1 -out ${F/.bam/} && realSFS ${F/.bam/}.saf.idx >${F/.bam/}.ml | awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ${F/.bam/}.ml >>mybams.het">>hets;
done
# execute all commands in hets; it will produce a tab-delimited table mybams.het:

>mybams.het
for F in `cat bams`; do
eval "awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ${F/.bam/}.ml" >>mybams.het
done
# [bam filename]   [total number of sites passing filters]   [heterozygosity]





#####-------------- Inbreeding ------------###########
export GENOME_REF=GCA_024679065.1_ASM2467906v1_genomic.fna
FILTERS='-minInd 1 -setMinDepthInd 10 -uniqueOnly 1 -minMapQ 20 -minQ 20'
>inb
for F in `cat bams`; do
echo "angsd -i $F -anc $GENOME_REF $FILTERS -GL 1 -HWE_pval_F 1 -doMajorMinor -doGlf 3 -out ${F/.bam/}">>inb;
done

# Rerunning angsd on all bams one more time to produce a glf file for ngsF (must use -doGlf 3)
FILTERS1='-minInd 1 -setMinDepthInd 10 -uniqueOnly 1 -minMapQ 20 -minQ 20'
TODO1='-doMajorMinor 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 3'
echo "angsd -b bams -GL 1 $FILTERS1 $TODO1 -P 12 -out inb">inb
ls6_launcher_creator.py -j inb -n inb -a IBN21018 -e kblack@utexas.edu -t 0:30:00 -w 1 
sbatch inb.slurm


# install ngsF (and pkg-config dependency through conda)
conda activate samtools
#conda install pkg-config
cdw
git clone https://github.com/fgvieira/ngsF.git
cd ngsF
nano Makefile
# add -I${TACC_GSL_INC}  to CC and CXX macros (CFLAGS= ...);
# and -L${TACC_GSL_LIB} to the 'LIB = ...' line.

module load gsl

export NIND=` cat bams.qc | wc -l `
export NS=$((`zcat inb.mafs.gz | wc -l`-1)) # 37192
echo "zcat inb.glf.gz | ngsF --n_ind $NIND --n_sites $NS --glf - --out inbreeding"> ngsf
#alternatively
gunzip inb.glf.gz 
echo "ngsF --n_ind $NIND --n_sites $NS --glf inb.glf --out inbreeding"> ngsf

# The output file is a text file with the per-individual inbreeding coefficients, one per line. 
# The parameters file is a binary file with the final parameters, namely global log-likelihood (double), 
# per-individual log-likelihood (N_IND doubles), per-individual inbreeding coefficients (N_IND doubles), 
# and per-site minor allele frequencies (N_SITES doubles)






#####-------------- Fst: global and per site ------------###########
# writing down 2d-SFS priors
realSFS Tampa.saf.idx Biscayne.saf.idx -P 24 > popTB.sfs ; realSFS fst index Tampa.saf.idx Biscayne.saf.idx -sfs popTB.sfs -fstout popTB

# global Fst between populations
realSFS fst stats popTB.fst.idx
# FST.Unweight[nObs:4971503]:0.002495 Fst.Weight:0.007394

# per-site Fst
realSFS fst print popTB.fst.idx > popTB.fst
# scaffold center-of-window fst-unweight, fst-weight
# sort file 
sort -r -k 3 popTB.fst | head -100 > popTB.fst.sort
