# Bonnethead_shark_demographic_modeling
Scripts and data to accompany the manuscript "Evidence of unidirectional gene flow in bonnethead sharks from the Gulf coast to the Atlantic coast of Florida"

Bonnethead_2bRAD.sh contains instructions to trim and filter reads, genotype samples, run model selection for demographic inference, and calculate Fst and inbreeding coefficients. Data input for these scripts include: d2.sfs (site frequency spectrum for GADMA) and "sfs" (folder): contains 100 bootstrapped sfs for custom model selection procedure. All genetic data can be accessed from the Sequence Read Archive under Bioproject: PRJNA1082259 (Accession nos.: SAMN40202348- SAMN40202454).

Use Bonnethead_popstructure.R to reproduce figures from the manuscript- including dendrograms and PCoA of IBS and relatedness genetic distances, Admixture barplots, and RDA biplots. Data input for these scripts include: bams.qc (list of bam files or individual samples), Bonnet2.ibsMat (Identity-by-State genetic distance matrix), four *qopt files (for admixture scenarios k=1-4), and Bonnets_meta.txt (containing site and morphological information for each sample).
