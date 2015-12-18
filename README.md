
quantitative-gene-expression

process data come from quantitative real-time PCR

1. the purpose is to display relative mRNA level of the target genes
2. computation based on delta-delta-Ct method
3. important parameters:
(1) gene used for internal control, for example, 'Cyclophinlin B'
(2) mapping the individual samples to the corresponding groups

currently, the program is based on two treatments with interactions, here they are called 'conditions' and 'treatments'.

So, each sample belong to one group within the condition groups, and also belong to one group within the treatment groups.

results are exported to several .CSV files and three .PNG figures

raw data file - this is the file exported from the machine, may need a little work before puting into the program.

relative expression levels for each sample with standard deviation from the triplicates

averaged expression levels for each group with standard deviation

p-values from pairwise student t-tests

p-values from two-way ANOVA analysis

figure 1 display bar chart from each sample, grouped by condition and treatment groups

figure 2 dispaly box-plot for each group, to help figure out the outliers

figure 3 display the bar chart of averge relative expression lever for each group

outlieres handling:
triplicates PCR for each sample and for each gene, including the internal control. Outliers from triplicates are removed based on big_d=max(max-median, median-min), small_d=min(max-median, median-min) big_d > 2.5 fold of small_d; 2) big_d is larger than 1/50 of the median value of the triplicates

outliers from each group are removed based on absolute deviation > 2 * std
