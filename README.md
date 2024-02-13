# CRISPRFACS_simulation
R code generating simulated CRISPR pool screen coupled with FACS sorting using different gate size.
Verious parameters can be changed in order to evaluate case scenarious.
--> PAramenters are assembled into a grid to generate as many test cases with variable gRNA counts per targets, replicates, sorted bins etc..
--> read counts per gRNA are modelled assa negative binomial distribution giving positive starting probability to those elements withih  the theoretical bin and negative to those excluded. This strategy ensures what is actually observed in CRISPR pool screen with FACS i.e. every gRNA is detected in every bin even opposite. The sorting in bins is solely contributing to generate an enrichments as reported elsewhere: see "A Genome-wide CRISPR Screen in Primary Immune Cells to Dissect Regulatory Networks" by Oren Parnas et al. 2015
