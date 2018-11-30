# APSIS: Analysis of Power for Sequencing-and-Imputation Studies

APSIS calculates statistical power for genome-wide association studies (GWAS), accounting for imperfect genotype imputation coverage and accuracy. In particular, APSIS can calculate power for studies in which a subset of GWAS participants are whole-genome sequenced, and remaining participants are array-genotyped and imputed using an augmented imputation reference panel. 

Users can compare power and total genotyping cost as functions of the number of GWAS participants sequenced (and included in the imputation reference panel) and imputed based on data from two admixed populations (African Americans and Latino Americans) and two European isolate populations (Sardinians and Finns) using three Illumina genotyping arrays (Infinium Core, OmniExpress, and Omni2.5).

#### Getting Started 

Load APSIS source code and data:

```
## If 'data.table' library is not already installed:
intall.packages("data.table")

## Note: this step involves reading ~28MB of data, which may take a moment to load:
source('apsis.r')

```

#### Calculating and Visualizing Power

Use diseaseModel() to specify effect size, minor allele frequency (MAF), and other parameters

```
my_model <- diseaseModel(
	maf=0.002, # minor allele frequency
	prev=0.1,  # disease prevalence
	rr=2.5     # use 'rr' for relative risk, or 'beta' for liability threshold model
)                  # for non-additive models, specify 'type' = "recessive" or "dominant"

```

Calculate power for the specified disease model, varying the numbers of participants sequenced and imputed:

```
## apsisPower() returns a data.table containing power calculations and parameters

my_pcal_dt_1 <- apsisPower(my_model, 
	n_imputed = seq(5000,80000,1000), 
	n_sequenced = c(0,250,1000),
	population = c("Sardinians", "African Americans"),
	snp_array = c("Infinium Core", "OmniExpress")
)

```

Visualize power calculations across genotyping arrays and populations:

```
library(ggplot2)
ggplot(my_pcal_dt_1, aes(y=Power, x=N_Imputed, linetype = Population, colour=factor(N_Sequenced))) + 
	geom_line() + scale_x_log10() + facet_grid(~Array)
```

![alt text](https://github.com/corbinq/APSIS/blob/master/example_plots/plot1.png)


#### Comparing Across Disease Models

Use diseaseModelList() to specify multiple disease models and compare results across parameter values:

```
## diseaseModelList() generates all MAF x RR x ... combinations across the specified values:

my_model_list <- diseaseModelList(
	maf = 10^seq(-3,-2, 0.02), 
	rr= c(1.5,1.7,2.5), 
	prev = 0.01
)

my_pcal_2 <- apsisPower(my_model_list, 
	n_imputed = 80000, 
	n_sequenced = c(0,1000),
	population = "African Americans",
	snp_array = "Infinium Core"
)
```

Visualize power as a function of MAF and effect size parameters:

```
library(ggplot2)
ggplot(my_pcal_2, aes(y=Power, x=MAF, colour=factor(EffSize), linetype=factor(N_Sequenced))) + 
	geom_line() + scale_x_log10() 
```

![alt text](https://github.com/corbinq/APSIS/blob/master/example_plots/plot2.png)
