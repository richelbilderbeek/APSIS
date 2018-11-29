# APSIS: Analysis of Power for Sequencing-and-Imputation Studies

APSIS calculates statistical power for genome-wide association studies (GWAS) using augmented reference panels.  Users can calculate power and total genotyping cost as a function of the number of GWAS participants sequenced (and included in the imputation reference panel) and imputed based on data from two admixed populations (African Americans and Latino Americans) and two European isolate populations (Sardinians and Finns), and three Illumina genotyping arrays.

#### Getting Started 
```
Load the source code and data (~28MB, which may take a moment):

## note: apsis requires the 'data.table' library
source('apsis.r')
```

#### Calculating and Visualizing Power

Use diseaseModel() to specify effect size, minor allele frequency (MAF), and other parameters

```
my_model <- diseaseModel(maf=0.002, rr=2.5, prev=0.1)
```

Effect size can be specified as 'rr' for relative risk, or 'beta' for the liability threshold model. Other required parameters are population MAF ('maf') and disease prevalence ('prev'). Additionally, you can specify 'type' = "recessive" or "dominant" for non-additive/multiplicative models. 

Calculate power for the specified disease model, varying the numbers of participants sequenced and imputed:

```
my_pcal_1 <- apsisPower(my_model, 
	n_imputed = seq(5000,80000,1000), 
	n_sequenced = c(0,250,1000),
	population = c("Sardinians", "African Americans"),
	snp_array = c("Infinium Core", "OmniExpress")
)
```

Visualize power calculations across genotyping arrays and populations:

```
library(ggplot2)
ggplot(my_pcal_1, aes(y=Power, x=N_Imputed, linetype = Population, colour=factor(N_Sequenced))) + 
	geom_line() + scale_x_log10() + facet_grid(~Array)
```

![alt text](https://raw.githubusercontent.com/corbinq/APSIS/example_plots/plot1.png)


#### Comparing Across Disease Models

Use diseaseModelList to compare results across disease model parameters:

```
my_model_list <- diseaseModelList(maf = (10^seq(-3,-2, 0.02)), rr= c(1.5,1.7,2.5), prev = 0.01)

my_pcal_2 <- apsisPower(my_model_list, 
	n_imputed = 80000, 
	n_sequenced = c(0,1000),
	population = "African Americans",
	snp_array = "Infinium Core"
)
```

Visualize power as a function of disease model parameters:

```
library(ggplot2)
ggplot(my_pcal_2, aes(y=Power, x=MAF, colour=factor(EffSize), linetype=factor(N_Sequenced))) + 
	geom_line() + scale_x_log10() 
```

![alt text](https://raw.githubusercontent.com/corbinq/APSIS/example_plots/plot2.png)
