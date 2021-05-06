source("apsis.r")

my_model <- diseaseModel(
	maf = 0.002, # minor allele frequency
	prev = 0.1,  # disease prevalence
	rr = 2.5     # use 'rr' for relative risk, or 'beta' for liability threshold model
)                  # for non-additive models, specify 'type' = "recessive" or "dominant"

my_model

my_pcal_dt_1 <- apsisPower(my_model,
	n_imputed = seq(5000, 80000, 1000),
	n_sequenced = c(0, 250, 1000),
	population = c("Sardinians", "African Americans", "Finns", "Latino Americans"),
	snp_array = c("Infinium Core", "OmniExpress", "Infinium Omni2.5")
)

my_pcal_dt_1

t <- dplyr::filter(my_pcal_dt_1, N_Imputed == 5000 & N_Sequenced == 1000)
t$power_per_dollar <- t$Power / t$Cost
t
t$Array <- factor(t$Array, levels = c("Infinium Core", "OmniExpress", "Infinium Omni2.5"))

# From the article:
#Illumina Infinium Core 	307K
#Illumina Infinium OmniExpress 	710K
#Illumina Infinium Omni2.5 	2.5M
t$n_snps <- NA
t$n_snps[t$Array == "Infinium Core"] <- 307000
t$n_snps[t$Array == "OmniExpress"] <- 710000
t$n_snps[t$Array == "Infinium Omni2.5"] <- 2500000


ggplot2::ggplot(t, ggplot2::aes(x = Array, y = power_per_dollar, fill = Population)) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(. ~ Population) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = ggplot2::element_text(size = 16)
  ) + ggplot2::ggsave("power_per_dollar_per_array.png", width = 7, height = 7)



ggplot2::ggplot(t, ggplot2::aes(x = n_snps, y = power_per_dollar, color = Population)) +
  ggplot2::geom_point(size = 4) +
  ggplot2::geom_line(size = 2) +
  ggplot2::scale_x_continuous(labels = scales::scientific) +
  ggplot2::facet_grid(. ~ Population) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = ggplot2::element_text(size = 16)
  ) + ggplot2::ggsave("power_per_dollar_per_n_snps.png", width = 7, height = 7)

