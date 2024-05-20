library(ggplot2)

#Make plot of peptide discovery across runs
good_peptide_discovery <- read.csv("/work/Bachelor_project/DA/Final_analysis/good_peptide_discovery.csv", 
                                   row.names=1)
ggplot(good_peptide_discovery) +
  geom_point(aes(x = Index, y = X1),
             color = "lightblue",
             alpha = 0.2,
             shape = 1) +
  geom_point(aes(x = Index, y = X2),
             color = "lightblue",
             alpha = 0.2,
             shape = 1) +
  geom_point(aes(x = Index, y = X3),
             color = "lightblue",
             alpha = 0.2,
             shape = 1) +
  geom_point(aes(x = Index, y = Mean),
             shape = 1) +
  labs(x = "Number of Runs", y = "Unique Peptides Discovered") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Peptide Discovery Across Runs")

#Makes the pepetides per run plot with a linear model
linear_model <- lm(good_peptide_discovery$Mean[1000:1419] ~ c(1000:1419))

ggplot(good_peptide_discovery) +
  geom_point(aes(x = Index, y = Mean),
             shape = 1) +
  geom_abline(slope = 36.76, intercept = (155487.03),
              color = "red") +
  labs(x = "Number of Runs", y = "Unique Peptides Discovered") +
  annotate("text",x = Inf, y = Inf, 
           label = paste("Slope", 36.76, sep = " = "),
           vjust = 4, hjust = 6.35) +
  annotate("text",x = Inf, y = Inf, 
           label = paste("R-Squared", round(summary(linear_model)$r.squared,2),
                         sep = " = "),
           vjust = 2, hjust = 4.93) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Peptide Discovery Across Runs, Linear Model")

ggplot(good_peptide_discovery[1000:1419,]) +
  geom_point(aes(Index, Mean), shape = 1) +
  geom_abline(slope = 36.76, 
              intercept = (155487.03), 
              color = "red",
              lwd = 1) +
  annotate("text",x = Inf, y = Inf, 
           label = paste("Slope", 36.76, sep = " = "),
           vjust = 4, hjust = 6.35) +
  annotate("text",x = Inf, y = Inf, 
           label = paste("R-Squared", round(summary(linear_model)$r.squared,2),
                         sep = " = "),
           vjust = 2, hjust = 4.93) +
  labs(x = "Number of Runs", y = "Unique Peptides Discovered") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Peptide Discovery Across Runs, Linear Model")

#Makes peptides per runs plot vs a logarithmic model
log_model <- lm(good_peptide_discovery$Mean[1:500] ~ log(good_peptide_discovery$Index[1:500]))
log_model
log_prediction <- c()
for (x in c(1:1419)){
  log_prediction <- append(log_prediction, log(x)*29023 - 23445)
}

good_peptide_discovery["log"] <- log_prediction

ggplot(data = good_peptide_discovery) +
  geom_point(aes(x= Index, y = Mean), color = "black",
             shape = 1) +
  geom_point(aes(x = Index, y = log), color = "red",
             shape = 1) +
  labs(x = "Number of Runs", y = "Unique Peptides Discovered") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, NA) +
  ggtitle("Peptide Discovery Across Runs, Logarithmic Model")

ggplot(good_peptide_discovery) +
  geom_point(aes(log(Index), Mean), shape = 1) +
  labs(x = "log(Number of runs)", y = "Unique Peptides Discovered") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Peptide Discovery Across Runs, Logarithmic Model") +
  geom_abline(slope = 29023, 
              intercept = (-23445), 
              color = "red",
              lwd = 1)

#Makes a plot of sum normalized peptide abundances with median line
peptides_sum_normalized <- read.csv("/work/Bachelor_project/DA/Final_analysis/sum_normalized_peptide_intensities.csv")
peptides_sum_normalized_median_list <- read.csv("/work/Bachelor_project/DA/Final_analysis/sum_normalized_peptide_intensities_median_list.csv")

ggplot(peptides_sum_normalized, aes(sample_coverage_n,
                               Summed.Intensity)) + 
  geom_point(alpha = 0.2) +
  geom_line(data = peptides_sum_normalized_median_list, 
            aes(sample_coverage,
                median_list_summed, 
                color = "Median"), 
            color = "red") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "Summed Intensity") +
  ggtitle("Sum Normalized Peptide Intensities")


#Can also be made with log on y-axis
ggplot(peptides_sum_normalized, aes(sample_coverage_n,
                                    log(Summed.Intensity))) + 
  geom_point(alpha = 0.2) +
  geom_line(data = peptides_sum_normalized_median_list, 
            aes(sample_coverage,
                log(median_list_summed), 
                color = "Median"), 
            color = "red") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "log(Summed Intensity)") +
  ggtitle("Summed Peptide Intensities, log Transformed")

#Make sum normalized as a heatmap
ggplot(peptides_sum_normalized, aes(sample_coverage_n,
                                    Summed.Intensity)) + 
  geom_bin2d(bins = 100) +
  scale_fill_continuous(name = "Count", low = "blue", high = "red") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "log(Summed Intensity)") +
  ggtitle("Heatmap of Summed Peptide Intensities")+
  scale_fill_continuous(trans = "log",
                        low = "blue", high = "red",
                        breaks = c(1, 5, 25, 125, 625, 3250),
                        name = "Count")

#Make log sum normalized as a heatmap
ggplot(peptides_sum_normalized, aes(sample_coverage_n,
                               log(Summed.Intensity))) + 
  geom_bin2d(bins = 100) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "log(Summed Intensity)") +
  ggtitle("Heatmap of Summed Peptide Intensities")+
  scale_fill_continuous(trans = "log",
                        type = "viridis",
                        breaks = c(1, 5, 25, 125, 625, 3250),
                        name = "Count")

#Make a plot of median normalized sum coverage peptides
peptides_median_normalized <- read.csv("/work/Bachelor_project/DA/Final_analysis/peptides_median_normalized.csv", 
                                       row.names=1)
peptides_median_normalized_median_list <- read.csv("/work/Bachelor_project/DA/Final_analysis/peptides_median_normalized_median_list.csv", 
                                                   row.names=1)
ggplot(peptides_median_normalized, aes(sample_coverage_n,
                      Median.Intensity)) + 
  geom_point(alpha = 0.2) +
  geom_line(data = peptides_median_normalized_median_list, 
            aes(Sample_coverage,
                short_median, 
                color = "Median"), 
            color = "red") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "Median Intensity") +
  ggtitle("Median Normalized Peptide Intensities")

#Make a heatmap of median normalized peptides
ggplot(peptides_median_normalized, 
       aes(sample_coverage_n, Median.Intensity)) +
  geom_bin2d(bins = 100) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "Median Intensity") +
  ggtitle("Heatmap of Median Normalized Peptide Intensities")+
  scale_fill_continuous(trans = "log",
                        type = "viridis",
                        breaks = c(1, 5, 25, 125, 625, 3250),
                        name = "Count")

#Make a plot for comparing Z-scores across sample coverage
z_plot_df <- read.csv("/work/Bachelor_project/DA/Final_analysis/new_z_scores.csv")

ggplot(z_plot_df) +
  geom_point(shape = 1, aes(unique_coverages, Z.score.summed,
                            color = "Summed Normalized")) +
  geom_point(shape = 1, aes(unique_coverages, Z.score.median,
                            color = "Median Normalized")) +
  labs(x = "Coverage", y = "Z-score") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        legend.position = "bottom") +
  ggtitle("Median Normalized vs Summed Normalized")

#Make peptides per coverage
peptides_per_coverage_df <- read.csv("/work/Bachelor_project/DA/Final_analysis/peptides_per_coverage_df.csv", 
                                     row.names=1)
ggplot(peptides_per_coverage_df, aes(Var1,Percent)) +
  geom_point(shape = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "Percentage of Total Peptides") +
  ggtitle("Peptides per Coverage")

ggplot(peptides_per_coverage_df, aes(Var1,log(Percent))) +
  geom_point(shape = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Sample Coverage", y = "log(Percentage of Total Peptides)") +
  ggtitle("Log Transformed Peptides per Coverage")

#Makes a plot of number of proteins with n peptides mapped to it
peptides_per_protein_per_run <- read.csv("/work/Bachelor_project/DA/Final_analysis/total_peptides_per_protein_per_run_10_good.csv", 
                                         row.names=1)
ggplot(peptides_per_protein_per_run) +
  geom_line(aes(Index, X2, colour = "1"), lwd = 1) + 
  geom_line(aes(Index, X3, colour = "2"), lwd = 1) +
  geom_line(aes(Index, X4, colour = "3"), lwd = 1) +
  geom_line(aes(Index, X5, colour = "4"), lwd = 1) +
  geom_line(aes(Index, X6, colour = "5"), lwd = 1) +
  geom_line(aes(Index, X7, colour = "6"), lwd = 1) +
  geom_line(aes(Index, X8, colour = "7"), lwd = 1) +
  geom_line(aes(Index, X9, colour = "8"), lwd = 1) +
  geom_line(aes(Index, X10, colour = "9"), lwd = 1) +
  geom_line(aes(Index, X11, colour = "10"), lwd = 1) +
  labs(x = "Number of Runs", y = "Number of Proteins") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  ggtitle("Total Peptides per Protein per Run") +
  guides(color = guide_legend(title = "Total peptides mapped to protein"))

#Make a plot of newly discovered proteins per run
peptides_per_protein_per_run <- read.csv("/work/Bachelor_project/DA/Final_analysis/total_peptides_per_protein_per_run_10_good.csv", 
                                         row.names=1)
ggplot(peptides_per_protein_per_run)+
  geom_point(aes(Index, 11888- X1), shape = 1) +
  labs(x = "Number of Runs", y = "Number of Unique Proteins") +
  scale_y_continuous(n.breaks = 6, labels = c(0, 2500, 5000, 7500, 10000, 12500)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Protein Discovery Across Runs")

#Make plot of peptides per protein as function of sample coverage of the protein
protein_sample_coverage_total_peptides <- read.csv("/work/Bachelor_project/DA/Final_analysis/protein_sample_coverage_total_peptides.csv", 
                                                   row.names=1)
ggplot(protein_sample_coverage_total_peptides) +
  geom_point(aes(x = protein_sample_coverage,
                 y = total_peptides_per_protein)) + 
  labs(x = "Protein Sample Coverage", 
       y = "Peptides per Protein") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Peptides per Protein by Sample Coverage")

ggplot(protein_sample_coverage_total_peptides) +
  geom_point(aes(x = protein_sample_coverage,
                 y = log(total_peptides_per_protein))) + 
  labs(x = "Protein Sample Coverage", 
       y = "log(Peptides per Protein)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Peptides per Protein by Sample Coverage")

#Make protein discovery plot with fitted linear model from run 1000 and onward
protein_discovery <- read.csv("/work/Bachelor_project/DA/Final_analysis/protein_discovery.csv")
protein_discovery_lm <- lm(protein_discovery$Mean[1000:1419] ~ protein_discovery$Index[1000:1419])

protein_discovery_lm
ggplot(protein_discovery) +
  geom_point(aes(x= Index, y = Mean)) +
  labs(x = "Number of Runs", y = "Unique Proteins Discovered") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Unique Proteins Discovered Across Runs") +
  geom_abline(slope = 1.234, intercept = 9953.063,
              color = "red") +
  annotate("text",x = Inf, y = Inf, 
           label = paste("Slope", 1.23, sep = " = "),
           vjust = 4, hjust = 7) +
  annotate("text",x = Inf, y = Inf, 
           label = paste("R-Squared", round(summary(protein_discovery_lm)$r.squared,4),
                         sep = " = "),
           vjust = 2, hjust = 4.35)

#Different PTM percentage plot as function of total PSM per run
PTMs_and_PSMs <- read.csv("/work/Bachelor_project/DA/Final_analysis/PTMs_and_PSMs.csv")
PTMs_and_PSMs[PTMs_and_PSMs == 0] <- NA

ggplot(PTMs_and_PSMs) + 
  geom_point(aes(total_psm, log(Methylations/total_psm), 
                 color = "Methylation"), alpha = 0.6) +
  geom_point(aes(total_psm, log(Oxidation/total_psm), 
                 color = "Oxidation"), alpha = 0.6) +
  geom_point(aes(total_psm, log(Acetyl/total_psm), 
                 color = "Acetyl"), alpha = 0.6) +
  geom_point(aes(total_psm, log(Carbamyl/total_psm), 
                 color = "Carbamyl"), alpha = 0.6) +
  geom_point(aes(total_psm, log(Phosphorylation/total_psm), 
                 color = "Phosphorylation"), alpha = 0.6) +
  geom_point(aes(total_psm, log(Carboxylation/total_psm), 
                 color = "Carboxylation"), alpha = 0.6) +
  geom_point(aes(total_psm, log(Formylation/total_psm), 
                 color = "Formylation"), alpha = 0.6) +
  labs(x = "Total PSM", y = "log(PTM percentage)") + 
  scale_color_manual(name='PTM type',
                     breaks=c('Methylation', 'Oxidation', 'Acetyl', 'Carbamyl', 
                              'Phosphorylation', 'Carboxylation', 'Formylation'),
                     values=c('Methylation'='blue','Oxidation' = 'green', 
                              'Acetyl'='red', 'Carbamyl' = 'orange',
                              'Phosphorylation'='yellow', 'Carboxylation' = 'purple',
                              'Formylation' = 'darkgreen')) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  ggtitle("PTM Percentage per PSM")
