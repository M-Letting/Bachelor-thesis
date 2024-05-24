#Creates a data frame with n rows and m columns with all elements as NA
create_big_df <- function(n, m){
  df <- data.frame(matrix(NA, 
                          nrow = n,
                          ncol = m))
  return(df)
}

#Creates a data frame with the number of peptides discovered across runs in a df
#Repeats this n times 
peptide_discovery <- function(df, n){
  result_df <- create_big_df(length(df), n)
  for (x in seq(1:n)){
    shuffled_df <- df[,sample(ncol(df))]
    result_df[1, x] <- sum(!is.na(shuffled_df[1]))
    for (i in c(2:ncol(shuffled_df))) {
      y <- rowSums(shuffled_df[, 1:i], na.rm = TRUE)
      y[y == 0] <- NA
      result_df[i, x] <- sum(!is.na(y))
      if (i %% 100 == 0){
        print(paste(x,i, sep = "."))
      }
    }
  }
  result_df["Index"] <- c(1:nrow(result_df))
  return(result_df)
}

#Normalizes a data frame by taking log of all values and subtracting the
#median of each column from all values of that column
median_normalize_dataframe <- function(df){
  df <- apply(df, 2, function(x) log(x))
  df <- as.data.frame(df)
  df <- apply(df, 2, function(x) x - median(x, na.rm = TRUE))
  result <- as.data.frame(df)
  return(result)
}

#Creates a the a dataframe for making summed intensities plot
sum_normalized <- function(df){
  require(dplyr)
  sample_coverage_n <- apply(df, 1,function(x) sum(!is.na(x)))
  sample_coverage_df <- as.data.frame(sample_coverage_n)
  df_normalized <- mutate_all(df, ~./median(., na.rm = TRUE))
  summed_intensity <- rowSums(df_normalized, na.rm =TRUE)
  sample_coverage_df["Summed.Intensity"] <- summed_intensity
  result_df <- sample_coverage_df[order(-sample_coverage_df$sample_coverage_n),]
  return(result_df)
}


#Make median list for summed intensities divided by median intensity
library("dplyr")
make_median_summed_list <- function(df, n){
  median_list_summed <- c()
  for (x in rev(unique(df$sample_coverage_n))){
    coverage_filtered <- filter(df,sample_coverage_n == x)
    median_list_summed <- append(median_list_summed,
                                 median(coverage_filtered$Summed.Intensity,
                                        na.rm = TRUE))
  }
  medians_summed <- data.frame(rev(unique(df$sample_coverage_n)), 
                               median_list_summed)
  medians_summed_short <- medians_summed[1:n,]
  colnames(medians_summed_short)[1] <- "sample_coverage"
  return(medians_summed_short)
}

#Creates a dataframe for making median normalized plot
median_normalized <- function(df){
  sample_coverage_n <- apply(df, 1,function(x) sum(!is.na(x)))
  median_df <- as.data.frame(sample_coverage_n)
  median_intensity <- apply(df, 
                            1, function(x) median(x, na.rm = TRUE))
  median_df["Median.Intensity"] <- median_intensity
  median_df <- median_df[order(median_df$sample_coverage_n,
                               decreasing = TRUE),]
  return(median_df)
}

#Creates a list with the medians of the first n unique sample coverages
#Works on a data frame created by sum_normalized()
library("dplyr")
make_median_sum_list <- function(df, n){
  median_list_summed <- c()
  for (x in rev(unique(df$sample_coverage_n))){
    coverage_filtered <- filter(df,sample_coverage_n == x)
    median_list_summed <- append(median_list_summed,
                                 median(coverage_filtered$Summed.Intensity,
                                        na.rm = TRUE))
  }
  medians_summed <- data.frame(rev(unique(df$sample_coverage_n)), 
                               median_list_summed)
  medians_summed_short <- medians_summed[1:n,]
  colnames(medians_summed_short)[1] <- "sample_coverage"
  return(medians_summed_short)
}

#Creates a list with the medians of the first n unique sample coverages
#Works on a data frame created by median_normalized()
library("dplyr")
make_median_median_list <- function(df, n){
  median_list_medians <- c()
  for (x in rev(unique(df$sample_coverage_n))){
    coverage_filtered <- filter(df, sample_coverage_n == x)
    median_list_medians <- append(median_list_medians,
                                  median(coverage_filtered$Median.Intensity,
                                         na.rm = TRUE))
  }
  short_median <- median_list_medians[1:n]
  result_df <- as.data.frame(short_median)
  result_df["Sample_coverage"] <- c(1:n)
  return(result_df)
}

#Make z-score plot
#Coverts a data frame with median or sum normalized data to z-scores
make_z_scores <- function(median_df, summed_df){
  z_score_df <- median_df
  z_score_df["Summed.Intensity"] <- summed_df[,2]
  unique_coverages <- unique(z_score_df$sample_coverage_n)
  z_plot_df <- as.data.frame(unique_coverages)
  
  median_of_unique_coverages_medians <- c()
  for (x in unique_coverages){
    to_be_added <- median(z_score_df[grep(x, z_score_df$sample_coverage_n), 2])
    median_of_unique_coverages_medians <- append(median_of_unique_coverages_medians, 
                                                 to_be_added)
  }
  
  median_of_unique_coverages_summed <- c()
  for (x in unique_coverages){
    to_be_added <- median(z_score_df[grep(x, z_score_df$sample_coverage_n), 3])
    median_of_unique_coverages_summed <- append(median_of_unique_coverages_summed, 
                                                to_be_added)
  }
  z_plot_df["Median.of.summed.unique.coverages"] <- median_of_unique_coverages_summed
  z_plot_df["Median.of.median.unique.coverages"] <- median_of_unique_coverages_medians
  
  mean_summed_unique <- mean(z_plot_df$Median.of.summed.unique.coverages)
  sd_summed_unique <- sd(z_plot_df$Median.of.summed.unique.coverages)
  z_score_summed <- (z_plot_df$Median.of.summed.unique.coverages - mean_summed_unique)/sd_summed_unique
  z_plot_df["Z.score.summed"] <- z_score_summed
  
  mean_median_unique <- mean(z_plot_df$Median.of.median.unique.coverages, na.rm = TRUE)
  sd_median_unique <- sd(z_plot_df$Median.of.median.unique.coverages, na.rm = TRUE)
  z_score_median <- (z_plot_df$Median.of.median.unique.coverages - mean_median_unique)/sd_median_unique
  z_plot_df["Z.score.median"] <- z_score_median
  
  return(z_plot_df)
}

#Petides per coverage
#Makes the data for creating peptides per coverage plot
make_peptides_per_coverage <- function(df){
  sample_coverage_n <- apply(df, 1,function(x) sum(!is.na(x)))
  sample_coverage_df <- as.data.frame(sample_coverage_n)
  unique_coverages <- unique(sample_coverage_df$sample_coverage_n)
  peptides_per_coverage <- as.data.frame(table(sample_coverage_df$sample_coverage_n))
  percent_per_coverage <- peptides_per_coverage$Freq / sum(peptides_per_coverage$Freq)
  peptides_per_coverage["Percent"] <- percent_per_coverage
  peptides_per_coverage["Var2"] <- as.double(peptides_per_coverage$Var1)
  return(peptides_per_coverage)
}

#Finds the number of proteins with up to n total peptides mapped to them
#repeats this 3 times before returning the average
make_total_peptides_per_protein_per_run <- function(df, n) {
  res1 <- create_big_df(ncol(df), n)
  sample1 <- df[,sample(1:ncol(df))]
  for (i in c(2:ncol(sample1))) {
    sum_of_rows <- rowSums(sample1[, 1:i], na.rm = TRUE)
    counts <- table(sum_of_rows)
    res1[i, ] <- counts[seq_len(n)]
  }
  print("Iteration 1 done")
  res2 <- create_big_df(ncol(df), n)
  sample2 <- df[,sample(1:ncol(df))]
  for (i in c(2:ncol(sample2))) {
    sum_of_rows <- rowSums(sample2[, 1:i], na.rm = TRUE)
    counts <- table(sum_of_rows)
    res2[i, ] <- counts[seq_len(n)]
  }
  print("Iteration 2 done")
  res3 <- create_big_df(ncol(df), n)
  sample2 <- df[,sample(1:ncol(df))]
  for (i in c(2:ncol(df))) {
    sum_of_rows <- rowSums(df[, 1:i], na.rm = TRUE)
    counts <- table(sum_of_rows)
    res3[i, ] <- counts[seq_len(n)]
  }
  print("Iteration 3 done")
  result_df <- create_big_df(ncol(df), n)
  for (x in seq_along(result_df)){
    result_df[x] <- rowMeans(as.data.frame(c(res1[x], res2[x],res3[x]),
                                           row.names = c(1:ncol(df))), 
                             na.rm = TRUE)
  }
  result_df["Index"] <- c(1:ncol(df))
  return(result_df)
}

#makes a correlation plot of n random samples from df, with the title string
make_cor_plot <- function(df, n, min_cov = 1, max_cov = nrow(df), string){
  to_sample <- df[df[,1] >= min_cov & df[,1] <= max_cov, ]
  sample1 <- sample(to_sample[,2], n)
  sample2 <- sample(to_sample[,2], n)
  samples <- data.frame(sort(sample1), sort(sample2))
  pearson <- round(cor(sort(sample1), sort(sample2)),4)
  r_squared <- round(summary(lm(sort(sample2)~sort(sample1)))$r.squared, 4)
  require(ggplot2)
  ggplot(samples) +
    geom_point(aes(sort.sample1., sort.sample2.),
               shape = 1) +
    geom_abline(slope = 1, color = "red") +
    labs(x = paste(n, "random samples"), y = paste(n, "random samples"))+
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text",x = Inf, y = Inf, 
             label = paste("Pearson", pearson, sep = "="),
             vjust = 2, hjust = 5.3) +
    annotate("text",x = Inf, y = Inf, 
             label = paste("R-squared", r_squared, sep = "="),
             vjust = 4, hjust = 5) +
    ggtitle(string)
}

#Makes a Venn diagram to compare the peptide discovery across machines
machine_comparrison <- function(df){
  exp_df <- df[, grep("EXP", names(df), value = TRUE)]
  exp_df <- exp_df[rowSums(is.na(exp_df)) != ncol(exp_df), ]
  exp_peptides <- rownames(exp_df)
  
  lum_df <- df[, grep("lum", names(df), value = TRUE,
                      ignore.case = TRUE)]
  lum_df <- lum_df[rowSums(is.na(lum_df)) != ncol(lum_df),]
  lum_peptides <- rownames(lum_df)
  
  qehf_df <- df[, grep("qehf", names(df), value = TRUE,
                       ignore.case = TRUE)]
  qehf_df <- qehf_df[rowSums(is.na(qehf_df)) != ncol(qehf_df), ]
  qehf_peptides <- rownames(qehf_df)
  
  eclip_df <- df[, grep("eclip", names(df), value = TRUE,
                        ignore.case = TRUE)]
  eclip_df <- eclip_df[rowSums(is.na(eclip_df)) != ncol(eclip_df), ]
  eclip_peptides <- rownames(eclip_df)
  
  require(tidyverse)
  require(VennDiagram)
  venn_plot <- venn.diagram(x = list(Exploris = exp_peptides, 
                                     Lumos = lum_peptides,
                                     QExactiveHF = qehf_peptides, 
                                     Eclipse = eclip_peptides),
                            filename = NULL,
                            category.names = c("Exploris", "Lumos",
                                               "Q-exactiveHF", "Ecplise"),
                            output = TRUE,
                            fill = c(alpha("#FF0000", 0.5), alpha("#00FF00", 0.5),
                                     alpha("#0000FF", 0.5),alpha("#FFFF00", 0.5)),
                            col = c("Red", "Green","Blue","Yellow"))
  grid.draw(venn_plot)
}
