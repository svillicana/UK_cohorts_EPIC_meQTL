library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

cis_meqtl_enrichment_bs <- fread(paste0(path,"enrichment/enrichment_top_cis_meQTLs_resamp_1000.txt"),
                                 header = TRUE, data.table = FALSE)
trans_meqtl_enrichment_bs <- fread(paste0(path,"enrichment/enrichment_top_trans_meQTLs_resamp_1000.txt"),
                                   header = TRUE, data.table = FALSE)

# Order factors  --------------------------------------------------------------

cat <- cis_meqtl_enrichment_bs %>%
  select(description, collection) %>%
  distinct() %>%
  arrange(collection,description)

cis_meqtl_enrichment_bs <- cis_meqtl_enrichment_bs %>%
  mutate(description = factor(description, levels = cat$description))

trans_meqtl_enrichment_bs <- trans_meqtl_enrichment_bs %>%
  mutate(description = factor(description, levels = cat$description))

# Plot distribution of OR -----------------------------------------------------

plot_dist_or_cis_meqtl <- cis_meqtl_enrichment_bs %>%
  ggplot(aes(x = oddsRatio, y = description, fill = collection)) +
  geom_density_ridges() +
  scale_x_continuous(trans = "log2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.y = element_blank())

plot_dist_or_trans_meqtl <- trans_meqtl_enrichment_bs %>%
  ggplot(aes(x = oddsRatio, y = description, fill = collection)) +
  geom_density_ridges() +
  scale_x_continuous(trans = "log2") +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.title.y = element_blank())

# Summary results -------------------------------------------------------------
  
cis_meqtl_enrichment_summary <- cis_meqtl_enrichment_bs %>%
  mutate(userSet = "cis-meQTLs") %>%
  group_by(across(c(userSet:collection,description:size))) %>%
  summarise(oddsRatio_mean = mean(oddsRatio),
            oddsRatio_lower = quantile(oddsRatio, 0.025),
            oddsRatio_upper = quantile(oddsRatio, 0.975)) %>%
  ungroup()

trans_meqtl_enrichment_summary <- trans_meqtl_enrichment_bs %>%
  mutate(userSet = "trans-meQTLs") %>%
  group_by(across(c(userSet:collection,description:size))) %>%
  summarise(oddsRatio_mean = mean(oddsRatio),
            oddsRatio_lower = quantile(oddsRatio, 0.025),
            oddsRatio_upper = quantile(oddsRatio, 0.975)) %>%
  ungroup()

meqtl_enrichment_summary <- bind_rows(cis_meqtl_enrichment_summary,
                                      trans_meqtl_enrichment_summary) %>%
  rename(oddsRatio = oddsRatio_mean) %>%
  arrange(desc(oddsRatio))

# Save -------------------------------------------------------------

# Plots
saveRDS(plot_dist_or_cis_meqtl, file = paste0(path,"enrichment/enrichment_top_cis_meQTLs_resamp_1000_dist.rds"))
saveRDS(plot_dist_or_trans_meqtl, file = paste0(path,"enrichment/enrichment_top_trans_meQTLs_resamp_1000_dist.rds"))

# Summary table
write.table(meqtl_enrichment_summary, file = paste0(path,"enrichment/enrichment_top_meQTLs_resamp_1000_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
