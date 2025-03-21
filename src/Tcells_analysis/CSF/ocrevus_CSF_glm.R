library(tidyverse)
library(ggbeeswarm)
library(betareg)
library(broom)
library(forestplot)
library(forcats)
library(RColorBrewer)
library(patchwork)
library(ggsignif)


setwd('~/autoimmune_10x/')

project.name <- "ocrevus_CSF_Yale"

df.metadata <- read_csv(paste0('data/', project.name,'/', project.name,'_metadata.csv'))
# df.metadata$condition <- df.metadata$condition %>% recode("0M" = 'baseline',
#                                  "6M" = 'post',
#                                  "12M" = 'post',
#                                  '18M' = 'post',
#                                  "CSF" = "Healthy")
df.output <- read_csv(paste0('output/', project.name, '/',  project.name, '_queryL2_Reference_Mapping.csv'))

df.output <- df.output %>% 
  left_join(df.metadata, by = "sample") %>%
  group_by(sample) %>%
  mutate(freq = n / sum(n))

df.output$condition <- recode(df.output$condition,
                              baseline = "Baseline",
                              healthy = "Healthy",
                              post = "Follow_up")
df.output$condition <- factor(df.output$condition, levels = c("Baseline", "Follow_up", "Healthy"))

df.ntotal <- df.output %>% group_by(sample) %>% summarise(n_total = sum(n))
cols <- c("donor", "sample", "disease", "disease_duration", "condition", "age", "race", 
          "sex", "immunosuppressant", "immunosuppressant_duration")
df.metadata <- df.output[,cols] %>% distinct()


df.output.wide <- df.output %>% 
  pivot_wider(values_from = freq, names_from = clusterL2,
              id_cols = c(donor, sample, disease, disease_duration, condition, age, race, sex,
                          immunosuppressant, immunosuppressant_duration), values_fill = 0) %>%
  left_join(df.ntotal, by='sample')

cells <- c("Tnaive", "TnaiveAct", "TnaiveMX1", "TnaiveSOX4","TcmTh0","TcmTh0Act", 
           "TcmTfh", "TcmPHLDA3", "TcmTh17", "TcmTh2",
           "TemTh1pre", "TemTh1", "TemTh117", "TemTph", "TemraTh1",
           "TregNaive", "TregAct", "TregEff")

## GLM for baseline vs post
df_MS <- df.output.wide %>% 
  filter(condition != "Healthy")
df_MS$condition <- factor(df_MS$condition, levels = c("Baseline", "Follow_up"))

list.res <- list()

for (cell in cells) {
  df_MS <- df_MS %>% 
    mutate(!!cell := !!sym(cell) + 0.001)
  formula_str <- as.formula(paste0(cell, " ~ condition + donor"))
  res <- betareg(formula_str, data = df_MS)
  tidy_res <- broom::tidy(res)
  
  
  converged <- res$converged
  
  if (!converged) {
    tidy_res$p.value <- NA
  }
  
  d.res <- tidy_res %>%
    filter(term == "conditionFollow_up") %>%
    mutate(clusterL2 = cell)
  list.res[[cell]] <- d.res
}

df.res <- bind_rows(list.res) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

df.forest <- df.res %>%
  select(clusterL2, estimate, `std.error`) %>%
  mutate(
    lower = estimate - 1.96 * `std.error`,
    upper = estimate + 1.96 * `std.error`
  )

# Label preparations
labeltext <- cbind(
  df.forest$cell,
  formatC(df.forest$Estimate, digits=2),
  formatC(df.forest$lower, digits=2),
  formatC(df.forest$upper, digits=2)
)

# # Plot
# pdf(file=paste0('output/', project.name, '/',  project.name, "_cellfreq_forest_clusterL2.base_post.glm.betareg.pdf"), width = 5, height = 5)
# forestplot(
#   # labeltext = labeltext,
#   labeltext = cells,
#   # graph.pos = 3,
#   mean = df.forest$estimate,
#   lower = df.forest$lower,
#   upper = df.forest$upper,
#   xlab = "Coefficient",
#   zero = 0,
#   line.margin = .15,
#   col=fpColors(box="royalblue", lines="darkblue", zero = "gray50")
# )
# dev.off()

## GLM for baseline vs healthy
df_MS <- df.output.wide %>% 
  filter(condition != "Follow_up")
df_MS$condition <- factor(df_MS$condition, levels = c("Healthy", "Baseline"))

list.res <- list()

for (cell in cells) {
  df_MS <- df_MS %>% 
    mutate(!!cell := !!sym(cell) + 0.001)
  formula_str <- as.formula(paste0(cell, " ~ condition"))
  res <- betareg(formula_str, data = df_MS)
  tidy_res <- broom::tidy(res)
  
  
  converged <- res$converged
  
  if (!converged) {
    tidy_res$p.value <- NA
  }
  
  d.res <- tidy_res %>%
    filter(term == "conditionBaseline") %>%
    mutate(clusterL2 = cell)
  list.res[[cell]] <- d.res
}

df.res.disease <- bind_rows(list.res) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

df.forest <- df.res.disease %>%
  select(clusterL2, estimate, `std.error`) %>%
  mutate(
    lower = estimate - 1.96 * `std.error`,
    upper = estimate + 1.96 * `std.error`
  )

# Label preparations
labeltext <- cbind(
  df.forest$cell,
  formatC(df.forest$Estimate, digits=2),
  formatC(df.forest$lower, digits=2),
  formatC(df.forest$upper, digits=2)
)

# Plot
# pdf(file=paste0('output/', project.name, '/',  project.name, "_cellfreq_forest_clusterL2.healthy_base.glm.betareg.pdf"), width = 5, height = 5)
# forestplot(
#   # labeltext = labeltext,
#   labeltext = cells,
#   # graph.pos = 3,
#   mean = df.forest$estimate,
#   lower = df.forest$lower,
#   upper = df.forest$upper,
#   xlab = "Coefficient",
#   zero = 0,
#   line.margin = .15,
#   col=fpColors(box="royalblue", lines="darkblue", zero = "gray50")
# )
# dev.off()


write_csv(df.res, paste0('output/', project.name, '/',  project.name, '_queryL2_base_post_glm_betareg.csv'))
write_csv(df.res.disease, paste0('output/', project.name, '/',  project.name, '_queryL2_healthy_base_glm_betareg.csv'))
write_csv(df.output.wide, paste0('output/', project.name, '/',  project.name, '_queryL2_output.csv'))


max_freq_data <- df.output %>%
  group_by(clusterL2) %>%
  summarise(y_max = max(freq*110))

signif_data <- df.res %>%
  left_join(max_freq_data, by = "clusterL2") %>%
  mutate(p.value = round(p.value, 3),
         p.value = ifelse(p.value == 0, "<1e-3", p.value)) %>%
  mutate(fdr = round(fdr, 3),
         fdr = ifelse(fdr == 0, "<1e-3", fdr)) %>%
  mutate(x_min = "Baseline", x_max = "Follow_up")

max_freq_data_disease <- df.output %>%
  group_by(clusterL2) %>%
  summarise(y_max = max(freq*135))

signif_data_disease <- df.res.disease %>%
  left_join(max_freq_data_disease, by = "clusterL2") %>%
  mutate(p.value = round(p.value, 3),
         p.value = ifelse(p.value == 0, "<1e-3", p.value)) %>%
  mutate(fdr = round(fdr, 3),
         fdr = ifelse(fdr == 0, "<1e-3", fdr)) %>%
  mutate(x_min = "Baseline", x_max = "Healthy")

df.output %>%
  mutate(pct_freq = freq*100) %>%
  ggplot(aes(condition, pct_freq)) +  
  geom_line(aes(group = donor)) +
  geom_point(aes(size = n)) +  
  facet_wrap(~clusterL2, scales = "free_y", nrow = 3) +  
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  theme_minimal() +
  theme(legend.position = "left",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size(range = c(2,12), breaks = 10^c(1,2,3,4)) + 
  geom_signif(data = as.data.frame(signif_data), 
                        aes(xmin = x_min, xmax = x_max, annotations = p.value, y_position = log10(y_max)),
                        textsize = 3, vjust = -0.2,
                        manual = TRUE,
                        tip_length = 0) + 
  geom_signif(data = as.data.frame(signif_data_disease), 
              aes(xmin = x_min, xmax = x_max, annotations = p.value, y_position = log10(y_max)),
              textsize = 3, vjust = -0.2,
              manual = TRUE,
              tip_length = 0)
ggsave(paste0('output/', project.name, '/',  project.name, '_queryL2_freq.pvalue.pdf'), width = 8, height = 10)

df.output %>%
  mutate(pct_freq = freq*100) %>%
  ggplot(aes(condition, pct_freq)) +  
  geom_line(aes(group = donor)) +
  geom_point(aes(size = n)) +  
  facet_wrap(~clusterL2, scales = "free_y", nrow = 3) +  
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  theme_minimal() +
  theme(legend.position = "left",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size(range = c(2,12), breaks = 10^c(1,2,3,4)) + 
  geom_signif(data = as.data.frame(signif_data), 
              aes(xmin = x_min, xmax = x_max, annotations = fdr, y_position = log10(y_max)),
              textsize = 3, vjust = -0.2,
              manual = TRUE,
              tip_length = 0) + 
  geom_signif(data = as.data.frame(signif_data_disease), 
              aes(xmin = x_min, xmax = x_max, annotations = fdr, y_position = log10(y_max)),
              textsize = 3, vjust = -0.2,
              manual = TRUE,
              tip_length = 0)
ggsave(paste0('output/', project.name, '/',  project.name, '_queryL2_freq.fdr.pdf'), width = 8, height = 10)

#### NMF

df.metadata <- read_csv(paste0('data/', project.name,'/', project.name,'_metadata.csv'))
df.output <- read_csv(paste0('output/', project.name, '/',  project.name, '_CD4T_NMF_arranged.csv'))

df.output <- df.output %>% 
  left_join(df.metadata, by = "sample", suffix = c("", ""))

df_MS <- df.output %>% 
  filter(condition != "Healthy")

list.res <- list()
for (cell in unique(df.output$clusterL2)) {
  for (comp in colnames(df.output) %>% str_subset("NMF_")){
    res <- glm(as.formula(paste0(comp, " ~ condition + donor")),
               data = df_MS %>% filter(clusterL2 == cell))
    # Use lmerTest's summary to get p-values
    res_summary <- summary(res)
    list.res[[paste(cell, comp)]] <- as_tibble(res_summary$coefficients, rownames = "var") %>%
      mutate(cell = cell, component = comp)
  }
}

df.res <- purrr::map_dfr(list.res, ~ .x %>% dplyr::filter(var == "conditionpost")) %>%
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|t|)`))


df.plot <- df.res %>%
  mutate(
    component = factor(component, levels = c(
      "NMF_0", "NMF_1", "NMF_2", "NMF_3", "NMF_4", "NMF_5", "NMF_6", "NMF_7", "NMF_8", "NMF_9",
      "NMF_10", "NMF_11"
    ))
  ) %>%
  mutate(component=fct_recode(component, 
                              'NMF0 Cytotoxic-F'='NMF_0', 'NMF1 Treg-F'='NMF_1', 
                              'NMF2 Th17-F'='NMF_2', 'NMF3 Naive-F'='NMF_3', 
                              'NMF4 Act-F'='NMF_4', 'NMF5 TregEff/Th2-F'='NMF_5', 
                              'NMF6 Tfh-F'='NMF_6', 'NMF7 IFN-F'='NMF_7', 
                              'NMF8 Cent.Mem.-F'='NMF_8', 
                              'NMF9 Thy.Emi.-F'='NMF_9', 
                              'NMF10 Tissue-F'='NMF_10',
                              'NMF11 Th1-F'='NMF_11'))


lim.est <- max(abs(df.plot$Estimate))
lim.padj <- 10^-10

# df.plot <- df.plot %>%
#   mutate(
#     Estimate = if_else(Estimate > lim.est, lim.est, if_else(Estimate < -lim.est, -lim.est, Estimate)),
#     # padj = if_else(padj > lim.padj, lim.padj, if_else(padj > 0.1, 0, padj))
#   )

dotplot <- ggplot(df.plot, aes(x = component, y = cell, size = -log10(padj), color = Estimate)) +
  geom_point() +
  scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(-lim.est, lim.est)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_size(range = c(0, 4), limits = c(0, 10)) +
  xlab("Component") +
  ylab("Cell Type") +
  labs(color = "Coef.") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

print(dotplot)
ggsave(paste0('output/', project.name, '/',  project.name, "dotplot.NMF.clusterL2.base_post.GLM.pdf"), dotplot, width = 5, height = 5)


write_csv(df.res, paste0('output/', project.name, '/',  project.name, '_NMF_queryL2_base_post_glm_withproj.csv'))


### disease effect
df_MS <- df.output %>% 
  filter(condition != "post")
df_MS$condition <- factor(df_MS$condition, levels = c("Healthy", "baseline"))

list.res <- list()
for (cell in unique(df.output$clusterL2)) {
  for (comp in colnames(df.output) %>% str_subset("NMF_")){
    res <- glm(as.formula(paste0(comp, " ~ condition")),
               data = df_MS %>% filter(clusterL2 == cell))
    # Use lmerTest's summary to get p-values
    res_summary <- summary(res)
    list.res[[paste(cell, comp)]] <- as_tibble(res_summary$coefficients, rownames = "var") %>%
      mutate(cell = cell, component = comp)
  }
}

df.res <- purrr::map_dfr(list.res, ~ .x %>% dplyr::filter(var == "conditionbaseline")) %>%
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|t|)`))


df.plot <- df.res %>%
  mutate(
    component = factor(component, levels = c(
      "NMF_0", "NMF_1", "NMF_2", "NMF_3", "NMF_4", "NMF_5", "NMF_6", "NMF_7", "NMF_8", "NMF_9",
      "NMF_10", "NMF_11"
    ))
  ) %>%
  mutate(component=fct_recode(component, 
                              'NMF0 Cytotoxic-F'='NMF_0', 'NMF1 Treg-F'='NMF_1', 
                              'NMF2 Th17-F'='NMF_2', 'NMF3 Naive-F'='NMF_3', 
                              'NMF4 Act-F'='NMF_4', 'NMF5 TregEff/Th2-F'='NMF_5', 
                              'NMF6 Tfh-F'='NMF_6', 'NMF7 IFN-F'='NMF_7', 
                              'NMF8 Cent.Mem.-F'='NMF_8', 
                              'NMF9 Thy.Emi.-F'='NMF_9', 
                              'NMF10 Tissue-F'='NMF_10',
                              'NMF11 Th1-F'='NMF_11'))


lim.est <- max(abs(df.plot$Estimate))
lim.padj <- 10^-10

# df.plot <- df.plot %>%
#   mutate(
#     Estimate = if_else(Estimate > lim.est, lim.est, if_else(Estimate < -lim.est, -lim.est, Estimate)),
#     # padj = if_else(padj > lim.padj, lim.padj, if_else(padj > 0.1, 0, padj))
#   )

dotplot <- ggplot(df.plot, aes(x = component, y = cell, size = -log10(padj), color = Estimate)) +
  geom_point() +
  scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(-lim.est, lim.est)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_size(range = c(0, 4), limits = c(0, 10)) +
  xlab("Component") +
  ylab("Cell Type") +
  labs(color = "Coef.") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, colour = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

print(dotplot)
ggsave(paste0('output/', project.name, '/',  project.name, "dotplot.NMF.clusterL2.haelty_base.GLM.pdf"), dotplot, width = 5, height = 5)


write_csv(df.res, paste0('output/', project.name, '/',  project.name, '_NMF_queryL2_healthy_base_glm_withproj.csv'))
