# R code to draw Fig S5, Fig S6 and Fig S7

library("here")
library("tidyverse")

# get command line arguments

args <- commandArgs(trailingOnly = TRUE)

sessioninfo.filename <- args[1]

########################################################################################
# get data

load(file=here("rdata","Ailaoshan_OTU_table.rdata"))

load(file=here("rdata","Ailaoshan_IUCNdata.rdata"))

# model summaries generated with Ailaoshan_model_summary.R
model.summary <- list(
  LSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_LSU_dt_groups.rds")),
  SSU = readRDS(file = here("rds", "Ailaoshan_model_summary_200_final_SSU_dt_groups.rds"))
)
beta.output <- lapply(model.summary, FUN = function (X) lapply(X$beta.output, as.data.frame))
rm(model.summary)

# relabel (non-avian) reptiles as squamates

taxa <- taxa %>% mutate(class = ifelse(class == "Reptiles", "Squamates", class))

########################################################################################
# species ordered by environmental coefficients

# LSU elevation

LSU.plot.elev <- beta.output$LSU$elev %>%
  mutate(`95% BCI excludes zero` = ifelse((`2.5%` > 0) | (`97.5%` < 0), TRUE, FALSE)) %>%
  inner_join(taxa, by = "OTU") %>% mutate(class = tolower(class)) %>%
  left_join(taxa.iucn %>% select(OTU, `IUCN category`), by = "OTU") %>%
  left_join(leech %>% select(OTU, AdultBodyMass_g, domestic) %>% distinct(), by = "OTU") %>%
  # add dummy reptile rows to ensure squamates facet is wide enough
  bind_rows(tibble(consensus.short = paste0("dummy", 1:4), class = rep("squamates",4), mean = c(-2,-2,2,2))) %>%
  ggplot(aes(x=reorder(consensus.short,mean), color = class)) +
  geom_abline(slope = 0, intercept = 0, color = "darkgrey") +
  geom_point(aes(y=mean)) +
  geom_linerange(aes(ymin = `25%`, ymax = `75%`), size = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.5) +
  labs(y = expression(estimated~beta[elev])) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank()) +
  geom_text(aes(y=`97.5%`, label = `IUCN category`), hjust = 0.5, vjust = 0, nudge_y = 0.1, size = 2, colour = "black") +
  facet_grid(~class, scales = "free_x", space = "free")

#

LSU.plot.reserve <- beta.output$LSU$reserve %>%
  mutate(`95% BCI excludes zero` = ifelse((`2.5%` > 0) | (`97.5%` < 0), TRUE, FALSE)) %>%
  inner_join(taxa, by = "OTU") %>% mutate(class = tolower(class)) %>%
  left_join(taxa.iucn %>% select(OTU, `IUCN category`), by = "OTU") %>%
  left_join(leech %>% select(OTU, AdultBodyMass_g, domestic) %>% distinct(), by = "OTU") %>%
  # add dummy reptile rows to ensure reptiles facet is wide enough
  bind_rows(tibble(consensus.short = paste0("dummy", 1:4), class = rep("squamates",4), mean = c(-2,-2,2,2))) %>%
  ggplot(aes(x=reorder(consensus.short, mean), color = class)) +
  geom_abline(slope = 0, intercept = 0, color = "darkgrey") +
  geom_point(aes(y=mean)) +
  geom_linerange(aes(ymin = `25%`, ymax = `75%`), size = 1, show.legend = FALSE) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.5) +
  labs(y = expression(estimated~beta[reserve])) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45), axis.title.x = element_blank()) +
  geom_text(aes(y=`97.5%`, label = `IUCN category`), hjust = 0.5, vjust = 0, nudge_y = 0.1, size = 2, colour = "black") +
  facet_grid(~class, scales = "free_x", space = "free")

# plot for publication
SSU.plot.elev
ggsave(here("figures","FigS7_SSU_elev.pdf"), width = 12, height = 6, useDingbats = FALSE)

# same plot with 10kg mammals to assist annotation
SSU.plot.elev + aes(color = (AdultBodyMass_g > 10000))
ggsave(here("figures","FigS7_SSU_elev_10kg.pdf"), width = 12, height = 6, useDingbats = FALSE)

# same plot with domestic mammals to assist annotation
SSU.plot.elev + aes(color = domestic)
ggsave(here("figures","FigS7_SSU_elev_domestic.pdf"), width = 12, height = 6, useDingbats = FALSE)

# same plot highlighting species whose 95% BCI excludes zero to assist annotation
SSU.plot.elev + aes(color = `95% BCI excludes zero`)
ggsave(here("figures","FigS7_SSU_elev_BCI.pdf"), width = 12, height = 6, useDingbats = FALSE)

########################################################################################

# session info
writeLines(capture.output(sessionInfo()), sessioninfo.filename)