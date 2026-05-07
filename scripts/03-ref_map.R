# 0. Init ----------------------------------------------------------------------
library(here); library(data.table); library(tidyverse)
library(giscoR); library(sf)

setwd(str_remove(here(), '/scripts'))



# 1. Import & prepare reference population data --------------------------------
ref <- fread('reference_composition.txt') %>%
  rename_all(tolower) %>%
  filter(used == 'Yes', !ancestry %in% c('MEA','OCE'))

ref_coord <- fread('igsr_populations.tsv') %>%
  filter(!grepl('Simons', `Data collections`)) %>%
  transmute(data = case_when(`Data collections` == 'Human Genome Diversity Project' ~ 'HGDP',
                             T ~ '1KG'),
            ancestry = case_when(`Superpopulation name` %in% c('Africa (HGDP)', 'African Ancestry') ~ 'AFR',
                                 `Superpopulation name` %in% c('Europe (HGDP)', 'European Ancestry') ~ 'EUR',
                                 `Superpopulation name` %in% c('East Asia (HGDP)', 'East Asian Ancestry') ~ 'EAS',
                                 `Superpopulation name` %in% c('America (HGDP)', 'American Ancestry') ~ 'IAM',
                                 T ~ NA),
            population = case_when(data == '1KG' ~ `Population code`,
                                   T ~ str_remove_all(`Population name`, ' ')),
            latitude = `Population latitude`,
            longitude = `Population longitude`) %>%
  drop_na()

ref <- ref %>%
  left_join(ref_coord) %>%
  mutate(ancestry = factor(ancestry, levels = c('AFR', 'EUR', 'EAS', 'IAM')))



# 2. Reference Population Map --------------------------------------------------
gisco_get_countries() %>%
  filter(ISO3_CODE != 'ATA') %>%
  ggplot() +
  geom_sf(fill = 'white') +
  geom_point(data = ref,
             aes(x = longitude, y = latitude, color = ancestry, size = n)) +
  scale_color_manual(name = 'Reference population', drop = F,
                     values = c('#ffc300', '#219ebc', '#a7c957', '#e76f51'),
                     labels = c(expression(italic(P)^AFR),
                                expression(italic(P)^EUR),
                                expression(italic(P)^EAS),
                                expression(italic(P)^IAM))) +
  scale_size_continuous(name = 'Sample size') +
  theme_void() +
  theme(legend.position = 'bottom', legend.box = 'vertical')

ggsave('figures/appendix_ref_location.pdf', width = 16, height = 12)
