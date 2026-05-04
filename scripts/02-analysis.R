# 0. Init ----------------------------------------------------------------------
library(here); library(tidyverse); library(cowplot); library(fixest)
library(flextable); library(modelsummary); library(ftExtra); library(caret)
library(scales); library(lemon); library(lme4); library(marginaleffects)
library(psych); library(lavaan); library(ggpubr); library(ggh4x)
library(gtsummary); library(ggtern); library(colorjam)
library(xgboost); library(survey); library(car)

available_cores <- parallel::detectCores(logical = T)
n_cores <- if (is.na(available_cores) || available_cores < 2) 1L else available_cores - 1L

setwd(str_remove(here(), '/scripts'))

race_label <- c('NHW', 'NHB', 'AINA', 'AAPI', 'Hispanic')

theme_custom <- function() {
  theme_classic() %+replace%
    theme(legend.position = 'top',
          plot.margin = unit(c(5.5,10,5.5,5.5), 'pt'),
          plot.title.position = 'plot', 
          plot.caption.position = 'plot',
          strip.background = element_blank(), 
          strip.text = element_text(face='bold', size=rel(0.8), margin=margin(0,0,2,0)),
          plot.title = element_text(face='bold', size=rel(1.2), hjust=0,
                                    margin=margin(0,0,5.5,0)))
}



# 1. Import & prepare data -----------------------------------------------------
genanc <- readRDS('all_stats.rds')
pheno <- readRDS('all_pheno.rds')

s4 <- genanc %>%
  filter(admixture=='s4') %>% 
  select(aid, matches('X[1-4]_b')) %>%
  rename_with(function(x) x = str_remove(x,'_b'), matches('X[1-4]_b')) %>%
  left_join(pheno)

# get the national means and SDs for standardization purpose (using the full phenotypic sample)
td <- svydesign(data = pheno %>% drop_na(wt4, clusterid4, strataid4, attract4),
                id = ~ clusterid4, strata = ~ strataid4, weights = ~ wt4)
nmean <- svymean(~ attract4, design = td)[[1]]
nsd <- sqrt(svyvar(~ attract4, design = td))[[1]]

fd <- s4 %>%
  rename(intrace3 = i3racec, intrace4 = i4racec, intbyear3 = i3byear, intbyear4 = i4byear) %>%
  rename_with(~ gsub('^(i)(\\d)(.*)$', 'int\\3\\2', .), 
              matches('i[34](female|edu|byear|spksp|ahexp)')) %>%
  pivot_longer(cols = matches('^age|^intrace|^intfemale|^intbyear|^intedu|^intbyear|^intspksp|^intahexp|^intid|^attract|^irace'),
               names_to = c('.value', 'wave'),
               names_pattern = '^(age|intrace|intfemale|intedu|intbyear|intspksp|intahexp|intid|attract|irace)(\\d+)$',
               values_drop_na = TRUE,
               values_transform = list(age = as.numeric)) %>%
  arrange(aid) %>%
  mutate(female = as.numeric(sex=='Female'),
         age2 = age*age, agey = round(age/12),
         agey = case_when(agey<=13 ~ 13, agey>=32 ~ 32, T ~ agey),
         # dichotomous variables
         attract2 = as.numeric(attract==5),
         attract = (attract-nmean)/nsd,
         intrace = relevel(factor(intrace), ref = 'NHW'),
         same_race = as.numeric(as.character(intrace)==as.character(best_race)),
         same_race2 = case_when(intrace=='NHW' ~ 'White', same_race==1 ~ 'Target',
                                same_race==0 ~ 'Other', T ~ NA),
         int2gender = case_when(intfemale==1 & female==1 ~ 'FF',
                                intfemale==1 & female==0 ~ 'FM',
                                intfemale==0 & female==1 ~ 'MF',
                                intfemale==0 & female==0 ~ 'MM',
                                T ~ NA),
         irace_wave = ifelse(is.na(irace)|is.na(wave), NA, paste(irace,wave)) %>% as.factor()) %>%
  group_by(wave) %>%
  mutate(intage = ifelse(intbyear>=median(intbyear, na.rm = T), 0, 1),
         intage = case_when(wave=='3' ~ iyear3 - intbyear,
                            wave=='4' ~ iyear4 - intbyear),
         intage = ifelse(intage<55, 0, 1)) %>%
  ungroup()

fd <- fd %>%
  drop_na(attract, female, best_race, agey, intid, wave) %>% # n = 35611
  group_by(intid, wave) %>%
  filter(n()>1) %>%
  group_by(agey) %>%
  filter(n()>1) %>%
  group_by(intid, wave) %>%
  filter(n()>1) %>%
  ungroup() # 35508



# 2. Descriptive tables & figures ----------------------------------------------
# 2.1 Respondents' descriptive statistics
set_gtsummary_theme(list('style_number-arg:big.mark' = ''))

td <- fd %>%
  mutate(best_race = case_when(best_race %in% c('NHB','Hispanic') ~ best_race,
                               T ~ 'Other Races'),
         best_race = factor(best_race, levels = c('NHB','Hispanic','Other Races'))) %>%
  select(aid, best_race, sex, byear, edu, attract, attract2, skin, hair, eye, 
         fam_ses, nhood_ses, paste0('X',1:4))

tbl <- lapply(c('NHB','Hispanic','Other Races'), function(r) {
  td %>%
    select(-c(attract, attract2)) %>%
    unique() %>%
    filter(best_race==r) %>%
    select(paste0('X',1:4), sex, byear, edu, skin, hair, eye, fam_ses, nhood_ses) %>%
    tbl_summary(
      statistic = list(all_continuous() ~ '{mean}--{sd}',
                       all_categorical() ~ '{p}%-- '),
      digits = everything() ~ 2,
      label = list(sex ~ 'Gender',
                   byear ~ 'Birth year',
                   edu ~ 'Educational attainment',
                   skin ~ 'Skin tone',
                   hair ~ 'Hair color',
                   eye ~ 'Eye color',
                   fam_ses ~ 'W1 social origins score',
                   nhood_ses ~ 'W1 neighborhood disadvantage',
                   X1 ~ '*P*^AFR^',
                   X2 ~ '*P*^EUR^',
                   X3 ~ '*P*^EAS^',
                   X4 ~ '*P*^IAM^')                            
    ) %>%
    # do not show missingness rows
    remove_row_type(type = 'missing') %>%
    # add sample sizes
    add_n() %>%
    # remove footnote
    modify_footnote(c(all_stat_cols()) ~ NA) %>%
    modify_header(stat_0 ~ '**Mean (SD) / %**')
})

tbl_add <- lapply(c('NHB','Hispanic','Other Races'), function(r) {
  nind <- td %>%
    filter(best_race==r) %>%
    drop_na(attract, attract2) %>%
    select(aid) %>%
    unique() %>% 
    nrow()
  
  t <- td %>%
    filter(best_race==r) %>%
    select(attract, attract2) %>%
    tbl_summary(
      type = list(attract ~ 'continuous'),
      statistic = list(all_continuous() ~ '{mean}--{sd}',
                       all_categorical() ~ '{p}%-- '),
      digits = everything() ~ 2,
      label = list(attract ~ 'Attrativeness',
                   attract2 ~ 'Very attractivess')                            
    ) %>%
    # do not show missingness rows
    remove_row_type(type = 'missing') %>%
    # add sample sizes
    add_n()
  
  t$table_body %>%
    mutate(n = paste0(n, '\n(N = ', nind, ')')) %>%
    return()
})

tbl <- map2(tbl, tbl_add, function(x,y) {
  x %>%
    modify_table_body(
      ~.x[1:4,] %>%
        rbind(y) %>%
        rbind(.x[-c(1:4),])
    )
})

wtbl <- tbl[[1]] %>%
  modify_table_body(
    ~.x %>%
      full_join(tbl[[2]]$table_body %>% rename(n1 = n, stat_1 = stat_0)) %>%
      full_join(tbl[[3]]$table_body %>% rename(n2 = n, stat_2 = stat_0)) %>%
      separate(stat_0, c('stat_0a','stat_0b'), '--') %>%
      separate(stat_1, c('stat_1a','stat_1b'), '--') %>%
      separate(stat_2, c('stat_2a','stat_2b'), '--')
  ) %>% 
  modify_header(n1 = '**N**', stat_0a = '**Mean / %**', stat_0b = '**SD**',
                stat_1a = '**Mean / %**', stat_1b = '**SD**', n2 = '**N**',
                stat_2a = '**Mean / %**', stat_2b = '**SD**') %>%
  modify_spanning_header(c(n, stat_0a, stat_0b) ~ '**NHB**', 
                         c(n1, stat_1a, stat_1b) ~ '**Hispanic**',
                         c(n2, stat_2a, stat_2b) ~ '**Other Races**') %>%
  as_flex_table() %>%
  colformat_md() %>%
  autofit()

save_as_html(wtbl, path = 'tables/appendix_desc_respondents.html')


# 2.2 Interviewers' descriptive statistics
td <- fd %>%
  transmute(best_race, intid = paste0('wave',wave,'_',intid),
            intage = case_when(wave=='3' ~ iyear3-intbyear, 
                               wave=='4' ~ iyear4-intbyear),
            intfemale, intrace, intedu, n = 1, 
            n2 = as.numeric(best_race=='NHB')) %>%
  drop_na() %>%
  group_by(intid) %>%
  mutate(intrace = factor(intrace, levels = c('NHW','NHB','AINA','AAPI','Hispanic','Other')),
         n = sum(n), n2 = sum(n2)) %>%
  ungroup()

tbl_unwt <- td %>%
  select(-best_race) %>%
  group_by(intid) %>%
  mutate(intage = mean(intage)) %>%
  unique() %>%
  tbl_summary(
    include = !intid,
    statistic = list(all_categorical() ~ '{p}%',
                     all_continuous() ~ '{mean} ({sd})'),
    digits = everything() ~ 0,
    label = list(n ~ 'Number of interviewees',
                 n2 ~ 'Number of NHB interviewees',
                 intage ~ "Interviewer's age",
                 intfemale ~ 'Female interviewer',
                 intrace ~ "Interviewer's race",
                 intedu ~ "Interviewer's educational attainment")
  )

tbl_wt <- td %>%
  filter(best_race=='NHB') %>%
  tbl_summary(
    include = !c(intid, best_race, n, n2), # wt
    statistic = list(all_categorical() ~ '{p}%',
                     all_continuous() ~ '{mean} ({sd})'),
    digits = everything() ~ 0,
    label = list(intage ~ "Interviewer's age",
                 intfemale ~ "Female interviewer",
                 intrace ~ "Interviewer's race",
                 intedu ~ "Interviewer's educational attainment")
  )

tbl <- tbl_unwt %>%
  modify_table_body(
    ~.x %>% 
      full_join(tbl_wt$table_body %>%
                  rename(stat_1 = stat_0)) %>%
      arrange(variable!='n', variable!='n2')
  ) %>%
  modify_footnote(c(all_stat_cols()) ~ NA) %>%
  modify_header(stat_0 = 'Unweighted', stat_1 = 'Weighted') %>%
  as_flex_table() %>%
  add_header_row(values = c('Characteristic','Mean (SD) / %'), colwidths = c(1,2)) %>%
  add_header_row(values = c('Characteristic',paste0('N = ',length(unique(td$intid)))),
                 colwidths = c(1,2)) %>%
  merge_at(i = 1:3, j = 1, part = 'header') %>%
  autofit()

save_as_html(tbl, path = 'tables/appendix_desc_interviewers.html')


# 2.3 Histograms of interview gaps (days between interviews per interviewer)
td <- pheno %>%
  select(aid, best_race, matches('i(ntid|year|month|day)[1-4]'), 
         intrace3 = i3racec, intrace4 = i4racec) %>%
  pivot_longer(cols = !c(aid, best_race),
               names_to = c('.value', 'wave'),
               names_pattern = '^(intid|intrace|iyear|imonth|iday)(\\d+)$',
               values_drop_na = TRUE) %>%
  mutate(wave = paste('Wave', wave),
         best_race = case_when(best_race %in% c('AINA','AAPI') ~ 'Other', T ~ best_race),
         best_race = factor(best_race, levels = c('NHW','NHB','Hispanic','Other')),
         intrace = factor(intrace, levels = c('NHW','NHB','Hispanic','Other')))

pd <- td %>%
  mutate(itime = iday + imonth*100 + iyear*10000,
         itime = as.Date(as.character(itime), format = '%Y%m%d')) %>%
  drop_na(intid) %>%
  group_by(wave, intid, intrace) %>%
  arrange(wave, intid, intrace, itime) %>%
  mutate(gap = as.numeric(difftime(itime, lag(itime), units = 'days'))) %>%
  drop_na(gap) # N = 64210
pd <- pd %>%
  mutate(gap = ifelse(gap>=10, 10, gap)) %>%
  group_by(wave, gap) %>%
  summarise(n = n()) %>%
  group_by(wave) %>%
  mutate(frac = n/sum(n)) %>%
  ungroup()

pd %>%
  ggplot(aes(x = gap, y = frac, fill = wave)) +
  geom_bar(stat = 'identity', color = 'black', width = 1, show.legend = F) +
  facet_wrap(~ wave, scales = 'free') +
  scale_x_continuous('Days between interviews', breaks = 0:10, labels = c(0:9,'10+'),
                     expand = c(0,0)) +
  scale_y_continuous('Fraction', limits = c(0,0.4), breaks = seq(0,0.4,0.1),
                     expand = c(0,0)) +
  scale_fill_viridis_d('') +
  theme_custom()
ggsave('figures/appendix_interviewgap_histogram.pdf', width = 6, height = 6)


# 2.4 Supplementary figure: histogram of attractiveness
pd <- pheno %>%
  select(aid, paste0('attract',1:4)) %>%
  gather(wave, attract, attract1:attract4) %>%
  mutate(wave = str_replace_all(wave, 'attract', 'Wave ')) %>%
  drop_na()
pd <- pd %>%
  mutate(sample = 'Full') %>%
  rbind(pd %>% filter(aid %in% s4$aid) %>% mutate(sample = 'Genotyped')) %>%
  group_by(sample, wave, attract) %>%
  summarise(n = n()) %>%
  group_by(sample, wave) %>%
  mutate(frac = n/sum(n)) %>%
  ungroup()

pd %>%
  ggplot(aes(x = attract, y = frac, fill = sample)) +
  geom_bar(stat = 'identity', position = 'identity', width = 1, color = 'black', alpha = 0.6) +
  facet_wrap(~ wave, scales = 'free') +
  coord_cartesian(expand = F) +
  scale_fill_manual('Sample', values = c('#c1121f','#669bbc')) +
  labs(x = 'Attractiveness scores', y = 'Fraction') +
  theme_custom()
ggsave('figures/appendix_attract_histogram.pdf', width = 5, height = 5)


# 2.5 Figure 1 & supplementary figure: racial disparity in attractiveness
m1 <- feols(attract ~ best_race2 | intid^wave, cluster = ~aid, fd)
m2 <- feols(attract2 ~ best_race2 | intid^wave, cluster = ~aid, fd)
td <- fd %>% 
  drop_na(best_race2) %>%
  mutate(
    race_effect = case_when(
      best_race2=='NHW' ~ 0,
      best_race2=='NHB' ~ m1$coefficients['best_race2NHB'],
      best_race2=='AINA' ~ m1$coefficients['best_race2AINA'],
      best_race2=='AAPI' ~ m1$coefficients['best_race2AAPI'],
      best_race2=='Hispanic' ~ m1$coefficients['best_race2Hispanic']),
    race_effect2 = case_when(
      best_race2=='NHW' ~ 0,
      best_race2=='NHB' ~ m2$coefficients['best_race2NHB'],
      best_race2=='AINA' ~ m2$coefficients['best_race2AINA'],
      best_race2=='AAPI' ~ m2$coefficients['best_race2AAPI'],
      best_race2=='Hispanic' ~ m2$coefficients['best_race2Hispanic']),
    # add the race effects and mean back
    attract3 = m1$residuals + race_effect + mean(attract),
    attract4 = m2$residuals + race_effect2 + mean(attract2)
  )
  
my_race_plot <- function(data, yvar, y_label, p_pos, tip_length, title_label, 
                         discrete, sec_y_label = NULL, ylim = c(0,0.18), label_digits = 2,
                         comparisons = list(c('NHW','NHB'), c('AAPI','NHB'), c('Hispanic','NHB'))) {
  data$y <- data[[yvar]]
  label_fmt <- paste0('%.', label_digits, 'f')
  value <- data %>% 
    filter(best_race!='AINA') %>%
    group_by(best_race) %>%
    summarise(n = n(), se = sd(y, na.rm = T)/sqrt(n), y = mean(y, na.rm = T)) %>%
    mutate(yp = y + se + 0.005, yl = sprintf(label_fmt, y))
  ref <- value %>% filter(best_race=='NHW') %>% pull(y)
  
  p <- data %>% 
    filter(best_race!='AINA') %>%
    ggbarplot(x = 'best_race', y = yvar, add = 'mean_se', fill = 'best_race') +
    geom_text(data = value, aes(x = best_race, y = yp, label = yl)) +
    geom_hline(yintercept = ref, linetype = 'dashed', color = 'gray') +
    stat_compare_means(method = 't.test', label = 'p.signif', 
                       label.y = p_pos, tip.length = tip_length,
                       symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,Inf),
                                          symbols = c('***','**','*','NS')),
                       comparisons = comparisons) +
    coord_cartesian(ylim = ylim)

  if (discrete) {
    p <- p +
      scale_y_continuous(y_label, breaks = seq(0,0.15,0.05), expand = c(0,0))
  } else {
    p <- p +
      scale_y_continuous(y_label, breaks = seq(0,0.15,0.05), expand = c(0,0),
                         sec.axis = sec_axis(~ .*nsd + nmean, sec_y_label))
  }
  
  p <- p +
    scale_fill_manual(values = c(NHW = '#007ea7', NHB = '#f9a620',
                                 AAPI = '#006400', Hispanic = '#dd0426')) +
    labs(x = '', title = title_label) +
    theme_custom() +
    theme(legend.position = 'none', plot.title.position = 'plot',
          plot.title = element_text(hjust = 0.6, size = 10))
  
  return(p)
}

plot_grid(my_race_plot(td, 'attract', 'Average attractiveness score\n(Standardized)', 
                       c(-0.1,-0.09,-0.08), 0.001, 'Unadjusted', F),
          my_race_plot(td, 'attract3', '', c(-0.32,-0.31,-0.3), 0.0005,
                       'Residualized', F, '(Unstandardized)'),
          my_race_plot(td, 'attract2', 'Fraction "very attractive"', c(0.1,0.11,0.12), 0.005, 'Unadjusted', T),
          my_race_plot(td, 'attract4', '', c(0.05,0.06,0.07), 0.0025, 'Residualized', T),
          align = 'v', nrow = 2, labels = c('A','B','C','D'))
ggsave('figures/figure1.pdf', width = 10, height = 9)

my_samerace_plot <- function(data, yvar, y_label, p_pos, discrete) {
  data$y <- data[[yvar]]
  value <- data %>% 
    filter(best_race!='AINA') %>%
    filter(best_race %in% c('NHW','NHB','Hispanic'), !is.na(same_race)) %>%
    mutate(same_race = factor(same_race, labels = c('Different', 'Same'))) %>%
    group_by(best_race, same_race) %>%
    summarise(n = n(), se = sd(y, na.rm = T)/sqrt(n), y = mean(y, na.rm = T)) %>%
    mutate(yp = ifelse(y < 0, y - se - 0.005, y + se + 0.005), yl = sprintf('%.3f', y))
  
  p <- data %>% 
    filter(best_race %in% c('NHW','NHB','Hispanic'), !is.na(same_race)) %>%
    mutate(same_race = factor(same_race, labels = c('Different', 'Same'))) %>%
    ggbarplot(x = 'same_race', y = yvar,  add = 'mean_se', fill = 'same_race',
              position = position_dodge(0.7), width = 1,
              facet.by = 'best_race', strip.position = 'bottom') +
    geom_text(data = value, aes(x = same_race, y = yp, label = yl)) +
    stat_compare_means(method = 't.test', label = 'p.signif',
                       label.y = p_pos, label.x = 1.5,
                       symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,Inf),
                                          symbols = c('***','**','*','NS')))
  
  if (discrete) {
    p <- p +
      scale_y_continuous(y_label)
  } else {
    p <- p +
      scale_y_continuous(y_label, sec.axis = sec_axis(~ .*nsd + nmean, '(Unstandardized)'))
  }
  
  p <- p +
    scale_fill_manual('Interviewer-interviewee race', values = c('#da627d','#450920')) +
    labs(x = '') +
    theme_custom() +
    theme(strip.placement = 'outside', 
          panel.spacing = unit(0,'lines'),
          strip.text = element_text(face = 'plain'),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_grid(
  plot_grid(my_race_plot(td %>% 
                           filter(best_race %in% c('NHW','NHB','Hispanic'),
                                  !is.na(same_race)),
                         'attract', 'Average attractiveness score\n(Standardized)', 
                         c(-0.1,-0.09,-0.08), 0.001, NULL, F, '(Unstandardized)',
                         comparisons = list(c('NHW','NHB'), c('Hispanic','NHB')),
                         ylim = c(-0.03,0.18), label_digits = 3),
            my_race_plot(td %>% 
                            filter(best_race %in% c('NHW','NHB','Hispanic'),
                                   !is.na(same_race)),
                          'attract2', 'Fraction "very attractive"', c(0.1,0.11,0.12),
                         0.005, NULL, T, comparisons = list(c('NHW','NHB'), c('Hispanic','NHB')),
                         ylim = c(-0.03,0.18), label_digits = 3),
            labels = c('A','B')),
  plot_grid(
    get_plot_component(my_samerace_plot(td, 'attract2', 'Fraction "very attractive"', 0.115, T),
                       'guide-box-top') %>% ggdraw(),
    plot_grid(my_samerace_plot(td, 'attract', 'Average attractiveness score\n(Standardized)', 0.08, F) + 
                theme(legend.position = 'none'),
              my_samerace_plot(td, 'attract2', 'Fraction "very attractive"', 0.125, T) + 
                theme(legend.position = 'none'),
              nrow = 1, labels = c('C','D')),
    nrow = 2, rel_heights = c(1,5)
  ),
  align = 'v', nrow = 2
)
ggsave('figures/appendix_attract_by_raceconcord.pdf', width = 10, height = 9)


# 2.6 Supplementary figure: distribution of GSPs by race
s4 %>%
  filter(best_race %in% c('NHW','NHB','Hispanic'), X3 <= 0.05) %>%
  ggtern(aes(x = X1, y = X2, z = X4, color = best_race)) +
  geom_point(alpha = 0.2, size = 3, show.legend = F) +
  facet_wrap(~ best_race) +
  scale_color_manual(name = '', values = c(NHW = '#007ea7', NHB = '#f9a620',
                                           Hispanic = '#dd0426')) +
  labs(x = expression(italic(P)^AFR), y = expression(italic(P)^EUR),
       z = expression(italic(P)^IAM)) +
  theme_custom()
ggsave('figures/appendix_gsp_ternary_easfilter.pdf', width = 15, height = 5)

ref <- s4 %>%
  filter(best_race=='NHW') %>%
  select(X1, X2, X4) %>%
  summarise_all(mean) %>%
  gather(anc, X, c(X1,X2,X4)) %>%
  mutate(anc = factor(anc, labels = c('AFR','EUR','IAM')))
s4 %>%
  filter(best_race %in% c('NHB','Hispanic')) %>%
  gather(anc, X, c(X1,X2,X4)) %>%
  mutate(anc = factor(anc, labels = c('AFR','EUR','IAM'))) %>%
  ggplot() +
  geom_histogram(aes(x = X, y = 0.05*..density.., fill = best_race),  alpha = 0.5,
                 position = 'identity', binwidth = 0.05, color = 'black', center = 0.025) +
  geom_vline(data = ref, aes(xintercept = X), color = '#007ea7', linetype = 'dashed') +
  facet_wrap(~ anc, scales = 'free',
             labeller = label_bquote(cols = bold(italic(P)^.(as.character(anc))))) +
  scale_x_continuous('Genetic similarity proportions', expand = c(0,0),
                     limit = c(0,1), breaks = seq(0,1,0.1)) +
  facetted_pos_scales(
    y = list(
      anc=='AFR' ~ scale_y_continuous('Fraction of sample', expand = c(0,0),
                                      limits = c(0,0.7), breaks = seq(0,0.7,0.1)),
      anc=='EUR' ~ scale_y_continuous('Fraction of sample', expand = c(0,0),
                                      limits = c(0,0.3), breaks = seq(0,0.3,0.05)),
      T ~ scale_y_continuous('Fraction of sample', expand = c(0,0),
                             limits = c(0,1), breaks = seq(0,1,0.1))
    )
  ) +
  scale_fill_manual(name = '', values = c(NHB = '#f9a620', Hispanic = '#dd0426')) +
  labs(y = 'Fraction of sample') +
  theme_custom() +
  theme(legend.position = 'top', panel.spacing.y = unit(0.25, 'cm'),
        strip.text.y = element_blank(), strip.placement = 'outside')
ggsave('figures/appendix_gsp_histo.pdf', width = 10, height = 5)


# 2.7 Figure 2 & supplementary figure: GSP-attractiveness binned scatter plots
# define a function for preparing data for attractiveness plots
my_attract_plot_data <- function(data, nhw_focalX, race_grp, discrete) {
  # underlying unbinned data
  if (!discrete) {
    pd <- data %>%
      filter(best_race==race_grp) %>%
      transmute(aid = paste0(aid,'_',wave), best_race, focalX,
                y = attract,
                ytype = 'Continuous attractiveness')
  } else {
    pd <- data %>%
      filter(best_race==race_grp) %>%
      transmute(aid = paste0(aid,'_',wave), best_race, focalX,
                attract = factor(6 - attract*nsd + nmean, 
                                 labels = c('Very attractive','Attractive',
                                            'About average','Unattractive', 'Very unattractive')))
    pd <- pd %>%
      cbind(predict(dummyVars('~ attract', pd, sep = '_'), newdata = pd) %>% 
              as.data.frame()) %>%
      gather(ytype, y, starts_with('attract_')) %>% 
      mutate(ytype = str_remove_all(ytype, 'attract_'),
             ytype = factor(ytype, levels = levels(pd$attract))) %>%
      filter(ytype!='Very unattractive')
  }
  
  # binned data
  pd2 <- pd %>%
    group_by(ytype) %>%
    arrange(ytype, focalX) %>%
    group_modify(~{
      n <- round(nrow(.x)/100)
      breaks <- quantile(.x$focalX, probs = seq(0,1,length.out=n+1))
      .x$bin <- cut(.x$focalX, breaks = breaks, labels = 1:n)
      .x$bin <- ifelse(is.na(.x$bin),1,.x$bin)
      return(.x)
    }) %>%
    group_by(ytype, bin) %>%
    summarise(across(c(focalX,y), mean)) %>%
    ungroup()
  if (discrete) {
    pd2 <- pd2 %>%
      group_by(bin) %>%
      arrange(bin, ytype) %>%
      mutate(y = cumsum(y)) %>%
      ungroup()
  }
  
  # fitted smoothed data
  pd3 <- pd %>%
    mutate(x_end = max(pd2$focalX))
  if (discrete) {
    pd3 <- pd3 %>%
      group_by(aid) %>%
      arrange(aid, ytype) %>%
      mutate(y = cumsum(y)) %>%
      ungroup()
  }
  if (min(pd3$focalX) > 0.25) {
    pd3 <- pd3 %>%
      add_row(expand_grid(focalX = seq(0,0.5,0.02), ytype = unique(pd3$ytype)))
  }
  pd3 <- pd3 %>%
    mutate(nhw = 0) %>% 
    add_row(
      expand_grid(focalX = nhw_focalX, ytype = unique(pd3$ytype)) %>%
        mutate(nhw = 1)
    ) %>%
    group_by(ytype) %>%
    group_modify(~{
      m1 <- lm(y ~ focalX, data = .x %>% filter(nhw==0))
      m2 <- nls(y ~ b0 + b1 * exp(focalX), data = .x %>% filter(nhw==0))
      .x$linear <- predict(m1, newdata = .x)
      .x$expo <- predict(m2, newdata = .x)
      .x <- .x %>% filter(focalX<=x_end | is.na(x_end))
    }) %>%
    ungroup() %>%
    gather(fit_method, hat_y, c(linear,expo))
  
  # predicted NHW values
  pd4 <- pd3 %>% filter(nhw==1) %>% select(-nhw)
  pd3 <- pd3 %>% filter(nhw==0) %>% select(-nhw) 
  
  return(list(pd1 = pd, pd2 = pd2, pd3 = pd3, pd4 = pd4))
}

# data for reference points
my_attract_plot_ref <- function(data, race_grp) {
  data %>%
    mutate(attract4 = as.numeric(attract*nsd + nmean>=4),
           attract3 = as.numeric(attract*nsd + nmean>=3)) %>%
    filter(best_race %in% c('NHW',race_grp)) %>%
    group_by(best_race) %>%
    summarise_at(c('focalX','attract',paste0('attract',2:4)), mean) %>%
    gather(ytype, y, c(attract,attract2,attract4,attract3)) %>%
    mutate(ytype = case_when(ytype=='attract2' ~ 'Very attractive',
                             ytype=='attract4' ~ 'Attractive',
                             ytype=='attract3' ~ 'About average',
                             T ~ 'Continuous attractiveness'))
}

wbdot <- my_attract_plot_ref(fd %>% rename(focalX = X1), 'NHB')
nhw_X1 <- wbdot %>%
  filter(best_race=='NHW') %>%
  pull(focalX) %>%
  unique()
pd <- my_attract_plot_data(fd %>% rename(focalX = X1), nhw_X1, 'NHB', F)
pd2 <- my_attract_plot_data(fd %>% rename(focalX = X1), nhw_X1, 'NHB', T) %>%
  lapply(function(d) {
    d %>% 
      filter(ytype %in% c('Very attractive','Attractive','About average'))
  })

# define a function for creating attractiveness plots
my_attract_plot <- function(pdlist, refdt, discrete, x_label, anc_color,
                            rug_y_override = NULL, rug_yend_override = NULL,
                            yref_override = NULL) {
  if (discrete) {
    tdod <- refdt %>% 
      filter(ytype!='Continuous attractiveness') %>%
      mutate(ytype = paste(ytype, '- ref'))
    rug_y <- 1.008; rug_yend <- 1.04
    yref <- seq(0,1,0.2)
    ylimit <- c(0,rug_yend)
    y_label <- 'Fraction with a given attractiveness level (or higher)'
    legend_pos <- c(1,0.7)
  } else {
    tdod <- refdt %>% 
      filter(ytype=='Continuous attractiveness') %>%
      mutate(ytype = paste(ytype, '- ref'))
    rug_y <- 0.76; rug_yend <- 0.8
    yref <- seq(-0.5,0.75,0.25)
    ylimit <- c(-0.5,rug_yend)
    y_label <- 'Attractiveness score\n(Standardized)'
    legend_pos <- 'none'
  }
  if (!is.null(rug_y_override)) rug_y <- rug_y_override
  if (!is.null(rug_yend_override)) rug_yend <- rug_yend_override
  if (!is.null(yref_override)) yref <- yref_override
  ylimit <- c(ylimit[[1]], rug_yend)
  ybreaks <- yref
  
  ytype_levels <- c(as.character(pdlist$pd2$ytype), tdod$ytype) %>%
    unique() %>%
    sort()
  p <- pdlist$pd2 %>%
    ggplot(aes(x = focalX, y = y)) +
    geom_point(aes(shape = ytype), size = 2, alpha = 0.4, show.legend = F) +
    geom_line(data = pdlist$pd3, aes(y = hat_y, linetype = ytype, alpha = fit_method), 
              color = anc_color, linewidth = 1) +
    geom_point(data = tdod, aes(shape = ytype, color = best_race), size = 3) +
    geom_point(data = pdlist$pd4, aes(y = hat_y, shape = ytype, alpha = fit_method), 
               size = 2.8, color = '#007ea7', stroke = 1.2) +
    geom_hline(data = tdod %>% filter(best_race=='NHW'), aes(yintercept = y),
               linetype = 'dashed', color = 'gray55') +
    geom_segment(data = pdlist$pd1 %>% 
                   unique() %>% 
                   mutate(rugy = rug_y, rugyend = rug_yend), 
                 aes(x = focalX, xend = focalX, y = rugy, yend = rugyend), 
                 color = anc_color, alpha = 0.1) + 
    geom_hline(yintercept = yref, color = 'gray', linetype = 'dotted', linewidth = 0.3) +
    geom_vline(xintercept = seq(0,1,0.25), color = 'gray', linetype = 'dotted', linewidth = 0.3) +
    scale_x_continuous(x_label, limit = c(-0.01,1), breaks = seq(0,1,0.25), expand = c(0,0))
  
  if (discrete) {
    p <- p +
      scale_y_continuous(y_label, limit = ylimit, breaks = ybreaks, expand = c(0,0))
  } else {
    p <- p +
      scale_y_continuous(y_label, limit = ylimit, breaks = ybreaks, expand = c(0,0),
                         sec.axis = sec_axis(~ .*nsd + nmean, '(Unstandardized)'))
  }
  
  p <- p +
    scale_shape_manual('Interviewer-Rated Attractiveness', 
                       values = c('Continuous attractiveness' = 1,
                                  'Continuous attractiveness - ref' = 16,
                                  'Very attractive' = 1,
                                  'Very attractive - ref' = 16,
                                  'Attractive' = 0,
                                  'Attractive - ref' = 15,
                                  'About average' = 2,
                                  'About average - ref' = 17),
                       labels = c('Very attractive' = 'Very attractive (5)',
                                  'Attractive' = 'Attractive (4)',
                                  'About average' = 'About average (3)'),
                       breaks = c('Very attractive','Attractive','About average')) +
    scale_linetype_manual('Interviewer-Rated Attractiveness', values = c(1,3,5),
                          labels = c('Very attractive (5)','Attractive (4)','About average (3)')) +
    scale_color_manual(values = c(NHB = 'black', Hispanic = 'black', NHW = '#007ea7')) +
    scale_alpha_manual(values = c(linear = 1, expo = 0.4), guide = 'none') +
    theme_custom() +
    theme(legend.position = legend_pos, legend.background = element_blank(),
          legend.title = element_text(hjust = 1), legend.text = element_text(hjust = 1),
          legend.box = 'vertical', legend.justification = 'right',
          legend.key.height = unit(0.4,'cm'), legend.key.width = unit(1.2,'cm'),
          strip.text.y = element_blank(), strip.placement = 'outside') +
    guides(color = 'none', linetype = 'none',
           shape = guide_legend(override.aes = list(color = 'black', size = 1.5, stroke = 1.05))) 

  return(p)
}

attract_plot_legend <- get_plot_component(my_attract_plot(pd, wbdot, F, expression(italic(P)^AFR), '#ffc300') +
                                            scale_alpha_manual('', values = c(linear = 1, expo = 0.4), 
                                                               labels = c(linear = 'Linear', expo = 'Exponential')) +
                                            theme(legend.position = 'top',
                                                  legend.justification = 'center') +
                                            guides(alpha = guide_legend(override.aes = list(shape = 5))),
                                          'guide-box-top') %>% 
  ggdraw()

# figure 2
plot_grid(attract_plot_legend,
          plot_grid(my_attract_plot(pd, wbdot, F, expression(italic(P)^AFR), '#ffc300'),
                    my_attract_plot(pd2, wbdot, T, expression(italic(P)^AFR), '#ffc300'),
                    rel_widths = c(1.1,1), labels = c('A','B')),
          nrow = 2, rel_heights = c(1,7))
ggsave('figures/figure2.pdf', width = 12, height = 6)

nhw_X4 <- my_attract_plot_ref(fd %>% rename(focalX = X4), 'Hispanic') %>%
  filter(best_race=='NHW') %>%
  pull(focalX) %>%
  unique()

attract_plot_legend2 <- get_plot_component(my_attract_plot(pd, wbdot, F, expression(italic(P)^AFR), 'gray30') +
                                            scale_alpha_manual('', values = c(linear = 1, expo = 0.4), 
                                                               labels = c(linear = 'Linear', expo = 'Exponential')) +
                                            theme(legend.position = 'top',
                                                  legend.justification = 'center') +
                                            guides(alpha = guide_legend(override.aes = list(shape = 5))),
                                          'guide-box-top') %>% 
  ggdraw()

# supplementary figure for Hispanic Americans
plot_grid(
  attract_plot_legend2,
  plot_grid(
    my_attract_plot(my_attract_plot_data(fd %>% rename(focalX = X1), 
                                         nhw_X1, 'Hispanic', F),
                    my_attract_plot_ref(fd %>% rename(focalX = X1), 'Hispanic'),
                    F, expression(italic(P)^AFR), '#ffc300'),
    my_attract_plot(my_attract_plot_data(fd %>% rename(focalX = X1), 
                                         nhw_X1, 'Hispanic', T) %>%
                      lapply(function(d) {
                        d %>% 
                          filter(ytype %in% c('Very attractive','Attractive','About average'))
                      }),
                    my_attract_plot_ref(fd %>% rename(focalX = X1), 'Hispanic'),
                    T, expression(italic(P)^AFR), '#ffc300'),
    my_attract_plot(my_attract_plot_data(fd %>% rename(focalX = X4),
                                         nhw_X4, 'Hispanic', F),
                    my_attract_plot_ref(fd %>% rename(focalX = X4), 'Hispanic'),
                    F, expression(italic(P)^IAM), '#e76f51'),
    my_attract_plot(my_attract_plot_data(fd %>% rename(focalX = X4), 
                                         nhw_X4, 'Hispanic', T) %>%
                      lapply(function(d) {
                        d %>% 
                          filter(ytype %in% c('Very attractive','Attractive','About average'))
                      }),
                    my_attract_plot_ref(fd %>% rename(focalX = X4), 'Hispanic'),
                    T, expression(italic(P)^IAM), '#e76f51'),
    nrow = 2, rel_widths = c(1.1,1,1.1,1), labels = c('A','B','C','D')
  ),
  nrow = 2, rel_heights = c(1,14)
)
ggsave('figures/appendix_attract_iam_hispanic.pdf', width = 12, height = 12)


# 2.8 Supplementary figure: attractiveness plot by monoracial vs other NHB
td <- fd %>%
  filter(contain_black, best_race=='NHB')

wbdot <- td %>%
  filter(partial_black==0) %>%
  mutate(best_race = 'NHB', panel = 'Monoracial NHB') %>%
  rbind(td %>% 
          filter(partial_black==1) %>%
          mutate(best_race = 'NHB', panel = 'Other NHB')) %>%
  rbind(fd %>% 
          filter(best_race=='NHW') %>%
          mutate(panel = 'Monoracial NHB')) %>%
  rbind(fd %>% 
          filter(best_race=='NHW') %>%
          mutate(panel = 'Other NHB')) %>%
  rename(focalX = X1) %>%
  mutate(attract4 = as.numeric(attract*nsd + nmean>=4), 
         attract3 = as.numeric(attract*nsd + nmean>=3)) %>%
  group_by(panel, best_race) %>%
  summarise_at(c('focalX','attract',paste0('attract',2:4)), mean) %>%
  gather(ytype, y, c(attract,attract2,attract4,attract3)) %>%
  mutate(ytype = case_when(ytype=='attract' ~ 'Continuous attractiveness',
                           ytype=='attract2' ~ 'Very attractive',
                           ytype=='attract4' ~ 'Attractive',
                           ytype=='attract3' ~ 'About average'))

nhw_X1 <- wbdot %>%
  filter(best_race=='NHW') %>%
  pull(focalX) %>%
  unique()

pd <- map2(my_attract_plot_data(td %>% filter(partial_black==0) %>% 
                                  rename(focalX = X1), nhw_X1, 'NHB', F),
           my_attract_plot_data(td %>% filter(partial_black==1) %>% 
                                  rename(focalX = X1), nhw_X1, 'NHB', F),
           function(x,y) {
             x %>% 
               mutate(panel = 'Monoracial NHB') %>%
               rbind(y %>% mutate(panel = 'Other NHB'))
           })

pd2 <- map2(my_attract_plot_data(td %>% filter(partial_black==0) %>% 
                                   rename(focalX = X1), nhw_X1, 'NHB', T),
            my_attract_plot_data(td %>% filter(partial_black==1) %>% 
                                   rename(focalX = X1), nhw_X1, 'NHB', T),
            function(x,y) {
              x %>% 
                mutate(panel = 'Monoracial NHB') %>%
                rbind(y %>% mutate(panel = 'Other NHB'))
            }) %>%
  lapply(function(d) {
    d %>% 
      filter(ytype %in% c('Very attractive','Attractive','About average'))
  })

textn <- pd2$pd1 %>% 
  select(aid, panel) %>% 
  unique() %>% 
  group_by(panel) %>% 
  summarise(nlabel = paste('NHB =', n()))

plot_grid(
  attract_plot_legend,
  plot_grid(my_attract_plot(pd, wbdot, F, '', '#ffc300', rug_y_override = 1.012,
                            rug_yend_override = 1.06, yref_override = seq(-0.5,1,0.25)) +
              geom_text(data = textn, aes(x = 0.03, y = Inf, label = nlabel),
                        hjust = 0, vjust = 3.2) +
              facet_wrap(~ panel, scales = 'free'),
            my_attract_plot(pd2, wbdot, T, NULL, '#ffc300') +
              geom_text(data = textn, aes(x = 0.03, y = Inf, label = nlabel),
                        hjust = 0, vjust = 3.2) +
              facet_wrap(~ panel, scales = 'free',
                         strip.position = 'bottom',
                         labeller = as_labeller(c('Monoracial NHB' = 'italic(P)^AFR',
                                                  'Other NHB' = 'italic(P)^AFR'),
                                                default = label_parsed)) +
              theme(strip.placement = 'outside', 
                    strip.text = element_text(face='plain',size=rel(1))),
            nrow = 2, align = 'v', labels = c('A','B')),
  nrow = 2, rel_heights = c(1,14)
)
ggsave('figures/appendix_attract_mono_multirace_black.pdf', width = 12, height = 12)


# 2.9 Supplementary figure: attractiveness plot by Black vs White interviewers
wbdot <- my_attract_plot_ref(fd %>% filter(intrace=='NHB') %>% rename(focalX = X1), 'NHB') %>%
  mutate(panel = 'Black Interviewers') %>%
  rbind(my_attract_plot_ref(fd %>% filter(intrace=='NHW') %>% rename(focalX = X1), 'NHB') %>%
          mutate(panel = 'White Interviewers'))

nhw_X1 <- wbdot %>%
  filter(best_race=='NHW') %>%
  pull(focalX) %>%
  unique()

pd <- map2(my_attract_plot_data(fd %>% filter(intrace=='NHB') %>% 
                                  rename(focalX = X1), nhw_X1, 'NHB', F),
           my_attract_plot_data(fd %>% filter(intrace=='NHW') %>% 
                                  rename(focalX = X1), nhw_X1, 'NHB', F),
           function(x,y) {
             x %>% 
               mutate(panel = 'Black Interviewers') %>%
               rbind(y %>% mutate(panel = 'White Interviewers'))
           })

pd2 <- map2(my_attract_plot_data(fd %>% filter(intrace=='NHB') %>% 
                                   rename(focalX = X1), nhw_X1, 'NHB', T),
            my_attract_plot_data(fd %>% filter(intrace=='NHW') %>% 
                                   rename(focalX = X1), nhw_X1, 'NHB', T),
            function(x,y) {
              x %>% 
                mutate(panel = 'Black Interviewers') %>%
                rbind(y %>% mutate(panel = 'White Interviewers'))
            }) %>%
  lapply(function(d) {
    d %>% 
      filter(ytype %in% c('Very attractive','Attractive','About average'))
  })

textn <- pd2$pd1 %>% 
  select(aid, panel) %>% 
  unique() %>% 
  group_by(panel) %>% 
  summarise(nlabel = paste('NHB =', n()))

plot_grid(
  attract_plot_legend,
  plot_grid(my_attract_plot(pd, wbdot, F, '', '#ffc300', rug_y_override = 1.516,
                             rug_yend_override = 1.58, yref_override = seq(-0.5,1.5,0.25)) +
              geom_text(data = textn, aes(x = 0.03, y = Inf, label = nlabel),
                        hjust = 0, vjust = 3.2) +
              facet_wrap(~ panel, scales = 'free'),
            my_attract_plot(pd2, wbdot, T, NULL, '#ffc300', rug_y_override = 1.0584,
                            rug_yend_override = 1.092) +
              geom_text(data = textn, aes(x = 0.03, y = Inf, label = nlabel),
                        hjust = 0, vjust = 3.2) +
              facet_wrap(~ panel, scales = 'free',
                         strip.position = 'bottom',
                         labeller = as_labeller(c('Black Interviewers' = 'italic(P)^AFR',
                                                  'White Interviewers' = 'italic(P)^AFR'),
                                                default = label_parsed)) +
              theme(strip.placement = 'outside', 
                    strip.text = element_text(face='plain',size=rel(1))),
            nrow = 2, align = 'v', labels = c('A','B')),
  nrow = 2, rel_heights = c(1,14)
)
ggsave('figures/appendix_attract_by_interviewer_race.pdf', width = 12, height = 12)



# 3. Main analyses -------------------------------------------------------------
# 3.1 Main regression
my_fe_dataprep <- function(data, y, covariates) {
  data %>%
    drop_na(all_of(covariates)) %>%
    group_by(wave, intid) %>%
    filter(if (y=='attract') n() > 1
           else length(unique(.data[[y]])) > 1) %>%
    ungroup() %>%
    as.data.frame()
}

my_fe_exp <- function(y, covariates, race_interaction = T) {
  rhs <- if (race_interaction) {
    paste(c('best_race', paste0(covariates, ':best_race')), collapse = ' + ')
  } else {
    paste(covariates, collapse = ' + ')
  }
  as.formula(paste(y, '~', rhs, '-1 | agey + wave^intid'))
}

my_feols <- function(data, y, covariates, race_interaction = T) {
  d <- my_fe_dataprep(data, y, covariates)
  exp <- my_fe_exp(y, covariates, race_interaction)
  reg <- feols(exp, cluster = ~aid, d)
  reg$nmdata <- d
  reg$fe_nid <- d[,'aid'] %>% unique() %>% length()
  return(reg)
}

my_feglm <- function(data, y, covariates, slope_vars = c('X1','X4'), 
                     slope_newdata = NULL, slope_by = 'best_race', race_interaction = T) {
  d <- my_fe_dataprep(data, y, covariates)
  exp <- my_fe_exp(y, covariates, race_interaction)
  reg <- feglm(exp, cluster = ~aid, family = binomial, data = d, data.save = T)
  reg$call$fml <- exp
  reg$call$data <- d
  reg$nmdata <- d
  reg$pseudo_r2 <- r2(reg, 'pr2')
  reg$fe_nid <- d[,'aid'] %>% unique() %>% length()
  slope_data <- if (is.null(slope_newdata)) reg$nmdata else slope_newdata(reg$nmdata)
  slope_args <- list(model = reg, newdata = slope_data, variables = slope_vars, vcov = F)
  if (!is.null(slope_by)) {
    slope_args$by <- slope_by
  }
  reg$me <- do.call(avg_slopes, slope_args)
  return(reg)
}

bc <- c('female', 'X1', 'X3', 'X4')
pc <- c('skin', 'hair', 'eye')
sc <- c('nhood_ses', 'fam_ses')

cov_list <- list(
  m1 = bc,
  m2 = c(bc,'irace_wave'),
  m3 = c(bc,'irace_wave',pc),
  m4 = c(bc,'irace_wave',pc,sc)
)

reg <- list(
  attract = lapply(cov_list, function(cov) {
    my_feols(my_fe_dataprep(fd %>% mutate(across(paste0('X',1:4), ~.*10)),
                            'attract',
                            cov_list[[4]]),
             'attract', cov)
  }),
  attract2 = lapply(cov_list, function(cov) {
    my_feglm(my_fe_dataprep(fd %>% mutate(across(paste0('X',1:4), ~.*10)),
                            'attract2',
                            cov_list[[4]]),
             'attract2', cov)
  })
)

# add AMEs by racial group and derive SEs through bootstrapping
for (i in seq_along(cov_list)) {
  set.seed(218)
  reg$attract2[[i]]$bme <- reg$attract2[[i]]$me %>%
    inferences(method = 'fwb') %>%
    mutate(p.value = 2*(1-pnorm(abs(estimate/std.error))))
}

reg_outcome <- c('Continuous Attractiveness','Very attractive')

rows <- c('Basic Controls', rep('X',4),
          'Racial Classification Controls', rep(' ',1), rep('X',3),
          'Physical Feature Controls', rep(' ',2), rep('X',2),
          'Socioeconomic Controls', rep(' ',3), 'X',
          'Age and Interviewer x Wave FEs', rep('X',4)) %>%
  matrix(nrow = 5, byrow = T) %>%
  as.data.frame()

my_coef_map <- function(var_names, include_x3 = F) {
  coef_omit <- if (include_x3) { 
    'NHB$|Hispanic$|NHW|AINA|AAPI|female|ses|irace|skin|hair|eye'
  } else {
    'NHB$|Hispanic$|NHW|AINA|AAPI|female|X3|ses|irace|skin|hair|eye'
  }
  
  var_names[!grepl(coef_omit, var_names)] %>%
    str_replace('X1', '*P*^AFR^ (10 pp)') %>%
    str_replace('X3', '*P*^EAS^ (10 pp)') %>%
    str_replace('X4', '*P*^IAM^ (10 pp)') %>%
    str_replace('same_race', 'Same race') %>%
    str_remove_all('best_race') %>%
    str_replace_all(':', ' x ')
}

my_main_table <- function(mlist, label, crows, AME = F, include_x3 = F, ame_by_race = T) {
  trows <- lapply(mlist, function(m) {
    d <- data.frame(V1 = c('Num.Inds','Num.Age','Num.Inter x Waves'), 
               V2 = c(m$fe_nid, m$fixef_sizes['agey'], m$fixef_sizes['wave^intid'])) %>%
      rbind(
        m$nmdata %>%
          group_by(best_race) %>%
          summarise(V2 = mean(.data[[m$fml[[2]]]])) %>%
          ungroup %>%
          filter(best_race %in% c('NHW','NHB','Hispanic')) %>%
          transmute(V1 = paste0('Mean outcome: ', best_race), V2 = sprintf('%.3f', V2))
      )
    if (!is.null(m$pseudo_r2)) {
      d <- d %>%
        rbind(data.frame(V1 = 'Pseudo R2', V2 = sprintf('%.3f', m$pseudo_r2)))
    }
    return(d)
  }) %>%
    do.call(cbind, .) %>%
    select(`m1.V1`, ends_with('V2'))
  names(trows) <- names(crows)
  rows <- rbind(crows, trows)
  names(rows) <- c('term', paste0('m',1:length(mlist)))
  rows <- rows %>%
    mutate(part = 'manual', statistic = '')
  
  vars <- lapply(mlist, function(r) {
    names(coef(r))
  }) %>%
    unlist() %>%
    unique() %>%
    my_coef_map(include_x3 = include_x3)
  
  coef_omit <- if (include_x3) {
    'NHB$|Hispanic$|NHW|AINA|AAPI|female|ses|irace|skin|hair|eye'
  } else {
    'NHB$|Hispanic$|NHW|AINA|AAPI|female|X3|ses|irace|skin|hair|eye'
  }
  
  if (AME & mlist[[1]]$fml[[2]]!='attract') {
    tblmlist <- lapply(mlist, function(m) {
      if (ame_by_race) {
        m$bme %>%
          filter(best_race %in% c('NHB','Hispanic')) %>%
          mutate(term = paste0(best_race, term))
      } else {
        m$bme
      }
    })
    rows <- rows %>%
      rbind(lapply(mlist, function(m) m$nobs) %>%
              bind_rows() %>%
              mutate(term = 'Num.Obs.', part = 'manual', statistic = ''))
  } else {
    tblmlist <- mlist
  }
  
  tbl <- modelsummary(tblmlist,
                      estimate = '{estimate}{stars}',
                      coef_omit = coef_omit,
                      coef_rename = vars, 
                      add_rows = rows,
                      gof_omit = 'IC|Std.Errors|FE|RMSE|Adj',
                      output = 'data.frame') %>%
    mutate(across(paste0('m',1:length(mlist)), 
                  ~ case_when(grepl('e',.) & statistic=='estimate' ~ 
                                paste0(sprintf('%.3f', parse_number(.)), gsub('.*[0-9]','',.)),
                              grepl('e',.) & statistic=='std.error' ~ 
                                paste0('(',sprintf('%.3f', parse_number(.)),')'),
                              T ~ .))) %>%
    mutate(part = case_when(grepl('Mean outcome:',term) ~ 'manual2',
                            grepl('Num.Obs',term) ~ 'manual3',
                            grepl('Num.',term) ~ 'manual4',
                            grepl('Pseudo',term) ~ 'manual5',
                            T ~ part)) %>%
    arrange(!(part=='estimates'), !(part=='manual'), !(part=='manual2'), !(part=='manual3'),
            !(part=='manual4'), !grepl('NHB',term), !grepl('Hispanic',term)) %>%
    mutate(term = ifelse(statistic=='std.error','',term)) %>%
    select(-c(part, statistic)) %>%
    flextable() %>%
    set_header_labels(values = c('', paste('Model',1:length(mlist)))) %>%
    add_header_row(values = c('', label), colwidths = c(1,length(mlist)), top = T) %>%
    align(align = 'center', part = 'header') %>%
    hline(i = c(length(vars)*2,length(vars)*2+nrow(crows))) %>%
    colformat_md() %>%
    autofit()
  
  return(tbl)
}

my_combine_tbl <- function(tbl_list, add_var_hline = T) {
  num_models <- tbl_list[[1]]$body$dataset %>%
    select(matches('^m[0-9]+')) %>%
    length()
  num_vars <- tbl_list[[1]]$body$dataset %>%
    filter(grepl('(NHB|Hispanic) ', term)) %>%
    nrow()
  
  tbl <- tbl_list[[1]]$body$dataset %>%
    rename_with(~ paste0('attract_',.), matches('^m[0-9]+')) %>%
    mutate(term = ifelse(term=='', paste('row', row_number()), term)) %>%
    full_join(tbl_list[[2]]$body$dataset %>%
                rename_with(~ paste0('attrac2_',.), matches('^m[0-9]+')) %>%
                mutate(term = ifelse(term=='', paste('row', row_number()), term))) %>%
    mutate(term = ifelse(grepl('row',term), '', term)) %>%
    flextable() %>%
    set_header_labels(values = c('', paste('Model',rep(1:num_models,2)))) %>%
    add_header_row(values = c('', reg_outcome[1:2]), 
                   colwidths = c(1,num_models,num_models), top = T) %>%
    align(align = 'center', part = 'header')
  
  if (add_var_hline) {
    tbl <- tbl %>%
      hline(i = c(num_vars*2, num_vars*2+5))
  }
   
    tbl %>%
    vline(j = num_models+1) %>%
    colformat_md() %>%
    autofit()
}
  
wtbl <- my_combine_tbl(
  list(my_main_table(reg$attract, reg_outcome[[1]], rows, F),
       my_main_table(reg$attract2, reg_outcome[[2]], rows, T))
)

save_as_html(wtbl, path = 'tables/table1.html')


# 3.2 Interaction/decomposition analyses
td <- fd %>%
  mutate(wave = case_when(wave %in% as.character(1:2) ~ '1', T ~ wave))
td <- td %>%
  cbind(predict(dummyVars('~ same_race2 + factor(female) + factor(intfemale) +
                          factor(strataid4) + factor(intage) + wave',
                          fd), newdata = td) %>% as.data.frame()) %>%
  rename_all(~ str_replace(.,'same_race2','samerace2')) %>%
  rename_all(~ str_replace(.,'factor\\(female\\)','female')) %>%
  rename_all(~ str_replace(.,'factor\\(intfemale\\)','intfemale')) %>%
  rename_all(~ str_replace(.,'factor\\(strataid4\\)','strataid4')) %>%
  rename_all(~ str_replace(.,'factor\\(intage\\)','intage')) %>%
  mutate(across(paste0('X',1:4), ~.*10),
         across(matches('(samerace2|female|intfemale|strataid4|intage|wave).+'), 
                ~ .*X1, .names = '{.col}_X1'),
         across(matches('(samerace2|female|intfemale|strataid4|intage|wave).+'), 
                ~ .*X4, .names = '{.col}_X4')
  ) %>%
  select(-ends_with('_X1_X4'))

intvars <- list(m1 = 'samerace2',
                m2 = 'female',
                m3 = 'intfemale',
                m4 = 'strataid4',
                m5 = 'intage',
                m6 = 'wave')
intvar_list <- intvars %>%
  lapply(function(intvar) td %>% select(matches(paste0('^',intvar,'.+'))) %>% names())
intvar_list <- lapply(intvar_list, function(v) v[-1])
intvar_list$m2 <- intvar_list$m2[grepl('_X',intvar_list$m2)]

int <- list(
  attract = lapply(intvar_list, function(cov) {
    my_feols(td, 'attract', c('female','X3',cov))
  }),
  attract2 = lapply(intvar_list, function(cov) {
    my_feglm(td, 'attract2', c('female','X3',cov),
             slope_vars = cov[grepl('_X', cov)],
             slope_newdata = function(d) d %>% filter(best_race %in% c('NHB','Hispanic')))
  })
)

intvar_main_list <- lapply(intvar_list, function(v) {
  c('X1','X4', v[!grepl('X1|X4', v)])
})
intvar_main_list$m2 <- c('female', intvar_main_list$m2)

int_main <- list(
  attract = map2(int$attract, intvar_main_list, function(x,y) {
    my_feols(x$nmdata, 'attract', y)
  }),
  attract2 = map2(int$attract2, intvar_main_list, function(x,y) {
    my_feglm(x$nmdata, 'attract2', y)
  })
)

my_feglm_boot <- function(data, indices, y, covariates, slope_vars = c('X1','X4'),
                          slope_newdata = NULL, slope_by = 'best_race') {
  my_feglm(data[indices,], y, covariates, slope_vars = slope_vars, 
           slope_newdata = slope_newdata, slope_by = slope_by)
}

boot_exports <- c('my_fe_dataprep', 'my_fe_exp', 'my_feglm', 'my_feglm_boot')
boot_cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(boot_cl, {
  library(tidyverse)
  library(fixest)
  library(marginaleffects)
  library(boot)
})
parallel::clusterExport(boot_cl, boot_exports, envir = environment())

for (i in seq_along(intvar_list)) {
  set.seed(218)
  parallel::clusterSetRNGStream(boot_cl, 218)
  
  b <- boot::boot(td,
                  statistic = function(data, indices, covariates, slope_vars) 
                    my_feglm_boot(data, indices, 'attract2', 
                                  covariates, slope_vars = slope_vars,
                                  slope_newdata = function(d) 
                                    d %>% filter(best_race %in% c('NHB','Hispanic')))$me$estimate,
                  R = 1000, parallel = 'snow', cl = boot_cl, 
                  covariates = c('female','X3',intvar_list[[i]]),
                  slope_vars = intvar_list[[i]][grepl('_X', intvar_list[[i]])])
  
  int$attract2[[i]]$bme <- int$attract2[[i]]$me %>%
    mutate(std.error = apply(b$t, 2, sd),
           p.value = 2*(1-pnorm(abs(estimate/std.error))))
}

for (i in seq_along(int_main$attract2)) {
  set.seed(218)
  parallel::clusterSetRNGStream(boot_cl, 218)
  
  b <- boot::boot(int_main$attract2[[i]]$nmdata,
                  statistic = function(data, indices, covariates)
                    my_feglm_boot(data, indices, 'attract2', covariates)$me$estimate,
                  R = 1000, parallel = 'snow', cl = boot_cl, covariates = intvar_main_list[[i]])
  
  int_main$attract2[[i]]$bme <- int_main$attract2[[i]]$me %>%
    mutate(std.error = apply(b$t, 2, sd),
           p.value = 2*(1-pnorm(abs(estimate/std.error))))
}

intvar_levels <- c('samerace2White', 'samerace2Target', 'samerace2Other',
                   'female0', 'female1', 'intfemale0', 'intfemale1',
                   paste0('strataid4',1:4), paste0('intage',0:1),
                   paste0('wave',c(1,3,4)))
intvar_levels <- paste0(rep(c('NHB','Hispanic'),each = length(intvar_levels)),
                        '_', rep(intvar_levels, 2))

subtitle_labels <- list(
  m1 = c(paste(c('White','Black','Other'),'interviewer'), 'P-value',
         paste(c('White','Hispanic','Other'),'interviewer'), 'P-value'),
  m2 = rep(c('Male respondent','Female respondent','P-value'), 2),
  m3 = rep(c('Male interviewer','Female interviewer','P-value'), 2),
  m4 = rep(c('West', 'Midwest', 'South', 'Northeast','P-value'), 2),
  m5 = rep(c('Interviewer under 55', 'Interviewer aged 55+', 'P-value'), 2),
  m6 = rep(c(paste('Age', int$attract$m6$nmdata %>%
             group_by(wave) %>%
             summarise(age = round(mean(age)/12)) %>% pull(age)),'P-value'), 2)
)

my_subgroup_ns <- function(data, int_prefix) {
  data %>%
    group_by(best_race, grp_var) %>%
    summarise(nobs = n(),
              nind = length(unique(aid)),
              nint = length(unique(intid))) %>%
    ungroup() %>%
    filter(best_race %in% c('NHB','Hispanic')) %>%
    gather(term, n, nobs:nint) %>%
    transmute(term, n, group = paste0(best_race,'_',int_prefix,grp_var),
              group = factor(group, levels = intvar_levels)) %>%
    spread(group, n)
}

my_test_intcoef <- function(mlist, intvarname_list, yvar, reg_label, race_label, anc_level) {
  tlevels <- unique(na.omit(td[[intvarname_list[[reg_label]]]]))
  flevels <- paste0('best_race', race_label, ':', intvarname_list[[reg_label]],
                    tlevels, '_', anc_level)
  exp <- paste(flevels[1], '=', flevels[-1])
  linearHypothesis(mlist[[yvar]][[reg_label]], exp) %>% 
    drop_na(`Pr(>Chisq)`) %>% 
    pull(`Pr(>Chisq)`)
}

my_int_table <- function(mlist, subtitle_list, intvarname_list, yvar, reg_label) {
  ndata <- mlist[[yvar]][[reg_label]]$nmdata
  ndata$grp_var <- ndata[[intvarname_list[[reg_label]]]]
  add_n <- my_subgroup_ns(ndata, intvarname_list[[reg_label]])
  m <- if (yvar=='attract2') {
    mlist[[yvar]][[reg_label]]$bme %>%
      filter(best_race %in% c('NHB','Hispanic')) %>%
      mutate(term = paste0(best_race, ' × ', term))
  } else {
    mlist[[yvar]][[reg_label]]
  }
  
  int_pvar <- lapply(c('attract','attract2'), function(y) {
    lapply(paste0('m',1:6), function(m) {
      lapply(c('NHB','Hispanic'), function(r) {
        lapply(c('X1','X4'), function(a) {
          data.frame(outcome = y, model = m, best_race = r, ancestry = a,
                     pvalue = my_test_intcoef(mlist, intvarname_list, y, m, r, a))
        }) %>%
          bind_rows()
      }) %>%
        bind_rows()
    })
  }) %>%
    bind_rows() %>%
    group_by(best_race, ancestry) %>%
    mutate(pvalue_adj = p.adjust(pvalue, method = 'fdr')) %>%
    ungroup()
  
  add_p <- int_pvar %>%
    filter(outcome==yvar, model==reg_label) %>%
    transmute(term = ancestry, best_race,
              p = paste0(sprintf('%.3f',pvalue), ' (', sprintf('%.3f',pvalue_adj),')')) %>%
    spread(best_race, p) %>%
    select(term, p_NHB = NHB, p_Hispanic = Hispanic)
  
  modelsummary(m,
               estimate = '{estimate}{stars}',
               coef_omit = '^(?!.*(NHB|Hispanic).*(X1|X4))',
               gof_omit = 'IC|Std.Errors|FE|RMSE|Adj|R|Num',
               output = 'data.frame') %>%
    mutate(term = str_remove_all(term, 'best_race')) %>%
    separate(term, c('group', 'term'), '_') %>%
    mutate(group = str_replace_all(group, ' × ', '_'),
           group = factor(group, levels = intvar_levels)) %>%
    spread(group, `(1)`) %>%
    rbind(add_n %>% mutate(part = '', statistic = '')) %>%
    left_join(add_p %>% mutate(statistic = 'estimate')) %>%
    mutate(term = case_when(term=='X1' ~ '*P*^AFR^ (10 pp)',
                            term=='X4' ~ '*P*^IAM^ (10 pp)',
                            term=='nobs' ~ 'N Observations',
                            term=='nind' ~ 'N Respondents',
                            term=='nint' ~ 'N Interviewers'),
           term = ifelse(statistic=='std.error','',term)) %>%
    select(-c(part, statistic)) %>%
    relocate(p_NHB, .before = matches('Hispanic')) %>%
    flextable() %>%
    set_header_labels(values = c('',subtitle_list[[reg_label]])) %>%
    add_header_row(values = c('','NHB','Hispanic'), 
                   colwidths = c(1,length(subtitle_list[[reg_label]])/2,
                                 length(subtitle_list[[reg_label]])/2)) %>%
    add_header_row(values = c('', ifelse(yvar=='attract',
                                         'Attractiveness Rating (Continuous)',
                                         'Pr("Very Attractive")')), 
                   colwidths = c(1,length(subtitle_list[[reg_label]]))) %>%
    align(align = 'center', part = 'header') %>%
    vline(j = 1+length(subtitle_list[[reg_label]])/2) %>%
    colformat_md() %>%
    autofit()
}

do.call(save_as_html, 
        c(
          flatten(lapply(paste0('m',1:6), function(r) {
            lapply(c('attract','attract2'), function(y) {
              my_int_table(int, subtitle_labels, intvars, y, r)
            })
          })),
          list(path = 'tables/appendix_decompose_table.html')
        ))

# create panel A for Figure 3: decomposition of overall effects
pd <- lapply(int_main, function(y) {
  lapply(y, function(m) {
    if (is.null(m$bme)) {
      td <- coeftable(m) %>%
        as.data.frame() %>%
        rownames_to_column('term') %>%
        rename_all(tolower) %>%
        rename(std.error = `std. error`, p.value = `pr(>|t|)`)
    } else {
      td <- m$bme %>%
        mutate(term = paste0(best_race,':',term)) %>%
        rename_all(tolower)
    }
    td %>%
      filter(grepl('X1', term), grepl('NHB', term)) %>%
      select(estimate, std.error, p.value)
  }) %>%
    bind_rows(.id = 'model') %>%
    mutate(subgroup = model, term = 'X1')
})

pd <- map2(int, pd, function(x, y) {
  lapply(x, function(m) {
    if (is.null(m$bme)) {
      td <- coeftable(m) %>%
        as.data.frame() %>%
        rownames_to_column('term') %>%
        rename_all(tolower) %>%
        rename(std.error = `std. error`, p.value = `pr(>|t|)`)
    } else {
      td <- m$bme %>%
        mutate(term = paste0(best_race,':',term)) %>%
        rename_all(tolower)
    }
    td %>%
      filter(grepl('X1', term), grepl('NHB', term)) %>%
      select(term, estimate, std.error, p.value)
  }) %>%
    bind_rows(.id = 'model') %>%
    separate(term, c('best_race','term'), ':') %>%
    mutate(subgroup = str_remove_all(term, '_X1')) %>%
    select(-best_race) %>%
    rbind(y) %>%
    add_row(model = paste0('m',1:6), subgroup = paste0('fill',1:6), estimate = 0) %>%
    add_row(model = 'm1', subgroup = 'fill1a', estimate = 0) %>%
    mutate(upper = estimate + 1.96*std.error, lower = estimate - 1.96*std.error)
})

add_p <- lapply(c('attract','attract2'), function(y) {
  lapply(paste0('m',1:6), function(m) {
    data.frame(outcome = y, model = m,
               pvalue = my_test_intcoef(int, intvars, y, m, 'NHB', 'X1'))
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  mutate(pvalue_adj = p.adjust(pvalue, method = 'fdr')) %>%
  ungroup() %>%
  mutate(pvalue = case_when(pvalue<0.001 ~ '***',
                            pvalue<0.01 ~ '**',
                            pvalue<0.05 ~ '*',
                            T ~ 'NS'),
         adj_sig = ifelse(pvalue_adj < 0.05, 1, 0)) %>%
  group_by(outcome) %>%
  group_split() %>%
  map(~{
    .x <- .x %>%
      left_join(data.frame(model = paste0('m',1:6),
                           x = c(3.5,8,12,17,22,27)))
  }) %>%
  set_names(c('attract','attract2'))

sglevel <- c('fill1a', 'm1', 'samerace2White', 'samerace2Target', 'samerace2Other',
             'm2', 'female0', 'female1',
             'm3', 'intfemale0', 'intfemale1',
             'm4', paste0('strataid4',1:4),
             'm5', paste0('intage',0:1),
             'm6', paste0('wave',c(1,3,4)))

sglabel <- c('fill1a', 'm1', 'White', 'Black', 'Other',
             'm2', 'Male', 'Female',
             'm3', 'Male ', 'Female ',
             'm4', 'West', 'Midwest', 'South', 'Northeast',
             'm5', 'Under 55', 'Age 55+',
             'm6', paste('Age', int$attract$m6$nmdata %>%
                            group_by(wave) %>%
                            summarise(age = round(mean(age)/12)) %>% pull(age)))

add_p2 <- lapply(add_p, function(d) {
  d %>%
    mutate(x = case_when(model=='m1' ~ 2.5,
                         model=='m2' ~ 2,
                         model=='m3' ~ 2,
                         model=='m4' ~ 3,
                         model=='m5' ~ 2,
                         model=='m6' ~ 2.5))
})

fig3a <- pd %>%
  bind_rows(.id = 'panel') %>%
  filter(!is.na(term) | subgroup=='fill1a') %>%
  mutate(across(c(estimate,upper,lower), ~ ifelse(panel=='attract2',-3*.,.)),
         subgroup = factor(subgroup, levels = rev(sglevel), labels = rev(sglabel))) %>%
  ggplot(aes(x = subgroup, y = estimate)) +
  geom_hline(yintercept = setdiff(seq(-0.18,0.18,0.06),0), linetype = 'dashed', color = 'gray') +
  geom_bar(aes(fill = subgroup), stat = 'identity', position = 'dodge', 
           color = 'black', show.legend = F) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = panel), width = 0.2, 
                show.legend = F) +
  geom_hline(yintercept = 0) +
  geom_text(data = add_p2 %>%
              bind_rows(.id = 'panel') %>%
              mutate(y = ifelse(panel=='attract', -0.24, 0.24)),
            aes(x = x, y = y, label = pvalue, color = factor(adj_sig)), show.legend = F) +
  geom_text(data = data.frame(y = c(-0.12,0.12), model = 'm1',
                              label = c('Attractiveness score','Very attractive')),
            aes(x = Inf, y = y, label = label), vjust = 1, fontface = 'bold') +
  facet_grid(model ~ ., scales = 'free_y', space = 'free', switch = 'y',
             labeller = as_labeller(c('m1' = 'Interviewer race',
                                      'm2' = 'Respondent gender',
                                      'm3' = 'Interviewer gender',
                                      'm4' = 'Region',
                                      'm5' = 'Interviewer age',
                                      'm6' = 'Respondent age'))) +
  coord_flip() +
  scale_y_continuous('Coefficients                                 AMEs    ',
                     breaks = seq(-0.24,0.24,0.06), 
                     labels = c(seq(-0.24,0,0.06), (seq(0,-0.24,-0.06)/3)[-1])) +
  scale_x_discrete('', labels = setdiff(rev(sglabel), 'fill1a') %>%
                     str_replace_all('m[1-6]', 'Overall'),
                   breaks = setdiff(rev(sglabel), 'fill1a')) +
  scale_fill_manual('', values = rev(c('white', 'white', '#ffe699', '#ffa726', '#f57c00',
                                       'white', '#ffe699', '#ffa726',
                                       'white', '#ffe699', '#ffa726',
                                       'white', '#ffe699', '#ffa726', '#f57c00', '#E65100',
                                       'white', '#ffe699', '#ffa726',
                                       'white', '#ffe699', '#ffa726', '#f57c00'))) +
  scale_color_manual('', values = c('black', 'red', 'black', 'blue')) +
  theme_custom() +
  theme(strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, vjust = 1))


# 3.3 GSPs complement skin tone & other physical measures in understanding attractiveness among NHB
td <- s4 %>%
  filter(best_race=='NHB') %>%
  mutate(skin2 = 6 - as.numeric(skin)) %>%
  drop_na(skin, X1)

tdis <- td %>%
  arrange(skin) %>%
  mutate(n = n() + 1 - row_number()) %>%
  select(pseudo_skin = skin, n) %>%
  arrange(n)

skin_prop <- td %>% 
  group_by(skin) %>%
  summarise(n = n()) %>%
  ungroup %>%
  mutate(prop = n/sum(n),
         label = paste0(skin, '\n(', sprintf('%.2f',prop*100), '%)'))

pd <- td %>%
  arrange(X1) %>%
  mutate(n = row_number()) %>%
  left_join(tdis) %>%
  pivot_longer(cols = matches('^attract'),
               names_to = c('.value', 'wave'),
               names_pattern = '^(attract)(\\d+)$',
               values_drop_na = T) %>%
  mutate(attract2 = as.numeric(attract==5), attract = (attract - nmean)/nsd)

# integrate information using the machine learning approach
mld <- fd %>% 
  cbind(predict(dummyVars('~ hair + eye', fd, sep = '_'), newdata = fd) %>%
          as.data.frame()) %>%
  rename_with(~str_replace(.,' ','_'), matches('^(hair|eye)_')) %>%
  filter(best_race=='NHB') %>% 
  drop_na(attract, skin, hair, eye, irace) %>%
  mutate(ib = as.numeric(irace=='Black'), attract = attract*nsd + nmean,
         skin = 6 - as.numeric(skin),
         hair_Other = hair_Blond + hair_Red + hair_Grey + hair_Other) %>%
  select(aid, X1, skin, matches('^(hair|eye)_'), ib, attract, attract2, 
         -c(hair_Blond, hair_Red, hair_Grey, hair_No_hair, eye_Hazel, eye_Blue, 
            eye_Green, eye_Other)) %>%
  # collapse to the individual-level
  group_by(aid) %>%
  summarise_all(mean) %>%
  ungroup()

set.seed(218)
train_index <- createDataPartition(mld$attract, p = 0.8, list = F)

# function to find optimal alpha using cross-validation
my_cv_glmnet <- function(data, yvar, alpha_seq = seq(0, 1, 0.01)) {
  set.seed(218)
  cv_foldid <- createFolds(data[train_index,][[yvar]]) %>%
    lapply(function(v) data.frame(index = v)) %>%
    bind_rows(.id = 'fold') %>%
    arrange(index) %>%
    mutate(fold = parse_number(fold)) %>%
    pull(fold)
  cv_errors <- list()
  cv_models <- list()
  for (i in seq_along(alpha_seq)) {
    # cross-validation to evaluate this alpha
    cv_models[[i]] <- glmnet::cv.glmnet(
      data[train_index,] %>% 
        select(-c(aid, attract, attract2)) %>%
        as.matrix(),
      data[train_index,yvar] %>% as.matrix(),
      foldid = cv_foldid,
      alpha = alpha_seq[i]
    )
    # store the mean cross-validation error
    cv_errors[[i]] <- data.frame(alpha = alpha_seq[[i]],
                               lambda = cv_models[[i]]$lambda,
                               cvm = cv_models[[i]]$cvm)
  }
  cv_errors <- bind_rows(cv_errors)
  return(list(model = cv_models, mean_cvm = cv_errors,
              best_tune = cv_errors %>% filter(cvm==min(cvm))))
}

my_find_optalpha <- function(data, yvar, alpha_seq = seq(0, 1, 0.01)) {
  round1 <- my_cv_glmnet(data, yvar)
  
  fitControl <- trainControl(method = 'repeatedcv',
                             number = 10,
                             repeats = 10)
  set.seed(218)
  round2 <- train(as.formula(paste(yvar, '~ .')),
                  data = data[train_index,] %>% 
                    select(-aid, -setdiff(c('attract','attract2'), yvar)),
                  method = 'glmnet',
                  preProc = c('center', 'scale'),
                  tuneLength = 25,
                  trControl = fitControl)
  
  fitGrid <-  round2$results %>%
    select(alpha, lambda) %>%
    rbind(round1$best_tune %>%
            select(alpha, lambda)) %>%
    rename(.alpha = alpha, .lambda = lambda)
  set.seed(218)
  round3 <- train(as.formula(paste(yvar, '~ .')),
                  data = data[train_index,] %>% 
                    select(-aid, -setdiff(c('attract','attract2'), yvar)),
                  method = 'glmnet',
                  preProc = c('center', 'scale'),
                  tuneGrid = fitGrid,
                  trControl = fitControl)
  
  return(round3)
}

mlm <- lapply(c('attract','attract2'), function(y) {
  set.seed(218)
  my_find_optalpha(mld, y)
}) %>%
  set_names(c('attract','attract2'))

glmnet_predict <- function(object, newdata) {
  # use the optimal lambda that minimizes cross-validation error
  if (grepl('train', class(object))[[1]]) {
    predict(object, newdata = newdata)
  } else {
    predict(object, newx = newdata, s = 'lambda.min')[,1]
  }
}

new_index <- round(skin_prop$prop*nrow(mld[-train_index,]))
tdis2 <- lapply(seq_along(new_index), function(i) {
  data.frame(pseudo_skin = rep(skin_prop$skin[[i]], new_index[[i]]))
}) %>%
  bind_rows() %>%
  mutate(n = rev(row_number())) %>%
  arrange(n)

my_attract_pseudo_skin <- function(y_var, y_label = NULL) {
  pdat <- pd %>%
    drop_na(attract) %>%
    gather(type, skin_grp, c(skin,pseudo_skin)) %>%
    group_by(skin_grp, type)%>%
    summarise(attract = mean(attract), attract2 = mean(attract2)) %>%
    ungroup() %>%
    mutate(skin_grp = factor(skin_grp, levels = skin_prop$skin),
           type = factor(type, levels = c('skin','pseudo_skin'))) %>%
    filter(skin_grp!='White')
  pdat$y <- pdat[[y_var]]
  
  p <- pdat %>%
    ggplot(aes(x = type, y = y, group = skin_grp, fill = type)) +
    geom_bar(stat = 'identity', position = 'dodge', alpha = 0.9, color = 'black') +
    scale_x_discrete('', labels = c('Observed', expression(italic(P)^AFR)))
  
  if (y_var=='attract') {
    p <- p +
      geom_hline(yintercept = nmean, linetype = 'dashed', color = 'gray40') +
      coord_cartesian(ylim = c(-0.25,0.35)) +
      scale_y_continuous('Average attractiveness score\n(Standardized)', expand = c(0,0), 
                         breaks = seq(-0.25,0.35,0.1),
                         sec.axis = sec_axis(~ .*nsd + nmean, '(Unstandardized)'))
  } else {
    p <- p +
      scale_y_continuous(y_label, expand = c(0,0))
  }
   p <- p +
    scale_fill_manual(values = c('#8ecae6','#ffc300','#d6ccc2'),
                      labels = c('Observed', expression(italic(P)^AFR), 'Composite')) + 
    labs(fill = '') +
    theme_custom() +
    theme(legend.position = 'none', strip.placement = 'outside')
   return(p)
}

my_attract_skin_byvar <- function(y_var, y_label = NULL) {
  pdat <- pd %>%
    drop_na(attract) %>%
    group_by(skin) %>%
    arrange(skin, X1) %>%
    mutate(higher = row_number() > n()/2) %>%
    group_by(skin, higher) %>%
    summarise(across(c(attract,attract2), mean)) %>%
    ungroup() %>%
    filter(skin!='White')
  
  pdat$y <- pdat[[y_var]]
  
  p <- pdat %>%
    ggplot(aes(x = skin, y = y, alpha = higher)) +
    geom_bar(stat = 'identity', position = 'dodge', fill = '#ffc300') +
    scale_x_discrete('', labels = skin_prop$label)
  
  if (y_var=='attract') {
    p <- p +
      geom_hline(yintercept = nmean, linetype = 'dashed', color = 'gray40') +
      coord_cartesian(ylim = c(-0.25,0.35)) +
      scale_y_continuous('Average attractiveness score\n(Standardized)', expand = c(0,0), 
                         breaks = seq(-0.25,0.35,0.1),
                         sec.axis = sec_axis(~ .*nsd + nmean, '(Unstandardized)'))
  } else {
    p <- p +
      scale_y_continuous(y_label, expand = c(0,0))
  }
  p <- p +
    scale_alpha_manual(values = c(0.4,0.9), labels = c('< median','> median')) +
    labs(y = '', alpha = expression(italic(P)^AFR)) +
    theme_custom() +
    theme(legend.position = 'none', strip.placement = 'outside')
  return(p)
}

my_ml_shap <- function(y_var, flip_coord) {
  set.seed(218)
  ex <- fastshap::explain(mlm[[y_var]], X = mld[train_index,] %>% 
                            select(-c(aid,attract,attract2)) %>% as.matrix(), 
                          pred_wrapper = glmnet_predict,
                          newdata = mld[-train_index,] %>% 
                            select(-c(aid,attract,attract2)) %>% as.matrix(),
                          nsim = 100, adjust = T, shap_only = F)
  
  p <- shapviz::shapviz(ex) %>%
    shapviz::sv_importance() +
    scale_y_discrete(labels = function(x) parse(text = case_when(
      x=='X1' ~ 'italic(P)^AFR',
      x=='skin' ~ 'Skin~tone',
      x=='ib' ~ 'Classified~Black',
      T ~ str_replace(sub('^(hair|eye)_(\\w+)$', '\\2~\\1', x), '_', '~')
    ))) +
    labs(subtitle = ifelse(y_var=='attract', 'Attractiveness score', 'Very attractive')) +
    theme_custom() +
    theme(plot.subtitle = element_text(hjust = 0.5))
  
  if (flip_coord) {
    p <- p +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  return(p)
}

plot_grid(
  plot_grid(my_attract_skin_byvar('attract')+theme(legend.position='top'),
            my_attract_pseudo_skin('attract'),
            align = 'v', rel_heights = c(1.2,1), nrow = 2, labels = c('A','C')), # rel_heights = c(1.2,1,1), nrow = 3
  plot_grid(my_attract_skin_byvar('attract2','Fraction "very attractive"')+
              theme(legend.position='top'),
            my_attract_pseudo_skin('attract2','Fraction "very attractive"'),
            align = 'v', rel_heights = c(1.2,1), nrow = 2, labels = c('B','D')),
  nrow = 1, rel_widths = c(1.05,1)
)
ggsave('figures/appendix_attract_by_skinafr.pdf', height = 9, width = 10)

# create panel B for Figure 3: SHAP values from machine learning models
fig3b <- lapply(c('attract','attract2'), function(y) my_ml_shap(y, F)) %>%
  plot_grid(plotlist = ., nrow = 2)

plot_grid(fig3a, fig3b, ncol = 2, rel_widths = c(3,2), labels = c('A','B'))
ggsave('figures/figure3.pdf', width = 12, height = 8)



# 4. Supplementary analyses ----------------------------------------------------
# 4.1 GSP - attractiveness association without interactions with race
m <- list(
  attract = lapply(cov_list, function(cov) {
    my_feols(fd %>% mutate(across(paste0('X',1:4), ~.*10)),
             'attract', cov, race_interaction = F)
  }),
  attract2 = lapply(cov_list, function(cov) {
    my_feglm(fd %>% mutate(across(paste0('X',1:4), ~.*10)), 'attract2', cov,
             race_interaction = F, slope_vars = c('X1','X3','X4'), slope_by = NULL)
  })
)

for (i in seq_along(cov_list)) {
  set.seed(218)
  m$attract2[[i]]$bme <- m$attract2[[i]]$me %>%
    inferences(method = 'fwb') %>%
    mutate(p.value = 2*(1-pnorm(abs(estimate/std.error))))
}

my_combine_tbl(list(
  my_main_table(m$attract, reg_outcome[[1]], rows, F, T, F),
  my_main_table(m$attract2, reg_outcome[[2]], rows, T, T, F)
), add_var_hline = F) %>%
  save_as_html(path = 'tables/appendix_attract_gsp.html')


# 4.2 Models with race concordance * self-reported race terms
m <- list(
  attract = lapply(cov_list, function(cov) {
    my_feols(fd %>% mutate(across(paste0('X',1:4), ~.*10)), 'attract', c(cov,'same_race'))
  }),
  attract2 = lapply(cov_list, function(cov) {
    my_feglm(fd %>% mutate(across(paste0('X',1:4), ~.*10)), 'attract2', c(cov,'same_race'),
             slope_vars = c('X1','X4','same_race'), slope_newdata = function(d) 
               d %>% filter(best_race %in% c('NHW','NHB','Hispanic')))
  })
)

for (i in seq_along(cov_list)) {
  set.seed(218)
  parallel::clusterSetRNGStream(boot_cl, 218)
  
  covariates_i <- c(cov_list[[i]],'same_race')
  
  b <- boot::boot(fd %>% mutate(across(paste0('X',1:4), ~.*10)),
                  statistic = function(data, indices, covariates) 
                    my_feglm_boot(data, indices, 'attract2', covariates,
                                  slope_vars = c('X1', 'X4', 'same_race'),
                                  slope_newdata = function(d) 
                                    d %>% filter(best_race %in% c('NHW','NHB','Hispanic')))$me$estimate,
                  R = 1000, parallel = 'snow', cl = boot_cl, 
                  covariates = c(cov_list[[i]], 'same_race'))
  
  m$attract2[[i]]$bme <- m$attract2[[i]]$me %>%
    mutate(std.error = apply(b$t, 2, sd),
           p.value = 2*(1-pnorm(abs(estimate/std.error))))
}
parallel::stopCluster(boot_cl)

my_combine_tbl(
  list(my_main_table(m$attract, reg_outcome[[1]], rows, F),
       my_main_table(m$attract2, reg_outcome[[2]], rows, T))
) %>%
  save_as_html(path = 'tables/appendix_attract_racebyconcordance.html')

