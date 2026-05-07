# 0. Init ----------------------------------------------------------------------
library(here); library(haven); library(tidyverse)

setwd(str_remove(here(), '/scripts'))



# 1. Import data ---------------------------------------------------------------
# 1.1 survey data
ah <- list()
ah$w1 <- read_xpt('allwave1.xpt') %>% rename_all(tolower)
for (i in 2:5) {
  ah[[paste0('w',i)]] <- read_xpt(paste0('wave', i, '.xpt')) %>% rename_all(tolower)
}
ah$wt4 <- read_xpt('weights4.xpt') %>% rename_all(tolower)
ah$int4 <- read_xpt('intid4.xpt') %>% rename_all(tolower)
ah$ses <- read_xpt('conses.xpt') %>% rename_all(tolower)


# 1.2 linkage data
link <- read_xpt('GID_link.xpt') %>% rename_all(tolower)



# 2. Prepare data --------------------------------------------------------------
# key variables from core survey data
ppheno <- ah$w4 %>%
  transmute(
    aid,
    # demographics
    sex = bio_sex4 %>% factor(labels = c('Male', 'Female')),
    byear = h4od1y, bmonth = h4od1m, age4 = (iyear4-byear)*12 + imonth4-bmonth,
    edu = case_when(
      h4ed2 %in% c(1,2) ~ 'Less than high school',
      h4ed2==3 ~ 'High school',
      h4ed2 %in% 4:6 ~ 'Some college',
      h4ed2>=7 & h4ed2<=13 ~ "Bachelor's degree or above",
      T ~ NA
    ) %>% factor(levels = c('Less than high school', 'High school', 'Some college',
                            "Bachelor's degree or above")),
    # interviewer-reported attractiveness and race
    h4ir1, h4ir4,
    # interview time
    iyear4, imonth4, iday4
  ) %>%
  full_join(
    ah$w3 %>% 
      mutate(age3 = (iyear3-h3od1y)*12 + imonth3-h3od1m,
             across(c('h3od2', starts_with('h3od4')), ~ ifelse(.>1, NA, .)),
             # physical features
             skin = h3ir17 %>% factor(labels = c('Black', 'Dark brown', 'Medium brown',
                                                 'Light brown', 'White')),
             hair = h3ir18 %>% factor(labels = c('No hair', 'Black', 'Brown', 'Blond',
                                                 'Red', 'Grey', 'Other')),
             eye = h3ir19 %>% factor(labels = c('Black', 'Brown', 'Hazel', 'Blue',
                                                'Green', 'Other')),
             # W3 interviewer's info
             i3female = as.numeric(figender==2),
             i3racec = case_when(fihisp==1 ~ 'Hispanic', 
                                 fihisp==0 & firace==1 ~ 'NHW',
                                 fihisp==0 & firace==2 ~ 'NHB',
                                 # fihisp==0 & !is.na(firace) ~ 'Other',
                                 fihisp==0 & firace==3 ~ 'AAPI',
                                 fihisp==0 & firace==4 ~ 'AINA',
                                 fihisp==0 & firace==5 ~ 'Other',
                                 T ~ NA),
             i3byear = ifelse(fibyear==9999, NA, fibyear),
             i3edu = case_when(
               fiedu==1 ~ 'High school or less',
               fiedu==2 ~ 'Some college',
               fiedu %in% 3:4 ~ "Bachelor's degree or above",
               T ~ NA
             ) %>% factor(levels = c('High school or less', 'Some college', 
                                     "Bachelor's degree or above")),
             i3spksp = ifelse(fispksp>1, NA, fispksp))
    ) %>%
  full_join(ah$w1) %>%
  full_join(ah$w2) %>%
  mutate(intid1 = as.numeric(intid), intid2 = as.numeric(intid2), 
         across(matches('h[12]gi1[ym]'), ~ ifelse(.>=96, NA, .)),
         age1 = (iyear-h1gi1y)*12 + imonth-h1gi1m,
         age2 = (iyear2-h2gi1y)*12 + imonth2-h2gi1m,
         # interviewer-rated attractiveness
         across(matches('h[1-4]ir1$'), ~ ifelse(.>5, NA, .)),
         # interviewer-classified race
         irace1 = ifelse(h1gi9>=6, NA, h1gi9),
         irace1 = factor(irace1, labels = c('White', 'Black', 'Native American',
                                            'Asian & Pacific Islander', 'Other')),
         irace3 = factor(h3ir4, labels = c('White', 'Black', 'Native American',
                                           'Asian & Pacific Islander')),
         irace4 = factor(h4ir4, labels = c('White', 'Black', 'Native American',
                                           'Asian & Pacific Islander')),
         # interview time
         iyear1 = iyear + 1900, iyear2 = iyear2 + 1900, imonth1 = imonth, iday1 = iday) %>%
  # survey design: cluster variable - school ID; strata variable - region
  full_join(ah$wt4 %>% select(aid, wt4 = gswgt4_2, lwt4 = gswgt4,
                              clusterid4 = psuscid, strataid4 = region)) %>%
  full_join(ah$int4 %>% 
              rename_all(tolower) %>%
              # W4 interviewer's info
              transmute(aid, intid4, i4female = case_when(i4gender==2 ~ 1,
                                                          i4gender==6 ~ NA, T ~ 0),
                        i4race = ifelse(i4race>5, NA, i4race), 
                        i4hisp = ifelse(i4hisp>1, NA, i4hisp),
                        i4racec = case_when(i4hisp==1 ~ 'Hispanic', 
                                            i4hisp==0 & i4race==1 ~ 'NHW',
                                            i4hisp==0 & i4race==2 ~ 'NHB',
                                            i4hisp==0 & i4race==3 ~ 'AAPI',
                                            i4hisp==0 & i4race==4 ~ 'AINA',
                                            i4hisp==0 & i4race==5 ~ 'Other',
                                            T ~ NA),
                        i4byear = ifelse(i4byear==9996, NA, i4byear),
                        i4edu = case_when(
                          i4edu %in% 1:2 ~ 'High school or less',
                          i4edu %in% 3:4 ~ 'Some college',
                          i4edu %in% 5:7 ~ "Bachelor's degree or above",
                          T ~ NA
                        ) %>% factor(levels = c('High school or less', 'Some college',
                                                "Bachelor's degree or above")),
                        i4spksp = ifelse(i4spksp>1, NA, i4spksp),
                        i4ahexp = case_when(i4ahexp %in% 1:2 ~ 1, i4ahexp==6 ~ NA, T ~ i4ahexp))) %>%
  full_join(ah$w5 %>% select(aid, matches('h5od4[abcdefg]'), h5od8)) %>%
  mutate(
    w3nhw = factor(h3od4a, labels = c('NA','NHW')),
    w3nhb = factor(h3od4b, labels = c('NA','NHB')),
    w3aina = factor(h3od4c, labels = c('NA','AINA')),
    w3aapi = factor(h3od4d, labels = c('NA','AAPI')),
    w3all = paste(w3nhw, w3nhb, w3aina, w3aapi, sep = ', '),
    w3all = case_when(w3all=='NA, NA, NA, NA' ~ NA, 
                      T ~ str_remove_all(w3all, '^(NA, )*|(, NA)*')),
    h3od6 = ifelse(h3od6<=4, h3od6, NA),
    h3od6 = factor(h3od6, labels = c('NHW', 'NHB', 'AINA', 'AAPI')),
    w5nhw = factor(h5od4a==1, labels = c('NA','NHW')),
    w5nhb = factor(h5od4b==1, labels = c('NA','NHB')),
    w5aina = factor(h5od4f==1, labels = c('NA','AINA')),
    w5aapi = factor(h5od4d==1|h5od4e==1, labels = c('NA','AAPI')),
    w5all = paste(w5nhw, w5nhb, w5aina, w5aapi, sep = ', '),
    w5all = case_when(w5all=='NA, NA, NA, NA' ~ NA, T ~ str_remove_all(w5all, '^(NA, )*|(, NA)*')),
    h5od8 = case_when(
      h5od8==1 ~ 'NHW',
      h5od8==2 ~ 'NHB',
      h5od8==21 ~ 'AINA',
      h5od8 %in% c(9:15,17:18,20) ~ 'AAPI',
      h5od8 %in% c(3:7) ~ 'Hispanic',
      h5od8==22 ~ NA
    ),
    w1nhw = factor(h1gi6a==1, labels = c('NA','NHW')),
    w1nhb = factor(h1gi6b==1, labels = c('NA','NHB')),
    w1aina = factor(h1gi6c==1, labels = c('NA','AINA')),
    w1aapi = factor(h1gi6d==1, labels = c('NA','AAPI')),
    w1all = paste(w1nhw, w1nhb, w1aina, w1aapi, sep = ', '),
    w1all = case_when(w1all=='NA, NA, NA, NA' ~ NA, 
                      T ~ str_remove_all(w1all, '^(NA, )*|(, NA)*')),
    h1gi8 = ifelse(h1gi8<=4, h1gi8, NA),
    h1gi8 = factor(h1gi8, labels = c('NHW', 'NHB', 'AINA', 'AAPI')),
    # best race
    best_race2 = case_when(
      h3od2==1 ~ 'Hispanic',
      h3od2==0 & !is.na(w3all) & !grepl(',', w3all) ~ w3all,
      h3od2==0 & !is.na(w3all) & grepl(',', w3all) & !is.na(h3od6) ~ h3od6,
      T ~ NA
    ),
    best_race = case_when(
      h3od2==0 & !is.na(w3all) & str_count(w3all, ',')==1 & grepl('NHW', w3all) ~ 
        str_remove(w3all, '^(NHW, )*|(, NHW)*'),
      T ~ best_race2
    ),
    best_race2 = case_when(
      is.na(best_race2) & h5od4c==1 ~ 'Hispanic',
      is.na(best_race2) & h5od4c==0 & !is.na(w5all) & !grepl(',', w5all) ~ w5all,
      is.na(best_race2) & h5od4c==0 & !is.na(w5all) & grepl(',', w5all) & !is.na(h5od8) ~ h5od8,
      T ~ best_race2
    ),
    best_race = case_when(
      is.na(best_race) & h5od4c==0 & !is.na(w5all) & str_count(w5all, ',')==1 & grepl('NHW', w5all) ~ 
        str_remove(w5all, '^(NHW, )*|(, NHW)*'),
      T ~ best_race2
    ),
    best_race2 = case_when(
      is.na(best_race2) & h1gi4==1 ~ 'Hispanic',
      is.na(best_race2) & h1gi4==0 & !is.na(w1all) & !grepl(',', w1all) ~ w1all,
      is.na(best_race2) & h1gi4==0 & !is.na(w1all) & grepl(',', w1all) & !is.na(h1gi8) ~ h1gi8,
      T ~ best_race2
    ),
    best_race = case_when(
      is.na(best_race) & h1gi4==0 & !is.na(w1all) & str_count(w1all, ',')==1 & grepl('NHW', w1all) ~
        str_remove(w1all, '^(NHW, )*|(, NHW)*'),
      T ~ best_race2
    ),
    across(matches('best_race'), ~ factor(., levels = c('NHW', 'NHB', 'AINA', 'AAPI', 'Hispanic'))),
    # always Black vs once non-Black
    w1all = paste0(w1all, ', ', factor(h1gi6e==1,labels=c('NA','Other'))),
    w5all = paste0(w5all, ', ', factor(h5od4g==1,labels=c('NA','Other'))),
    across(c(w1all,w5all), ~ str_remove_all(., '^(NA, )*|(, NA)*')),
    all = paste(w3all, w1all, w5all, sep = ', '),
    all = str_remove_all(all, '^(NA, )*|(, NA)*'),
    all = case_when(all=='NA' ~ NA, T ~ all),
    contain_black = grepl('NHB', all),
    partial_black = case_when(
      contain_black & grepl('NHW|AAPI|AINA|Other', all) ~ 1,
      contain_black & h3od2==1 ~ 1,
      contain_black & h5od4c==1 ~ 1,
      contain_black & h1gi4==1 ~ 1,
      T ~ 0
    )
  ) %>%
  # merge in W1 constructed SES variables
  full_join(ah$ses %>% transmute(aid, nhood_ses = nhood1_d, fam_ses = sespc_al)) %>%
  select(aid, sex, byear, bmonth, paste0('age',1:4), edu, matches('race'),
         matches('_black'), skin, hair, eye, attract1 = h1ir1, attract2 = h2ir1, 
         attract3 = h3ir1, attract4 = h4ir1, wt4, lwt4, clusterid4, strataid4, 
         matches('intid[1-4]|i[34](female|edu|byear|spksp)'), i4ahexp, nhood_ses,
         fam_ses, matches('i(year|month|day)[1-4]'), all_race = all)

pheno <- link %>% 
  transmute(aid, FID = as.character(gfam), IID = as.character(gid)) %>%
  left_join(ppheno)
summary(pheno)



# 3. Export data ---------------------------------------------------------------
saveRDS(pheno, 'pheno.rds')
saveRDS(ppheno, 'all_pheno.rds')

