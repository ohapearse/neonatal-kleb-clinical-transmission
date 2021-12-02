#just playing around and making a change


library(msm)
library(tidyverse)
library(here)
library(lubridate)
library(survival)

delr <- function(id){
  df <- msm_kleb2 %>%
    filter(Neo == id) %>%
    arrange(-interval) 
  num <- df %>%
    filter(Klebsiella.pneumoniae == 2) %>%
    pull(time_point)
  num2 <- c(num, "dis")
  num2 <- num2[1]
    #slice_max(Klebsiella.pneumoniae)
  if(num2 == "dis"){
    end <- df
  } else if(num2 == "int"){
    end <- df %>%
      filter(time_point == "int" |
               time_point == "adm" |
               time_point == "birth")
  } else{
    end <- df %>%
      filter(time_point == "adm" |
               time_point == "birth")
  }
  return(end)
}

#loading data ----
kleb <- read.csv(here("data", "neo_kleb_msm_data.csv"))
kleb_baseline <- read.csv(here("data", "kleb_tind_surv.csv"))
neo_data <- read.csv(here("data", "neonatal_crf_data_qual_2_use_this.csv"))
mat_data <- read.csv(here("data", "mat_crf_qual_update.csv"))
labour_data <- read.csv(here("data", "labour_qual_updated.csv"))
dyads <- read.csv(here("data", "full_list_dyads_qual.csv"))
cot_nums <- read.csv(here("data", "cot_nums_neonates.csv"))
cot_cont <- read.csv(here("data", "cots_contamination_tmerge.csv"))
mat_cont <- read.csv(here("data", "maternal_contamination_tmerge.csv"))
abx <- read.csv(here("data", "neo_abx_data_qual.csv"))
mat_hand <- read.csv(here("data", "maternal_hand_tmerge.csv"))
location <- read.csv(here("data", "neonates_on_chatinkha_tmerge.csv"))

#prepping labour data
labour_data2 <- labour_data %>%
  select(-X.1, -X)

#prepping maternal data
mat_data2 <- mat_data  %>%
  mutate(matabxall = case_when(matabx == 1 |
                                 matcpt == 1 |
                                 matipt == 1 ~ 1,
                               TRUE ~ 0)) %>%
  select(-X.1, -X)

#prepping cot number data
cot_nums <- cot_nums %>%
  rename(neoid = Neo)

#prepping cot contamination data
cot_cont <- cot_cont %>%
  rename(neoid = Neo)

#prepping maternal contamination data
mat_cont2 <- dyads %>%
  left_join(mat_cont, by = "Mat") %>%
  select(Mat, Neo, Klebsiella.pneumoniae, interval) %>%
  rename(neoid = Neo)

#hands data
mat_hand2 <- dyads %>%
  left_join(mat_hand, by = "Mat") %>%
  select(Mat, Neo, Klebsiella.pneumoniae, interval) %>%
  rename(neoid = Neo)

#adding some data as there is a mistake with the crf
mat_hand2[mat_hand2 == -27] <- 4

#prepping antibiotic data
abx2 <- abx %>%
  mutate(cef_all = case_when(cef == 1 |
                               cef_met == 1 |
                               cef_gent == 1 ~ 1,
                             TRUE ~ 0))

#prepping location data
location2 <- location %>%
  select(-X) %>%
  rename(neoid = Neo)

#adding birth timepoint
birth <- data.frame(
  Neo = unique(kleb$Neo),
  Klebsiella.pneumoniae = 0,
  interval = 0,
  time_point = "birth")

kleb2 <- kleb %>%
  mutate(interval = case_when(interval == 0 ~ 0.5,
                              TRUE ~ as.numeric(interval))) %>%
  select(-X)

kleb3 <- rbind(birth, kleb2)

#combining data ----
kleb_msm_neo <- kleb3 %>%
  left_join(dyads, by = "Neo") %>%
  rename(neoid = Neo, pid = Mat) %>%
  left_join(neo_data, by = "neoid") %>%
  left_join(labour_data2, by = "neoid") %>%
  left_join(mat_data2, by = "pid") %>%
  select(-X.x, -X.y) %>%
  mutate(wt = wt / 1000) %>%
  mutate(outborn = case_when(qech == 1 ~ 0,
                             qech == 0 ~ 1)) %>%
  mutate(notbreast = case_when(breast == 1 ~ 0,
                               breast == 0 ~ 1)) %>%
  mutate(preterm = case_when(term == 1 ~ 0,
                             term == 0 ~ 1)) %>%
  mutate(lbw = case_when(wt <= 2.5 ~ 1,
                         TRUE ~ 0))

#prepping a data frame to test tmerge
surv_test <- kleb_msm_neo %>%
  #select(neoid, pid, status, time, wt, sex) %>%
  filter(!is.na(interval))

kleb4 <- kleb2 %>%
  rename(neoid = Neo,
         status = Klebsiella.pneumoniae,
         time = interval) %>%
  filter(!is.na(neoid),
         neoid != "RNKZH")


#saving time independent data frame
#write.csv(surv_test, here("outputs", "kleb_tind_surv.csv"))

#creating the baseline tmerge dataframe
surv_td <- tmerge(data1 = kleb_baseline, data2 = kleb_baseline, id = neoid, colonisation = event(time, status))
surv_td <- tmerge(data1 = surv_td, data2 = kleb4, id = neoid, colonisation = event(time, status))
surv_td <- tmerge(data1 = surv_td, data2 = mat_cont2, id = neoid, matcont = tdc(interval, Klebsiella.pneumoniae), options=list(na.rm=FALSE))

#adding in the antibiotic data
surv_td <- tmerge(data1 = surv_test, data2 = abx, id = neoid, abx = tdc(time, abx))
surv_td <- tmerge(data1 = surv_td, data2 = abx, id = neoid, xpen_gent = tdc(time, xpen_gent))
surv_td <- tmerge(data1 = surv_td, data2 = abx2, id = neoid, cef_all = tdc(time, cef_all))
surv_td <- surv_td %>% mutate(abx = replace_na(abx, 0),
                              xpen_gent = replace_na(xpen_gent, 0),
                              cef_all = replace_na(cef_all, 0))
#adding number of neonates per cot as a time dependent variable
surv_td <- tmerge(data1 = surv_td, data2 = cot_nums, id = neoid, cot_nums = tdc(interval, cotneonum), options=list(na.rm=FALSE))
#adding the contamination of the cot as a time dependent variable
surv_td <- tmerge(data1 = surv_td, data2 = cot_cont, id = neoid, cot_cont = tdc(interval, `Klebsiella.pneumoniae`))
#adding the maternal contamination as a time dependent variable
surv_td <- tmerge(data1 = surv_td, data2 = mat_cont2, id = neoid, mat_cont = tdc(interval, `Klebsiella.pneumoniae`))
#adding maternal hand contamination as a time dependent variable
surv_td <- tmerge(data1 = surv_td, data2 = mat_hand2, id = neoid, mat_hand = tdc(interval, `Klebsiella.pneumoniae`))
#adding the location
surv_td <- tmerge(data1 = surv_td, data2 = location2, id = neoid, onward = tdc(wardtime, wardstat))
#adding the enumerator 
surv_td <- tmerge(data1 = surv_td, data2 = surv_td, id = neoid, enum = cumtdc(tstart))

#this is now replacing some of the maternal contamination NAs with 0s as it is likely that 
#if a mother is negative in the first swab they were also negative before that
surv_td_2 <- surv_td %>%
  mutate(mat_cont = case_when(neoid == "DHCKN" |
                                neoid == "VNKBR" |
                                neoid == "VNKRD" |
                                neoid == "ZRKND" |
                                neoid == "ZVDSC" ~ as.integer(0),
                              TRUE ~ mat_cont))

msm_kleb2 <- surv_test %>%
  mutate(Klebsiella.pneumoniae = case_when(Klebsiella.pneumoniae == 0 ~ 1,
                                           Klebsiella.pneumoniae == 1 ~ 2)) %>%
  filter(!is.na(interval)) %>%
  rename(Neo = neoid) %>%
  arrange(Neo, interval)

#creating a data frame which adds the assumption of negativity at birth

msm_kleb3 <- do.call("rbind", lapply(X = unique(msm_kleb2$Neo), FUN = delr)) %>%
  arrange(Neo, interval)

#contingency table of state changes ----
con_tab <- statetable.msm(Klebsiella.pneumoniae, Neo, msm_kleb3)

#no transitions back to 0
Q <- matrix(c(-0.13, 0, 0.13, 0), ncol = 2)
Q2 <- matrix(c(-1, 0, 1, 0), ncol = 2)
Q3 <- matrix(c(-0.1, 0, 0.1, 0), ncol = 2)


#transitions back to 0
Qback <- matrix(c(-0.2, 0.1, 0.2, -0.1), ncol = 2)

#starting values for state changes ----
Q.crude <- crudeinits.msm(formula = Klebsiella.pneumoniae ~ interval, 
                          subject = Neo, 
                          data = msm_kleb3,
                          qmatrix = Q)

#running msm ----
#this one allows transitions back to uncolonised
kleb_msmback <- msm(formula = Klebsiella.pneumoniae ~ interval,
                 subject = Neo, 
                 data = msm_kleb2,
                 qmatrix = Qback)

#this one does not allow transitions back to uncolonised
kleb_msm3 <- msm(formula = Klebsiella.pneumoniae ~ interval,
                 subject = Neo, 
                 data = msm_kleb3,
                 qmatrix = Q,
                 covariates = ~ notbreast + 
                   o2 + matabxall + dlivmod +
                   sex + lbw + outborn + preterm +
                   outborn + intabx + kmc + mathiv)

estimates <- exp(kleb_msm3$estimates)
names(estimates) <- c("Baseline", "notbreast", "o2", "matabxall",
                         "dlivmod", "sex", "lbw", "outborn", "preterm",
                         "intabx", "kmc", "mathiv")

#Warning message:
#  In msm.check.model(mf.trans$"(fromstate)", mf.trans$"(tostate)",  :
#                       Absorbing - absorbing transition at observation 16
#this warning indicates that the participant has transitioned from colonised 
#to colonised which 

#- 2 * log-likelihood - lower is better
#taking out the transitions back to uncolonised reduces the log likelihood
#does this mean it fits better than having transitions back to uncolonised

#this one has some transition points
kleb_msm2 <- msm(formula = Klebsiella.pneumoniae ~ interval,
                subject = Neo, 
                data = msm_kleb3,
                qmatrix = Q.crude)

#estimating the sojourn time in each state ----
kleb_time <- sojourn.msm(kleb_msm3)

#estimating the transition probabilities in each state ----
oneday <- pmatrix.msm(kleb_msm3, t = 1)
fiveday <- pmatrix.msm(kleb_msm3, t = 5)

#hazard rates for each transition ----
kleb_haz <- hazard.msm(kleb_msm3)

#plotting the prevalence ----
plot.prevalence.msm(kleb_msm3)


#pci allows you to model transition intensities as a piecewise constant
#pci = c(t1 of change, t2 of change ...)
