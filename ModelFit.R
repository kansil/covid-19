library(dplyr, warn.conflicts = FALSE)
library(readr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(rstanarm)

args <- commandArgs(trailingOnly = TRUE)


csse_cases_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
csse_deaths_raw <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")

countries <- c("Turkey","Italy","Germany","Belgium", "Spain", "United Kingdom", "France", "Sweden")
populations <- data.frame(Country=countries, Population=c(84230075,60471893 ,83753431, 11583871, 46752725, 67847297, 65257889, 10092420 ))



# Data processing

cases_temp <- pivot_longer(csse_cases_raw, cols = contains("/2"), names_to="datestr", values_to="CumCases") %>%
  select(datestr=datestr, Country="Country/Region", Region="Province/State", CumCases=CumCases) %>%
  filter(Country %in% countries) %>%
  filter(is.na(Region))%>%
  select(-Region) %>%
  mutate(Date = mdy(datestr)) %>%
  full_join(populations) %>% 
  mutate(CumCasesPerCap = CumCases/Population) %>%
  group_by(Country) %>%
  mutate(Cases = diff(c(0,CumCases))) %>% 
  mutate(CasesPerCap = Cases/Population) 

deaths_temp <- pivot_longer(csse_deaths_raw, cols = contains("/2"), names_to="datestr", values_to="CumDeaths") %>%
  select(datestr=datestr, Country="Country/Region", Region="Province/State", CumDeaths=CumDeaths) %>%
  filter(Country %in% countries) %>%
  filter(is.na(Region))%>%
  select(-Region) %>%
  mutate(Date = mdy(datestr)) %>%
  full_join(populations) %>% 
  mutate(CumDeathsPerCap = CumDeaths/Population) %>%
  group_by(Country) %>%
  mutate(Deaths = diff(c(0,CumDeaths))) %>% 
  mutate(DeathsPerCap = Deaths/Population) 

cases_datum <- cases_temp %>%
  group_by(Country) %>%
  summarize(EpidemicDatum = Date[min(which(CumCasesPerCap >= 0.3/1e6))])

deaths_datum <- deaths_temp %>%
  group_by(Country) %>%
  summarize(EpidemicDatum = Date[min(which(CumDeathsPerCap >= 0.3/1e6))])

if (length(args) > 0) {
  lastknownday <- ymd(args[1])  
} else {
  lastknownday <- ymd((cases_temp %>% filter(Country == "Turkey") %>% arrange(Date) %>% tail(1))$Date)
}

if (length(args) > 1) {
  horizon <- ymd(args[2])  
} else {
  horizon <- 21
}


cases_final <- cases_temp %>% merge(cases_datum) %>% 
  mutate(EpidemicDay = as.numeric(Date - EpidemicDatum)) %>% 
  select(-EpidemicDatum) %>%
  filter(CumCases >= 0, EpidemicDay >= 0) %>% 
  filter(Cases >= 0) %>%
  filter(Date <= lastknownday) %>%
  select(Country, Date, CumCases, CumCasesPerCap, Cases, CasesPerCap, Population, EpidemicDay)

deaths_final <- deaths_temp %>% merge(deaths_datum) %>% 
  mutate(EpidemicDay = as.numeric(Date - EpidemicDatum)) %>% 
  select(-EpidemicDatum) %>%
  filter(CumDeaths >= 0, EpidemicDay >= 0) %>% 
  filter(Deaths >= 0) %>%
  filter(Date <= lastknownday) %>%
  select(Country, Date, CumDeaths, CumDeathsPerCap, Deaths, DeathsPerCap, Population, EpidemicDay)

# MCMC Sampling

options(mc.cores = parallel::detectCores())

burnin <- 2000
totaliter <- 5000

case_model <- stan_glmer.nb(
  Cases ~ (poly(EpidemicDay, 2) + poly(EpidemicDay, 2) || Country),
  offset = log(Population), data=cases_final,
  chains = 12, iter = totaliter, warmup = burnin, adapt_delta = 0.99)

case_posterior <- posterior_predict(case_model, newdata = cases_final)

death_model <- stan_glmer.nb(
  Deaths ~ (poly(EpidemicDay, 2) + poly(EpidemicDay, 2) || Country),
  offset = log(Population), data=deaths_final,
  chains = 12, iter = totaliter, warmup = burnin, adapt_delta = 0.99)

death_posterior <- posterior_predict(death_model, newdata = deaths_final)


# Generate projections and PI bands

case_projections <- cases_final %>% group_by(Country) %>%
         complete(Date = full_seq(c(min(Date),lastknownday+horizon),1)) %>%
         ungroup() %>%
         arrange(Country, Date) %>%
         fill(Country, CumCases, CumCasesPerCap, Cases, CasesPerCap, Population) %>% 
         group_by(Country) %>%
         arrange(Date) %>%
         mutate(EpidemicDay= 0:(n()-1)) %>%
         ungroup %>%
         arrange(Country,Date)
         
case_pp <- posterior_predict(case_model, newdata = case_projections)

case_med <- apply(case_pp, 2, quantile, probs = 0.5)
case_lo_95 <- apply(case_pp, 2, quantile, probs = 0.025)
case_hi_95 <- apply(case_pp, 2, quantile, probs = 0.975)
case_lo_99 <- apply(case_pp, 2, quantile, probs = 0.005)
case_hi_99 <- apply(case_pp, 2, quantile, probs = 0.995)
case_lo_80 <- apply(case_pp, 2, quantile, probs = 0.10)
case_hi_80 <- apply(case_pp, 2, quantile, probs = 0.90)
case_lo_60 <- apply(case_pp, 2, quantile, probs = 0.20)
case_hi_60 <- apply(case_pp, 2, quantile, probs = 0.80)

case_bands <- cbind(case_projections, case_med,
                    case_lo_60, case_hi_60,
                    case_lo_80,case_hi_80,
                    case_lo_95, case_hi_95,
                    case_lo_99,case_hi_99)

death_projection_template <- deaths_final %>% group_by(Country) %>%
         complete(Date = full_seq(c(min(Date),lastknownday+horizon),1)) %>%
         ungroup() %>%
         arrange(Country, Date) %>%
         fill(Country, CumDeaths, CumDeathsPerCap, Deaths, DeathsPerCap, Population) %>% 
         group_by(Country) %>%
         arrange(Date) %>%
         mutate(EpidemicDay= 0:(n()-1)) %>%
         ungroup %>%
         arrange(Country,Date)

death_pp <- posterior_predict(death_model, newdata = death_projections)

death_med <- apply(death_pp, 2, quantile, probs = 0.5)
death_lo_95 <- apply(death_pp, 2, quantile, probs = 0.025)
death_hi_95 <- apply(death_pp, 2, quantile, probs = 0.975)
death_lo_99 <- apply(death_pp, 2, quantile, probs = 0.005)
death_hi_99 <- apply(death_pp, 2, quantile, probs = 0.995)
death_lo_80 <- apply(death_pp, 2, quantile, probs = 0.10)
death_hi_80 <- apply(death_pp, 2, quantile, probs = 0.90)
death_lo_60 <- apply(death_pp, 2, quantile, probs = 0.20)
death_hi_60 <- apply(death_pp, 2, quantile, probs = 0.80)

death_bands <- cbind(death_projections, death_med,
                     death_lo_60, death_hi_60,
                     death_lo_80, death_hi_80,
                     death_lo_95, death_hi_95,
                     death_lo_99, death_hi_99)

# Save Model Data

save(list=c("horizon",
            "lastknownday",
            "death_model",
            "case_model",
            "case_bands",
            "death_bands",
            "deaths_datum",
            "cases_datum",
            "populations",
            "countries",
            "cases_final",
            "deaths_final",
            "death_pp",
            "case_pp",
            "latest_cumcase",
            "latest_cumdeath",
            "case_posterior",
            "death_posterior"), 
     file = sprintf("models/model-data-%s.Rdata", lastknownday))

# Write Projection Tables

tr_proj_case <- case_bands %>% filter(Country == "Turkey", Date > lastknownday) %>% 
  select(Date, EpidemicDay, case_med, case_lo_60, case_hi_60,
         case_lo_80, case_hi_80, case_lo_95, case_hi_95, 
         case_lo_99, case_hi_99) %>%
  mutate(cum_case_lo_60 = cumsum(case_lo_60) + latest_cumcase,
         cum_case_hi_60 = cumsum(case_hi_60) + latest_cumcase,
         cum_case_lo_80 = cumsum(case_lo_80) + latest_cumcase,
         cum_case_hi_80 = cumsum(case_hi_80) + latest_cumcase,
         cum_case_lo_95 = cumsum(case_lo_95) + latest_cumcase,
         cum_case_hi_95 = cumsum(case_hi_95) + latest_cumcase,
         cum_case_lo_99 = cumsum(case_lo_99) + latest_cumcase,
         cum_case_hi_99 = cumsum(case_hi_99) + latest_cumcase)

write_csv(tr_proj_case, sprintf("tables/case-projections-%s.csv", lastknownday))

tr_proj_death <- death_bands %>% filter(Country == "Turkey", Date > lastknownday) %>% 
  select(Date, EpidemicDay, death_med, death_lo_60, death_hi_60,
         death_lo_80, death_hi_80, death_lo_95, death_hi_95, 
         death_lo_99, death_hi_99) %>%
  mutate(cum_death_lo_60 = cumsum(death_lo_60) + latest_cumdeath,
         cum_death_hi_60 = cumsum(death_hi_60) + latest_cumdeath,
         cum_death_lo_80 = cumsum(death_lo_80) + latest_cumdeath,
         cum_death_hi_80 = cumsum(death_hi_80) + latest_cumdeath,
         cum_death_lo_95 = cumsum(death_lo_95) + latest_cumdeath,
         cum_death_hi_95 = cumsum(death_hi_95) + latest_cumdeath,
         cum_death_lo_99 = cumsum(death_lo_99) + latest_cumdeath,
         cum_death_hi_99 = cumsum(death_hi_99) + latest_cumdeath)

write_csv(tr_proj_death, sprintf("tables/death-projections-%s.csv", lastknownday))


# Plotting

## Daily Plots

theme_set(theme_classic())

casedatum <- ymd((cases_datum %>% filter(Country == "Turkey"))$EpidemicDatum)

ggplot(case_bands %>% filter(Country == "Turkey", Date > lastknownday)) + 
  theme(legend.position = "right") +
  geom_ribbon(aes(x=Date, ymin=case_lo_99, ymax=case_hi_99, fill="99% PI")) + 
  geom_ribbon(aes(x=Date, ymin=case_lo_95, ymax=case_hi_95, fill="95% PI")) + 
  geom_ribbon(aes(x=Date, ymin=case_lo_80, ymax=case_hi_80, fill="80% PI")) + 
  geom_point(data= cases_final %>% filter(Country == "Turkey") , aes(x=Date, y=Cases), col="black", shape=3) +
  geom_line(aes(x=Date, y=case_med, col="MLE")) +
  geom_vline(xintercept=seq(casedatum, ymd(lastknownday) +  horizon, by=1),alpha=0.1) +
  geom_hline(yintercept=seq(0,6000, by=100),alpha=0.1) +
  labs(y= "Cases/Day", x="Date", fill="", colour="") +
  scale_color_manual(values=c("MLE" = "black")) +
  scale_fill_manual(values=c("80% PI" = "#dcedc1", "95% PI" = "#ffd3b6", "99% PI" = "#ffaaa5")) +
  xlim(casedatum,lastknownday+horizon) +
  ggsave(sprintf("plots/cases-daily-%s.png", lastknownday), dpi=300, dev='png', height=5, width=10, units = "in")

deathdatum <- ymd((deaths_datum %>% filter(Country == "Turkey"))$EpidemicDatum)

ggplot(death_bands %>% filter(Country == "Turkey", Date > lastknownday )) + 
  theme(legend.position = "right") +
  geom_ribbon(aes(x=Date, ymin=death_lo_99, ymax=death_hi_99, fill="99% PI")) + 
  geom_ribbon(aes(x=Date, ymin=death_lo_95, ymax=death_hi_95, fill="95% PI")) + 
  geom_ribbon(aes(x=Date, ymin=death_lo_80, ymax=death_hi_80, fill="80% PI")) + 
  geom_point(data= deaths_final %>% filter(Country == "Turkey") , aes(x=Date, y=Deaths), col="black", shape=3) +
  geom_line(aes(x=Date, y=death_med, col="MLE")) +
  geom_vline(xintercept=seq(deathdatum, ymd(lastknownday) +  horizon, by=1),alpha=0.1) +
  geom_hline(yintercept=seq(0,130, by=5),alpha=0.1) +
  labs(y= "Deaths/Day", x="Date", fill="", colour="") +
  scale_color_manual(values=c("MLE" = "black")) +
  scale_fill_manual(values=c("80% PI" = "#dcedc1", "95% PI" = "#ffd3b6", "99% PI" = "#ffaaa5")) +
  ggsave(sprintf("plots/deaths-daily-%s.png", lastknownday), dpi=300, dev='png', height=5, width=10, units = "in")


## Cumulative Plots

latest_cumcase <- (cases_final %>% filter(Country == "Turkey") %>% arrange(Date) %>% tail(1))$CumCases

ylimhi = max(cumsum(case_bands %>% filter(Country == "Turkey", Date > lastknownday) %>% select(case_hi_99))+latest_cumcase)
ylimlo = (cases_final %>% filter(Country == "Turkey", Date == (lastknownday - 15)) %>% tail(1))$CumCases

ggplot(case_bands %>% filter(Country == "Turkey", Date > lastknownday))+ 
  geom_point(data= cases_final%>% filter(Country == "Turkey") , aes(x=Date, y=CumCases), col="black", shape=3) +
  geom_ribbon(aes(x=Date, ymin=cumsum(case_lo_99)+latest_cumcase, ymax=cumsum(case_hi_99)+latest_cumcase, fill="99% PI" )) + 
  geom_ribbon(aes(x=Date, ymin=cumsum(case_lo_95)+latest_cumcase, ymax=cumsum(case_hi_95)+latest_cumcase, fill="95% PI")) + 
  geom_ribbon(aes(x=Date, ymin=cumsum(case_lo_80)+latest_cumcase, ymax=cumsum(case_hi_80)+latest_cumcase, fill="80% PI")) + 
  geom_line(aes(x=Date, y=cumsum(case_med)+latest_cumcase, col="MLE")) + 
  theme(legend.position = "right")+
  xlim(lastknownday - 15,lastknownday+horizon) +
  ylim(ylimlo, ylimhi) +
  labs(y = "Cumulative Cases to Date",
       x = "",
       fill = "",
       colour = "") +
  scale_color_manual(values = c("MLE" = "black")) +
  scale_fill_manual(values=c("80% PI" = "#dcedc1", "95% PI" = "#ffd3b6", "99% PI" = "#ffaaa5")) +
  geom_vline(xintercept=seq(lastknownday - 15, lastknownday + horizon, by=1),alpha=0.1) +
  geom_hline(yintercept=seq(ylimlo,ylimhi, by=1000),alpha=0.1) +
  ggsave(sprintf("plots/cases-cumulative-%s.png", lastknownday), dpi=300, dev='png', height=5, width=10, units = "in")

latest_cumdeath <- (deaths_final %>% filter(Country == "Turkey") %>% arrange(Date) %>% tail(1))$CumDeaths

ylimhi = max(cumsum(death_bands %>% filter(Country == "Turkey", Date > lastknownday) %>% select(death_hi_99))+latest_cumdeath)
ylimlo = (deaths_final %>% filter(Country == "Turkey", Date == (lastknownday - 15)) %>% tail(1))$CumDeaths

ggplot(death_bands %>% filter(Country == "Turkey", Date > lastknownday))+ 
  geom_point(data= deaths_final%>% filter(Country == "Turkey") , aes(x=Date, y=CumDeaths), col="black", shape=3) +
  geom_ribbon(aes(x=Date, ymin=cumsum(death_lo_99)+latest_cumdeath, ymax=cumsum(death_hi_99)+latest_cumdeath, fill="99% PI" )) + 
  geom_ribbon(aes(x=Date, ymin=cumsum(death_lo_95)+latest_cumdeath, ymax=cumsum(death_hi_95)+latest_cumdeath, fill="95% PI")) + 
  geom_ribbon(aes(x=Date, ymin=cumsum(death_lo_80)+latest_cumdeath, ymax=cumsum(death_hi_80)+latest_cumdeath, fill="80% PI")) + 
  geom_line(aes(x=Date, y=cumsum(death_med)+latest_cumdeath, col="MLE")) + 
  theme(legend.position = "right")+
  xlim(lastknownday - 15,lastknownday+horizon) +
  ylim(ylimlo, ylimhi) +
  labs(y = "Cumulative Deaths to Date",
       x = "",
       fill = "",
       colour = "") +
  scale_color_manual(values = c("MLE" = "black")) +
  scale_fill_manual(values=c("80% PI" = "#dcedc1", "95% PI" = "#ffd3b6", "99% PI" = "#ffaaa5")) +
  geom_vline(xintercept=seq(lastknownday - 15, lastknownday + horizon, by=1),alpha=0.1) +
  geom_hline(yintercept=seq(ylimlo,ylimhi, by=20),alpha=0.1) +
  ggsave(sprintf("plots/deaths-cumulative-%s.png", lastknownday), dpi=300, dev='png', height=5, width=10, units = "in")

