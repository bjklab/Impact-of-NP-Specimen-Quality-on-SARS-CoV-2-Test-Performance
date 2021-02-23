#' #############################################
#' 
#' load libraries and set seed
#' 
#' #############################################
library(tidyverse)
library(ggforce)
library(patchwork)
library(tidybayes)
library(brms)
set.seed(16)



#' #############################################
#' 
#' load data - create versions with missing and imputed missing
#' 
#' #############################################

bd_dat <- read_csv("./data/np_ct_plus_chest_ct_o2.csv")
bd_dat

bd_dat %>%
  # filter to MS2 target range
  filter(ms2_ct >= 20 & ms2_ct <= 25) %>%
  identity() -> bd_dat
bd_dat


bd_dat %>%
  count(!is.na(sarscov2_ct))


bd_dat %>%
  gt::gt_preview()


bd_dat %>%
  ggplot(data = ., aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.8, shape = 16, size = 0.5) + 
  facet_matrix(vars(ms2_ct, bactin_ct, sarscov2_ct))



bd_dat %>%
  summarise(`Number of Specimens` = n(),
            ms2_med = signif(median(ms2_ct, na.rm = TRUE),2),
            ms2_iqr = signif(IQR(ms2_ct, na.rm = TRUE),2),
            ms2_q1 = signif(quantile(ms2_ct, probs = c(0.25), na.rm = TRUE), 2),
            ms2_q3 = signif(quantile(ms2_ct, probs = c(0.75), na.rm = TRUE), 2),
            bactin_med = signif(median(bactin_ct, na.rm = TRUE),2),
            bactin_iqr = signif(IQR(bactin_ct, na.rm = TRUE),2),
            bactin_q1 = signif(quantile(bactin_ct, probs = c(0.25), na.rm = TRUE), 2),
            bactin_q3 = signif(quantile(bactin_ct, probs = c(0.75), na.rm = TRUE), 2),
            sarscov2_med = signif(median(sarscov2_ct, na.rm = TRUE),2),
            sarscov2_iqr = signif(IQR(sarscov2_ct, na.rm = TRUE),2),
            sarscov2_q1 = signif(quantile(sarscov2_ct, probs = c(0.25), na.rm = TRUE), 2),
            sarscov2_q3 = signif(quantile(sarscov2_ct, probs = c(0.75), na.rm = TRUE), 2),
            ms2_sum = glue::glue("{ms2_med} ({ms2_q1} - {ms2_q3})"),
            bactin_sum = glue::glue("{bactin_med} ({bactin_q1} - {bactin_q3})"),
            sarscov2_sum = glue::glue("{sarscov2_med} ({sarscov2_q1} - {sarscov2_q3})")) %>%
  rename(`MS2 Ct` = ms2_sum, `Beta-actin Ct` = bactin_sum, `SARS-CoV-2 Ct` = sarscov2_sum) %>%
  select(-contains("med"),-contains("iqr"), -contains("q1"), -contains("q3")) -> bd_dat_summary
bd_dat_summary


bd_dat_summary %>%
  gt::gt()



#' #############################################
#' 
#' IMPUTE MISSING DATA
#' 
#' #############################################

bd_dat %>%
  # impute missing data as maximum observed Ct Value (for that analyte, among actual clinical specimens)
  # mutate_at(.vars = vars(ms2_ct, bactin_ct, sarscov2_ct), .funs = ~ replace(.x, is.na(.x), max(.x, na.rm = TRUE))) %>%
  mutate_at(.vars = vars(ms2_ct, bactin_ct, sarscov2_ct), .funs = ~ replace(.x, is.na(.x), 40)) %>%
  identity() -> bd_impute

bd_impute

bd_impute %>%
  gt::gt_preview()



bd_impute %>%
  summarise(`Number of Specimens` = n(),
            ms2_med = signif(median(ms2_ct, na.rm = TRUE),2),
            ms2_iqr = signif(IQR(ms2_ct, na.rm = TRUE),2),
            bactin_med = signif(median(bactin_ct, na.rm = TRUE),2),
            bactin_iqr = signif(IQR(bactin_ct, na.rm = TRUE),2),
            sarscov2_med = signif(median(sarscov2_ct, na.rm = TRUE),2),
            sarscov2_iqr = signif(IQR(sarscov2_ct, na.rm = TRUE),2),
            ms2_sum = glue::glue("{ms2_med} ({ms2_iqr})"),
            bactin_sum = glue::glue("{bactin_med} ({bactin_iqr})"),
            sarscov2_sum = glue::glue("{sarscov2_med} ({sarscov2_iqr})")) %>%
  rename(`MS2 Ct` = ms2_sum, `Beta-actin Ct` = bactin_sum, `SARS-CoV-2 Ct` = sarscov2_sum) %>%
  select(-contains("med"),-contains("iqr")) -> bd_impute_summary
bd_impute_summary


bd_dat_summary %>%
  bind_rows(bd_impute_summary) %>%
  gt::gt()



#' plot comparison of PCR analytes
bd_impute %>%
  dplyr::rename(`MS2 Ct` = ms2_ct, `β-actin Ct` = bactin_ct, `SARS-CoV-2 Ct` = sarscov2_ct) %>%
  ggplot(data = ., aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.8, shape = 21, size = 0.5) + 
  geom_autodensity() +
  facet_matrix(vars(`MS2 Ct`, `β-actin Ct`, `SARS-CoV-2 Ct`), layer.diag = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        #strip.text = ggtext::element_markdown(color = "black"),
        #strip.text.x = ggtext::element_markdown(color = "black"),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black")) -> p_rtpcr
p_rtpcr

# p_rtpcr %>%
#   ggsave(filename = "./figs/p_rtpcr_analyte_comparison.pdf", height = 4, width = 5, units = "in", device = cairo_pdf)
# p_rtpcr %>%
#   ggsave(filename = "./figs/p_rtpcr_analyte_comparison.svg", height = 4, width = 5, units = "in")
# p_rtpcr %>%
#   ggsave(filename = "./figs/p_rtpcr_analyte_comparison.png", height = 4, width = 5, units = "in", dpi = 600)


bd_impute %>%
  dplyr::rename(`MS2 Ct` = ms2_ct, `β-actin Ct` = bactin_ct, `SARS-CoV-2 Ct` = sarscov2_ct) %>%
  ggplot(data = ., aes(x = `β-actin Ct`, y =`SARS-CoV-2 Ct`)) + 
  geom_point(alpha = 0.8, shape = 21, size = 0.5) + 
  theme_bw() +
  theme(strip.background = element_blank(),
        #strip.text = ggtext::element_markdown(color = "black"),
        #strip.text.x = ggtext::element_markdown(color = "black"),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black")) -> p_rtpcr2
p_rtpcr2




bd_impute %>%
  select(ms2_ct, bactin_ct, sarscov2_ct) %>%
  summarise(`Number of Specimens` = n(),
            ms2_med = signif(median(ms2_ct, na.rm = TRUE),2),
            ms2_iqr = signif(IQR(ms2_ct, na.rm = TRUE),2),
            bactin_med = signif(median(bactin_ct, na.rm = TRUE),2),
            bactin_iqr = signif(IQR(bactin_ct, na.rm = TRUE),2),
            sarscov2_med = signif(median(sarscov2_ct, na.rm = TRUE),2),
            sarscov2_iqr = signif(IQR(sarscov2_ct, na.rm = TRUE),2),
            ms2_sum = glue::glue("{ms2_med} ({ms2_iqr})"),
            bactin_sum = glue::glue("{bactin_med} ({bactin_iqr})"),
            sarscov2_sum = glue::glue("{sarscov2_med} ({sarscov2_iqr})")) %>%
  rename(`MS2 Ct` = ms2_sum, `Beta-actin Ct` = bactin_sum, `SARS-CoV-2 Ct` = sarscov2_sum) %>%
  select(-contains("med"),-contains("iqr")) -> bd_impute_summary
bd_impute_summary


bd_impute_summary %>%
  gt::gt()



#' #############################################
#' 
#' SCALED & RAW DATA MODELS
#' 
#' #############################################

#' #############################################
#' 
#' * model - binomial with Ct scaled - impute
#' 
#' #############################################


# bd_impute %>%
#   mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
#   select(covid_positive, ms2_ct, bactin_ct, sarscov2_ct) %>%
#   mutate_at(.vars = vars(-covid_positive), .funs = list("flip" = ~ -log(.x)),) %>%
#   mutate_at(.vars = vars(-covid_positive, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
#   brm(data = ., family = bernoulli,
#       covid_positive ~ 1 + bactin_ct_scale,
#       prior = c(prior(student_t(3, 0, 2.5), class = Intercept),
#                 prior(student_t(3, 0, 2.5), class = b)
#       ),
#       iter = 2000,
#       warmup = 1000,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.99, max_treedepth = 16),
#       backend = "cmdstanr",
#       seed = 16) -> m_logit_scale_brms
# 
# 
# m_logit_scale_brms %>% write_rds(file = "./models/binomial/m_logit_scale_brms.rds.gz", compress = "gz")
# m_logit_scale_brms$fit %>% write_rds(file = "./models/binomial/m_logit_scale_brms_stanfit.rds.gz", compress = "gz")
m_logit_scale_brms <- read_rds(file = "./models/binomial/m_logit_scale_brms.rds.gz")

m_logit_scale_brms$formula

pp_check(m_logit_scale_brms)
m_logit_scale_brms$fit -> m_logit_scale_stan
rstan::check_hmc_diagnostics(m_logit_scale_stan)


# binomial model results on OR scale
m_logit_scale_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  # odds ratio scale
  mutate_if(.predicate = ~ is.numeric(.x), .funs = list("OR" = ~ exp(.x)))


m_logit_scale_brms %>%
  #get_variables()
  tidybayes::spread_draws(b_bactin_ct_scale) %>%
  mutate_at(.vars = vars(b_bactin_ct_scale), .funs = list("OR" = ~ exp(.x))) %>%
  median_hdi()

m_logit_scale_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_at(.vars = vars(-param), .funs = list("OR" = ~ exp(.x))) %>%
  gt::gt() %>%
  gt::fmt_scientific(columns = 6:9, decimals = 3) %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3) %>%
  gt::data_color(
    columns = vars(Estimate, Q2.5, Q97.5),
    colors = scales::col_bin(
      palette = c("#FB6467FF","white","#B7E4F9FF"),
      domain = NULL,
      na.color = NA,
      bins = c(-Inf,0,Inf)
    )
  ) %>%
  gt::data_color(
    columns = vars(Estimate_OR, Q2.5_OR, Q97.5_OR),
    colors = scales::col_bin(
      palette = c("#FB6467FF","white","#B7E4F9FF"),
      domain = NULL,
      na.color = NA,
      bins = c(-Inf,1,Inf)
    )
  )





# binomial model results on probability scale

m_logit_scale_brms$data %>%
  tibble() %>%
  expand(bactin_ct_scale = modelr::seq_range(bactin_ct_scale, n = 100),
  ) %>% 
  mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  tidybayes::add_fitted_draws(model = m_logit_scale_brms) %>%
  # mutate(bactin_ct = exp(-bactin_ct_scale),
  #        .value = exp(.value * -1)) %>%
  ggplot(data = ., aes(x = bactin_ct, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  theme_bw() +
  theme(legend.position = c(0.88,0.70),
        strip.background = element_blank(),
        legend.title = ggtext::element_markdown(),
        legend.background = element_rect(color = "black", size = 0.25, fill = "white"),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(color = "black")
  ) +
  labs(x = "\u03B2-actin Ct Value",
    #x = expression(paste(beta,"-actin Ct Value")),
    y = "Probability SARS-CoV-2 Positive",
    #title = "Impact of Specimen Quality on Test Sensitivity",
    fill = "Posterior<br>Credible<br>Interval"
  ) -> p_logit_scale_brms
p_logit_scale_brms


# p_logit_scale_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_brms.pdf", device = cairo_pdf, height = 3, width = 4, units = "in")
# p_logit_scale_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_brms.svg", height = 3, width = 4, units = "in")
# p_logit_scale_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_brms.png", height = 3, width = 4, units = "in", dpi = 600)



# posterior contrast: B-actin Ct increases from 28 to 32
m_logit_scale_brms$data %>%
  tibble() %>%
  expand(bactin_ct_scale = c((28 - mean(bd_impute$bactin_ct))/sd(bd_impute$bactin_ct),(32 - mean(bd_impute$bactin_ct))/sd(bd_impute$bactin_ct)),
  ) %>% 
  mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  tidybayes::add_fitted_draws(model = m_logit_scale_brms) %>%
  ungroup() %>%
  select(bactin_ct, .value, .draw) %>%
  mutate(bactin_ct = factor(bactin_ct)) %>%
  spread(key = bactin_ct, value = .value) %>%
  rename_all(~ paste0("bactin_ct_", .x)) %>%
  mutate(contrast = bactin_ct_32 - bactin_ct_28) %>%
  median_hdi() %>%
  select(contains("contrast"), .interval)



# A tibble: 1 x 4
# contrast contrast.lower contrast.upper .interval
# <dbl>          <dbl>          <dbl> <chr>    
#   1  -0.0634        -0.0917        -0.0335 hdi    


m_logit_scale_brms$data %>%
  tibble() %>%
  expand(bactin_ct_scale = c((28 - mean(bd_impute$bactin_ct))/sd(bd_impute$bactin_ct),(32 - mean(bd_impute$bactin_ct))/sd(bd_impute$bactin_ct)),
  ) %>% 
  mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  tidybayes::add_fitted_draws(model = m_logit_scale_brms) %>%
  ungroup() %>%
  select(bactin_ct, .value, .draw) %>%
  mutate(bactin_ct = factor(bactin_ct)) %>%
  spread(key = bactin_ct, value = .value) %>%
  rename_all(~ paste0("bactin_ct_", .x)) %>%
  mutate(contrast = bactin_ct_32 - bactin_ct_28) %>%
  median_hdi() %>%
  select(contains("contrast"), .interval) %>%
  gt::gt() %>%
  gt::fmt_percent(columns = 1:3)




#' #############################################
#' 
#' * model - linear with Ct raw - impute
#' 
#' #############################################

#' brms
# bd_impute %>%
#   select(ms2_ct, bactin_ct, sarscov2_ct) %>%
#   mutate_all(.funs = list("flip" = ~ -log(.x)),) %>%
#   mutate_at(.vars = vars(-contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
#   #select(contains("scale")) %>%
#   brm(data = ., family = gaussian,
#       sarscov2_ct ~ 1 + bactin_ct,
#       #prior = c(prior(student_t(3,-4, 2.5), class = Intercept),
#       #          prior(student_t(3,0, 2.5), class = b)
#       #),
#       iter = 2000,
#       warmup = 1000,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.99, max_treedepth = 16),
#       backend = "cmdstanr",
#       seed = 16) -> m_lin_raw_brms
# 
# 
# m_lin_raw_brms %>% write_rds(file = "./models/binomial/m_lin_raw_brms.rds.gz", compress = "gz")
# m_lin_raw_brms$fit %>% write_rds(file = "./models/binomial/m_lin_raw_brms_stanfit.rds.gz", compress = "gz")
m_lin_raw_brms <- read_rds(file = "./models/binomial/m_lin_raw_brms.rds.gz")

m_lin_raw_brms$formula

pp_check(m_lin_raw_brms)
m_lin_raw_brms$fit -> m_lin_raw_stan
rstan::check_hmc_diagnostics(m_lin_raw_stan)


# linear model results
m_lin_raw_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param")

m_lin_raw_brms %>%
  #get_variables()
  tidybayes::spread_draws(b_bactin_ct) %>%
  median_hdi()

m_lin_raw_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3) %>%
  gt::data_color(
    columns = vars(Estimate, Q2.5, Q97.5),
    colors = scales::col_bin(
      palette = c("#FB6467FF","white","#B7E4F9FF"),
      domain = NULL,
      na.color = NA,
      bins = c(-Inf,0,Inf)
    )
  )



m_lin_raw_brms$data %>%
  tibble() %>%
  expand(bactin_ct = modelr::seq_range(bactin_ct, n = 100),
  ) %>% 
  #mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  tidybayes::add_fitted_draws(model = m_lin_raw_brms) %>%
  ggplot(data = ., aes(x = bactin_ct, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  theme_bw() +
  theme(legend.position = c(0.12,0.70),
        strip.background = element_blank(),
        legend.title = ggtext::element_markdown(),
        legend.background = element_rect(color = "black", size = 0.25, fill = "white"),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        #axis.title.x = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(color = "black")
  ) +
  labs(#x = "\u03B2-actin Ct Value",
    x = expression(paste(beta,"-actin Ct Value")),
    y = "SARS-CoV-2 Ct Value",
    #title = "Impact of Specimen Quality on Test Sensitivity",
    fill = "Posterior<br>Credible<br>Interval"
  ) -> p_lin_raw_brms
p_lin_raw_brms


# p_lin_raw_brms %>%
#   ggsave(plot = ., filename = "./figs/p_lin_raw_brms.pdf", height = 3, width = 4, units = "in")
# p_lin_raw_brms %>%
#   ggsave(plot = ., filename = "./figs/p_lin_raw_brms.svg", height = 3, width = 4, units = "in")
# p_lin_raw_brms %>%
#   ggsave(plot = ., filename = "./figs/p_lin_raw_brms.png", height = 3, width = 4, units = "in", dpi = 600)



#' update linear model without imputed data
bd_dat %>%
  select(ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_all(.funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  update(object = m_lin_raw_brms, newdata = .,
         iter = 2000,
         warmup = 1000,
         chains = 4,
         cores = 4,
         control = list("adapt_delta" = 0.99, max_treedepth = 16),
         backend = "cmdstanr",
         seed = 16) -> m_lin_raw_without_impute_brms



m_lin_raw_without_impute_brms$formula

pp_check(m_lin_raw_without_impute_brms)
m_lin_raw_without_impute_brms$fit -> m_lin_raw_without_impute_stan
rstan::check_hmc_diagnostics(m_lin_raw_without_impute_stan)


# linear model results
m_lin_raw_without_impute_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param")

m_lin_raw_without_impute_brms %>%
  #get_variables()
  tidybayes::spread_draws(b_bactin_ct) %>%
  median_hdi()


m_lin_raw_without_impute_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3) %>%
  gt::data_color(
    columns = vars(Estimate, Q2.5, Q97.5),
    colors = scales::col_bin(
      palette = c("#FB6467FF","white","#B7E4F9FF"),
      domain = NULL,
      na.color = NA,
      bins = c(-Inf,0,Inf)
    )
  )





m_lin_raw_without_impute_brms$data %>%
  tibble() %>%
  expand(bactin_ct = modelr::seq_range(bactin_ct, n = 100),
  ) %>% 
  tidybayes::add_fitted_draws(model = m_lin_raw_without_impute_brms) %>%
  ggplot(data = ., aes(x = bactin_ct, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  geom_point(data = bd_dat, aes(x = bactin_ct, y = sarscov2_ct), shape = 21) +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        legend.title = ggtext::element_markdown(),
        #legend.background = element_rect(color = "black", size = 0.25, fill = "white"),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        #axis.title.x = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(color = "black")
  ) +
  labs(#x = "\u03B2-actin Ct Value",
    x = expression(paste(beta,"-actin Ct Value")),
    y = "SARS-CoV-2 Ct Value",
    #title = "Impact of Specimen Quality on Test Sensitivity",
    fill = "Posterior<br>Credible<br>Interval"
  ) -> p_lin_raw_without_impute_brms
p_lin_raw_without_impute_brms












