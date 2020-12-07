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
  gt::gt_preview()


bd_dat %>%
  ggplot(data = ., aes(x = .panel_x, y = .panel_y)) + 
  geom_point(alpha = 0.8, shape = 16, size = 0.5) + 
  facet_matrix(vars(ms2_ct, bactin_ct, sarscov2_ct))



bd_dat %>%
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
  select(-contains("med"),-contains("iqr")) -> bd_dat_summary
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
        strip.text.x = ggtext::element_markdown(color = "black"),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black")) -> p_rtpcr
p_rtpcr





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
#' B-ACTIN vs CHEST CT INFILTRATES
#' 
#' #############################################

bd_impute %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  select(covid_positive, chest_ct_infiltrates, contains("o2"), ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94), .funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  # mutate(chest_ct_infiltrates = ifelse(is.na(chest_ct_infiltrates), "Unknown", ifelse(chest_ct_infiltrates == TRUE, "Present", "Absent"))) %>%
  # mutate(o2_below_94 = ifelse(is.na(o2_below_94), "O2 Not Measured",
  #                             ifelse(o2_below_94 == TRUE, "O2 Saturation < 94%", "O2 Saturation \u2265 94%"))) %>%
  filter(!is.na(chest_ct_infiltrates)) %>%
  brm(data = ., family = gaussian,
      bactin_ct ~ 1 + chest_ct_infiltrates,
      # prior = c(prior(student_t(3, 0, 2.5), class = Intercept),
      #            prior(student_t(3, 0, 2.5), class = b)
      # ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.9999, max_treedepth = 22),
      backend = "cmdstanr",
      seed = 16) -> m_bactin_vs_chestct_brms


m_bactin_vs_chestct_brms$fit %>% rstan::check_hmc_diagnostics()


m_bactin_vs_chestct_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  #mutate_at(.vars = vars(-param), .funs = list("OR" = ~ exp(.x))) %>%
  gt::gt() %>%
  #gt::fmt_scientific(columns = 6:9, decimals = 3) %>%
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
  identity() #%>%
# gt::data_color(
#   columns = vars(Estimate_OR, Q2.5_OR, Q97.5_OR),
#   colors = scales::col_bin(
#     palette = c("#FB6467FF","white","#B7E4F9FF"),
#     domain = NULL,
#     na.color = NA,
#     bins = c(-Inf,1,Inf)
#   )
# )









#' #############################################
#' 
#' SCALED & RAW DATA MODELS
#' 
#' #############################################

#' #############################################
#' 
#' * model - binomial with Ct scaled MIXED radiology effect - impute
#' 
#' #############################################

#' #' brms
#' bd_impute %>%
#'   mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
#'   select(covid_positive, chest_ct_infiltrates, ms2_ct, bactin_ct, sarscov2_ct) %>%
#'   mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates), .funs = list("flip" = ~ -log(.x)),) %>%
#'   mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
#'   mutate(chest_ct_infiltrates = ifelse(is.na(chest_ct_infiltrates), "Unknown", ifelse(chest_ct_infiltrates == TRUE, "Present", "Absent"))) %>%
#'   #filter(!is.na(chest_ct_infiltrates)) %>%
#'   brm(data = ., family = bernoulli,
#'       formula = covid_positive ~ (bactin_ct_scale + 1 | chest_ct_infiltrates),
#'       # prior = c(prior(student_t(3, 0, 2.5), class = Intercept),
#'       #            prior(student_t(3, 0, 2.5), class = b)
#'       # ),
#'       iter = 2000,
#'       warmup = 1000,
#'       chains = 4,
#'       cores = 4,
#'       control = list("adapt_delta" = 0.999, max_treedepth = 22),
#'       backend = "cmdstanr",
#'       seed = 16) -> m_logit_scale_radiology_brms
#' 
#' 
#' m_logit_scale_radiology_brms %>% write_rds(file = "./models/binomial/m_logit_scale_radiology_brms.rds.gz", compress = "gz")
#' m_logit_scale_radiology_brms$fit %>% write_rds(file = "./models/binomial/m_logit_scale_radiology_brms_stanfit.rds.gz", compress = "gz")
m_logit_scale_radiology_brms <- read_rds(file = "./models/binomial/m_logit_scale_radiology_brms.rds.gz")

m_logit_scale_radiology_brms$formula


pp_check(m_logit_scale_radiology_brms)
m_logit_scale_radiology_brms$fit -> m_logit_scale_radiology_stan
rstan::check_hmc_diagnostics(m_logit_scale_radiology_stan)



# binomial model results on OR scale
m_logit_scale_radiology_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  # odds ratio scale
  mutate_if(.predicate = ~ is.numeric(.x), .funs = list("OR" = ~ exp(.x)))


m_logit_scale_radiology_brms %>%
  #get_variables()
  tidybayes::spread_draws(b_bactin_ct_scale) %>%
  mutate_at(.vars = vars(b_bactin_ct_scale), .funs = list("OR" = ~ exp(.x))) %>%
  median_hdi()


# binomial model results on probability scale

m_logit_scale_radiology_brms$data %>%
  count(chest_ct_infiltrates) %>%
  gt::gt_preview()


m_logit_scale_radiology_brms$data %>%
  tibble() %>%
  expand(bactin_ct_scale = modelr::seq_range(bactin_ct_scale, n = 100),
         chest_ct_infiltrates = unique(chest_ct_infiltrates)) %>% 
  tidybayes::add_fitted_draws(model = m_logit_scale_radiology_brms) %>%
  mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  mutate(chest_ct_number = case_when(chest_ct_infiltrates == "Absent" ~ 67,
                                     chest_ct_infiltrates == "Present" ~ 45,
                                     chest_ct_infiltrates == "Unknown" ~ 1199),
         chest_ct_infiltrates = glue::glue("Chest CT Infiltrate:<br>{chest_ct_infiltrates} (n={chest_ct_number})")) %>%
  mutate(chest_ct_infiltrates_f = factor(chest_ct_infiltrates, levels = rev(sort(unique(chest_ct_infiltrates))))) %>%
  mutate(chest_ct_infiltrates_f = fct_rev(chest_ct_infiltrates_f)) %>%
  ggplot(data = ., aes(x = bactin_ct, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(facets = ~ chest_ct_infiltrates_f, scales = "free") +
  theme_bw() +
  theme(legend.position = c(0.23,0.82),
        strip.background = element_blank(),
        legend.title = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(color = "black"),
        legend.background = element_rect(color = "black", size = 0.25, fill = "white")
  ) +
  labs(x = "\u03B2-actin Ct Value",
    #x = expression(paste(beta,"-actin Ct Value")),
    y = "Probability SARS-CoV-2 Positive",
    #title = "Impact of Specimen Quality on Test Sensitivity",
    fill = "Posterior<br>Credible<br>Interval"
  ) -> p_logit_scale_radiology_brms
p_logit_scale_radiology_brms



# p_logit_scale_radiology_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_radiology_brms.pdf", device = cairo_pdf, height = 5, width = 7, units = "in")
# p_logit_scale_radiology_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_radiology_brms.svg", height = 5, width = 7, units = "in")
# p_logit_scale_radiology_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_radiology_brms.png", height = 5, width = 7, units = "in", dpi = 600)






#' #############################################
#' 
#' * model - binomial with Ct scaled ADJUST for radiology - impute
#' 
#' #############################################

#' #' brms
#' bd_impute %>%
#'   mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
#'   select(covid_positive, chest_ct_infiltrates, ms2_ct, bactin_ct, sarscov2_ct) %>%
#'   mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates), .funs = list("flip" = ~ -log(.x)),) %>%
#'   mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
#'   mutate(chest_ct_infiltrates = ifelse(is.na(chest_ct_infiltrates), "Unknown", ifelse(chest_ct_infiltrates == TRUE, "Present", "Absent"))) %>%
#'   filter(chest_ct_infiltrates != "Unknown") %>%
#'   mutate(chest_ct_infiltrates = chest_ct_infiltrates == "Present") %>%
#'   brm(data = ., family = bernoulli,
#'       covid_positive ~ (1 + bactin_ct_scale + chest_ct_infiltrates),
#'       # prior = c(prior(student_t(3, 0, 2.5), class = Intercept),
#'       #           prior(student_t(3, 0, 2.5), class = b)
#'       # ),
#'       iter = 2000,
#'       warmup = 1000,
#'       chains = 4,
#'       cores = 4,
#'       control = list("adapt_delta" = 0.999, max_treedepth = 22),
#'       backend = "cmdstanr",
#'       seed = 16) -> m_logit_scale_radiology_adjust_brms
#' 
#' 
#' m_logit_scale_radiology_adjust_brms %>% write_rds(file = "./models/binomial/m_logit_scale_radiology_adjust_brms.rds.gz", compress = "gz")
#' m_logit_scale_radiology_adjust_brms$fit %>% write_rds(file = "./models/binomial/m_logit_scale_radiology_adjust_brms_stanfit.rds.gz", compress = "gz")
m_logit_scale_radiology_adjust_brms <- read_rds(file = "./models/binomial/m_logit_scale_radiology_adjust_brms.rds.gz")

m_logit_scale_radiology_adjust_brms$formula


pp_check(m_logit_scale_radiology_adjust_brms)
m_logit_scale_radiology_adjust_brms$fit -> m_logit_scale_radiology_adjust_stan
rstan::check_hmc_diagnostics(m_logit_scale_radiology_adjust_stan)



# binomial model results on OR scale
m_logit_scale_radiology_adjust_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  # odds ratio scale
  mutate_if(.predicate = ~ is.numeric(.x), .funs = list("OR" = ~ exp(.x)))


m_logit_scale_radiology_adjust_brms %>%
  #get_variables()
  tidybayes::spread_draws(b_bactin_ct_scale) %>%
  mutate_at(.vars = vars(b_bactin_ct_scale), .funs = list("OR" = ~ exp(.x))) %>%
  median_hdi()

#OR 0.685 (95%CI 0.133 - 1.55)



# binomial model results on probability scale

m_logit_scale_radiology_adjust_brms$data %>%
  count(chest_ct_infiltrates) %>%
  gt::gt_preview()


m_logit_scale_radiology_adjust_brms$data %>%
  tibble() %>%
  expand(bactin_ct_scale = modelr::seq_range(bactin_ct_scale, n = 100),
         chest_ct_infiltrates = c(TRUE, FALSE)
  ) %>% 
  mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  tidybayes::add_fitted_draws(model = m_logit_scale_radiology_adjust_brms) %>%
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
        #axis.title.x = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(color = "black")
  ) +
  labs(#x = "\u03B2-actin Ct Value",
    x = expression(paste(beta,"-actin Ct Value")),
    y = "Probability SARS-CoV-2 Positive",
    #title = "Impact of Specimen Quality on Test Sensitivity",
    fill = "Posterior<br>Credible<br>Interval"
  ) -> p_logit_scale_radiology_adjust_brms
p_logit_scale_radiology_adjust_brms



# p_logit_scale_radiology_adjust_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_radiology_adjust_brms.pdf", height = 3, width = 4, units = "in")
# p_logit_scale_radiology_adjust_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_radiology_adjust_brms.svg", height = 3, width = 4, units = "in")
# p_logit_scale_radiology_adjust_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_radiology_adjust_brms.png", height = 3, width = 4, units = "in", dpi = 600)




