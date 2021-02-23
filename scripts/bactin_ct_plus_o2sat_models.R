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
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  select(covid_positive, chest_ct_infiltrates, contains("o2"), ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94), .funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  count(chest_ct_infiltrates)


bd_impute %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  select(covid_positive, chest_ct_infiltrates, contains("o2"), ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94), .funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  count(o2_below_94)



#' #############################################
#' 
#' B-ACTIN vs O2 SATURATION
#' 
#' #############################################

bd_impute %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  select(covid_positive, chest_ct_infiltrates, contains("o2"), ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates), .funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  mutate(chest_ct_infiltrates = ifelse(is.na(chest_ct_infiltrates), "Unknown", ifelse(chest_ct_infiltrates == TRUE, "Present", "Absent"))) %>%
  mutate(o2_below_94 = ifelse(is.na(o2_below_94), "O2 Not Measured",
                              ifelse(o2_below_94 == TRUE, "O2 Saturation < 94%", "O2 Saturation \u2265 94%"))) %>%
  brm(data = ., family = gaussian,
      bactin_ct ~ 1 + o2_sat_scale,
      # prior = c(prior(student_t(3, 0, 2.5), class = Intercept),
      #            prior(student_t(3, 0, 2.5), class = b)
      # ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.9999, max_treedepth = 22),
      backend = "cmdstanr",
      seed = 16) -> m_bactin_vs_o2sat_brms


m_bactin_vs_o2sat_brms$fit %>% rstan::check_hmc_diagnostics()


m_bactin_vs_o2sat_brms %>%
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



bd_impute %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  select(covid_positive, chest_ct_infiltrates, contains("o2"), ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94), .funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -o2_below_94, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  # mutate(chest_ct_infiltrates = ifelse(is.na(chest_ct_infiltrates), "Unknown", ifelse(chest_ct_infiltrates == TRUE, "Present", "Absent"))) %>%
  # mutate(o2_below_94 = ifelse(is.na(o2_below_94), "O2 Not Measured",
  #                             ifelse(o2_below_94 == TRUE, "O2 Saturation < 94%", "O2 Saturation \u2265 94%"))) %>%
  filter(!is.na(o2_below_94)) %>%
  brm(data = ., family = gaussian,
      bactin_ct ~ 1 + o2_below_94,
      # prior = c(prior(student_t(3, 0, 2.5), class = Intercept),
      #            prior(student_t(3, 0, 2.5), class = b)
      # ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.9999, max_treedepth = 22),
      backend = "cmdstanr",
      seed = 16) -> m_bactin_vs_o2cat_brms


m_bactin_vs_o2cat_brms$fit %>% rstan::check_hmc_diagnostics()


m_bactin_vs_o2cat_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  # odds ratio scale
  mutate_at(.vars = vars(-param), .funs = list("OR" = ~ exp(.x))) %>%
  select(param, contains("OR")) %>%
  filter(param == "b_o2_below_94TRUE") %>%
  gt::gt() %>%
  #gt::fmt_scientific(columns = 6:9, decimals = 3) %>%
  gt::fmt_number(columns = 2:6, n_sigfig = 3) %>%
  # gt::data_color(
  #   columns = vars(Estimate, Q2.5, Q97.5),
  #   colors = scales::col_bin(
  #     palette = c("#FB6467FF","white","#B7E4F9FF"),
  #     domain = NULL,
  #     na.color = NA,
  #     bins = c(-Inf,0,Inf)
  #   )
  # ) %>%
  # identity() #%>%
gt::data_color(
  columns = vars(Estimate_OR, Q2.5_OR, Q97.5_OR),
  colors = scales::col_bin(
    palette = c("#FB6467FF","white","#B7E4F9FF"),
    domain = NULL,
    na.color = NA,
    bins = c(-Inf,1,Inf)
  )
)





#' #############################################
#' 
#' SCALED & RAW DATA MODELS
#' 
#' #############################################

#' #############################################
#' 
#' * model - binomial with Ct scaled MIXED O2 sat effect - impute
#' 
#' #############################################

bd_impute %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  filter(!is.na(o2_below_94)) %>%
  count(covid_positive, o2_below_94)


#' brms
bd_impute %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  select(covid_positive, chest_ct_infiltrates, contains("o2"), ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates), .funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  mutate(chest_ct_infiltrates = ifelse(is.na(chest_ct_infiltrates), "Unknown", ifelse(chest_ct_infiltrates == TRUE, "Present", "Absent"))) %>%
  mutate(o2_below_94 = ifelse(is.na(o2_below_94), "O2 Not Measured",
                              ifelse(o2_below_94 == TRUE, "O2 Saturation < 94%", "O2 Saturation \u2265 94%"))) %>%
  #filter(!is.na(o2_below_94)) %>%
  brm(data = ., family = bernoulli,
      covid_positive ~ (bactin_ct_scale + 1 | o2_below_94),
      prior = c(prior(normal(0, 0.5), class = Intercept),
                prior(normal(0, 0.5), class = sd),
                prior(lkj(1), class = cor)
      ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.9999, max_treedepth = 22),
      backend = "cmdstanr",
      seed = 16) -> m_logit_scale_o2sat_brms


m_logit_scale_o2sat_brms %>% write_rds(file = "./models/binomial/m_logit_scale_o2sat_brms.rds.gz", compress = "gz")
m_logit_scale_o2sat_brms$fit %>% write_rds(file = "./models/binomial/m_logit_scale_o2sat_brms_stanfit.rds.gz", compress = "gz")
m_logit_scale_o2sat_brms <- read_rds(file = "./models/binomial/m_logit_scale_o2sat_brms.rds.gz")

m_logit_scale_o2sat_brms$formula


pp_check(m_logit_scale_o2sat_brms)
m_logit_scale_o2sat_brms$fit -> m_logit_scale_o2sat_stan
rstan::check_hmc_diagnostics(m_logit_scale_o2sat_stan)



# binomial model results on OR scale
m_logit_scale_o2sat_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  # odds ratio scale
  mutate_if(.predicate = ~ is.numeric(.x), .funs = list("OR" = ~ exp(.x)))



# binomial model results on probability scale

m_logit_scale_o2sat_brms$data %>%
  count(o2_below_94) %>%
  mutate(proportion = n / sum(n, na.rm = TRUE)) %>%
  gt::gt_preview() %>%
  gt::fmt_percent(columns = 4, decimals = 1)


m_logit_scale_o2sat_brms$data %>%
  tibble() %>%
  expand(bactin_ct_scale = modelr::seq_range(bactin_ct_scale, n = 100),
         o2_below_94 = unique(o2_below_94)) %>% 
  tidybayes::add_fitted_draws(model = m_logit_scale_o2sat_brms) %>%
  mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  mutate(o2_cat_number = case_when(o2_below_94 == "O2 Saturation < 94%" ~ 93,
                                   o2_below_94 == "O2 Saturation \u2265 94%" ~ 332,
                                   o2_below_94 == "O2 Not Measured" ~ 857,
                                   ),
         o2_below_94 = glue::glue("{o2_below_94} (n={o2_cat_number})")
         ) %>%
  mutate(o2_below_94 = gsub("O2","O<sub>2</sub>", o2_below_94)) %>%
  mutate(o2_below_94 = factor(o2_below_94), levels = rev(sort(unique(o2_below_94)))) %>%
  ggplot(data = ., aes(x = bactin_ct, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(facets = ~ o2_below_94, scales = "free") +
  theme_bw() +
  theme(legend.position = c(0.93,0.81),
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
  ) -> p_logit_scale_o2sat_brms
p_logit_scale_o2sat_brms



# p_logit_scale_o2sat_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_o2sat_brms.pdf", device = cairo_pdf, height = 5, width = 7, units = "in")
# p_logit_scale_o2sat_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_o2sat_brms.svg", height = 5, width = 7, units = "in")
# p_logit_scale_o2sat_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_o2sat_brms.png", height = 5, width = 7, units = "in", dpi = 600)






#' #############################################
#' 
#' * model - binomial with Ct scaled ADJUST for o2sat - impute
#' 
#' #############################################

#' brms
bd_impute %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(covid_positive = as.numeric(grepl("POS",result))) %>%
  select(covid_positive, contains("o2"), chest_ct_infiltrates, ms2_ct, bactin_ct, sarscov2_ct) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -contains("o2")), .funs = list("flip" = ~ -log(.x)),) %>%
  mutate_at(.vars = vars(-covid_positive, -chest_ct_infiltrates, -contains("o2"), -contains("flip")), .funs = list("scale" = ~ scale(.x)[,1])) %>%
  mutate(chest_ct_infiltrates = ifelse(is.na(chest_ct_infiltrates), "Unknown", ifelse(chest_ct_infiltrates == TRUE, "Present", "Absent"))) %>%
  #filter(chest_ct_infiltrates != "Unknown") %>%
  mutate(chest_ct_infiltrates = chest_ct_infiltrates == "Present") %>%
  filter(!is.na(o2_below_94)) %>%
  brm(data = ., family = bernoulli,
      covid_positive ~ (1 + bactin_ct_scale + o2_below_94),
      # prior = c(prior(student_t(3, 0, 2.5), class = Intercept),
      #           prior(student_t(3, 0, 2.5), class = b)
      # ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 22),
      backend = "cmdstanr",
      seed = 16) -> m_logit_scale_o2sat_adjust_brms


m_logit_scale_o2sat_adjust_brms %>% write_rds(file = "./models/binomial/m_logit_scale_o2sat_adjust_brms.rds.gz", compress = "gz")
m_logit_scale_o2sat_adjust_brms$fit %>% write_rds(file = "./models/binomial/m_logit_scale_o2sat_adjust_brms_stanfit.rds.gz", compress = "gz")
m_logit_scale_o2sat_adjust_brms <- read_rds(file = "./models/binomial/m_logit_scale_o2sat_adjust_brms.rds.gz")

m_logit_scale_o2sat_adjust_brms$formula


pp_check(m_logit_scale_o2sat_adjust_brms)
m_logit_scale_o2sat_adjust_brms$fit -> m_logit_scale_o2sat_adjust_stan
rstan::check_hmc_diagnostics(m_logit_scale_o2sat_adjust_stan)



# binomial model results on OR scale
m_logit_scale_o2sat_adjust_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  # odds ratio scale
  mutate_if(.predicate = ~ is.numeric(.x), .funs = list("OR" = ~ exp(.x)))


m_logit_scale_o2sat_adjust_brms %>%
  #get_variables()
  tidybayes::spread_draws(b_bactin_ct_scale) %>%
  mutate_at(.vars = vars(b_bactin_ct_scale), .funs = list("OR" = ~ exp(.x))) %>%
  median_hdi()

#OR 0.8 (95%CI 0.458 - 1.18)



# binomial model results on probability scale

m_logit_scale_o2sat_adjust_brms$data %>%
  count(o2_below_94) %>%
  gt::gt_preview()


m_logit_scale_o2sat_adjust_brms$data %>%
  tibble() %>%
  expand(bactin_ct_scale = modelr::seq_range(bactin_ct_scale, n = 100),
         o2_below_94 = c(TRUE,FALSE)
  ) %>% 
  mutate(bactin_ct = bactin_ct_scale * sd(bd_impute$bactin_ct) + mean(bd_impute$bactin_ct)) %>%
  tidybayes::add_fitted_draws(model = m_logit_scale_o2sat_adjust_brms) %>%
  # mutate(bactin_ct = exp(-bactin_ct_scale),
  #        .value = exp(.value * -1)) %>%
  ggplot(data = ., aes(x = bactin_ct, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  theme_bw() +
  theme(legend.position = c(0.87,0.70),
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
  ) -> p_logit_scale_o2sat_adjust_brms
p_logit_scale_o2sat_adjust_brms



# p_logit_scale_o2sat_adjust_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_o2sat_adjust_brms.pdf", device = cairo_pdf, height = 3, width = 4, units = "in")
# p_logit_scale_o2sat_adjust_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_o2sat_adjust_brms.svg", height = 3, width = 4, units = "in")
# p_logit_scale_o2sat_adjust_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_o2sat_adjust_brms.png", height = 3, width = 4, units = "in", dpi = 600)








#' #############################################
#' 
#' COMBINE O2 Sat & CHEST CT PLOTS
#' 
#' #############################################

library(patchwork)

(p_logit_scale_radiology_brms + theme(legend.position = "none")) /
  (p_logit_scale_o2sat_brms + theme(legend.position = c(0.22,0.69))) +
  plot_annotation(tag_levels = "A") -> p_logit_scale_ct_o2_brms
p_logit_scale_ct_o2_brms



# p_logit_scale_ct_o2_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_ct_o2_brms.pdf", device = cairo_pdf, height = 7, width = 7, units = "in")
# p_logit_scale_ct_o2_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_ct_o2_brms.svg", height = 7, width = 7, units = "in")
# p_logit_scale_ct_o2_brms %>%
#   ggsave(plot = ., filename = "./figs/p_logit_scale_ct_o2_brms.png", height = 7, width = 7, units = "in", dpi = 600)






#' #############################################
#' 
#' LINEAR MODEL ADJUSTED FOR O2
#' 
#' #############################################

#' #' brms
bd_dat %>%
  select(bactin_ct, sarscov2_ct, o2_sat) %>%
  mutate_all(.funs = list("scale" = ~ scale(.x)[,1])) %>%
  mutate(o2_below_94 = o2_sat < 0.94) %>%
  mutate(o2_below_94 = ifelse(is.na(o2_below_94), "O2 Not Measured",
                              ifelse(o2_below_94 == TRUE, "O2 Saturation < 94%", "O2 Saturation \u2265 94%"))) %>%
  brm(data = ., family = gaussian,
      sarscov2_ct ~ 1 + bactin_ct_scale + o2_sat_scale,
      #prior = c(prior(student_t(3,-4, 2.5), class = Intercept),
      #          prior(student_t(3,0, 2.5), class = b)
      #),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_lin_raw_o2_brms


m_lin_raw_o2_brms %>% write_rds(file = "./models/binomial/m_lin_raw_o2_brms.rds.gz", compress = "gz")
m_lin_raw_o2_brms$fit %>% write_rds(file = "./models/binomial/m_lin_raw_o2_brms_stanfit.rds.gz", compress = "gz")
m_lin_raw_o2_brms <- read_rds(file = "./models/binomial/m_lin_raw_o2_brms.rds.gz")

m_lin_raw_o2_brms$formula

pp_check(m_lin_raw_o2_brms)
m_lin_raw_o2_brms$fit -> m_lin_raw_o2_stan
rstan::check_hmc_diagnostics(m_lin_raw_o2_stan)


paletteer::paletteer_d(palette = "ggsci::red_material",n=5)
paletteer::paletteer_d(palette = "nord::frost",n=4)
paletteer::paletteer_d(palette = "ggsci::schwifty_rickandmorty",n=5)

m_lin_raw_o2_brms %>%
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
  # gt::data_color(
  #   columns = vars(Estimate, Q2.5, Q97.5),
  #   colors = scales::col_numeric(
  #     # custom defined values - notice that order matters!
  #     #palette = c("#ffffff", "#f2fbd2", "#c9ecb4", "#93d3ab", "#35b0ab"),
  #     palette = as.character(paletteer::paletteer_d("ggsci::red_material", n = 5)),
  #     domain = NULL
  #   )
  # )


