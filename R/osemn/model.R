library("tidyverse")
library("ggpubr")
library("ggRandomForests")
library("GGally")
library("sjPlot")
library("survival")
library("survminer")
library("coxme")
library("flexsurv")
#library("coxme")
library("hermite")
#library("ARTool")
library("randomForestSRC")
library("lme4")
library("lmerTest")
library("glmulti")
library("rstatix")
library("car")
library("sjstats")
library("modelr")
library("stargazer")

(cores <- parallel::detectCores())
doParallel::registerDoParallel(cores - 1 )  ## registerDoMC( detectCores()-1 ) in Linux

options(rf.cores = cores - 1, 
        mc.cores = cores - 1)

dat_effects <- list(
  protocol =  c("classic", "pq"),
  arch = c("x86", "arm"),
  scheme = c("no", "single", "double"),
  speed = c("slow", "average", "fast")
)

dat <- read_csv("./data/output/pqtls.csv") %>%
  mutate(
    protocol = factor(protocol, levels = dat_effects$protocol),
    arch = factor(arch, levels = dat_effects$arch),
    scheme = factor(scheme, levels = dat_effects$scheme),
    speed = factor(speed, levels = dat_effects$speed),
    status = 1
  ) %>%
  as.data.frame


dat_train <- dat %>% filter(label == "train") %>% select(-label) %>% as.data.frame
dat_test  <- dat %>% filter(label == "train") %>% select(-label) %>% as.data.frame


fit_hermite <- function(data) {
  # La primera fase por máximo verosimilitud
  # equivale numéricamente a minimizar la logverosimilitud negativa
  maxlikelihood <- function(params) {
    -sum(log(dhermite(data, a = params[1], b = params[2])))
  }
  m <- mean(data)
  s <- sd(data)
  a1 <- abs(2*m - s^2)
  a2 <- (s^2 - m)/2
  partial_model <- optim(
    par = c(a = 0, b = 0), # start from max likelihood  hermite
    fn = maxlikelihood, 
    control = list(maxit = 1000000)
  )
  
  # El segunda fase por maximización de bondad de ajuste
  # Equivale numéricamente a minimizar el rechazo de la gamma
  pval <- 0
  goodness_of_fit <- function(params) {
    pval <- ad.test( # Usamos una prueba Anderson-Darling
      x = data, 
      null = "phermite",# con hipótesis nula la gamma
      a = params[1], 
      b = params[2], 
      estimated = TRUE # e hipótesis compuesta con ajuste Braun
    )$p.value
    (1 - pval) # pues la hipotesis nula es que SÍ es gamma
  }
  (final_model <- optim(
    par = partial_model$par,
    fn = goodness_of_fit, 
    method="L-BFGS-B",
    lower=c(0, 0),
    upper=c(Inf, Inf),
    control = list(maxit = 1000000)))
}
#dat_params <- dat$total %>% fit_hermite

dat_survdiff <- survdiff(Surv(total, status) ~ protocol, data = dat, rho = 1)  

dat_survfit <- survfit(Surv(total, status) ~ protocol, data = dat) 
pdf("./data/output/pqtls_survfit.pdf", 800, 400, 24)
p <- dat_survfit %>% ggsurvplot(data = dat,
    conf.int = TRUE,
    risk.table = TRUE,
    surv.median.line = "hv",
    fun = "pct",
    #break.x.by = 
    #pval.method = TRUE,
    pval = "Peto-Peto, p < 1e-03 ***",
    log.rank.weights = "S1", # peto peto
    #test.for.trend = TRUE,
    #facet.by = c("arch"),
    #color = scheme,
    #add.all = TRUE,
    ylim = c(0,51),
    break.x.by = 2e+10 / 5,
    ylab = "TLS batch proportion (%)",
    xlab = "CPU clock cycles",
    legend.title = "Protocol",
    #legend.labs = c("both", "classic", "pq"),
    #linetype = "strata",
    palette = c("#2E9FDF", "#0f6b4f"), #, "#E7B800", "#2E9FDF"),
    ggtheme = ggthemes::theme_tufte(base_size = 24),
    tables.theme = ggthemes::theme_tufte(base_size = 14))
p$plot %>% ggsave(
  filename = "./data/output/pqtls_survfit_plot.pdf",
  plot = .,
  width=600,
  height=200,
  units="mm",
  scale=1/3,
  pointsize = 24,
  dpi=1200,
)
p$table %>% ggsave(
  filename = "./data/output/pqtls_survfit_table.pdf",
  plot = .,
  width=600,
  height=100,
  units="mm",
  scale=1/3,
  pointsize=14,
  dpi=1200,
)
dev.off()

#png("./data/output/pqtls_cumhaz.png", 800, 400, units = "px", 24)
p <-dat_survfit %>% ggsurvplot(data = dat, 
  fun = "cumhaz",
  conf.int = TRUE,
  conf.int.km = TRUE,
  risk.table = TRUE, #"nrisk_cumevents", # "nrisk_cumevents" "nrisk_cumhazard"
  #pval = TRUE,
  #pval.method = "Hola",
  #log.rank.weights = "S1",
  #test.for.trend = TRUE,
  #facet.by = "arch",
  #group.by = "scheme",
  add.all = FALSE,
  ylab = "Expected TLS executions",
  xlab = "CPU clock cycles",
  legend.title = "Protocol",
  #legend.labs = c("classic", "pq"),
  #linetype = "strata",
  #palette = c("gray", "#E7B800", "#2E9FDF"),
  palette = c("#2E9FDF", "#0f6b4f"),# c("gray", "#0f6b4f", "#2E9FDF"),
  ylim = c(0,3),
  ggtheme = ggthemes::theme_tufte(base_size = 24),
  tables.theme = ggthemes::theme_tufte(14))
p$plot %>% ggsave(
  filename = "./data/output/pqtls_cumhaz_plot.pdf",
  plot = .,
  width=600,
  height=300,
  units="mm",
  scale=1/4,
  dpi=1200,
)
p$table %>% ggsave(
  filename = "./data/output/pqtls_cumhaz_table.pdf",
  plot = .,
  width=600,
  height=100,
  units="mm",
  scale=1/3,
  pointsize=14,
  dpi=1200,
)
dev.off()

p <- NULL


#srv <- survreg(Surv(total, status) ~ protocol*arch*scheme, data = dat, dist = "weibull")

dat_survreg <- flexsurvreg(Surv(total, status) ~ protocol, data = dat, dist = "weibull")
png("./data/output/pqtls_survreg.png", 400, 450, units = "px", 12)
dat_survreg %>% ggflexsurvplot(data = dat, 
  fun = "survival",
  conf.int = TRUE,
  conf.int.km = TRUE,
  risk.table = TRUE,
  #pval = TRUE,
  #pval.method = "Hola",
  #log.rank.weights = "S1",
  #test.for.trend = TRUE,
  facet.by = "arch",
  group.by = "scheme",
  add.all = FALSE,
  ylab = "Batch proportion (%)",
  xlab = "Clock cycles",
  legend.title = "Protocol",
  #legend.labs = c("classic", "pq"),
  #linetype = "strata",
  #palette = c("gray", "#E7B800", "#2E9FDF"),
  ggtheme = ggthemes::theme_tufte(base_size = 12),
  tables.theme = ggthemes::theme_tufte(base_size = 12))
dev.off()




dat_model_mixef <-  lmer(
  total ~ iteration * protocol + 
    (iteration | arch:scheme) + 
    (iteration * protocol || arch), 
  data = dat
)


(dat_model <- glm(
  total ~ protocol*arch*scheme,
  data = dat
))
(dat_model %>% mape(dat_test))

(dat_model_mixef %>% mape(dat_test))


sink("./data/output/pqtls_model.txt")
dat_model %>% stargazer(model.names = TRUE) 
sink()

png("./data/output/pqtls_model.png", 400, 450, units = "px", 12)
dat_model %>% sjPlot::plot_model( 
  axis.title = "CPU clock cycles",
  colors = c("#6D9EC1", "#5da68a"), 
  show.values = TRUE,
  value.offset = .1,
  value.size = 4,
  dot.size = 1,
  #line.size = 1.5,
  vline.color = "gray",
  width = 0.5#1.5
) + ggthemes::theme_tufte(base_size = 12)
dev.off()

png("./data/output/pqtls_residuals.png", 400, 450, units = "px", 12)
dat_model %>% ggnostic(
    mapping = ggplot2::aes(color = scheme),
    columnsY = c("total", ".fitted", ".se.fit", ".resid", ".std.resid", ".hat", ".sigma", ".cooksd"),
    continuous = list(default = ggally_smooth, .fitted = ggally_smooth), 
  ) 
#p + ggthemes::theme_tufte(base_size = 36) + scale_fill_manual()
dev.off()

############
# ANOVA
# -----

sink("./data/output/pqtls_anova.txt")
(dat_anova <- dat_model %>% Anova(
  type = 3, 
  white.adjust = "hc4",
  test.statistic = "F") %>%
  anova_stats()) %>% left_join(
    dat_model %>% Anova(
    type = 3, 
    white.adjust = "hc4",
    test.statistic = "Wald") %>% 
    tidy %>% dplyr::select(-df),
    by = "term"
  ) %>%
  dplyr::select(term, df, sumsq, meansq, statistic, omegasq) %>%
  knitr::kable(format = "latex")
sink()


png("./data/output/pqtls_arm_anova.png", 400, 475, units = "px", 36)
dat %>% filter(arch == "arm") %>% 
  ggboxplot(
    x = "protocol",
    y = "total",
    color = "scheme",
    combine = TRUE,
    #dotsize = 0.3,
    #binwidth = 0.1,
    notch = TRUE,
    alpha = 0,
    palette = c("#0f6b4f", "#388a66", "#5da68a"),
    #palette = "locus",
    add = c(
      "jitter", 
      #"boxplot",
      "median_iqr"
    ),
    add.params = list(
      alpha = 0.1,
      #notch = TRUE,
      #combine = TRUE
      #size = 0.5
      #binwidth = 0.1,
      #dotsize = 3,
      jitter = 0.1
    ),
    ylim = c(0, 2.8e+07),
    ggtheme = ggthemes::theme_tufte(base_size = 36)
    #shape = "arch"
  ) + 
  stat_compare_means(
    label.y = 2.8e+7,
    label.x = 1,
    method = "anova"
  ) +
  stat_compare_means(
    label.y = 2.7e+7, # arm is 2 orders of magnitude lesser
    #method = "anova",
    label = "p.signif",
    comparisons = list(
      c("classic", "pq")
    )
  ) + 
  stat_compare_means(
    aes(group = scheme),
    method = "anova",
    label = "p.signif",
    #ref.group = "no",
    label.y = 2.3e+7
  )
dev.off()

png("./data/output/pqtls_x86_anova.png", 400, 475, units = "px", 36)
dat %>% filter(arch == "x86") %>% 
  ggboxplot(
    x = "protocol",
    y = "total",
    color = "scheme",
    combine = TRUE,
    #dotsize = 0.3,
    #binwidth = 0.1,
    notch = TRUE,
    alpha = 0,
    palette = c("#0f6b4f", "#388a66", "#5da68a"),
    #palette = "locus",
    add = c(
      "jitter", 
      #"boxplot",
      "median_iqr"
    ),
    add.params = list(
      alpha = 0.2,
      #notch = TRUE,
      #combine = TRUE
      #size = 0.5
      #binwidth = 0.1,
      #dotsize = 3,
      jitter = 0.1
    ),
    #ylim = c(0, 2.8e+07),
    ggtheme = ggthemes::theme_tufte(base_size = 36)
    #shape = "arch"
  ) + 
  stat_compare_means(
    label.y = 2.2e+10,
    label.x = 1,
    method = "anova"
  ) +
  stat_compare_means(
    label.y = 2.1e+10, # arm is 2 orders of magnitude lesser
    #method = "anova",
    label = "p.signif",
    comparisons = list(
      c("classic", "pq")
    )
  ) + 
  stat_compare_means(
    aes(group = scheme),
    method = "anova",
    label = "p.signif",
    #ref.group = "no",
    label.y = 2.0e+10
  )
dev.off()

