library(tidyverse)

files <- list(
  pq = list(
    x86 = list(
      no = "./data/input/new/pq/x86/No-Dilitium-r2",
      single = "./data/input/new/pq/x86/Single-Dilitium-r2",
      double = "./data/input/new/pq/x86/Double-Dilitium-r2"  
    ),
    arm = list(
      no = "./data/input/new/pq/arm/No-Dilitium",
      single = "./data/input/new/pq/arm/Single-Dilitium",
      double = "./data/input/new/pq/arm/Double-Dilitium"
    )
  ),
  classic = list(
    x86 = list(
      no = "./data/input/new/classic/x86/Classic-No-Sign.csv",
      single = "./data/input/new/classic/x86/Classic-Single-Sign.csv",
      double = "./data/input/new/classic/x86/Classic-Double-Sign.csv"  
    ),
    arm = list(
      no = "./data/input/new/classic/arm/Classic-No-Sign-Raspberry.csv",
      single = "./data/input/new/classic/arm/Classic-Single-Sign-Raspberry.csv",
      double = "./data/input/new/classic/arm/Classic-Double-Sign-Raspberry.csv"
    )
  )
)

# Predictable sampling seed
set.seed(7)
dat <- files %>% imap_dfr(function(archs, protocol) {
  archs %>% imap_dfr(function(schemes, arch) {
    schemes %>% imap_dfr(function(file, scheme) {
      read_csv(file, 
        col_names = c(
          "auth",
          "exchange",
          "cypher",
          "total"
        )) %>% 
        bind_rows(sample_n(., 50)) %>% 
        slice(1:1000) %>%
        mutate(iteration = 1:1000, scheme = scheme) #%>% slice(1:950)
    }) %>% mutate(arch = arch)
  }) %>% mutate(protocol = protocol) %>%
  mutate(label = if_else(iteration %in% sample(1:1000, size = 800, replace = FALSE), "train", "test"))
}) #%>% mutate(status = 1)

dat_summarised <- dat %>% group_by(protocol, arch, scheme) %>%
  summarise(
    m_ = mean(total),
    s_ = sd(total),
    n_ = n()
  )


dat <- dat %>% left_join(dat_summarised) %>%
  mutate(stotal = ((total - m_)/s_)^2   ) %>%
  mutate(speed = factor(if_else(
    dt(abs(total - m_)/(s_/sqrt(n_)), df = n_ - 1) < 0.05,
      "average", if_else(
        total - m_ < 0,
        "slow",
        "fast"
      )), levels = c("slow", "average", "fast"))) %>% 
  select(-m_, -s_, -n_, -stotal)


# Wilcox test for train-test homogeneity
# 0.9724
(wilcox.test(total ~ -1 + label, data = dat))


dat %>% write_csv("./data/output/pqtls.csv")
