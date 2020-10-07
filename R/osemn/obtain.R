library(tidyverse)

files <- list(
  pq = list(
    x86 = list(
      no = "./data/input/pq/x86/no.csv",
      single = "./data/input/pq/x86/single.csv",
      double = "./data/input/pq/x86/double.csv"  
    ),
    arm = list(
      no = "./data/input/pq/arm/no.csv",
      single = "./data/input/pq/arm/single.csv",
      double = "./data/input/pq/arm/double.csv"
    )
  ),
  classic = list(
    x86 = list(
      no = "./data/input/classic/x86/no.csv",
      single = "./data/input/classic/x86/single.csv",
      double = "./data/input/classic/x86/double.csv"  
    ),
    arm = list(
      no = "./data/input/classic/arm/no.csv",
      single = "./data/input/classic/arm/single.csv",
      double = "./data/input/classic/arm/double.csv"
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
test <- wilcox.test(total ~ -1 + label, data = dat)

if (test$p.value < 0.05) stop("Train-test is not homogeneous")

dat %>% write_csv("./data/output/pqtls.csv")
