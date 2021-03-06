# PQ-TLS Exploratory Data Analyser

require("tidyverse")
require("broom")
require("dlookr")
require("RColorBrewer")
require("ggiraph")
require("ggiraphExtra")
require("plyr")
require("crayon")
require("tinytex")

# Logarithmic scale to test for how speedy are clock times
axis_labels <- list(
  '1 cycle' = 1, 
  '10 cycles' = 10, 
  '100 cycles' = 100,
  '1000 cycles' = 1000,
  '10000 cycles' = 10000,
  '100000 cycles' = 100000,
  '1000000 cycles' = 1000000,
  '10000000 cycles' = 10000000,
  '100000000 cycles' = 100000000,
  '1000000000 cycles' = 1000000000
)

datas <- list()
data <- NULL

feature_engineer <- function(data) {
  data %>% mutate(
    iteration = 1:nrow(.),
    type = if_else(total %in% boxplot.stats(total)$out, "steady", "outlier"),
    slow = if_else(total < min(boxplot.stats(total)$conf), 0, 1),
    fast = if_else(total > max(boxplot.stats(total)$conf), 0, 1))
}

eda <- function(file, name, estimate_auth = FALSE) {
  dir <- getwd()
  data <- read_csv(file, 
    col_names = c(
      "auth",
      "exchange",
      "cypher",
      "total"
    )) %>% 
    mutate(
      iteration = 1:nrow(.),
      type = if_else(total %in% boxplot.stats(total)$out, "steady", "outlier"),
      slow = if_else(total < min(boxplot.stats(total)$conf), 0, 1),
      fast = if_else(total > max(boxplot.stats(total)$conf), 0, 1))
  if (estimate_auth) {
    data <- data %>% mutate(auth = total - exchange - cypher)
  }

  data %>% diagnose_report(
    target = total, 
    output_format = "html",
    output_file = paste0(name,"_diagnosis.html"),
    output_dir = "./data/output")

  data %>% select(-type, -slow, -fast) %>% eda_report(
    target = total, 
    output_format = "html",
    output_file = paste0(name,"_eda.html"),
    output_dir = "./data/output")
  data %>% select(-iteration, -type, -slow, -fast) %>% transformation_report(
    target = total,
    output_format = "html",
    output_file = paste0(name,"_transformation.html"),
    output_dir = "./data/output")

  # Transform the dataset ----
  data <- data %>% mutate(total_lvl = cut(exchange, 
    breaks = c(0, boxplot.stats(exchange)$conf, Inf),
    right = FALSE))
  datas[[name]] <- data



  try({
    setwd(dir)
    setwd("./data/output")
    before_after <- data %>%
      mutate(elapsed = cumsum(total)) %>%
      arrange(elapsed) %>% 
      mutate(before = as.numeric(elapsed - lag(elapsed)),
            after =  as.numeric(lead(elapsed) - elapsed)) %>%
      filter(!is.na(before) & !is.na(after))
      ggplot(before_after, aes(x = before, y = after)) +
          geom_hex() +
          scale_fill_gradientn(colours = rev(brewer.pal(5, "Spectral"))) +
          scale_y_continuous( minor_breaks = NULL) +
          scale_x_continuous( minor_breaks = NULL) +
          labs( x = "Cycles before TLS", y = "Cycles after TLS") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(filename = paste0(name,"_time_plot.png"),
        width = 3,
        height = 3)
  })

  try({
    setwd(dir)
    setwd("./data/output")
    linr <- data %>% lm(total ~ auth*cypher*total_lvl,
        data = .) 
    sink(paste0(name,"_linr_summary.txt"))
    linr %>% summary %>% print
    sink()

    linr %>% ggPredict(interactive = FALSE)
    ggsave(filename = paste0(name,"_linr_a.png"),
      width = 3,
      height = 3)
  })

  try({
    setwd(dir)
    setwd("./data/output")
    logr <- data %>% glm(slow ~ auth*cypher*total_lvl,
      data = .,
      family = "binomial")
    sink(paste0(name,"_logr_summary.txt"))
    logr %>% summary %>% print
    sink()

    logr %>% ggPredict(interactive = FALSE)
    ggsave(filename = paste0(name,"_logr.png"),
      width = 3,
      height = 3)
  })

  setwd(dir)
}



eda("./data/input/pq/x86/no.csv", "pq_x86_no_auth", TRUE)
eda("./data/input/pq/x86/single.csv", "pq_x86_single_auth", FALSE)
eda("./data/input/pq/x86/double.csv", "pq_x86_double_auth", FALSE)
eda("./data/input/pq/arm/no.csv", "pq_arm_no_auth", TRUE)
eda("./data/input/pq/arm/single.csv", "pq_arm_single_auth", FALSE)
eda("./data/input/pq/arm/double.csv", "pq_arm_double_auth", FALSE)


eda("./data/input/classic/x86/no.csv", "classic_x86_no_auth", TRUE)
eda("./data/input/classic/x86/single.csv", "classic_x86_single_auth", FALSE)
eda("./data/input/classic/x86/double.csv", "classic_x86_double_auth", FALSE)
eda("./data/input/classic/arm/no.csv", "classic_arm_no_auth", TRUE)
eda("./data/input/classic/arm/single.csv", "classic_arm_single_auth", FALSE)
eda("./data/input/classic/arm/double.csv", "classic_arm_double_auth", FALSE)
