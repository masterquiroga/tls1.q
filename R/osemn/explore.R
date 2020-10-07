library("tidyverse")
library("greybox")
library("ggcorrplot")
library("GGally")
library("ggiraphExtra")

dat_effects <- list(
  protocol =  c("pq", "classic"),
  arch = c("x86", "arm"),
  scheme = c("no", "single", "double"),
  speed = c("slow", "average", "fast")
)

dat <- read_csv("./data/output/pqtls.csv") %>%
  mutate(
    protocol = factor(protocol, levels = dat_effects$protocol),
    arch = factor(arch, levels = dat_effects$arch),
    scheme = factor(scheme, levels = dat_effects$scheme),
    speed = factor(speed, levels = dat_effects$speed)
  ) %>%
  as.data.frame

assocs <- dat %>% select(-iteration, -label) %>% greybox::assoc()
#pdf(file = "./data/output/pqtls_assoc.pdf", width = 400, height = 400, pointsize = 1024)
ggcorrplot(
  assocs$value, 
  type = "upper",
  colors = c("#6D9EC1", "white","#0f6b4f"), # "#E46726"),
  hc.order = TRUE, 
  show.legend = TRUE,
  #method = "circle",
  p.mat = assocs$p.value, 
  lab = TRUE,
  insig = "blank", 
  ggtheme=ggthemes::theme_tufte(base_size = 12)) %>% ggsave(
  filename = "./data/output/pqtls_assoc.pdf", plot = .,
  width=600,
  height=600,
  units="mm",
  scale=1/4,
  dpi=600,
)
dev.off()


#pdf("./data/output/pqtls_pairs.pdf", width = 400, height = 400, pointsize = 12)
(p <- dat %>% select(
    protocol,
    arch,
    scheme,
    speed,
    iteration,
    auth, 
    exchange,
    cypher,
    total
  ) %>% ggpairs(.,
    switch = "y",
    upper = list(
      continuous = function(data, mapping, ..., low = "white", high = "black") {
        ggplot(data = data, mapping = mapping) +
          geom_hex(...) +
          scale_fill_gradient(low = low, high = high)
      },
      combo = wrap("box", outlier.size = 0.005, outlier.alpha = 0.05)
    ),
    diag = list(
      #continuous = 'barDiag'
    ),
    lower = list(
      continuous = wrap("points", size = 0.005, alpha = 0.1),
      discrete = function(data, mapping, ...) {
        ggplot(data = data, mapping = mapping) +
          geom_count(...)
      }
      #combo = 'facetdensity',
      #discrete = 'ratio'
    )
  ) + ggthemes::theme_tufte()) %>% ggsave(
  filename = "./data/output/pqtls_pairs.pdf",
  plot = .,
  width=400,
  height=400,
  units="mm",
  scale=1/3,
  dpi=1200,
)
dev.off()
# dat %>% select(
#     protocol,
#     scheme,
#     iteration,
#     cypher,
#     arch,
#     auth, 
#     exchange,
#     total
#   ) %>%
#   greybox::spread(histograms = TRUE)


sink("./data/output/pqtls_detailed-summary.txt")
dat %>% group_by(protocol, arch, scheme) %>%
  summarize(
    min = min(total),
    mean = mean(total),
    median = median(total),
    max = max(total),
    sd = sd(total),
    iqr = IQR(total)
  ) %>% knitr::kable(format = 'latex', 
  caption = "Detailed Summary of Clock Times: The clock times are grouped by protocol tested, execution architecture and authentication (handshake) scheme.")
sink()

sink("./data/output/pqtls_summary.txt")
dat %>% group_by(protocol) %>%
  summarize(
    min = min(total),
    mean = mean(total),
    median = median(total),
    max = max(total),
    sd = sd(total),
    iqr = IQR(total)
  ) %>% knitr::kable(format = 'latex', 
  caption = "Summary of Clock Times: The clock times are grouped only by protocol tested regardless of execution architecture and authentication (handshake) scheme.")
sink()




