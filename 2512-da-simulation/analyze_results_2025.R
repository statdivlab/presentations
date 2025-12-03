library(tidyverse)

# load in inference results from smaller models
files <- list.files("new_t1e_res", full.names = T)
t1e_res1 <- readRDS(files[1])
t1e_res <- t1e_res1
for (file in files) {
  t1e_batch <- readRDS(file)
  t1e_res <- rbind(t1e_res, t1e_batch)
}

# get data into needed format 
plot_res <- t1e_res %>%
  mutate(method = ifelse(method == "linda", "LinDA", method),
         method = ifelse(method == "maaslin", "MaAsLin3", method),
         method = ifelse(method == "CLR t-test", "CLR linear\nmodel", method)) %>% 
  filter(method != "ANCOM-BC2 with SA") %>% 
  mutate(method = factor(method, 
                         levels = 
                           c("ANCOM-BC2", "radEmu", "DESeq2", "MaAsLin3", 
                             "CLR linear\nmodel", "ALDEx2", "LinDA")))
  

# drop ANCOM-BC2 with sensitivity analysis
# t1e plot
plot_res %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_qq(distribution = stats::qunif,
          aes(sample = p,
              color = method, 
          ), geom="line",
          linewidth = 0.5) +
  labs(color = "Test") +
  xlab("Theoretical p-value quantiles") +
  ylab("Empirical p-value quantiles") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 60, size=6,hjust= 1),
        axis.text.y = element_text(size = 6)) +
  scale_x_sqrt(breaks=c(0,0.01,0.05,0.25,0.5,1)) +
  scale_y_sqrt(breaks = c(0,0.01,0.05,0.25,0.5,1)) +
  coord_equal() +
  guides(color = guide_legend(position = "right", nrow = 7)) + 
  NULL
ggsave("t1e_figure_updated.png", width = 8, height = 6)

t1e_res %>%
  filter(method != "ANCOM-BC2 with SA") %>% 
  group_by(method) %>%
  summarise(prop_1 = mean(p <= 0.01, na.rm = TRUE),
            prop_5 = mean(p <= 0.05, na.rm = TRUE),
            prop_10 = mean(p <= 0.1, na.rm = TRUE))

# MSE's
plot_res %>%
  ggplot(aes(x = method, y = mse, color = method)) + 
  geom_boxplot() +
  ylim(c(0, 2.5)) + 
  labs(y = "Mean squared error", x = "Method") + 
  ggtitle("Error in estimates across methods") + 
  theme_bw(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = "none")
ggsave("mse_figure_updated.png", width = 8, height = 6)

# distribution of estimate of tested parameter 
plot_res %>%
  ggplot(aes(x = method, y = est, color = method)) + 
  geom_boxplot() +
  labs(y = "Estimate for null category", x = "Method") + 
  ggtitle("Null estimates across methods") + 
  theme_bw(base_size = 16) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = "none") + 
  geom_hline(aes(yintercept = 0), linetype = "dotted") 
ggsave("null_error_figure_updated.png", width = 8, height = 6)
