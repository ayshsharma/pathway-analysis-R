# Plotting boxplot
df_long$T <- factor(df_long$T)

boxplot(Z ~ T, data = df_long,
        outline = FALSE, las = 2,  # las=2 rotates x labels
        xlab = "T", ylab = "Z",
        main = "Z distribution by T")





# For readability, plotting every 25th T (still boxplot)
keep <- as.integer(as.character(df_long$T)) %% 50 == 0
boxplot(Z ~ T, data = df_long[keep, ],
        outline = FALSE, las = 2,
        xlab = "T (every 50)", ylab = "Z")






# Next plot
library(ggplot2)
df_long$T_num <- as.integer(as.character(df_long$T))
df_25 <- subset(df_long, T_num %% 50 == 0)

ggplot(df_25, aes(x = factor(T_num), y = Z)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "T (every 50)", y = "Z") +
  theme_bw()





# Next plot
df_long2 <- transform(df_long, T_num = as.numeric(T))

ggplot(df_long2, aes(x = T_num, y = Z)) +
  geom_point(alpha = 0.15, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "T", y = "Z", title = "All points with best-fit line") +
  theme_bw()




# Next plot
library(dplyr)
df_long2 <- df_long %>%
  mutate(T_num = as.numeric(T))

df_med <- df_long2 %>%
  group_by(T_num) %>%
  summarise(Z_med = median(Z, na.rm = TRUE), .groups = "drop")

ggplot(df_long2, aes(x = factor(T_num), y = Z)) +
  geom_boxplot(outlier.shape = NA) +
  geom_smooth(
    data = df_med,
    aes(x = T_num, y = Z_med, group = 1),
    method = "lm", se = TRUE, color = "red"
  ) +
  labs(x = "T", y = "Z", title = "Z by T with linear fit to median(Z)") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())






# Next plot (Best, imo)
df_trend <- df_long %>%
  mutate(T_num = as.numeric(T)) %>%
  group_by(T_num) %>%
  summarise(Z_med = median(Z, na.rm = TRUE), .groups = "drop")

ggplot(df_trend, aes(T_num, Z_med)) +
  geom_line(alpha = 0.4) +
  geom_smooth(method = "loess", span = 0.2, se = TRUE, color = "red") +
  labs(x = "T", y = "median(Z across KEGG IDs)",
       title = "Trend of Z vs T (smoothed median)") +
  theme_bw()



# Next plot
df2 <- df_long %>% mutate(T_num = as.numeric(T))

df2 <- df2 %>% mutate(T_bin = 25 * floor((T_num - 1) / 25) + 1) # Bin T into groups of 25

df_bin <- df2 %>%
  group_by(T_bin) %>%
  summarise(Z_med = median(Z, na.rm = TRUE), .groups = "drop")

ggplot(df2, aes(factor(T_bin), Z)) +
  geom_boxplot(outlier.shape = NA) +
  geom_smooth(data = df_bin, aes(x = T_bin, y = Z_med, group = 1),
              method = "loess", span = 0.4, se = TRUE, color = "red") +
  labs(x = "T (binned by 25)", y = "Z") +
  theme_bw()







# Next (hopefully last) plot
df2 <- df_long %>% mutate(T_num = as.numeric(T))

df2 <- df2 %>% mutate(T_bin = 25 * floor((T_num - 1) / 25) + 1) # Bin T into groups of 25 (1–25, 26–50, ...)

df_bin <- df2 %>%
  group_by(T_bin) %>%
  summarise(Z_med = median(Z, na.rm = TRUE), .groups = "drop")

ggplot(df2, aes(x = T_bin, y = Z)) +
  geom_boxplot(aes(group = T_bin), outlier.shape = NA) +
  geom_smooth(
    data = df_bin,
    aes(x = T_bin, y = Z_med),
    method = "loess", span = 0.4, se = TRUE, color = "red"
  ) +
  scale_x_continuous(breaks = seq(min(df2$T_bin), max(df2$T_bin), by = 100)) +
  labs(x = "T (binned, lower bound)", y = "Z",
       title = "Boxplots by T-bin + trend of median(Z)") +
  theme_bw()