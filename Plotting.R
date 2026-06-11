# Creating a dataframe with all the compounds and their values for easy plotting and access
df_long <- do.call(rbind, lapply(names(Z_T_relationship), function(k) {
  x <- Z_T_relationship[[k]]
  data.frame(kegg_id = k, T = x$T, Z = x$Z)
}))


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



library(dplyr)

df_400T <- df_long %>%
  mutate(T_num = as.numeric(T)) %>%
  filter(T_num %in% sort(unique(T_num), decreasing = FALSE)[1:400]) %>%
  group_by(T_num) %>%
  summarise(Z_med = mean(Z, na.rm = TRUE), .groups = "drop")



ggplot(df_400T, aes(T_num, Z_med)) +
  geom_line(alpha = 0.4) +
  geom_smooth(method = "loess", span = 0.2, se = TRUE, color = "red") +
  labs(x = "T", y = "mean(Z across KEGG IDs)",
       title = "Trend of Z vs T (smoothed mean)") +
  theme_bw()





# Time for some linear regression!
library(dplyr)

# 1) mean Z at each T (across all kegg_id), restricted to T in [0, 200]
df_mean <- df_long %>%
  mutate(T_num = as.numeric(T)) %>%
  filter(T_num >= 0, T_num <= 400) %>%
  group_by(T_num) %>%
  summarise(meanZ = mean(Z, na.rm = TRUE),
            medianZ = median(Z, na.rm = TRUE),
            n = dplyr::n(),           # optional: how many points went into the mean
            .groups = "drop") %>%
  arrange(T_num)

# 2) linear regression: meanZ = b0 + b1*T
fit_mean <- lm(meanZ ~ T_num, data = df_mean)
fit_median <- lm(medianZ ~ T_num, data = df_mean)
summary(fit_mean)
summary(fit_median)

anova(fit_mean, fit_median)



# 3) print the fitted function
b0_mean <- coef(fit_mean)[["(Intercept)"]]
b1_mean <- coef(fit_mean)[["T_num"]]
cat(sprintf("mean(Z) = %.6f + %.6f * T\n", b0_mean, b1_mean))

b0_median <- coef(fit_median)[["(Intercept)"]]
b1_median <- coef(fit_median)[["T_num"]]
cat(sprintf("median(Z) = %.6f + %.6f * T\n", b0_median, b1_median))

# Plotting residuals
plot(df_mean$T_num, resid(fit_mean),
     xlab = "T", ylab = "Residual",
     main = "Residuals vs T (0–400)")
abline(h = 0, col = "red")





df_mean_0_400 <- df_mean %>% filter(T_num >= 0, T_num <= 400)

fit_0_400 <- lm(meanZ ~ T_num, data = df_mean_0_400)

fit_no0 <- lm(meanZ ~ T_num, data = df_mean_0_400 %>% filter(T_num != 0))

summary(fit_0_400)
summary(fit_no0)
# T = 0 seems to be an extremely weird outlier (perhaps an issue with the API/my algorithms/something else)

# Regardless, it seems to be wise to exclude T=0 when performing the linear regression
fit_1_400 <- lm(meanZ ~ T_num, data = df_mean_0_400 %>% filter(T_num >= 1))
df_mean_1_400 <- df_mean_0_400[df_mean_0_400$T_num>=1, ]
plot(df_mean_1_400$T_num, resid(fit_1_400),
     xlab = "T", ylab = "Residual",
     main = "Residuals vs T (0–100)")
abline(h = 0, col = "red")
summary(fit_1_400)

# Trying a quadratic fit, just in case it is better
model_quad <- lm(meanZ ~ poly(T_num, 2, raw = TRUE), 
                 data = df_mean_0_400 %>% filter(T_num >= 1))
summary(model_quad)


anova(fit_1_400, model_quad) # Comparing the linear fit to the quadratic one

ggplot(df_mean_0_400 %>% filter(T_num >= 1),
       aes(T_num, meanZ)) +
  geom_point() +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 2, raw = TRUE),
              se = FALSE)
par(mfrow = c(1,1))
plot(model_quad)


model_cubic <- lm(meanZ ~ poly(T_num, 3, raw = TRUE),
                  data = df_mean_0_400 %>% filter(T_num >= 1))

summary(model_cubic)
anova(fit_1_400, model_quad, model_cubic)

par(mfrow = c(2,2))
plot(model_cubic)




model_quartic <- lm(meanZ ~ poly(T_num, 4, raw = TRUE),
                    data = df_mean_0_400 %>% filter(T_num >= 1))

summary(model_quartic)
anova(fit_1_400, model_quad, model_cubic, model_quartic)

par(mfrow = c(2,2))
plot(model_quartic)

# library(mgcv)
# 
# model_gam <- gam(meanZ ~ s(T_num),
#                  data = df_mean_0_400 %>% filter(T_num >= 1))
# 
# summary(model_gam)
# 
# par(mfrow = c(2,2))
# plot(model_gam)
# 
# gam.check(model_gam)
# 
# AIC(fit_1_400, model_quad, model_cubic,
#     model_quartic, model_gam)





