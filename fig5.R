####fig5c in situ CAR proportion vs mes infl
a <- read.csv("PATH/to/suppl_table_7")

v <- log(a$CART_IF_stain)
u <- log(a$Mes_infl + 1)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(insitu_CAR = v, Mes_infl = u)

library(viridis)
ggplot(uv, 
       aes(x = insitu_CAR, y = Mes_infl)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####fig5c blood CAR proportion vs mes infl
v <- log(a$X12W_postCAR_CD3)
u <- log(a$Mes_infl + 1)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(blood_CAR = v, Mes_infl = u)

library(viridis)
ggplot(uv, 
       aes(x = blood_CAR, y = Mes_infl)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()


