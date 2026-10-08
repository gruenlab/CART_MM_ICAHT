###ED3a clinical score vs CRS grade
a <- read.csv("PATH/to/suppl_table_7")
b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)

library(ggplot2)
library(dplyr)

df <- df %>%
  mutate(CRS_grade = as.factor(CRS_grade))

df <- data.frame(clinical_score = v, CRS_grade = a$CRS_grade, condition = a$condition)

kw <- kruskal.test(clinical_score ~ CRS_grade, data = df)
p_label <- paste0("p = ", signif(kw$p.value, 3))

ggplot(df, aes(x = factor(CRS_grade), y = clinical_score)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(
    aes(color = condition),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "CRS grade",
    y = "clinical score",
    color = "Condition",
    subtitle = p_label
  ) +
  theme_classic()

### mes infl vs CRS grade

df <- df %>%
  mutate(CRS_grade = as.factor(CRS_grade))

df <- data.frame(Mes_infl = a$Mes_infl, CRS_grade = a$CRS_grade, condition = a$condition)

kw <- kruskal.test(Mes_infl ~ CRS_grade, data = df)
p_label <- paste0("p = ", signif(kw$p.value, 3))

ggplot(df, aes(x = factor(CRS_grade), y = Mes_infl)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(
    aes(color = condition),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "CRS grade",
    y = "clinical score",
    color = "Condition",
    subtitle = p_label
  ) +
  theme_classic()

### ED3b treatment instensity vs conditions
df <- data.frame(Mes_infl = a$Mes_infl, treatment_intensity = a$treatment, condition = a$condition)

kw <- t.test(treatment_intensity ~ condition, data = df)

p_label <- paste0("t-test p = ", signif(kw$p.value, 3))

ggplot(df, aes(x = condition , y =  treatment_intensity)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(
    aes(color = condition),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    
    y = "treatment_intensity",
    color = "Condition",subtitle = p_label
  ) +
  theme_classic()

# ED3b treatment instensity vs mes infl

kw <- kruskal.test(Mes_infl ~ treatment_intensity, data = df)
p_label <- paste0("p = ", signif(kw$p.value, 3))

ggplot(df, aes(x = factor(treatment_intensity), y = Mes_infl)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(
    aes(color = condition),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  labs(
    x = "treatment_intensity",
    y = "Mes_infl score",
    color = "Condition",
    subtitle = p_label
  ) +
  theme_classic()

### ED3c baseline vs D90 clinical scores
z <- read.csv("PATH/to/suppl_table_6",header=TRUE)
y <- z[,c("Hb", "PLT",  "ANC")]
yn <- t(t(y)/colSums(y))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

yng <- apply(yn,1,gmean)


y <- z[,c("Hb_baseline", "PLT_baseline",  "ANC_baseline")]
yn <- t(t(y)/colSums(y))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(yn,1,gmean)

v <- log(bng)
u <- log(yng)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(baseline_clinical_score = v, D90_clinical_score = u)

ggplot(uv, aes(x = D90_clinical_score, y = baseline_clinical_score)) +
  geom_point() +
  theme_classic()+
  geom_smooth(method=lm , color="red", fill="#69b3a2") +
  sm_statCorr(color="red", fill="#69b3a2", text_size = 5) +sm_classic()

####ED3d age vs basline_clinical_score
u <- z$age
v <- log(bng)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data.frame(basline_clinical_score = v, age = u)


ggplot(uv, 
       aes(x = age, y = basline_clinical_score)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "red", fill = "red") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "red", fill = "red", 
              text_size = 5) +
  sm_classic()

####ED3d age vs D90_clinical_score
v <- log(yng)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data.frame(D90_clinical_score = v, age = u)


ggplot(uv, 
       aes(x = age, y = D90_clinical_score)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "red", fill = "red") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "red", fill = "red", 
              text_size = 5) +
  sm_classic()

####ED3d age vs mes_infl_score
v <- log(z$mes_infl_score)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data.frame(mes_infl_score = v, age = u)


ggplot(uv, 
       aes(x = age, y = mes_infl_score)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "red", fill = "red") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "red", fill = "red", 
              text_size = 5) +
  sm_classic()

####ED3e baseline IL6 level vs D90_clinical_score
b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)
u <- log(a$IL6_pre)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(D90_clinical_score = v, IL6_baseline = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_baseline, y = D90_clinical_score)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####ED3e baseline IL6 level vs mes_infl_score

v <- log(a$Mes_infl + 1)
u <- log(a$IL6_pre)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(Mes_infl = v, IL6_baseline = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_baseline, y = Mes_infl)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####ED3e baseline IL6 level vs HSC_sen
v <- log(a$HSC_sen + 1)
u <- log(a$IL6_pre)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(HSC_sen = v, IL6_baseline = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_baseline, y = HSC_sen)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####ED3f peak IL6 level vs clinical_score
v <- log(bng)
u <- log(a$IL6_peak)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(clinical_score = v, IL6_peak = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_peak, y = clinical_score)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()


####ED3f peak IL6 level vs mes_infl
v <- log(a$Mes_infl + 1)
u <- log(a$IL6_peak)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(mes_infl = v, IL6_peak = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_peak, y = mes_infl)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()


####ED3f peak IL6 level vs HSC_sen
v <- log(a$HSC_sen + 1)
u <- log(a$IL6_peak)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(HSC_sen = v, IL6_peak = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_peak, y = HSC_sen)) +
  geom_point(color = "grey30") +   # optional: change point color
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  # new line + shade colors
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()


