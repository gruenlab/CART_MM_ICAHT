####fig4f D90_clinical_score vs mes infl
a <- read.csv("PATH/to/suppl_table_6")
b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)
u <- log(a$Mes_infl + 1)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(D90_clinical_score = v, Mes_infl = u)

library(viridis)
ggplot(uv, 
       aes(x = Mes_infl, y = D90_clinical_score)) +
  geom_point(color = "grey30") +  
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

