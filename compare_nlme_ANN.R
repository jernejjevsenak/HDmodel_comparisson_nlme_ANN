library(dplyr)
library(brnn)
library(lmfor)
library(nlme)
library(reshape2)
library(MLmetrics)

source("fit_nlme_nested.R")
data_PIAB <- read.csv("Norway_spruce_height2018.csv", sep = ",")

data_PIAB_plot <- data_PIAB # data frame for plot level comparison
data_PIAB_nested <- data_PIAB # data frame for nested (plot ID + grouping) comparison
data_PIAB_grouped_plot <- data_PIAB # data frame for regional level comparison
data_PIAB_regional <- data_PIAB # data frame for regional level comparison

#################################
# 1 COMPARISSON ON A PLOT LEVEL #
#################################

# 1.1 ANN plotID
ANN1 <- brnn(height ~ DBH + plotID, data = data_PIAB_plot, neurons = 3)
data_PIAB_plot$h_predicted_ANN1 <- predict(ANN1)

# 1.2 ANN plotID + BAL 
ANN2 <- brnn(height ~ DBH + plotID + BAL, data = data_PIAB_plot, neurons = 3)
data_PIAB_plot$h_predicted_ANN2 <- predict(ANN2)

# 1.3 lmfor curtis
lmfor_curtis <- ImputeHeights(data_PIAB_plot$DBH, data_PIAB_plot$height, data_PIAB_plot$plotID, 
                              modelName = "curtis", nranp = 2, varf = TRUE,
                              addResidual = FALSE, makeplot = FALSE, level = 1,
                              start=NA, bh=1.3, control=list(maxIter = 1000, msMaxIter =400),random=NA)
data_PIAB_plot$h_predicted_curtis <- lmfor_curtis$hpred

# 1.4 lmfor naslund
lmfor_naslund <- ImputeHeights(data_PIAB_plot$DBH, data_PIAB_plot$height, data_PIAB_plot$plotID, 
                               modelName = "naslund", nranp = 2, varf = TRUE,
                               addResidual = FALSE, makeplot = FALSE, level = 1,
                               start=NA, bh=1.3, control=list(maxIter = 1000, msMaxIter =400),random=NA)
data_PIAB_plot$h_predicted_naslund <- lmfor_naslund$hpred

# 1.5 lmfor michailoff
lmfor_michailoff <- ImputeHeights(data_PIAB_plot$DBH, data_PIAB_plot$height, data_PIAB_plot$plotID, 
                                  modelName = "michailoff", nranp = 2, varf = TRUE,
                                  addResidual = FALSE, makeplot = FALSE, level = 1,
                                  start=NA, bh=1.3, control=list(maxIter = 1000, msMaxIter =400),random=NA)
data_PIAB_plot$h_predicted_michailoff <- lmfor_michailoff$hpred

# 1.6 lmfor wykoff
lmfor_wykoff <- ImputeHeights(data_PIAB_plot$DBH, data_PIAB_plot$height, data_PIAB_plot$plotID, 
                              modelName = "wykoff", nranp = 2, varf = TRUE,
                              addResidual = FALSE, makeplot = FALSE, level = 1,
                              start=NA, bh=1.3, control=list(maxIter = 1000, msMaxIter =400),random=NA)
data_PIAB_plot$h_predicted_wykoff <- lmfor_wykoff$hpred

# Evaluation at the plot level
DF_plot <- melt(data_PIAB_plot, id.vars = c("plotID", "DBH", "species", "height", "BAL", "PC1_cat", "plotID_PC1_cat"))

group_by(DF_plot, variable) %>% summarise(
  n = n(),
  RMSE = RMSE(height, value),
  rRMSE = RMSE/mean(value)*100,
  RRSE = MLmetrics::RRSE(height, value),
  R2_Score = MLmetrics::R2_Score(height, value),
  bias = mean(height, na.rm = TRUE) - mean(value, na.rm = TRUE)
)

#######################################################
# 2 COMPARISSON WITH NESTED GROUPS (PLOTID + PC1_cat) #
#######################################################

# 2.1 ANN plotID_PC1_cat
ANN1 <- brnn(height ~ DBH + plotID_PC1_cat, data = data_PIAB_nested, neurons = 3)
data_PIAB_nested$h_predicted_ANN1 <- predict(ANN1)

# 2.2 ANN plotID_PC1_cat + BAL
ANN2 <- brnn(height ~ DBH + plotID_PC1_cat + BAL, data = data_PIAB_nested, neurons = 3)
data_PIAB_nested$h_predicted_ANN2 <- predict(ANN2)

# 2.3 nlme curtis
nlme_curtis <- fit_nlme_two_levels(d = data_PIAB_nested$DBH, 
                                   h = data_PIAB_nested$height, 
                                   plot = data_PIAB_nested$plotID, 
                                   plot_PCA = data_PIAB_nested$PC1_cat,
                                   modelName = "curtis")
data_PIAB_nested$h_predicted_curtis <- predict(nlme_curtis)

# 2.4 nlme naslund
nlme_naslund <- fit_nlme_two_levels(d = data_PIAB_nested$DBH, 
                                    h = data_PIAB_nested$height, 
                                    plot = data_PIAB_nested$plotID, 
                                    plot_PCA = data_PIAB_nested$PC1_cat,
                                    modelName = "naslund")
data_PIAB_nested$h_predicted_naslund <- predict(nlme_naslund)

# 2.5 nlme michailoff
nlme_michailoff <- fit_nlme_two_levels(d = data_PIAB_nested$DBH, 
                                       h = data_PIAB_nested$height, 
                                       plot = data_PIAB_nested$plotID, 
                                       plot_PCA = data_PIAB_nested$PC1_cat,
                                       modelName = "michailoff")
data_PIAB_nested$h_predicted_michailoff <- predict(nlme_michailoff)

# 2.6 nlme wykoff
nlme_wykoff <- fit_nlme_two_levels(d = data_PIAB_nested$DBH, 
                                   h = data_PIAB_nested$height, 
                                   plot = data_PIAB_nested$plotID, 
                                   plot_PCA = data_PIAB_nested$PC1_cat,
                                   modelName = "wykoff")
data_PIAB_nested$h_predicted_wykoff <- predict(nlme_wykoff)

# Evaluation at the nested group and plot effects
DF_nested <- melt(data_PIAB_nested, id.vars = c("plotID", "DBH", "species", "height", "BAL", "PC1_cat", "plotID_PC1_cat"))

group_by(DF_nested, variable) %>% summarise(
  n = n(),
  RMSE = RMSE(height, value),
  rRMSE = RMSE/mean(value)*100,
  RRSE = MLmetrics::RRSE(height, value),
  R2_Score = MLmetrics::R2_Score(height, value),
  bias = mean(height, na.rm = TRUE) - mean(value, na.rm = TRUE)
)


###########################################
# 3 COMPARISSON AT THE GROUPED PLOT LEVEL #
###########################################

# 3.1 ANN grouped plotID
ANN1 <- brnn(height ~ DBH + PC1_cat, data = data_PIAB_grouped_plot, neurons = 3)
data_PIAB_grouped_plot$h_predicted_ANN1 <- predict(ANN1)

# 3.2 ANN regionalID + BAL 
ANN2 <- brnn(height ~ DBH + PC1_cat + BAL, data = data_PIAB_grouped_plot, neurons = 3)
data_PIAB_grouped_plot$h_predicted_ANN2 <- predict(ANN2)

# 3.3 lmfor curtis
lmfor_curtis <- ImputeHeights(data_PIAB_grouped_plot$DBH, data_PIAB_grouped_plot$height, data_PIAB_grouped_plot$PC1_cat, 
                       modelName = "curtis", nranp = 2, varf = TRUE, addResidual = FALSE, makeplot=TRUE, level = 1,
                       start=NA, bh=1.3, control=list(maxIter = 1000, msmaxIter = 1000), random=NA)
data_PIAB_grouped_plot$h_predicted_lmfor_curtis <- lmfor_curtis$hpred

# 3.4 Naslund
lmfor_naslund <- ImputeHeights(data_PIAB_grouped_plot$DBH, data_PIAB_grouped_plot$height, data_PIAB_grouped_plot$PC1_cat, 
                              modelName = "naslund", nranp = 2, varf = TRUE, addResidual = FALSE, makeplot=TRUE, level = 1,
                              start=NA, bh=1.3, control=list(maxIter = 1000, msmaxIter = 1000), random=NA)
data_PIAB_grouped_plot$h_predicted_lmfor_naslund <- lmfor_naslund$hpred

# 3.5 Michailoff
lmfor_michailoff <- ImputeHeights(data_PIAB_grouped_plot$DBH, data_PIAB_grouped_plot$height, data_PIAB_grouped_plot$PC1_cat, 
                              modelName = "michailoff", nranp = 2, varf = TRUE, addResidual = FALSE, makeplot=TRUE, level = 1,
                              start=NA, bh=1.3, control=list(maxIter = 1000, msmaxIter = 1000), random=NA)
data_PIAB_grouped_plot$h_predicted_lmfor_michailoff <- lmfor_michailoff$hpred

# 3.6 wykoff
lmfor_wykoff <- ImputeHeights(data_PIAB_grouped_plot$DBH, data_PIAB_grouped_plot$height, data_PIAB_grouped_plot$PC1_cat, 
                              modelName = "wykoff", nranp = 2, varf = TRUE, addResidual = FALSE, makeplot=TRUE, level = 1,
                              start=NA, bh=1.3, control=list(maxIter = 1000, msmaxIter = 1000), random=NA)
data_PIAB_grouped_plot$h_predicted_lmfor_wykoff <- lmfor_wykoff$hpred

# Evaluation at the grouped plot level
DF_grouped_plot <- melt(data_PIAB_grouped_plot, id.vars = c("plotID", "DBH", "species", "height", "BAL", "PC1_cat", "plotID_PC1_cat"))

group_by(DF_grouped_plot, variable) %>% summarise(
  n = n(),
  RMSE = RMSE(height, value),
  rRMSE = RMSE/mean(value)*100,
  RRSE = MLmetrics::RRSE(height, value),
  R2_Score = MLmetrics::R2_Score(height, value),
  bias = mean(height, na.rm = TRUE) - mean(value, na.rm = TRUE)
  )


#######################################
# 4 COMPARISSON AT THE REGIONAL LEVEL #
#######################################

# 4.1 ANN regional
ANN1 <- brnn(height ~ DBH, data = data_PIAB_regional, neurons = 3)
data_PIAB_regional$h_predicted_ANN1 <- predict(ANN1)

# 4.2 ANN regional + BAL 
ANN2 <- brnn(height ~ DBH + BAL, data = data_PIAB_regional, neurons = 3)
data_PIAB_regional$h_predicted_ANN2 <- predict(ANN2)

# 4.3 lmfor curtis
lmfor_curtis <- startHDcurtis(data_PIAB_regional$DBH, data_PIAB_regional$height, bh=1.3) 
data_PIAB_regional$h_predicted_lmfor_curtis <- HDcurtis(d = data_PIAB_regional$DBH, a = lmfor_curtis[1], b = lmfor_curtis[2])

# 4.4 Naslund
lmfor_naslund <- startHDnaslund(data_PIAB_regional$DBH, data_PIAB_regional$height, bh=1.3) 
data_PIAB_regional$h_predicted_lmfor_naslund <- HDnaslund(d = data_PIAB_regional$DBH, a = lmfor_naslund[1], b = lmfor_naslund[2])

# 4.5 Michailoff
lmfor_michailoff <- startHDmichailoff(data_PIAB_regional$DBH, data_PIAB_regional$height, bh=1.3) 
data_PIAB_regional$h_predicted_lmfor_michailoff <- HDmichailoff(d = data_PIAB_regional$DBH, a = lmfor_michailoff[1], b = lmfor_michailoff[2])

# 4.6 wykoff
lmfor_wykoff <- startHDwykoff(data_PIAB_regional$DBH, data_PIAB_regional$height, bh=1.3)
data_PIAB_regional$h_predicted_lmfor_wykoff <- HDwykoff(d = data_PIAB_regional$DBH, a = lmfor_wykoff[1], b = lmfor_wykoff[2])

# Evaluation at the regional level
DF_regional <- melt(data_PIAB_regional, id.vars = c("plotID", "DBH", "species", "height", "BAL", "PC1_cat", "plotID_PC1_cat"))

group_by(DF_regional, variable) %>% summarise(
  n = n(),
  n = n(),
  RMSE = RMSE(height, value),
  rRMSE = RMSE/mean(value)*100,
  RRSE = MLmetrics::RRSE(height, value),
  R2_Score = MLmetrics::R2_Score(height, value),
  bias = mean(height, na.rm = TRUE) - mean(value, na.rm = TRUE)
)