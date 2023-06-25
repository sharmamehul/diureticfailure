library(sjPlot)
library(sjmisc)
library(ggplot2)

storage <- list()
for(i in names(file)[-1]){
  model <- cor.test(~ logBE + get(i), data = file, method = "spearman")
  storage[[i]] <- model$p.value
}

library(ROCR)
for_ROC <- cbind(random$`All binary`, random$Homedosemultiplier)
colnames(for_ROC) <- c("outcome", "hdm")
for_ROC <- as.data.frame(for_ROC)

sample_size = floor(0.9*nrow(for_ROC))
set.seed(777)

picked = sample(seq_len(nrow(for_ROC)),size = sample_size)
train =for_ROC[picked,]
train <- as.data.frame(train)
holdout =for_ROC[-picked,]

fit <- glm(diureticdata_clean_forR_new_AS$Homedosemultiplier ~ diureticdata_clean_forR_new_AS$Creatinineatdischarge)
summary(fit)
glm.probs <- predict(fit, newdata = random, type = "response")
pred <- prediction(glm.probs, random$`All binary`)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b= 1)
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc

test <- glm(diureticdata_clean_forR$Homedosemultiplier ~ diureticdata_clean_forR$Age + diureticdata_clean_forR$Gender + 
              diureticdata_clean_forR$HRonadmission +
              diureticdata_clean_forR$EF + diureticdata_clean_forR$CAD + diureticdata_clean_forR$PriorPCICABG +
              diureticdata_clean_forR$HTN + diureticdata_clean_forR$HLD +
              diureticdata_clean_forR$PAD + diureticdata_clean_forR$AF +
              diureticdata_clean_forR$DM + diureticdata_clean_forR$Insulin +
              diureticdata_clean_forR$CKD + diureticdata_clean_forR$COPDAsthma +
              diureticdata_clean_forR$Cerebrovasculardisease + diureticdata_clean_forR$Pacemaker +
              diureticdata_clean_forR$ICD + diureticdata_clean_forR$Priorhistoryofsmoking +
              diureticdata_clean_forR$Activesmoker + diureticdata_clean_forR$Liverdisease +
              diureticdata_clean_forR$BUNatdischarge + diureticdata_clean_forR$Creatinineatdischarge +
              diureticdata_clean_forR$Homediureticdosing +
              diureticdata_clean_forR$Initialhospitaldiureticdose + diureticdata_clean_forR$Hospitallength + 
              diureticdata_clean_forR$Deltacreatinineat72hours_A +
              diureticdata_clean_forR$AKIat72hours + diureticdata_clean_forR$AKIondischarge +
              diureticdata_clean_forR$Deltacreatinineatdischarge + diureticdata_clean_forR$Baselineneutrophil +
              diureticdata_clean_forR$Baselinelymphocyte +
              diureticdata_clean_forR$`%deltaweight` + diureticdata_clean_forR$Baselinetroponin +
              diureticdata_clean_forR$SBPondischarge + diureticdata_clean_forR$Sodiumatdischarge +
              diureticdata_clean_forR$Potassiumatdischarge, family = "binomial")
summary(test)

apply(diureticdata_clean_forR[,5:6], 2, function(x) summary(glm(diureticdata_clean_forR$Homedosemultiplier ~ x)))



exp(cbind(coef(test), confint(test)))

coef <- test$coefficients[2]
pval <- anova(test)$'Pr(>F)'[1]
coef
pval



test <- glm(diuretic_data_csv_AS_forR$All ~ diuretic_data_csv_AS_forR$Initialhospitaldiureticdose, family = "binomial")
summary(test)
coef <- test$coefficients[2]
pval <- anova(test)$'Pr(>F)'[1]
coef
pval



test <- glm(diureticdata_clean_forR$`All binary` ~ diureticdata_clean_forR$Initialhospitaldiureticdose, family = "binomial")
summary(test)

test <- rlm(diureticdata_clean_forR$weight_change ~ diureticdata_clean_forR$Gender)
summary(test)
f.robftest(test)
test <- kruskal.test(diureticdata_clean_forR$weight_change ~ as.factor(diureticdata_clean_forR$Homedosemultiplier))
summary(test)

table(diureticdata_clean_forR$Homedosemultiplier, diureticdata_clean_forR$`All binary`)

ggplot(diureticdata_clean_forR, aes(as.factor(`All binary`), weight_change)) + geom_boxplot()

c <- glm(diureticdata_clean_forR$`All binary` ~ diureticdata_clean_forR$Homedosemultiplier +
           diureticdata_clean_forR$Deltacreatinineatdischarge + diureticdata_clean_forR$Baselinecreatinine +
           diureticdata_clean_forR$Deltacreatinineat72hours_A +
           diureticdata_clean_forR$AKIat72hours + diureticdata_clean_forR$AKIondischarge, family = "binomial")
summary(c)

q <- glm(diureticdata_clean_forR$Homedosemultiplier ~ diureticdata_clean_forR$`%deltaweight` + diureticdata_clean_forR$Deltacreatinineatdischarge)
summary(q)


# EF vs deltacreatinine

a <- lm(Deltacreatinineatdischarge ~ EF*Homedosemultiplier, data = diureticdata_clean_forR_new_AS)
plot_model(a, type = "int", mdrt.values = "quart")
summary(a)

data <- as.data.frame(diureticdata_clean_forR_new_AS)
data <- data[!is.na(data$EF), ]
data$EF_cat <- sapply(data$EF, function(x) {
  if(x < 40){
    y <- 1
  } else if(x < 50){
    y <- 2
  } else {
    y <- 3
  }
  return(y)
})
data$EF_cat <- as.factor(data$EF_cat)
ggplot(data[data$EF_cat!= 2, ], aes(Homedosemultiplier, Deltacreatinineatdischarge, color = EF_cat)) + geom_smooth(method = "lm", se = FALSE) + geom_point()

# can we see relationship between HDM and outcome , groupping by EF??
data.lowEF <- diureticdata_clean_forR_new_AS[diureticdata_clean_forR_new_AS$EF <= 40, ]
data.highEF <- diureticdata_clean_forR_new_AS[diureticdata_clean_forR_new_AS$EF >= 50, ]

m1 <- glm(data.lowEF$`All binary` ~ data.lowEF$HomeDoseMultiplierGroups)
m2 <- glm(data.highEF$`All binary` ~ data.highEF$HomeDoseMultiplierGroups)
summary(m1) # No , not sig
summary(m2) # no, not sig


# Other heart function indicators?
m3 <- lm(data$Deltacreatinineatdischarge ~ data$BaselineBNP * data$Homedosemultiplier)
m4 <- lm(data$Deltacreatinineatdischarge ~ data$Baselinepotassium * data$Homedosemultiplier)
m5 <- lm(Deltacreatinineatdischarge ~ EF * Homedosemultiplier + BaselineBNP + Baselinepotassium, data = diureticdata_clean_forR_new_AS)

# is weight change affected?
m6 <- lm(`%deltaweight` ~ EF * Homedosemultiplier, data = diureticdata_clean_forR_new_AS)

# Is EF*homedosemultiplier + weight chage associated with outcome?
m7 <- glm(`All binary` ~ Baselinecreatinine * Homedosemultiplier + `%deltaweight`, data = data) # did not work

# Is there subgroup of patients where relationship between HDM and deltacreatinine is different
m1.1 <- glm(`All binary` ~ Deltacreatinineatdischarge + `%deltaweight`, data = diureticdata_clean_forR_new_AS)
m2.1 <- lm(Deltacreatinineatdischarge ~ Homedosemultiplier, data = diureticdata_clean_forR_new_AS)
m1.2 <- glm(`All binary` ~ Deltacreatinineatdischarge + `%deltaweight`, data = data.highEF)
m2.2 <- lm(Deltacreatinineatdischarge ~ Homedosemultiplier, data = data.highEF)
m1.3 <- glm(`All binary` ~ Deltacreatinineatdischarge + `%deltaweight`, data = data.lowEF)
m2.3 <- lm(Deltacreatinineatdischarge ~ Homedosemultiplier, data = data.lowEF)

         