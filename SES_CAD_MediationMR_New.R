
###### X1 education ######
GWAS_education <- fread("Input/GWAS_education.csv")
GWAS_education <- as.data.frame(GWAS_education)
str(GWAS_education)
GWAS_education$SS <- as.numeric(GWAS_education$SS)
sum(is.na(GWAS_education$SNP))

###### X2 occupation ######
GWAS_occupation <- fread("Input/GWAS_occupation_New.csv")
GWAS_occupation <- as.data.frame(GWAS_occupation)
str(GWAS_occupation)
GWAS_occupation <- GWAS_occupation[,-1]
GWAS_occupation$N <- as.numeric(GWAS_occupation$N)

###### X3 income ######
GWAS_income <- fread("Input/GWAS_income.csv")
GWAS_income <- as.data.frame(GWAS_income)
str(GWAS_income)
GWAS_income <- GWAS_income[,-1]
GWAS_income$SS <- 397751
sum(is.na(GWAS_income$SNP))

###### X4 regional deprivation(TDI) ######
GWAS_TDI <- fread("Input/GWAS_TDI.csv")
GWAS_TDI <- as.data.frame(GWAS_TDI)
GWAS_TDI$SS <- 462464
str(GWAS_TDI)

###### Y CAD ######
GWAS_CAD <- fread("Input/GWAS_CAD.csv")
GWAS_CAD <- as.data.frame(GWAS_CAD)
str(GWAS_CAD)
GWAS_CAD$SS <- as.numeric(GWAS_CAD$SS)

###### M1.1 Diet ######
###### _____carbohydrate ######
GWAS_carbohydrate <- fread("Input/GWAS_carbohydrateNew.tsv.gz")
GWAS_carbohydrate <- as.data.frame(GWAS_carbohydrate)
str(GWAS_carbohydrate)
colnames(GWAS_carbohydrate) <- c("SNP","CHR","BP","other allele","effect allele","Freq","BETA","SE","P")
GWAS_carbohydrate$N <- 282271

###### _____protein ######
GWAS_protein <- fread("Input/GWAS_proteinNew.tsv.gz")
GWAS_protein <- as.data.frame(GWAS_protein)
str(GWAS_protein)
colnames(GWAS_protein) <- c("SNP","CHR","BP","other allele","effect allele","Freq","BETA","SE","P")
GWAS_protein$N <- 282271

###### _____fat ######
GWAS_fat <- fread("Input/GWAS_fatNew.tsv.gz")
GWAS_fat <- as.data.frame(GWAS_fat)
str(GWAS_fat)
colnames(GWAS_fat) <- c("SNP","CHR","BP","other allele","effect allele","Freq","BETA","SE","P")
GWAS_fat$N <- 282271

###### _____FreshFruit ######
GWAS_FreshFruit <- fread("Input/GWAS_FreshFruit.tsv.gz")
GWAS_FreshFruit <- as.data.frame(GWAS_FreshFruit)
str(GWAS_FreshFruit)
names(GWAS_FreshFruit)[1] <- "SNP"
GWAS_FreshFruit$SS <- 447401

###### _____DriedFruit ######
GWAS_DriedFruit <- fread("Input/GWAS_DriedFruit.tsv.gz")
GWAS_DriedFruit <- as.data.frame(GWAS_DriedFruit)
str(GWAS_DriedFruit)
names(GWAS_DriedFruit)[1] <- "SNP"
GWAS_DriedFruit$SS <- 444741

###### _____RawVegetables ######
GWAS_RawVegetables <- fread("Input/GWAS_RawVegetables.tsv.gz")
GWAS_RawVegetables <- as.data.frame(GWAS_RawVegetables)
str(GWAS_RawVegetables)
names(GWAS_RawVegetables)[1] <- "SNP"
GWAS_RawVegetables$SS <- 443633

###### _____CookedVegetables ######
GWAS_CookedVegetables <- fread("Input/GWAS_CookedVegetables.tsv.gz")
GWAS_CookedVegetables <- as.data.frame(GWAS_CookedVegetables)
str(GWAS_CookedVegetables)
names(GWAS_CookedVegetables)[1] <- "SNP"
GWAS_CookedVegetables$SS <- 444190

###### _____salt ######
GWAS_salt <- fread("Input/GWAS_salt.tsv.gz")
GWAS_salt <- as.data.frame(GWAS_salt)
str(GWAS_salt)
names(GWAS_salt)[1] <- "SNP"
GWAS_salt$SS <- 448890

###### _____processed ######
GWAS_processed <- fread("Input/GWAS_processed.tsv.gz")
GWAS_processed <- as.data.frame(GWAS_processed)
str(GWAS_processed)
names(GWAS_processed)[1] <- "SNP"
GWAS_processed$SS <- 448303

###### _____OilyFish ######
GWAS_OilyFish <- fread("Input/GWAS_OilyFish.tsv.gz")
GWAS_OilyFish <- as.data.frame(GWAS_OilyFish)
str(GWAS_OilyFish)
names(GWAS_OilyFish)[1] <- "SNP"
GWAS_OilyFish$SS <- 446854


###### _____wholegrain ######
GWAS_wholegrain <- fread("Input/GWAS_wholegrain.tsv.gz")
GWAS_wholegrain <- as.data.frame(GWAS_wholegrain)
str(GWAS_wholegrain)
names(GWAS_wholegrain)[1] <- "SNP"
GWAS_wholegrain$SS <- 434087

###### M1.2 Physical activity and sedentary behavior ######
###### _____DeviceOverallActivity ######
GWAS_DeviceOverallActivity <- fread("Input/GWAS_DeviceOverallActivity.csv.gz")
GWAS_DeviceOverallActivity <- as.data.frame(GWAS_DeviceOverallActivity)
str(GWAS_DeviceOverallActivity) # P_BOLT_LMM_INF
GWAS_DeviceOverallActivity$N <- 91105
x1 <- c("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","BETA","SE","P_BOLT_LMM_INF","N")
GWAS_DeviceOverallActivity <- GWAS_DeviceOverallActivity[x1]

###### _____PA ######
GWAS_PA <- fread("Input/GWAS_PAnew.tsv.gz")
GWAS_PA <- as.data.frame(GWAS_PA)
str(GWAS_PA)
head(GWAS_PA)
colnames(GWAS_PA) <- c("SNP","CHR","BP","other allele","effect allele","Freq","BETA","SE","P","N")
sum(is.na(GWAS_PA$SNP))

###### _____DeviceModerate ######
GWAS_DeviceModerate <- fread("Input/GWAS_DeviceModerate.csv.gz")
GWAS_DeviceModerate <- as.data.frame(GWAS_DeviceModerate)
str(GWAS_DeviceModerate)
GWAS_DeviceModerate$N <- 91105
x1 <- c("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","BETA","SE","P_BOLT_LMM_INF","N")
GWAS_DeviceModerate <- GWAS_DeviceModerate[x1]

###### _____DeviceSedentary ######
GWAS_DeviceSedentary <- fread("Input/GWAS_DeviceSedentary.csv.gz")
GWAS_DeviceSedentary <- as.data.frame(GWAS_DeviceSedentary)
str(GWAS_DeviceSedentary)
GWAS_DeviceSedentary$N <- 91105
x1 <- c("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","BETA","SE","P_BOLT_LMM_INF","N")
GWAS_DeviceSedentary <- GWAS_DeviceSedentary[x1]

###### _____LST ######
GWAS_LST <- fread("Input/GWAS_LSTnew.tsv.gz")
GWAS_LST <- as.data.frame(GWAS_LST)
str(GWAS_LST)
colnames(GWAS_LST) <- c("SNP","CHR","BP","other allele","effect allele","Freq","BETA","SE","P","N")
sum(is.na(GWAS_LST$SNP))

###### _____television ######
GWAS_television <- fread("Input/GWAS_televisionnew.tsv.gz")
GWAS_television <- as.data.frame(GWAS_television)
str(GWAS_television)
GWAS_television$N <- as.numeric(GWAS_television$N)
colnames(GWAS_television) <- c("SNP","CHR","BP","other allele","effect allele","Freq","BETA","SE","P","N")

###### _____computer ######
GWAS_computer <- fread("Input/Computer.all.gz")
GWAS_computer <- as.data.frame(GWAS_computer)
str(GWAS_computer)
sum(is.na(GWAS_computer$SNP)) # 0
x1 <- c("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","BETA","SE","P_BOLT_LMM_INF")
GWAS_computer <- GWAS_computer[x1]
GWAS_computer$N <- 408815

###### _____driving ######
GWAS_driving <- fread("Input/Driving.all.gz")
GWAS_driving <- as.data.frame(GWAS_driving)
str(GWAS_driving)
sum(is.na(GWAS_driving$SNP)) # 0
x1 <- c("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","BETA","SE","P_BOLT_LMM_INF")
GWAS_driving <- GWAS_driving[x1]
GWAS_driving$N <- 408815

###### M1.3 Smoking ######
###### _____Lifetime smoking index ######
GWAS_LifetimeSmoking <- fread("Input/GWAS_LifetimeSmoking.txt")
GWAS_LifetimeSmoking <- as.data.frame(GWAS_LifetimeSmoking)
str(GWAS_LifetimeSmoking)
sum(is.na(GWAS_LifetimeSmoking$SNP)) # 0
GWAS_LifetimeSmoking$SS <- 462690

###### M1.4 Sleep health ######
###### _____SleepDuration ######
GWAS_SleepDuration <- fread("Input/GWAS_SleepDuration.txt.gz")
GWAS_SleepDuration <- as.data.frame(GWAS_SleepDuration)
str(GWAS_SleepDuration)
sum(is.na(GWAS_SleepDuration$SNP)) # 0
GWAS_SleepDuration$N <- as.numeric(GWAS_SleepDuration$N)
range(GWAS_SleepDuration$OR) # [1] -0.6111  0.5308
names(GWAS_SleepDuration)[8] <- "BETA"

###### _____DeviceSleepDuration ######
GWAS_DeviceSleepDuration <- fread("Input/GWAS_DeviceSleepDuration.csv.gz")
GWAS_DeviceSleepDuration <- as.data.frame(GWAS_DeviceSleepDuration)
str(GWAS_DeviceSleepDuration)
sum(is.na(GWAS_DeviceSleepDuration$SNP)) # 0
GWAS_DeviceSleepDuration$N <- 91105
x1 <- c("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","BETA","SE","P_BOLT_LMM_INF","N")
GWAS_DeviceSleepDuration <- GWAS_DeviceSleepDuration[x1]

###### _____Insomnia ######
GWAS_Insomnia <- fread("Input/GWAS_Insomnia.txt.gz")
GWAS_Insomnia <- as.data.frame(GWAS_Insomnia)
str(GWAS_Insomnia)
sum(is.na(GWAS_Insomnia$SNP)) # 0
GWAS_Insomnia$N <- as.numeric(GWAS_Insomnia$N)
range(GWAS_Insomnia$OR) #[1] 0.2818 2.5620
GWAS_Insomnia$BETA <- log(GWAS_Insomnia$OR)

###### M2.1 WBS ######
GWAS_WBS <- fread("Input/GWAS_WBS.txt.gz")
GWAS_WBS <- as.data.frame(GWAS_WBS)
str(GWAS_WBS)
GWAS_WBS <- GWAS_WBS[,-1]
names(GWAS_WBS)[1] <- "SNP"
GWAS_WBS$N <- as.numeric(GWAS_WBS$N)
sum(is.na(GWAS_WBS$SNP)) # 0

###### M2.2 MentalProblems ######
GWAS_MentalProblems <- fread("Input/GWAS_MentalProblems.csv")
GWAS_MentalProblems <- as.data.frame(GWAS_MentalProblems)
str(GWAS_MentalProblems)
GWAS_MentalProblems <- GWAS_MentalProblems[,-1]
GWAS_MentalProblems$SS <- 322580

###### M2.3 depression ######
GWAS_depression <- fread("Input/GWAS_depression.txt")
GWAS_depression <- as.data.frame(GWAS_depression)
str(GWAS_depression)
GWAS_depression$SampleSize <- 500199
names(GWAS_depression)[1] <- "SNP"
sum(is.na(GWAS_depression$SNP)) # 0
head(GWAS_depression)

###### M3.1 BMI ######
GWAS_BMI <- fread("Input/GWAS_BMI.csv")
GWAS_BMI <- as.data.frame(GWAS_BMI)
str(GWAS_BMI)
GWAS_BMI <- GWAS_BMI[,-1]
sum(is.na(GWAS_BMI$SNP)) # 0

###### M3.2 lipid ######
GWAS_lipid <- fread("Input/GWAS_lipidnew.tsv.gz")
GWAS_lipid <- as.data.frame(GWAS_lipid)
str(GWAS_lipid)
head(GWAS_lipid)
colnames(GWAS_lipid) <- c("SNP","CHR","BP","other allele","effect allele","N","Freq","BETA","SE","P")
sum(is.na(GWAS_lipid$SNP))

###### M3.3 glucose ######
GWAS_glucose <- fread("Input/GWAS_glucosenew.tsv.gz")
GWAS_glucose <- as.data.frame(GWAS_glucose)
str(GWAS_glucose)
head(GWAS_glucose)
colnames(GWAS_glucose) <- c("SNP","CHR","BP","other allele","effect allele","Freq","BETA","SE","P","N")
sum(is.na(GWAS_glucose$SNP))

###### M3.4 SBP, DBP ######
GWAS_SBP <- fread("Input/GWAS_SBP.csv")
GWAS_SBP <- as.data.frame(GWAS_SBP)
str(GWAS_SBP)
GWAS_SBP <- GWAS_SBP[,-1]
sum(is.na(GWAS_SBP$SNP))# 0

GWAS_DBP <- fread("Input/GWAS_DBP.csv")
GWAS_DBP <- as.data.frame(GWAS_DBP)
str(GWAS_DBP)
GWAS_DBP <- GWAS_DBP[,-1]
sum(is.na(GWAS_DBP$SNP))# 0


###### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ######
###### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ######


###### 1.1.1 education & CAD ######
###### ______education on CAD → total effect of education on CAD【total causal effect of X on Y（β）】######
IV_education <- read.xlsx("Input/IV_education.xlsx")
str(IV_education)
IV_education <- format_data(IV_education,
                            type = "exposure", 
                            snp_col = "SNP", 
                            beta_col = "Effect.size",
                            se_col = "SE",
                            effect_allele_col = "Allele.1",
                            other_allele_col = "Allele2",
                            eaf_col = "Frequency.Allele.1",
                            pval_col = "P-value",
                            samplesize_col = "N")
str(IV_education)
IV_education$exposure <- "educational attainment"
IV_education$samplesize.exposure <- 766345

IV_education$r2 <- (2 * (IV_education$beta.exposure^2) * IV_education$eaf.exposure * (1 - IV_education$eaf.exposure)) /
  (2 * (IV_education$beta.exposure^2) * IV_education$eaf.exposure * (1 - IV_education$eaf.exposure) +
     2 * IV_education$samplesize.exposure * IV_education$eaf.exposure * (1 - IV_education$eaf.exposure) * (IV_education$se.exposure^2))

IV_education$F <- (IV_education$r2 * (IV_education$samplesize.exposure - 2))/(1 - IV_education$r2)
range(IV_education$F) #[1]  28.22258 373.77680
IV_education_meanF <- mean(IV_education$F) #[1] 46.63938

# write.xlsx(IV_education, file = "Output/UVMR_IV/IV_education.xlsx")

out_education_ON_CAD <- extract_outcome_data(
  snps = IV_education$SNP,
  outcomes = 'ieu-a-7')

UVdat_education_ON_CAD <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_CAD
)
range(UVdat_education_ON_CAD$F)#[1]  28.22258 373.77680
str(UVdat_education_ON_CAD)
range(UVdat_education_ON_CAD$pval.outcome) #[1] 3.91003e-06 9.99182e-01
range(UVdat_education_ON_CAD$F)
write.xlsx(UVdat_education_ON_CAD, file = "Output/UVMR_IV/UVdat_education_ON_CAD.xlsx")

result_education_ON_CAD <- mr(UVdat_education_ON_CAD)


scatter_plot_education_ON_CAD <- mr_scatter_plot(result_education_ON_CAD, UVdat_education_ON_CAD)
scatter_plot_education_ON_CAD[[1]]
ggsave(scatter_plot_education_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_CAD.pdf", width=7, height=7)


pleiotropy_education_ON_CAD <- mr_pleiotropy_test(UVdat_education_ON_CAD)
str(pleiotropy_education_ON_CAD)
pleiotropy_education_ON_CAD <- pleiotropy_education_ON_CAD[,-c(1,2)]
pleiotropy_education_ON_CAD$outcome <- "coronary heart disease"


heterogeneity_education_ON_CAD <- mr_heterogeneity(UVdat_education_ON_CAD)
heterogeneity_education_ON_CAD <- heterogeneity_education_ON_CAD[,-c(1,2)]
heterogeneity_education_ON_CAD$outcome <- "coronary heart disease"


singleSNP_education_ON_CAD <- mr_singlesnp(UVdat_education_ON_CAD)
funnel_plot_education_ON_CAD <- mr_funnel_plot(singleSNP_education_ON_CAD)
funnel_plot_education_ON_CAD[[1]]
ggsave(funnel_plot_education_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_CAD.pdf", width=7, height=7)


loo_education_ON_CAD <- mr_leaveoneout(UVdat_education_ON_CAD)
loo_plot_education_ON_CAD <- mr_leaveoneout_plot(loo_education_ON_CAD)
loo_plot_education_ON_CAD[[1]]
ggsave(loo_plot_education_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_education_ON_CAD <- mr_forest_plot(singleSNP_education_ON_CAD)
singleSNP_plot_education_ON_CAD[[1]]
ggsave(singleSNP_plot_education_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_CAD.pdf", width=7, height=21)

PRESSO_education_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_CAD, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_CAD <- generate_odds_ratios(result_education_ON_CAD)
result_education_ON_CAD <- result_education_ON_CAD[,-c(1,2)]
result_education_ON_CAD$outcome <- "coronary heart disease"


###### ______CAD on education ######
IV_CAD <- format_data(GWAS_CAD,
                      type = "exposure", 
                      snp_col = "SNP", 
                      beta_col = "ES",
                      se_col = "SE",
                      effect_allele_col = "ALT",
                      other_allele_col = "REF",
                      eaf_col = "AF",
                      pval_col = "P",
                      samplesize_col = "SS")
str(IV_CAD)
IV_CAD$exposure <- "coronary heart disease"

IV_CAD <- IV_CAD[IV_CAD$pval.exposure < 5e-8,]
IV_CAD <- clump_data(IV_CAD)

IV_CAD$r2 <- (2 * (IV_CAD$beta.exposure^2) * IV_CAD$eaf.exposure * (1 - IV_CAD$eaf.exposure)) /
  (2 * (IV_CAD$beta.exposure^2) * IV_CAD$eaf.exposure * (1 - IV_CAD$eaf.exposure) +
     2 * IV_CAD$samplesize.exposure * IV_CAD$eaf.exposure * (1 - IV_CAD$eaf.exposure) * (IV_CAD$se.exposure^2))

IV_CAD$F <- (IV_CAD$r2 * (IV_CAD$samplesize.exposure - 2))/(1 - IV_CAD$r2)
range(IV_CAD$F) #[1]  29.88125 443.09730
IV_CAD_meanF <- mean(IV_CAD$F) #[1] 60.60692

write.xlsx(IV_CAD, file = "Output/UVMR_IV/IV_CAD.xlsx")

out_CAD_ON_education <- extract_outcome_data(
  snps = IV_CAD$SNP,
  outcomes = 'ieu-a-1239')

UVdat_CAD_ON_education <- harmonise_data(
  exposure_dat =  IV_CAD, 
  outcome_dat = out_CAD_ON_education
)
str(UVdat_CAD_ON_education)
range(UVdat_CAD_ON_education$F) # [1]  29.88125 443.09730
tt <- mean(UVdat_CAD_ON_education$F)# [1] 60.60692
range(UVdat_CAD_ON_education$pval.outcome) # [1] 0.00471998 0.97300000

result_CAD_ON_education <- mr(UVdat_CAD_ON_education)

#散点图
scatter_plot_CAD_ON_education <- mr_scatter_plot(result_CAD_ON_education, UVdat_CAD_ON_education)
scatter_plot_CAD_ON_education[[1]]
ggsave(scatter_plot_CAD_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_CAD_ON_education.pdf", width=7, height=7)


pleiotropy_CAD_ON_education <- mr_pleiotropy_test(UVdat_CAD_ON_education)
str(pleiotropy_CAD_ON_education)
pleiotropy_CAD_ON_education <- pleiotropy_CAD_ON_education[,-c(1,2)]
pleiotropy_CAD_ON_education$outcome <- "educational attainment"


heterogeneity_CAD_ON_education <- mr_heterogeneity(UVdat_CAD_ON_education)
heterogeneity_CAD_ON_education <- heterogeneity_CAD_ON_education[,-c(1,2)]
heterogeneity_CAD_ON_education$outcome <- "educational attainment"


singleSNP_CAD_ON_education <- mr_singlesnp(UVdat_CAD_ON_education)
funnel_plot_CAD_ON_education <- mr_funnel_plot(singleSNP_CAD_ON_education)
funnel_plot_CAD_ON_education[[1]]
ggsave(funnel_plot_CAD_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_CAD_ON_education.pdf", width=7, height=7)


loo_CAD_ON_education <- mr_leaveoneout(UVdat_CAD_ON_education)
loo_plot_CAD_ON_education <- mr_leaveoneout_plot(loo_CAD_ON_education)
loo_plot_CAD_ON_education[[1]]
ggsave(loo_plot_CAD_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_CAD_ON_education.pdf", width=7, height=7)


singleSNP_plot_CAD_ON_education <- mr_forest_plot(singleSNP_CAD_ON_education)
singleSNP_plot_CAD_ON_education[[1]]
ggsave(singleSNP_plot_CAD_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_CAD_ON_education.pdf", width=7, height=7)

# MRPRESSO
PRESSO_CAD_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_CAD_ON_education, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_CAD_ON_education <- generate_odds_ratios(result_CAD_ON_education)
result_CAD_ON_education <- result_CAD_ON_education[,-c(1,2)]
result_CAD_ON_education$outcome <- "educational attainment"


###### 1.1.2 occupation & CAD ######
###### ______occupation on CAD ######
IV_occupation <- format_data(GWAS_occupation,
                         type = "exposure", 
                         snp_col = "SNP", 
                         beta_col = "beta",
                         se_col = "standard_error",
                         effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele",
                         eaf_col = "EAF",
                         pval_col = "p_value",
                         samplesize_col = "N")
str(IV_occupation)
IV_occupation$exposure <- "occupation"

IV_occupation <- IV_occupation[IV_occupation$pval.exposure < 5e-8,]
IV_occupation <- clump_data(IV_occupation)

IV_occupation$r2 <- (2 * (IV_occupation$beta.exposure^2) * IV_occupation$eaf.exposure * (1 - IV_occupation$eaf.exposure)) /
  (2 * (IV_occupation$beta.exposure^2) * IV_occupation$eaf.exposure * (1 - IV_occupation$eaf.exposure) +
     2 * IV_occupation$samplesize.exposure * IV_occupation$eaf.exposure * (1 - IV_occupation$eaf.exposure) * (IV_occupation$se.exposure^2))

IV_occupation$F <- (IV_occupation$r2 * (IV_occupation$samplesize.exposure - 2))/(1 - IV_occupation$r2)
range(IV_occupation$F) # [1] 29.80588 89.18191
IV_occupation_meanF <- mean(IV_occupation$F) # [1] 38.88091

write.xlsx(IV_occupation, file = "Output/UVMR_IV/IV_occupation.xlsx")

out_occupation_ON_CAD <- extract_outcome_data(
  snps = IV_occupation$SNP,
  outcomes = 'ieu-a-7')

UVdat_occupation_ON_CAD <- harmonise_data(
  exposure_dat =  IV_occupation, 
  outcome_dat = out_occupation_ON_CAD
)
str(UVdat_occupation_ON_CAD)
range(UVdat_occupation_ON_CAD$pval.outcome) 
range(UVdat_occupation_ON_CAD$F) 

write.xlsx(UVdat_occupation_ON_CAD, file = "Output/UVMR_IV/UVdat_occupation_ON_CAD.xlsx")

result_occupation_ON_CAD <- mr(UVdat_occupation_ON_CAD)

#散点图
scatter_plot_occupation_ON_CAD <- mr_scatter_plot(result_occupation_ON_CAD, UVdat_occupation_ON_CAD)
scatter_plot_occupation_ON_CAD[[1]]
ggsave(scatter_plot_occupation_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_occupation_ON_CAD.pdf", width=7, height=7)


pleiotropy_occupation_ON_CAD <- mr_pleiotropy_test(UVdat_occupation_ON_CAD)
str(pleiotropy_occupation_ON_CAD)
pleiotropy_occupation_ON_CAD <- pleiotropy_occupation_ON_CAD[,-c(1,2)]
pleiotropy_occupation_ON_CAD$outcome <- "coronary heart disease"


heterogeneity_occupation_ON_CAD <- mr_heterogeneity(UVdat_occupation_ON_CAD)
heterogeneity_occupation_ON_CAD <- heterogeneity_occupation_ON_CAD[,-c(1,2)]
heterogeneity_occupation_ON_CAD$outcome <- "coronary heart disease"


singleSNP_occupation_ON_CAD <- mr_singlesnp(UVdat_occupation_ON_CAD)
funnel_plot_occupation_ON_CAD <- mr_funnel_plot(singleSNP_occupation_ON_CAD)
funnel_plot_occupation_ON_CAD[[1]]
ggsave(funnel_plot_occupation_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_occupation_ON_CAD.pdf", width=7, height=7)


loo_occupation_ON_CAD <- mr_leaveoneout(UVdat_occupation_ON_CAD)
loo_plot_occupation_ON_CAD <- mr_leaveoneout_plot(loo_occupation_ON_CAD)
loo_plot_occupation_ON_CAD[[1]]
ggsave(loo_plot_occupation_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_occupation_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_occupation_ON_CAD <- mr_forest_plot(singleSNP_occupation_ON_CAD)
singleSNP_plot_occupation_ON_CAD[[1]]
ggsave(singleSNP_plot_occupation_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_occupation_ON_CAD.pdf", width=7, height=7)

# MRPRESSO
PRESSO_occupation_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_occupation_ON_CAD, NbDistribution = 10000,  
                                  SignifThreshold = 0.05)

result_occupation_ON_CAD <- generate_odds_ratios(result_occupation_ON_CAD)
result_occupation_ON_CAD <- result_occupation_ON_CAD[,-c(1,2)]
result_occupation_ON_CAD$outcome <- "coronary heart disease"


###### ______CAD on occupation ######
out_CAD_ON_occupation <- format_data(GWAS_occupation,
                                 type = "outcome", 
                                 snps = IV_CAD$SNP, 
                                 snp_col = "SNP",
                                 beta_col = "beta",
                                 se_col = "standard_error",
                                 effect_allele_col = "effect_allele",
                                 other_allele_col = "other_allele",
                                 eaf_col = "EAF",
                                 pval_col = "p_value")
out_CAD_ON_occupation$outcome <- "occupation"

UVdat_CAD_ON_occupation <- harmonise_data(
  exposure_dat =  IV_CAD, 
  outcome_dat = out_CAD_ON_occupation
)
str(UVdat_CAD_ON_occupation)
range(UVdat_CAD_ON_occupation$F)
tt <- mean(UVdat_CAD_ON_occupation$F)# [1] 62.19545
range(UVdat_CAD_ON_occupation$pval.outcome)

write.xlsx(UVdat_CAD_ON_occupation, file = "Output/UVMR_IV/UVdat_CAD_ON_occupation.xlsx")

result_CAD_ON_occupation <- mr(UVdat_CAD_ON_occupation)

#散点图
scatter_plot_CAD_ON_occupation <- mr_scatter_plot(result_CAD_ON_occupation, UVdat_CAD_ON_occupation)
scatter_plot_CAD_ON_occupation[[1]]
ggsave(scatter_plot_CAD_ON_occupation[[1]], file="Output/UVMR_Secondary Results/scatter_plot_CAD_ON_occupation.pdf", width=7, height=7)


pleiotropy_CAD_ON_occupation <- mr_pleiotropy_test(UVdat_CAD_ON_occupation)
str(pleiotropy_CAD_ON_occupation)
pleiotropy_CAD_ON_occupation <- pleiotropy_CAD_ON_occupation[,-c(1,2)]
pleiotropy_CAD_ON_occupation$outcome <- "occupation"


heterogeneity_CAD_ON_occupation <- mr_heterogeneity(UVdat_CAD_ON_occupation)
heterogeneity_CAD_ON_occupation <- heterogeneity_CAD_ON_occupation[,-c(1,2)]
heterogeneity_CAD_ON_occupation$outcome <- "occupation"


singleSNP_CAD_ON_occupation <- mr_singlesnp(UVdat_CAD_ON_occupation)
funnel_plot_CAD_ON_occupation <- mr_funnel_plot(singleSNP_CAD_ON_occupation)
funnel_plot_CAD_ON_occupation[[1]]
ggsave(funnel_plot_CAD_ON_occupation[[1]], file="Output/UVMR_Secondary Results/funnel_plot_CAD_ON_occupation.pdf", width=7, height=7)


loo_CAD_ON_occupation <- mr_leaveoneout(UVdat_CAD_ON_occupation)
loo_plot_CAD_ON_occupation <- mr_leaveoneout_plot(loo_CAD_ON_occupation)
loo_plot_CAD_ON_occupation[[1]]
ggsave(loo_plot_CAD_ON_occupation[[1]], file="Output/UVMR_Secondary Results/loo_plot_CAD_ON_occupation.pdf", width=7, height=7)


singleSNP_plot_CAD_ON_occupation <- mr_forest_plot(singleSNP_CAD_ON_occupation)
singleSNP_plot_CAD_ON_occupation[[1]]
ggsave(singleSNP_plot_CAD_ON_occupation[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_CAD_ON_occupation.pdf", width=7, height=7)

# MRPRESSO
PRESSO_CAD_ON_occupation <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_CAD_ON_occupation, NbDistribution = 10000,  
                                  SignifThreshold = 0.05)

result_CAD_ON_occupation <- generate_odds_ratios(result_CAD_ON_occupation)
result_CAD_ON_occupation <- result_CAD_ON_occupation[,-c(1,2)]
result_CAD_ON_occupation$outcome <- "occupation"


###### 1.1.3 income & CAD ######
###### ______income on CAD ######
IV_income <- format_data(GWAS_income,
                      type = "exposure", 
                      snp_col = "SNP", 
                      beta_col = "ES",
                      se_col = "SE",
                      effect_allele_col = "ALT",
                      other_allele_col = "REF",
                      eaf_col = "AF",
                      pval_col = "P",
                      samplesize_col = "SS")
str(IV_income)
IV_income$exposure <- "household income"

IV_income <- IV_income[IV_income$pval.exposure < 5e-8,]
IV_income <- clump_data(IV_income)

IV_income$r2 <- (2 * (IV_income$beta.exposure^2) * IV_income$eaf.exposure * (1 - IV_income$eaf.exposure)) /
  (2 * (IV_income$beta.exposure^2) * IV_income$eaf.exposure * (1 - IV_income$eaf.exposure) +
     2 * IV_income$samplesize.exposure * IV_income$eaf.exposure * (1 - IV_income$eaf.exposure) * (IV_income$se.exposure^2))

IV_income$F <- (IV_income$r2 * (IV_income$samplesize.exposure - 2))/(1 - IV_income$r2)
range(IV_income$F) # [1]  29.98273 101.38726
IV_income_meanF <- mean(IV_income$F) # [1] 40.98244

write.xlsx(IV_income, file = "Output/UVMR_IV/IV_income.xlsx")

out_income_ON_CAD <- extract_outcome_data(
  snps = IV_income$SNP,
  outcomes = 'ieu-a-7')

UVdat_income_ON_CAD <- harmonise_data(
  exposure_dat =  IV_income, 
  outcome_dat = out_income_ON_CAD
)
str(UVdat_income_ON_CAD)
range(UVdat_income_ON_CAD$pval.outcome) 
range(UVdat_income_ON_CAD$F) 

write.xlsx(UVdat_income_ON_CAD, file = "Output/UVMR_IV/UVdat_income_ON_CAD.xlsx")

result_income_ON_CAD <- mr(UVdat_income_ON_CAD)

#散点图
scatter_plot_income_ON_CAD <- mr_scatter_plot(result_income_ON_CAD, UVdat_income_ON_CAD)
scatter_plot_income_ON_CAD[[1]]
ggsave(scatter_plot_income_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_income_ON_CAD.pdf", width=7, height=7)


pleiotropy_income_ON_CAD <- mr_pleiotropy_test(UVdat_income_ON_CAD)
str(pleiotropy_income_ON_CAD)
pleiotropy_income_ON_CAD <- pleiotropy_income_ON_CAD[,-c(1,2)]
pleiotropy_income_ON_CAD$outcome <- "coronary heart disease"


heterogeneity_income_ON_CAD <- mr_heterogeneity(UVdat_income_ON_CAD)
heterogeneity_income_ON_CAD <- heterogeneity_income_ON_CAD[,-c(1,2)]
heterogeneity_income_ON_CAD$outcome <- "coronary heart disease"


singleSNP_income_ON_CAD <- mr_singlesnp(UVdat_income_ON_CAD)
funnel_plot_income_ON_CAD <- mr_funnel_plot(singleSNP_income_ON_CAD)
funnel_plot_income_ON_CAD[[1]]
ggsave(funnel_plot_income_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_income_ON_CAD.pdf", width=7, height=7)


loo_income_ON_CAD <- mr_leaveoneout(UVdat_income_ON_CAD)
loo_plot_income_ON_CAD <- mr_leaveoneout_plot(loo_income_ON_CAD)
loo_plot_income_ON_CAD[[1]]
ggsave(loo_plot_income_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_income_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_income_ON_CAD <- mr_forest_plot(singleSNP_income_ON_CAD)
singleSNP_plot_income_ON_CAD[[1]]
ggsave(singleSNP_plot_income_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_income_ON_CAD.pdf", width=7, height=7)

# MRPRESSO
PRESSO_income_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_income_ON_CAD, NbDistribution = 10000,  
                                  SignifThreshold = 0.05)

result_income_ON_CAD <- generate_odds_ratios(result_income_ON_CAD)
result_income_ON_CAD <- result_income_ON_CAD[,-c(1,2)]
result_income_ON_CAD$outcome <- "coronary heart disease"

###### ______CAD on income ######
out_CAD_ON_income <- format_data(GWAS_income,
                      type = "outcome", 
                      snps = IV_CAD$SNP, 
                      snp_col = "SNP",
                      beta_col = "ES",
                      se_col = "SE",
                      effect_allele_col = "ALT",
                      other_allele_col = "REF",
                      eaf_col = "AF",
                      pval_col = "P")
out_CAD_ON_income$outcome <- "household income"

UVdat_CAD_ON_income <- harmonise_data(
  exposure_dat =  IV_CAD, 
  outcome_dat = out_CAD_ON_income
)
range(UVdat_CAD_ON_income$F)
tt <- mean(UVdat_CAD_ON_income$F)# [1] 61.8569
range(UVdat_CAD_ON_income$pval.outcome)

write.xlsx(UVdat_CAD_ON_income, file = "Output/UVMR_IV/UVdat_CAD_ON_income.xlsx")

result_CAD_ON_income <- mr(UVdat_CAD_ON_income)

#散点图
scatter_plot_CAD_ON_income <- mr_scatter_plot(result_CAD_ON_income, UVdat_CAD_ON_income)
scatter_plot_CAD_ON_income[[1]]
ggsave(scatter_plot_CAD_ON_income[[1]], file="Output/UVMR_Secondary Results/scatter_plot_CAD_ON_income.pdf", width=7, height=7)


pleiotropy_CAD_ON_income <- mr_pleiotropy_test(UVdat_CAD_ON_income)
str(pleiotropy_CAD_ON_income)
pleiotropy_CAD_ON_income <- pleiotropy_CAD_ON_income[,-c(1,2)]
pleiotropy_CAD_ON_income$outcome <- "household income"


heterogeneity_CAD_ON_income <- mr_heterogeneity(UVdat_CAD_ON_income)
heterogeneity_CAD_ON_income <- heterogeneity_CAD_ON_income[,-c(1,2)]
heterogeneity_CAD_ON_income$outcome <- "household income"


singleSNP_CAD_ON_income <- mr_singlesnp(UVdat_CAD_ON_income)
funnel_plot_CAD_ON_income <- mr_funnel_plot(singleSNP_CAD_ON_income)
funnel_plot_CAD_ON_income[[1]]
ggsave(funnel_plot_CAD_ON_income[[1]], file="Output/UVMR_Secondary Results/funnel_plot_CAD_ON_income.pdf", width=7, height=7)


loo_CAD_ON_income <- mr_leaveoneout(UVdat_CAD_ON_income)
loo_plot_CAD_ON_income <- mr_leaveoneout_plot(loo_CAD_ON_income)
loo_plot_CAD_ON_income[[1]]
ggsave(loo_plot_CAD_ON_income[[1]], file="Output/UVMR_Secondary Results/loo_plot_CAD_ON_income.pdf", width=7, height=7)


singleSNP_plot_CAD_ON_income <- mr_forest_plot(singleSNP_CAD_ON_income)
singleSNP_plot_CAD_ON_income[[1]]
ggsave(singleSNP_plot_CAD_ON_income[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_CAD_ON_income.pdf", width=7, height=7)

# MRPRESSO
PRESSO_CAD_ON_income <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_CAD_ON_income, NbDistribution = 10000,  
                                  SignifThreshold = 0.05)

result_CAD_ON_income <- generate_odds_ratios(result_CAD_ON_income)
result_CAD_ON_income <- result_CAD_ON_income[,-c(1,2)]
result_CAD_ON_income$outcome <- "household income"


###### 1.1.4 TDI & CAD ######
###### ______TDI on CAD ######
IV_TDI <- fread("Input/IV_TDI.csv")
IV_TDI <- as.data.frame(IV_TDI)
str(IV_TDI)
IV_TDI <- IV_TDI[,-1]
IV_TDI$samplesize.exposure <- as.numeric(IV_TDI$samplesize.exposure)
IV_TDI<-format_data(IV_TDI,
                    type = "exposure", 
                    snp_col = "SNP", 
                    beta_col = "beta.exposure",
                    se_col = "se.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure",
                    eaf_col = "eaf.exposure",
                    pval_col = "pval.exposure",
                    samplesize_col = "samplesize.exposure")
IV_TDI$exposure <- "regional deprivation"

IV_TDI$r2 <- (2 * (IV_TDI$beta.exposure^2) * IV_TDI$eaf.exposure * (1 - IV_TDI$eaf.exposure)) /
  (2 * (IV_TDI$beta.exposure^2) * IV_TDI$eaf.exposure * (1 - IV_TDI$eaf.exposure) +
     2 * IV_TDI$samplesize.exposure * IV_TDI$eaf.exposure * (1 - IV_TDI$eaf.exposure) * (IV_TDI$se.exposure^2))

IV_TDI$F <- (IV_TDI$r2 * (IV_TDI$samplesize.exposure - 2))/(1 - IV_TDI$r2)
range(IV_TDI$F) #[1] 30.19642 49.76711
IV_TDI_meanF <- mean(IV_TDI$F) #[1] 35.43213

write.xlsx(IV_TDI, file = "Output/UVMR_IV/IV_TDI.xlsx")

out_TDI_ON_CAD <- extract_outcome_data(
  snps = IV_TDI$SNP,
  outcomes = 'ieu-a-7')
write.xlsx(out_TDI_ON_CAD, file = "Output/UVMR_IV/out_TDI_ON_CAD.xlsx")

UVdat_TDI_ON_CAD <- harmonise_data(
  exposure_dat =  IV_TDI, 
  outcome_dat = out_TDI_ON_CAD
)
str(UVdat_TDI_ON_CAD)
range(UVdat_TDI_ON_CAD$pval.outcome)
range(UVdat_TDI_ON_CAD$F)
write.xlsx(UVdat_TDI_ON_CAD, file = "Output/UVMR_IV/UVdat_TDI_ON_CAD.xlsx")

result_TDI_ON_CAD <- mr(UVdat_TDI_ON_CAD)

#散点图
scatter_plot_TDI_ON_CAD <- mr_scatter_plot(result_TDI_ON_CAD, UVdat_TDI_ON_CAD)
scatter_plot_TDI_ON_CAD[[1]]
ggsave(scatter_plot_TDI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_TDI_ON_CAD.pdf", width=7, height=7)


pleiotropy_TDI_ON_CAD <- mr_pleiotropy_test(UVdat_TDI_ON_CAD)
str(pleiotropy_TDI_ON_CAD)
pleiotropy_TDI_ON_CAD <- pleiotropy_TDI_ON_CAD[,-c(1,2)]
pleiotropy_TDI_ON_CAD$outcome <- "coronary heart disease"


heterogeneity_TDI_ON_CAD <- mr_heterogeneity(UVdat_TDI_ON_CAD)
heterogeneity_TDI_ON_CAD <- heterogeneity_TDI_ON_CAD[,-c(1,2)]
heterogeneity_TDI_ON_CAD$outcome <- "coronary heart disease"


singleSNP_TDI_ON_CAD <- mr_singlesnp(UVdat_TDI_ON_CAD)
funnel_plot_TDI_ON_CAD <- mr_funnel_plot(singleSNP_TDI_ON_CAD)
funnel_plot_TDI_ON_CAD[[1]]
ggsave(funnel_plot_TDI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_TDI_ON_CAD.pdf", width=7, height=7)


loo_TDI_ON_CAD <- mr_leaveoneout(UVdat_TDI_ON_CAD)
loo_plot_TDI_ON_CAD <- mr_leaveoneout_plot(loo_TDI_ON_CAD)
loo_plot_TDI_ON_CAD[[1]]
ggsave(loo_plot_TDI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_TDI_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_TDI_ON_CAD <- mr_forest_plot(singleSNP_TDI_ON_CAD)
singleSNP_plot_TDI_ON_CAD[[1]]
ggsave(singleSNP_plot_TDI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_TDI_ON_CAD.pdf", width=7, height=7)

# MRPRESSO
PRESSO_TDI_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_TDI_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_TDI_ON_CAD <- generate_odds_ratios(result_TDI_ON_CAD)
result_TDI_ON_CAD <- result_TDI_ON_CAD[,-c(1,2)]
result_TDI_ON_CAD$outcome <- "coronary heart disease"


###### ______CAD on TDI ######
out_CAD_ON_TDI <- extract_outcome_data(
  snps = IV_CAD$SNP,
  outcomes = 'ukb-b-10011')

UVdat_CAD_ON_TDI <- harmonise_data(
  exposure_dat =  IV_CAD, 
  outcome_dat = out_CAD_ON_TDI
)
range(UVdat_CAD_ON_TDI$F)
tt <- mean(UVdat_CAD_ON_TDI$F)#[1] 60.60692
range(UVdat_CAD_ON_TDI$pval.outcome)

write.xlsx(out_CAD_ON_TDI,file="Output/UVMR_IV/out_CAD_ON_TDI.xlsx")
write.xlsx(UVdat_CAD_ON_TDI,file="Output/UVMR_IV/UVdat_CAD_ON_TDI.xlsx")

result_CAD_ON_TDI <- mr(UVdat_CAD_ON_TDI)

#散点图
scatter_plot_CAD_ON_TDI <- mr_scatter_plot(result_CAD_ON_TDI, UVdat_CAD_ON_TDI)
scatter_plot_CAD_ON_TDI[[1]]
ggsave(scatter_plot_CAD_ON_TDI[[1]], file="Output/UVMR_Secondary Results/scatter_plot_CAD_ON_TDI.pdf", width=7, height=7)


pleiotropy_CAD_ON_TDI <- mr_pleiotropy_test(UVdat_CAD_ON_TDI)
str(pleiotropy_CAD_ON_TDI)
pleiotropy_CAD_ON_TDI <- pleiotropy_CAD_ON_TDI[,-c(1,2)]
pleiotropy_CAD_ON_TDI$outcome <- "regional deprivation"


heterogeneity_CAD_ON_TDI <- mr_heterogeneity(UVdat_CAD_ON_TDI)
str(heterogeneity_CAD_ON_TDI)
heterogeneity_CAD_ON_TDI <- heterogeneity_CAD_ON_TDI[,-c(1,2)]
heterogeneity_CAD_ON_TDI$outcome <- "regional deprivation"


singleSNP_CAD_ON_TDI <- mr_singlesnp(UVdat_CAD_ON_TDI)
funnel_plot_CAD_ON_TDI <- mr_funnel_plot(singleSNP_CAD_ON_TDI)
funnel_plot_CAD_ON_TDI[[1]]
ggsave(funnel_plot_CAD_ON_TDI[[1]], file="Output/UVMR_Secondary Results/funnel_plot_CAD_ON_TDI.pdf", width=7, height=7)


loo_CAD_ON_TDI <- mr_leaveoneout(UVdat_CAD_ON_TDI)
loo_plot_CAD_ON_TDI <- mr_leaveoneout_plot(loo_CAD_ON_TDI)
loo_plot_CAD_ON_TDI[[1]]
ggsave(loo_plot_CAD_ON_TDI[[1]], file="Output/UVMR_Secondary Results/loo_plot_CAD_ON_TDI.pdf", width=7, height=7)


singleSNP_plot_CAD_ON_TDI <- mr_forest_plot(singleSNP_CAD_ON_TDI)
singleSNP_plot_CAD_ON_TDI[[1]]
ggsave(singleSNP_plot_CAD_ON_TDI[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_CAD_ON_TDI.pdf", width=7, height=7)

# MRPRESSO
PRESSO_CAD_ON_TDI <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_CAD_ON_TDI, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_CAD_ON_TDI <- generate_odds_ratios(result_CAD_ON_TDI)
result_CAD_ON_TDI <- result_CAD_ON_TDI[,-c(1,2)]
result_CAD_ON_TDI$outcome <- "regional deprivation"


###### Main Results
SummaryResultsUVMR_SES_ON_CAD <- rbind(result_education_ON_CAD,
                                       result_occupation_ON_CAD,
                                       result_income_ON_CAD,
                                       result_TDI_ON_CAD)
write.xlsx(SummaryResultsUVMR_SES_ON_CAD, file = "Output/UVMR_Main Results/SummaryResultsUVMR_SES_ON_CAD.xlsx")

SummaryResultsUVMR_CAD_ON_SES <- rbind(result_CAD_ON_education,
                                       result_CAD_ON_occupation,
                                       result_CAD_ON_income,
                                       result_CAD_ON_TDI)
write.xlsx(SummaryResultsUVMR_CAD_ON_SES, file = "Output/UVMR_Main Results/SummaryResultsUVMR_CAD_ON_SES.xlsx")


###### Secondary Results

Summary_UVPleiotropy_SES_ON_CAD <- rbind(pleiotropy_education_ON_CAD,
                                         pleiotropy_occupation_ON_CAD,
                                         pleiotropy_income_ON_CAD,
                                         pleiotropy_TDI_ON_CAD)
write.xlsx(Summary_UVPleiotropy_SES_ON_CAD,"Output/UVMR_Secondary Results/Summary_UVPleiotropy_SES_ON_CAD.xlsx")

Summary_UVPleiotropy_CAD_ON_SES <- rbind(pleiotropy_CAD_ON_education,
                                         pleiotropy_CAD_ON_occupation,
                                         pleiotropy_CAD_ON_income,
                                         pleiotropy_CAD_ON_TDI)
write.xlsx(Summary_UVPleiotropy_CAD_ON_SES,"Output/UVMR_Secondary Results/Summary_UVPleiotropy_CAD_ON_SES.xlsx")


Summary_UVHeterogeneity_SES_ON_CAD <- rbind(heterogeneity_education_ON_CAD,
                                            heterogeneity_occupation_ON_CAD,
                                            heterogeneity_income_ON_CAD,
                                            heterogeneity_TDI_ON_CAD)
write.xlsx(Summary_UVHeterogeneity_SES_ON_CAD,"Output/UVMR_Secondary Results/Summary_UVHeterogeneity_SES_ON_CAD.xlsx")

Summary_UVHeterogeneity_CAD_ON_SES <- rbind(heterogeneity_CAD_ON_education,
                                            heterogeneity_CAD_ON_occupation,
                                            heterogeneity_CAD_ON_income,
                                            heterogeneity_CAD_ON_TDI)
write.xlsx(Summary_UVHeterogeneity_CAD_ON_SES,"Output/UVMR_Secondary Results/Summary_UVHeterogeneity_CAD_ON_SES.xlsx")


# the direct effect of education, occupation, income, TDI on CAD
IV_educationStringent <- subset(GWAS_education,P<5e-8)
IV_educationStringent <- clump_data(IV_educationStringent)

IV_occupation <- subset(GWAS_occupation,p_value<5e-8)
IV_occupation <- clump_data(IV_occupation)

IV_income <- subset(GWAS_income,P<5e-8)
IV_income <- clump_data(IV_income)

IV_TDI <- subset(GWAS_TDI,P<5e-8)
IV_TDI <- clump_data(IV_TDI)

###### _____CAD ~ education + occupation ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_occupation <- IV_occupation["SNP"]

SNP_IV_education.occupation <- rbind(SNP_IV_educationStringent,SNP_IV_occupation)# 377SNP
SNP_IV_education.occupation <- unique(SNP_IV_education.occupation) # 1SNP
SNP_IV_education.occupation <- clump_data(SNP_IV_education.occupation) # 352SNP
write.xlsx(SNP_IV_education.occupation, file = "Output/MVMR_IV/SNP_IV_education.occupation.xlsx")

exp_education_MVwith.occupation <- format_data(
  GWAS_education,
  snps = SNP_IV_education.occupation$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.occupation$exposure <- "education"

exp_occupation_MVwith.education <- format_data(
  GWAS_occupation,
  snps = SNP_IV_education.occupation$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value") 
exp_occupation_MVwith.education$outcome <- "occupation"

education.occupation <- harmonise_data(exposure_dat = exp_education_MVwith.occupation, outcome_dat = exp_occupation_MVwith.education, action = 1)

out_CAD_MV.education.occupation <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.occupation$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.occupation$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.occupation, outcome_dat = out_CAD_MV.education.occupation, action = 1)

education.occupation.dat <- education.occupation[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.occupation.dat) <- c("SNP", "beta.education", "beta.occupation", "se.education", "se.occupation")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.occupation_ON_CAD <- merge(education.occupation.dat, education.CAD.dat, by = "SNP") # 334SNP
write.xlsx(MVdat_education.occupation_ON_CAD, file = "Output/MVMR_IV/MVdat_education.occupation_ON_CAD.xlsx")
MVdat_education.occupation_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_education.occupation_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.occupation_ON_CAD[,c("beta.education","beta.occupation")])
bxse = as.matrix(MVdat_education.occupation_ON_CAD[,c("se.education","se.occupation")])
by = as.vector(MVdat_education.occupation_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.occupation_ON_CAD$se.CAD)

MVdatForm_education.occupation_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","occupation"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.occupation_ON_CAD <- mr_mvivw(MVdatForm_education.occupation_ON_CAD)
tt <- mvmr(MVdatForm_education.occupation_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.occupation_ON_CAD <- mr_mvegger(MVdatForm_education.occupation_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.occupation_ON_CAD <- mr_mvlasso(MVdatForm_education.occupation_ON_CAD)

# MV Median
result_MV.Median_education.occupation_ON_CAD <- mr_mvmedian(MVdatForm_education.occupation_ON_CAD)



F_education.occupation <- strength_mvmr(r_input = MVdatForm_education.occupation_ON_CAD, gencov = 0)


mv_hete_education.occupation <- pleiotropy_mvmr(r_input = MVdatForm_education.occupation_ON_CAD, gencov = 0)


PRESSO_education.occupation_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                         BetaExposure = c("beta.education", "beta.occupation"), 
                                         SdOutcome = "se.CAD", 
                                         SdExposure = c("se.education", "se.occupation"),
                                         OUTLIERtest = TRUE, 
                                         DISTORTIONtest = TRUE, 
                                         data = MVdat_education.occupation_ON_CAD,
                                         NbDistribution = 1000, 
                                         SignifThreshold = 0.05)


###### _____CAD ~ education + income ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_income <- IV_income["SNP"]

SNP_IV_education.income <- rbind(SNP_IV_educationStringent,SNP_IV_income)# 401SNP
SNP_IV_education.income <- unique(SNP_IV_education.income) # 5SNP
SNP_IV_education.income <- clump_data(SNP_IV_education.income) # 351SNP
write.xlsx(SNP_IV_education.income, file = "Output/MVMR_IV/SNP_IV_education.income.xlsx")

exp_education_MVwith.income <- format_data(
  GWAS_education,
  snps = SNP_IV_education.income$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.income$exposure <- "education"

exp_income_MVwith.education <- format_data(
  GWAS_income,
  snps = SNP_IV_education.income$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_income_MVwith.education$outcome <- "income"

education.income <- harmonise_data(exposure_dat = exp_education_MVwith.income, outcome_dat = exp_income_MVwith.education, action = 1)

out_CAD_MV.education.income <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.income$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.income$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.income, outcome_dat = out_CAD_MV.education.income, action = 1)

education.income.dat <- education.income[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.income.dat) <- c("SNP", "beta.education", "beta.income", "se.education", "se.income")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.income_ON_CAD <- merge(education.income.dat, education.CAD.dat, by = "SNP") # 350SNP
write.xlsx(MVdat_education.income_ON_CAD, file = "Output/MVMR_IV/MVdat_education.income_ON_CAD.xlsx")
MVdat_education.income_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_education.income_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.income_ON_CAD[,c("beta.education","beta.income")])
bxse = as.matrix(MVdat_education.income_ON_CAD[,c("se.education","se.income")])
by = as.vector(MVdat_education.income_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.income_ON_CAD$se.CAD)

MVdatForm_education.income_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","income"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.income_ON_CAD <- mr_mvivw(MVdatForm_education.income_ON_CAD)
tt <- mvmr(MVdatForm_education.income_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.income_ON_CAD <- mr_mvegger(MVdatForm_education.income_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.income_ON_CAD <- mr_mvlasso(MVdatForm_education.income_ON_CAD)

# MV Median
result_MV.Median_education.income_ON_CAD <- mr_mvmedian(MVdatForm_education.income_ON_CAD)



F_education.income <- strength_mvmr(r_input = MVdatForm_education.income_ON_CAD, gencov = 0)


mv_hete_education.income <- pleiotropy_mvmr(r_input = MVdatForm_education.income_ON_CAD, gencov = 0)


PRESSO_education.income_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.education", "beta.income"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.education", "se.income"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_education.income_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)


###### _____CAD ~ education + TDI ######
str(IV_educationStringent)
str(IV_TDI)

SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_TDI <- IV_TDI["SNP"]

SNP_IV_education.TDI <- rbind(SNP_IV_educationStringent,SNP_IV_TDI)# 366SNP
SNP_IV_education.TDI <- unique(SNP_IV_education.TDI) # 0SNP
SNP_IV_education.TDI <- clump_data(SNP_IV_education.TDI) # 348SNP
write.xlsx(SNP_IV_education.TDI, file = "Output/MVMR_IV/SNP_IV_education.TDI.xlsx")

exp_education_MVwith.TDI <- format_data(
  GWAS_education,
  snps = SNP_IV_education.TDI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.TDI$exposure <- "education"

exp_TDI_MVwith.education <- format_data(
  GWAS_TDI,
  snps = SNP_IV_education.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_TDI_MVwith.education$outcome <- "TDI"

education.TDI <- harmonise_data(exposure_dat = exp_education_MVwith.TDI, outcome_dat = exp_TDI_MVwith.education, action = 1)

out_CAD_MV.education.TDI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.TDI$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.TDI, outcome_dat = out_CAD_MV.education.TDI, action = 1)

education.TDI.dat <- education.TDI[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.TDI.dat) <- c("SNP", "beta.education", "beta.TDI", "se.education", "se.TDI")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.TDI_ON_CAD <- merge(education.TDI.dat, education.CAD.dat, by = "SNP") # 347SNP
write.xlsx(MVdat_education.TDI_ON_CAD, file = "Output/MVMR_IV/MVdat_education.TDI_ON_CAD.xlsx")
MVdat_education.TDI_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_education.TDI_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.TDI_ON_CAD[,c("beta.education","beta.TDI")])
bxse = as.matrix(MVdat_education.TDI_ON_CAD[,c("se.education","se.TDI")])
by = as.vector(MVdat_education.TDI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.TDI_ON_CAD$se.CAD)

MVdatForm_education.TDI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","TDI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.TDI_ON_CAD <- mr_mvivw(MVdatForm_education.TDI_ON_CAD)
tt <- mvmr(MVdatForm_education.TDI_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.TDI_ON_CAD <- mr_mvegger(MVdatForm_education.TDI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.TDI_ON_CAD <- mr_mvlasso(MVdatForm_education.TDI_ON_CAD)

# MV Median
result_MV.Median_education.TDI_ON_CAD <- mr_mvmedian(MVdatForm_education.TDI_ON_CAD)



F_education.TDI <- strength_mvmr(r_input = MVdatForm_education.TDI_ON_CAD, gencov = 0)


mv_hete_education.TDI <- pleiotropy_mvmr(r_input = MVdatForm_education.TDI_ON_CAD, gencov = 0)


PRESSO_education.TDI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                            BetaExposure = c("beta.education", "beta.TDI"), 
                                            SdOutcome = "se.CAD", 
                                            SdExposure = c("se.education", "se.TDI"),
                                            OUTLIERtest = TRUE, 
                                            DISTORTIONtest = TRUE, 
                                            data = MVdat_education.TDI_ON_CAD,
                                            NbDistribution = 1000, 
                                            SignifThreshold = 0.05)


###### _____CAD ~ occupation + income ######
str(IV_occupation)
str(IV_income)

SNP_IV_occupation <- IV_occupation["SNP"]
SNP_IV_income <- IV_income["SNP"]

SNP_IV_occupation.income <- rbind(SNP_IV_occupation,SNP_IV_income)# 82SNP
SNP_IV_occupation.income <- unique(SNP_IV_occupation.income) # 1SNP
SNP_IV_occupation.income <- clump_data(SNP_IV_occupation.income) # 66SNP
write.xlsx(SNP_IV_occupation.income, file = "Output/MVMR_IV/SNP_IV_occupation.income.xlsx")

exp_occupation_MVwith.income <- format_data(
  GWAS_occupation,
  snps = SNP_IV_occupation.income$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value") 
exp_occupation_MVwith.income$exposure <- "occupation"

exp_income_MVwith.occupation <- format_data(
  GWAS_income,
  snps = SNP_IV_occupation.income$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_income_MVwith.occupation$outcome <- "income"

occupation.income <- harmonise_data(exposure_dat = exp_occupation_MVwith.income, outcome_dat = exp_income_MVwith.occupation, action = 1)

out_CAD_MV.occupation.income <- format_data(
  GWAS_CAD,
  snps = SNP_IV_occupation.income$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.occupation.income$outcome <- "CAD"

occupation.CAD <- harmonise_data(exposure_dat = exp_occupation_MVwith.income, outcome_dat = out_CAD_MV.occupation.income, action = 1)

occupation.income.dat <- occupation.income[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(occupation.income.dat) <- c("SNP", "beta.occupation", "beta.income", "se.occupation", "se.income")

occupation.CAD.dat <- occupation.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(occupation.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_occupation.income_ON_CAD <- merge(occupation.income.dat, occupation.CAD.dat, by = "SNP") # 59SNP
write.xlsx(MVdat_occupation.income_ON_CAD, file = "Output/MVMR_IV/MVdat_occupation.income_ON_CAD.xlsx")
MVdat_occupation.income_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_occupation.income_ON_CAD.xlsx")

bx = as.matrix(MVdat_occupation.income_ON_CAD[,c("beta.occupation","beta.income")])
bxse = as.matrix(MVdat_occupation.income_ON_CAD[,c("se.occupation","se.income")])
by = as.vector(MVdat_occupation.income_ON_CAD$beta.CAD)
byse = as.vector(MVdat_occupation.income_ON_CAD$se.CAD)

MVdatForm_occupation.income_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("occupation","income"),outcome = "CAD")

# MV-IVW
result_MV.IVW_occupation.income_ON_CAD <- mr_mvivw(MVdatForm_occupation.income_ON_CAD)
tt <- mvmr(MVdatForm_occupation.income_ON_CAD)

# MV MR-Egger
result_MV.Egger_occupation.income_ON_CAD <- mr_mvegger(MVdatForm_occupation.income_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_occupation.income_ON_CAD <- mr_mvlasso(MVdatForm_occupation.income_ON_CAD)

# MV Median
result_MV.Median_occupation.income_ON_CAD <- mr_mvmedian(MVdatForm_occupation.income_ON_CAD)



F_occupation.income <- strength_mvmr(r_input = MVdatForm_occupation.income_ON_CAD, gencov = 0)


mv_hete_occupation.income <- pleiotropy_mvmr(r_input = MVdatForm_occupation.income_ON_CAD, gencov = 0)


PRESSO_occupation.income_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                            BetaExposure = c("beta.occupation", "beta.income"), 
                                            SdOutcome = "se.CAD", 
                                            SdExposure = c("se.occupation", "se.income"),
                                            OUTLIERtest = TRUE, 
                                            DISTORTIONtest = TRUE, 
                                            data = MVdat_occupation.income_ON_CAD,
                                            NbDistribution = 1000, 
                                            SignifThreshold = 0.05)


###### _____CAD ~ occupation + TDI ######
str(IV_occupation)
str(IV_TDI)

SNP_IV_occupation <- IV_occupation["SNP"]
SNP_IV_TDI <- IV_TDI["SNP"]

SNP_IV_occupation.TDI <- rbind(SNP_IV_occupation,SNP_IV_TDI)# 47SNP
SNP_IV_occupation.TDI <- unique(SNP_IV_occupation.TDI) # 0SNP
SNP_IV_occupation.TDI <- clump_data(SNP_IV_occupation.TDI) # 46SNP
write.xlsx(SNP_IV_occupation.TDI, file = "Output/MVMR_IV/SNP_IV_occupation.TDI.xlsx")

exp_occupation_MVwith.TDI <- format_data(
  GWAS_occupation,
  snps = SNP_IV_occupation.TDI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value") 
exp_occupation_MVwith.TDI$exposure <- "occupation"

exp_TDI_MVwith.occupation <- format_data(
  GWAS_TDI,
  snps = SNP_IV_occupation.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_TDI_MVwith.occupation$outcome <- "TDI"

occupation.TDI <- harmonise_data(exposure_dat = exp_occupation_MVwith.TDI, outcome_dat = exp_TDI_MVwith.occupation, action = 1)

out_CAD_MV.occupation.TDI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_occupation.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.occupation.TDI$outcome <- "CAD"

occupation.CAD <- harmonise_data(exposure_dat = exp_occupation_MVwith.TDI, outcome_dat = out_CAD_MV.occupation.TDI, action = 1)

occupation.TDI.dat <- occupation.TDI[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(occupation.TDI.dat) <- c("SNP", "beta.occupation", "beta.TDI", "se.occupation", "se.TDI")

occupation.CAD.dat <- occupation.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(occupation.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_occupation.TDI_ON_CAD <- merge(occupation.TDI.dat, occupation.CAD.dat, by = "SNP") # 40SNP
write.xlsx(MVdat_occupation.TDI_ON_CAD, file = "Output/MVMR_IV/MVdat_occupation.TDI_ON_CAD.xlsx")
MVdat_occupation.TDI_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_occupation.TDI_ON_CAD.xlsx")

bx = as.matrix(MVdat_occupation.TDI_ON_CAD[,c("beta.occupation","beta.TDI")])
bxse = as.matrix(MVdat_occupation.TDI_ON_CAD[,c("se.occupation","se.TDI")])
by = as.vector(MVdat_occupation.TDI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_occupation.TDI_ON_CAD$se.CAD)

MVdatForm_occupation.TDI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("occupation","TDI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_occupation.TDI_ON_CAD <- mr_mvivw(MVdatForm_occupation.TDI_ON_CAD)
tt <- mvmr(MVdatForm_occupation.TDI_ON_CAD)

# MV MR-Egger
result_MV.Egger_occupation.TDI_ON_CAD <- mr_mvegger(MVdatForm_occupation.TDI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_occupation.TDI_ON_CAD <- mr_mvlasso(MVdatForm_occupation.TDI_ON_CAD)

# MV Median
result_MV.Median_occupation.TDI_ON_CAD <- mr_mvmedian(MVdatForm_occupation.TDI_ON_CAD)



F_occupation.TDI <- strength_mvmr(r_input = MVdatForm_occupation.TDI_ON_CAD, gencov = 0)


mv_hete_occupation.TDI <- pleiotropy_mvmr(r_input = MVdatForm_occupation.TDI_ON_CAD, gencov = 0)


PRESSO_occupation.TDI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                             BetaExposure = c("beta.occupation", "beta.TDI"), 
                                             SdOutcome = "se.CAD", 
                                             SdExposure = c("se.occupation", "se.TDI"),
                                             OUTLIERtest = TRUE, 
                                             DISTORTIONtest = TRUE, 
                                             data = MVdat_occupation.TDI_ON_CAD,
                                             NbDistribution = 1000, 
                                             SignifThreshold = 0.05)


###### _____CAD ~ income + TDI ######
str(IV_income)
str(IV_TDI)

SNP_IV_income <- IV_income["SNP"]
SNP_IV_TDI <- IV_TDI["SNP"]

SNP_IV_income.TDI <- rbind(SNP_IV_income,SNP_IV_TDI)# 71SNP
SNP_IV_income.TDI <- unique(SNP_IV_income.TDI) # 1SNP
SNP_IV_income.TDI <- clump_data(SNP_IV_income.TDI) # 67SNP
write.xlsx(SNP_IV_income.TDI, file = "Output/MVMR_IV/SNP_IV_income.TDI.xlsx")

exp_income_MVwith.TDI <- format_data(
  GWAS_income,
  snps = SNP_IV_income.TDI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_income_MVwith.TDI$exposure <- "income"

exp_TDI_MVwith.income <- format_data(
  GWAS_TDI,
  snps = SNP_IV_income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_TDI_MVwith.income$outcome <- "TDI"

income.TDI <- harmonise_data(exposure_dat = exp_income_MVwith.TDI, outcome_dat = exp_TDI_MVwith.income, action = 1)

out_CAD_MV.income.TDI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.income.TDI$outcome <- "CAD"

income.CAD <- harmonise_data(exposure_dat = exp_income_MVwith.TDI, outcome_dat = out_CAD_MV.income.TDI, action = 1)

income.TDI.dat <- income.TDI[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(income.TDI.dat) <- c("SNP", "beta.income", "beta.TDI", "se.income", "se.TDI")

income.CAD.dat <- income.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(income.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_income.TDI_ON_CAD <- merge(income.TDI.dat, income.CAD.dat, by = "SNP") # 66SNP
write.xlsx(MVdat_income.TDI_ON_CAD, file = "Output/MVMR_IV/MVdat_income.TDI_ON_CAD.xlsx")
MVdat_income.TDI_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_income.TDI_ON_CAD.xlsx")

bx = as.matrix(MVdat_income.TDI_ON_CAD[,c("beta.income","beta.TDI")])
bxse = as.matrix(MVdat_income.TDI_ON_CAD[,c("se.income","se.TDI")])
by = as.vector(MVdat_income.TDI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_income.TDI_ON_CAD$se.CAD)

MVdatForm_income.TDI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("income","TDI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_income.TDI_ON_CAD <- mr_mvivw(MVdatForm_income.TDI_ON_CAD)
tt <- mvmr(MVdatForm_income.TDI_ON_CAD)

# MV MR-Egger
result_MV.Egger_income.TDI_ON_CAD <- mr_mvegger(MVdatForm_income.TDI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_income.TDI_ON_CAD <- mr_mvlasso(MVdatForm_income.TDI_ON_CAD)

# MV Median
result_MV.Median_income.TDI_ON_CAD <- mr_mvmedian(MVdatForm_income.TDI_ON_CAD)




F_income.TDI <- strength_mvmr(r_input = MVdatForm_income.TDI_ON_CAD, gencov = 0)


mv_hete_income.TDI <- pleiotropy_mvmr(r_input = MVdatForm_income.TDI_ON_CAD, gencov = 0)


PRESSO_income.TDI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                          BetaExposure = c("beta.income", "beta.TDI"), 
                                          SdOutcome = "se.CAD", 
                                          SdExposure = c("se.income", "se.TDI"),
                                          OUTLIERtest = TRUE, 
                                          DISTORTIONtest = TRUE, 
                                          data = MVdat_income.TDI_ON_CAD,
                                          NbDistribution = 1000, 
                                          SignifThreshold = 0.05)


###### _____CAD ~ education + occupation + income ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_occupation <- IV_occupation["SNP"]
SNP_IV_income <- IV_income["SNP"]

SNP_IV_education.occupation.income <- rbind(SNP_IV_educationStringent,SNP_IV_occupation,SNP_IV_income)# 430SNP
SNP_IV_education.occupation.income <- unique(SNP_IV_education.occupation.income) # 7SNP
SNP_IV_education.occupation.income <- clump_data(SNP_IV_education.occupation.income) # 355SNP
write.xlsx(SNP_IV_education.occupation.income, file = "Output/MVMR_IV/SNP_IV_education.occupation.income.xlsx")

exp_education_MVwith.occupation.income <- format_data(
  GWAS_education,
  snps = SNP_IV_education.occupation.income$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.occupation.income$exposure <- "education"

exp_occupation_MVwith.education.income <- format_data(
  GWAS_occupation,
  snps = SNP_IV_education.occupation.income$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value") 
exp_occupation_MVwith.education.income$outcome <- "occupation"

exp_income_MVwith.education.occupation <- format_data(
  GWAS_income,
  snps = SNP_IV_education.occupation.income$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_income_MVwith.education.occupation$outcome <- "income"

education.occupation <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.income, outcome_dat = exp_occupation_MVwith.education.income, action = 1)
education.income <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.income, outcome_dat = exp_income_MVwith.education.occupation, action = 1)

out_CAD_MV.education.occupation.income <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.occupation.income$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.occupation.income$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.income, outcome_dat = out_CAD_MV.education.occupation.income, action = 1)

education.occupation.dat <- education.occupation[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.occupation.dat) <- c("SNP", "beta.education", "beta.occupation", "se.education", "se.occupation")

education.income.dat <- education.income[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.income.dat) <- c("SNP", "beta.income", "se.income")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.occupation.dat, education.income.dat, by = "SNP")
MVdat_education.occupation.income_ON_CAD <- merge(mvmr.dat1, education.CAD.dat, by = "SNP") # 336SNP
write.xlsx(MVdat_education.occupation.income_ON_CAD, file = "Output/MVMR_IV/MVdat_education.occupation.income_ON_CAD.xlsx")
MVdat_education.occupation.income_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_education.occupation.income_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.occupation.income_ON_CAD[,c("beta.education","beta.occupation","beta.income")])
bxse = as.matrix(MVdat_education.occupation.income_ON_CAD[,c("se.education","se.occupation","se.income")])
by = as.vector(MVdat_education.occupation.income_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.occupation.income_ON_CAD$se.CAD)

MVdatForm_education.occupation.income_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","occupation","income"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.occupation.income_ON_CAD <- mr_mvivw(MVdatForm_education.occupation.income_ON_CAD)
tt <- mvmr(MVdatForm_education.occupation.income_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.occupation.income_ON_CAD <- mr_mvegger(MVdatForm_education.occupation.income_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.occupation.income_ON_CAD <- mr_mvlasso(MVdatForm_education.occupation.income_ON_CAD)

# MV Median
result_MV.Median_education.occupation.income_ON_CAD <- mr_mvmedian(MVdatForm_education.occupation.income_ON_CAD)



F_education.occupation.income <- strength_mvmr(r_input = MVdatForm_education.occupation.income_ON_CAD, gencov = 0)


mv_hete_education.occupation.income <- pleiotropy_mvmr(r_input = MVdatForm_education.occupation.income_ON_CAD, gencov = 0)


PRESSO_education.occupation.income_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.education", "beta.occupation","beta.income"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.education", "se.occupation","se.income"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_education.occupation.income_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)


###### _____CAD ~ education + occupation + TDI ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_occupation <- IV_occupation["SNP"]
SNP_IV_TDI <- IV_TDI["SNP"]

SNP_IV_education.occupation.TDI <- rbind(SNP_IV_educationStringent,SNP_IV_occupation,SNP_IV_TDI)# 395SNP
SNP_IV_education.occupation.TDI <- unique(SNP_IV_education.occupation.TDI) # 1SNP
SNP_IV_education.occupation.TDI <- clump_data(SNP_IV_education.occupation.TDI) # 352SNP
write.xlsx(SNP_IV_education.occupation.TDI, file = "Output/MVMR_IV/SNP_IV_education.occupation.TDI.xlsx")

exp_education_MVwith.occupation.TDI <- format_data(
  GWAS_education,
  snps = SNP_IV_education.occupation.TDI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.occupation.TDI$exposure <- "education"

exp_occupation_MVwith.education.TDI <- format_data(
  GWAS_occupation,
  snps = SNP_IV_education.occupation.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value") 
exp_occupation_MVwith.education.TDI$outcome <- "occupation"

exp_TDI_MVwith.education.occupation <- format_data(
  GWAS_TDI,
  snps = SNP_IV_education.occupation.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_TDI_MVwith.education.occupation$outcome <- "TDI"

education.occupation <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.TDI, outcome_dat = exp_occupation_MVwith.education.TDI, action = 1)
education.TDI <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.TDI, outcome_dat = exp_TDI_MVwith.education.occupation, action = 1)

out_CAD_MV.education.occupation.TDI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.occupation.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.occupation.TDI$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.TDI, outcome_dat = out_CAD_MV.education.occupation.TDI, action = 1)

education.occupation.dat <- education.occupation[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.occupation.dat) <- c("SNP", "beta.education", "beta.occupation", "se.education", "se.occupation")

education.TDI.dat <- education.TDI[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.TDI.dat) <- c("SNP", "beta.TDI", "se.TDI")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.occupation.dat, education.TDI.dat, by = "SNP")
MVdat_education.occupation.TDI_ON_CAD <- merge(mvmr.dat1, education.CAD.dat, by = "SNP") # 333SNP
write.xlsx(MVdat_education.occupation.TDI_ON_CAD, file = "Output/MVMR_IV/MVdat_education.occupation.TDI_ON_CAD.xlsx")
MVdat_education.occupation.TDI_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_education.occupation.TDI_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.occupation.TDI_ON_CAD[,c("beta.education","beta.occupation","beta.TDI")])
bxse = as.matrix(MVdat_education.occupation.TDI_ON_CAD[,c("se.education","se.occupation","se.TDI")])
by = as.vector(MVdat_education.occupation.TDI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.occupation.TDI_ON_CAD$se.CAD)

MVdatForm_education.occupation.TDI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","occupation","TDI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.occupation.TDI_ON_CAD <- mr_mvivw(MVdatForm_education.occupation.TDI_ON_CAD)
tt <- mvmr(MVdatForm_education.occupation.TDI_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.occupation.TDI_ON_CAD <- mr_mvegger(MVdatForm_education.occupation.TDI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.occupation.TDI_ON_CAD <- mr_mvlasso(MVdatForm_education.occupation.TDI_ON_CAD)

# MV Median
result_MV.Median_education.occupation.TDI_ON_CAD <- mr_mvmedian(MVdatForm_education.occupation.TDI_ON_CAD)



F_education.occupation.TDI <- strength_mvmr(r_input = MVdatForm_education.occupation.TDI_ON_CAD, gencov = 0)


mv_hete_education.occupation.TDI <- pleiotropy_mvmr(r_input = MVdatForm_education.occupation.TDI_ON_CAD, gencov = 0)


PRESSO_education.occupation.TDI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                       BetaExposure = c("beta.education", "beta.occupation","beta.TDI"), 
                                                       SdOutcome = "se.CAD", 
                                                       SdExposure = c("se.education", "se.occupation","se.TDI"),
                                                       OUTLIERtest = TRUE, 
                                                       DISTORTIONtest = TRUE, 
                                                       data = MVdat_education.occupation.TDI_ON_CAD,
                                                       NbDistribution = 1000, 
                                                       SignifThreshold = 0.05)


###### _____CAD ~ education + income + TDI ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_income <- IV_income["SNP"]
SNP_IV_TDI <- IV_TDI["SNP"]

SNP_IV_education.income.TDI <- rbind(SNP_IV_educationStringent,SNP_IV_income,SNP_IV_TDI)# 419SNP
SNP_IV_education.income.TDI <- unique(SNP_IV_education.income.TDI) # 6SNP
SNP_IV_education.income.TDI <- clump_data(SNP_IV_education.income.TDI) # 352SNP
write.xlsx(SNP_IV_education.income.TDI, file = "Output/MVMR_IV/SNP_IV_education.income.TDI.xlsx")

exp_education_MVwith.income.TDI <- format_data(
  GWAS_education,
  snps = SNP_IV_education.income.TDI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.income.TDI$exposure <- "education"

exp_income_MVwith.education.TDI <- format_data(
  GWAS_income,
  snps = SNP_IV_education.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_income_MVwith.education.TDI$outcome <- "income"

exp_TDI_MVwith.education.income <- format_data(
  GWAS_TDI,
  snps = SNP_IV_education.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_TDI_MVwith.education.income$outcome <- "TDI"

education.income <- harmonise_data(exposure_dat = exp_education_MVwith.income.TDI, outcome_dat = exp_income_MVwith.education.TDI, action = 2)
education.TDI <- harmonise_data(exposure_dat = exp_education_MVwith.income.TDI, outcome_dat = exp_TDI_MVwith.education.income, action = 2)

out_CAD_MV.education.income.TDI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.income.TDI$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.income.TDI, outcome_dat = out_CAD_MV.education.income.TDI, action = 2)

education.income.dat <- education.income[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.income.dat) <- c("SNP", "beta.education", "beta.income", "se.education", "se.income")

education.TDI.dat <- education.TDI[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.TDI.dat) <- c("SNP", "beta.TDI", "se.TDI")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.income.dat, education.TDI.dat, by = "SNP")
MVdat_education.income.TDI_ON_CAD <- merge(mvmr.dat1, education.CAD.dat, by = "SNP") # 351SNP
write.xlsx(MVdat_education.income.TDI_ON_CAD, file = "Output/MVMR_IV/MVdat_education.income.TDI_ON_CAD.xlsx")
MVdat_education.income.TDI_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_education.income.TDI_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.income.TDI_ON_CAD[,c("beta.education","beta.income","beta.TDI")])
bxse = as.matrix(MVdat_education.income.TDI_ON_CAD[,c("se.education","se.income","se.TDI")])
by = as.vector(MVdat_education.income.TDI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.income.TDI_ON_CAD$se.CAD)

MVdatForm_education.income.TDI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","income","TDI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.income.TDI_ON_CAD <- mr_mvivw(MVdatForm_education.income.TDI_ON_CAD)
tt <- mvmr(MVdatForm_education.income.TDI_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.income.TDI_ON_CAD <- mr_mvegger(MVdatForm_education.income.TDI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.income.TDI_ON_CAD <- mr_mvlasso(MVdatForm_education.income.TDI_ON_CAD)

# MV Median
result_MV.Median_education.income.TDI_ON_CAD <- mr_mvmedian(MVdatForm_education.income.TDI_ON_CAD)



F_education.income.TDI <- strength_mvmr(r_input = MVdatForm_education.income.TDI_ON_CAD, gencov = 0)


mv_hete_education.income.TDI <- pleiotropy_mvmr(r_input = MVdatForm_education.income.TDI_ON_CAD, gencov = 0)


PRESSO_education.income.TDI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                    BetaExposure = c("beta.education", "beta.income","beta.TDI"), 
                                                    SdOutcome = "se.CAD", 
                                                    SdExposure = c("se.education", "se.income","se.TDI"),
                                                    OUTLIERtest = TRUE, 
                                                    DISTORTIONtest = TRUE, 
                                                    data = MVdat_education.income.TDI_ON_CAD,
                                                    NbDistribution = 1000, 
                                                    SignifThreshold = 0.05)


###### _____CAD ~ occupation + income + TDI ######
SNP_IV_occupation <- IV_occupation["SNP"]
SNP_IV_income <- IV_income["SNP"]
SNP_IV_TDI <- IV_TDI["SNP"]

SNP_IV_occupation.income.TDI <- rbind(SNP_IV_occupation,SNP_IV_income,SNP_IV_TDI)# 100SNP
SNP_IV_occupation.income.TDI <- unique(SNP_IV_occupation.income.TDI) # 2SNP
SNP_IV_occupation.income.TDI <- clump_data(SNP_IV_occupation.income.TDI) # 80SNP
write.xlsx(SNP_IV_occupation.income.TDI, file = "Output/MVMR_IV/SNP_IV_occupation.income.TDI.xlsx")

exp_occupation_MVwith.income.TDI <- format_data(
  GWAS_occupation,
  snps = SNP_IV_occupation.income.TDI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value") 
exp_occupation_MVwith.income.TDI$exposure <- "occupation"

exp_income_MVwith.occupation.TDI <- format_data(
  GWAS_income,
  snps = SNP_IV_occupation.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_income_MVwith.occupation.TDI$outcome <- "income"

exp_TDI_MVwith.occupation.income <- format_data(
  GWAS_TDI,
  snps = SNP_IV_occupation.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_TDI_MVwith.occupation.income$outcome <- "TDI"

occupation.income <- harmonise_data(exposure_dat = exp_occupation_MVwith.income.TDI, outcome_dat = exp_income_MVwith.occupation.TDI, action = 1)
occupation.TDI <- harmonise_data(exposure_dat = exp_occupation_MVwith.income.TDI, outcome_dat = exp_TDI_MVwith.occupation.income, action = 1)

out_CAD_MV.occupation.income.TDI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_occupation.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.occupation.income.TDI$outcome <- "CAD"

occupation.CAD <- harmonise_data(exposure_dat = exp_occupation_MVwith.income.TDI, outcome_dat = out_CAD_MV.occupation.income.TDI, action = 1)

occupation.income.dat <- occupation.income[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(occupation.income.dat) <- c("SNP", "beta.occupation", "beta.income", "se.occupation", "se.income")

occupation.TDI.dat <- occupation.TDI[c("SNP", "beta.outcome", "se.outcome")]
colnames(occupation.TDI.dat) <- c("SNP", "beta.TDI", "se.TDI")

occupation.CAD.dat <- occupation.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(occupation.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(occupation.income.dat, occupation.TDI.dat, by = "SNP")
MVdat_occupation.income.TDI_ON_CAD <- merge(mvmr.dat1, occupation.CAD.dat, by = "SNP") # 72SNP
write.xlsx(MVdat_occupation.income.TDI_ON_CAD, file = "Output/MVMR_IV/MVdat_occupation.income.TDI_ON_CAD.xlsx")
MVdat_occupation.income.TDI_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_occupation.income.TDI_ON_CAD.xlsx")

bx = as.matrix(MVdat_occupation.income.TDI_ON_CAD[,c("beta.occupation","beta.income","beta.TDI")])
bxse = as.matrix(MVdat_occupation.income.TDI_ON_CAD[,c("se.occupation","se.income","se.TDI")])
by = as.vector(MVdat_occupation.income.TDI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_occupation.income.TDI_ON_CAD$se.CAD)

MVdatForm_occupation.income.TDI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("occupation","income","TDI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_occupation.income.TDI_ON_CAD <- mr_mvivw(MVdatForm_occupation.income.TDI_ON_CAD)
tt <- mvmr(MVdatForm_occupation.income.TDI_ON_CAD)

# MV MR-Egger
result_MV.Egger_occupation.income.TDI_ON_CAD <- mr_mvegger(MVdatForm_occupation.income.TDI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_occupation.income.TDI_ON_CAD <- mr_mvlasso(MVdatForm_occupation.income.TDI_ON_CAD)

# MV Median
result_MV.Median_occupation.income.TDI_ON_CAD <- mr_mvmedian(MVdatForm_occupation.income.TDI_ON_CAD)



F_occupation.income.TDI <- strength_mvmr(r_input = MVdatForm_occupation.income.TDI_ON_CAD, gencov = 0)


mv_hete_occupation.income.TDI <- pleiotropy_mvmr(r_input = MVdatForm_occupation.income.TDI_ON_CAD, gencov = 0)


PRESSO_occupation.income.TDI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.occupation", "beta.income","beta.TDI"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.occupation", "se.income","se.TDI"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_occupation.income.TDI_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)


###### _____CAD ~ education + occupation + income +TDI ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_occupation <- IV_occupation["SNP"]
SNP_IV_income <- IV_income["SNP"]
SNP_IV_TDI <- IV_TDI["SNP"]

SNP_IV_education.occupation.income.TDI <- rbind(SNP_IV_educationStringent,SNP_IV_occupation,SNP_IV_income,SNP_IV_TDI)# 448SNP
SNP_IV_education.occupation.income.TDI <- unique(SNP_IV_education.occupation.income.TDI) # 8SNP
SNP_IV_education.occupation.income.TDI <- clump_data(SNP_IV_education.occupation.income.TDI) # 356SNP
write.xlsx(SNP_IV_education.occupation.income.TDI, file = "Output/MVMR_IV/SNP_IV_education.occupation.income.TDI.xlsx")

exp_education_MVwith.occupation.income.TDI <- format_data(
  GWAS_education,
  snps = SNP_IV_education.occupation.income.TDI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.occupation.income.TDI$exposure <- "education"

exp_occupation_MVwith.education.income.TDI <- format_data(
  GWAS_occupation,
  snps = SNP_IV_education.occupation.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value") 
exp_occupation_MVwith.education.income.TDI$outcome <- "occupation"

exp_income_MVwith.education.occupation.TDI <- format_data(
  GWAS_income,
  snps = SNP_IV_education.occupation.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_income_MVwith.education.occupation.TDI$outcome <- "income"

exp_TDI_MVwith.education.occupation.income <- format_data(
  GWAS_TDI,
  snps = SNP_IV_education.occupation.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_TDI_MVwith.education.occupation.income$outcome <- "TDI"

education.occupation <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.income.TDI, outcome_dat = exp_occupation_MVwith.education.income.TDI, action = 1)
education.income <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.income.TDI, outcome_dat = exp_income_MVwith.education.occupation.TDI, action = 1)
education.TDI <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.income.TDI, outcome_dat = exp_TDI_MVwith.education.occupation.income, action = 1)

out_CAD_MV.education.occupation.income.TDI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.occupation.income.TDI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.occupation.income.TDI$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.occupation.income.TDI, outcome_dat = out_CAD_MV.education.occupation.income.TDI, action = 1)

education.occupation.dat <- education.occupation[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.occupation.dat) <- c("SNP", "beta.education", "beta.occupation", "se.education", "se.occupation")

education.income.dat <- education.income[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.income.dat) <- c("SNP", "beta.income", "se.income")

education.TDI.dat <- education.TDI[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.TDI.dat) <- c("SNP", "beta.TDI", "se.TDI")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.occupation.dat, education.income.dat, by = "SNP")
mvmr.dat2 <- merge(mvmr.dat1, education.TDI.dat, by = "SNP")
MVdat_education.occupation.income.TDI_ON_CAD <- merge(mvmr.dat2, education.CAD.dat, by = "SNP") # 336SNP
write.xlsx(MVdat_education.occupation.income.TDI_ON_CAD, file = "Output/MVMR_IV/MVdat_education.occupation.income.TDI_ON_CAD.xlsx")
MVdat_education.occupation.income.TDI_ON_CAD <- read.xlsx("Output/MVMR_IV/MVdat_education.occupation.income.TDI_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.occupation.income.TDI_ON_CAD[,c("beta.education","beta.occupation","beta.income","beta.TDI")])
bxse = as.matrix(MVdat_education.occupation.income.TDI_ON_CAD[,c("se.education","se.occupation","se.income","se.TDI")])
by = as.vector(MVdat_education.occupation.income.TDI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.occupation.income.TDI_ON_CAD$se.CAD)

MVdatForm_education.occupation.income.TDI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","occupation","income","TDI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.occupation.income.TDI_ON_CAD <- mr_mvivw(MVdatForm_education.occupation.income.TDI_ON_CAD)
tt <- mvmr(MVdatForm_education.occupation.income.TDI_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.occupation.income.TDI_ON_CAD <- mr_mvegger(MVdatForm_education.occupation.income.TDI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.occupation.income.TDI_ON_CAD <- mr_mvlasso(MVdatForm_education.occupation.income.TDI_ON_CAD)

# MV Median
result_MV.Median_education.occupation.income.TDI_ON_CAD <- mr_mvmedian(MVdatForm_education.occupation.income.TDI_ON_CAD)



F_education.occupation.income.TDI <- strength_mvmr(r_input = MVdatForm_education.occupation.income.TDI_ON_CAD, gencov = 0)


mv_hete_education.occupation.income.TDI <- pleiotropy_mvmr(r_input = MVdatForm_education.occupation.income.TDI_ON_CAD, gencov = 0)


PRESSO_education.occupation.income.TDI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                       BetaExposure = c("beta.education", "beta.occupation","beta.income","beta.TDI"), 
                                                       SdOutcome = "se.CAD", 
                                                       SdExposure = c("se.education", "se.occupation","se.income","se.TDI"),
                                                       OUTLIERtest = TRUE, 
                                                       DISTORTIONtest = TRUE, 
                                                       data = MVdat_education.occupation.income.TDI_ON_CAD,
                                                       NbDistribution = 1000, 
                                                       SignifThreshold = 0.05)


###### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ######
###### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ######

# mediation screening criteria I: there is a causal relationship between X and the mediator, and the impact of X on the mediator should be unidirectional
###### ___1.1）X on each M UVMR→effect of X on each M（β1） ######

#### ________M1.1 Diet_carbohydrate ####
str(IV_education)

out_education_ON_carbohydrate <- format_data(GWAS_carbohydrate,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "BETA",
                                           se_col = "SE",
                                           effect_allele_col = "effect allele",
                                           other_allele_col = "other allele",
                                           eaf_col = "Freq",
                                           pval_col = "P")
str(out_education_ON_carbohydrate)
out_education_ON_carbohydrate$outcome <- "carbohydrate"

UVdat_education_ON_carbohydrate <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_carbohydrate
)
str(UVdat_education_ON_carbohydrate)
range(UVdat_education_ON_carbohydrate$F) 
range(UVdat_education_ON_carbohydrate$pval.outcome) 
write.xlsx(UVdat_education_ON_carbohydrate, file = "Output/UVMR_IV/UVdat_education_ON_carbohydrate.xlsx")

result_education_ON_carbohydrate <- mr(UVdat_education_ON_carbohydrate)


scatter_plot_education_ON_carbohydrate <- mr_scatter_plot(result_education_ON_carbohydrate, UVdat_education_ON_carbohydrate)
scatter_plot_education_ON_carbohydrate[[1]]
ggsave(scatter_plot_education_ON_carbohydrate[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_carbohydrate.pdf", width=7, height=7)


pleiotropy_education_ON_carbohydrate <- mr_pleiotropy_test(UVdat_education_ON_carbohydrate)
str(pleiotropy_education_ON_carbohydrate)
pleiotropy_education_ON_carbohydrate <- pleiotropy_education_ON_carbohydrate[,-c(1,2)]


heterogeneity_education_ON_carbohydrate <- mr_heterogeneity(UVdat_education_ON_carbohydrate)
str(heterogeneity_education_ON_carbohydrate)
heterogeneity_education_ON_carbohydrate <- heterogeneity_education_ON_carbohydrate[,-c(1,2)]


singleSNP_education_ON_carbohydrate <- mr_singlesnp(UVdat_education_ON_carbohydrate)
funnel_plot_education_ON_carbohydrate <- mr_funnel_plot(singleSNP_education_ON_carbohydrate)
funnel_plot_education_ON_carbohydrate[[1]]
ggsave(funnel_plot_education_ON_carbohydrate[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_carbohydrate.pdf", width=7, height=7)


loo_education_ON_carbohydrate <- mr_leaveoneout(UVdat_education_ON_carbohydrate)
loo_plot_education_ON_carbohydrate <- mr_leaveoneout_plot(loo_education_ON_carbohydrate)
loo_plot_education_ON_carbohydrate[[1]]
ggsave(loo_plot_education_ON_carbohydrate[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_carbohydrate.pdf", width=7, height=21)


singleSNP_plot_education_ON_carbohydrate <- mr_forest_plot(singleSNP_education_ON_carbohydrate)
singleSNP_plot_education_ON_carbohydrate[[1]]
ggsave(singleSNP_plot_education_ON_carbohydrate[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_carbohydrate.pdf", width=7, height=21)

PRESSO_education_ON_carbohydrate <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_carbohydrate, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_carbohydrate <- generate_odds_ratios(result_education_ON_carbohydrate)
result_education_ON_carbohydrate <- result_education_ON_carbohydrate[,-c(1,2)]


#### ________M1.1 Diet_protein ####
str(IV_education)

out_education_ON_protein <- format_data(GWAS_protein,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "BETA",
                                           se_col = "SE",
                                           effect_allele_col = "effect allele",
                                           other_allele_col = "other allele",
                                           eaf_col = "Freq",
                                           pval_col = "P")
str(out_education_ON_protein)
out_education_ON_protein$outcome <- "protein"

UVdat_education_ON_protein <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_protein
)
str(UVdat_education_ON_protein)
range(UVdat_education_ON_protein$F) 
range(UVdat_education_ON_protein$pval.outcome) 
write.xlsx(UVdat_education_ON_protein, file = "Output/UVMR_IV/UVdat_education_ON_protein.xlsx")

result_education_ON_protein <- mr(UVdat_education_ON_protein)


scatter_plot_education_ON_protein <- mr_scatter_plot(result_education_ON_protein, UVdat_education_ON_protein)
scatter_plot_education_ON_protein[[1]]
ggsave(scatter_plot_education_ON_protein[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_protein.pdf", width=7, height=7)


pleiotropy_education_ON_protein <- mr_pleiotropy_test(UVdat_education_ON_protein)
str(pleiotropy_education_ON_protein)
pleiotropy_education_ON_protein <- pleiotropy_education_ON_protein[,-c(1,2)]


heterogeneity_education_ON_protein <- mr_heterogeneity(UVdat_education_ON_protein)
str(heterogeneity_education_ON_protein)
heterogeneity_education_ON_protein <- heterogeneity_education_ON_protein[,-c(1,2)]


singleSNP_education_ON_protein <- mr_singlesnp(UVdat_education_ON_protein)
funnel_plot_education_ON_protein <- mr_funnel_plot(singleSNP_education_ON_protein)
funnel_plot_education_ON_protein[[1]]
ggsave(funnel_plot_education_ON_protein[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_protein.pdf", width=7, height=7)


loo_education_ON_protein <- mr_leaveoneout(UVdat_education_ON_protein)
loo_plot_education_ON_protein <- mr_leaveoneout_plot(loo_education_ON_protein)
loo_plot_education_ON_protein[[1]]
ggsave(loo_plot_education_ON_protein[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_protein.pdf", width=7, height=21)


singleSNP_plot_education_ON_protein <- mr_forest_plot(singleSNP_education_ON_protein)
singleSNP_plot_education_ON_protein[[1]]
ggsave(singleSNP_plot_education_ON_protein[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_protein.pdf", width=7, height=21)

PRESSO_education_ON_protein <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_protein, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_protein <- generate_odds_ratios(result_education_ON_protein)
result_education_ON_protein <- result_education_ON_protein[,-c(1,2)]


#### ________M1.1 Diet_fat ####
str(IV_education)

out_education_ON_fat <- format_data(GWAS_fat,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "BETA",
                                           se_col = "SE",
                                           effect_allele_col = "effect allele",
                                           other_allele_col = "other allele",
                                           eaf_col = "Freq",
                                           pval_col = "P")
str(out_education_ON_fat)
out_education_ON_fat$outcome <- "fat"

UVdat_education_ON_fat <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_fat
)
str(UVdat_education_ON_fat)
range(UVdat_education_ON_fat$F) 
range(UVdat_education_ON_fat$pval.outcome) 
write.xlsx(UVdat_education_ON_fat, file = "Output/UVMR_IV/UVdat_education_ON_fat.xlsx")

result_education_ON_fat <- mr(UVdat_education_ON_fat)


scatter_plot_education_ON_fat <- mr_scatter_plot(result_education_ON_fat, UVdat_education_ON_fat)
scatter_plot_education_ON_fat[[1]]
ggsave(scatter_plot_education_ON_fat[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_fat.pdf", width=7, height=7)


pleiotropy_education_ON_fat <- mr_pleiotropy_test(UVdat_education_ON_fat)
str(pleiotropy_education_ON_fat)
pleiotropy_education_ON_fat <- pleiotropy_education_ON_fat[,-c(1,2)]


heterogeneity_education_ON_fat <- mr_heterogeneity(UVdat_education_ON_fat)
str(heterogeneity_education_ON_fat)
heterogeneity_education_ON_fat <- heterogeneity_education_ON_fat[,-c(1,2)]


singleSNP_education_ON_fat <- mr_singlesnp(UVdat_education_ON_fat)
funnel_plot_education_ON_fat <- mr_funnel_plot(singleSNP_education_ON_fat)
funnel_plot_education_ON_fat[[1]]
ggsave(funnel_plot_education_ON_fat[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_fat.pdf", width=7, height=7)


loo_education_ON_fat <- mr_leaveoneout(UVdat_education_ON_fat)
loo_plot_education_ON_fat <- mr_leaveoneout_plot(loo_education_ON_fat)
loo_plot_education_ON_fat[[1]]
ggsave(loo_plot_education_ON_fat[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_fat.pdf", width=7, height=21)


singleSNP_plot_education_ON_fat <- mr_forest_plot(singleSNP_education_ON_fat)
singleSNP_plot_education_ON_fat[[1]]
ggsave(singleSNP_plot_education_ON_fat[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_fat.pdf", width=7, height=21)

PRESSO_education_ON_fat <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_fat, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_fat <- generate_odds_ratios(result_education_ON_fat)
result_education_ON_fat <- result_education_ON_fat[,-c(1,2)]


#### ________M1.1 Diet_FreshFruit ####
str(IV_education)

out_education_ON_FreshFruit <- format_data(GWAS_FreshFruit,
                                    type = "outcome", 
                                    snps = IV_education$SNP, 
                                    snp_col = "SNP",
                                    beta_col = "beta",
                                    se_col = "standard_error",
                                    effect_allele_col = "effect_allele",
                                    other_allele_col = "other_allele",
                                    eaf_col = "effect_allele_frequency",
                                    pval_col = "p_value")
str(out_education_ON_FreshFruit)
out_education_ON_FreshFruit$outcome <- "FreshFruit"

UVdat_education_ON_FreshFruit <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_FreshFruit
)
str(UVdat_education_ON_FreshFruit)
range(UVdat_education_ON_FreshFruit$F) 
range(UVdat_education_ON_FreshFruit$pval.outcome) 
UVdat_education_ON_FreshFruit <- UVdat_education_ON_FreshFruit[UVdat_education_ON_FreshFruit$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_FreshFruit, file = "Output/UVMR_IV/UVdat_education_ON_FreshFruit.xlsx")

result_education_ON_FreshFruit <- mr(UVdat_education_ON_FreshFruit)


scatter_plot_education_ON_FreshFruit <- mr_scatter_plot(result_education_ON_FreshFruit, UVdat_education_ON_FreshFruit)
scatter_plot_education_ON_FreshFruit[[1]]
ggsave(scatter_plot_education_ON_FreshFruit[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_FreshFruit.pdf", width=7, height=7)


pleiotropy_education_ON_FreshFruit <- mr_pleiotropy_test(UVdat_education_ON_FreshFruit)
str(pleiotropy_education_ON_FreshFruit)
pleiotropy_education_ON_FreshFruit <- pleiotropy_education_ON_FreshFruit[,-c(1,2)]


heterogeneity_education_ON_FreshFruit <- mr_heterogeneity(UVdat_education_ON_FreshFruit)
str(heterogeneity_education_ON_FreshFruit)
heterogeneity_education_ON_FreshFruit <- heterogeneity_education_ON_FreshFruit[,-c(1,2)]


singleSNP_education_ON_FreshFruit <- mr_singlesnp(UVdat_education_ON_FreshFruit)
funnel_plot_education_ON_FreshFruit <- mr_funnel_plot(singleSNP_education_ON_FreshFruit)
funnel_plot_education_ON_FreshFruit[[1]]
ggsave(funnel_plot_education_ON_FreshFruit[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_FreshFruit.pdf", width=7, height=7)


loo_education_ON_FreshFruit <- mr_leaveoneout(UVdat_education_ON_FreshFruit)
loo_plot_education_ON_FreshFruit <- mr_leaveoneout_plot(loo_education_ON_FreshFruit)
loo_plot_education_ON_FreshFruit[[1]]
ggsave(loo_plot_education_ON_FreshFruit[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_FreshFruit.pdf", width=7, height=21)


singleSNP_plot_education_ON_FreshFruit <- mr_forest_plot(singleSNP_education_ON_FreshFruit)
singleSNP_plot_education_ON_FreshFruit[[1]]
ggsave(singleSNP_plot_education_ON_FreshFruit[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_FreshFruit.pdf", width=7, height=21)

PRESSO_education_ON_FreshFruit <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_FreshFruit, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_FreshFruit <- generate_odds_ratios(result_education_ON_FreshFruit)
result_education_ON_FreshFruit <- result_education_ON_FreshFruit[,-c(1,2)]


#### ________M1.1 Diet_DriedFruit ####
str(IV_education)

out_education_ON_DriedFruit <- format_data(GWAS_DriedFruit,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "beta",
                                           se_col = "standard_error",
                                           effect_allele_col = "effect_allele",
                                           other_allele_col = "other_allele",
                                           eaf_col = "effect_allele_frequency",
                                           pval_col = "p_value")
str(out_education_ON_DriedFruit)
out_education_ON_DriedFruit$outcome <- "DriedFruit"

UVdat_education_ON_DriedFruit <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_DriedFruit
)
str(UVdat_education_ON_DriedFruit)
range(UVdat_education_ON_DriedFruit$F) 
range(UVdat_education_ON_DriedFruit$pval.outcome) 
UVdat_education_ON_DriedFruit <- UVdat_education_ON_DriedFruit[UVdat_education_ON_DriedFruit$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_DriedFruit, file = "Output/UVMR_IV/UVdat_education_ON_DriedFruit.xlsx")

result_education_ON_DriedFruit <- mr(UVdat_education_ON_DriedFruit)


scatter_plot_education_ON_DriedFruit <- mr_scatter_plot(result_education_ON_DriedFruit, UVdat_education_ON_DriedFruit)
scatter_plot_education_ON_DriedFruit[[1]]
ggsave(scatter_plot_education_ON_DriedFruit[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_DriedFruit.pdf", width=7, height=7)


pleiotropy_education_ON_DriedFruit <- mr_pleiotropy_test(UVdat_education_ON_DriedFruit)
str(pleiotropy_education_ON_DriedFruit)
pleiotropy_education_ON_DriedFruit <- pleiotropy_education_ON_DriedFruit[,-c(1,2)]


heterogeneity_education_ON_DriedFruit <- mr_heterogeneity(UVdat_education_ON_DriedFruit)
str(heterogeneity_education_ON_DriedFruit)
heterogeneity_education_ON_DriedFruit <- heterogeneity_education_ON_DriedFruit[,-c(1,2)]


singleSNP_education_ON_DriedFruit <- mr_singlesnp(UVdat_education_ON_DriedFruit)
funnel_plot_education_ON_DriedFruit <- mr_funnel_plot(singleSNP_education_ON_DriedFruit)
funnel_plot_education_ON_DriedFruit[[1]]
ggsave(funnel_plot_education_ON_DriedFruit[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_DriedFruit.pdf", width=7, height=7)


loo_education_ON_DriedFruit <- mr_leaveoneout(UVdat_education_ON_DriedFruit)
loo_plot_education_ON_DriedFruit <- mr_leaveoneout_plot(loo_education_ON_DriedFruit)
loo_plot_education_ON_DriedFruit[[1]]
ggsave(loo_plot_education_ON_DriedFruit[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_DriedFruit.pdf", width=7, height=21)


singleSNP_plot_education_ON_DriedFruit <- mr_forest_plot(singleSNP_education_ON_DriedFruit)
singleSNP_plot_education_ON_DriedFruit[[1]]
ggsave(singleSNP_plot_education_ON_DriedFruit[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_DriedFruit.pdf", width=7, height=21)

PRESSO_education_ON_DriedFruit <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_DriedFruit, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_DriedFruit <- generate_odds_ratios(result_education_ON_DriedFruit)
result_education_ON_DriedFruit <- result_education_ON_DriedFruit[,-c(1,2)]


#### ________M1.1 Diet_RawVegetables ####
str(IV_education)

out_education_ON_RawVegetables <- format_data(GWAS_RawVegetables,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "beta",
                                           se_col = "standard_error",
                                           effect_allele_col = "effect_allele",
                                           other_allele_col = "other_allele",
                                           eaf_col = "effect_allele_frequency",
                                           pval_col = "p_value")
str(out_education_ON_RawVegetables)
out_education_ON_RawVegetables$outcome <- "RawVegetables"

UVdat_education_ON_RawVegetables <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_RawVegetables
)
str(UVdat_education_ON_RawVegetables)
range(UVdat_education_ON_RawVegetables$F) 
range(UVdat_education_ON_RawVegetables$pval.outcome) 
UVdat_education_ON_RawVegetables <- UVdat_education_ON_RawVegetables[UVdat_education_ON_RawVegetables$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_RawVegetables, file = "Output/UVMR_IV/UVdat_education_ON_RawVegetables.xlsx")

result_education_ON_RawVegetables <- mr(UVdat_education_ON_RawVegetables)


scatter_plot_education_ON_RawVegetables <- mr_scatter_plot(result_education_ON_RawVegetables, UVdat_education_ON_RawVegetables)
scatter_plot_education_ON_RawVegetables[[1]]
ggsave(scatter_plot_education_ON_RawVegetables[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_RawVegetables.pdf", width=7, height=7)


pleiotropy_education_ON_RawVegetables <- mr_pleiotropy_test(UVdat_education_ON_RawVegetables)
str(pleiotropy_education_ON_RawVegetables)
pleiotropy_education_ON_RawVegetables <- pleiotropy_education_ON_RawVegetables[,-c(1,2)]


heterogeneity_education_ON_RawVegetables <- mr_heterogeneity(UVdat_education_ON_RawVegetables)
str(heterogeneity_education_ON_RawVegetables)
heterogeneity_education_ON_RawVegetables <- heterogeneity_education_ON_RawVegetables[,-c(1,2)]


singleSNP_education_ON_RawVegetables <- mr_singlesnp(UVdat_education_ON_RawVegetables)
funnel_plot_education_ON_RawVegetables <- mr_funnel_plot(singleSNP_education_ON_RawVegetables)
funnel_plot_education_ON_RawVegetables[[1]]
ggsave(funnel_plot_education_ON_RawVegetables[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_RawVegetables.pdf", width=7, height=7)


loo_education_ON_RawVegetables <- mr_leaveoneout(UVdat_education_ON_RawVegetables)
loo_plot_education_ON_RawVegetables <- mr_leaveoneout_plot(loo_education_ON_RawVegetables)
loo_plot_education_ON_RawVegetables[[1]]
ggsave(loo_plot_education_ON_RawVegetables[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_RawVegetables.pdf", width=7, height=21)


singleSNP_plot_education_ON_RawVegetables <- mr_forest_plot(singleSNP_education_ON_RawVegetables)
singleSNP_plot_education_ON_RawVegetables[[1]]
ggsave(singleSNP_plot_education_ON_RawVegetables[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_RawVegetables.pdf", width=7, height=21)

PRESSO_education_ON_RawVegetables <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_RawVegetables, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_RawVegetables <- generate_odds_ratios(result_education_ON_RawVegetables)
result_education_ON_RawVegetables <- result_education_ON_RawVegetables[,-c(1,2)]


#### ________M1.1 Diet_CookedVegetables ####
str(IV_education)

out_education_ON_CookedVegetables <- format_data(GWAS_CookedVegetables,
                                              type = "outcome", 
                                              snps = IV_education$SNP, 
                                              snp_col = "SNP",
                                              beta_col = "beta",
                                              se_col = "standard_error",
                                              effect_allele_col = "effect_allele",
                                              other_allele_col = "other_allele",
                                              eaf_col = "effect_allele_frequency",
                                              pval_col = "p_value")
str(out_education_ON_CookedVegetables)
out_education_ON_CookedVegetables$outcome <- "CookedVegetables"

UVdat_education_ON_CookedVegetables <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_CookedVegetables
)
str(UVdat_education_ON_CookedVegetables)
range(UVdat_education_ON_CookedVegetables$F) 
range(UVdat_education_ON_CookedVegetables$pval.outcome) 
UVdat_education_ON_CookedVegetables <- UVdat_education_ON_CookedVegetables[UVdat_education_ON_CookedVegetables$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_CookedVegetables, file = "Output/UVMR_IV/UVdat_education_ON_CookedVegetables.xlsx")

result_education_ON_CookedVegetables <- mr(UVdat_education_ON_CookedVegetables)


scatter_plot_education_ON_CookedVegetables <- mr_scatter_plot(result_education_ON_CookedVegetables, UVdat_education_ON_CookedVegetables)
scatter_plot_education_ON_CookedVegetables[[1]]
ggsave(scatter_plot_education_ON_CookedVegetables[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_CookedVegetables.pdf", width=7, height=7)


pleiotropy_education_ON_CookedVegetables <- mr_pleiotropy_test(UVdat_education_ON_CookedVegetables)
str(pleiotropy_education_ON_CookedVegetables)
pleiotropy_education_ON_CookedVegetables <- pleiotropy_education_ON_CookedVegetables[,-c(1,2)]


heterogeneity_education_ON_CookedVegetables <- mr_heterogeneity(UVdat_education_ON_CookedVegetables)
str(heterogeneity_education_ON_CookedVegetables)
heterogeneity_education_ON_CookedVegetables <- heterogeneity_education_ON_CookedVegetables[,-c(1,2)]


singleSNP_education_ON_CookedVegetables <- mr_singlesnp(UVdat_education_ON_CookedVegetables)
funnel_plot_education_ON_CookedVegetables <- mr_funnel_plot(singleSNP_education_ON_CookedVegetables)
funnel_plot_education_ON_CookedVegetables[[1]]
ggsave(funnel_plot_education_ON_CookedVegetables[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_CookedVegetables.pdf", width=7, height=7)


loo_education_ON_CookedVegetables <- mr_leaveoneout(UVdat_education_ON_CookedVegetables)
loo_plot_education_ON_CookedVegetables <- mr_leaveoneout_plot(loo_education_ON_CookedVegetables)
loo_plot_education_ON_CookedVegetables[[1]]
ggsave(loo_plot_education_ON_CookedVegetables[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_CookedVegetables.pdf", width=7, height=21)


singleSNP_plot_education_ON_CookedVegetables <- mr_forest_plot(singleSNP_education_ON_CookedVegetables)
singleSNP_plot_education_ON_CookedVegetables[[1]]
ggsave(singleSNP_plot_education_ON_CookedVegetables[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_CookedVegetables.pdf", width=7, height=21)

PRESSO_education_ON_CookedVegetables <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_CookedVegetables, NbDistribution = 10000,  
                                               SignifThreshold = 0.05)

result_education_ON_CookedVegetables <- generate_odds_ratios(result_education_ON_CookedVegetables)
result_education_ON_CookedVegetables <- result_education_ON_CookedVegetables[,-c(1,2)]


#### ________M1.1 Diet_salt ####
str(IV_education)

out_education_ON_salt <- format_data(GWAS_salt,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "beta",
                                           se_col = "standard_error",
                                           effect_allele_col = "effect_allele",
                                           other_allele_col = "other_allele",
                                           eaf_col = "effect_allele_frequency",
                                           pval_col = "p_value")
str(out_education_ON_salt)
out_education_ON_salt$outcome <- "salt"

UVdat_education_ON_salt <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_salt
)
str(UVdat_education_ON_salt)
range(UVdat_education_ON_salt$F) 
range(UVdat_education_ON_salt$pval.outcome) 
UVdat_education_ON_salt <- UVdat_education_ON_salt[UVdat_education_ON_salt$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_salt, file = "Output/UVMR_IV/UVdat_education_ON_salt.xlsx")

result_education_ON_salt <- mr(UVdat_education_ON_salt)


scatter_plot_education_ON_salt <- mr_scatter_plot(result_education_ON_salt, UVdat_education_ON_salt)
scatter_plot_education_ON_salt[[1]]
ggsave(scatter_plot_education_ON_salt[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_salt.pdf", width=7, height=7)


pleiotropy_education_ON_salt <- mr_pleiotropy_test(UVdat_education_ON_salt)
str(pleiotropy_education_ON_salt)
pleiotropy_education_ON_salt <- pleiotropy_education_ON_salt[,-c(1,2)]


heterogeneity_education_ON_salt <- mr_heterogeneity(UVdat_education_ON_salt)
str(heterogeneity_education_ON_salt)
heterogeneity_education_ON_salt <- heterogeneity_education_ON_salt[,-c(1,2)]


singleSNP_education_ON_salt <- mr_singlesnp(UVdat_education_ON_salt)
funnel_plot_education_ON_salt <- mr_funnel_plot(singleSNP_education_ON_salt)
funnel_plot_education_ON_salt[[1]]
ggsave(funnel_plot_education_ON_salt[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_salt.pdf", width=7, height=7)


loo_education_ON_salt <- mr_leaveoneout(UVdat_education_ON_salt)
loo_plot_education_ON_salt <- mr_leaveoneout_plot(loo_education_ON_salt)
loo_plot_education_ON_salt[[1]]
ggsave(loo_plot_education_ON_salt[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_salt.pdf", width=7, height=21)


singleSNP_plot_education_ON_salt <- mr_forest_plot(singleSNP_education_ON_salt)
singleSNP_plot_education_ON_salt[[1]]
ggsave(singleSNP_plot_education_ON_salt[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_salt.pdf", width=7, height=21)

PRESSO_education_ON_salt <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_salt, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_salt <- generate_odds_ratios(result_education_ON_salt)
result_education_ON_salt <- result_education_ON_salt[,-c(1,2)]


#### ________M1.1 Diet_processed ####
str(IV_education)

out_education_ON_processed <- format_data(GWAS_processed,
                                     type = "outcome", 
                                     snps = IV_education$SNP, 
                                     snp_col = "SNP",
                                     beta_col = "beta",
                                     se_col = "standard_error",
                                     effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",
                                     eaf_col = "effect_allele_frequency",
                                     pval_col = "p_value")
str(out_education_ON_processed)
out_education_ON_processed$outcome <- "processed"

UVdat_education_ON_processed <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_processed
)
str(UVdat_education_ON_processed)
range(UVdat_education_ON_processed$F) 
range(UVdat_education_ON_processed$pval.outcome) 
UVdat_education_ON_processed <- UVdat_education_ON_processed[UVdat_education_ON_processed$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_processed, file = "Output/UVMR_IV/UVdat_education_ON_processed.xlsx")

result_education_ON_processed <- mr(UVdat_education_ON_processed)


scatter_plot_education_ON_processed <- mr_scatter_plot(result_education_ON_processed, UVdat_education_ON_processed)
scatter_plot_education_ON_processed[[1]]
ggsave(scatter_plot_education_ON_processed[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_processed.pdf", width=7, height=7)


pleiotropy_education_ON_processed <- mr_pleiotropy_test(UVdat_education_ON_processed)
str(pleiotropy_education_ON_processed)
pleiotropy_education_ON_processed <- pleiotropy_education_ON_processed[,-c(1,2)]


heterogeneity_education_ON_processed <- mr_heterogeneity(UVdat_education_ON_processed)
str(heterogeneity_education_ON_processed)
heterogeneity_education_ON_processed <- heterogeneity_education_ON_processed[,-c(1,2)]


singleSNP_education_ON_processed <- mr_singlesnp(UVdat_education_ON_processed)
funnel_plot_education_ON_processed <- mr_funnel_plot(singleSNP_education_ON_processed)
funnel_plot_education_ON_processed[[1]]
ggsave(funnel_plot_education_ON_processed[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_processed.pdf", width=7, height=7)


loo_education_ON_processed <- mr_leaveoneout(UVdat_education_ON_processed)
loo_plot_education_ON_processed <- mr_leaveoneout_plot(loo_education_ON_processed)
loo_plot_education_ON_processed[[1]]
ggsave(loo_plot_education_ON_processed[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_processed.pdf", width=7, height=21)


singleSNP_plot_education_ON_processed <- mr_forest_plot(singleSNP_education_ON_processed)
singleSNP_plot_education_ON_processed[[1]]
ggsave(singleSNP_plot_education_ON_processed[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_processed.pdf", width=7, height=21)

PRESSO_education_ON_processed <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_processed, NbDistribution = 10000,  
                                      SignifThreshold = 0.05)

result_education_ON_processed <- generate_odds_ratios(result_education_ON_processed)
result_education_ON_processed <- result_education_ON_processed[,-c(1,2)]


#### ________M1.1 Diet_OilyFish ####
str(IV_education)

out_education_ON_OilyFish <- format_data(GWAS_OilyFish,
                                          type = "outcome", 
                                          snps = IV_education$SNP, 
                                          snp_col = "SNP",
                                          beta_col = "beta",
                                          se_col = "standard_error",
                                          effect_allele_col = "effect_allele",
                                          other_allele_col = "other_allele",
                                          eaf_col = "effect_allele_frequency",
                                          pval_col = "p_value")
str(out_education_ON_OilyFish)
out_education_ON_OilyFish$outcome <- "OilyFish"

UVdat_education_ON_OilyFish <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_OilyFish
)
str(UVdat_education_ON_OilyFish)
range(UVdat_education_ON_OilyFish$F) 
range(UVdat_education_ON_OilyFish$pval.outcome) 
UVdat_education_ON_OilyFish <- UVdat_education_ON_OilyFish[UVdat_education_ON_OilyFish$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_OilyFish, file = "Output/UVMR_IV/UVdat_education_ON_OilyFish.xlsx")

result_education_ON_OilyFish <- mr(UVdat_education_ON_OilyFish)


scatter_plot_education_ON_OilyFish <- mr_scatter_plot(result_education_ON_OilyFish, UVdat_education_ON_OilyFish)
scatter_plot_education_ON_OilyFish[[1]]
ggsave(scatter_plot_education_ON_OilyFish[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_OilyFish.pdf", width=7, height=7)


pleiotropy_education_ON_OilyFish <- mr_pleiotropy_test(UVdat_education_ON_OilyFish)
str(pleiotropy_education_ON_OilyFish)
pleiotropy_education_ON_OilyFish <- pleiotropy_education_ON_OilyFish[,-c(1,2)]


heterogeneity_education_ON_OilyFish <- mr_heterogeneity(UVdat_education_ON_OilyFish)
str(heterogeneity_education_ON_OilyFish)
heterogeneity_education_ON_OilyFish <- heterogeneity_education_ON_OilyFish[,-c(1,2)]


singleSNP_education_ON_OilyFish <- mr_singlesnp(UVdat_education_ON_OilyFish)
funnel_plot_education_ON_OilyFish <- mr_funnel_plot(singleSNP_education_ON_OilyFish)
funnel_plot_education_ON_OilyFish[[1]]
ggsave(funnel_plot_education_ON_OilyFish[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_OilyFish.pdf", width=7, height=7)


loo_education_ON_OilyFish <- mr_leaveoneout(UVdat_education_ON_OilyFish)
loo_plot_education_ON_OilyFish <- mr_leaveoneout_plot(loo_education_ON_OilyFish)
loo_plot_education_ON_OilyFish[[1]]
ggsave(loo_plot_education_ON_OilyFish[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_OilyFish.pdf", width=7, height=21)


singleSNP_plot_education_ON_OilyFish <- mr_forest_plot(singleSNP_education_ON_OilyFish)
singleSNP_plot_education_ON_OilyFish[[1]]
ggsave(singleSNP_plot_education_ON_OilyFish[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_OilyFish.pdf", width=7, height=21)

PRESSO_education_ON_OilyFish <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                           OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_OilyFish, NbDistribution = 10000,  
                                           SignifThreshold = 0.05)

result_education_ON_OilyFish <- generate_odds_ratios(result_education_ON_OilyFish)
result_education_ON_OilyFish <- result_education_ON_OilyFish[,-c(1,2)]


#### ________M1.1 Diet_wholegrain ####
str(IV_education)

out_education_ON_wholegrain <- format_data(GWAS_wholegrain,
                                         type = "outcome", 
                                         snps = IV_education$SNP, 
                                         snp_col = "SNP",
                                         beta_col = "beta",
                                         se_col = "standard_error",
                                         effect_allele_col = "effect_allele",
                                         other_allele_col = "other_allele",
                                         eaf_col = "effect_allele_frequency",
                                         pval_col = "p_value")
str(out_education_ON_wholegrain)
out_education_ON_wholegrain$outcome <- "wholegrain"

UVdat_education_ON_wholegrain <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_wholegrain
)
str(UVdat_education_ON_wholegrain)
range(UVdat_education_ON_wholegrain$F) 
range(UVdat_education_ON_wholegrain$pval.outcome) 
UVdat_education_ON_wholegrain <- UVdat_education_ON_wholegrain[UVdat_education_ON_wholegrain$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_wholegrain, file = "Output/UVMR_IV/UVdat_education_ON_wholegrain.xlsx")

result_education_ON_wholegrain <- mr(UVdat_education_ON_wholegrain)


scatter_plot_education_ON_wholegrain <- mr_scatter_plot(result_education_ON_wholegrain, UVdat_education_ON_wholegrain)
scatter_plot_education_ON_wholegrain[[1]]
ggsave(scatter_plot_education_ON_wholegrain[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_wholegrain.pdf", width=7, height=7)


pleiotropy_education_ON_wholegrain <- mr_pleiotropy_test(UVdat_education_ON_wholegrain)
str(pleiotropy_education_ON_wholegrain)
pleiotropy_education_ON_wholegrain <- pleiotropy_education_ON_wholegrain[,-c(1,2)]


heterogeneity_education_ON_wholegrain <- mr_heterogeneity(UVdat_education_ON_wholegrain)
str(heterogeneity_education_ON_wholegrain)
heterogeneity_education_ON_wholegrain <- heterogeneity_education_ON_wholegrain[,-c(1,2)]


singleSNP_education_ON_wholegrain <- mr_singlesnp(UVdat_education_ON_wholegrain)
funnel_plot_education_ON_wholegrain <- mr_funnel_plot(singleSNP_education_ON_wholegrain)
funnel_plot_education_ON_wholegrain[[1]]
ggsave(funnel_plot_education_ON_wholegrain[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_wholegrain.pdf", width=7, height=7)


loo_education_ON_wholegrain <- mr_leaveoneout(UVdat_education_ON_wholegrain)
loo_plot_education_ON_wholegrain <- mr_leaveoneout_plot(loo_education_ON_wholegrain)
loo_plot_education_ON_wholegrain[[1]]
ggsave(loo_plot_education_ON_wholegrain[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_wholegrain.pdf", width=7, height=21)


singleSNP_plot_education_ON_wholegrain <- mr_forest_plot(singleSNP_education_ON_wholegrain)
singleSNP_plot_education_ON_wholegrain[[1]]
ggsave(singleSNP_plot_education_ON_wholegrain[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_wholegrain.pdf", width=7, height=21)

PRESSO_education_ON_wholegrain <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_wholegrain, NbDistribution = 10000,  
                                          SignifThreshold = 0.05)

result_education_ON_wholegrain <- generate_odds_ratios(result_education_ON_wholegrain)
result_education_ON_wholegrain <- result_education_ON_wholegrain[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceOverallActivity ####
str(IV_education)

out_education_ON_DeviceOverallActivity <- format_data(GWAS_DeviceOverallActivity,
                                   type = "outcome", 
                                   snps = IV_education$SNP, 
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "ALLELE1",
                                   other_allele_col = "ALLELE0",
                                   eaf_col = "A1FREQ",
                                   pval_col = "P_BOLT_LMM_INF")
str(out_education_ON_DeviceOverallActivity)
out_education_ON_DeviceOverallActivity$outcome <- "DeviceOverallActivity"

UVdat_education_ON_DeviceOverallActivity <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_DeviceOverallActivity
)
str(UVdat_education_ON_DeviceOverallActivity)
range(UVdat_education_ON_DeviceOverallActivity$pval.outcome)
range(UVdat_education_ON_DeviceOverallActivity$F)
UVdat_education_ON_DeviceOverallActivity <- UVdat_education_ON_DeviceOverallActivity[UVdat_education_ON_DeviceOverallActivity$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_DeviceOverallActivity, file = "Output/UVMR_IV/UVdat_education_ON_DeviceOverallActivity.xlsx")

result_education_ON_DeviceOverallActivity <- mr(UVdat_education_ON_DeviceOverallActivity)


scatter_plot_education_ON_DeviceOverallActivity <- mr_scatter_plot(result_education_ON_DeviceOverallActivity, UVdat_education_ON_DeviceOverallActivity)
scatter_plot_education_ON_DeviceOverallActivity[[1]]
ggsave(scatter_plot_education_ON_DeviceOverallActivity[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_DeviceOverallActivity.pdf", width=7, height=7)


pleiotropy_education_ON_DeviceOverallActivity <- mr_pleiotropy_test(UVdat_education_ON_DeviceOverallActivity)
str(pleiotropy_education_ON_DeviceOverallActivity)
pleiotropy_education_ON_DeviceOverallActivity <- pleiotropy_education_ON_DeviceOverallActivity[,-c(1,2)]


heterogeneity_education_ON_DeviceOverallActivity <- mr_heterogeneity(UVdat_education_ON_DeviceOverallActivity)
str(heterogeneity_education_ON_DeviceOverallActivity)
heterogeneity_education_ON_DeviceOverallActivity <- heterogeneity_education_ON_DeviceOverallActivity[,-c(1,2)]


singleSNP_education_ON_DeviceOverallActivity <- mr_singlesnp(UVdat_education_ON_DeviceOverallActivity)
funnel_plot_education_ON_DeviceOverallActivity <- mr_funnel_plot(singleSNP_education_ON_DeviceOverallActivity)
funnel_plot_education_ON_DeviceOverallActivity[[1]]
ggsave(funnel_plot_education_ON_DeviceOverallActivity[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_DeviceOverallActivity.pdf", width=7, height=7)


loo_education_ON_DeviceOverallActivity <- mr_leaveoneout(UVdat_education_ON_DeviceOverallActivity)
loo_plot_education_ON_DeviceOverallActivity <- mr_leaveoneout_plot(loo_education_ON_DeviceOverallActivity)
loo_plot_education_ON_DeviceOverallActivity[[1]]
ggsave(loo_plot_education_ON_DeviceOverallActivity[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_DeviceOverallActivity.pdf", width=7, height=21)


singleSNP_plot_education_ON_DeviceOverallActivity <- mr_forest_plot(singleSNP_education_ON_DeviceOverallActivity)
singleSNP_plot_education_ON_DeviceOverallActivity[[1]]
ggsave(singleSNP_plot_education_ON_DeviceOverallActivity[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_DeviceOverallActivity.pdf", width=7, height=21)

PRESSO_education_ON_DeviceOverallActivity <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_DeviceOverallActivity, NbDistribution = 10000,  
                                    SignifThreshold = 0.05)

result_education_ON_DeviceOverallActivity <- generate_odds_ratios(result_education_ON_DeviceOverallActivity)
result_education_ON_DeviceOverallActivity <- result_education_ON_DeviceOverallActivity[,-c(1,2)]


#### ________M1.2 Physical activity_PA ####
str(IV_education)

out_education_ON_PA <- format_data(GWAS_PA,
                                   type = "outcome", 
                                   snps = IV_education$SNP, 
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "effect allele",
                                   other_allele_col = "other allele",
                                   eaf_col = "Freq",
                                   pval_col = "P")
str(out_education_ON_PA)
out_education_ON_PA$outcome <- "MVPA"

UVdat_education_ON_PA <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_PA
)
str(UVdat_education_ON_PA)
range(UVdat_education_ON_PA$pval.outcome)
range(UVdat_education_ON_PA$F)
UVdat_education_ON_PA <- UVdat_education_ON_PA[UVdat_education_ON_PA$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_PA, file = "Output/UVMR_IV/UVdat_education_ON_PA.xlsx")

result_education_ON_PA <- mr(UVdat_education_ON_PA)


scatter_plot_education_ON_PA <- mr_scatter_plot(result_education_ON_PA, UVdat_education_ON_PA)
scatter_plot_education_ON_PA[[1]]
ggsave(scatter_plot_education_ON_PA[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_PA.pdf", width=7, height=7)


pleiotropy_education_ON_PA <- mr_pleiotropy_test(UVdat_education_ON_PA)
str(pleiotropy_education_ON_PA)
pleiotropy_education_ON_PA <- pleiotropy_education_ON_PA[,-c(1,2)]


heterogeneity_education_ON_PA <- mr_heterogeneity(UVdat_education_ON_PA)
str(heterogeneity_education_ON_PA)
heterogeneity_education_ON_PA <- heterogeneity_education_ON_PA[,-c(1,2)]


singleSNP_education_ON_PA <- mr_singlesnp(UVdat_education_ON_PA)
funnel_plot_education_ON_PA <- mr_funnel_plot(singleSNP_education_ON_PA)
funnel_plot_education_ON_PA[[1]]
ggsave(funnel_plot_education_ON_PA[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_PA.pdf", width=7, height=7)


loo_education_ON_PA <- mr_leaveoneout(UVdat_education_ON_PA)
loo_plot_education_ON_PA <- mr_leaveoneout_plot(loo_education_ON_PA)
loo_plot_education_ON_PA[[1]]
ggsave(loo_plot_education_ON_PA[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_PA.pdf", width=7, height=21)


singleSNP_plot_education_ON_PA <- mr_forest_plot(singleSNP_education_ON_PA)
singleSNP_plot_education_ON_PA[[1]]
ggsave(singleSNP_plot_education_ON_PA[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_PA.pdf", width=7, height=21)

PRESSO_education_ON_PA <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_PA, NbDistribution = 10000,  
                                    SignifThreshold = 0.05)

result_education_ON_PA <- generate_odds_ratios(result_education_ON_PA)
result_education_ON_PA <- result_education_ON_PA[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceModerate ####
str(IV_education)

out_education_ON_DeviceModerate <- format_data(GWAS_DeviceModerate,
                                                      type = "outcome", 
                                                      snps = IV_education$SNP, 
                                                      snp_col = "SNP",
                                                      beta_col = "BETA",
                                                      se_col = "SE",
                                                      effect_allele_col = "ALLELE1",
                                                      other_allele_col = "ALLELE0",
                                                      eaf_col = "A1FREQ",
                                                      pval_col = "P_BOLT_LMM_INF")
str(out_education_ON_DeviceModerate)
out_education_ON_DeviceModerate$outcome <- "DeviceModerate"

UVdat_education_ON_DeviceModerate <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_DeviceModerate
)
str(UVdat_education_ON_DeviceModerate)
range(UVdat_education_ON_DeviceModerate$pval.outcome)
range(UVdat_education_ON_DeviceModerate$F)
write.xlsx(UVdat_education_ON_DeviceModerate, file = "Output/UVMR_IV/UVdat_education_ON_DeviceModerate.xlsx")

result_education_ON_DeviceModerate <- mr(UVdat_education_ON_DeviceModerate)


scatter_plot_education_ON_DeviceModerate <- mr_scatter_plot(result_education_ON_DeviceModerate, UVdat_education_ON_DeviceModerate)
scatter_plot_education_ON_DeviceModerate[[1]]
ggsave(scatter_plot_education_ON_DeviceModerate[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_DeviceModerate.pdf", width=7, height=7)


pleiotropy_education_ON_DeviceModerate <- mr_pleiotropy_test(UVdat_education_ON_DeviceModerate)
str(pleiotropy_education_ON_DeviceModerate)
pleiotropy_education_ON_DeviceModerate <- pleiotropy_education_ON_DeviceModerate[,-c(1,2)]


heterogeneity_education_ON_DeviceModerate <- mr_heterogeneity(UVdat_education_ON_DeviceModerate)
str(heterogeneity_education_ON_DeviceModerate)
heterogeneity_education_ON_DeviceModerate <- heterogeneity_education_ON_DeviceModerate[,-c(1,2)]


singleSNP_education_ON_DeviceModerate <- mr_singlesnp(UVdat_education_ON_DeviceModerate)
funnel_plot_education_ON_DeviceModerate <- mr_funnel_plot(singleSNP_education_ON_DeviceModerate)
funnel_plot_education_ON_DeviceModerate[[1]]
ggsave(funnel_plot_education_ON_DeviceModerate[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_DeviceModerate.pdf", width=7, height=7)


loo_education_ON_DeviceModerate <- mr_leaveoneout(UVdat_education_ON_DeviceModerate)
loo_plot_education_ON_DeviceModerate <- mr_leaveoneout_plot(loo_education_ON_DeviceModerate)
loo_plot_education_ON_DeviceModerate[[1]]
ggsave(loo_plot_education_ON_DeviceModerate[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_DeviceModerate.pdf", width=7, height=21)


singleSNP_plot_education_ON_DeviceModerate <- mr_forest_plot(singleSNP_education_ON_DeviceModerate)
singleSNP_plot_education_ON_DeviceModerate[[1]]
ggsave(singleSNP_plot_education_ON_DeviceModerate[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_DeviceModerate.pdf", width=7, height=21)

PRESSO_education_ON_DeviceModerate <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_DeviceModerate, NbDistribution = 10000,  
                                                       SignifThreshold = 0.05)

result_education_ON_DeviceModerate <- generate_odds_ratios(result_education_ON_DeviceModerate)
result_education_ON_DeviceModerate <- result_education_ON_DeviceModerate[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceSedentary ####
str(IV_education)

out_education_ON_DeviceSedentary <- format_data(GWAS_DeviceSedentary,
                                                      type = "outcome", 
                                                      snps = IV_education$SNP, 
                                                      snp_col = "SNP",
                                                      beta_col = "BETA",
                                                      se_col = "SE",
                                                      effect_allele_col = "ALLELE1",
                                                      other_allele_col = "ALLELE0",
                                                      eaf_col = "A1FREQ",
                                                      pval_col = "P_BOLT_LMM_INF")
str(out_education_ON_DeviceSedentary)
out_education_ON_DeviceSedentary$outcome <- "DeviceSedentary"

UVdat_education_ON_DeviceSedentary <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_DeviceSedentary
)
str(UVdat_education_ON_DeviceSedentary)
range(UVdat_education_ON_DeviceSedentary$pval.outcome)
range(UVdat_education_ON_DeviceSedentary$F)
write.xlsx(UVdat_education_ON_DeviceSedentary, file = "Output/UVMR_IV/UVdat_education_ON_DeviceSedentary.xlsx")

result_education_ON_DeviceSedentary <- mr(UVdat_education_ON_DeviceSedentary)


scatter_plot_education_ON_DeviceSedentary <- mr_scatter_plot(result_education_ON_DeviceSedentary, UVdat_education_ON_DeviceSedentary)
scatter_plot_education_ON_DeviceSedentary[[1]]
ggsave(scatter_plot_education_ON_DeviceSedentary[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_DeviceSedentary.pdf", width=7, height=7)


pleiotropy_education_ON_DeviceSedentary <- mr_pleiotropy_test(UVdat_education_ON_DeviceSedentary)
str(pleiotropy_education_ON_DeviceSedentary)
pleiotropy_education_ON_DeviceSedentary <- pleiotropy_education_ON_DeviceSedentary[,-c(1,2)]


heterogeneity_education_ON_DeviceSedentary <- mr_heterogeneity(UVdat_education_ON_DeviceSedentary)
str(heterogeneity_education_ON_DeviceSedentary)
heterogeneity_education_ON_DeviceSedentary <- heterogeneity_education_ON_DeviceSedentary[,-c(1,2)]


singleSNP_education_ON_DeviceSedentary <- mr_singlesnp(UVdat_education_ON_DeviceSedentary)
funnel_plot_education_ON_DeviceSedentary <- mr_funnel_plot(singleSNP_education_ON_DeviceSedentary)
funnel_plot_education_ON_DeviceSedentary[[1]]
ggsave(funnel_plot_education_ON_DeviceSedentary[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_DeviceSedentary.pdf", width=7, height=7)


loo_education_ON_DeviceSedentary <- mr_leaveoneout(UVdat_education_ON_DeviceSedentary)
loo_plot_education_ON_DeviceSedentary <- mr_leaveoneout_plot(loo_education_ON_DeviceSedentary)
loo_plot_education_ON_DeviceSedentary[[1]]
ggsave(loo_plot_education_ON_DeviceSedentary[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_DeviceSedentary.pdf", width=7, height=21)


singleSNP_plot_education_ON_DeviceSedentary <- mr_forest_plot(singleSNP_education_ON_DeviceSedentary)
singleSNP_plot_education_ON_DeviceSedentary[[1]]
ggsave(singleSNP_plot_education_ON_DeviceSedentary[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_DeviceSedentary.pdf", width=7, height=21)

PRESSO_education_ON_DeviceSedentary <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_DeviceSedentary, NbDistribution = 10000,  
                                                       SignifThreshold = 0.05)

result_education_ON_DeviceSedentary <- generate_odds_ratios(result_education_ON_DeviceSedentary)
result_education_ON_DeviceSedentary <- result_education_ON_DeviceSedentary[,-c(1,2)]


#### ________M1.2 Physical activity_LST ####
str(IV_education)

out_education_ON_LST <- format_data(GWAS_LST,
                                    type = "outcome", 
                                    snps = IV_education$SNP, 
                                    snp_col = "SNP",
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    effect_allele_col = "effect allele",
                                    other_allele_col = "other allele",
                                    eaf_col = "Freq",
                                    pval_col = "P")
str(out_education_ON_LST)
out_education_ON_LST$outcome <- "LST"

UVdat_education_ON_LST <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_LST
)
str(UVdat_education_ON_LST)
range(UVdat_education_ON_LST$F) # [1]  28.22258 373.77680
range(UVdat_education_ON_LST$pval.outcome) # [1] 2.727e-19 9.997e-01
UVdat_education_ON_LST <- UVdat_education_ON_LST[UVdat_education_ON_LST$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_LST, file = "Output/UVMR_IV/UVdat_education_ON_LST.xlsx")

result_education_ON_LST <- mr(UVdat_education_ON_LST)


scatter_plot_education_ON_LST <- mr_scatter_plot(result_education_ON_LST, UVdat_education_ON_LST)
scatter_plot_education_ON_LST[[1]]
ggsave(scatter_plot_education_ON_LST[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_LST.pdf", width=7, height=7)


pleiotropy_education_ON_LST <- mr_pleiotropy_test(UVdat_education_ON_LST)
str(pleiotropy_education_ON_LST)
pleiotropy_education_ON_LST <- pleiotropy_education_ON_LST[,-c(1,2)]


heterogeneity_education_ON_LST <- mr_heterogeneity(UVdat_education_ON_LST)
str(heterogeneity_education_ON_LST)
heterogeneity_education_ON_LST <- heterogeneity_education_ON_LST[,-c(1,2)]


singleSNP_education_ON_LST <- mr_singlesnp(UVdat_education_ON_LST)
funnel_plot_education_ON_LST <- mr_funnel_plot(singleSNP_education_ON_LST)
funnel_plot_education_ON_LST[[1]]
ggsave(funnel_plot_education_ON_LST[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_LST.pdf", width=7, height=7)


loo_education_ON_LST <- mr_leaveoneout(UVdat_education_ON_LST)
loo_plot_education_ON_LST <- mr_leaveoneout_plot(loo_education_ON_LST)
loo_plot_education_ON_LST[[1]]
ggsave(loo_plot_education_ON_LST[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_LST.pdf", width=7, height=21)


singleSNP_plot_education_ON_LST <- mr_forest_plot(singleSNP_education_ON_LST)
singleSNP_plot_education_ON_LST[[1]]
ggsave(singleSNP_plot_education_ON_LST[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_LST.pdf", width=7, height=21)

PRESSO_education_ON_LST <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_LST, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_LST <- generate_odds_ratios(result_education_ON_LST)
result_education_ON_LST <- result_education_ON_LST[,-c(1,2)]


#### ________M1.2 Physical activity_television ####
str(IV_education)

out_education_ON_television <- format_data(GWAS_television,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "BETA",
                                           se_col = "SE",
                                           effect_allele_col = "effect allele",
                                           other_allele_col = "other allele",
                                           eaf_col = "Freq",
                                           pval_col = "P")
str(out_education_ON_television)
out_education_ON_television$outcome <- "television"

UVdat_education_ON_television <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_television
)
str(UVdat_education_ON_television)
range(UVdat_education_ON_television$F)
range(UVdat_education_ON_television$pval.outcome)
UVdat_education_ON_television <- UVdat_education_ON_television[UVdat_education_ON_television$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_television, file = "Output/UVMR_IV/UVdat_education_ON_television.xlsx")

result_education_ON_television <- mr(UVdat_education_ON_television)


scatter_plot_education_ON_television <- mr_scatter_plot(result_education_ON_television, UVdat_education_ON_television)
scatter_plot_education_ON_television[[1]]
ggsave(scatter_plot_education_ON_television[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_television.pdf", width=7, height=7)


pleiotropy_education_ON_television <- mr_pleiotropy_test(UVdat_education_ON_television)
str(pleiotropy_education_ON_television)
pleiotropy_education_ON_television <- pleiotropy_education_ON_television[,-c(1,2)]


heterogeneity_education_ON_television <- mr_heterogeneity(UVdat_education_ON_television)
str(heterogeneity_education_ON_television)
heterogeneity_education_ON_television <- heterogeneity_education_ON_television[,-c(1,2)]


singleSNP_education_ON_television <- mr_singlesnp(UVdat_education_ON_television)
funnel_plot_education_ON_television <- mr_funnel_plot(singleSNP_education_ON_television)
funnel_plot_education_ON_television[[1]]
ggsave(funnel_plot_education_ON_television[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_television.pdf", width=7, height=7)


loo_education_ON_television <- mr_leaveoneout(UVdat_education_ON_television)
loo_plot_education_ON_television <- mr_leaveoneout_plot(loo_education_ON_television)
loo_plot_education_ON_television[[1]]
ggsave(loo_plot_education_ON_television[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_television.pdf", width=7, height=21)


singleSNP_plot_education_ON_television <- mr_forest_plot(singleSNP_education_ON_television)
singleSNP_plot_education_ON_television[[1]]
ggsave(singleSNP_plot_education_ON_television[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_television.pdf", width=7, height=21)

PRESSO_education_ON_television <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_television, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_television <- generate_odds_ratios(result_education_ON_television)
result_education_ON_television <- result_education_ON_television[,-c(1,2)]


#### ________M1.2 Physical activity_computer ####
str(IV_education)

out_education_ON_computer <- format_data(GWAS_computer,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "BETA",
                                           se_col = "SE",
                                           effect_allele_col = "ALLELE1",
                                           other_allele_col = "ALLELE0",
                                           eaf_col = "A1FREQ",
                                           pval_col = "P_BOLT_LMM_INF")
str(out_education_ON_computer)
out_education_ON_computer$outcome <- "computer"

UVdat_education_ON_computer <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_computer
)
str(UVdat_education_ON_computer)
range(UVdat_education_ON_computer$F)
range(UVdat_education_ON_computer$pval.outcome)
UVdat_education_ON_computer <- UVdat_education_ON_computer[UVdat_education_ON_computer$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_computer, file = "Output/UVMR_IV/UVdat_education_ON_computer.xlsx")

result_education_ON_computer <- mr(UVdat_education_ON_computer)


scatter_plot_education_ON_computer <- mr_scatter_plot(result_education_ON_computer, UVdat_education_ON_computer)
scatter_plot_education_ON_computer[[1]]
ggsave(scatter_plot_education_ON_computer[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_computer.pdf", width=7, height=7)


pleiotropy_education_ON_computer <- mr_pleiotropy_test(UVdat_education_ON_computer)
str(pleiotropy_education_ON_computer)
pleiotropy_education_ON_computer <- pleiotropy_education_ON_computer[,-c(1,2)]


heterogeneity_education_ON_computer <- mr_heterogeneity(UVdat_education_ON_computer)
str(heterogeneity_education_ON_computer)
heterogeneity_education_ON_computer <- heterogeneity_education_ON_computer[,-c(1,2)]


singleSNP_education_ON_computer <- mr_singlesnp(UVdat_education_ON_computer)
funnel_plot_education_ON_computer <- mr_funnel_plot(singleSNP_education_ON_computer)
funnel_plot_education_ON_computer[[1]]
ggsave(funnel_plot_education_ON_computer[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_computer.pdf", width=7, height=7)


loo_education_ON_computer <- mr_leaveoneout(UVdat_education_ON_computer)
loo_plot_education_ON_computer <- mr_leaveoneout_plot(loo_education_ON_computer)
loo_plot_education_ON_computer[[1]]
ggsave(loo_plot_education_ON_computer[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_computer.pdf", width=7, height=21)


singleSNP_plot_education_ON_computer <- mr_forest_plot(singleSNP_education_ON_computer)
singleSNP_plot_education_ON_computer[[1]]
ggsave(singleSNP_plot_education_ON_computer[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_computer.pdf", width=7, height=21)

PRESSO_education_ON_computer <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_computer, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_computer <- generate_odds_ratios(result_education_ON_computer)
result_education_ON_computer <- result_education_ON_computer[,-c(1,2)]


#### ________M1.2 Physical activity_driving ####
str(IV_education)

out_education_ON_driving <- format_data(GWAS_driving,
                                         type = "outcome", 
                                         snps = IV_education$SNP, 
                                         snp_col = "SNP",
                                         beta_col = "BETA",
                                         se_col = "SE",
                                         effect_allele_col = "ALLELE1",
                                         other_allele_col = "ALLELE0",
                                         eaf_col = "A1FREQ",
                                         pval_col = "P_BOLT_LMM_INF")
str(out_education_ON_driving)
out_education_ON_driving$outcome <- "driving"

UVdat_education_ON_driving <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_driving
)
str(UVdat_education_ON_driving)
range(UVdat_education_ON_driving$F)
range(UVdat_education_ON_driving$pval.outcome)
UVdat_education_ON_driving <- UVdat_education_ON_driving[UVdat_education_ON_driving$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_driving, file = "Output/UVMR_IV/UVdat_education_ON_driving.xlsx")

result_education_ON_driving <- mr(UVdat_education_ON_driving)


scatter_plot_education_ON_driving <- mr_scatter_plot(result_education_ON_driving, UVdat_education_ON_driving)
scatter_plot_education_ON_driving[[1]]
ggsave(scatter_plot_education_ON_driving[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_driving.pdf", width=7, height=7)


pleiotropy_education_ON_driving <- mr_pleiotropy_test(UVdat_education_ON_driving)
str(pleiotropy_education_ON_driving)
pleiotropy_education_ON_driving <- pleiotropy_education_ON_driving[,-c(1,2)]


heterogeneity_education_ON_driving <- mr_heterogeneity(UVdat_education_ON_driving)
str(heterogeneity_education_ON_driving)
heterogeneity_education_ON_driving <- heterogeneity_education_ON_driving[,-c(1,2)]


singleSNP_education_ON_driving <- mr_singlesnp(UVdat_education_ON_driving)
funnel_plot_education_ON_driving <- mr_funnel_plot(singleSNP_education_ON_driving)
funnel_plot_education_ON_driving[[1]]
ggsave(funnel_plot_education_ON_driving[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_driving.pdf", width=7, height=7)


loo_education_ON_driving <- mr_leaveoneout(UVdat_education_ON_driving)
loo_plot_education_ON_driving <- mr_leaveoneout_plot(loo_education_ON_driving)
loo_plot_education_ON_driving[[1]]
ggsave(loo_plot_education_ON_driving[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_driving.pdf", width=7, height=21)


singleSNP_plot_education_ON_driving <- mr_forest_plot(singleSNP_education_ON_driving)
singleSNP_plot_education_ON_driving[[1]]
ggsave(singleSNP_plot_education_ON_driving[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_driving.pdf", width=7, height=21)

PRESSO_education_ON_driving <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_driving, NbDistribution = 10000,  
                                          SignifThreshold = 0.05)

result_education_ON_driving <- generate_odds_ratios(result_education_ON_driving)
result_education_ON_driving <- result_education_ON_driving[,-c(1,2)]


#### ________M1.3 Smoking_LifetimeSmoking ####
str(IV_education)

out_education_ON_LifetimeSmoking <- format_data(GWAS_LifetimeSmoking,
                                        type = "outcome", 
                                        snps = IV_education$SNP, 
                                        snp_col = "SNP",
                                        beta_col = "BETA",
                                        se_col = "SE",
                                        effect_allele_col = "EFFECT_ALLELE",
                                        other_allele_col = "OTHER_ALLELE",
                                        eaf_col = "EAF",
                                        pval_col = "P")
str(out_education_ON_LifetimeSmoking)
out_education_ON_LifetimeSmoking$outcome <- "Lifetime smoking index"

UVdat_education_ON_LifetimeSmoking <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_LifetimeSmoking
)
str(UVdat_education_ON_LifetimeSmoking)
range(UVdat_education_ON_LifetimeSmoking$F)
range(UVdat_education_ON_LifetimeSmoking$pval.outcome)
UVdat_education_ON_LifetimeSmoking <- UVdat_education_ON_LifetimeSmoking[UVdat_education_ON_LifetimeSmoking$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_LifetimeSmoking, file = "Output/UVMR_IV/UVdat_education_ON_LifetimeSmoking.xlsx")

result_education_ON_LifetimeSmoking <- mr(UVdat_education_ON_LifetimeSmoking)


scatter_plot_education_ON_LifetimeSmoking <- mr_scatter_plot(result_education_ON_LifetimeSmoking, UVdat_education_ON_LifetimeSmoking)
scatter_plot_education_ON_LifetimeSmoking[[1]]
ggsave(scatter_plot_education_ON_LifetimeSmoking[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_LifetimeSmoking.pdf", width=7, height=7)


pleiotropy_education_ON_LifetimeSmoking <- mr_pleiotropy_test(UVdat_education_ON_LifetimeSmoking)
str(pleiotropy_education_ON_LifetimeSmoking)
pleiotropy_education_ON_LifetimeSmoking <- pleiotropy_education_ON_LifetimeSmoking[,-c(1,2)]


heterogeneity_education_ON_LifetimeSmoking <- mr_heterogeneity(UVdat_education_ON_LifetimeSmoking)
str(heterogeneity_education_ON_LifetimeSmoking)
heterogeneity_education_ON_LifetimeSmoking <- heterogeneity_education_ON_LifetimeSmoking[,-c(1,2)]


singleSNP_education_ON_LifetimeSmoking <- mr_singlesnp(UVdat_education_ON_LifetimeSmoking)
funnel_plot_education_ON_LifetimeSmoking <- mr_funnel_plot(singleSNP_education_ON_LifetimeSmoking)
funnel_plot_education_ON_LifetimeSmoking[[1]]
ggsave(funnel_plot_education_ON_LifetimeSmoking[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_LifetimeSmoking.pdf", width=7, height=7)


loo_education_ON_LifetimeSmoking <- mr_leaveoneout(UVdat_education_ON_LifetimeSmoking)
loo_plot_education_ON_LifetimeSmoking <- mr_leaveoneout_plot(loo_education_ON_LifetimeSmoking)
loo_plot_education_ON_LifetimeSmoking[[1]]
ggsave(loo_plot_education_ON_LifetimeSmoking[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_LifetimeSmoking.pdf", width=7, height=21)


singleSNP_plot_education_ON_LifetimeSmoking <- mr_forest_plot(singleSNP_education_ON_LifetimeSmoking)
singleSNP_plot_education_ON_LifetimeSmoking[[1]]
ggsave(singleSNP_plot_education_ON_LifetimeSmoking[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_LifetimeSmoking.pdf", width=7, height=21)

PRESSO_education_ON_LifetimeSmoking <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_LifetimeSmoking, NbDistribution = 10000,  
                                         SignifThreshold = 0.05)

result_education_ON_LifetimeSmoking <- generate_odds_ratios(result_education_ON_LifetimeSmoking)
str(result_education_ON_LifetimeSmoking)
result_education_ON_LifetimeSmoking <- result_education_ON_LifetimeSmoking[,-c(1,2)]


#### ________M1.4 Sleep health_SleepDuration ####
str(IV_education)

out_education_ON_SleepDuration <- format_data(GWAS_SleepDuration,
                                      type = "outcome", 
                                      snps = IV_education$SNP, 
                                      snp_col = "SNP",
                                      beta_col = "BETA",
                                      se_col = "SE",
                                      effect_allele_col = "A1",
                                      other_allele_col = "A2",
                                      eaf_col = "MAF",
                                      pval_col = "P")
str(out_education_ON_SleepDuration)
out_education_ON_SleepDuration$outcome <- "SleepDuration"

UVdat_education_ON_SleepDuration <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_SleepDuration
)
str(UVdat_education_ON_SleepDuration)
range(UVdat_education_ON_SleepDuration$F)
range(UVdat_education_ON_SleepDuration$pval.outcome)
UVdat_education_ON_SleepDuration <- UVdat_education_ON_SleepDuration[UVdat_education_ON_SleepDuration$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_SleepDuration, file = "Output/UVMR_IV/UVdat_education_ON_SleepDuration.xlsx")

result_education_ON_SleepDuration <- mr(UVdat_education_ON_SleepDuration)


scatter_plot_education_ON_SleepDuration <- mr_scatter_plot(result_education_ON_SleepDuration, UVdat_education_ON_SleepDuration)
scatter_plot_education_ON_SleepDuration[[1]]
ggsave(scatter_plot_education_ON_SleepDuration[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_SleepDuration.pdf", width=7, height=7)


pleiotropy_education_ON_SleepDuration <- mr_pleiotropy_test(UVdat_education_ON_SleepDuration)
str(pleiotropy_education_ON_SleepDuration)
pleiotropy_education_ON_SleepDuration <- pleiotropy_education_ON_SleepDuration[,-c(1,2)]


heterogeneity_education_ON_SleepDuration <- mr_heterogeneity(UVdat_education_ON_SleepDuration)
str(heterogeneity_education_ON_SleepDuration)
heterogeneity_education_ON_SleepDuration <- heterogeneity_education_ON_SleepDuration[,-c(1,2)]


singleSNP_education_ON_SleepDuration <- mr_singlesnp(UVdat_education_ON_SleepDuration)
funnel_plot_education_ON_SleepDuration <- mr_funnel_plot(singleSNP_education_ON_SleepDuration)
funnel_plot_education_ON_SleepDuration[[1]]
ggsave(funnel_plot_education_ON_SleepDuration[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_SleepDuration.pdf", width=7, height=7)


loo_education_ON_SleepDuration <- mr_leaveoneout(UVdat_education_ON_SleepDuration)
loo_plot_education_ON_SleepDuration <- mr_leaveoneout_plot(loo_education_ON_SleepDuration)
loo_plot_education_ON_SleepDuration[[1]]
ggsave(loo_plot_education_ON_SleepDuration[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_SleepDuration.pdf", width=7, height=21)


singleSNP_plot_education_ON_SleepDuration <- mr_forest_plot(singleSNP_education_ON_SleepDuration)
singleSNP_plot_education_ON_SleepDuration[[1]]
ggsave(singleSNP_plot_education_ON_SleepDuration[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_SleepDuration.pdf", width=7, height=21)

PRESSO_education_ON_SleepDuration <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_SleepDuration, NbDistribution = 10000,  
                                       SignifThreshold = 0.05)

result_education_ON_SleepDuration <- generate_odds_ratios(result_education_ON_SleepDuration)
str(result_education_ON_SleepDuration)
result_education_ON_SleepDuration <- result_education_ON_SleepDuration[,-c(1,2)]


#### ________M1.4 Sleep health_DeviceSleepDuration ####
str(IV_education)

out_education_ON_DeviceSleepDuration <- format_data(GWAS_DeviceSleepDuration,
                                                type = "outcome", 
                                                snps = IV_education$SNP, 
                                                snp_col = "SNP",
                                                beta_col = "BETA",
                                                se_col = "SE",
                                                effect_allele_col = "ALLELE1",
                                                other_allele_col = "ALLELE0",
                                                eaf_col = "A1FREQ",
                                                pval_col = "P_BOLT_LMM_INF")
str(out_education_ON_DeviceSleepDuration)
out_education_ON_DeviceSleepDuration$outcome <- "DeviceSleepDuration"

UVdat_education_ON_DeviceSleepDuration <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_DeviceSleepDuration
)
str(UVdat_education_ON_DeviceSleepDuration)
range(UVdat_education_ON_DeviceSleepDuration$pval.outcome)
range(UVdat_education_ON_DeviceSleepDuration$F)
UVdat_education_ON_DeviceSleepDuration <- UVdat_education_ON_DeviceSleepDuration[UVdat_education_ON_DeviceSleepDuration$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_DeviceSleepDuration, file = "Output/UVMR_IV/UVdat_education_ON_DeviceSleepDuration.xlsx")

result_education_ON_DeviceSleepDuration <- mr(UVdat_education_ON_DeviceSleepDuration)


scatter_plot_education_ON_DeviceSleepDuration <- mr_scatter_plot(result_education_ON_DeviceSleepDuration, UVdat_education_ON_DeviceSleepDuration)
scatter_plot_education_ON_DeviceSleepDuration[[1]]
ggsave(scatter_plot_education_ON_DeviceSleepDuration[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_DeviceSleepDuration.pdf", width=7, height=7)


pleiotropy_education_ON_DeviceSleepDuration <- mr_pleiotropy_test(UVdat_education_ON_DeviceSleepDuration)
str(pleiotropy_education_ON_DeviceSleepDuration)
pleiotropy_education_ON_DeviceSleepDuration <- pleiotropy_education_ON_DeviceSleepDuration[,-c(1,2)]


heterogeneity_education_ON_DeviceSleepDuration <- mr_heterogeneity(UVdat_education_ON_DeviceSleepDuration)
str(heterogeneity_education_ON_DeviceSleepDuration)
heterogeneity_education_ON_DeviceSleepDuration <- heterogeneity_education_ON_DeviceSleepDuration[,-c(1,2)]


singleSNP_education_ON_DeviceSleepDuration <- mr_singlesnp(UVdat_education_ON_DeviceSleepDuration)
funnel_plot_education_ON_DeviceSleepDuration <- mr_funnel_plot(singleSNP_education_ON_DeviceSleepDuration)
funnel_plot_education_ON_DeviceSleepDuration[[1]]
ggsave(funnel_plot_education_ON_DeviceSleepDuration[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_DeviceSleepDuration.pdf", width=7, height=7)


loo_education_ON_DeviceSleepDuration <- mr_leaveoneout(UVdat_education_ON_DeviceSleepDuration)
loo_plot_education_ON_DeviceSleepDuration <- mr_leaveoneout_plot(loo_education_ON_DeviceSleepDuration)
loo_plot_education_ON_DeviceSleepDuration[[1]]
ggsave(loo_plot_education_ON_DeviceSleepDuration[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_DeviceSleepDuration.pdf", width=7, height=21)


singleSNP_plot_education_ON_DeviceSleepDuration <- mr_forest_plot(singleSNP_education_ON_DeviceSleepDuration)
singleSNP_plot_education_ON_DeviceSleepDuration[[1]]
ggsave(singleSNP_plot_education_ON_DeviceSleepDuration[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_DeviceSleepDuration.pdf", width=7, height=21)

PRESSO_education_ON_DeviceSleepDuration <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                                 OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_DeviceSleepDuration, NbDistribution = 10000,  
                                                 SignifThreshold = 0.05)

result_education_ON_DeviceSleepDuration <- generate_odds_ratios(result_education_ON_DeviceSleepDuration)
result_education_ON_DeviceSleepDuration <- result_education_ON_DeviceSleepDuration[,-c(1,2)]


#### ________M1.4 Sleep health_Insomnia ####
str(IV_education)

out_education_ON_Insomnia <- format_data(GWAS_Insomnia,
                                         type = "outcome", 
                                         snps = IV_education$SNP, 
                                         snp_col = "SNP",
                                         beta_col = "BETA",
                                         se_col = "SE",
                                         effect_allele_col = "A1",
                                         other_allele_col = "A2",
                                         eaf_col = "MAF",
                                         pval_col = "P")
str(out_education_ON_Insomnia)
out_education_ON_Insomnia$outcome <- "Insomnia"

UVdat_education_ON_Insomnia <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_Insomnia
)
str(UVdat_education_ON_Insomnia)
range(UVdat_education_ON_Insomnia$F)
range(UVdat_education_ON_Insomnia$pval.outcome)
UVdat_education_ON_Insomnia <- UVdat_education_ON_Insomnia[UVdat_education_ON_Insomnia$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_Insomnia, file = "Output/UVMR_IV/UVdat_education_ON_Insomnia.xlsx")

result_education_ON_Insomnia <- mr(UVdat_education_ON_Insomnia)


scatter_plot_education_ON_Insomnia <- mr_scatter_plot(result_education_ON_Insomnia, UVdat_education_ON_Insomnia)
scatter_plot_education_ON_Insomnia[[1]]
ggsave(scatter_plot_education_ON_Insomnia[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_Insomnia.pdf", width=7, height=7)


pleiotropy_education_ON_Insomnia <- mr_pleiotropy_test(UVdat_education_ON_Insomnia)
str(pleiotropy_education_ON_Insomnia)
pleiotropy_education_ON_Insomnia <- pleiotropy_education_ON_Insomnia[,-c(1,2)]


heterogeneity_education_ON_Insomnia <- mr_heterogeneity(UVdat_education_ON_Insomnia)
str(heterogeneity_education_ON_Insomnia)
heterogeneity_education_ON_Insomnia <- heterogeneity_education_ON_Insomnia[,-c(1,2)]


singleSNP_education_ON_Insomnia <- mr_singlesnp(UVdat_education_ON_Insomnia)
funnel_plot_education_ON_Insomnia <- mr_funnel_plot(singleSNP_education_ON_Insomnia)
funnel_plot_education_ON_Insomnia[[1]]
ggsave(funnel_plot_education_ON_Insomnia[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_Insomnia.pdf", width=7, height=7)


loo_education_ON_Insomnia <- mr_leaveoneout(UVdat_education_ON_Insomnia)
loo_plot_education_ON_Insomnia <- mr_leaveoneout_plot(loo_education_ON_Insomnia)
loo_plot_education_ON_Insomnia[[1]]
ggsave(loo_plot_education_ON_Insomnia[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_Insomnia.pdf", width=7, height=21)


singleSNP_plot_education_ON_Insomnia <- mr_forest_plot(singleSNP_education_ON_Insomnia)
singleSNP_plot_education_ON_Insomnia[[1]]
ggsave(singleSNP_plot_education_ON_Insomnia[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_Insomnia.pdf", width=7, height=21)

PRESSO_education_ON_Insomnia <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_Insomnia, NbDistribution = 10000,  
                                          SignifThreshold = 0.05)

result_education_ON_Insomnia <- generate_odds_ratios(result_education_ON_Insomnia)
result_education_ON_Insomnia <- result_education_ON_Insomnia[,-c(1,2)]


#### ________M2.1 Well-being spectrum(WBS) ####
str(IV_education)

out_education_ON_WBS <- format_data(GWAS_WBS,
                                    type = "outcome", 
                                    snps = IV_education$SNP, 
                                    snp_col = "SNP",
                                    beta_col = "BETA",
                                    se_col = "SE",
                                    effect_allele_col = "A1",
                                    other_allele_col = "A2",
                                    eaf_col = "EAF",
                                    pval_col = "PVAL")
str(out_education_ON_WBS)
out_education_ON_WBS$outcome <- "WBS"

UVdat_education_ON_WBS <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_WBS
)
str(UVdat_education_ON_WBS)
range(UVdat_education_ON_WBS$F) 
range(UVdat_education_ON_WBS$pval.outcome) 
UVdat_education_ON_WBS <- UVdat_education_ON_WBS[UVdat_education_ON_WBS$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_WBS, file = "Output/UVMR_IV/UVdat_education_ON_WBS.xlsx")

result_education_ON_WBS <- mr(UVdat_education_ON_WBS)


scatter_plot_education_ON_WBS <- mr_scatter_plot(result_education_ON_WBS, UVdat_education_ON_WBS)
scatter_plot_education_ON_WBS[[1]]
ggsave(scatter_plot_education_ON_WBS[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_WBS.pdf", width=7, height=7)


pleiotropy_education_ON_WBS <- mr_pleiotropy_test(UVdat_education_ON_WBS)
str(pleiotropy_education_ON_WBS)
pleiotropy_education_ON_WBS <- pleiotropy_education_ON_WBS[,-c(1,2)]


heterogeneity_education_ON_WBS <- mr_heterogeneity(UVdat_education_ON_WBS)
str(heterogeneity_education_ON_WBS)
heterogeneity_education_ON_WBS <- heterogeneity_education_ON_WBS[,-c(1,2)]


singleSNP_education_ON_WBS <- mr_singlesnp(UVdat_education_ON_WBS)
funnel_plot_education_ON_WBS <- mr_funnel_plot(singleSNP_education_ON_WBS)
funnel_plot_education_ON_WBS[[1]]
ggsave(funnel_plot_education_ON_WBS[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_WBS.pdf", width=7, height=7)


loo_education_ON_WBS <- mr_leaveoneout(UVdat_education_ON_WBS)
loo_plot_education_ON_WBS <- mr_leaveoneout_plot(loo_education_ON_WBS)
loo_plot_education_ON_WBS[[1]]
ggsave(loo_plot_education_ON_WBS[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_WBS.pdf", width=7, height=21)


singleSNP_plot_education_ON_WBS <- mr_forest_plot(singleSNP_education_ON_WBS)
singleSNP_plot_education_ON_WBS[[1]]
ggsave(singleSNP_plot_education_ON_WBS[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_WBS.pdf", width=7, height=21)

PRESSO_education_ON_WBS <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_WBS, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_WBS <- generate_odds_ratios(result_education_ON_WBS)
result_education_ON_WBS <- result_education_ON_WBS[,-c(1,2)]


#### ________M2.2 MentalProblems ####
str(IV_education)

out_education_ON_MentalProblems <- format_data(GWAS_MentalProblems,
                                    type = "outcome", 
                                    snps = IV_education$SNP, 
                                    snp_col = "SNP",
                                    beta_col = "ES",
                                    se_col = "SE",
                                    effect_allele_col = "ALT",
                                    other_allele_col = "REF",
                                    eaf_col = "AF",
                                    pval_col = "P")
str(out_education_ON_MentalProblems)
out_education_ON_MentalProblems$outcome <- "MentalProblems"

UVdat_education_ON_MentalProblems <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_MentalProblems
)
str(UVdat_education_ON_MentalProblems)
range(UVdat_education_ON_MentalProblems$F) 
range(UVdat_education_ON_MentalProblems$pval.outcome) 
UVdat_education_ON_MentalProblems <- UVdat_education_ON_MentalProblems[UVdat_education_ON_MentalProblems$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_MentalProblems, file = "Output/UVMR_IV/UVdat_education_ON_MentalProblems.xlsx")

result_education_ON_MentalProblems <- mr(UVdat_education_ON_MentalProblems)


scatter_plot_education_ON_MentalProblems <- mr_scatter_plot(result_education_ON_MentalProblems, UVdat_education_ON_MentalProblems)
scatter_plot_education_ON_MentalProblems[[1]]
ggsave(scatter_plot_education_ON_MentalProblems[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_MentalProblems.pdf", width=7, height=7)


pleiotropy_education_ON_MentalProblems <- mr_pleiotropy_test(UVdat_education_ON_MentalProblems)
str(pleiotropy_education_ON_MentalProblems)
pleiotropy_education_ON_MentalProblems <- pleiotropy_education_ON_MentalProblems[,-c(1,2)]


heterogeneity_education_ON_MentalProblems <- mr_heterogeneity(UVdat_education_ON_MentalProblems)
str(heterogeneity_education_ON_MentalProblems)
heterogeneity_education_ON_MentalProblems <- heterogeneity_education_ON_MentalProblems[,-c(1,2)]


singleSNP_education_ON_MentalProblems <- mr_singlesnp(UVdat_education_ON_MentalProblems)
funnel_plot_education_ON_MentalProblems <- mr_funnel_plot(singleSNP_education_ON_MentalProblems)
funnel_plot_education_ON_MentalProblems[[1]]
ggsave(funnel_plot_education_ON_MentalProblems[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_MentalProblems.pdf", width=7, height=7)


loo_education_ON_MentalProblems <- mr_leaveoneout(UVdat_education_ON_MentalProblems)
loo_plot_education_ON_MentalProblems <- mr_leaveoneout_plot(loo_education_ON_MentalProblems)
loo_plot_education_ON_MentalProblems[[1]]
ggsave(loo_plot_education_ON_MentalProblems[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_MentalProblems.pdf", width=7, height=21)


singleSNP_plot_education_ON_MentalProblems <- mr_forest_plot(singleSNP_education_ON_MentalProblems)
singleSNP_plot_education_ON_MentalProblems[[1]]
ggsave(singleSNP_plot_education_ON_MentalProblems[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_MentalProblems.pdf", width=7, height=21)

PRESSO_education_ON_MentalProblems <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_MentalProblems, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_MentalProblems <- generate_odds_ratios(result_education_ON_MentalProblems)
result_education_ON_MentalProblems <- result_education_ON_MentalProblems[,-c(1,2)]


#### ________M2.3 Depression ####
str(IV_education)

out_education_ON_depression <- format_data(GWAS_depression,
                                           type = "outcome", 
                                           snps = IV_education$SNP, 
                                           snp_col = "SNP",
                                           beta_col = "LogOR",
                                           se_col = "StdErrLogOR",
                                           effect_allele_col = "A1",
                                           other_allele_col = "A2",
                                           eaf_col = "Freq",
                                           pval_col = "P")
str(out_education_ON_depression)
out_education_ON_depression$outcome <- "depression"

UVdat_education_ON_depression <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_depression
)
str(UVdat_education_ON_depression)
range(UVdat_education_ON_depression$F)
range(UVdat_education_ON_depression$pval.outcome)
UVdat_education_ON_depression <- UVdat_education_ON_depression[UVdat_education_ON_depression$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_depression, file = "Output/UVMR_IV/UVdat_education_ON_depression.xlsx")

result_education_ON_depression <- mr(UVdat_education_ON_depression)


scatter_plot_education_ON_depression <- mr_scatter_plot(result_education_ON_depression, UVdat_education_ON_depression)
scatter_plot_education_ON_depression[[1]]
ggsave(scatter_plot_education_ON_depression[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_depression.pdf", width=7, height=7)


pleiotropy_education_ON_depression <- mr_pleiotropy_test(UVdat_education_ON_depression)
str(pleiotropy_education_ON_depression)
pleiotropy_education_ON_depression <- pleiotropy_education_ON_depression[,-c(1,2)]


heterogeneity_education_ON_depression <- mr_heterogeneity(UVdat_education_ON_depression)
str(heterogeneity_education_ON_depression)
heterogeneity_education_ON_depression <- heterogeneity_education_ON_depression[,-c(1,2)]


singleSNP_education_ON_depression <- mr_singlesnp(UVdat_education_ON_depression)
funnel_plot_education_ON_depression <- mr_funnel_plot(singleSNP_education_ON_depression)
funnel_plot_education_ON_depression[[1]]
ggsave(funnel_plot_education_ON_depression[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_depression.pdf", width=7, height=7)


loo_education_ON_depression <- mr_leaveoneout(UVdat_education_ON_depression)
loo_plot_education_ON_depression <- mr_leaveoneout_plot(loo_education_ON_depression)
loo_plot_education_ON_depression[[1]]
ggsave(loo_plot_education_ON_depression[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_depression.pdf", width=7, height=21)


singleSNP_plot_education_ON_depression <- mr_forest_plot(singleSNP_education_ON_depression)
singleSNP_plot_education_ON_depression[[1]]
ggsave(singleSNP_plot_education_ON_depression[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_depression.pdf", width=7, height=21)

PRESSO_education_ON_depression <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_depression, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_education_ON_depression <- generate_odds_ratios(result_education_ON_depression)
str(result_education_ON_depression)
result_education_ON_depression <- result_education_ON_depression[,-c(1,2)]
result_education_ON_depression$outcome <- "depression"


#### ________M3.1 BMI ####
str(IV_education)

out_education_ON_BMI <- extract_outcome_data(
  snps = IV_education$SNP,
  outcomes = 'ieu-b-40')
str(out_education_ON_BMI)
out_education_ON_BMI$outcome <- "BMI"

UVdat_education_ON_BMI <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_BMI
)
str(UVdat_education_ON_BMI)
range(UVdat_education_ON_BMI$pval.outcome)
range(UVdat_education_ON_BMI$F)
UVdat_education_ON_BMI <- UVdat_education_ON_BMI[UVdat_education_ON_BMI$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_BMI, file = "Output/UVMR_IV/UVdat_education_ON_BMI.xlsx")

result_education_ON_BMI <- mr(UVdat_education_ON_BMI)


scatter_plot_education_ON_BMI <- mr_scatter_plot(result_education_ON_BMI, UVdat_education_ON_BMI)
scatter_plot_education_ON_BMI[[1]]
ggsave(scatter_plot_education_ON_BMI[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_BMI.pdf", width=7, height=7)


pleiotropy_education_ON_BMI <- mr_pleiotropy_test(UVdat_education_ON_BMI)
str(pleiotropy_education_ON_BMI)
pleiotropy_education_ON_BMI <- pleiotropy_education_ON_BMI[,-c(1,2)]


heterogeneity_education_ON_BMI <- mr_heterogeneity(UVdat_education_ON_BMI)
str(heterogeneity_education_ON_BMI)
heterogeneity_education_ON_BMI <- heterogeneity_education_ON_BMI[,-c(1,2)]


singleSNP_education_ON_BMI <- mr_singlesnp(UVdat_education_ON_BMI)
funnel_plot_education_ON_BMI <- mr_funnel_plot(singleSNP_education_ON_BMI)
funnel_plot_education_ON_BMI[[1]]
ggsave(funnel_plot_education_ON_BMI[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_BMI.pdf", width=7, height=7)


loo_education_ON_BMI <- mr_leaveoneout(UVdat_education_ON_BMI)
loo_plot_education_ON_BMI <- mr_leaveoneout_plot(loo_education_ON_BMI)
loo_plot_education_ON_BMI[[1]]
ggsave(loo_plot_education_ON_BMI[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_BMI.pdf", width=7, height=21)


singleSNP_plot_education_ON_BMI <- mr_forest_plot(singleSNP_education_ON_BMI)
singleSNP_plot_education_ON_BMI[[1]]
ggsave(singleSNP_plot_education_ON_BMI[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_BMI.pdf", width=7, height=21)

PRESSO_education_ON_BMI <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_BMI, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_BMI <- generate_odds_ratios(result_education_ON_BMI)
result_education_ON_BMI <- result_education_ON_BMI[,-c(1,2)]


#### ________M3.2 lipid ####
str(IV_education)

out_education_ON_lipid <- format_data(GWAS_lipid,
                                      type = "outcome", 
                                      snps = IV_education$SNP, 
                                      snp_col = "SNP",
                                      beta_col = "BETA",
                                      se_col = "SE",
                                      effect_allele_col = "effect allele",
                                      other_allele_col = "other allele",
                                      eaf_col = "Freq",
                                      pval_col = "P")
str(out_education_ON_lipid)
out_education_ON_lipid$outcome <- "blood lipids"

UVdat_education_ON_lipid <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_lipid
)
range(UVdat_education_ON_lipid$pval.outcome)
range(UVdat_education_ON_lipid$F)
str(UVdat_education_ON_lipid)
UVdat_education_ON_lipid <- UVdat_education_ON_lipid[UVdat_education_ON_lipid$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_lipid, file = "Output/UVMR_IV/UVdat_education_ON_lipid.xlsx")

result_education_ON_lipid <- mr(UVdat_education_ON_lipid)


scatter_plot_education_ON_lipid <- mr_scatter_plot(result_education_ON_lipid, UVdat_education_ON_lipid)
scatter_plot_education_ON_lipid[[1]]
ggsave(scatter_plot_education_ON_lipid[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_lipid.pdf", width=7, height=7)


pleiotropy_education_ON_lipid <- mr_pleiotropy_test(UVdat_education_ON_lipid)
str(pleiotropy_education_ON_lipid)
pleiotropy_education_ON_lipid <- pleiotropy_education_ON_lipid[,-c(1,2)]


heterogeneity_education_ON_lipid <- mr_heterogeneity(UVdat_education_ON_lipid)
str(heterogeneity_education_ON_lipid)
heterogeneity_education_ON_lipid <- heterogeneity_education_ON_lipid[,-c(1,2)]


singleSNP_education_ON_lipid <- mr_singlesnp(UVdat_education_ON_lipid)
funnel_plot_education_ON_lipid <- mr_funnel_plot(singleSNP_education_ON_lipid)
funnel_plot_education_ON_lipid[[1]]
ggsave(funnel_plot_education_ON_lipid[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_lipid.pdf", width=7, height=7)


loo_education_ON_lipid <- mr_leaveoneout(UVdat_education_ON_lipid)
loo_plot_education_ON_lipid <- mr_leaveoneout_plot(loo_education_ON_lipid)
loo_plot_education_ON_lipid[[1]]
ggsave(loo_plot_education_ON_lipid[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_lipid.pdf", width=7, height=21)


singleSNP_plot_education_ON_lipid <- mr_forest_plot(singleSNP_education_ON_lipid)
singleSNP_plot_education_ON_lipid[[1]]
ggsave(singleSNP_plot_education_ON_lipid[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_lipid.pdf", width=7, height=21)

PRESSO_education_ON_lipid <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_lipid, NbDistribution = 10000,  
                                       SignifThreshold = 0.05)

result_education_ON_lipid <- generate_odds_ratios(result_education_ON_lipid)
result_education_ON_lipid <- result_education_ON_lipid[,-c(1,2)]


#### ________M3.3 glucose ####
str(IV_education)

out_education_ON_glucose <- format_data(GWAS_glucose,
                                        type = "outcome", 
                                        snps = IV_education$SNP, 
                                        snp_col = "SNP",
                                        beta_col = "BETA",
                                        se_col = "SE",
                                        effect_allele_col = "effect allele",
                                        other_allele_col = "other allele",
                                        eaf_col = "Freq",
                                        pval_col = "P")
str(out_education_ON_glucose)
out_education_ON_glucose$outcome <- "blood glucose"

UVdat_education_ON_glucose <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_glucose
)
range(UVdat_education_ON_glucose$pval.outcome)
range(UVdat_education_ON_glucose$F)
str(UVdat_education_ON_glucose)
UVdat_education_ON_glucose <- UVdat_education_ON_glucose[UVdat_education_ON_glucose$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_glucose, file = "Output/UVMR_IV/UVdat_education_ON_glucose.xlsx")

result_education_ON_glucose <- mr(UVdat_education_ON_glucose)


scatter_plot_education_ON_glucose <- mr_scatter_plot(result_education_ON_glucose, UVdat_education_ON_glucose)
scatter_plot_education_ON_glucose[[1]]
ggsave(scatter_plot_education_ON_glucose[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_glucose.pdf", width=7, height=7)


pleiotropy_education_ON_glucose <- mr_pleiotropy_test(UVdat_education_ON_glucose)
str(pleiotropy_education_ON_glucose)
pleiotropy_education_ON_glucose <- pleiotropy_education_ON_glucose[,-c(1,2)]


heterogeneity_education_ON_glucose <- mr_heterogeneity(UVdat_education_ON_glucose)
str(heterogeneity_education_ON_glucose)
heterogeneity_education_ON_glucose <- heterogeneity_education_ON_glucose[,-c(1,2)]


singleSNP_education_ON_glucose <- mr_singlesnp(UVdat_education_ON_glucose)
funnel_plot_education_ON_glucose <- mr_funnel_plot(singleSNP_education_ON_glucose)
funnel_plot_education_ON_glucose[[1]]
ggsave(funnel_plot_education_ON_glucose[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_glucose.pdf", width=7, height=7)


loo_education_ON_glucose <- mr_leaveoneout(UVdat_education_ON_glucose)
loo_plot_education_ON_glucose <- mr_leaveoneout_plot(loo_education_ON_glucose)
loo_plot_education_ON_glucose[[1]]
ggsave(loo_plot_education_ON_glucose[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_glucose.pdf", width=7, height=21)


singleSNP_plot_education_ON_glucose <- mr_forest_plot(singleSNP_education_ON_glucose)
singleSNP_plot_education_ON_glucose[[1]]
ggsave(singleSNP_plot_education_ON_glucose[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_glucose.pdf", width=7, height=21)

PRESSO_education_ON_glucose <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_glucose, NbDistribution = 10000,  
                                         SignifThreshold = 0.05)

result_education_ON_glucose <- generate_odds_ratios(result_education_ON_glucose)
result_education_ON_glucose <- result_education_ON_glucose[,-c(1,2)]


#### ________M3.4 SBP, DBP ####

## SBP
str(IV_education)

out_education_ON_SBP <- extract_outcome_data(
  snps = IV_education$SNP,
  outcomes = 'ieu-b-38')

str(out_education_ON_SBP)
out_education_ON_SBP$outcome <- "SBP"

UVdat_education_ON_SBP <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_SBP
)
range(UVdat_education_ON_SBP$pval.outcome)
range(UVdat_education_ON_SBP$F)
str(UVdat_education_ON_SBP)
UVdat_education_ON_SBP <- UVdat_education_ON_SBP[UVdat_education_ON_SBP$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_SBP, file = "Output/UVMR_IV/UVdat_education_ON_SBP.xlsx")

result_education_ON_SBP <- mr(UVdat_education_ON_SBP)


scatter_plot_education_ON_SBP <- mr_scatter_plot(result_education_ON_SBP, UVdat_education_ON_SBP)
scatter_plot_education_ON_SBP[[1]]
ggsave(scatter_plot_education_ON_SBP[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_SBP.pdf", width=7, height=7)


pleiotropy_education_ON_SBP <- mr_pleiotropy_test(UVdat_education_ON_SBP)
str(pleiotropy_education_ON_SBP)
pleiotropy_education_ON_SBP <- pleiotropy_education_ON_SBP[,-c(1,2)]


heterogeneity_education_ON_SBP <- mr_heterogeneity(UVdat_education_ON_SBP)
str(heterogeneity_education_ON_SBP)
heterogeneity_education_ON_SBP <- heterogeneity_education_ON_SBP[,-c(1,2)]


singleSNP_education_ON_SBP <- mr_singlesnp(UVdat_education_ON_SBP)
funnel_plot_education_ON_SBP <- mr_funnel_plot(singleSNP_education_ON_SBP)
funnel_plot_education_ON_SBP[[1]]
ggsave(funnel_plot_education_ON_SBP[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_SBP.pdf", width=7, height=7)


loo_education_ON_SBP <- mr_leaveoneout(UVdat_education_ON_SBP)
loo_plot_education_ON_SBP <- mr_leaveoneout_plot(loo_education_ON_SBP)
loo_plot_education_ON_SBP[[1]]
ggsave(loo_plot_education_ON_SBP[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_SBP.pdf", width=7, height=21)


singleSNP_plot_education_ON_SBP <- mr_forest_plot(singleSNP_education_ON_SBP)
singleSNP_plot_education_ON_SBP[[1]]
ggsave(singleSNP_plot_education_ON_SBP[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_SBP.pdf", width=7, height=21)

PRESSO_education_ON_SBP <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_SBP, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_SBP <- generate_odds_ratios(result_education_ON_SBP)
result_education_ON_SBP <- result_education_ON_SBP[,-c(1,2)]


## DBP
str(IV_education)

out_education_ON_DBP <- extract_outcome_data(
  snps = IV_education$SNP,
  outcomes = 'ieu-b-39')

str(out_education_ON_DBP)
out_education_ON_DBP$outcome <- "DBP"

UVdat_education_ON_DBP <- harmonise_data(
  exposure_dat =  IV_education, 
  outcome_dat = out_education_ON_DBP
)
range(UVdat_education_ON_DBP$pval.outcome)
range(UVdat_education_ON_DBP$F)
str(UVdat_education_ON_DBP)
UVdat_education_ON_DBP <- UVdat_education_ON_DBP[UVdat_education_ON_DBP$pval.outcome>5e-8,]
write.xlsx(UVdat_education_ON_DBP, file = "Output/UVMR_IV/UVdat_education_ON_DBP.xlsx")

result_education_ON_DBP <- mr(UVdat_education_ON_DBP)


scatter_plot_education_ON_DBP <- mr_scatter_plot(result_education_ON_DBP, UVdat_education_ON_DBP)
scatter_plot_education_ON_DBP[[1]]
ggsave(scatter_plot_education_ON_DBP[[1]], file="Output/UVMR_Secondary Results/scatter_plot_education_ON_DBP.pdf", width=7, height=7)


pleiotropy_education_ON_DBP <- mr_pleiotropy_test(UVdat_education_ON_DBP)
str(pleiotropy_education_ON_DBP)
pleiotropy_education_ON_DBP <- pleiotropy_education_ON_DBP[,-c(1,2)]


heterogeneity_education_ON_DBP <- mr_heterogeneity(UVdat_education_ON_DBP)
str(heterogeneity_education_ON_DBP)
heterogeneity_education_ON_DBP <- heterogeneity_education_ON_DBP[,-c(1,2)]


singleSNP_education_ON_DBP <- mr_singlesnp(UVdat_education_ON_DBP)
funnel_plot_education_ON_DBP <- mr_funnel_plot(singleSNP_education_ON_DBP)
funnel_plot_education_ON_DBP[[1]]
ggsave(funnel_plot_education_ON_DBP[[1]], file="Output/UVMR_Secondary Results/funnel_plot_education_ON_DBP.pdf", width=7, height=7)


loo_education_ON_DBP <- mr_leaveoneout(UVdat_education_ON_DBP)
loo_plot_education_ON_DBP <- mr_leaveoneout_plot(loo_education_ON_DBP)
loo_plot_education_ON_DBP[[1]]
ggsave(loo_plot_education_ON_DBP[[1]], file="Output/UVMR_Secondary Results/loo_plot_education_ON_DBP.pdf", width=7, height=21)


singleSNP_plot_education_ON_DBP <- mr_forest_plot(singleSNP_education_ON_DBP)
singleSNP_plot_education_ON_DBP[[1]]
ggsave(singleSNP_plot_education_ON_DBP[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_education_ON_DBP.pdf", width=7, height=21)

PRESSO_education_ON_DBP <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_education_ON_DBP, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_education_ON_DBP <- generate_odds_ratios(result_education_ON_DBP)
result_education_ON_DBP <- result_education_ON_DBP[,-c(1,2)]



#### ________M1.1 Diet_protein ####
IV_protein <-format_data(GWAS_protein,
                            type = "exposure", 
                            snp_col = "SNP", 
                            beta_col = "BETA",
                            se_col = "SE",
                            effect_allele_col = "effect allele",
                            other_allele_col = "other allele",
                            eaf_col = "Freq",
                            pval_col = "P",
                            samplesize_col = "N")

str(IV_protein)
IV_protein$exposure <- "protein"

IV_protein <- IV_protein[IV_protein$pval.exposure < 5e-8,]
IV_protein <- clump_data(IV_protein)

IV_protein$r2 <- (2 * (IV_protein$beta.exposure^2) * IV_protein$eaf.exposure * (1 - IV_protein$eaf.exposure)) /
  (2 * (IV_protein$beta.exposure^2) * IV_protein$eaf.exposure * (1 - IV_protein$eaf.exposure) +
     2 * IV_protein$samplesize.exposure * IV_protein$eaf.exposure * (1 - IV_protein$eaf.exposure) * (IV_protein$se.exposure^2))

IV_protein$F <- (IV_protein$r2 * (IV_protein$samplesize.exposure - 2))/(1 - IV_protein$r2)
range(IV_protein$F) # [1]  34.84109 118.37356
IV_protein_meanF <- mean(IV_protein$F) # [1] 63.10193

write.xlsx(IV_protein, file = "Output/UVMR_IV/IV_protein.xlsx")

out_protein_ON_education <- extract_outcome_data(
  snps = IV_protein$SNP,
  outcomes = 'ieu-a-1239')
str(out_protein_ON_education)
out_protein_ON_education$outcome <- "educational attainment"

UVdat_protein_ON_education <- harmonise_data(
  exposure_dat =  IV_protein, 
  outcome_dat = out_protein_ON_education
)
str(UVdat_protein_ON_education)
range(UVdat_protein_ON_education$F)
range(UVdat_protein_ON_education$pval.outcome)
write.xlsx(UVdat_protein_ON_education,file = "Output/UVMR_IV/UVdat_protein_ON_education.xlsx")

result_protein_ON_education <- mr(UVdat_protein_ON_education)


scatter_plot_protein_ON_education <- mr_scatter_plot(result_protein_ON_education, UVdat_protein_ON_education)
scatter_plot_protein_ON_education[[1]]
ggsave(scatter_plot_protein_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_protein_ON_education.pdf", width=7, height=7)


pleiotropy_protein_ON_education <- mr_pleiotropy_test(UVdat_protein_ON_education)
str(pleiotropy_protein_ON_education)
pleiotropy_protein_ON_education <- pleiotropy_protein_ON_education[,-c(1,2)]


heterogeneity_protein_ON_education <- mr_heterogeneity(UVdat_protein_ON_education)
str(heterogeneity_protein_ON_education)
heterogeneity_protein_ON_education <- heterogeneity_protein_ON_education[,-c(1,2)]


singleSNP_protein_ON_education <- mr_singlesnp(UVdat_protein_ON_education)
funnel_plot_protein_ON_education <- mr_funnel_plot(singleSNP_protein_ON_education)
funnel_plot_protein_ON_education[[1]]
ggsave(funnel_plot_protein_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_protein_ON_education.pdf", width=7, height=7)


loo_protein_ON_education <- mr_leaveoneout(UVdat_protein_ON_education)
loo_plot_protein_ON_education <- mr_leaveoneout_plot(loo_protein_ON_education)
loo_plot_protein_ON_education[[1]]
ggsave(loo_plot_protein_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_protein_ON_education.pdf", width=7, height=21)


singleSNP_plot_protein_ON_education <- mr_forest_plot(singleSNP_protein_ON_education)
singleSNP_plot_protein_ON_education[[1]]
ggsave(singleSNP_plot_protein_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_protein_ON_education.pdf", width=7, height=7)

PRESSO_protein_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_protein_ON_education, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_protein_ON_education <- generate_odds_ratios(result_protein_ON_education)
result_protein_ON_education <- result_protein_ON_education[,-c(1,2)]


#### ________M1.1 Diet_fat ####
IV_fat <-format_data(GWAS_fat,
                         type = "exposure", 
                         snp_col = "SNP", 
                         beta_col = "BETA",
                         se_col = "SE",
                         effect_allele_col = "effect allele",
                         other_allele_col = "other allele",
                         eaf_col = "Freq",
                         pval_col = "P",
                         samplesize_col = "N")

str(IV_fat)
IV_fat$exposure <- "fat"

IV_fat <- IV_fat[IV_fat$pval.exposure < 5e-8,]
IV_fat <- clump_data(IV_fat)

IV_fat$r2 <- (2 * (IV_fat$beta.exposure^2) * IV_fat$eaf.exposure * (1 - IV_fat$eaf.exposure)) /
  (2 * (IV_fat$beta.exposure^2) * IV_fat$eaf.exposure * (1 - IV_fat$eaf.exposure) +
     2 * IV_fat$samplesize.exposure * IV_fat$eaf.exposure * (1 - IV_fat$eaf.exposure) * (IV_fat$se.exposure^2))

IV_fat$F <- (IV_fat$r2 * (IV_fat$samplesize.exposure - 2))/(1 - IV_fat$r2)
range(IV_fat$F) # [1] 35.12706 99.89831
IV_fat_meanF <- mean(IV_fat$F) # [1] 69.0989

write.xlsx(IV_fat, file = "Output/UVMR_IV/IV_fat.xlsx")

out_fat_ON_education <- extract_outcome_data(
  snps = IV_fat$SNP,
  outcomes = 'ieu-a-1239')
str(out_fat_ON_education)
out_fat_ON_education$outcome <- "educational attainment"

UVdat_fat_ON_education <- harmonise_data(
  exposure_dat =  IV_fat, 
  outcome_dat = out_fat_ON_education
)
str(UVdat_fat_ON_education)
range(UVdat_fat_ON_education$F)
range(UVdat_fat_ON_education$pval.outcome)
write.xlsx(UVdat_fat_ON_education,file = "Output/UVMR_IV/UVdat_fat_ON_education.xlsx")

result_fat_ON_education <- mr(UVdat_fat_ON_education)


scatter_plot_fat_ON_education <- mr_scatter_plot(result_fat_ON_education, UVdat_fat_ON_education)
scatter_plot_fat_ON_education[[1]]
ggsave(scatter_plot_fat_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_fat_ON_education.pdf", width=7, height=7)


pleiotropy_fat_ON_education <- mr_pleiotropy_test(UVdat_fat_ON_education)
str(pleiotropy_fat_ON_education)
pleiotropy_fat_ON_education <- pleiotropy_fat_ON_education[,-c(1,2)]


heterogeneity_fat_ON_education <- mr_heterogeneity(UVdat_fat_ON_education)
str(heterogeneity_fat_ON_education)
heterogeneity_fat_ON_education <- heterogeneity_fat_ON_education[,-c(1,2)]


singleSNP_fat_ON_education <- mr_singlesnp(UVdat_fat_ON_education)
funnel_plot_fat_ON_education <- mr_funnel_plot(singleSNP_fat_ON_education)
funnel_plot_fat_ON_education[[1]]
ggsave(funnel_plot_fat_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_fat_ON_education.pdf", width=7, height=7)


loo_fat_ON_education <- mr_leaveoneout(UVdat_fat_ON_education)
loo_plot_fat_ON_education <- mr_leaveoneout_plot(loo_fat_ON_education)
loo_plot_fat_ON_education[[1]]
ggsave(loo_plot_fat_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_fat_ON_education.pdf", width=7, height=21)


singleSNP_plot_fat_ON_education <- mr_forest_plot(singleSNP_fat_ON_education)
singleSNP_plot_fat_ON_education[[1]]
ggsave(singleSNP_plot_fat_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_fat_ON_education.pdf", width=7, height=7)

PRESSO_fat_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_fat_ON_education, NbDistribution = 10000,  
                                         SignifThreshold = 0.05)

result_fat_ON_education <- generate_odds_ratios(result_fat_ON_education)
result_fat_ON_education <- result_fat_ON_education[,-c(1,2)]


#### ________M1.1 Diet_FreshFruit ####
IV_FreshFruit <-format_data(GWAS_FreshFruit,
                     type = "exposure", 
                     snp_col = "SNP", 
                     beta_col = "beta",
                     se_col = "standard_error",
                     effect_allele_col = "effect_allele",
                     other_allele_col = "other_allele",
                     eaf_col = "effect_allele_frequency",
                     pval_col = "p_value",
                     samplesize_col = "SS")

str(IV_FreshFruit)
IV_FreshFruit$exposure <- "FreshFruit"

IV_FreshFruit <- IV_FreshFruit[IV_FreshFruit$pval.exposure < 5e-8,]
IV_FreshFruit <- clump_data(IV_FreshFruit)

IV_FreshFruit$r2 <- (2 * (IV_FreshFruit$beta.exposure^2) * IV_FreshFruit$eaf.exposure * (1 - IV_FreshFruit$eaf.exposure)) /
  (2 * (IV_FreshFruit$beta.exposure^2) * IV_FreshFruit$eaf.exposure * (1 - IV_FreshFruit$eaf.exposure) +
     2 * IV_FreshFruit$samplesize.exposure * IV_FreshFruit$eaf.exposure * (1 - IV_FreshFruit$eaf.exposure) * (IV_FreshFruit$se.exposure^2))

IV_FreshFruit$F <- (IV_FreshFruit$r2 * (IV_FreshFruit$samplesize.exposure - 2))/(1 - IV_FreshFruit$r2)
range(IV_FreshFruit$F) # [1]  29.07646 290.02630
IV_FreshFruit_meanF <- mean(IV_FreshFruit$F) # [1] 47.20986

write.xlsx(IV_FreshFruit, file = "Output/UVMR_IV/IV_FreshFruit.xlsx")

out_FreshFruit_ON_education <- extract_outcome_data(
  snps = IV_FreshFruit$SNP,
  outcomes = 'ieu-a-1239')
str(out_FreshFruit_ON_education)
out_FreshFruit_ON_education$outcome <- "educational attainment"

UVdat_FreshFruit_ON_education <- harmonise_data(
  exposure_dat =  IV_FreshFruit, 
  outcome_dat = out_FreshFruit_ON_education
)
str(UVdat_FreshFruit_ON_education)
range(UVdat_FreshFruit_ON_education$F)
range(UVdat_FreshFruit_ON_education$pval.outcome)
UVdat_FreshFruit_ON_education <- UVdat_FreshFruit_ON_education[UVdat_FreshFruit_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_FreshFruit_ON_education,file = "Output/UVMR_IV/UVdat_FreshFruit_ON_education.xlsx")

result_FreshFruit_ON_education <- mr(UVdat_FreshFruit_ON_education)


scatter_plot_FreshFruit_ON_education <- mr_scatter_plot(result_FreshFruit_ON_education, UVdat_FreshFruit_ON_education)
scatter_plot_FreshFruit_ON_education[[1]]
ggsave(scatter_plot_FreshFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_FreshFruit_ON_education.pdf", width=7, height=7)


pleiotropy_FreshFruit_ON_education <- mr_pleiotropy_test(UVdat_FreshFruit_ON_education)
str(pleiotropy_FreshFruit_ON_education)
pleiotropy_FreshFruit_ON_education <- pleiotropy_FreshFruit_ON_education[,-c(1,2)]


heterogeneity_FreshFruit_ON_education <- mr_heterogeneity(UVdat_FreshFruit_ON_education)
str(heterogeneity_FreshFruit_ON_education)
heterogeneity_FreshFruit_ON_education <- heterogeneity_FreshFruit_ON_education[,-c(1,2)]


singleSNP_FreshFruit_ON_education <- mr_singlesnp(UVdat_FreshFruit_ON_education)
funnel_plot_FreshFruit_ON_education <- mr_funnel_plot(singleSNP_FreshFruit_ON_education)
funnel_plot_FreshFruit_ON_education[[1]]
ggsave(funnel_plot_FreshFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_FreshFruit_ON_education.pdf", width=7, height=7)


loo_FreshFruit_ON_education <- mr_leaveoneout(UVdat_FreshFruit_ON_education)
loo_plot_FreshFruit_ON_education <- mr_leaveoneout_plot(loo_FreshFruit_ON_education)
loo_plot_FreshFruit_ON_education[[1]]
ggsave(loo_plot_FreshFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_FreshFruit_ON_education.pdf", width=7, height=21)


singleSNP_plot_FreshFruit_ON_education <- mr_forest_plot(singleSNP_FreshFruit_ON_education)
singleSNP_plot_FreshFruit_ON_education[[1]]
ggsave(singleSNP_plot_FreshFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_FreshFruit_ON_education.pdf", width=7, height=7)

PRESSO_FreshFruit_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_FreshFruit_ON_education, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_FreshFruit_ON_education <- generate_odds_ratios(result_FreshFruit_ON_education)
result_FreshFruit_ON_education <- result_FreshFruit_ON_education[,-c(1,2)]


#### ________M1.1 Diet_DriedFruit ####
IV_DriedFruit <-format_data(GWAS_DriedFruit,
                            type = "exposure", 
                            snp_col = "SNP", 
                            beta_col = "beta",
                            se_col = "standard_error",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele",
                            eaf_col = "effect_allele_frequency",
                            pval_col = "p_value",
                            samplesize_col = "SS")

str(IV_DriedFruit)
IV_DriedFruit$exposure <- "DriedFruit"

IV_DriedFruit <- IV_DriedFruit[IV_DriedFruit$pval.exposure < 5e-8,]
IV_DriedFruit <- clump_data(IV_DriedFruit)

IV_DriedFruit$r2 <- (2 * (IV_DriedFruit$beta.exposure^2) * IV_DriedFruit$eaf.exposure * (1 - IV_DriedFruit$eaf.exposure)) /
  (2 * (IV_DriedFruit$beta.exposure^2) * IV_DriedFruit$eaf.exposure * (1 - IV_DriedFruit$eaf.exposure) +
     2 * IV_DriedFruit$samplesize.exposure * IV_DriedFruit$eaf.exposure * (1 - IV_DriedFruit$eaf.exposure) * (IV_DriedFruit$se.exposure^2))

IV_DriedFruit$F <- (IV_DriedFruit$r2 * (IV_DriedFruit$samplesize.exposure - 2))/(1 - IV_DriedFruit$r2)
range(IV_DriedFruit$F) #[1]  29.65542 137.01605
IV_DriedFruit_meanF <- mean(IV_DriedFruit$F) #[1] 40.86215

write.xlsx(IV_DriedFruit, file = "Output/UVMR_IV/IV_DriedFruit.xlsx")

out_DriedFruit_ON_education <- extract_outcome_data(
  snps = IV_DriedFruit$SNP,
  outcomes = 'ieu-a-1239')
str(out_DriedFruit_ON_education)
out_DriedFruit_ON_education$outcome <- "educational attainment"

UVdat_DriedFruit_ON_education <- harmonise_data(
  exposure_dat =  IV_DriedFruit, 
  outcome_dat = out_DriedFruit_ON_education
)
str(UVdat_DriedFruit_ON_education)
range(UVdat_DriedFruit_ON_education$F)
range(UVdat_DriedFruit_ON_education$pval.outcome)
UVdat_DriedFruit_ON_education <- UVdat_DriedFruit_ON_education[UVdat_DriedFruit_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_DriedFruit_ON_education,file = "Output/UVMR_IV/UVdat_DriedFruit_ON_education.xlsx")

result_DriedFruit_ON_education <- mr(UVdat_DriedFruit_ON_education)


scatter_plot_DriedFruit_ON_education <- mr_scatter_plot(result_DriedFruit_ON_education, UVdat_DriedFruit_ON_education)
scatter_plot_DriedFruit_ON_education[[1]]
ggsave(scatter_plot_DriedFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DriedFruit_ON_education.pdf", width=7, height=7)


pleiotropy_DriedFruit_ON_education <- mr_pleiotropy_test(UVdat_DriedFruit_ON_education)
str(pleiotropy_DriedFruit_ON_education)
pleiotropy_DriedFruit_ON_education <- pleiotropy_DriedFruit_ON_education[,-c(1,2)]


heterogeneity_DriedFruit_ON_education <- mr_heterogeneity(UVdat_DriedFruit_ON_education)
str(heterogeneity_DriedFruit_ON_education)
heterogeneity_DriedFruit_ON_education <- heterogeneity_DriedFruit_ON_education[,-c(1,2)]


singleSNP_DriedFruit_ON_education <- mr_singlesnp(UVdat_DriedFruit_ON_education)
funnel_plot_DriedFruit_ON_education <- mr_funnel_plot(singleSNP_DriedFruit_ON_education)
funnel_plot_DriedFruit_ON_education[[1]]
ggsave(funnel_plot_DriedFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DriedFruit_ON_education.pdf", width=7, height=7)


loo_DriedFruit_ON_education <- mr_leaveoneout(UVdat_DriedFruit_ON_education)
loo_plot_DriedFruit_ON_education <- mr_leaveoneout_plot(loo_DriedFruit_ON_education)
loo_plot_DriedFruit_ON_education[[1]]
ggsave(loo_plot_DriedFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_DriedFruit_ON_education.pdf", width=7, height=7)


singleSNP_plot_DriedFruit_ON_education <- mr_forest_plot(singleSNP_DriedFruit_ON_education)
singleSNP_plot_DriedFruit_ON_education[[1]]
ggsave(singleSNP_plot_DriedFruit_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DriedFruit_ON_education.pdf", width=7, height=7)

PRESSO_DriedFruit_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DriedFruit_ON_education, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_DriedFruit_ON_education <- generate_odds_ratios(result_DriedFruit_ON_education)
result_DriedFruit_ON_education <- result_DriedFruit_ON_education[,-c(1,2)]


#### ________M1.1 Diet_RawVegetables ####
IV_RawVegetables <-format_data(GWAS_RawVegetables,
                            type = "exposure", 
                            snp_col = "SNP", 
                            beta_col = "beta",
                            se_col = "standard_error",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele",
                            eaf_col = "effect_allele_frequency",
                            pval_col = "p_value",
                            samplesize_col = "SS")

str(IV_RawVegetables)
IV_RawVegetables$exposure <- "RawVegetables"

IV_RawVegetables <- IV_RawVegetables[IV_RawVegetables$pval.exposure < 5e-8,]
IV_RawVegetables <- clump_data(IV_RawVegetables)

IV_RawVegetables$r2 <- (2 * (IV_RawVegetables$beta.exposure^2) * IV_RawVegetables$eaf.exposure * (1 - IV_RawVegetables$eaf.exposure)) /
  (2 * (IV_RawVegetables$beta.exposure^2) * IV_RawVegetables$eaf.exposure * (1 - IV_RawVegetables$eaf.exposure) +
     2 * IV_RawVegetables$samplesize.exposure * IV_RawVegetables$eaf.exposure * (1 - IV_RawVegetables$eaf.exposure) * (IV_RawVegetables$se.exposure^2))

IV_RawVegetables$F <- (IV_RawVegetables$r2 * (IV_RawVegetables$samplesize.exposure - 2))/(1 - IV_RawVegetables$r2)
range(IV_RawVegetables$F) #[1] 29.65811 93.01331
IV_RawVegetables_meanF <- mean(IV_RawVegetables$F) #[1] 41.02934

write.xlsx(IV_RawVegetables, file = "Output/UVMR_IV/IV_RawVegetables.xlsx")

out_RawVegetables_ON_education <- extract_outcome_data(
  snps = IV_RawVegetables$SNP,
  outcomes = 'ieu-a-1239')
str(out_RawVegetables_ON_education)
out_RawVegetables_ON_education$outcome <- "educational attainment"

UVdat_RawVegetables_ON_education <- harmonise_data(
  exposure_dat =  IV_RawVegetables, 
  outcome_dat = out_RawVegetables_ON_education
)
str(UVdat_RawVegetables_ON_education)
range(UVdat_RawVegetables_ON_education$F)
range(UVdat_RawVegetables_ON_education$pval.outcome)
UVdat_RawVegetables_ON_education <- UVdat_RawVegetables_ON_education[UVdat_RawVegetables_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_RawVegetables_ON_education,file = "Output/UVMR_IV/UVdat_RawVegetables_ON_education.xlsx")

result_RawVegetables_ON_education <- mr(UVdat_RawVegetables_ON_education)


scatter_plot_RawVegetables_ON_education <- mr_scatter_plot(result_RawVegetables_ON_education, UVdat_RawVegetables_ON_education)
scatter_plot_RawVegetables_ON_education[[1]]
ggsave(scatter_plot_RawVegetables_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_RawVegetables_ON_education.pdf", width=7, height=7)


pleiotropy_RawVegetables_ON_education <- mr_pleiotropy_test(UVdat_RawVegetables_ON_education)
str(pleiotropy_RawVegetables_ON_education)
pleiotropy_RawVegetables_ON_education <- pleiotropy_RawVegetables_ON_education[,-c(1,2)]


heterogeneity_RawVegetables_ON_education <- mr_heterogeneity(UVdat_RawVegetables_ON_education)
str(heterogeneity_RawVegetables_ON_education)
heterogeneity_RawVegetables_ON_education <- heterogeneity_RawVegetables_ON_education[,-c(1,2)]


singleSNP_RawVegetables_ON_education <- mr_singlesnp(UVdat_RawVegetables_ON_education)
funnel_plot_RawVegetables_ON_education <- mr_funnel_plot(singleSNP_RawVegetables_ON_education)
funnel_plot_RawVegetables_ON_education[[1]]
ggsave(funnel_plot_RawVegetables_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_RawVegetables_ON_education.pdf", width=7, height=7)


loo_RawVegetables_ON_education <- mr_leaveoneout(UVdat_RawVegetables_ON_education)
loo_plot_RawVegetables_ON_education <- mr_leaveoneout_plot(loo_RawVegetables_ON_education)
loo_plot_RawVegetables_ON_education[[1]]
ggsave(loo_plot_RawVegetables_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_RawVegetables_ON_education.pdf", width=7, height=7)


singleSNP_plot_RawVegetables_ON_education <- mr_forest_plot(singleSNP_RawVegetables_ON_education)
singleSNP_plot_RawVegetables_ON_education[[1]]
ggsave(singleSNP_plot_RawVegetables_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_RawVegetables_ON_education.pdf", width=7, height=7)

PRESSO_RawVegetables_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_RawVegetables_ON_education, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_RawVegetables_ON_education <- generate_odds_ratios(result_RawVegetables_ON_education)
result_RawVegetables_ON_education <- result_RawVegetables_ON_education[,-c(1,2)]


#### ________M1.1 Diet_salt ####
IV_salt <-format_data(GWAS_salt,
                               type = "exposure", 
                               snp_col = "SNP", 
                               beta_col = "beta",
                               se_col = "standard_error",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               eaf_col = "effect_allele_frequency",
                               pval_col = "p_value",
                               samplesize_col = "SS")

str(IV_salt)
IV_salt$exposure <- "salt"

IV_salt <- IV_salt[IV_salt$pval.exposure < 5e-8,]
IV_salt <- clump_data(IV_salt)

IV_salt$r2 <- (2 * (IV_salt$beta.exposure^2) * IV_salt$eaf.exposure * (1 - IV_salt$eaf.exposure)) /
  (2 * (IV_salt$beta.exposure^2) * IV_salt$eaf.exposure * (1 - IV_salt$eaf.exposure) +
     2 * IV_salt$samplesize.exposure * IV_salt$eaf.exposure * (1 - IV_salt$eaf.exposure) * (IV_salt$se.exposure^2))

IV_salt$F <- (IV_salt$r2 * (IV_salt$samplesize.exposure - 2))/(1 - IV_salt$r2)
range(IV_salt$F) #[1]  29.40288 194.82851
IV_salt_meanF <- mean(IV_salt$F) #[1] 48.83646

write.xlsx(IV_salt, file = "Output/UVMR_IV/IV_salt.xlsx")

out_salt_ON_education <- extract_outcome_data(
  snps = IV_salt$SNP,
  outcomes = 'ieu-a-1239')
str(out_salt_ON_education)
out_salt_ON_education$outcome <- "educational attainment"

UVdat_salt_ON_education <- harmonise_data(
  exposure_dat =  IV_salt, 
  outcome_dat = out_salt_ON_education
)
str(UVdat_salt_ON_education)
range(UVdat_salt_ON_education$F)
range(UVdat_salt_ON_education$pval.outcome)
UVdat_salt_ON_education <- UVdat_salt_ON_education[UVdat_salt_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_salt_ON_education,file = "Output/UVMR_IV/UVdat_salt_ON_education.xlsx")

result_salt_ON_education <- mr(UVdat_salt_ON_education)


scatter_plot_salt_ON_education <- mr_scatter_plot(result_salt_ON_education, UVdat_salt_ON_education)
scatter_plot_salt_ON_education[[1]]
ggsave(scatter_plot_salt_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_salt_ON_education.pdf", width=7, height=7)


pleiotropy_salt_ON_education <- mr_pleiotropy_test(UVdat_salt_ON_education)
str(pleiotropy_salt_ON_education)
pleiotropy_salt_ON_education <- pleiotropy_salt_ON_education[,-c(1,2)]


heterogeneity_salt_ON_education <- mr_heterogeneity(UVdat_salt_ON_education)
str(heterogeneity_salt_ON_education)
heterogeneity_salt_ON_education <- heterogeneity_salt_ON_education[,-c(1,2)]


singleSNP_salt_ON_education <- mr_singlesnp(UVdat_salt_ON_education)
funnel_plot_salt_ON_education <- mr_funnel_plot(singleSNP_salt_ON_education)
funnel_plot_salt_ON_education[[1]]
ggsave(funnel_plot_salt_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_salt_ON_education.pdf", width=7, height=7)


loo_salt_ON_education <- mr_leaveoneout(UVdat_salt_ON_education)
loo_plot_salt_ON_education <- mr_leaveoneout_plot(loo_salt_ON_education)
loo_plot_salt_ON_education[[1]]
ggsave(loo_plot_salt_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_salt_ON_education.pdf", width=7, height=7)


singleSNP_plot_salt_ON_education <- mr_forest_plot(singleSNP_salt_ON_education)
singleSNP_plot_salt_ON_education[[1]]
ggsave(singleSNP_plot_salt_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_salt_ON_education.pdf", width=7, height=7)

PRESSO_salt_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_salt_ON_education, NbDistribution = 10000,  
                                               SignifThreshold = 0.05)

result_salt_ON_education <- generate_odds_ratios(result_salt_ON_education)
result_salt_ON_education <- result_salt_ON_education[,-c(1,2)]


#### ________M1.1 Diet_processed ####
IV_processed <-format_data(GWAS_processed,
                      type = "exposure", 
                      snp_col = "SNP", 
                      beta_col = "beta",
                      se_col = "standard_error",
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      eaf_col = "effect_allele_frequency",
                      pval_col = "p_value",
                      samplesize_col = "SS")

str(IV_processed)
IV_processed$exposure <- "processed"

IV_processed <- IV_processed[IV_processed$pval.exposure < 5e-8,]
IV_processed <- clump_data(IV_processed)

IV_processed$r2 <- (2 * (IV_processed$beta.exposure^2) * IV_processed$eaf.exposure * (1 - IV_processed$eaf.exposure)) /
  (2 * (IV_processed$beta.exposure^2) * IV_processed$eaf.exposure * (1 - IV_processed$eaf.exposure) +
     2 * IV_processed$samplesize.exposure * IV_processed$eaf.exposure * (1 - IV_processed$eaf.exposure) * (IV_processed$se.exposure^2))

IV_processed$F <- (IV_processed$r2 * (IV_processed$samplesize.exposure - 2))/(1 - IV_processed$r2)
range(IV_processed$F) #[1] 30.15430 52.72489
IV_processed_meanF <- mean(IV_processed$F) #[1] 35.84384

write.xlsx(IV_processed, file = "Output/UVMR_IV/IV_processed.xlsx")

out_processed_ON_education <- extract_outcome_data(
  snps = IV_processed$SNP,
  outcomes = 'ieu-a-1239')
str(out_processed_ON_education)
out_processed_ON_education$outcome <- "educational attainment"

UVdat_processed_ON_education <- harmonise_data(
  exposure_dat =  IV_processed, 
  outcome_dat = out_processed_ON_education
)
str(UVdat_processed_ON_education)
range(UVdat_processed_ON_education$F)
range(UVdat_processed_ON_education$pval.outcome)
UVdat_processed_ON_education <- UVdat_processed_ON_education[UVdat_processed_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_processed_ON_education,file = "Output/UVMR_IV/UVdat_processed_ON_education.xlsx")

result_processed_ON_education <- mr(UVdat_processed_ON_education)


scatter_plot_processed_ON_education <- mr_scatter_plot(result_processed_ON_education, UVdat_processed_ON_education)
scatter_plot_processed_ON_education[[1]]
ggsave(scatter_plot_processed_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_processed_ON_education.pdf", width=7, height=7)


pleiotropy_processed_ON_education <- mr_pleiotropy_test(UVdat_processed_ON_education)
str(pleiotropy_processed_ON_education)
pleiotropy_processed_ON_education <- pleiotropy_processed_ON_education[,-c(1,2)]


heterogeneity_processed_ON_education <- mr_heterogeneity(UVdat_processed_ON_education)
str(heterogeneity_processed_ON_education)
heterogeneity_processed_ON_education <- heterogeneity_processed_ON_education[,-c(1,2)]


singleSNP_processed_ON_education <- mr_singlesnp(UVdat_processed_ON_education)
funnel_plot_processed_ON_education <- mr_funnel_plot(singleSNP_processed_ON_education)
funnel_plot_processed_ON_education[[1]]
ggsave(funnel_plot_processed_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_processed_ON_education.pdf", width=7, height=7)


loo_processed_ON_education <- mr_leaveoneout(UVdat_processed_ON_education)
loo_plot_processed_ON_education <- mr_leaveoneout_plot(loo_processed_ON_education)
loo_plot_processed_ON_education[[1]]
ggsave(loo_plot_processed_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_processed_ON_education.pdf", width=7, height=7)


singleSNP_plot_processed_ON_education <- mr_forest_plot(singleSNP_processed_ON_education)
singleSNP_plot_processed_ON_education[[1]]
ggsave(singleSNP_plot_processed_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_processed_ON_education.pdf", width=7, height=7)

PRESSO_processed_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_processed_ON_education, NbDistribution = 10000,  
                                      SignifThreshold = 0.05)

result_processed_ON_education <- generate_odds_ratios(result_processed_ON_education)
result_processed_ON_education <- result_processed_ON_education[,-c(1,2)]


#### ________M1.1 Diet_OilyFish ####
IV_OilyFish <-format_data(GWAS_OilyFish,
                           type = "exposure", 
                           snp_col = "SNP", 
                           beta_col = "beta",
                           se_col = "standard_error",
                           effect_allele_col = "effect_allele",
                           other_allele_col = "other_allele",
                           eaf_col = "effect_allele_frequency",
                           pval_col = "p_value",
                           samplesize_col = "SS")

str(IV_OilyFish)
IV_OilyFish$exposure <- "OilyFish"

IV_OilyFish <- IV_OilyFish[IV_OilyFish$pval.exposure < 5e-8,]
IV_OilyFish <- clump_data(IV_OilyFish)

IV_OilyFish$r2 <- (2 * (IV_OilyFish$beta.exposure^2) * IV_OilyFish$eaf.exposure * (1 - IV_OilyFish$eaf.exposure)) /
  (2 * (IV_OilyFish$beta.exposure^2) * IV_OilyFish$eaf.exposure * (1 - IV_OilyFish$eaf.exposure) +
     2 * IV_OilyFish$samplesize.exposure * IV_OilyFish$eaf.exposure * (1 - IV_OilyFish$eaf.exposure) * (IV_OilyFish$se.exposure^2))

IV_OilyFish$F <- (IV_OilyFish$r2 * (IV_OilyFish$samplesize.exposure - 2))/(1 - IV_OilyFish$r2)
range(IV_OilyFish$F) #[1] 29.54126 98.56254
IV_OilyFish_meanF <- mean(IV_OilyFish$F) #[1] 45.04969

write.xlsx(IV_OilyFish, file = "Output/UVMR_IV/IV_OilyFish.xlsx")

out_OilyFish_ON_education <- extract_outcome_data(
  snps = IV_OilyFish$SNP,
  outcomes = 'ieu-a-1239')
str(out_OilyFish_ON_education)
out_OilyFish_ON_education$outcome <- "educational attainment"

UVdat_OilyFish_ON_education <- harmonise_data(
  exposure_dat =  IV_OilyFish, 
  outcome_dat = out_OilyFish_ON_education
)
str(UVdat_OilyFish_ON_education)
range(UVdat_OilyFish_ON_education$F)
range(UVdat_OilyFish_ON_education$pval.outcome)
UVdat_OilyFish_ON_education <- UVdat_OilyFish_ON_education[UVdat_OilyFish_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_OilyFish_ON_education,file = "Output/UVMR_IV/UVdat_OilyFish_ON_education.xlsx")

result_OilyFish_ON_education <- mr(UVdat_OilyFish_ON_education)


scatter_plot_OilyFish_ON_education <- mr_scatter_plot(result_OilyFish_ON_education, UVdat_OilyFish_ON_education)
scatter_plot_OilyFish_ON_education[[1]]
ggsave(scatter_plot_OilyFish_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_OilyFish_ON_education.pdf", width=7, height=7)


pleiotropy_OilyFish_ON_education <- mr_pleiotropy_test(UVdat_OilyFish_ON_education)
str(pleiotropy_OilyFish_ON_education)
pleiotropy_OilyFish_ON_education <- pleiotropy_OilyFish_ON_education[,-c(1,2)]


heterogeneity_OilyFish_ON_education <- mr_heterogeneity(UVdat_OilyFish_ON_education)
str(heterogeneity_OilyFish_ON_education)
heterogeneity_OilyFish_ON_education <- heterogeneity_OilyFish_ON_education[,-c(1,2)]


singleSNP_OilyFish_ON_education <- mr_singlesnp(UVdat_OilyFish_ON_education)
funnel_plot_OilyFish_ON_education <- mr_funnel_plot(singleSNP_OilyFish_ON_education)
funnel_plot_OilyFish_ON_education[[1]]
ggsave(funnel_plot_OilyFish_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_OilyFish_ON_education.pdf", width=7, height=7)


loo_OilyFish_ON_education <- mr_leaveoneout(UVdat_OilyFish_ON_education)
loo_plot_OilyFish_ON_education <- mr_leaveoneout_plot(loo_OilyFish_ON_education)
loo_plot_OilyFish_ON_education[[1]]
ggsave(loo_plot_OilyFish_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_OilyFish_ON_education.pdf", width=7, height=7)


singleSNP_plot_OilyFish_ON_education <- mr_forest_plot(singleSNP_OilyFish_ON_education)
singleSNP_plot_OilyFish_ON_education[[1]]
ggsave(singleSNP_plot_OilyFish_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_OilyFish_ON_education.pdf", width=7, height=7)

PRESSO_OilyFish_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                           OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_OilyFish_ON_education, NbDistribution = 10000,  
                                           SignifThreshold = 0.05)

result_OilyFish_ON_education <- generate_odds_ratios(result_OilyFish_ON_education)
result_OilyFish_ON_education <- result_OilyFish_ON_education[,-c(1,2)]


#### ________M1.1 Diet_wholegrain ####
IV_wholegrain <-format_data(GWAS_wholegrain,
                          type = "exposure", 
                          snp_col = "SNP", 
                          beta_col = "beta",
                          se_col = "standard_error",
                          effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele",
                          eaf_col = "effect_allele_frequency",
                          pval_col = "p_value",
                          samplesize_col = "SS")

str(IV_wholegrain)
IV_wholegrain$exposure <- "wholegrain"

IV_wholegrain <- IV_wholegrain[IV_wholegrain$pval.exposure < 5e-8,]
IV_wholegrain <- clump_data(IV_wholegrain)

IV_wholegrain$r2 <- (2 * (IV_wholegrain$beta.exposure^2) * IV_wholegrain$eaf.exposure * (1 - IV_wholegrain$eaf.exposure)) /
  (2 * (IV_wholegrain$beta.exposure^2) * IV_wholegrain$eaf.exposure * (1 - IV_wholegrain$eaf.exposure) +
     2 * IV_wholegrain$samplesize.exposure * IV_wholegrain$eaf.exposure * (1 - IV_wholegrain$eaf.exposure) * (IV_wholegrain$se.exposure^2))

IV_wholegrain$F <- (IV_wholegrain$r2 * (IV_wholegrain$samplesize.exposure - 2))/(1 - IV_wholegrain$r2)
range(IV_wholegrain$F) #[1] 29.51506 84.78115
IV_wholegrain_meanF <- mean(IV_wholegrain$F) #[1] 37.34682

write.xlsx(IV_wholegrain, file = "Output/UVMR_IV/IV_wholegrain.xlsx")

out_wholegrain_ON_education <- extract_outcome_data(
  snps = IV_wholegrain$SNP,
  outcomes = 'ieu-a-1239')
str(out_wholegrain_ON_education)
out_wholegrain_ON_education$outcome <- "educational attainment"

UVdat_wholegrain_ON_education <- harmonise_data(
  exposure_dat =  IV_wholegrain, 
  outcome_dat = out_wholegrain_ON_education
)
str(UVdat_wholegrain_ON_education)
range(UVdat_wholegrain_ON_education$F)
range(UVdat_wholegrain_ON_education$pval.outcome)
UVdat_wholegrain_ON_education <- UVdat_wholegrain_ON_education[UVdat_wholegrain_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_wholegrain_ON_education,file = "Output/UVMR_IV/UVdat_wholegrain_ON_education.xlsx")

result_wholegrain_ON_education <- mr(UVdat_wholegrain_ON_education)


scatter_plot_wholegrain_ON_education <- mr_scatter_plot(result_wholegrain_ON_education, UVdat_wholegrain_ON_education)
scatter_plot_wholegrain_ON_education[[1]]
ggsave(scatter_plot_wholegrain_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_wholegrain_ON_education.pdf", width=7, height=7)


pleiotropy_wholegrain_ON_education <- mr_pleiotropy_test(UVdat_wholegrain_ON_education)
str(pleiotropy_wholegrain_ON_education)
pleiotropy_wholegrain_ON_education <- pleiotropy_wholegrain_ON_education[,-c(1,2)]


heterogeneity_wholegrain_ON_education <- mr_heterogeneity(UVdat_wholegrain_ON_education)
str(heterogeneity_wholegrain_ON_education)
heterogeneity_wholegrain_ON_education <- heterogeneity_wholegrain_ON_education[,-c(1,2)]


singleSNP_wholegrain_ON_education <- mr_singlesnp(UVdat_wholegrain_ON_education)
funnel_plot_wholegrain_ON_education <- mr_funnel_plot(singleSNP_wholegrain_ON_education)
funnel_plot_wholegrain_ON_education[[1]]
ggsave(funnel_plot_wholegrain_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_wholegrain_ON_education.pdf", width=7, height=7)


loo_wholegrain_ON_education <- mr_leaveoneout(UVdat_wholegrain_ON_education)
loo_plot_wholegrain_ON_education <- mr_leaveoneout_plot(loo_wholegrain_ON_education)
loo_plot_wholegrain_ON_education[[1]]
ggsave(loo_plot_wholegrain_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_wholegrain_ON_education.pdf", width=7, height=7)


singleSNP_plot_wholegrain_ON_education <- mr_forest_plot(singleSNP_wholegrain_ON_education)
singleSNP_plot_wholegrain_ON_education[[1]]
ggsave(singleSNP_plot_wholegrain_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_wholegrain_ON_education.pdf", width=7, height=7)

PRESSO_wholegrain_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_wholegrain_ON_education, NbDistribution = 10000,  
                                          SignifThreshold = 0.05)

result_wholegrain_ON_education <- generate_odds_ratios(result_wholegrain_ON_education)
result_wholegrain_ON_education <- result_wholegrain_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceOverallActivity ####
IV_DeviceOverallActivity <-format_data(GWAS_DeviceOverallActivity,
                          type = "exposure", 
                          snp_col = "SNP", 
                          beta_col = "BETA",
                          se_col = "SE",
                          effect_allele_col = "ALLELE1",
                          other_allele_col = "ALLELE0",
                          eaf_col = "A1FREQ",
                          pval_col = "P_BOLT_LMM_INF",
                          samplesize_col = "N")

str(IV_DeviceOverallActivity)
IV_DeviceOverallActivity$exposure <- "DeviceOverallActivity"

IV_DeviceOverallActivity <- IV_DeviceOverallActivity[IV_DeviceOverallActivity$pval.exposure < 5e-8,]
IV_DeviceOverallActivity <- clump_data(IV_DeviceOverallActivity)

IV_DeviceOverallActivity$r2 <- (2 * (IV_DeviceOverallActivity$beta.exposure^2) * IV_DeviceOverallActivity$eaf.exposure * (1 - IV_DeviceOverallActivity$eaf.exposure)) /
  (2 * (IV_DeviceOverallActivity$beta.exposure^2) * IV_DeviceOverallActivity$eaf.exposure * (1 - IV_DeviceOverallActivity$eaf.exposure) +
     2 * IV_DeviceOverallActivity$samplesize.exposure * IV_DeviceOverallActivity$eaf.exposure * (1 - IV_DeviceOverallActivity$eaf.exposure) * (IV_DeviceOverallActivity$se.exposure^2))

IV_DeviceOverallActivity$F <- (IV_DeviceOverallActivity$r2 * (IV_DeviceOverallActivity$samplesize.exposure - 2))/(1 - IV_DeviceOverallActivity$r2)
range(IV_DeviceOverallActivity$F) #[1] 29.75679 47.25188
IV_DeviceOverallActivity_meanF <- mean(IV_DeviceOverallActivity$F) #[1] 34.85534

write.xlsx(IV_DeviceOverallActivity, file = "Output/UVMR_IV/IV_DeviceOverallActivity.xlsx")

out_DeviceOverallActivity_ON_education <- extract_outcome_data(
  snps = IV_DeviceOverallActivity$SNP,
  outcomes = 'ieu-a-1239')
str(out_DeviceOverallActivity_ON_education)
out_DeviceOverallActivity_ON_education$outcome <- "educational attainment"

UVdat_DeviceOverallActivity_ON_education <- harmonise_data(
  exposure_dat =  IV_DeviceOverallActivity, 
  outcome_dat = out_DeviceOverallActivity_ON_education
)
str(UVdat_DeviceOverallActivity_ON_education)
range(UVdat_DeviceOverallActivity_ON_education$F)
range(UVdat_DeviceOverallActivity_ON_education$pval.outcome)
UVdat_DeviceOverallActivity_ON_education <- UVdat_DeviceOverallActivity_ON_education[UVdat_DeviceOverallActivity_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_DeviceOverallActivity_ON_education,file = "Output/UVMR_IV/UVdat_DeviceOverallActivity_ON_education.xlsx")

result_DeviceOverallActivity_ON_education <- mr(UVdat_DeviceOverallActivity_ON_education)


scatter_plot_DeviceOverallActivity_ON_education <- mr_scatter_plot(result_DeviceOverallActivity_ON_education, UVdat_DeviceOverallActivity_ON_education)
scatter_plot_DeviceOverallActivity_ON_education[[1]]
ggsave(scatter_plot_DeviceOverallActivity_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DeviceOverallActivity_ON_education.pdf", width=7, height=7)


pleiotropy_DeviceOverallActivity_ON_education <- mr_pleiotropy_test(UVdat_DeviceOverallActivity_ON_education)
str(pleiotropy_DeviceOverallActivity_ON_education)
pleiotropy_DeviceOverallActivity_ON_education <- pleiotropy_DeviceOverallActivity_ON_education[,-c(1,2)]


heterogeneity_DeviceOverallActivity_ON_education <- mr_heterogeneity(UVdat_DeviceOverallActivity_ON_education)
str(heterogeneity_DeviceOverallActivity_ON_education)
heterogeneity_DeviceOverallActivity_ON_education <- heterogeneity_DeviceOverallActivity_ON_education[,-c(1,2)]


singleSNP_DeviceOverallActivity_ON_education <- mr_singlesnp(UVdat_DeviceOverallActivity_ON_education)
funnel_plot_DeviceOverallActivity_ON_education <- mr_funnel_plot(singleSNP_DeviceOverallActivity_ON_education)
funnel_plot_DeviceOverallActivity_ON_education[[1]]
ggsave(funnel_plot_DeviceOverallActivity_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DeviceOverallActivity_ON_education.pdf", width=7, height=7)


loo_DeviceOverallActivity_ON_education <- mr_leaveoneout(UVdat_DeviceOverallActivity_ON_education)
loo_plot_DeviceOverallActivity_ON_education <- mr_leaveoneout_plot(loo_DeviceOverallActivity_ON_education)
loo_plot_DeviceOverallActivity_ON_education[[1]]
ggsave(loo_plot_DeviceOverallActivity_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_DeviceOverallActivity_ON_education.pdf", width=7, height=7)


singleSNP_plot_DeviceOverallActivity_ON_education <- mr_forest_plot(singleSNP_DeviceOverallActivity_ON_education)
singleSNP_plot_DeviceOverallActivity_ON_education[[1]]
ggsave(singleSNP_plot_DeviceOverallActivity_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DeviceOverallActivity_ON_education.pdf", width=7, height=7)

PRESSO_DeviceOverallActivity_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DeviceOverallActivity_ON_education, NbDistribution = 10000,  
                                          SignifThreshold = 0.05)

result_DeviceOverallActivity_ON_education <- generate_odds_ratios(result_DeviceOverallActivity_ON_education)
result_DeviceOverallActivity_ON_education <- result_DeviceOverallActivity_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_PA ####
IV_PA <-format_data(GWAS_PA,
                    type = "exposure", 
                    snp_col = "SNP", 
                    beta_col = "BETA",
                    se_col = "SE",
                    effect_allele_col = "effect allele",
                    other_allele_col = "other allele",
                    eaf_col = "Freq",
                    pval_col = "P",
                    samplesize_col = "N")

str(IV_PA)
IV_PA$exposure <- "physical activity"

IV_PA <- IV_PA[IV_PA$pval.exposure < 5e-8,]
IV_PA <- clump_data(IV_PA)

IV_PA$r2 <- (2 * (IV_PA$beta.exposure^2) * IV_PA$eaf.exposure * (1 - IV_PA$eaf.exposure)) /
  (2 * (IV_PA$beta.exposure^2) * IV_PA$eaf.exposure * (1 - IV_PA$eaf.exposure) +
     2 * IV_PA$samplesize.exposure * IV_PA$eaf.exposure * (1 - IV_PA$eaf.exposure) * (IV_PA$se.exposure^2))

IV_PA$F <- (IV_PA$r2 * (IV_PA$samplesize.exposure - 2))/(1 - IV_PA$r2)
range(IV_PA$F) # [1] 30.11588 81.42880
IV_PA_meanF <- mean(IV_PA$F) # [1] 37.19831

write.xlsx(IV_PA, file = "Output/UVMR_IV/IV_PA.xlsx")

out_PA_ON_education <- extract_outcome_data(
  snps = IV_PA$SNP,
  outcomes = 'ieu-a-1239')
str(out_PA_ON_education)
out_PA_ON_education$outcome <- "educational attainment"

UVdat_PA_ON_education <- harmonise_data(
  exposure_dat =  IV_PA, 
  outcome_dat = out_PA_ON_education
)
range(UVdat_PA_ON_education$pval.outcome)
range(UVdat_PA_ON_education$F)
str(UVdat_PA_ON_education)
UVdat_PA_ON_education <- UVdat_PA_ON_education[UVdat_PA_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_PA_ON_education,file = "Output/UVMR_IV/UVdat_PA_ON_education.xlsx")

result_PA_ON_education <- mr(UVdat_PA_ON_education)


scatter_plot_PA_ON_education <- mr_scatter_plot(result_PA_ON_education, UVdat_PA_ON_education)
scatter_plot_PA_ON_education[[1]]
ggsave(scatter_plot_PA_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_PA_ON_education.pdf", width=7, height=7)


pleiotropy_PA_ON_education <- mr_pleiotropy_test(UVdat_PA_ON_education)
str(pleiotropy_PA_ON_education)
pleiotropy_PA_ON_education <- pleiotropy_PA_ON_education[,-c(1,2)]


heterogeneity_PA_ON_education <- mr_heterogeneity(UVdat_PA_ON_education)
str(heterogeneity_PA_ON_education)
heterogeneity_PA_ON_education <- heterogeneity_PA_ON_education[,-c(1,2)]


singleSNP_PA_ON_education <- mr_singlesnp(UVdat_PA_ON_education)
funnel_plot_PA_ON_education <- mr_funnel_plot(singleSNP_PA_ON_education)
funnel_plot_PA_ON_education[[1]]
ggsave(funnel_plot_PA_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_PA_ON_education.pdf", width=7, height=7)


loo_PA_ON_education <- mr_leaveoneout(UVdat_PA_ON_education)
loo_plot_PA_ON_education <- mr_leaveoneout_plot(loo_PA_ON_education)
loo_plot_PA_ON_education[[1]]
ggsave(loo_plot_PA_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_PA_ON_education.pdf", width=7, height=7)


singleSNP_plot_PA_ON_education <- mr_forest_plot(singleSNP_PA_ON_education)
singleSNP_plot_PA_ON_education[[1]]
ggsave(singleSNP_plot_PA_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_PA_ON_education.pdf", width=7, height=7)

PRESSO_PA_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_PA_ON_education, NbDistribution = 10000,  
                                    SignifThreshold = 0.05)

result_PA_ON_education <- generate_odds_ratios(result_PA_ON_education)
result_PA_ON_education <- result_PA_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceModerate ####
IV_DeviceModerate <-format_data(GWAS_DeviceModerate,
                                       type = "exposure", 
                                       snp_col = "SNP", 
                                       beta_col = "BETA",
                                       se_col = "SE",
                                       effect_allele_col = "ALLELE1",
                                       other_allele_col = "ALLELE0",
                                       eaf_col = "A1FREQ",
                                       pval_col = "P_BOLT_LMM_INF",
                                       samplesize_col = "N")

str(IV_DeviceModerate)
IV_DeviceModerate$exposure <- "DeviceModerate"

IV_DeviceModerate <- IV_DeviceModerate[IV_DeviceModerate$pval.exposure < 5e-8,]
IV_DeviceModerate <- clump_data(IV_DeviceModerate)

IV_DeviceModerate$r2 <- (2 * (IV_DeviceModerate$beta.exposure^2) * IV_DeviceModerate$eaf.exposure * (1 - IV_DeviceModerate$eaf.exposure)) /
  (2 * (IV_DeviceModerate$beta.exposure^2) * IV_DeviceModerate$eaf.exposure * (1 - IV_DeviceModerate$eaf.exposure) +
     2 * IV_DeviceModerate$samplesize.exposure * IV_DeviceModerate$eaf.exposure * (1 - IV_DeviceModerate$eaf.exposure) * (IV_DeviceModerate$se.exposure^2))

IV_DeviceModerate$F <- (IV_DeviceModerate$r2 * (IV_DeviceModerate$samplesize.exposure - 2))/(1 - IV_DeviceModerate$r2)
range(IV_DeviceModerate$F) #[1] 32.09066
IV_DeviceModerate_meanF <- mean(IV_DeviceModerate$F) #[1] 32.09066

write.xlsx(IV_DeviceModerate, file = "Output/UVMR_IV/IV_DeviceModerate.xlsx")

out_DeviceModerate_ON_education <- extract_outcome_data(
  snps = IV_DeviceModerate$SNP,
  outcomes = 'ieu-a-1239')
str(out_DeviceModerate_ON_education)
out_DeviceModerate_ON_education$outcome <- "educational attainment"

UVdat_DeviceModerate_ON_education <- harmonise_data(
  exposure_dat =  IV_DeviceModerate, 
  outcome_dat = out_DeviceModerate_ON_education
)
str(UVdat_DeviceModerate_ON_education)
range(UVdat_DeviceModerate_ON_education$F)
range(UVdat_DeviceModerate_ON_education$pval.outcome)
write.xlsx(UVdat_DeviceModerate_ON_education,file = "Output/UVMR_IV/UVdat_DeviceModerate_ON_education.xlsx")

result_DeviceModerate_ON_education <- mr(UVdat_DeviceModerate_ON_education)


scatter_plot_DeviceModerate_ON_education <- mr_scatter_plot(result_DeviceModerate_ON_education, UVdat_DeviceModerate_ON_education)
scatter_plot_DeviceModerate_ON_education[[1]]
ggsave(scatter_plot_DeviceModerate_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DeviceModerate_ON_education.pdf", width=7, height=7)

result_DeviceModerate_ON_education <- generate_odds_ratios(result_DeviceModerate_ON_education)
result_DeviceModerate_ON_education <- result_DeviceModerate_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceSedentary ####
IV_DeviceSedentary <-format_data(GWAS_DeviceSedentary,
                                       type = "exposure", 
                                       snp_col = "SNP", 
                                       beta_col = "BETA",
                                       se_col = "SE",
                                       effect_allele_col = "ALLELE1",
                                       other_allele_col = "ALLELE0",
                                       eaf_col = "A1FREQ",
                                       pval_col = "P_BOLT_LMM_INF",
                                       samplesize_col = "N")

str(IV_DeviceSedentary)
IV_DeviceSedentary$exposure <- "DeviceSedentary"

IV_DeviceSedentary <- IV_DeviceSedentary[IV_DeviceSedentary$pval.exposure < 5e-8,]
IV_DeviceSedentary <- clump_data(IV_DeviceSedentary)

IV_DeviceSedentary$r2 <- (2 * (IV_DeviceSedentary$beta.exposure^2) * IV_DeviceSedentary$eaf.exposure * (1 - IV_DeviceSedentary$eaf.exposure)) /
  (2 * (IV_DeviceSedentary$beta.exposure^2) * IV_DeviceSedentary$eaf.exposure * (1 - IV_DeviceSedentary$eaf.exposure) +
     2 * IV_DeviceSedentary$samplesize.exposure * IV_DeviceSedentary$eaf.exposure * (1 - IV_DeviceSedentary$eaf.exposure) * (IV_DeviceSedentary$se.exposure^2))

IV_DeviceSedentary$F <- (IV_DeviceSedentary$r2 * (IV_DeviceSedentary$samplesize.exposure - 2))/(1 - IV_DeviceSedentary$r2)
range(IV_DeviceSedentary$F) #[1] 30.18161 35.44605
IV_DeviceSedentary_meanF <- mean(IV_DeviceSedentary$F) #[1] 33.54843

write.xlsx(IV_DeviceSedentary, file = "Output/UVMR_IV/IV_DeviceSedentary.xlsx")

out_DeviceSedentary_ON_education <- extract_outcome_data(
  snps = IV_DeviceSedentary$SNP,
  outcomes = 'ieu-a-1239')
str(out_DeviceSedentary_ON_education)
out_DeviceSedentary_ON_education$outcome <- "educational attainment"

UVdat_DeviceSedentary_ON_education <- harmonise_data(
  exposure_dat =  IV_DeviceSedentary, 
  outcome_dat = out_DeviceSedentary_ON_education
)
str(UVdat_DeviceSedentary_ON_education)
range(UVdat_DeviceSedentary_ON_education$F)
range(UVdat_DeviceSedentary_ON_education$pval.outcome)
UVdat_DeviceSedentary_ON_education <- UVdat_DeviceSedentary_ON_education[UVdat_DeviceSedentary_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_DeviceSedentary_ON_education,file = "Output/UVMR_IV/UVdat_DeviceSedentary_ON_education.xlsx")

result_DeviceSedentary_ON_education <- mr(UVdat_DeviceSedentary_ON_education)


scatter_plot_DeviceSedentary_ON_education <- mr_scatter_plot(result_DeviceSedentary_ON_education, UVdat_DeviceSedentary_ON_education)
scatter_plot_DeviceSedentary_ON_education[[1]]
ggsave(scatter_plot_DeviceSedentary_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DeviceSedentary_ON_education.pdf", width=7, height=7)


pleiotropy_DeviceSedentary_ON_education <- mr_pleiotropy_test(UVdat_DeviceSedentary_ON_education)
str(pleiotropy_DeviceSedentary_ON_education)
pleiotropy_DeviceSedentary_ON_education <- pleiotropy_DeviceSedentary_ON_education[,-c(1,2)]


heterogeneity_DeviceSedentary_ON_education <- mr_heterogeneity(UVdat_DeviceSedentary_ON_education)
str(heterogeneity_DeviceSedentary_ON_education)
heterogeneity_DeviceSedentary_ON_education <- heterogeneity_DeviceSedentary_ON_education[,-c(1,2)]


singleSNP_DeviceSedentary_ON_education <- mr_singlesnp(UVdat_DeviceSedentary_ON_education)
funnel_plot_DeviceSedentary_ON_education <- mr_funnel_plot(singleSNP_DeviceSedentary_ON_education)
funnel_plot_DeviceSedentary_ON_education[[1]]
ggsave(funnel_plot_DeviceSedentary_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DeviceSedentary_ON_education.pdf", width=7, height=7)


loo_DeviceSedentary_ON_education <- mr_leaveoneout(UVdat_DeviceSedentary_ON_education)
loo_plot_DeviceSedentary_ON_education <- mr_leaveoneout_plot(loo_DeviceSedentary_ON_education)
loo_plot_DeviceSedentary_ON_education[[1]]
ggsave(loo_plot_DeviceSedentary_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_DeviceSedentary_ON_education.pdf", width=7, height=7)


singleSNP_plot_DeviceSedentary_ON_education <- mr_forest_plot(singleSNP_DeviceSedentary_ON_education)
singleSNP_plot_DeviceSedentary_ON_education[[1]]
ggsave(singleSNP_plot_DeviceSedentary_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DeviceSedentary_ON_education.pdf", width=7, height=7)

PRESSO_DeviceSedentary_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DeviceSedentary_ON_education, NbDistribution = 10000,  
                                                       SignifThreshold = 0.05)

result_DeviceSedentary_ON_education <- generate_odds_ratios(result_DeviceSedentary_ON_education)
result_DeviceSedentary_ON_education <- result_DeviceSedentary_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_LST ####
IV_LST <-format_data(GWAS_LST,
                    type = "exposure", 
                    snp_col = "SNP", 
                    beta_col = "BETA",
                    se_col = "SE",
                    effect_allele_col = "effect allele",
                    other_allele_col = "other allele",
                    eaf_col = "Freq",
                    pval_col = "P",
                    samplesize_col = "N")

str(IV_LST)
IV_LST$exposure <- "LST"

IV_LST <- IV_LST[IV_LST$pval.exposure < 5e-8,]
IV_LST <- clump_data(IV_LST)

IV_LST$r2 <- (2 * (IV_LST$beta.exposure^2) * IV_LST$eaf.exposure * (1 - IV_LST$eaf.exposure)) /
  (2 * (IV_LST$beta.exposure^2) * IV_LST$eaf.exposure * (1 - IV_LST$eaf.exposure) +
     2 * IV_LST$samplesize.exposure * IV_LST$eaf.exposure * (1 - IV_LST$eaf.exposure) * (IV_LST$se.exposure^2))

IV_LST$F <- (IV_LST$r2 * (IV_LST$samplesize.exposure - 2))/(1 - IV_LST$r2)
range(IV_LST$F) # [1] 29.34017 97.16289
IV_LST_meanF <- mean(IV_LST$F) # [1] 39.27694

write.xlsx(IV_LST, file = "Output/UVMR_IV/IV_LST.xlsx")

out_LST_ON_education <- extract_outcome_data(
  snps = IV_LST$SNP,
  outcomes = 'ieu-a-1239')
str(out_LST_ON_education)
out_LST_ON_education$outcome <- "educational attainment"

UVdat_LST_ON_education <- harmonise_data(
  exposure_dat =  IV_LST, 
  outcome_dat = out_LST_ON_education
)
str(UVdat_LST_ON_education)
range(UVdat_LST_ON_education$pval.outcome)
range(UVdat_LST_ON_education$F)
UVdat_LST_ON_education <- UVdat_LST_ON_education[UVdat_LST_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_LST_ON_education,file = "Output/UVMR_IV/UVdat_LST_ON_education.xlsx")

result_LST_ON_education <- mr(UVdat_LST_ON_education)


scatter_plot_LST_ON_education <- mr_scatter_plot(result_LST_ON_education, UVdat_LST_ON_education)
scatter_plot_LST_ON_education[[1]]
ggsave(scatter_plot_LST_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_LST_ON_education.pdf", width=7, height=7)


pleiotropy_LST_ON_education <- mr_pleiotropy_test(UVdat_LST_ON_education)
str(pleiotropy_LST_ON_education)
pleiotropy_LST_ON_education <- pleiotropy_LST_ON_education[,-c(1,2)]


heterogeneity_LST_ON_education <- mr_heterogeneity(UVdat_LST_ON_education)
str(heterogeneity_LST_ON_education)
heterogeneity_LST_ON_education <- heterogeneity_LST_ON_education[,-c(1,2)]


singleSNP_LST_ON_education <- mr_singlesnp(UVdat_LST_ON_education)
funnel_plot_LST_ON_education <- mr_funnel_plot(singleSNP_LST_ON_education)
funnel_plot_LST_ON_education[[1]]
ggsave(funnel_plot_LST_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_LST_ON_education.pdf", width=7, height=7)


loo_LST_ON_education <- mr_leaveoneout(UVdat_LST_ON_education)
loo_plot_LST_ON_education <- mr_leaveoneout_plot(loo_LST_ON_education)
loo_plot_LST_ON_education[[1]]
ggsave(loo_plot_LST_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_LST_ON_education.pdf", width=7, height=21)


singleSNP_plot_LST_ON_education <- mr_forest_plot(singleSNP_LST_ON_education)
singleSNP_plot_LST_ON_education[[1]]
ggsave(singleSNP_plot_LST_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_LST_ON_education.pdf", width=7, height=7)

PRESSO_LST_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_LST_ON_education, NbDistribution = 10000,  
                                    SignifThreshold = 0.05)

result_LST_ON_education <- generate_odds_ratios(result_LST_ON_education)
result_LST_ON_education <- result_LST_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_television ####
IV_television <-format_data(GWAS_television,
                            type = "exposure", 
                            snp_col = "SNP", 
                            beta_col = "BETA",
                            se_col = "SE",
                            effect_allele_col = "effect allele",
                            other_allele_col = "other allele",
                            eaf_col = "Freq",
                            pval_col = "P",
                            samplesize_col = "N")

str(IV_television)
IV_television$exposure <- "television"

IV_television <- IV_television[IV_television$pval.exposure < 5e-8,]
IV_television <- clump_data(IV_television)

IV_television$r2 <- (2 * (IV_television$beta.exposure^2) * IV_television$eaf.exposure * (1 - IV_television$eaf.exposure)) /
  (2 * (IV_television$beta.exposure^2) * IV_television$eaf.exposure * (1 - IV_television$eaf.exposure) +
     2 * IV_television$samplesize.exposure * IV_television$eaf.exposure * (1 - IV_television$eaf.exposure) * (IV_television$se.exposure^2))

IV_television$F <- (IV_television$r2 * (IV_television$samplesize.exposure - 2))/(1 - IV_television$r2)
range(IV_television$F) # [1]  30.02792 144.18876
IV_television_meanF <- mean(IV_television$F) # [1] 42.20146

write.xlsx(IV_television, file = "Output/UVMR_IV/IV_television.xlsx")

out_television_ON_education <- extract_outcome_data(
  snps = IV_television$SNP,
  outcomes = 'ieu-a-1239')
str(out_television_ON_education)
out_television_ON_education$outcome <- "educational attainment"

UVdat_television_ON_education <- harmonise_data(
  exposure_dat =  IV_television, 
  outcome_dat = out_television_ON_education
)
range(UVdat_television_ON_education$pval.outcome)
str(UVdat_television_ON_education)
UVdat_television_ON_education <- UVdat_television_ON_education[UVdat_television_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_television_ON_education,file = "Output/UVMR_IV/UVdat_television_ON_education.xlsx")

result_television_ON_education <- mr(UVdat_television_ON_education)


scatter_plot_television_ON_education <- mr_scatter_plot(result_television_ON_education, UVdat_television_ON_education)
scatter_plot_television_ON_education[[1]]
ggsave(scatter_plot_television_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_television_ON_education.pdf", width=7, height=7)


pleiotropy_television_ON_education <- mr_pleiotropy_test(UVdat_television_ON_education)
str(pleiotropy_television_ON_education)
pleiotropy_television_ON_education <- pleiotropy_television_ON_education[,-c(1,2)]


heterogeneity_television_ON_education <- mr_heterogeneity(UVdat_television_ON_education)
str(heterogeneity_television_ON_education)
heterogeneity_television_ON_education <- heterogeneity_television_ON_education[,-c(1,2)]


singleSNP_television_ON_education <- mr_singlesnp(UVdat_television_ON_education)
funnel_plot_television_ON_education <- mr_funnel_plot(singleSNP_television_ON_education)
funnel_plot_television_ON_education[[1]]
ggsave(funnel_plot_television_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_television_ON_education.pdf", width=7, height=7)


loo_television_ON_education <- mr_leaveoneout(UVdat_television_ON_education)
loo_plot_television_ON_education <- mr_leaveoneout_plot(loo_television_ON_education)
loo_plot_television_ON_education[[1]]
ggsave(loo_plot_television_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_television_ON_education.pdf", width=7, height=7)


singleSNP_plot_television_ON_education <- mr_forest_plot(singleSNP_television_ON_education)
singleSNP_plot_television_ON_education[[1]]
ggsave(singleSNP_plot_television_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_television_ON_education.pdf", width=7, height=7)

PRESSO_television_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_television_ON_education, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_television_ON_education <- generate_odds_ratios(result_television_ON_education)
result_television_ON_education <- result_television_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_computer ####
IV_computer <- subset(GWAS_computer,P_BOLT_LMM_INF<5e-8)
IV_computer <- clump_data(IV_computer)

IV_computer <-format_data(IV_computer,
                            type = "exposure", 
                            snp_col = "SNP", 
                            beta_col = "BETA",
                            se_col = "SE",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            eaf_col = "A1FREQ",
                            pval_col = "P_BOLT_LMM_INF",
                            samplesize_col = "N")

str(IV_computer)
IV_computer$exposure <- "computer"

IV_computer$r2 <- (2 * (IV_computer$beta.exposure^2) * IV_computer$eaf.exposure * (1 - IV_computer$eaf.exposure)) /
  (2 * (IV_computer$beta.exposure^2) * IV_computer$eaf.exposure * (1 - IV_computer$eaf.exposure) +
     2 * IV_computer$samplesize.exposure * IV_computer$eaf.exposure * (1 - IV_computer$eaf.exposure) * (IV_computer$se.exposure^2))

IV_computer$F <- (IV_computer$r2 * (IV_computer$samplesize.exposure - 2))/(1 - IV_computer$r2)
range(IV_computer$F) # [1] 29.89740 49.12252
IV_computer_meanF <- mean(IV_computer$F) # [1] 33.47493

write.xlsx(IV_computer, file = "Output/UVMR_IV/IV_computer.xlsx")

out_computer_ON_education <- extract_outcome_data(
  snps = IV_computer$SNP,
  outcomes = 'ieu-a-1239')
str(out_computer_ON_education)
out_computer_ON_education$outcome <- "educational attainment"

UVdat_computer_ON_education <- harmonise_data(
  exposure_dat =  IV_computer, 
  outcome_dat = out_computer_ON_education
)
str(UVdat_computer_ON_education)
range(UVdat_computer_ON_education$pval.outcome)
UVdat_computer_ON_education <- UVdat_computer_ON_education[UVdat_computer_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_computer_ON_education,file = "Output/UVMR_IV/UVdat_computer_ON_education.xlsx")

result_computer_ON_education <- mr(UVdat_computer_ON_education)


scatter_plot_computer_ON_education <- mr_scatter_plot(result_computer_ON_education, UVdat_computer_ON_education)
scatter_plot_computer_ON_education[[1]]
ggsave(scatter_plot_computer_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_computer_ON_education.pdf", width=7, height=7)


pleiotropy_computer_ON_education <- mr_pleiotropy_test(UVdat_computer_ON_education)
str(pleiotropy_computer_ON_education)
pleiotropy_computer_ON_education <- pleiotropy_computer_ON_education[,-c(1,2)]


heterogeneity_computer_ON_education <- mr_heterogeneity(UVdat_computer_ON_education)
str(heterogeneity_computer_ON_education)
heterogeneity_computer_ON_education <- heterogeneity_computer_ON_education[,-c(1,2)]


singleSNP_computer_ON_education <- mr_singlesnp(UVdat_computer_ON_education)
funnel_plot_computer_ON_education <- mr_funnel_plot(singleSNP_computer_ON_education)
funnel_plot_computer_ON_education[[1]]
ggsave(funnel_plot_computer_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_computer_ON_education.pdf", width=7, height=7)


loo_computer_ON_education <- mr_leaveoneout(UVdat_computer_ON_education)
loo_plot_computer_ON_education <- mr_leaveoneout_plot(loo_computer_ON_education)
loo_plot_computer_ON_education[[1]]
ggsave(loo_plot_computer_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_computer_ON_education.pdf", width=7, height=7)


singleSNP_plot_computer_ON_education <- mr_forest_plot(singleSNP_computer_ON_education)
singleSNP_plot_computer_ON_education[[1]]
ggsave(singleSNP_plot_computer_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_computer_ON_education.pdf", width=7, height=7)

PRESSO_computer_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_computer_ON_education, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_computer_ON_education <- generate_odds_ratios(result_computer_ON_education)
result_computer_ON_education <- result_computer_ON_education[,-c(1,2)]


#### ________M1.2 Physical activity_driving ####
IV_driving <-format_data(GWAS_driving,
                            type = "exposure", 
                            snp_col = "SNP", 
                            beta_col = "BETA",
                            se_col = "SE",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            eaf_col = "A1FREQ",
                            pval_col = "P_BOLT_LMM_INF",
                            samplesize_col = "N")

str(IV_driving)
IV_driving$exposure <- "driving"

IV_driving <- IV_driving[IV_driving$pval.exposure < 5e-8,]
IV_driving <- clump_data(IV_driving)

IV_driving$r2 <- (2 * (IV_driving$beta.exposure^2) * IV_driving$eaf.exposure * (1 - IV_driving$eaf.exposure)) /
  (2 * (IV_driving$beta.exposure^2) * IV_driving$eaf.exposure * (1 - IV_driving$eaf.exposure) +
     2 * IV_driving$samplesize.exposure * IV_driving$eaf.exposure * (1 - IV_driving$eaf.exposure) * (IV_driving$se.exposure^2))

IV_driving$F <- (IV_driving$r2 * (IV_driving$samplesize.exposure - 2))/(1 - IV_driving$r2)
range(IV_driving$F) # [1] 29.91545 45.00637
IV_driving_meanF <- mean(IV_driving$F) # [1] 37.58354

write.xlsx(IV_driving, file = "Output/UVMR_IV/IV_driving.xlsx")

out_driving_ON_education <- extract_outcome_data(
  snps = IV_driving$SNP,
  outcomes = 'ieu-a-1239')
str(out_driving_ON_education)
out_driving_ON_education$outcome <- "educational attainment"

UVdat_driving_ON_education <- harmonise_data(
  exposure_dat =  IV_driving, 
  outcome_dat = out_driving_ON_education
)
str(UVdat_driving_ON_education)
range(UVdat_driving_ON_education$pval.outcome)
UVdat_driving_ON_education <- UVdat_driving_ON_education[UVdat_driving_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_driving_ON_education,file = "Output/UVMR_IV/UVdat_driving_ON_education.xlsx")

result_driving_ON_education <- mr(UVdat_driving_ON_education)


scatter_plot_driving_ON_education <- mr_scatter_plot(result_driving_ON_education, UVdat_driving_ON_education)
scatter_plot_driving_ON_education[[1]]
ggsave(scatter_plot_driving_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_driving_ON_education.pdf", width=7, height=7)


pleiotropy_driving_ON_education <- mr_pleiotropy_test(UVdat_driving_ON_education)
str(pleiotropy_driving_ON_education)
pleiotropy_driving_ON_education <- pleiotropy_driving_ON_education[,-c(1,2)]


heterogeneity_driving_ON_education <- mr_heterogeneity(UVdat_driving_ON_education)
str(heterogeneity_driving_ON_education)
heterogeneity_driving_ON_education <- heterogeneity_driving_ON_education[,-c(1,2)]


singleSNP_driving_ON_education <- mr_singlesnp(UVdat_driving_ON_education)
funnel_plot_driving_ON_education <- mr_funnel_plot(singleSNP_driving_ON_education)
funnel_plot_driving_ON_education[[1]]
ggsave(funnel_plot_driving_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_driving_ON_education.pdf", width=7, height=7)


loo_driving_ON_education <- mr_leaveoneout(UVdat_driving_ON_education)
loo_plot_driving_ON_education <- mr_leaveoneout_plot(loo_driving_ON_education)
loo_plot_driving_ON_education[[1]]
ggsave(loo_plot_driving_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_driving_ON_education.pdf", width=7, height=7)


singleSNP_plot_driving_ON_education <- mr_forest_plot(singleSNP_driving_ON_education)
singleSNP_plot_driving_ON_education[[1]]
ggsave(singleSNP_plot_driving_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_driving_ON_education.pdf", width=7, height=7)

PRESSO_driving_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_driving_ON_education, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_driving_ON_education <- generate_odds_ratios(result_driving_ON_education)
result_driving_ON_education <- result_driving_ON_education[,-c(1,2)]


#### ________M1.3 Smoking_LifetimeSmoking ####
IV_LifetimeSmoking <-format_data(GWAS_LifetimeSmoking,
                         type = "exposure", 
                         snp_col = "SNP", 
                         beta_col = "BETA",
                         se_col = "SE",
                         effect_allele_col = "EFFECT_ALLELE",
                         other_allele_col = "OTHER_ALLELE",
                         eaf_col = "EAF",
                         pval_col = "P",
                         samplesize_col = "SS")

str(IV_LifetimeSmoking)
IV_LifetimeSmoking$exposure <- "Lifetime smoking index"

IV_LifetimeSmoking <- IV_LifetimeSmoking[IV_LifetimeSmoking$pval.exposure < 5e-8,]
IV_LifetimeSmoking <- clump_data(IV_LifetimeSmoking)

IV_LifetimeSmoking$r2 <- (2 * (IV_LifetimeSmoking$beta.exposure^2) * IV_LifetimeSmoking$eaf.exposure * (1 - IV_LifetimeSmoking$eaf.exposure)) /
  (2 * (IV_LifetimeSmoking$beta.exposure^2) * IV_LifetimeSmoking$eaf.exposure * (1 - IV_LifetimeSmoking$eaf.exposure) +
     2 * IV_LifetimeSmoking$samplesize.exposure * IV_LifetimeSmoking$eaf.exposure * (1 - IV_LifetimeSmoking$eaf.exposure) * (IV_LifetimeSmoking$se.exposure^2))

IV_LifetimeSmoking$F <- (IV_LifetimeSmoking$r2 * (IV_LifetimeSmoking$samplesize.exposure - 2))/(1 - IV_LifetimeSmoking$r2)
range(IV_LifetimeSmoking$F) # [1]  30.04769 172.80187
IV_LifetimeSmoking_meanF <- mean(IV_LifetimeSmoking$F) # [1] 44.08217

write.xlsx(IV_LifetimeSmoking, file = "Output/UVMR_IV/IV_LifetimeSmoking.xlsx")

out_LifetimeSmoking_ON_education <- extract_outcome_data(
  snps = IV_LifetimeSmoking$SNP,
  outcomes = 'ieu-a-1239')
str(out_LifetimeSmoking_ON_education)
out_LifetimeSmoking_ON_education$outcome <- "educational attainment"

UVdat_LifetimeSmoking_ON_education <- harmonise_data(
  exposure_dat =  IV_LifetimeSmoking, 
  outcome_dat = out_LifetimeSmoking_ON_education
)
str(UVdat_LifetimeSmoking_ON_education)
range(UVdat_LifetimeSmoking_ON_education$pval.outcome)
range(UVdat_LifetimeSmoking_ON_education$F)
UVdat_LifetimeSmoking_ON_education <- UVdat_LifetimeSmoking_ON_education[UVdat_LifetimeSmoking_ON_education$pval.outcome>5e-8,]

write.xlsx(UVdat_LifetimeSmoking_ON_education,file = "Output/UVMR_IV/UVdat_LifetimeSmoking_ON_education.xlsx")

result_LifetimeSmoking_ON_education <- mr(UVdat_LifetimeSmoking_ON_education)


scatter_plot_LifetimeSmoking_ON_education <- mr_scatter_plot(result_LifetimeSmoking_ON_education, UVdat_LifetimeSmoking_ON_education)
scatter_plot_LifetimeSmoking_ON_education[[1]]
ggsave(scatter_plot_LifetimeSmoking_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_LifetimeSmoking_ON_education.pdf", width=7, height=7)


pleiotropy_LifetimeSmoking_ON_education <- mr_pleiotropy_test(UVdat_LifetimeSmoking_ON_education)
str(pleiotropy_LifetimeSmoking_ON_education)
pleiotropy_LifetimeSmoking_ON_education <- pleiotropy_LifetimeSmoking_ON_education[,-c(1,2)]


heterogeneity_LifetimeSmoking_ON_education <- mr_heterogeneity(UVdat_LifetimeSmoking_ON_education)
str(heterogeneity_LifetimeSmoking_ON_education)
heterogeneity_LifetimeSmoking_ON_education <- heterogeneity_LifetimeSmoking_ON_education[,-c(1,2)]


singleSNP_LifetimeSmoking_ON_education <- mr_singlesnp(UVdat_LifetimeSmoking_ON_education)
funnel_plot_LifetimeSmoking_ON_education <- mr_funnel_plot(singleSNP_LifetimeSmoking_ON_education)
funnel_plot_LifetimeSmoking_ON_education[[1]]
ggsave(funnel_plot_LifetimeSmoking_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_LifetimeSmoking_ON_education.pdf", width=7, height=7)


loo_LifetimeSmoking_ON_education <- mr_leaveoneout(UVdat_LifetimeSmoking_ON_education)
loo_plot_LifetimeSmoking_ON_education <- mr_leaveoneout_plot(loo_LifetimeSmoking_ON_education)
loo_plot_LifetimeSmoking_ON_education[[1]]
ggsave(loo_plot_LifetimeSmoking_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_LifetimeSmoking_ON_education.pdf", width=7, height=21)


singleSNP_plot_LifetimeSmoking_ON_education <- mr_forest_plot(singleSNP_LifetimeSmoking_ON_education)
singleSNP_plot_LifetimeSmoking_ON_education[[1]]
ggsave(singleSNP_plot_LifetimeSmoking_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_LifetimeSmoking_ON_education.pdf", width=7, height=7)

PRESSO_LifetimeSmoking_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_LifetimeSmoking_ON_education, NbDistribution = 10000,  
                                         SignifThreshold = 0.05)

result_LifetimeSmoking_ON_education <- generate_odds_ratios(result_LifetimeSmoking_ON_education)
result_LifetimeSmoking_ON_education <- result_LifetimeSmoking_ON_education[,-c(1,2)]


#### ________M1.4 Sleep health_SleepDuration ####
IV_SleepDuration <-format_data(GWAS_SleepDuration,
                       type = "exposure", 
                       snp_col = "SNP", 
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "A1",
                       other_allele_col = "A2",
                       eaf_col = "MAF",
                       pval_col = "P",
                       samplesize_col = "N")

str(IV_SleepDuration)
IV_SleepDuration$exposure <- "SleepDuration"

IV_SleepDuration <- IV_SleepDuration[IV_SleepDuration$pval.exposure < 5e-8,]
IV_SleepDuration <- clump_data(IV_SleepDuration)

IV_SleepDuration$r2 <- (2 * (IV_SleepDuration$beta.exposure^2) * IV_SleepDuration$eaf.exposure * (1 - IV_SleepDuration$eaf.exposure)) /
  (2 * (IV_SleepDuration$beta.exposure^2) * IV_SleepDuration$eaf.exposure * (1 - IV_SleepDuration$eaf.exposure) +
     2 * IV_SleepDuration$samplesize.exposure * IV_SleepDuration$eaf.exposure * (1 - IV_SleepDuration$eaf.exposure) * (IV_SleepDuration$se.exposure^2))

IV_SleepDuration$F <- (IV_SleepDuration$r2 * (IV_SleepDuration$samplesize.exposure - 2))/(1 - IV_SleepDuration$r2)
range(IV_SleepDuration$F) #[1]  29.80872 190.17499
IV_SleepDuration_meanF <- mean(IV_SleepDuration$F) #[1] 39.78967

write.xlsx(IV_SleepDuration, file = "Output/UVMR_IV/IV_SleepDuration.xlsx")

out_SleepDuration_ON_education <- extract_outcome_data(
  snps = IV_SleepDuration$SNP,
  outcomes = 'ieu-a-1239')
str(out_SleepDuration_ON_education)
out_SleepDuration_ON_education$outcome <- "educational attainment"

UVdat_SleepDuration_ON_education <- harmonise_data(
  exposure_dat =  IV_SleepDuration, 
  outcome_dat = out_SleepDuration_ON_education
)
str(UVdat_SleepDuration_ON_education)
range(UVdat_SleepDuration_ON_education$F)
range(UVdat_SleepDuration_ON_education$pval.outcome)
UVdat_SleepDuration_ON_education <- UVdat_SleepDuration_ON_education[UVdat_SleepDuration_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_SleepDuration_ON_education,file = "Output/UVMR_IV/UVdat_SleepDuration_ON_education.xlsx")

result_SleepDuration_ON_education <- mr(UVdat_SleepDuration_ON_education)


scatter_plot_SleepDuration_ON_education <- mr_scatter_plot(result_SleepDuration_ON_education, UVdat_SleepDuration_ON_education)
scatter_plot_SleepDuration_ON_education[[1]]
ggsave(scatter_plot_SleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_SleepDuration_ON_education.pdf", width=7, height=7)


pleiotropy_SleepDuration_ON_education <- mr_pleiotropy_test(UVdat_SleepDuration_ON_education)
str(pleiotropy_SleepDuration_ON_education)
pleiotropy_SleepDuration_ON_education <- pleiotropy_SleepDuration_ON_education[,-c(1,2)]


heterogeneity_SleepDuration_ON_education <- mr_heterogeneity(UVdat_SleepDuration_ON_education)
str(heterogeneity_SleepDuration_ON_education)
heterogeneity_SleepDuration_ON_education <- heterogeneity_SleepDuration_ON_education[,-c(1,2)]


singleSNP_SleepDuration_ON_education <- mr_singlesnp(UVdat_SleepDuration_ON_education)
funnel_plot_SleepDuration_ON_education <- mr_funnel_plot(singleSNP_SleepDuration_ON_education)
funnel_plot_SleepDuration_ON_education[[1]]
ggsave(funnel_plot_SleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_SleepDuration_ON_education.pdf", width=7, height=7)


loo_SleepDuration_ON_education <- mr_leaveoneout(UVdat_SleepDuration_ON_education)
loo_plot_SleepDuration_ON_education <- mr_leaveoneout_plot(loo_SleepDuration_ON_education)
loo_plot_SleepDuration_ON_education[[1]]
ggsave(loo_plot_SleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_SleepDuration_ON_education.pdf", width=7, height=7)


singleSNP_plot_SleepDuration_ON_education <- mr_forest_plot(singleSNP_SleepDuration_ON_education)
singleSNP_plot_SleepDuration_ON_education[[1]]
ggsave(singleSNP_plot_SleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_SleepDuration_ON_education.pdf", width=7, height=7)

PRESSO_SleepDuration_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_SleepDuration_ON_education, NbDistribution = 10000,  
                                       SignifThreshold = 0.05)

result_SleepDuration_ON_education <- generate_odds_ratios(result_SleepDuration_ON_education)
result_SleepDuration_ON_education <- result_SleepDuration_ON_education[,-c(1,2)]


#### ________M1.4 Sleep health_DeviceSleepDuration ####
IV_DeviceSleepDuration <-format_data(GWAS_DeviceSleepDuration,
                               type = "exposure", 
                               snp_col = "SNP", 
                               beta_col = "BETA",
                               se_col = "SE",
                               effect_allele_col = "ALLELE1",
                               other_allele_col = "ALLELE0",
                               eaf_col = "A1FREQ",
                               pval_col = "P_BOLT_LMM_INF",
                               samplesize_col = "N")

str(IV_DeviceSleepDuration)
IV_DeviceSleepDuration$exposure <- "DeviceSleepDuration"

IV_DeviceSleepDuration <- IV_DeviceSleepDuration[IV_DeviceSleepDuration$pval.exposure < 5e-8,]
IV_DeviceSleepDuration <- clump_data(IV_DeviceSleepDuration)

IV_DeviceSleepDuration$r2 <- (2 * (IV_DeviceSleepDuration$beta.exposure^2) * IV_DeviceSleepDuration$eaf.exposure * (1 - IV_DeviceSleepDuration$eaf.exposure)) /
  (2 * (IV_DeviceSleepDuration$beta.exposure^2) * IV_DeviceSleepDuration$eaf.exposure * (1 - IV_DeviceSleepDuration$eaf.exposure) +
     2 * IV_DeviceSleepDuration$samplesize.exposure * IV_DeviceSleepDuration$eaf.exposure * (1 - IV_DeviceSleepDuration$eaf.exposure) * (IV_DeviceSleepDuration$se.exposure^2))

IV_DeviceSleepDuration$F <- (IV_DeviceSleepDuration$r2 * (IV_DeviceSleepDuration$samplesize.exposure - 2))/(1 - IV_DeviceSleepDuration$r2)
range(IV_DeviceSleepDuration$F) #[1] 30.21700 83.41078
IV_DeviceSleepDuration_meanF <- mean(IV_DeviceSleepDuration$F) #[1] 42.67557

write.xlsx(IV_DeviceSleepDuration, file = "Output/UVMR_IV/IV_DeviceSleepDuration.xlsx")

out_DeviceSleepDuration_ON_education <- extract_outcome_data(
  snps = IV_DeviceSleepDuration$SNP,
  outcomes = 'ieu-a-1239')
str(out_DeviceSleepDuration_ON_education)
out_DeviceSleepDuration_ON_education$outcome <- "educational attainment"

UVdat_DeviceSleepDuration_ON_education <- harmonise_data(
  exposure_dat =  IV_DeviceSleepDuration, 
  outcome_dat = out_DeviceSleepDuration_ON_education
)
str(UVdat_DeviceSleepDuration_ON_education)
range(UVdat_DeviceSleepDuration_ON_education$F)
range(UVdat_DeviceSleepDuration_ON_education$pval.outcome)
UVdat_DeviceSleepDuration_ON_education <- UVdat_DeviceSleepDuration_ON_education[UVdat_DeviceSleepDuration_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_DeviceSleepDuration_ON_education,file = "Output/UVMR_IV/UVdat_DeviceSleepDuration_ON_education.xlsx")

result_DeviceSleepDuration_ON_education <- mr(UVdat_DeviceSleepDuration_ON_education)


scatter_plot_DeviceSleepDuration_ON_education <- mr_scatter_plot(result_DeviceSleepDuration_ON_education, UVdat_DeviceSleepDuration_ON_education)
scatter_plot_DeviceSleepDuration_ON_education[[1]]
ggsave(scatter_plot_DeviceSleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DeviceSleepDuration_ON_education.pdf", width=7, height=7)


pleiotropy_DeviceSleepDuration_ON_education <- mr_pleiotropy_test(UVdat_DeviceSleepDuration_ON_education)
str(pleiotropy_DeviceSleepDuration_ON_education)
pleiotropy_DeviceSleepDuration_ON_education <- pleiotropy_DeviceSleepDuration_ON_education[,-c(1,2)]


heterogeneity_DeviceSleepDuration_ON_education <- mr_heterogeneity(UVdat_DeviceSleepDuration_ON_education)
str(heterogeneity_DeviceSleepDuration_ON_education)
heterogeneity_DeviceSleepDuration_ON_education <- heterogeneity_DeviceSleepDuration_ON_education[,-c(1,2)]


singleSNP_DeviceSleepDuration_ON_education <- mr_singlesnp(UVdat_DeviceSleepDuration_ON_education)
funnel_plot_DeviceSleepDuration_ON_education <- mr_funnel_plot(singleSNP_DeviceSleepDuration_ON_education)
funnel_plot_DeviceSleepDuration_ON_education[[1]]
ggsave(funnel_plot_DeviceSleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DeviceSleepDuration_ON_education.pdf", width=7, height=7)


loo_DeviceSleepDuration_ON_education <- mr_leaveoneout(UVdat_DeviceSleepDuration_ON_education)
loo_plot_DeviceSleepDuration_ON_education <- mr_leaveoneout_plot(loo_DeviceSleepDuration_ON_education)
loo_plot_DeviceSleepDuration_ON_education[[1]]
ggsave(loo_plot_DeviceSleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_DeviceSleepDuration_ON_education.pdf", width=7, height=7)


singleSNP_plot_DeviceSleepDuration_ON_education <- mr_forest_plot(singleSNP_DeviceSleepDuration_ON_education)
singleSNP_plot_DeviceSleepDuration_ON_education[[1]]
ggsave(singleSNP_plot_DeviceSleepDuration_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DeviceSleepDuration_ON_education.pdf", width=7, height=7)

PRESSO_DeviceSleepDuration_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DeviceSleepDuration_ON_education, NbDistribution = 10000,  
                                               SignifThreshold = 0.05)

result_DeviceSleepDuration_ON_education <- generate_odds_ratios(result_DeviceSleepDuration_ON_education)
result_DeviceSleepDuration_ON_education <- result_DeviceSleepDuration_ON_education[,-c(1,2)]


#### ________M1.4 Sleep health_Insomnia ####
IV_Insomnia <-format_data(GWAS_Insomnia,
                               type = "exposure", 
                               snp_col = "SNP", 
                               beta_col = "BETA",
                               se_col = "SE",
                               effect_allele_col = "A1",
                               other_allele_col = "A2",
                               eaf_col = "MAF",
                               pval_col = "P",
                               samplesize_col = "N")

str(IV_Insomnia)
IV_Insomnia$exposure <- "Insomnia"

IV_Insomnia <- IV_Insomnia[IV_Insomnia$pval.exposure < 5e-8,]
IV_Insomnia <- clump_data(IV_Insomnia)

IV_Insomnia$r2 <- (2 * (IV_Insomnia$beta.exposure^2) * IV_Insomnia$eaf.exposure * (1 - IV_Insomnia$eaf.exposure)) /
  (2 * (IV_Insomnia$beta.exposure^2) * IV_Insomnia$eaf.exposure * (1 - IV_Insomnia$eaf.exposure) +
     2 * IV_Insomnia$samplesize.exposure * IV_Insomnia$eaf.exposure * (1 - IV_Insomnia$eaf.exposure) * (IV_Insomnia$se.exposure^2))

IV_Insomnia$F <- (IV_Insomnia$r2 * (IV_Insomnia$samplesize.exposure - 2))/(1 - IV_Insomnia$r2)
range(IV_Insomnia$F) #[1] 30.36327 94.69962
IV_Insomnia_meanF <- mean(IV_Insomnia$F) #[1] 41.19055

write.xlsx(IV_Insomnia, file = "Output/UVMR_IV/IV_Insomnia.xlsx")

out_Insomnia_ON_education <- extract_outcome_data(
  snps = IV_Insomnia$SNP,
  outcomes = 'ieu-a-1239')
str(out_Insomnia_ON_education)
out_Insomnia_ON_education$outcome <- "educational attainment"

UVdat_Insomnia_ON_education <- harmonise_data(
  exposure_dat =  IV_Insomnia, 
  outcome_dat = out_Insomnia_ON_education
)
str(UVdat_Insomnia_ON_education)
range(UVdat_Insomnia_ON_education$F)
range(UVdat_Insomnia_ON_education$pval.outcome)
UVdat_Insomnia_ON_education <- UVdat_Insomnia_ON_education[UVdat_Insomnia_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_Insomnia_ON_education,file = "Output/UVMR_IV/UVdat_Insomnia_ON_education.xlsx")

result_Insomnia_ON_education <- mr(UVdat_Insomnia_ON_education)


scatter_plot_Insomnia_ON_education <- mr_scatter_plot(result_Insomnia_ON_education, UVdat_Insomnia_ON_education)
scatter_plot_Insomnia_ON_education[[1]]
ggsave(scatter_plot_Insomnia_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_Insomnia_ON_education.pdf", width=7, height=7)


pleiotropy_Insomnia_ON_education <- mr_pleiotropy_test(UVdat_Insomnia_ON_education)
str(pleiotropy_Insomnia_ON_education)
pleiotropy_Insomnia_ON_education <- pleiotropy_Insomnia_ON_education[,-c(1,2)]


heterogeneity_Insomnia_ON_education <- mr_heterogeneity(UVdat_Insomnia_ON_education)
str(heterogeneity_Insomnia_ON_education)
heterogeneity_Insomnia_ON_education <- heterogeneity_Insomnia_ON_education[,-c(1,2)]


singleSNP_Insomnia_ON_education <- mr_singlesnp(UVdat_Insomnia_ON_education)
funnel_plot_Insomnia_ON_education <- mr_funnel_plot(singleSNP_Insomnia_ON_education)
funnel_plot_Insomnia_ON_education[[1]]
ggsave(funnel_plot_Insomnia_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_Insomnia_ON_education.pdf", width=7, height=7)


loo_Insomnia_ON_education <- mr_leaveoneout(UVdat_Insomnia_ON_education)
loo_plot_Insomnia_ON_education <- mr_leaveoneout_plot(loo_Insomnia_ON_education)
loo_plot_Insomnia_ON_education[[1]]
ggsave(loo_plot_Insomnia_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_Insomnia_ON_education.pdf", width=7, height=7)


singleSNP_plot_Insomnia_ON_education <- mr_forest_plot(singleSNP_Insomnia_ON_education)
singleSNP_plot_Insomnia_ON_education[[1]]
ggsave(singleSNP_plot_Insomnia_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_Insomnia_ON_education.pdf", width=7, height=7)

PRESSO_Insomnia_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_Insomnia_ON_education, NbDistribution = 10000,  
                                               SignifThreshold = 0.05)

result_Insomnia_ON_education <- generate_odds_ratios(result_Insomnia_ON_education)
result_Insomnia_ON_education <- result_Insomnia_ON_education[,-c(1,2)]


#### ________M2.1 Well-being spectrum(WBS) ####
IV_WBS <- GWAS_WBS[GWAS_WBS$PVAL < 5e-8,]
str(IV_WBS)

IV_WBS <-format_data(IV_WBS,
                     type = "exposure", 
                     snp_col = "SNP", 
                     beta_col = "BETA",
                     se_col = "SE",
                     effect_allele_col = "A1",
                     other_allele_col = "A2",
                     eaf_col = "EAF",
                     pval_col = "PVAL",
                     samplesize_col = "N")

str(IV_WBS)
IV_WBS$exposure <- "WBS"

IV_WBS <- clump_data(IV_WBS)

IV_WBS$r2 <- (2 * (IV_WBS$beta.exposure^2) * IV_WBS$eaf.exposure * (1 - IV_WBS$eaf.exposure)) /
  (2 * (IV_WBS$beta.exposure^2) * IV_WBS$eaf.exposure * (1 - IV_WBS$eaf.exposure) +
     2 * IV_WBS$samplesize.exposure * IV_WBS$eaf.exposure * (1 - IV_WBS$eaf.exposure) * (IV_WBS$se.exposure^2))

IV_WBS$F <- (IV_WBS$r2 * (IV_WBS$samplesize.exposure - 2))/(1 - IV_WBS$r2)
range(IV_WBS$F) #[1]  29.7274 109.2183
IV_WBS_meanF <- mean(IV_WBS$F) #[1] 42.78211

write.xlsx(IV_WBS, file = "Output/UVMR_IV/IV_WBS.xlsx")

out_WBS_ON_education <- extract_outcome_data(
  snps = IV_WBS$SNP,
  outcomes = 'ieu-a-1239')
str(out_WBS_ON_education)
out_WBS_ON_education$outcome <- "educational attainment"

UVdat_WBS_ON_education <- harmonise_data(
  exposure_dat =  IV_WBS, 
  outcome_dat = out_WBS_ON_education
)
str(UVdat_WBS_ON_education)
range(UVdat_WBS_ON_education$F)

write.xlsx(UVdat_WBS_ON_education,file = "Output/UVMR_IV/UVdat_WBS_ON_education.xlsx")

result_WBS_ON_education <- mr(UVdat_WBS_ON_education)


scatter_plot_WBS_ON_education <- mr_scatter_plot(result_WBS_ON_education, UVdat_WBS_ON_education)
scatter_plot_WBS_ON_education[[1]]
ggsave(scatter_plot_WBS_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_WBS_ON_education.pdf", width=7, height=7)


pleiotropy_WBS_ON_education <- mr_pleiotropy_test(UVdat_WBS_ON_education)
str(pleiotropy_WBS_ON_education)
pleiotropy_WBS_ON_education <- pleiotropy_WBS_ON_education[,-c(1,2)]


heterogeneity_WBS_ON_education <- mr_heterogeneity(UVdat_WBS_ON_education)
str(heterogeneity_WBS_ON_education)
heterogeneity_WBS_ON_education <- heterogeneity_WBS_ON_education[,-c(1,2)]


singleSNP_WBS_ON_education <- mr_singlesnp(UVdat_WBS_ON_education)
funnel_plot_WBS_ON_education <- mr_funnel_plot(singleSNP_WBS_ON_education)
funnel_plot_WBS_ON_education[[1]]
ggsave(funnel_plot_WBS_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_WBS_ON_education.pdf", width=7, height=7)


loo_WBS_ON_education <- mr_leaveoneout(UVdat_WBS_ON_education)
loo_plot_WBS_ON_education <- mr_leaveoneout_plot(loo_WBS_ON_education)
loo_plot_WBS_ON_education[[1]]
ggsave(loo_plot_WBS_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_WBS_ON_education.pdf", width=7, height=21)


singleSNP_plot_WBS_ON_education <- mr_forest_plot(singleSNP_WBS_ON_education)
singleSNP_plot_WBS_ON_education[[1]]
ggsave(singleSNP_plot_WBS_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_WBS_ON_education.pdf", width=7, height=7)

PRESSO_WBS_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_WBS_ON_education, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_WBS_ON_education <- generate_odds_ratios(result_WBS_ON_education)
result_WBS_ON_education <- result_WBS_ON_education[,-c(1,2)]


#### ________M2.2 MentalProblems ####
IV_MentalProblems <-format_data(GWAS_MentalProblems,
                     type = "exposure", 
                     snp_col = "SNP", 
                     beta_col = "ES",
                     se_col = "SE",
                     effect_allele_col = "ALT",
                     other_allele_col = "REF",
                     eaf_col = "AF",
                     pval_col = "P",
                     samplesize_col = "SS")

str(IV_MentalProblems)
IV_MentalProblems$exposure <- "MentalProblems"

IV_MentalProblems <- IV_MentalProblems[IV_MentalProblems$pval.exposure < 5e-8,]
IV_MentalProblems <- clump_data(IV_MentalProblems)

IV_MentalProblems$r2 <- (2 * (IV_MentalProblems$beta.exposure^2) * IV_MentalProblems$eaf.exposure * (1 - IV_MentalProblems$eaf.exposure)) /
  (2 * (IV_MentalProblems$beta.exposure^2) * IV_MentalProblems$eaf.exposure * (1 - IV_MentalProblems$eaf.exposure) +
     2 * IV_MentalProblems$samplesize.exposure * IV_MentalProblems$eaf.exposure * (1 - IV_MentalProblems$eaf.exposure) * (IV_MentalProblems$se.exposure^2))

IV_MentalProblems$F <- (IV_MentalProblems$r2 * (IV_MentalProblems$samplesize.exposure - 2))/(1 - IV_MentalProblems$r2)
range(IV_MentalProblems$F) #[1] 30.85997 53.59228
IV_MentalProblems_meanF <- mean(IV_MentalProblems$F) #[1] 37.92488

write.xlsx(IV_MentalProblems, file = "Output/UVMR_IV/IV_MentalProblems.xlsx")

out_MentalProblems_ON_education <- extract_outcome_data(
  snps = IV_MentalProblems$SNP,
  outcomes = 'ieu-a-1239')
str(out_MentalProblems_ON_education)
out_MentalProblems_ON_education$outcome <- "educational attainment"

UVdat_MentalProblems_ON_education <- harmonise_data(
  exposure_dat =  IV_MentalProblems, 
  outcome_dat = out_MentalProblems_ON_education
)
str(UVdat_MentalProblems_ON_education)
range(UVdat_MentalProblems_ON_education$F)
range(UVdat_MentalProblems_ON_education$pval.outcome)
UVdat_MentalProblems_ON_education <- UVdat_MentalProblems_ON_education[UVdat_MentalProblems_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_MentalProblems_ON_education,file = "Output/UVMR_IV/UVdat_MentalProblems_ON_education.xlsx")

result_MentalProblems_ON_education <- mr(UVdat_MentalProblems_ON_education)


scatter_plot_MentalProblems_ON_education <- mr_scatter_plot(result_MentalProblems_ON_education, UVdat_MentalProblems_ON_education)
scatter_plot_MentalProblems_ON_education[[1]]
ggsave(scatter_plot_MentalProblems_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_MentalProblems_ON_education.pdf", width=7, height=7)


pleiotropy_MentalProblems_ON_education <- mr_pleiotropy_test(UVdat_MentalProblems_ON_education)
str(pleiotropy_MentalProblems_ON_education)
pleiotropy_MentalProblems_ON_education <- pleiotropy_MentalProblems_ON_education[,-c(1,2)]


heterogeneity_MentalProblems_ON_education <- mr_heterogeneity(UVdat_MentalProblems_ON_education)
str(heterogeneity_MentalProblems_ON_education)
heterogeneity_MentalProblems_ON_education <- heterogeneity_MentalProblems_ON_education[,-c(1,2)]


singleSNP_MentalProblems_ON_education <- mr_singlesnp(UVdat_MentalProblems_ON_education)
funnel_plot_MentalProblems_ON_education <- mr_funnel_plot(singleSNP_MentalProblems_ON_education)
funnel_plot_MentalProblems_ON_education[[1]]
ggsave(funnel_plot_MentalProblems_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_MentalProblems_ON_education.pdf", width=7, height=7)


loo_MentalProblems_ON_education <- mr_leaveoneout(UVdat_MentalProblems_ON_education)
loo_plot_MentalProblems_ON_education <- mr_leaveoneout_plot(loo_MentalProblems_ON_education)
loo_plot_MentalProblems_ON_education[[1]]
ggsave(loo_plot_MentalProblems_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_MentalProblems_ON_education.pdf", width=7, height=21)


singleSNP_plot_MentalProblems_ON_education <- mr_forest_plot(singleSNP_MentalProblems_ON_education)
singleSNP_plot_MentalProblems_ON_education[[1]]
ggsave(singleSNP_plot_MentalProblems_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_MentalProblems_ON_education.pdf", width=7, height=7)

PRESSO_MentalProblems_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_MentalProblems_ON_education, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_MentalProblems_ON_education <- generate_odds_ratios(result_MentalProblems_ON_education)
result_MentalProblems_ON_education <- result_MentalProblems_ON_education[,-c(1,2)]


#### ________M2.3 Depression ####
IV_depression <- read.xlsx("Input/IV_depression.xlsx")
str(IV_depression)

IV_depression <-format_data(IV_depression,
                            type = "exposure", 
                            snp_col = "Marker.Name", 
                            beta_col = "Log(Odds.Ratio)",
                            se_col = "Standard.error.of.the.Log(Odds.Ratio)",
                            effect_allele_col = "A1",
                            other_allele_col = "A2",
                            eaf_col = "Allele.Frequency.in.meta-analysis",
                            pval_col = "P-value")

str(IV_depression)
IV_depression$exposure <- "depression"
IV_depression$samplesize.exposure <- 807553

IV_depression$r2 <- (2 * (IV_depression$beta.exposure^2) * IV_depression$eaf.exposure * (1 - IV_depression$eaf.exposure)) /
  (2 * (IV_depression$beta.exposure^2) * IV_depression$eaf.exposure * (1 - IV_depression$eaf.exposure) +
     2 * IV_depression$samplesize.exposure * IV_depression$eaf.exposure * (1 - IV_depression$eaf.exposure) * (IV_depression$se.exposure^2))

IV_depression$F <- (IV_depression$r2 * (IV_depression$samplesize.exposure - 2))/(1 - IV_depression$r2)
range(IV_depression$F) #[1]  30.24993 107.35236
IV_depression_meanF <- mean(IV_depression$F) #[1] 42.72652

write.xlsx(IV_depression, file = "Output/UVMR_IV/IV_depression.xlsx")

out_depression_ON_education <- extract_outcome_data(
  snps = IV_depression$SNP,
  outcomes = 'ieu-a-1239')
str(out_depression_ON_education)
out_depression_ON_education$outcome <- "educational attainment"

UVdat_depression_ON_education <- harmonise_data(
  exposure_dat =  IV_depression, 
  outcome_dat = out_depression_ON_education
)
range(UVdat_depression_ON_education$pval.outcome)
range(UVdat_depression_ON_education$F)
str(UVdat_depression_ON_education)
UVdat_depression_ON_education <- UVdat_depression_ON_education[UVdat_depression_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_depression_ON_education,file = "Output/UVMR_IV/UVdat_depression_ON_education.xlsx")

result_depression_ON_education <- mr(UVdat_depression_ON_education)


scatter_plot_depression_ON_education <- mr_scatter_plot(result_depression_ON_education, UVdat_depression_ON_education)
scatter_plot_depression_ON_education[[1]]
ggsave(scatter_plot_depression_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_depression_ON_education.pdf", width=7, height=7)


pleiotropy_depression_ON_education <- mr_pleiotropy_test(UVdat_depression_ON_education)
str(pleiotropy_depression_ON_education)
pleiotropy_depression_ON_education <- pleiotropy_depression_ON_education[,-c(1,2)]


heterogeneity_depression_ON_education <- mr_heterogeneity(UVdat_depression_ON_education)
str(heterogeneity_depression_ON_education)
heterogeneity_depression_ON_education <- heterogeneity_depression_ON_education[,-c(1,2)]


singleSNP_depression_ON_education <- mr_singlesnp(UVdat_depression_ON_education)
funnel_plot_depression_ON_education <- mr_funnel_plot(singleSNP_depression_ON_education)
funnel_plot_depression_ON_education[[1]]
ggsave(funnel_plot_depression_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_depression_ON_education.pdf", width=7, height=7)


loo_depression_ON_education <- mr_leaveoneout(UVdat_depression_ON_education)
loo_plot_depression_ON_education <- mr_leaveoneout_plot(loo_depression_ON_education)
loo_plot_depression_ON_education[[1]]
ggsave(loo_plot_depression_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_depression_ON_education.pdf", width=7, height=21)


singleSNP_plot_depression_ON_education <- mr_forest_plot(singleSNP_depression_ON_education)
singleSNP_plot_depression_ON_education[[1]]
ggsave(singleSNP_plot_depression_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_depression_ON_education.pdf", width=7, height=21)

PRESSO_depression_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                            OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_depression_ON_education, NbDistribution = 10000,  
                                            SignifThreshold = 0.05)

result_depression_ON_education <- generate_odds_ratios(result_depression_ON_education)
result_depression_ON_education <- result_depression_ON_education[,-c(1,2)]


#### ________M3.1 BMI ####
IV_BMI <-format_data(GWAS_BMI,
                     type = "exposure", 
                     snp_col = "SNP", 
                     beta_col = "ES",
                     se_col = "SE",
                     effect_allele_col = "ALT",
                     other_allele_col = "REF",
                     eaf_col = "AF",
                     pval_col = "P",
                     samplesize_col = "SS")

str(IV_BMI)
IV_BMI$exposure <- "BMI"

IV_BMI <- IV_BMI[IV_BMI$pval.exposure < 5e-8,]
IV_BMI <- clump_data(IV_BMI)

IV_BMI$r2 <- (2 * (IV_BMI$beta.exposure^2) * IV_BMI$eaf.exposure * (1 - IV_BMI$eaf.exposure)) /
  (2 * (IV_BMI$beta.exposure^2) * IV_BMI$eaf.exposure * (1 - IV_BMI$eaf.exposure) +
     2 * IV_BMI$samplesize.exposure * IV_BMI$eaf.exposure * (1 - IV_BMI$eaf.exposure) * (IV_BMI$se.exposure^2))

IV_BMI$F <- (IV_BMI$r2 * (IV_BMI$samplesize.exposure - 2))/(1 - IV_BMI$r2)
range(IV_BMI$F) #[1]   28.62242 1360.30446
IV_BMI_meanF <- mean(IV_BMI$F)#[1] 72.47668

write.xlsx(IV_BMI, file = "Output/UVMR_IV/IV_BMI.xlsx")

out_BMI_ON_education <- extract_outcome_data(
  snps = IV_BMI$SNP,
  outcomes = 'ieu-a-1239')
str(out_BMI_ON_education)
out_BMI_ON_education$outcome <- "educational attainment"

UVdat_BMI_ON_education <- harmonise_data(
  exposure_dat =  IV_BMI, 
  outcome_dat = out_BMI_ON_education
)
range(UVdat_BMI_ON_education$pval.outcome)
str(UVdat_BMI_ON_education)
UVdat_BMI_ON_education <- UVdat_BMI_ON_education[UVdat_BMI_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_BMI_ON_education,file = "Output/UVMR_IV/UVdat_BMI_ON_education.xlsx")

result_BMI_ON_education <- mr(UVdat_BMI_ON_education)


scatter_plot_BMI_ON_education <- mr_scatter_plot(result_BMI_ON_education, UVdat_BMI_ON_education)
scatter_plot_BMI_ON_education[[1]]
ggsave(scatter_plot_BMI_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_BMI_ON_education.pdf", width=7, height=7)


pleiotropy_BMI_ON_education <- mr_pleiotropy_test(UVdat_BMI_ON_education)
str(pleiotropy_BMI_ON_education)
pleiotropy_BMI_ON_education <- pleiotropy_BMI_ON_education[,-c(1,2)]


heterogeneity_BMI_ON_education <- mr_heterogeneity(UVdat_BMI_ON_education)
str(heterogeneity_BMI_ON_education)
heterogeneity_BMI_ON_education <- heterogeneity_BMI_ON_education[,-c(1,2)]


singleSNP_BMI_ON_education <- mr_singlesnp(UVdat_BMI_ON_education)
funnel_plot_BMI_ON_education <- mr_funnel_plot(singleSNP_BMI_ON_education)
funnel_plot_BMI_ON_education[[1]]
ggsave(funnel_plot_BMI_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_BMI_ON_education.pdf", width=7, height=7)


loo_BMI_ON_education <- mr_leaveoneout(UVdat_BMI_ON_education)
loo_plot_BMI_ON_education <- mr_leaveoneout_plot(loo_BMI_ON_education)
loo_plot_BMI_ON_education[[1]]
ggsave(loo_plot_BMI_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_BMI_ON_education.pdf", width=7, height=21)


singleSNP_plot_BMI_ON_education <- mr_forest_plot(singleSNP_BMI_ON_education)
singleSNP_plot_BMI_ON_education[[1]]
ggsave(singleSNP_plot_BMI_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_BMI_ON_education.pdf", width=7, height=21)

PRESSO_BMI_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_BMI_ON_education, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_BMI_ON_education <- generate_odds_ratios(result_BMI_ON_education)
result_BMI_ON_education <- result_BMI_ON_education[,-c(1,2)]


#### ________M3.2 lipid ####
IV_lipid <-format_data(GWAS_lipid,
                       type = "exposure", 
                       snp_col = "SNP", 
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col = "effect allele",
                       other_allele_col = "other allele",
                       eaf_col = "Freq",
                       pval_col = "P",
                       samplesize_col = "N")

str(IV_lipid)
IV_lipid$exposure <- "blood lipids"

IV_lipid <- IV_lipid[IV_lipid$pval.exposure < 5e-8,]
IV_lipid <- clump_data(IV_lipid)

IV_lipid$r2 <- (2 * (IV_lipid$beta.exposure^2) * IV_lipid$eaf.exposure * (1 - IV_lipid$eaf.exposure)) /
  (2 * (IV_lipid$beta.exposure^2) * IV_lipid$eaf.exposure * (1 - IV_lipid$eaf.exposure) +
     2 * IV_lipid$samplesize.exposure * IV_lipid$eaf.exposure * (1 - IV_lipid$eaf.exposure) * (IV_lipid$se.exposure^2))

IV_lipid$F <- (IV_lipid$r2 * (IV_lipid$samplesize.exposure - 2))/(1 - IV_lipid$r2)
range(IV_lipid$F)#[1]   29.76042 2837.00157
IV_lipid_meanF <- mean(IV_lipid$F)#[1] 170.7524

write.xlsx(IV_lipid, file = "Output/UVMR_IV/IV_lipid.xlsx")

out_lipid_ON_education <- extract_outcome_data(
  snps = IV_lipid$SNP,
  outcomes = 'ieu-a-1239')
str(out_lipid_ON_education)
out_lipid_ON_education$outcome <- "educational attainment"

UVdat_lipid_ON_education <- harmonise_data(
  exposure_dat =  IV_lipid, 
  outcome_dat = out_lipid_ON_education
)
range(UVdat_lipid_ON_education$pval.outcome)
str(UVdat_lipid_ON_education)
UVdat_lipid_ON_education <- UVdat_lipid_ON_education[UVdat_lipid_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_lipid_ON_education,file = "Output/UVMR_IV/UVdat_lipid_ON_education.xlsx")

result_lipid_ON_education <- mr(UVdat_lipid_ON_education)


scatter_plot_lipid_ON_education <- mr_scatter_plot(result_lipid_ON_education, UVdat_lipid_ON_education)
scatter_plot_lipid_ON_education[[1]]
ggsave(scatter_plot_lipid_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_lipid_ON_education.pdf", width=7, height=7)


pleiotropy_lipid_ON_education <- mr_pleiotropy_test(UVdat_lipid_ON_education)
str(pleiotropy_lipid_ON_education)
pleiotropy_lipid_ON_education <- pleiotropy_lipid_ON_education[,-c(1,2)]


heterogeneity_lipid_ON_education <- mr_heterogeneity(UVdat_lipid_ON_education)
str(heterogeneity_lipid_ON_education)
heterogeneity_lipid_ON_education <- heterogeneity_lipid_ON_education[,-c(1,2)]


singleSNP_lipid_ON_education <- mr_singlesnp(UVdat_lipid_ON_education)
funnel_plot_lipid_ON_education <- mr_funnel_plot(singleSNP_lipid_ON_education)
funnel_plot_lipid_ON_education[[1]]
ggsave(funnel_plot_lipid_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_lipid_ON_education.pdf", width=7, height=7)


loo_lipid_ON_education <- mr_leaveoneout(UVdat_lipid_ON_education)
loo_plot_lipid_ON_education <- mr_leaveoneout_plot(loo_lipid_ON_education)
loo_plot_lipid_ON_education[[1]]
ggsave(loo_plot_lipid_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_lipid_ON_education.pdf", width=7, height=21)


singleSNP_plot_lipid_ON_education <- mr_forest_plot(singleSNP_lipid_ON_education)
singleSNP_plot_lipid_ON_education[[1]]
ggsave(singleSNP_plot_lipid_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_lipid_ON_education.pdf", width=7, height=21)

PRESSO_lipid_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                       OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_lipid_ON_education, NbDistribution = 10000,  
                                       SignifThreshold = 0.05)

result_lipid_ON_education <- generate_odds_ratios(result_lipid_ON_education)
result_lipid_ON_education <- result_lipid_ON_education[,-c(1,2)]


#### ________M3.3 glucose ####
IV_glucose <-format_data(GWAS_glucose,
                         type = "exposure", 
                         snp_col = "SNP", 
                         beta_col = "BETA",
                         se_col = "SE",
                         effect_allele_col = "effect allele",
                         other_allele_col = "other allele",
                         eaf_col = "Freq",
                         pval_col = "P",
                         samplesize_col = "N")

str(IV_glucose)
IV_glucose$exposure <- "blood glucose"

IV_glucose <- IV_glucose[IV_glucose$pval.exposure < 5e-8,]
IV_glucose <- clump_data(IV_glucose)

IV_glucose$r2 <- (2 * (IV_glucose$beta.exposure^2) * IV_glucose$eaf.exposure * (1 - IV_glucose$eaf.exposure)) /
  (2 * (IV_glucose$beta.exposure^2) * IV_glucose$eaf.exposure * (1 - IV_glucose$eaf.exposure) +
     2 * IV_glucose$samplesize.exposure * IV_glucose$eaf.exposure * (1 - IV_glucose$eaf.exposure) * (IV_glucose$se.exposure^2))

IV_glucose$F <- (IV_glucose$r2 * (IV_glucose$samplesize.exposure - 2))/(1 - IV_glucose$r2)
range(IV_glucose$F)#[1]   24.99961 1198.46069
IV_glucose_meanF <- mean(IV_glucose$F)#[1] 100.4452

write.xlsx(IV_glucose, file = "Output/UVMR_IV/IV_glucose.xlsx")

out_glucose_ON_education <- extract_outcome_data(
  snps = IV_glucose$SNP,
  outcomes = 'ieu-a-1239')
str(out_glucose_ON_education)
out_glucose_ON_education$outcome <- "educational attainment"

UVdat_glucose_ON_education <- harmonise_data(
  exposure_dat =  IV_glucose, 
  outcome_dat = out_glucose_ON_education
)
range(UVdat_glucose_ON_education$pval.outcome)
str(UVdat_glucose_ON_education)
UVdat_glucose_ON_education <- UVdat_glucose_ON_education[UVdat_glucose_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_glucose_ON_education,file = "Output/UVMR_IV/UVdat_glucose_ON_education.xlsx")

result_glucose_ON_education <- mr(UVdat_glucose_ON_education)


scatter_plot_glucose_ON_education <- mr_scatter_plot(result_glucose_ON_education, UVdat_glucose_ON_education)
scatter_plot_glucose_ON_education[[1]]
ggsave(scatter_plot_glucose_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_glucose_ON_education.pdf", width=7, height=7)


pleiotropy_glucose_ON_education <- mr_pleiotropy_test(UVdat_glucose_ON_education)
str(pleiotropy_glucose_ON_education)
pleiotropy_glucose_ON_education <- pleiotropy_glucose_ON_education[,-c(1,2)]


heterogeneity_glucose_ON_education <- mr_heterogeneity(UVdat_glucose_ON_education)
str(heterogeneity_glucose_ON_education)
heterogeneity_glucose_ON_education <- heterogeneity_glucose_ON_education[,-c(1,2)]


singleSNP_glucose_ON_education <- mr_singlesnp(UVdat_glucose_ON_education)
funnel_plot_glucose_ON_education <- mr_funnel_plot(singleSNP_glucose_ON_education)
funnel_plot_glucose_ON_education[[1]]
ggsave(funnel_plot_glucose_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_glucose_ON_education.pdf", width=7, height=7)


loo_glucose_ON_education <- mr_leaveoneout(UVdat_glucose_ON_education)
loo_plot_glucose_ON_education <- mr_leaveoneout_plot(loo_glucose_ON_education)
loo_plot_glucose_ON_education[[1]]
ggsave(loo_plot_glucose_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_glucose_ON_education.pdf", width=7, height=21)


singleSNP_plot_glucose_ON_education <- mr_forest_plot(singleSNP_glucose_ON_education)
singleSNP_plot_glucose_ON_education[[1]]
ggsave(singleSNP_plot_glucose_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_glucose_ON_education.pdf", width=7, height=21)

PRESSO_glucose_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_glucose_ON_education, NbDistribution = 10000,  
                                         SignifThreshold = 0.05)

result_glucose_ON_education <- generate_odds_ratios(result_glucose_ON_education)
result_glucose_ON_education <- result_glucose_ON_education[,-c(1,2)]


#### ________M3.4 SBP, DBP ####

## SBP
IV_SBP <-format_data(GWAS_SBP,
                     type = "exposure", 
                     snp_col = "SNP", 
                     beta_col = "ES",
                     se_col = "SE",
                     effect_allele_col = "ALT",
                     other_allele_col = "REF",
                     eaf_col = "AF",
                     pval_col = "P",
                     samplesize_col = "SS")

str(IV_SBP)
IV_SBP$exposure <- "SBP"

IV_SBP <- IV_SBP[IV_SBP$pval.exposure < 5e-8,]
tt <- clump_data(IV_SBP)
IV_SBP <- clump_data(IV_SBP)

IV_SBP$r2 <- (2 * (IV_SBP$beta.exposure^2) * IV_SBP$eaf.exposure * (1 - IV_SBP$eaf.exposure)) /
  (2 * (IV_SBP$beta.exposure^2) * IV_SBP$eaf.exposure * (1 - IV_SBP$eaf.exposure) +
     2 * IV_SBP$samplesize.exposure * IV_SBP$eaf.exposure * (1 - IV_SBP$eaf.exposure) * (IV_SBP$se.exposure^2))

IV_SBP$F <- (IV_SBP$r2 * (IV_SBP$samplesize.exposure - 2))/(1 - IV_SBP$r2)
range(IV_SBP$F)#[1]  29.72041 627.54577
IV_SBP_meanF <- mean(IV_SBP$F)#[1] 75.12852

write.xlsx(IV_SBP, file = "Output/UVMR_IV/IV_SBP.xlsx")

out_SBP_ON_education <- extract_outcome_data(
  snps = IV_SBP$SNP,
  outcomes = 'ieu-a-1239')
str(out_SBP_ON_education)
out_SBP_ON_education$outcome <- "educational attainment"

UVdat_SBP_ON_education <- harmonise_data(
  exposure_dat =  IV_SBP, 
  outcome_dat = out_SBP_ON_education
)
range(UVdat_SBP_ON_education$pval.outcome)
str(UVdat_SBP_ON_education)
UVdat_SBP_ON_education <- UVdat_SBP_ON_education[UVdat_SBP_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_SBP_ON_education,file = "Output/UVMR_IV/UVdat_SBP_ON_education.xlsx")

result_SBP_ON_education <- mr(UVdat_SBP_ON_education)


scatter_plot_SBP_ON_education <- mr_scatter_plot(result_SBP_ON_education, UVdat_SBP_ON_education)
scatter_plot_SBP_ON_education[[1]]
ggsave(scatter_plot_SBP_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_SBP_ON_education.pdf", width=7, height=7)


pleiotropy_SBP_ON_education <- mr_pleiotropy_test(UVdat_SBP_ON_education)
str(pleiotropy_SBP_ON_education)
pleiotropy_SBP_ON_education <- pleiotropy_SBP_ON_education[,-c(1,2)]


heterogeneity_SBP_ON_education <- mr_heterogeneity(UVdat_SBP_ON_education)
str(heterogeneity_SBP_ON_education)
heterogeneity_SBP_ON_education <- heterogeneity_SBP_ON_education[,-c(1,2)]


singleSNP_SBP_ON_education <- mr_singlesnp(UVdat_SBP_ON_education)
funnel_plot_SBP_ON_education <- mr_funnel_plot(singleSNP_SBP_ON_education)
funnel_plot_SBP_ON_education[[1]]
ggsave(funnel_plot_SBP_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_SBP_ON_education.pdf", width=7, height=7)


loo_SBP_ON_education <- mr_leaveoneout(UVdat_SBP_ON_education)
loo_plot_SBP_ON_education <- mr_leaveoneout_plot(loo_SBP_ON_education)
loo_plot_SBP_ON_education[[1]]
ggsave(loo_plot_SBP_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_SBP_ON_education.pdf", width=7, height=21)


singleSNP_plot_SBP_ON_education <- mr_forest_plot(singleSNP_SBP_ON_education)
singleSNP_plot_SBP_ON_education[[1]]
ggsave(singleSNP_plot_SBP_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_SBP_ON_education.pdf", width=7, height=21)

PRESSO_SBP_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_SBP_ON_education, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_SBP_ON_education <- generate_odds_ratios(result_SBP_ON_education)
result_SBP_ON_education <- result_SBP_ON_education[,-c(1,2)]


## DBP
IV_DBP <-format_data(GWAS_DBP,
                     type = "exposure", 
                     snp_col = "SNP", 
                     beta_col = "ES",
                     se_col = "SE",
                     effect_allele_col = "ALT",
                     other_allele_col = "REF",
                     eaf_col = "AF",
                     pval_col = "P",
                     samplesize_col = "SS")

str(IV_DBP)
IV_DBP$exposure <- "DBP"

IV_DBP <- IV_DBP[IV_DBP$pval.exposure < 5e-8,]
IV_DBP <- clump_data(IV_DBP)

IV_DBP$r2 <- (2 * (IV_DBP$beta.exposure^2) * IV_DBP$eaf.exposure * (1 - IV_DBP$eaf.exposure)) /
  (2 * (IV_DBP$beta.exposure^2) * IV_DBP$eaf.exposure * (1 - IV_DBP$eaf.exposure) +
     2 * IV_DBP$samplesize.exposure * IV_DBP$eaf.exposure * (1 - IV_DBP$eaf.exposure) * (IV_DBP$se.exposure^2))

IV_DBP$F <- (IV_DBP$r2 * (IV_DBP$samplesize.exposure - 2))/(1 - IV_DBP$r2)
range(IV_DBP$F)#[1]  29.64768 815.81417
IV_DBP_meanF <- mean(IV_DBP$F)#[1] 79.41997

write.xlsx(IV_DBP, file = "Output/UVMR_IV/IV_DBP.xlsx")

out_DBP_ON_education <- extract_outcome_data(
  snps = IV_DBP$SNP,
  outcomes = 'ieu-a-1239')
str(out_DBP_ON_education)
out_DBP_ON_education$outcome <- "educational attainment"

UVdat_DBP_ON_education <- harmonise_data(
  exposure_dat =  IV_DBP, 
  outcome_dat = out_DBP_ON_education
)
range(UVdat_DBP_ON_education$pval.outcome)
str(UVdat_DBP_ON_education)
UVdat_DBP_ON_education <- UVdat_DBP_ON_education[UVdat_DBP_ON_education$pval.outcome>5e-8,]
write.xlsx(UVdat_DBP_ON_education,file = "Output/UVMR_IV/UVdat_DBP_ON_education.xlsx")

result_DBP_ON_education <- mr(UVdat_DBP_ON_education)


scatter_plot_DBP_ON_education <- mr_scatter_plot(result_DBP_ON_education, UVdat_DBP_ON_education)
scatter_plot_DBP_ON_education[[1]]
ggsave(scatter_plot_DBP_ON_education[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DBP_ON_education.pdf", width=7, height=7)


pleiotropy_DBP_ON_education <- mr_pleiotropy_test(UVdat_DBP_ON_education)
str(pleiotropy_DBP_ON_education)
pleiotropy_DBP_ON_education <- pleiotropy_DBP_ON_education[,-c(1,2)]


heterogeneity_DBP_ON_education <- mr_heterogeneity(UVdat_DBP_ON_education)
str(heterogeneity_DBP_ON_education)
heterogeneity_DBP_ON_education <- heterogeneity_DBP_ON_education[,-c(1,2)]


singleSNP_DBP_ON_education <- mr_singlesnp(UVdat_DBP_ON_education)
funnel_plot_DBP_ON_education <- mr_funnel_plot(singleSNP_DBP_ON_education)
funnel_plot_DBP_ON_education[[1]]
ggsave(funnel_plot_DBP_ON_education[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DBP_ON_education.pdf", width=7, height=7)


loo_DBP_ON_education <- mr_leaveoneout(UVdat_DBP_ON_education)
loo_plot_DBP_ON_education <- mr_leaveoneout_plot(loo_DBP_ON_education)
loo_plot_DBP_ON_education[[1]]
ggsave(loo_plot_DBP_ON_education[[1]], file="Output/UVMR_Secondary Results/loo_plot_DBP_ON_education.pdf", width=7, height=21)


singleSNP_plot_DBP_ON_education <- mr_forest_plot(singleSNP_DBP_ON_education)
singleSNP_plot_DBP_ON_education[[1]]
ggsave(singleSNP_plot_DBP_ON_education[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DBP_ON_education.pdf", width=7, height=21)

PRESSO_DBP_ON_education <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DBP_ON_education, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_DBP_ON_education <- generate_odds_ratios(result_DBP_ON_education)
result_DBP_ON_education <- result_DBP_ON_education[,-c(1,2)]

###### Main Results
SummaryResultsUVMR_education_ON_Mediator <- rbind(result_education_ON_carbohydrate,
                                                  result_education_ON_protein,
                                                  result_education_ON_fat,
                                                  result_education_ON_FreshFruit,
                                                  result_education_ON_DriedFruit,
                                                  result_education_ON_RawVegetables,
                                                  result_education_ON_CookedVegetables,
                                                  result_education_ON_salt,
                                                  result_education_ON_processed,
                                                  result_education_ON_OilyFish,
                                                  result_education_ON_wholegrain,
                                                  result_education_ON_DeviceOverallActivity,
                                                  result_education_ON_PA,
                                                  result_education_ON_DeviceModerate,
                                                  result_education_ON_DeviceSedentary,
                                                  result_education_ON_LST,
                                                  result_education_ON_television,
                                                  result_education_ON_computer,
                                                  result_education_ON_driving,
                                                  result_education_ON_LifetimeSmoking,
                                                  result_education_ON_SleepDuration,
                                                  result_education_ON_DeviceSleepDuration,
                                                  result_education_ON_Insomnia,
                                                  result_education_ON_WBS,
                                                  result_education_ON_MentalProblems,
                                                  result_education_ON_depression,
                                                  result_education_ON_BMI,
                                                  result_education_ON_lipid,
                                                  result_education_ON_glucose,
                                                  result_education_ON_SBP,
                                                  result_education_ON_DBP)
write.xlsx(SummaryResultsUVMR_education_ON_Mediator, file = "Output/UVMR_Main Results/SummaryResultsUVMR_education_ON_Mediator.xlsx")

SummaryResultsUVMR_Mediator_ON_education <- rbind(result_protein_ON_education,
                                                  result_fat_ON_education,
                                                  result_FreshFruit_ON_education,
                                                  result_DriedFruit_ON_education,
                                                  result_RawVegetables_ON_education,
                                                  result_salt_ON_education,
                                                  result_processed_ON_education,
                                                  result_OilyFish_ON_education,
                                                  result_wholegrain_ON_education,
                                                  result_DeviceOverallActivity_ON_education,
                                                  result_PA_ON_education,
                                                  result_DeviceModerate_ON_education,
                                                  result_DeviceSedentary_ON_education,
                                                  result_LST_ON_education,
                                                  result_television_ON_education,
                                                  result_computer_ON_education,
                                                  result_driving_ON_education,
                                                  result_LifetimeSmoking_ON_education,
                                                  result_SleepDuration_ON_education,
                                                  result_DeviceSleepDuration_ON_education,
                                                  result_Insomnia_ON_education,
                                                  result_WBS_ON_education,
                                                  result_MentalProblems_ON_education,
                                                  result_depression_ON_education,
                                                  result_BMI_ON_education,
                                                  result_lipid_ON_education,
                                                  result_glucose_ON_education,
                                                  result_SBP_ON_education,
                                                  result_DBP_ON_education)
write.xlsx(SummaryResultsUVMR_Mediator_ON_education, file = "Output/UVMR_Main Results/SummaryResultsUVMR_Mediator_ON_education.xlsx")


###### Secondary Results

Summary_UVPleiotropy_education_ON_Mediator <- rbind(pleiotropy_education_ON_carbohydrate,
                                                    pleiotropy_education_ON_protein,
                                                    pleiotropy_education_ON_fat,
                                                    pleiotropy_education_ON_FreshFruit,
                                                    pleiotropy_education_ON_DriedFruit,
                                                    pleiotropy_education_ON_RawVegetables,
                                                    pleiotropy_education_ON_CookedVegetables,
                                                    pleiotropy_education_ON_salt,
                                                    pleiotropy_education_ON_processed,
                                                    pleiotropy_education_ON_OilyFish,
                                                    pleiotropy_education_ON_wholegrain,
                                                    pleiotropy_education_ON_DeviceOverallActivity,
                                                    pleiotropy_education_ON_PA,
                                                    pleiotropy_education_ON_DeviceModerate,
                                                    pleiotropy_education_ON_DeviceSedentary,
                                                    pleiotropy_education_ON_LST,
                                                    pleiotropy_education_ON_television,
                                                    pleiotropy_education_ON_computer,
                                                    pleiotropy_education_ON_driving,
                                                    pleiotropy_education_ON_LifetimeSmoking,
                                                    pleiotropy_education_ON_SleepDuration,
                                                    pleiotropy_education_ON_DeviceSleepDuration,
                                                    pleiotropy_education_ON_Insomnia,
                                                    pleiotropy_education_ON_WBS,
                                                    pleiotropy_education_ON_MentalProblems,
                                                    pleiotropy_education_ON_depression,
                                                    pleiotropy_education_ON_BMI,
                                                    pleiotropy_education_ON_lipid,
                                                    pleiotropy_education_ON_glucose,
                                                    pleiotropy_education_ON_SBP,
                                                    pleiotropy_education_ON_DBP)
write.xlsx(Summary_UVPleiotropy_education_ON_Mediator,"Output/UVMR_Secondary Results/Summary_UVPleiotropy_education_ON_Mediator.xlsx")

Summary_UVPleiotropy_Mediator_ON_education <- rbind(pleiotropy_protein_ON_education,
                                                    pleiotropy_fat_ON_education,
                                                    pleiotropy_FreshFruit_ON_education,
                                                    pleiotropy_DriedFruit_ON_education,
                                                    pleiotropy_RawVegetables_ON_education,
                                                    pleiotropy_salt_ON_education,
                                                    pleiotropy_processed_ON_education,
                                                    pleiotropy_OilyFish_ON_education,
                                                    pleiotropy_wholegrain_ON_education,
                                                    pleiotropy_DeviceOverallActivity_ON_education,
                                                    pleiotropy_PA_ON_education,
                                                    pleiotropy_DeviceSedentary_ON_education,
                                                    pleiotropy_LST_ON_education,
                                                    pleiotropy_television_ON_education,
                                                    pleiotropy_computer_ON_education,
                                                    pleiotropy_driving_ON_education,
                                                    pleiotropy_LifetimeSmoking_ON_education,
                                                    pleiotropy_SleepDuration_ON_education,
                                                    pleiotropy_DeviceSleepDuration_ON_education,
                                                    pleiotropy_Insomnia_ON_education,
                                                    pleiotropy_WBS_ON_education,
                                                    pleiotropy_MentalProblems_ON_education,
                                                    pleiotropy_depression_ON_education,
                                                    pleiotropy_BMI_ON_education,
                                                    pleiotropy_lipid_ON_education,
                                                    pleiotropy_glucose_ON_education,
                                                    pleiotropy_SBP_ON_education,
                                                    pleiotropy_DBP_ON_education)
write.xlsx(Summary_UVPleiotropy_Mediator_ON_education,"Output/UVMR_Secondary Results/Summary_UVPleiotropy_Mediator_ON_education.xlsx")


Summary_UVHeterogeneity_education_ON_Mediator <- rbind(heterogeneity_education_ON_carbohydrate,
                                                       heterogeneity_education_ON_protein,
                                                       heterogeneity_education_ON_fat,
                                                       heterogeneity_education_ON_FreshFruit,
                                                       heterogeneity_education_ON_DriedFruit,
                                                       heterogeneity_education_ON_RawVegetables,
                                                       heterogeneity_education_ON_CookedVegetables,
                                                       heterogeneity_education_ON_salt,
                                                       heterogeneity_education_ON_processed,
                                                       heterogeneity_education_ON_OilyFish,
                                                       heterogeneity_education_ON_wholegrain,
                                                       heterogeneity_education_ON_DeviceOverallActivity,
                                                       heterogeneity_education_ON_PA,
                                                       heterogeneity_education_ON_DeviceModerate,
                                                       heterogeneity_education_ON_DeviceSedentary,
                                                       heterogeneity_education_ON_LST,
                                                       heterogeneity_education_ON_television,
                                                       heterogeneity_education_ON_computer,
                                                       heterogeneity_education_ON_driving,
                                                       heterogeneity_education_ON_LifetimeSmoking,
                                                       heterogeneity_education_ON_SleepDuration,
                                                       heterogeneity_education_ON_DeviceSleepDuration,
                                                       heterogeneity_education_ON_Insomnia,
                                                       heterogeneity_education_ON_WBS,
                                                       heterogeneity_education_ON_MentalProblems,
                                                       heterogeneity_education_ON_depression,
                                                       heterogeneity_education_ON_BMI,
                                                       heterogeneity_education_ON_lipid,
                                                       heterogeneity_education_ON_glucose,
                                                       heterogeneity_education_ON_SBP,
                                                       heterogeneity_education_ON_DBP)
write.xlsx(Summary_UVHeterogeneity_education_ON_Mediator,"Output/UVMR_Secondary Results/Summary_UVHeterogeneity_education_ON_Mediator.xlsx")

Summary_UVHeterogeneity_Mediator_ON_education <- rbind(heterogeneity_protein_ON_education,
                                                       heterogeneity_fat_ON_education,
                                                       heterogeneity_FreshFruit_ON_education,
                                                       heterogeneity_DriedFruit_ON_education,
                                                       heterogeneity_RawVegetables_ON_education,
                                                       heterogeneity_salt_ON_education,
                                                       heterogeneity_processed_ON_education,
                                                       heterogeneity_OilyFish_ON_education,
                                                       heterogeneity_wholegrain_ON_education,
                                                       heterogeneity_DeviceOverallActivity_ON_education,
                                                       heterogeneity_PA_ON_education,
                                                       heterogeneity_DeviceSedentary_ON_education,
                                                       heterogeneity_LST_ON_education,
                                                       heterogeneity_television_ON_education,
                                                       heterogeneity_computer_ON_education,
                                                       heterogeneity_driving_ON_education,
                                                       heterogeneity_LifetimeSmoking_ON_education,
                                                       heterogeneity_SleepDuration_ON_education,
                                                       heterogeneity_DeviceSleepDuration_ON_education,
                                                       heterogeneity_Insomnia_ON_education,
                                                       heterogeneity_WBS_ON_education,
                                                       heterogeneity_MentalProblems_ON_education,
                                                       heterogeneity_depression_ON_education,
                                                       heterogeneity_BMI_ON_education,
                                                       heterogeneity_lipid_ON_education,
                                                       heterogeneity_glucose_ON_education,
                                                       heterogeneity_SBP_ON_education,
                                                       heterogeneity_DBP_ON_education)
write.xlsx(Summary_UVHeterogeneity_Mediator_ON_education,"Output/UVMR_Secondary Results/Summary_UVHeterogeneity_Mediator_ON_education.xlsx")



## mediation screening criteria II: the causal relationship between the mediator and Y is consistent whether X is adjusted or not
# the mediator should be causally associated with the outcome
# +the mediator should have a direct causal effect on the outcome independently of education

## the mediator should be causally associated with the outcome in UVMR

#### ________M1.1 Diet_protein ####
str(IV_protein)

out_protein_ON_CAD <- extract_outcome_data(
  snps = IV_protein$SNP,
  outcomes = 'ieu-a-7')
str(out_protein_ON_CAD)
out_protein_ON_CAD$outcome <- "coronary heart disease"

UVdat_protein_ON_CAD <- harmonise_data(
  exposure_dat =  IV_protein, 
  outcome_dat = out_protein_ON_CAD
)
str(UVdat_protein_ON_CAD)
range(UVdat_protein_ON_CAD$F)
range(UVdat_protein_ON_CAD$pval.outcome)
UVdat_protein_ON_CAD <- UVdat_protein_ON_CAD[UVdat_protein_ON_CAD$pval.outcome>5e-8,]

write.xlsx(UVdat_protein_ON_CAD,file = "Output/UVMR_IV/UVdat_protein_ON_CAD.xlsx")

result_protein_ON_CAD <- mr(UVdat_protein_ON_CAD)


scatter_plot_protein_ON_CAD <- mr_scatter_plot(result_protein_ON_CAD, UVdat_protein_ON_CAD)
scatter_plot_protein_ON_CAD[[1]]
ggsave(scatter_plot_protein_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_protein_ON_CAD.pdf", width=7, height=7)


pleiotropy_protein_ON_CAD <- mr_pleiotropy_test(UVdat_protein_ON_CAD)
str(pleiotropy_protein_ON_CAD)
pleiotropy_protein_ON_CAD <- pleiotropy_protein_ON_CAD[,-c(1,2)]


heterogeneity_protein_ON_CAD <- mr_heterogeneity(UVdat_protein_ON_CAD)
str(heterogeneity_protein_ON_CAD)
heterogeneity_protein_ON_CAD <- heterogeneity_protein_ON_CAD[,-c(1,2)]


singleSNP_protein_ON_CAD <- mr_singlesnp(UVdat_protein_ON_CAD)
funnel_plot_protein_ON_CAD <- mr_funnel_plot(singleSNP_protein_ON_CAD)
funnel_plot_protein_ON_CAD[[1]]
ggsave(funnel_plot_protein_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_protein_ON_CAD.pdf", width=7, height=7)


loo_protein_ON_CAD <- mr_leaveoneout(UVdat_protein_ON_CAD)
loo_plot_protein_ON_CAD <- mr_leaveoneout_plot(loo_protein_ON_CAD)
loo_plot_protein_ON_CAD[[1]]
ggsave(loo_plot_protein_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_protein_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_protein_ON_CAD <- mr_forest_plot(singleSNP_protein_ON_CAD)
singleSNP_plot_protein_ON_CAD[[1]]
ggsave(singleSNP_plot_protein_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_protein_ON_CAD.pdf", width=7, height=7)

PRESSO_protein_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_protein_ON_CAD, NbDistribution = 10000,  
                                      SignifThreshold = 0.05)

result_protein_ON_CAD <- generate_odds_ratios(result_protein_ON_CAD)
result_protein_ON_CAD <- result_protein_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_fat ####
str(IV_fat)

out_fat_ON_CAD <- extract_outcome_data(
  snps = IV_fat$SNP,
  outcomes = 'ieu-a-7')
str(out_fat_ON_CAD)
out_fat_ON_CAD$outcome <- "coronary heart disease"

UVdat_fat_ON_CAD <- harmonise_data(
  exposure_dat =  IV_fat, 
  outcome_dat = out_fat_ON_CAD
)
str(UVdat_fat_ON_CAD)
range(UVdat_fat_ON_CAD$F)
range(UVdat_fat_ON_CAD$pval.outcome)
UVdat_fat_ON_CAD <- UVdat_fat_ON_CAD[UVdat_fat_ON_CAD$pval.outcome>5e-8,]

write.xlsx(UVdat_fat_ON_CAD,file = "Output/UVMR_IV/UVdat_fat_ON_CAD.xlsx")

result_fat_ON_CAD <- mr(UVdat_fat_ON_CAD)


scatter_plot_fat_ON_CAD <- mr_scatter_plot(result_fat_ON_CAD, UVdat_fat_ON_CAD)
scatter_plot_fat_ON_CAD[[1]]
ggsave(scatter_plot_fat_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_fat_ON_CAD.pdf", width=7, height=7)


pleiotropy_fat_ON_CAD <- mr_pleiotropy_test(UVdat_fat_ON_CAD)
str(pleiotropy_fat_ON_CAD)
pleiotropy_fat_ON_CAD <- pleiotropy_fat_ON_CAD[,-c(1,2)]


heterogeneity_fat_ON_CAD <- mr_heterogeneity(UVdat_fat_ON_CAD)
str(heterogeneity_fat_ON_CAD)
heterogeneity_fat_ON_CAD <- heterogeneity_fat_ON_CAD[,-c(1,2)]


singleSNP_fat_ON_CAD <- mr_singlesnp(UVdat_fat_ON_CAD)
funnel_plot_fat_ON_CAD <- mr_funnel_plot(singleSNP_fat_ON_CAD)
funnel_plot_fat_ON_CAD[[1]]
ggsave(funnel_plot_fat_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_fat_ON_CAD.pdf", width=7, height=7)


loo_fat_ON_CAD <- mr_leaveoneout(UVdat_fat_ON_CAD)
loo_plot_fat_ON_CAD <- mr_leaveoneout_plot(loo_fat_ON_CAD)
loo_plot_fat_ON_CAD[[1]]
ggsave(loo_plot_fat_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_fat_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_fat_ON_CAD <- mr_forest_plot(singleSNP_fat_ON_CAD)
singleSNP_plot_fat_ON_CAD[[1]]
ggsave(singleSNP_plot_fat_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_fat_ON_CAD.pdf", width=7, height=7)

PRESSO_fat_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_fat_ON_CAD, NbDistribution = 10000,  
                                   SignifThreshold = 0.05)

result_fat_ON_CAD <- generate_odds_ratios(result_fat_ON_CAD)
result_fat_ON_CAD <- result_fat_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_FreshFruit ####
str(IV_FreshFruit)

out_FreshFruit_ON_CAD <- extract_outcome_data(
  snps = IV_FreshFruit$SNP,
  outcomes = 'ieu-a-7')
str(out_FreshFruit_ON_CAD)
out_FreshFruit_ON_CAD$outcome <- "coronary heart disease"

UVdat_FreshFruit_ON_CAD <- harmonise_data(
  exposure_dat =  IV_FreshFruit, 
  outcome_dat = out_FreshFruit_ON_CAD
)
str(UVdat_FreshFruit_ON_CAD)
range(UVdat_FreshFruit_ON_CAD$F)
range(UVdat_FreshFruit_ON_CAD$pval.outcome)

write.xlsx(UVdat_FreshFruit_ON_CAD,file = "Output/UVMR_IV/UVdat_FreshFruit_ON_CAD.xlsx")

result_FreshFruit_ON_CAD <- mr(UVdat_FreshFruit_ON_CAD)


scatter_plot_FreshFruit_ON_CAD <- mr_scatter_plot(result_FreshFruit_ON_CAD, UVdat_FreshFruit_ON_CAD)
scatter_plot_FreshFruit_ON_CAD[[1]]
ggsave(scatter_plot_FreshFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_FreshFruit_ON_CAD.pdf", width=7, height=7)


pleiotropy_FreshFruit_ON_CAD <- mr_pleiotropy_test(UVdat_FreshFruit_ON_CAD)
str(pleiotropy_FreshFruit_ON_CAD)
pleiotropy_FreshFruit_ON_CAD <- pleiotropy_FreshFruit_ON_CAD[,-c(1,2)]


heterogeneity_FreshFruit_ON_CAD <- mr_heterogeneity(UVdat_FreshFruit_ON_CAD)
str(heterogeneity_FreshFruit_ON_CAD)
heterogeneity_FreshFruit_ON_CAD <- heterogeneity_FreshFruit_ON_CAD[,-c(1,2)]


singleSNP_FreshFruit_ON_CAD <- mr_singlesnp(UVdat_FreshFruit_ON_CAD)
funnel_plot_FreshFruit_ON_CAD <- mr_funnel_plot(singleSNP_FreshFruit_ON_CAD)
funnel_plot_FreshFruit_ON_CAD[[1]]
ggsave(funnel_plot_FreshFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_FreshFruit_ON_CAD.pdf", width=7, height=7)


loo_FreshFruit_ON_CAD <- mr_leaveoneout(UVdat_FreshFruit_ON_CAD)
loo_plot_FreshFruit_ON_CAD <- mr_leaveoneout_plot(loo_FreshFruit_ON_CAD)
loo_plot_FreshFruit_ON_CAD[[1]]
ggsave(loo_plot_FreshFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_FreshFruit_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_FreshFruit_ON_CAD <- mr_forest_plot(singleSNP_FreshFruit_ON_CAD)
singleSNP_plot_FreshFruit_ON_CAD[[1]]
ggsave(singleSNP_plot_FreshFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_FreshFruit_ON_CAD.pdf", width=7, height=7)

PRESSO_FreshFruit_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_FreshFruit_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_FreshFruit_ON_CAD <- generate_odds_ratios(result_FreshFruit_ON_CAD)
result_FreshFruit_ON_CAD <- result_FreshFruit_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_DriedFruit ####
str(IV_DriedFruit)

out_DriedFruit_ON_CAD <- extract_outcome_data(
  snps = IV_DriedFruit$SNP,
  outcomes = 'ieu-a-7')
str(out_DriedFruit_ON_CAD)
out_DriedFruit_ON_CAD$outcome <- "coronary heart disease"

UVdat_DriedFruit_ON_CAD <- harmonise_data(
  exposure_dat =  IV_DriedFruit, 
  outcome_dat = out_DriedFruit_ON_CAD
)
str(UVdat_DriedFruit_ON_CAD)
range(UVdat_DriedFruit_ON_CAD$F)
range(UVdat_DriedFruit_ON_CAD$pval.outcome)
UVdat_DriedFruit_ON_CAD <- UVdat_DriedFruit_ON_CAD[UVdat_DriedFruit_ON_CAD$pval.outcome>5e-8,]
write.xlsx(UVdat_DriedFruit_ON_CAD,file = "Output/UVMR_IV/UVdat_DriedFruit_ON_CAD.xlsx")

result_DriedFruit_ON_CAD <- mr(UVdat_DriedFruit_ON_CAD)


scatter_plot_DriedFruit_ON_CAD <- mr_scatter_plot(result_DriedFruit_ON_CAD, UVdat_DriedFruit_ON_CAD)
scatter_plot_DriedFruit_ON_CAD[[1]]
ggsave(scatter_plot_DriedFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DriedFruit_ON_CAD.pdf", width=7, height=7)


pleiotropy_DriedFruit_ON_CAD <- mr_pleiotropy_test(UVdat_DriedFruit_ON_CAD)
str(pleiotropy_DriedFruit_ON_CAD)
pleiotropy_DriedFruit_ON_CAD <- pleiotropy_DriedFruit_ON_CAD[,-c(1,2)]


heterogeneity_DriedFruit_ON_CAD <- mr_heterogeneity(UVdat_DriedFruit_ON_CAD)
str(heterogeneity_DriedFruit_ON_CAD)
heterogeneity_DriedFruit_ON_CAD <- heterogeneity_DriedFruit_ON_CAD[,-c(1,2)]


singleSNP_DriedFruit_ON_CAD <- mr_singlesnp(UVdat_DriedFruit_ON_CAD)
funnel_plot_DriedFruit_ON_CAD <- mr_funnel_plot(singleSNP_DriedFruit_ON_CAD)
funnel_plot_DriedFruit_ON_CAD[[1]]
ggsave(funnel_plot_DriedFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DriedFruit_ON_CAD.pdf", width=7, height=7)


loo_DriedFruit_ON_CAD <- mr_leaveoneout(UVdat_DriedFruit_ON_CAD)
loo_plot_DriedFruit_ON_CAD <- mr_leaveoneout_plot(loo_DriedFruit_ON_CAD)
loo_plot_DriedFruit_ON_CAD[[1]]
ggsave(loo_plot_DriedFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_DriedFruit_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_DriedFruit_ON_CAD <- mr_forest_plot(singleSNP_DriedFruit_ON_CAD)
singleSNP_plot_DriedFruit_ON_CAD[[1]]
ggsave(singleSNP_plot_DriedFruit_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DriedFruit_ON_CAD.pdf", width=7, height=7)

PRESSO_DriedFruit_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DriedFruit_ON_CAD, NbDistribution = 10000,  
                              SignifThreshold = 0.05)

result_DriedFruit_ON_CAD <- generate_odds_ratios(result_DriedFruit_ON_CAD)
result_DriedFruit_ON_CAD <- result_DriedFruit_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_RawVegetables ####
str(IV_RawVegetables)

out_RawVegetables_ON_CAD <- extract_outcome_data(
  snps = IV_RawVegetables$SNP,
  outcomes = 'ieu-a-7')
str(out_RawVegetables_ON_CAD)
out_RawVegetables_ON_CAD$outcome <- "coronary heart disease"

UVdat_RawVegetables_ON_CAD <- harmonise_data(
  exposure_dat =  IV_RawVegetables, 
  outcome_dat = out_RawVegetables_ON_CAD
)
str(UVdat_RawVegetables_ON_CAD)
range(UVdat_RawVegetables_ON_CAD$F)
range(UVdat_RawVegetables_ON_CAD$pval.outcome)

write.xlsx(UVdat_RawVegetables_ON_CAD,file = "Output/UVMR_IV/UVdat_RawVegetables_ON_CAD.xlsx")

result_RawVegetables_ON_CAD <- mr(UVdat_RawVegetables_ON_CAD)


scatter_plot_RawVegetables_ON_CAD <- mr_scatter_plot(result_RawVegetables_ON_CAD, UVdat_RawVegetables_ON_CAD)
scatter_plot_RawVegetables_ON_CAD[[1]]
ggsave(scatter_plot_RawVegetables_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_RawVegetables_ON_CAD.pdf", width=7, height=7)


pleiotropy_RawVegetables_ON_CAD <- mr_pleiotropy_test(UVdat_RawVegetables_ON_CAD)
str(pleiotropy_RawVegetables_ON_CAD)
pleiotropy_RawVegetables_ON_CAD <- pleiotropy_RawVegetables_ON_CAD[,-c(1,2)]


heterogeneity_RawVegetables_ON_CAD <- mr_heterogeneity(UVdat_RawVegetables_ON_CAD)
str(heterogeneity_RawVegetables_ON_CAD)
heterogeneity_RawVegetables_ON_CAD <- heterogeneity_RawVegetables_ON_CAD[,-c(1,2)]


singleSNP_RawVegetables_ON_CAD <- mr_singlesnp(UVdat_RawVegetables_ON_CAD)
funnel_plot_RawVegetables_ON_CAD <- mr_funnel_plot(singleSNP_RawVegetables_ON_CAD)
funnel_plot_RawVegetables_ON_CAD[[1]]
ggsave(funnel_plot_RawVegetables_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_RawVegetables_ON_CAD.pdf", width=7, height=7)


loo_RawVegetables_ON_CAD <- mr_leaveoneout(UVdat_RawVegetables_ON_CAD)
loo_plot_RawVegetables_ON_CAD <- mr_leaveoneout_plot(loo_RawVegetables_ON_CAD)
loo_plot_RawVegetables_ON_CAD[[1]]
ggsave(loo_plot_RawVegetables_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_RawVegetables_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_RawVegetables_ON_CAD <- mr_forest_plot(singleSNP_RawVegetables_ON_CAD)
singleSNP_plot_RawVegetables_ON_CAD[[1]]
ggsave(singleSNP_plot_RawVegetables_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_RawVegetables_ON_CAD.pdf", width=7, height=7)

PRESSO_RawVegetables_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_RawVegetables_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_RawVegetables_ON_CAD <- generate_odds_ratios(result_RawVegetables_ON_CAD)
result_RawVegetables_ON_CAD <- result_RawVegetables_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_salt ####
str(IV_salt)

out_salt_ON_CAD <- extract_outcome_data(
  snps = IV_salt$SNP,
  outcomes = 'ieu-a-7')

str(out_salt_ON_CAD)
out_salt_ON_CAD$outcome <- "coronary heart disease"

UVdat_salt_ON_CAD <- harmonise_data(
  exposure_dat =  IV_salt, 
  outcome_dat = out_salt_ON_CAD
)
str(UVdat_salt_ON_CAD)
range(UVdat_salt_ON_CAD$F)
range(UVdat_salt_ON_CAD$pval.outcome)

write.xlsx(UVdat_salt_ON_CAD,file = "Output/UVMR_IV/UVdat_salt_ON_CAD.xlsx")

result_salt_ON_CAD <- mr(UVdat_salt_ON_CAD)


scatter_plot_salt_ON_CAD <- mr_scatter_plot(result_salt_ON_CAD, UVdat_salt_ON_CAD)
scatter_plot_salt_ON_CAD[[1]]
ggsave(scatter_plot_salt_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_salt_ON_CAD.pdf", width=7, height=7)


pleiotropy_salt_ON_CAD <- mr_pleiotropy_test(UVdat_salt_ON_CAD)
str(pleiotropy_salt_ON_CAD)
pleiotropy_salt_ON_CAD <- pleiotropy_salt_ON_CAD[,-c(1,2)]


heterogeneity_salt_ON_CAD <- mr_heterogeneity(UVdat_salt_ON_CAD)
str(heterogeneity_salt_ON_CAD)
heterogeneity_salt_ON_CAD <- heterogeneity_salt_ON_CAD[,-c(1,2)]


singleSNP_salt_ON_CAD <- mr_singlesnp(UVdat_salt_ON_CAD)
funnel_plot_salt_ON_CAD <- mr_funnel_plot(singleSNP_salt_ON_CAD)
funnel_plot_salt_ON_CAD[[1]]
ggsave(funnel_plot_salt_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_salt_ON_CAD.pdf", width=7, height=7)


loo_salt_ON_CAD <- mr_leaveoneout(UVdat_salt_ON_CAD)
loo_plot_salt_ON_CAD <- mr_leaveoneout_plot(loo_salt_ON_CAD)
loo_plot_salt_ON_CAD[[1]]
ggsave(loo_plot_salt_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_salt_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_salt_ON_CAD <- mr_forest_plot(singleSNP_salt_ON_CAD)
singleSNP_plot_salt_ON_CAD[[1]]
ggsave(singleSNP_plot_salt_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_salt_ON_CAD.pdf", width=7, height=7)

PRESSO_salt_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_salt_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_salt_ON_CAD <- generate_odds_ratios(result_salt_ON_CAD)
result_salt_ON_CAD <- result_salt_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_processed ####
str(IV_processed)

out_processed_ON_CAD <- extract_outcome_data(
  snps = IV_processed$SNP,
  outcomes = 'ieu-a-7')

str(out_processed_ON_CAD)
out_processed_ON_CAD$outcome <- "coronary heart disease"

UVdat_processed_ON_CAD <- harmonise_data(
  exposure_dat =  IV_processed, 
  outcome_dat = out_processed_ON_CAD
)
str(UVdat_processed_ON_CAD)
range(UVdat_processed_ON_CAD$F)
range(UVdat_processed_ON_CAD$pval.outcome)

write.xlsx(UVdat_processed_ON_CAD,file = "Output/UVMR_IV/UVdat_processed_ON_CAD.xlsx")

result_processed_ON_CAD <- mr(UVdat_processed_ON_CAD)


scatter_plot_processed_ON_CAD <- mr_scatter_plot(result_processed_ON_CAD, UVdat_processed_ON_CAD)
scatter_plot_processed_ON_CAD[[1]]
ggsave(scatter_plot_processed_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_processed_ON_CAD.pdf", width=7, height=7)


pleiotropy_processed_ON_CAD <- mr_pleiotropy_test(UVdat_processed_ON_CAD)
str(pleiotropy_processed_ON_CAD)
pleiotropy_processed_ON_CAD <- pleiotropy_processed_ON_CAD[,-c(1,2)]


heterogeneity_processed_ON_CAD <- mr_heterogeneity(UVdat_processed_ON_CAD)
str(heterogeneity_processed_ON_CAD)
heterogeneity_processed_ON_CAD <- heterogeneity_processed_ON_CAD[,-c(1,2)]


singleSNP_processed_ON_CAD <- mr_singlesnp(UVdat_processed_ON_CAD)
funnel_plot_processed_ON_CAD <- mr_funnel_plot(singleSNP_processed_ON_CAD)
funnel_plot_processed_ON_CAD[[1]]
ggsave(funnel_plot_processed_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_processed_ON_CAD.pdf", width=7, height=7)


loo_processed_ON_CAD <- mr_leaveoneout(UVdat_processed_ON_CAD)
loo_plot_processed_ON_CAD <- mr_leaveoneout_plot(loo_processed_ON_CAD)
loo_plot_processed_ON_CAD[[1]]
ggsave(loo_plot_processed_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_processed_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_processed_ON_CAD <- mr_forest_plot(singleSNP_processed_ON_CAD)
singleSNP_plot_processed_ON_CAD[[1]]
ggsave(singleSNP_plot_processed_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_processed_ON_CAD.pdf", width=7, height=7)

PRESSO_processed_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_processed_ON_CAD, NbDistribution = 10000,  
                                SignifThreshold = 0.05)

result_processed_ON_CAD <- generate_odds_ratios(result_processed_ON_CAD)
result_processed_ON_CAD <- result_processed_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_OilyFish ####
str(IV_OilyFish)

out_OilyFish_ON_CAD <- extract_outcome_data(
  snps = IV_OilyFish$SNP,
  outcomes = 'ieu-a-7')
str(out_OilyFish_ON_CAD)
out_OilyFish_ON_CAD$outcome <- "coronary heart disease"

UVdat_OilyFish_ON_CAD <- harmonise_data(
  exposure_dat =  IV_OilyFish, 
  outcome_dat = out_OilyFish_ON_CAD
)
str(UVdat_OilyFish_ON_CAD)
range(UVdat_OilyFish_ON_CAD$F)
range(UVdat_OilyFish_ON_CAD$pval.outcome)

write.xlsx(UVdat_OilyFish_ON_CAD,file = "Output/UVMR_IV/UVdat_OilyFish_ON_CAD.xlsx")

result_OilyFish_ON_CAD <- mr(UVdat_OilyFish_ON_CAD)


scatter_plot_OilyFish_ON_CAD <- mr_scatter_plot(result_OilyFish_ON_CAD, UVdat_OilyFish_ON_CAD)
scatter_plot_OilyFish_ON_CAD[[1]]
ggsave(scatter_plot_OilyFish_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_OilyFish_ON_CAD.pdf", width=7, height=7)


pleiotropy_OilyFish_ON_CAD <- mr_pleiotropy_test(UVdat_OilyFish_ON_CAD)
str(pleiotropy_OilyFish_ON_CAD)
pleiotropy_OilyFish_ON_CAD <- pleiotropy_OilyFish_ON_CAD[,-c(1,2)]


heterogeneity_OilyFish_ON_CAD <- mr_heterogeneity(UVdat_OilyFish_ON_CAD)
str(heterogeneity_OilyFish_ON_CAD)
heterogeneity_OilyFish_ON_CAD <- heterogeneity_OilyFish_ON_CAD[,-c(1,2)]


singleSNP_OilyFish_ON_CAD <- mr_singlesnp(UVdat_OilyFish_ON_CAD)
funnel_plot_OilyFish_ON_CAD <- mr_funnel_plot(singleSNP_OilyFish_ON_CAD)
funnel_plot_OilyFish_ON_CAD[[1]]
ggsave(funnel_plot_OilyFish_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_OilyFish_ON_CAD.pdf", width=7, height=7)


loo_OilyFish_ON_CAD <- mr_leaveoneout(UVdat_OilyFish_ON_CAD)
loo_plot_OilyFish_ON_CAD <- mr_leaveoneout_plot(loo_OilyFish_ON_CAD)
loo_plot_OilyFish_ON_CAD[[1]]
ggsave(loo_plot_OilyFish_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_OilyFish_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_OilyFish_ON_CAD <- mr_forest_plot(singleSNP_OilyFish_ON_CAD)
singleSNP_plot_OilyFish_ON_CAD[[1]]
ggsave(singleSNP_plot_OilyFish_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_OilyFish_ON_CAD.pdf", width=7, height=7)

PRESSO_OilyFish_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_OilyFish_ON_CAD, NbDistribution = 10000,  
                                     SignifThreshold = 0.05)

result_OilyFish_ON_CAD <- generate_odds_ratios(result_OilyFish_ON_CAD)
result_OilyFish_ON_CAD <- result_OilyFish_ON_CAD[,-c(1,2)]


#### ________M1.1 Diet_wholegrain ####
str(IV_wholegrain)

out_wholegrain_ON_CAD <- extract_outcome_data(
  snps = IV_wholegrain$SNP,
  outcomes = 'ieu-a-7')
str(out_wholegrain_ON_CAD)
out_wholegrain_ON_CAD$outcome <- "coronary heart disease"

UVdat_wholegrain_ON_CAD <- harmonise_data(
  exposure_dat =  IV_wholegrain, 
  outcome_dat = out_wholegrain_ON_CAD
)
str(UVdat_wholegrain_ON_CAD)
range(UVdat_wholegrain_ON_CAD$F)
range(UVdat_wholegrain_ON_CAD$pval.outcome)

write.xlsx(UVdat_wholegrain_ON_CAD,file = "Output/UVMR_IV/UVdat_wholegrain_ON_CAD.xlsx")

result_wholegrain_ON_CAD <- mr(UVdat_wholegrain_ON_CAD)


scatter_plot_wholegrain_ON_CAD <- mr_scatter_plot(result_wholegrain_ON_CAD, UVdat_wholegrain_ON_CAD)
scatter_plot_wholegrain_ON_CAD[[1]]
ggsave(scatter_plot_wholegrain_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_wholegrain_ON_CAD.pdf", width=7, height=7)


pleiotropy_wholegrain_ON_CAD <- mr_pleiotropy_test(UVdat_wholegrain_ON_CAD)
str(pleiotropy_wholegrain_ON_CAD)
pleiotropy_wholegrain_ON_CAD <- pleiotropy_wholegrain_ON_CAD[,-c(1,2)]


heterogeneity_wholegrain_ON_CAD <- mr_heterogeneity(UVdat_wholegrain_ON_CAD)
str(heterogeneity_wholegrain_ON_CAD)
heterogeneity_wholegrain_ON_CAD <- heterogeneity_wholegrain_ON_CAD[,-c(1,2)]


singleSNP_wholegrain_ON_CAD <- mr_singlesnp(UVdat_wholegrain_ON_CAD)
funnel_plot_wholegrain_ON_CAD <- mr_funnel_plot(singleSNP_wholegrain_ON_CAD)
funnel_plot_wholegrain_ON_CAD[[1]]
ggsave(funnel_plot_wholegrain_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_wholegrain_ON_CAD.pdf", width=7, height=7)


loo_wholegrain_ON_CAD <- mr_leaveoneout(UVdat_wholegrain_ON_CAD)
loo_plot_wholegrain_ON_CAD <- mr_leaveoneout_plot(loo_wholegrain_ON_CAD)
loo_plot_wholegrain_ON_CAD[[1]]
ggsave(loo_plot_wholegrain_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_wholegrain_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_wholegrain_ON_CAD <- mr_forest_plot(singleSNP_wholegrain_ON_CAD)
singleSNP_plot_wholegrain_ON_CAD[[1]]
ggsave(singleSNP_plot_wholegrain_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_wholegrain_ON_CAD.pdf", width=7, height=7)

PRESSO_wholegrain_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_wholegrain_ON_CAD, NbDistribution = 10000,  
                                    SignifThreshold = 0.05)

result_wholegrain_ON_CAD <- generate_odds_ratios(result_wholegrain_ON_CAD)
result_wholegrain_ON_CAD <- result_wholegrain_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceOverallActivity ####
str(IV_DeviceOverallActivity)

out_DeviceOverallActivity_ON_CAD <- extract_outcome_data(
  snps = IV_DeviceOverallActivity$SNP,
  outcomes = 'ieu-a-7')
str(out_DeviceOverallActivity_ON_CAD)
out_DeviceOverallActivity_ON_CAD$outcome <- "coronary heart disease"

UVdat_DeviceOverallActivity_ON_CAD <- harmonise_data(
  exposure_dat =  IV_DeviceOverallActivity, 
  outcome_dat = out_DeviceOverallActivity_ON_CAD
)
str(UVdat_DeviceOverallActivity_ON_CAD)
range(UVdat_DeviceOverallActivity_ON_CAD$F)
range(UVdat_DeviceOverallActivity_ON_CAD$pval.outcome)
write.xlsx(UVdat_DeviceOverallActivity_ON_CAD,file = "Output/UVMR_IV/UVdat_DeviceOverallActivity_ON_CAD.xlsx")

result_DeviceOverallActivity_ON_CAD <- mr(UVdat_DeviceOverallActivity_ON_CAD)


scatter_plot_DeviceOverallActivity_ON_CAD <- mr_scatter_plot(result_DeviceOverallActivity_ON_CAD, UVdat_DeviceOverallActivity_ON_CAD)
scatter_plot_DeviceOverallActivity_ON_CAD[[1]]
ggsave(scatter_plot_DeviceOverallActivity_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DeviceOverallActivity_ON_CAD.pdf", width=7, height=7)


pleiotropy_DeviceOverallActivity_ON_CAD <- mr_pleiotropy_test(UVdat_DeviceOverallActivity_ON_CAD)
str(pleiotropy_DeviceOverallActivity_ON_CAD)
pleiotropy_DeviceOverallActivity_ON_CAD <- pleiotropy_DeviceOverallActivity_ON_CAD[,-c(1,2)]


heterogeneity_DeviceOverallActivity_ON_CAD <- mr_heterogeneity(UVdat_DeviceOverallActivity_ON_CAD)
str(heterogeneity_DeviceOverallActivity_ON_CAD)
heterogeneity_DeviceOverallActivity_ON_CAD <- heterogeneity_DeviceOverallActivity_ON_CAD[,-c(1,2)]


singleSNP_DeviceOverallActivity_ON_CAD <- mr_singlesnp(UVdat_DeviceOverallActivity_ON_CAD)
funnel_plot_DeviceOverallActivity_ON_CAD <- mr_funnel_plot(singleSNP_DeviceOverallActivity_ON_CAD)
funnel_plot_DeviceOverallActivity_ON_CAD[[1]]
ggsave(funnel_plot_DeviceOverallActivity_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DeviceOverallActivity_ON_CAD.pdf", width=7, height=7)


loo_DeviceOverallActivity_ON_CAD <- mr_leaveoneout(UVdat_DeviceOverallActivity_ON_CAD)
loo_plot_DeviceOverallActivity_ON_CAD <- mr_leaveoneout_plot(loo_DeviceOverallActivity_ON_CAD)
loo_plot_DeviceOverallActivity_ON_CAD[[1]]
ggsave(loo_plot_DeviceOverallActivity_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_DeviceOverallActivity_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_DeviceOverallActivity_ON_CAD <- mr_forest_plot(singleSNP_DeviceOverallActivity_ON_CAD)
singleSNP_plot_DeviceOverallActivity_ON_CAD[[1]]
ggsave(singleSNP_plot_DeviceOverallActivity_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DeviceOverallActivity_ON_CAD.pdf", width=7, height=7)

PRESSO_DeviceOverallActivity_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DeviceOverallActivity_ON_CAD, NbDistribution = 10000,  
                              SignifThreshold = 0.05)

result_DeviceOverallActivity_ON_CAD <- generate_odds_ratios(result_DeviceOverallActivity_ON_CAD)
result_DeviceOverallActivity_ON_CAD <- result_DeviceOverallActivity_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_PA ####
str(IV_PA)

out_PA_ON_CAD <- extract_outcome_data(
  snps = IV_PA$SNP,
  outcomes = 'ieu-a-7')
str(out_PA_ON_CAD)
out_PA_ON_CAD$outcome <- "coronary heart disease"

UVdat_PA_ON_CAD <- harmonise_data(
  exposure_dat =  IV_PA, 
  outcome_dat = out_PA_ON_CAD
)
str(UVdat_PA_ON_CAD)
range(UVdat_PA_ON_CAD$F)
range(UVdat_PA_ON_CAD$pval.outcome)
write.xlsx(UVdat_PA_ON_CAD,file = "Output/UVMR_IV/UVdat_PA_ON_CAD.xlsx")

result_PA_ON_CAD <- mr(UVdat_PA_ON_CAD)


scatter_plot_PA_ON_CAD <- mr_scatter_plot(result_PA_ON_CAD, UVdat_PA_ON_CAD)
scatter_plot_PA_ON_CAD[[1]]
ggsave(scatter_plot_PA_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_PA_ON_CAD.pdf", width=7, height=7)


pleiotropy_PA_ON_CAD <- mr_pleiotropy_test(UVdat_PA_ON_CAD)
str(pleiotropy_PA_ON_CAD)
pleiotropy_PA_ON_CAD <- pleiotropy_PA_ON_CAD[,-c(1,2)]


heterogeneity_PA_ON_CAD <- mr_heterogeneity(UVdat_PA_ON_CAD)
str(heterogeneity_PA_ON_CAD)
heterogeneity_PA_ON_CAD <- heterogeneity_PA_ON_CAD[,-c(1,2)]


singleSNP_PA_ON_CAD <- mr_singlesnp(UVdat_PA_ON_CAD)
funnel_plot_PA_ON_CAD <- mr_funnel_plot(singleSNP_PA_ON_CAD)
funnel_plot_PA_ON_CAD[[1]]
ggsave(funnel_plot_PA_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_PA_ON_CAD.pdf", width=7, height=7)


loo_PA_ON_CAD <- mr_leaveoneout(UVdat_PA_ON_CAD)
loo_plot_PA_ON_CAD <- mr_leaveoneout_plot(loo_PA_ON_CAD)
loo_plot_PA_ON_CAD[[1]]
ggsave(loo_plot_PA_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_PA_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_PA_ON_CAD <- mr_forest_plot(singleSNP_PA_ON_CAD)
singleSNP_plot_PA_ON_CAD[[1]]
ggsave(singleSNP_plot_PA_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_PA_ON_CAD.pdf", width=7, height=7)

PRESSO_PA_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_PA_ON_CAD, NbDistribution = 10000,  
                              SignifThreshold = 0.05)

result_PA_ON_CAD <- generate_odds_ratios(result_PA_ON_CAD)
result_PA_ON_CAD <- result_PA_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceModerate ####
str(IV_DeviceModerate)

out_DeviceModerate_ON_CAD <- extract_outcome_data(
  snps = IV_DeviceModerate$SNP,
  outcomes = 'ieu-a-7')
str(out_DeviceModerate_ON_CAD)
out_DeviceModerate_ON_CAD$outcome <- "coronary heart disease"

UVdat_DeviceModerate_ON_CAD <- harmonise_data(
  exposure_dat =  IV_DeviceModerate, 
  outcome_dat = out_DeviceModerate_ON_CAD
)
str(UVdat_DeviceModerate_ON_CAD)
range(UVdat_DeviceModerate_ON_CAD$F)
range(UVdat_DeviceModerate_ON_CAD$pval.outcome)
write.xlsx(UVdat_DeviceModerate_ON_CAD,file = "Output/UVMR_IV/UVdat_DeviceModerate_ON_CAD.xlsx")

result_DeviceModerate_ON_CAD <- mr(UVdat_DeviceModerate_ON_CAD)

result_DeviceModerate_ON_CAD <- generate_odds_ratios(result_DeviceModerate_ON_CAD)
result_DeviceModerate_ON_CAD <- result_DeviceModerate_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_DeviceSedentary ####
str(IV_DeviceSedentary)

out_DeviceSedentary_ON_CAD <- extract_outcome_data(
  snps = IV_DeviceSedentary$SNP,
  outcomes = 'ieu-a-7')
str(out_DeviceSedentary_ON_CAD)
out_DeviceSedentary_ON_CAD$outcome <- "coronary heart disease"

UVdat_DeviceSedentary_ON_CAD <- harmonise_data(
  exposure_dat =  IV_DeviceSedentary, 
  outcome_dat = out_DeviceSedentary_ON_CAD
)
str(UVdat_DeviceSedentary_ON_CAD)
range(UVdat_DeviceSedentary_ON_CAD$F)
range(UVdat_DeviceSedentary_ON_CAD$pval.outcome)
write.xlsx(UVdat_DeviceSedentary_ON_CAD,file = "Output/UVMR_IV/UVdat_DeviceSedentary_ON_CAD.xlsx")

result_DeviceSedentary_ON_CAD <- mr(UVdat_DeviceSedentary_ON_CAD)


scatter_plot_DeviceSedentary_ON_CAD <- mr_scatter_plot(result_DeviceSedentary_ON_CAD, UVdat_DeviceSedentary_ON_CAD)
scatter_plot_DeviceSedentary_ON_CAD[[1]]
ggsave(scatter_plot_DeviceSedentary_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DeviceSedentary_ON_CAD.pdf", width=7, height=7)


pleiotropy_DeviceSedentary_ON_CAD <- mr_pleiotropy_test(UVdat_DeviceSedentary_ON_CAD)
str(pleiotropy_DeviceSedentary_ON_CAD)
pleiotropy_DeviceSedentary_ON_CAD <- pleiotropy_DeviceSedentary_ON_CAD[,-c(1,2)]


heterogeneity_DeviceSedentary_ON_CAD <- mr_heterogeneity(UVdat_DeviceSedentary_ON_CAD)
str(heterogeneity_DeviceSedentary_ON_CAD)
heterogeneity_DeviceSedentary_ON_CAD <- heterogeneity_DeviceSedentary_ON_CAD[,-c(1,2)]


singleSNP_DeviceSedentary_ON_CAD <- mr_singlesnp(UVdat_DeviceSedentary_ON_CAD)
funnel_plot_DeviceSedentary_ON_CAD <- mr_funnel_plot(singleSNP_DeviceSedentary_ON_CAD)
funnel_plot_DeviceSedentary_ON_CAD[[1]]
ggsave(funnel_plot_DeviceSedentary_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DeviceSedentary_ON_CAD.pdf", width=7, height=7)


loo_DeviceSedentary_ON_CAD <- mr_leaveoneout(UVdat_DeviceSedentary_ON_CAD)
loo_plot_DeviceSedentary_ON_CAD <- mr_leaveoneout_plot(loo_DeviceSedentary_ON_CAD)
loo_plot_DeviceSedentary_ON_CAD[[1]]
ggsave(loo_plot_DeviceSedentary_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_DeviceSedentary_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_DeviceSedentary_ON_CAD <- mr_forest_plot(singleSNP_DeviceSedentary_ON_CAD)
singleSNP_plot_DeviceSedentary_ON_CAD[[1]]
ggsave(singleSNP_plot_DeviceSedentary_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DeviceSedentary_ON_CAD.pdf", width=7, height=7)

PRESSO_DeviceSedentary_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DeviceSedentary_ON_CAD, NbDistribution = 10000,  
                              SignifThreshold = 0.05)

result_DeviceSedentary_ON_CAD <- generate_odds_ratios(result_DeviceSedentary_ON_CAD)
result_DeviceSedentary_ON_CAD <- result_DeviceSedentary_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_LST ####
str(IV_LST)

out_LST_ON_CAD <- extract_outcome_data(
  snps = IV_LST$SNP,
  outcomes = 'ieu-a-7')
str(out_LST_ON_CAD)
out_LST_ON_CAD$outcome <- "coronary heart disease"

UVdat_LST_ON_CAD <- harmonise_data(
  exposure_dat =  IV_LST, 
  outcome_dat = out_LST_ON_CAD
)
str(UVdat_LST_ON_CAD)
range(UVdat_LST_ON_CAD$F)
range(UVdat_LST_ON_CAD$pval.outcome)
write.xlsx(UVdat_LST_ON_CAD,file = "Output/UVMR_IV/UVdat_LST_ON_CAD.xlsx")

result_LST_ON_CAD <- mr(UVdat_LST_ON_CAD)


scatter_plot_LST_ON_CAD <- mr_scatter_plot(result_LST_ON_CAD, UVdat_LST_ON_CAD)
scatter_plot_LST_ON_CAD[[1]]
ggsave(scatter_plot_LST_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_LST_ON_CAD.pdf", width=7, height=7)


pleiotropy_LST_ON_CAD <- mr_pleiotropy_test(UVdat_LST_ON_CAD)
str(pleiotropy_LST_ON_CAD)
pleiotropy_LST_ON_CAD <- pleiotropy_LST_ON_CAD[,-c(1,2)]


heterogeneity_LST_ON_CAD <- mr_heterogeneity(UVdat_LST_ON_CAD)
str(heterogeneity_LST_ON_CAD)
heterogeneity_LST_ON_CAD <- heterogeneity_LST_ON_CAD[,-c(1,2)]


singleSNP_LST_ON_CAD <- mr_singlesnp(UVdat_LST_ON_CAD)
funnel_plot_LST_ON_CAD <- mr_funnel_plot(singleSNP_LST_ON_CAD)
funnel_plot_LST_ON_CAD[[1]]
ggsave(funnel_plot_LST_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_LST_ON_CAD.pdf", width=7, height=7)


loo_LST_ON_CAD <- mr_leaveoneout(UVdat_LST_ON_CAD)
loo_plot_LST_ON_CAD <- mr_leaveoneout_plot(loo_LST_ON_CAD)
loo_plot_LST_ON_CAD[[1]]
ggsave(loo_plot_LST_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_LST_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_LST_ON_CAD <- mr_forest_plot(singleSNP_LST_ON_CAD)
singleSNP_plot_LST_ON_CAD[[1]]
ggsave(singleSNP_plot_LST_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_LST_ON_CAD.pdf", width=7, height=7)

PRESSO_LST_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_LST_ON_CAD, NbDistribution = 10000,  
                              SignifThreshold = 0.05)

result_LST_ON_CAD <- generate_odds_ratios(result_LST_ON_CAD)
result_LST_ON_CAD <- result_LST_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_television ####
str(IV_television)

out_television_ON_CAD <- extract_outcome_data(
  snps = IV_television$SNP,
  outcomes = 'ieu-a-7')
str(out_television_ON_CAD)
out_television_ON_CAD$outcome <- "coronary heart disease"

UVdat_television_ON_CAD <- harmonise_data(
  exposure_dat =  IV_television,
  outcome_dat = out_television_ON_CAD
)
str(UVdat_television_ON_CAD)
range(UVdat_television_ON_CAD$F)
range(UVdat_television_ON_CAD$pval.outcome)
write.xlsx(UVdat_television_ON_CAD,file = "Output/UVMR_IV/UVdat_television_ON_CAD.xlsx")

result_television_ON_CAD <- mr(UVdat_television_ON_CAD)


scatter_plot_television_ON_CAD <- mr_scatter_plot(result_television_ON_CAD, UVdat_television_ON_CAD)
scatter_plot_television_ON_CAD[[1]]
ggsave(scatter_plot_television_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_television_ON_CAD.pdf", width=7, height=7)


pleiotropy_television_ON_CAD <- mr_pleiotropy_test(UVdat_television_ON_CAD)
str(pleiotropy_television_ON_CAD)
pleiotropy_television_ON_CAD <- pleiotropy_television_ON_CAD[,-c(1,2)]


heterogeneity_television_ON_CAD <- mr_heterogeneity(UVdat_television_ON_CAD)
str(heterogeneity_television_ON_CAD)
heterogeneity_television_ON_CAD <- heterogeneity_television_ON_CAD[,-c(1,2)]


singleSNP_television_ON_CAD <- mr_singlesnp(UVdat_television_ON_CAD)
funnel_plot_television_ON_CAD <- mr_funnel_plot(singleSNP_television_ON_CAD)
funnel_plot_television_ON_CAD[[1]]
ggsave(funnel_plot_television_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_television_ON_CAD.pdf", width=7, height=7)


loo_television_ON_CAD <- mr_leaveoneout(UVdat_television_ON_CAD)
loo_plot_television_ON_CAD <- mr_leaveoneout_plot(loo_television_ON_CAD)
loo_plot_television_ON_CAD[[1]]
ggsave(loo_plot_television_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_television_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_television_ON_CAD <- mr_forest_plot(singleSNP_television_ON_CAD)
singleSNP_plot_television_ON_CAD[[1]]
ggsave(singleSNP_plot_television_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_television_ON_CAD.pdf", width=7, height=21)

PRESSO_television_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_television_ON_CAD, NbDistribution = 10000,  
                                      SignifThreshold = 0.05)

result_television_ON_CAD <- generate_odds_ratios(result_television_ON_CAD)
result_television_ON_CAD <- result_television_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_computer ####
str(IV_computer)

out_computer_ON_CAD <- extract_outcome_data(
  snps = IV_computer$SNP,
  outcomes = 'ieu-a-7')
str(out_computer_ON_CAD)
out_computer_ON_CAD$outcome <- "coronary heart disease"

UVdat_computer_ON_CAD <- harmonise_data(
  exposure_dat =  IV_computer,
  outcome_dat = out_computer_ON_CAD
)
str(UVdat_computer_ON_CAD)
range(UVdat_computer_ON_CAD$F)
range(UVdat_computer_ON_CAD$pval.outcome)
write.xlsx(UVdat_computer_ON_CAD,file = "Output/UVMR_IV/UVdat_computer_ON_CAD.xlsx")

result_computer_ON_CAD <- mr(UVdat_computer_ON_CAD)


scatter_plot_computer_ON_CAD <- mr_scatter_plot(result_computer_ON_CAD, UVdat_computer_ON_CAD)
scatter_plot_computer_ON_CAD[[1]]
ggsave(scatter_plot_computer_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_computer_ON_CAD.pdf", width=7, height=7)


pleiotropy_computer_ON_CAD <- mr_pleiotropy_test(UVdat_computer_ON_CAD)
str(pleiotropy_computer_ON_CAD)
pleiotropy_computer_ON_CAD <- pleiotropy_computer_ON_CAD[,-c(1,2)]


heterogeneity_computer_ON_CAD <- mr_heterogeneity(UVdat_computer_ON_CAD)
str(heterogeneity_computer_ON_CAD)
heterogeneity_computer_ON_CAD <- heterogeneity_computer_ON_CAD[,-c(1,2)]


singleSNP_computer_ON_CAD <- mr_singlesnp(UVdat_computer_ON_CAD)
funnel_plot_computer_ON_CAD <- mr_funnel_plot(singleSNP_computer_ON_CAD)
funnel_plot_computer_ON_CAD[[1]]
ggsave(funnel_plot_computer_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_computer_ON_CAD.pdf", width=7, height=7)


loo_computer_ON_CAD <- mr_leaveoneout(UVdat_computer_ON_CAD)
loo_plot_computer_ON_CAD <- mr_leaveoneout_plot(loo_computer_ON_CAD)
loo_plot_computer_ON_CAD[[1]]
ggsave(loo_plot_computer_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_computer_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_computer_ON_CAD <- mr_forest_plot(singleSNP_computer_ON_CAD)
singleSNP_plot_computer_ON_CAD[[1]]
ggsave(singleSNP_plot_computer_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_computer_ON_CAD.pdf", width=7, height=21)

PRESSO_computer_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_computer_ON_CAD, NbDistribution = 10000,  
                                      SignifThreshold = 0.05)

result_computer_ON_CAD <- generate_odds_ratios(result_computer_ON_CAD)
result_computer_ON_CAD <- result_computer_ON_CAD[,-c(1,2)]


#### ________M1.2 Physical activity_driving ####
str(IV_driving)

out_driving_ON_CAD <- extract_outcome_data(
  snps = IV_driving$SNP,
  outcomes = 'ieu-a-7')
str(out_driving_ON_CAD)
out_driving_ON_CAD$outcome <- "coronary heart disease"

UVdat_driving_ON_CAD <- harmonise_data(
  exposure_dat =  IV_driving, 
  outcome_dat = out_driving_ON_CAD
)
str(UVdat_driving_ON_CAD)
range(UVdat_driving_ON_CAD$F)
range(UVdat_driving_ON_CAD$pval.outcome)
write.xlsx(UVdat_driving_ON_CAD,file = "Output/UVMR_IV/UVdat_driving_ON_CAD.xlsx")

result_driving_ON_CAD <- mr(UVdat_driving_ON_CAD)


scatter_plot_driving_ON_CAD <- mr_scatter_plot(result_driving_ON_CAD, UVdat_driving_ON_CAD)
scatter_plot_driving_ON_CAD[[1]]
ggsave(scatter_plot_driving_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_driving_ON_CAD.pdf", width=7, height=7)


pleiotropy_driving_ON_CAD <- mr_pleiotropy_test(UVdat_driving_ON_CAD)
str(pleiotropy_driving_ON_CAD)
pleiotropy_driving_ON_CAD <- pleiotropy_driving_ON_CAD[,-c(1,2)]


heterogeneity_driving_ON_CAD <- mr_heterogeneity(UVdat_driving_ON_CAD)
str(heterogeneity_driving_ON_CAD)
heterogeneity_driving_ON_CAD <- heterogeneity_driving_ON_CAD[,-c(1,2)]


singleSNP_driving_ON_CAD <- mr_singlesnp(UVdat_driving_ON_CAD)
funnel_plot_driving_ON_CAD <- mr_funnel_plot(singleSNP_driving_ON_CAD)
funnel_plot_driving_ON_CAD[[1]]
ggsave(funnel_plot_driving_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_driving_ON_CAD.pdf", width=7, height=7)


loo_driving_ON_CAD <- mr_leaveoneout(UVdat_driving_ON_CAD)
loo_plot_driving_ON_CAD <- mr_leaveoneout_plot(loo_driving_ON_CAD)
loo_plot_driving_ON_CAD[[1]]
ggsave(loo_plot_driving_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_driving_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_driving_ON_CAD <- mr_forest_plot(singleSNP_driving_ON_CAD)
singleSNP_plot_driving_ON_CAD[[1]]
ggsave(singleSNP_plot_driving_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_driving_ON_CAD.pdf", width=7, height=7)

PRESSO_driving_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                              OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_driving_ON_CAD, NbDistribution = 10000,  
                              SignifThreshold = 0.05)

result_driving_ON_CAD <- generate_odds_ratios(result_driving_ON_CAD)
result_driving_ON_CAD <- result_driving_ON_CAD[,-c(1,2)]


#### ________M1.3 Smoking_LifetimeSmoking ####
str(IV_LifetimeSmoking)

out_LifetimeSmoking_ON_CAD <- extract_outcome_data(
  snps = IV_LifetimeSmoking$SNP,
  outcomes = 'ieu-a-7')
str(out_LifetimeSmoking_ON_CAD)
out_LifetimeSmoking_ON_CAD$outcome <- "coronary heart disease"

UVdat_LifetimeSmoking_ON_CAD <- harmonise_data(
  exposure_dat =  IV_LifetimeSmoking, 
  outcome_dat = out_LifetimeSmoking_ON_CAD
)
str(UVdat_LifetimeSmoking_ON_CAD)
range(UVdat_LifetimeSmoking_ON_CAD$F)
range(UVdat_LifetimeSmoking_ON_CAD$pval.outcome)
UVdat_LifetimeSmoking_ON_CAD <- UVdat_LifetimeSmoking_ON_CAD[UVdat_LifetimeSmoking_ON_CAD$pval.outcome>5e-6,]
write.xlsx(UVdat_LifetimeSmoking_ON_CAD,file = "Output/UVMR_IV/UVdat_LifetimeSmoking_ON_CAD.xlsx")

set.seed(12345)

result_LifetimeSmoking_ON_CAD <- mr(UVdat_LifetimeSmoking_ON_CAD)


scatter_plot_LifetimeSmoking_ON_CAD <- mr_scatter_plot(result_LifetimeSmoking_ON_CAD, UVdat_LifetimeSmoking_ON_CAD)
scatter_plot_LifetimeSmoking_ON_CAD[[1]]
ggsave(scatter_plot_LifetimeSmoking_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_LifetimeSmoking_ON_CAD.pdf", width=7, height=7)


pleiotropy_LifetimeSmoking_ON_CAD <- mr_pleiotropy_test(UVdat_LifetimeSmoking_ON_CAD)
str(pleiotropy_LifetimeSmoking_ON_CAD)
pleiotropy_LifetimeSmoking_ON_CAD <- pleiotropy_LifetimeSmoking_ON_CAD[,-c(1,2)]


heterogeneity_LifetimeSmoking_ON_CAD <- mr_heterogeneity(UVdat_LifetimeSmoking_ON_CAD)
str(heterogeneity_LifetimeSmoking_ON_CAD)
heterogeneity_LifetimeSmoking_ON_CAD <- heterogeneity_LifetimeSmoking_ON_CAD[,-c(1,2)]


singleSNP_LifetimeSmoking_ON_CAD <- mr_singlesnp(UVdat_LifetimeSmoking_ON_CAD)
funnel_plot_LifetimeSmoking_ON_CAD <- mr_funnel_plot(singleSNP_LifetimeSmoking_ON_CAD)
funnel_plot_LifetimeSmoking_ON_CAD[[1]]
ggsave(funnel_plot_LifetimeSmoking_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_LifetimeSmoking_ON_CAD.pdf", width=7, height=7)


loo_LifetimeSmoking_ON_CAD <- mr_leaveoneout(UVdat_LifetimeSmoking_ON_CAD)
loo_plot_LifetimeSmoking_ON_CAD <- mr_leaveoneout_plot(loo_LifetimeSmoking_ON_CAD)
loo_plot_LifetimeSmoking_ON_CAD[[1]]
ggsave(loo_plot_LifetimeSmoking_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_LifetimeSmoking_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_LifetimeSmoking_ON_CAD <- mr_forest_plot(singleSNP_LifetimeSmoking_ON_CAD)
singleSNP_plot_LifetimeSmoking_ON_CAD[[1]]
ggsave(singleSNP_plot_LifetimeSmoking_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_LifetimeSmoking_ON_CAD.pdf", width=7, height=7)

PRESSO_LifetimeSmoking_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_LifetimeSmoking_ON_CAD, NbDistribution = 10000,  
                                   SignifThreshold = 0.05)

result_LifetimeSmoking_ON_CAD <- generate_odds_ratios(result_LifetimeSmoking_ON_CAD)
result_LifetimeSmoking_ON_CAD <- result_LifetimeSmoking_ON_CAD[,-c(1,2)]


#### ________M1.4 Sleep health_SleepDuration ####
str(IV_SleepDuration)

out_SleepDuration_ON_CAD <- extract_outcome_data(
  snps = IV_SleepDuration$SNP,
  outcomes = 'ieu-a-7')
str(out_SleepDuration_ON_CAD)
out_SleepDuration_ON_CAD$outcome <- "coronary heart disease"

UVdat_SleepDuration_ON_CAD <- harmonise_data(
  exposure_dat =  IV_SleepDuration, 
  outcome_dat = out_SleepDuration_ON_CAD
)
str(UVdat_SleepDuration_ON_CAD)
range(UVdat_SleepDuration_ON_CAD$F)
range(UVdat_SleepDuration_ON_CAD$pval.outcome)
write.xlsx(UVdat_SleepDuration_ON_CAD,file = "Output/UVMR_IV/UVdat_SleepDuration_ON_CAD.xlsx")

result_SleepDuration_ON_CAD <- mr(UVdat_SleepDuration_ON_CAD)


scatter_plot_SleepDuration_ON_CAD <- mr_scatter_plot(result_SleepDuration_ON_CAD, UVdat_SleepDuration_ON_CAD)
scatter_plot_SleepDuration_ON_CAD[[1]]
ggsave(scatter_plot_SleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_SleepDuration_ON_CAD.pdf", width=7, height=7)


pleiotropy_SleepDuration_ON_CAD <- mr_pleiotropy_test(UVdat_SleepDuration_ON_CAD)
str(pleiotropy_SleepDuration_ON_CAD)
pleiotropy_SleepDuration_ON_CAD <- pleiotropy_SleepDuration_ON_CAD[,-c(1,2)]


heterogeneity_SleepDuration_ON_CAD <- mr_heterogeneity(UVdat_SleepDuration_ON_CAD)
str(heterogeneity_SleepDuration_ON_CAD)
heterogeneity_SleepDuration_ON_CAD <- heterogeneity_SleepDuration_ON_CAD[,-c(1,2)]


singleSNP_SleepDuration_ON_CAD <- mr_singlesnp(UVdat_SleepDuration_ON_CAD)
funnel_plot_SleepDuration_ON_CAD <- mr_funnel_plot(singleSNP_SleepDuration_ON_CAD)
funnel_plot_SleepDuration_ON_CAD[[1]]
ggsave(funnel_plot_SleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_SleepDuration_ON_CAD.pdf", width=7, height=7)


loo_SleepDuration_ON_CAD <- mr_leaveoneout(UVdat_SleepDuration_ON_CAD)
loo_plot_SleepDuration_ON_CAD <- mr_leaveoneout_plot(loo_SleepDuration_ON_CAD)
loo_plot_SleepDuration_ON_CAD[[1]]
ggsave(loo_plot_SleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_SleepDuration_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_SleepDuration_ON_CAD <- mr_forest_plot(singleSNP_SleepDuration_ON_CAD)
singleSNP_plot_SleepDuration_ON_CAD[[1]]
ggsave(singleSNP_plot_SleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_SleepDuration_ON_CAD.pdf", width=7, height=7)

PRESSO_SleepDuration_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                 OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_SleepDuration_ON_CAD, NbDistribution = 10000,  
                                 SignifThreshold = 0.05)

result_SleepDuration_ON_CAD <- generate_odds_ratios(result_SleepDuration_ON_CAD)
result_SleepDuration_ON_CAD <- result_SleepDuration_ON_CAD[,-c(1,2)]


#### ________M1.4 Sleep health_DeviceSleepDuration ####
str(IV_DeviceSleepDuration)

out_DeviceSleepDuration_ON_CAD <- extract_outcome_data(
  snps = IV_DeviceSleepDuration$SNP,
  outcomes = 'ieu-a-7')
str(out_DeviceSleepDuration_ON_CAD)
out_DeviceSleepDuration_ON_CAD$outcome <- "coronary heart disease"

UVdat_DeviceSleepDuration_ON_CAD <- harmonise_data(
  exposure_dat =  IV_DeviceSleepDuration, 
  outcome_dat = out_DeviceSleepDuration_ON_CAD
)
str(UVdat_DeviceSleepDuration_ON_CAD)
range(UVdat_DeviceSleepDuration_ON_CAD$F)
range(UVdat_DeviceSleepDuration_ON_CAD$pval.outcome)
write.xlsx(UVdat_DeviceSleepDuration_ON_CAD,file = "Output/UVMR_IV/UVdat_DeviceSleepDuration_ON_CAD.xlsx")

result_DeviceSleepDuration_ON_CAD <- mr(UVdat_DeviceSleepDuration_ON_CAD)


scatter_plot_DeviceSleepDuration_ON_CAD <- mr_scatter_plot(result_DeviceSleepDuration_ON_CAD, UVdat_DeviceSleepDuration_ON_CAD)
scatter_plot_DeviceSleepDuration_ON_CAD[[1]]
ggsave(scatter_plot_DeviceSleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DeviceSleepDuration_ON_CAD.pdf", width=7, height=7)


pleiotropy_DeviceSleepDuration_ON_CAD <- mr_pleiotropy_test(UVdat_DeviceSleepDuration_ON_CAD)
str(pleiotropy_DeviceSleepDuration_ON_CAD)
pleiotropy_DeviceSleepDuration_ON_CAD <- pleiotropy_DeviceSleepDuration_ON_CAD[,-c(1,2)]


heterogeneity_DeviceSleepDuration_ON_CAD <- mr_heterogeneity(UVdat_DeviceSleepDuration_ON_CAD)
str(heterogeneity_DeviceSleepDuration_ON_CAD)
heterogeneity_DeviceSleepDuration_ON_CAD <- heterogeneity_DeviceSleepDuration_ON_CAD[,-c(1,2)]


singleSNP_DeviceSleepDuration_ON_CAD <- mr_singlesnp(UVdat_DeviceSleepDuration_ON_CAD)
funnel_plot_DeviceSleepDuration_ON_CAD <- mr_funnel_plot(singleSNP_DeviceSleepDuration_ON_CAD)
funnel_plot_DeviceSleepDuration_ON_CAD[[1]]
ggsave(funnel_plot_DeviceSleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DeviceSleepDuration_ON_CAD.pdf", width=7, height=7)


loo_DeviceSleepDuration_ON_CAD <- mr_leaveoneout(UVdat_DeviceSleepDuration_ON_CAD)
loo_plot_DeviceSleepDuration_ON_CAD <- mr_leaveoneout_plot(loo_DeviceSleepDuration_ON_CAD)
loo_plot_DeviceSleepDuration_ON_CAD[[1]]
ggsave(loo_plot_DeviceSleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_DeviceSleepDuration_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_DeviceSleepDuration_ON_CAD <- mr_forest_plot(singleSNP_DeviceSleepDuration_ON_CAD)
singleSNP_plot_DeviceSleepDuration_ON_CAD[[1]]
ggsave(singleSNP_plot_DeviceSleepDuration_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DeviceSleepDuration_ON_CAD.pdf", width=7, height=7)

PRESSO_DeviceSleepDuration_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DeviceSleepDuration_ON_CAD, NbDistribution = 10000,  
                                    SignifThreshold = 0.05)

result_DeviceSleepDuration_ON_CAD <- generate_odds_ratios(result_DeviceSleepDuration_ON_CAD)
result_DeviceSleepDuration_ON_CAD <- result_DeviceSleepDuration_ON_CAD[,-c(1,2)]


#### ________M1.4 Sleep health_Insomnia ####
str(IV_Insomnia)

out_Insomnia_ON_CAD <- extract_outcome_data(
  snps = IV_Insomnia$SNP,
  outcomes = 'ieu-a-7')
str(out_Insomnia_ON_CAD)
out_Insomnia_ON_CAD$outcome <- "coronary heart disease"

UVdat_Insomnia_ON_CAD <- harmonise_data(
  exposure_dat =  IV_Insomnia, 
  outcome_dat = out_Insomnia_ON_CAD
)
str(UVdat_Insomnia_ON_CAD)
range(UVdat_Insomnia_ON_CAD$F)
range(UVdat_Insomnia_ON_CAD$pval.outcome)
write.xlsx(UVdat_Insomnia_ON_CAD,file = "Output/UVMR_IV/UVdat_Insomnia_ON_CAD.xlsx")

result_Insomnia_ON_CAD <- mr(UVdat_Insomnia_ON_CAD)


scatter_plot_Insomnia_ON_CAD <- mr_scatter_plot(result_Insomnia_ON_CAD, UVdat_Insomnia_ON_CAD)
scatter_plot_Insomnia_ON_CAD[[1]]
ggsave(scatter_plot_Insomnia_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_Insomnia_ON_CAD.pdf", width=7, height=7)


pleiotropy_Insomnia_ON_CAD <- mr_pleiotropy_test(UVdat_Insomnia_ON_CAD)
str(pleiotropy_Insomnia_ON_CAD)
pleiotropy_Insomnia_ON_CAD <- pleiotropy_Insomnia_ON_CAD[,-c(1,2)]


heterogeneity_Insomnia_ON_CAD <- mr_heterogeneity(UVdat_Insomnia_ON_CAD)
str(heterogeneity_Insomnia_ON_CAD)
heterogeneity_Insomnia_ON_CAD <- heterogeneity_Insomnia_ON_CAD[,-c(1,2)]


singleSNP_Insomnia_ON_CAD <- mr_singlesnp(UVdat_Insomnia_ON_CAD)
funnel_plot_Insomnia_ON_CAD <- mr_funnel_plot(singleSNP_Insomnia_ON_CAD)
funnel_plot_Insomnia_ON_CAD[[1]]
ggsave(funnel_plot_Insomnia_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_Insomnia_ON_CAD.pdf", width=7, height=7)


loo_Insomnia_ON_CAD <- mr_leaveoneout(UVdat_Insomnia_ON_CAD)
loo_plot_Insomnia_ON_CAD <- mr_leaveoneout_plot(loo_Insomnia_ON_CAD)
loo_plot_Insomnia_ON_CAD[[1]]
ggsave(loo_plot_Insomnia_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_Insomnia_ON_CAD.pdf", width=7, height=7)


singleSNP_plot_Insomnia_ON_CAD <- mr_forest_plot(singleSNP_Insomnia_ON_CAD)
singleSNP_plot_Insomnia_ON_CAD[[1]]
ggsave(singleSNP_plot_Insomnia_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_Insomnia_ON_CAD.pdf", width=7, height=7)

PRESSO_Insomnia_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                         OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_Insomnia_ON_CAD, NbDistribution = 10000,  
                                         SignifThreshold = 0.05)

result_Insomnia_ON_CAD <- generate_odds_ratios(result_Insomnia_ON_CAD)
result_Insomnia_ON_CAD <- result_Insomnia_ON_CAD[,-c(1,2)]


#### ________M2.1 Well-being spectrum(WBS) ####
str(IV_WBS)

out_WBS_ON_CAD <- extract_outcome_data(
  snps = IV_WBS$SNP,
  outcomes = 'ieu-a-7')
str(out_WBS_ON_CAD)
out_WBS_ON_CAD$outcome <- "coronary heart disease"

UVdat_WBS_ON_CAD <- harmonise_data(
  exposure_dat =  IV_WBS, 
  outcome_dat = out_WBS_ON_CAD
)
str(UVdat_WBS_ON_CAD)
range(UVdat_WBS_ON_CAD$F)
range(UVdat_WBS_ON_CAD$pval.outcome)
write.xlsx(UVdat_WBS_ON_CAD,file = "Output/UVMR_IV/UVdat_WBS_ON_CAD.xlsx")

result_WBS_ON_CAD <- mr(UVdat_WBS_ON_CAD)


scatter_plot_WBS_ON_CAD <- mr_scatter_plot(result_WBS_ON_CAD, UVdat_WBS_ON_CAD)
scatter_plot_WBS_ON_CAD[[1]]
ggsave(scatter_plot_WBS_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_WBS_ON_CAD.pdf", width=7, height=7)


pleiotropy_WBS_ON_CAD <- mr_pleiotropy_test(UVdat_WBS_ON_CAD)
str(pleiotropy_WBS_ON_CAD)
pleiotropy_WBS_ON_CAD <- pleiotropy_WBS_ON_CAD[,-c(1,2)]


heterogeneity_WBS_ON_CAD <- mr_heterogeneity(UVdat_WBS_ON_CAD)
str(heterogeneity_WBS_ON_CAD)
heterogeneity_WBS_ON_CAD <- heterogeneity_WBS_ON_CAD[,-c(1,2)]


singleSNP_WBS_ON_CAD <- mr_singlesnp(UVdat_WBS_ON_CAD)
funnel_plot_WBS_ON_CAD <- mr_funnel_plot(singleSNP_WBS_ON_CAD)
funnel_plot_WBS_ON_CAD[[1]]
ggsave(funnel_plot_WBS_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_WBS_ON_CAD.pdf", width=7, height=7)


loo_WBS_ON_CAD <- mr_leaveoneout(UVdat_WBS_ON_CAD)
loo_plot_WBS_ON_CAD <- mr_leaveoneout_plot(loo_WBS_ON_CAD)
loo_plot_WBS_ON_CAD[[1]]
ggsave(loo_plot_WBS_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_WBS_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_WBS_ON_CAD <- mr_forest_plot(singleSNP_WBS_ON_CAD)
singleSNP_plot_WBS_ON_CAD[[1]]
ggsave(singleSNP_plot_WBS_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_WBS_ON_CAD.pdf", width=7, height=7)

PRESSO_WBS_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_WBS_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_WBS_ON_CAD <- generate_odds_ratios(result_WBS_ON_CAD)
result_WBS_ON_CAD <- result_WBS_ON_CAD[,-c(1,2)]


#### ________M2.2 MentalProblems ####
str(IV_MentalProblems)

out_MentalProblems_ON_CAD <- extract_outcome_data(
  snps = IV_MentalProblems$SNP,
  outcomes = 'ieu-a-7')
str(out_MentalProblems_ON_CAD)
out_MentalProblems_ON_CAD$outcome <- "coronary heart disease"

UVdat_MentalProblems_ON_CAD <- harmonise_data(
  exposure_dat =  IV_MentalProblems, 
  outcome_dat = out_MentalProblems_ON_CAD
)
str(UVdat_MentalProblems_ON_CAD)
range(UVdat_MentalProblems_ON_CAD$F)
range(UVdat_MentalProblems_ON_CAD$pval.outcome)
write.xlsx(UVdat_MentalProblems_ON_CAD,file = "Output/UVMR_IV/UVdat_MentalProblems_ON_CAD.xlsx")

result_MentalProblems_ON_CAD <- mr(UVdat_MentalProblems_ON_CAD)


scatter_plot_MentalProblems_ON_CAD <- mr_scatter_plot(result_MentalProblems_ON_CAD, UVdat_MentalProblems_ON_CAD)
scatter_plot_MentalProblems_ON_CAD[[1]]
ggsave(scatter_plot_MentalProblems_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_MentalProblems_ON_CAD.pdf", width=7, height=7)


pleiotropy_MentalProblems_ON_CAD <- mr_pleiotropy_test(UVdat_MentalProblems_ON_CAD)
str(pleiotropy_MentalProblems_ON_CAD)
pleiotropy_MentalProblems_ON_CAD <- pleiotropy_MentalProblems_ON_CAD[,-c(1,2)]


heterogeneity_MentalProblems_ON_CAD <- mr_heterogeneity(UVdat_MentalProblems_ON_CAD)
str(heterogeneity_MentalProblems_ON_CAD)
heterogeneity_MentalProblems_ON_CAD <- heterogeneity_MentalProblems_ON_CAD[,-c(1,2)]


singleSNP_MentalProblems_ON_CAD <- mr_singlesnp(UVdat_MentalProblems_ON_CAD)
funnel_plot_MentalProblems_ON_CAD <- mr_funnel_plot(singleSNP_MentalProblems_ON_CAD)
funnel_plot_MentalProblems_ON_CAD[[1]]
ggsave(funnel_plot_MentalProblems_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_MentalProblems_ON_CAD.pdf", width=7, height=7)


loo_MentalProblems_ON_CAD <- mr_leaveoneout(UVdat_MentalProblems_ON_CAD)
loo_plot_MentalProblems_ON_CAD <- mr_leaveoneout_plot(loo_MentalProblems_ON_CAD)
loo_plot_MentalProblems_ON_CAD[[1]]
ggsave(loo_plot_MentalProblems_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_MentalProblems_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_MentalProblems_ON_CAD <- mr_forest_plot(singleSNP_MentalProblems_ON_CAD)
singleSNP_plot_MentalProblems_ON_CAD[[1]]
ggsave(singleSNP_plot_MentalProblems_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_MentalProblems_ON_CAD.pdf", width=7, height=21)

PRESSO_MentalProblems_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_MentalProblems_ON_CAD, NbDistribution = 10000,  
                                      SignifThreshold = 0.05)

result_MentalProblems_ON_CAD <- generate_odds_ratios(result_MentalProblems_ON_CAD)
result_MentalProblems_ON_CAD <- result_MentalProblems_ON_CAD[,-c(1,2)]


#### ________M2.3 Depression ####
str(IV_depression)

out_depression_ON_CAD <- extract_outcome_data(
  snps = IV_depression$SNP,
  outcomes = 'ieu-a-7')
str(out_depression_ON_CAD)
out_depression_ON_CAD$outcome <- "coronary heart disease"

UVdat_depression_ON_CAD <- harmonise_data(
  exposure_dat =  IV_depression, 
  outcome_dat = out_depression_ON_CAD
)
str(UVdat_depression_ON_CAD)
range(UVdat_depression_ON_CAD$F)
range(UVdat_depression_ON_CAD$pval.outcome)
write.xlsx(UVdat_depression_ON_CAD,file = "Output/UVMR_IV/UVdat_depression_ON_CAD.xlsx")

result_depression_ON_CAD <- mr(UVdat_depression_ON_CAD)


scatter_plot_depression_ON_CAD <- mr_scatter_plot(result_depression_ON_CAD, UVdat_depression_ON_CAD)
scatter_plot_depression_ON_CAD[[1]]
ggsave(scatter_plot_depression_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_depression_ON_CAD.pdf", width=7, height=7)


pleiotropy_depression_ON_CAD <- mr_pleiotropy_test(UVdat_depression_ON_CAD)
str(pleiotropy_depression_ON_CAD)
pleiotropy_depression_ON_CAD <- pleiotropy_depression_ON_CAD[,-c(1,2)]


heterogeneity_depression_ON_CAD <- mr_heterogeneity(UVdat_depression_ON_CAD)
str(heterogeneity_depression_ON_CAD)
heterogeneity_depression_ON_CAD <- heterogeneity_depression_ON_CAD[,-c(1,2)]


singleSNP_depression_ON_CAD <- mr_singlesnp(UVdat_depression_ON_CAD)
funnel_plot_depression_ON_CAD <- mr_funnel_plot(singleSNP_depression_ON_CAD)
funnel_plot_depression_ON_CAD[[1]]
ggsave(funnel_plot_depression_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_depression_ON_CAD.pdf", width=7, height=7)


loo_depression_ON_CAD <- mr_leaveoneout(UVdat_depression_ON_CAD)
loo_plot_depression_ON_CAD <- mr_leaveoneout_plot(loo_depression_ON_CAD)
loo_plot_depression_ON_CAD[[1]]
ggsave(loo_plot_depression_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_depression_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_depression_ON_CAD <- mr_forest_plot(singleSNP_depression_ON_CAD)
singleSNP_plot_depression_ON_CAD[[1]]
ggsave(singleSNP_plot_depression_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_depression_ON_CAD.pdf", width=7, height=21)

PRESSO_depression_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                      OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_depression_ON_CAD, NbDistribution = 10000,  
                                      SignifThreshold = 0.05)

result_depression_ON_CAD <- generate_odds_ratios(result_depression_ON_CAD)
result_depression_ON_CAD <- result_depression_ON_CAD[,-c(1,2)]


#### ________M3.1 BMI ####
str(IV_BMI)

out_BMI_ON_CAD <- extract_outcome_data(
  snps = IV_BMI$SNP,
  outcomes = 'ieu-a-7')
str(out_BMI_ON_CAD)
out_BMI_ON_CAD$outcome <- "coronary heart disease"

UVdat_BMI_ON_CAD <- harmonise_data(
  exposure_dat =  IV_BMI, 
  outcome_dat = out_BMI_ON_CAD
)
str(UVdat_BMI_ON_CAD)
range(UVdat_BMI_ON_CAD$F)
range(UVdat_BMI_ON_CAD$pval.outcome)
UVdat_BMI_ON_CAD <- UVdat_BMI_ON_CAD[UVdat_BMI_ON_CAD$pval.outcome>5e-8,]
write.xlsx(UVdat_BMI_ON_CAD,file = "Output/UVMR_IV/UVdat_BMI_ON_CAD.xlsx")

result_BMI_ON_CAD <- mr(UVdat_BMI_ON_CAD)


scatter_plot_BMI_ON_CAD <- mr_scatter_plot(result_BMI_ON_CAD, UVdat_BMI_ON_CAD)
scatter_plot_BMI_ON_CAD[[1]]
ggsave(scatter_plot_BMI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_BMI_ON_CAD.pdf", width=7, height=7)


pleiotropy_BMI_ON_CAD <- mr_pleiotropy_test(UVdat_BMI_ON_CAD)
str(pleiotropy_BMI_ON_CAD)
pleiotropy_BMI_ON_CAD <- pleiotropy_BMI_ON_CAD[,-c(1,2)]


heterogeneity_BMI_ON_CAD <- mr_heterogeneity(UVdat_BMI_ON_CAD)
str(heterogeneity_BMI_ON_CAD)
heterogeneity_BMI_ON_CAD <- heterogeneity_BMI_ON_CAD[,-c(1,2)]


singleSNP_BMI_ON_CAD <- mr_singlesnp(UVdat_BMI_ON_CAD)
funnel_plot_BMI_ON_CAD <- mr_funnel_plot(singleSNP_BMI_ON_CAD)
funnel_plot_BMI_ON_CAD[[1]]
ggsave(funnel_plot_BMI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_BMI_ON_CAD.pdf", width=7, height=7)


loo_BMI_ON_CAD <- mr_leaveoneout(UVdat_BMI_ON_CAD)
loo_plot_BMI_ON_CAD <- mr_leaveoneout_plot(loo_BMI_ON_CAD)
loo_plot_BMI_ON_CAD[[1]]
ggsave(loo_plot_BMI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_BMI_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_BMI_ON_CAD <- mr_forest_plot(singleSNP_BMI_ON_CAD)
singleSNP_plot_BMI_ON_CAD[[1]]
ggsave(singleSNP_plot_BMI_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_BMI_ON_CAD.pdf", width=7, height=21)

PRESSO_BMI_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_BMI_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_BMI_ON_CAD <- generate_odds_ratios(result_BMI_ON_CAD)
result_BMI_ON_CAD <- result_BMI_ON_CAD[,-c(1,2)]


#### ________M3.2 lipid ####
str(IV_lipid)

out_lipid_ON_CAD <- extract_outcome_data(
  snps = IV_lipid$SNP,
  outcomes = 'ieu-a-7')
str(out_lipid_ON_CAD)
out_lipid_ON_CAD$outcome <- "coronary heart disease"

UVdat_lipid_ON_CAD <- harmonise_data(
  exposure_dat =  IV_lipid, 
  outcome_dat = out_lipid_ON_CAD
)
str(UVdat_lipid_ON_CAD)
range(UVdat_lipid_ON_CAD$F)
range(UVdat_lipid_ON_CAD$pval.outcome)
UVdat_lipid_ON_CAD <- UVdat_lipid_ON_CAD[UVdat_lipid_ON_CAD$pval.outcome>5e-8,]
write.xlsx(UVdat_lipid_ON_CAD,file = "Output/UVMR_IV/UVdat_lipid_ON_CAD.xlsx")

result_lipid_ON_CAD <- mr(UVdat_lipid_ON_CAD)


scatter_plot_lipid_ON_CAD <- mr_scatter_plot(result_lipid_ON_CAD, UVdat_lipid_ON_CAD)
scatter_plot_lipid_ON_CAD[[1]]
ggsave(scatter_plot_lipid_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_lipid_ON_CAD.pdf", width=7, height=7)


pleiotropy_lipid_ON_CAD <- mr_pleiotropy_test(UVdat_lipid_ON_CAD)
str(pleiotropy_lipid_ON_CAD)
pleiotropy_lipid_ON_CAD <- pleiotropy_lipid_ON_CAD[,-c(1,2)]


heterogeneity_lipid_ON_CAD <- mr_heterogeneity(UVdat_lipid_ON_CAD)
str(heterogeneity_lipid_ON_CAD)
heterogeneity_lipid_ON_CAD <- heterogeneity_lipid_ON_CAD[,-c(1,2)]


singleSNP_lipid_ON_CAD <- mr_singlesnp(UVdat_lipid_ON_CAD)
funnel_plot_lipid_ON_CAD <- mr_funnel_plot(singleSNP_lipid_ON_CAD)
funnel_plot_lipid_ON_CAD[[1]]
ggsave(funnel_plot_lipid_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_lipid_ON_CAD.pdf", width=7, height=7)


loo_lipid_ON_CAD <- mr_leaveoneout(UVdat_lipid_ON_CAD)
loo_plot_lipid_ON_CAD <- mr_leaveoneout_plot(loo_lipid_ON_CAD)
loo_plot_lipid_ON_CAD[[1]]
ggsave(loo_plot_lipid_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_lipid_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_lipid_ON_CAD <- mr_forest_plot(singleSNP_lipid_ON_CAD)
singleSNP_plot_lipid_ON_CAD[[1]]
ggsave(singleSNP_plot_lipid_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_lipid_ON_CAD.pdf", width=7, height=21)

PRESSO_lipid_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                 OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_lipid_ON_CAD, NbDistribution = 10000,  
                                 SignifThreshold = 0.05)

result_lipid_ON_CAD <- generate_odds_ratios(result_lipid_ON_CAD)
result_lipid_ON_CAD <- result_lipid_ON_CAD[,-c(1,2)]


#### ________M3.3 glucose ####
str(IV_glucose)

out_glucose_ON_CAD <- extract_outcome_data(
  snps = IV_glucose$SNP,
  outcomes = 'ieu-a-7')
str(out_glucose_ON_CAD)
out_glucose_ON_CAD$outcome <- "coronary heart disease"

UVdat_glucose_ON_CAD <- harmonise_data(
  exposure_dat =  IV_glucose, 
  outcome_dat = out_glucose_ON_CAD
)
str(UVdat_glucose_ON_CAD)
range(UVdat_glucose_ON_CAD$F)
range(UVdat_glucose_ON_CAD$pval.outcome)
UVdat_glucose_ON_CAD <- UVdat_glucose_ON_CAD[UVdat_glucose_ON_CAD$pval.outcome>5e-8,]
write.xlsx(UVdat_glucose_ON_CAD,file = "Output/UVMR_IV/UVdat_glucose_ON_CAD.xlsx")

result_glucose_ON_CAD <- mr(UVdat_glucose_ON_CAD)


scatter_plot_glucose_ON_CAD <- mr_scatter_plot(result_glucose_ON_CAD, UVdat_glucose_ON_CAD)
scatter_plot_glucose_ON_CAD[[1]]
ggsave(scatter_plot_glucose_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_glucose_ON_CAD.pdf", width=7, height=7)


pleiotropy_glucose_ON_CAD <- mr_pleiotropy_test(UVdat_glucose_ON_CAD)
str(pleiotropy_glucose_ON_CAD)
pleiotropy_glucose_ON_CAD <- pleiotropy_glucose_ON_CAD[,-c(1,2)]


heterogeneity_glucose_ON_CAD <- mr_heterogeneity(UVdat_glucose_ON_CAD)
str(heterogeneity_glucose_ON_CAD)
heterogeneity_glucose_ON_CAD <- heterogeneity_glucose_ON_CAD[,-c(1,2)]


singleSNP_glucose_ON_CAD <- mr_singlesnp(UVdat_glucose_ON_CAD)
funnel_plot_glucose_ON_CAD <- mr_funnel_plot(singleSNP_glucose_ON_CAD)
funnel_plot_glucose_ON_CAD[[1]]
ggsave(funnel_plot_glucose_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_glucose_ON_CAD.pdf", width=7, height=7)


loo_glucose_ON_CAD <- mr_leaveoneout(UVdat_glucose_ON_CAD)
loo_plot_glucose_ON_CAD <- mr_leaveoneout_plot(loo_glucose_ON_CAD)
loo_plot_glucose_ON_CAD[[1]]
ggsave(loo_plot_glucose_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_glucose_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_glucose_ON_CAD <- mr_forest_plot(singleSNP_glucose_ON_CAD)
singleSNP_plot_glucose_ON_CAD[[1]]
ggsave(singleSNP_plot_glucose_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_glucose_ON_CAD.pdf", width=7, height=21)

PRESSO_glucose_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_glucose_ON_CAD, NbDistribution = 10000,  
                                   SignifThreshold = 0.05)

result_glucose_ON_CAD <- generate_odds_ratios(result_glucose_ON_CAD)
result_glucose_ON_CAD <- result_glucose_ON_CAD[,-c(1,2)]


#### ________M3.4 SBP, DBP ####
## SBP
str(IV_SBP)

out_SBP_ON_CAD <- extract_outcome_data(
  snps = IV_SBP$SNP,
  outcomes = 'ieu-a-7')
str(out_SBP_ON_CAD)
out_SBP_ON_CAD$outcome <- "coronary heart disease"

UVdat_SBP_ON_CAD <- harmonise_data(
  exposure_dat =  IV_SBP, 
  outcome_dat = out_SBP_ON_CAD
)
str(UVdat_SBP_ON_CAD)
range(UVdat_SBP_ON_CAD$F)
range(UVdat_SBP_ON_CAD$pval.outcome)
UVdat_SBP_ON_CAD <- UVdat_SBP_ON_CAD[UVdat_SBP_ON_CAD$pval.outcome>5e-8,]
write.xlsx(UVdat_SBP_ON_CAD,file = "Output/UVMR_IV/UVdat_SBP_ON_CAD.xlsx")

result_SBP_ON_CAD <- mr(UVdat_SBP_ON_CAD)


scatter_plot_SBP_ON_CAD <- mr_scatter_plot(result_SBP_ON_CAD, UVdat_SBP_ON_CAD)
scatter_plot_SBP_ON_CAD[[1]]
ggsave(scatter_plot_SBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_SBP_ON_CAD.pdf", width=7, height=7)


pleiotropy_SBP_ON_CAD <- mr_pleiotropy_test(UVdat_SBP_ON_CAD)
str(pleiotropy_SBP_ON_CAD)
pleiotropy_SBP_ON_CAD <- pleiotropy_SBP_ON_CAD[,-c(1,2)]


heterogeneity_SBP_ON_CAD <- mr_heterogeneity(UVdat_SBP_ON_CAD)
str(heterogeneity_SBP_ON_CAD)
heterogeneity_SBP_ON_CAD <- heterogeneity_SBP_ON_CAD[,-c(1,2)]


singleSNP_SBP_ON_CAD <- mr_singlesnp(UVdat_SBP_ON_CAD)
funnel_plot_SBP_ON_CAD <- mr_funnel_plot(singleSNP_SBP_ON_CAD)
funnel_plot_SBP_ON_CAD[[1]]
ggsave(funnel_plot_SBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_SBP_ON_CAD.pdf", width=7, height=7)


loo_SBP_ON_CAD <- mr_leaveoneout(UVdat_SBP_ON_CAD)
loo_plot_SBP_ON_CAD <- mr_leaveoneout_plot(loo_SBP_ON_CAD)
loo_plot_SBP_ON_CAD[[1]]
ggsave(loo_plot_SBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_SBP_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_SBP_ON_CAD <- mr_forest_plot(singleSNP_SBP_ON_CAD)
singleSNP_plot_SBP_ON_CAD[[1]]
ggsave(singleSNP_plot_SBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_SBP_ON_CAD.pdf", width=7, height=21)

PRESSO_SBP_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_SBP_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_SBP_ON_CAD <- generate_odds_ratios(result_SBP_ON_CAD)
result_SBP_ON_CAD <- result_SBP_ON_CAD[,-c(1,2)]


## DBP
str(IV_DBP)

out_DBP_ON_CAD <- extract_outcome_data(
  snps = IV_DBP$SNP,
  outcomes = 'ieu-a-7')
str(out_DBP_ON_CAD)
out_DBP_ON_CAD$outcome <- "coronary heart disease"

UVdat_DBP_ON_CAD <- harmonise_data(
  exposure_dat =  IV_DBP, 
  outcome_dat = out_DBP_ON_CAD
)
str(UVdat_DBP_ON_CAD)
range(UVdat_DBP_ON_CAD$F)
range(UVdat_DBP_ON_CAD$pval.outcome)
UVdat_DBP_ON_CAD <- UVdat_DBP_ON_CAD[UVdat_DBP_ON_CAD$pval.outcome>5e-8,]
write.xlsx(UVdat_DBP_ON_CAD,file = "Output/UVMR_IV/UVdat_DBP_ON_CAD.xlsx")

result_DBP_ON_CAD <- mr(UVdat_DBP_ON_CAD)


scatter_plot_DBP_ON_CAD <- mr_scatter_plot(result_DBP_ON_CAD, UVdat_DBP_ON_CAD)
scatter_plot_DBP_ON_CAD[[1]]
ggsave(scatter_plot_DBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/scatter_plot_DBP_ON_CAD.pdf", width=7, height=7)


pleiotropy_DBP_ON_CAD <- mr_pleiotropy_test(UVdat_DBP_ON_CAD)
str(pleiotropy_DBP_ON_CAD)
pleiotropy_DBP_ON_CAD <- pleiotropy_DBP_ON_CAD[,-c(1,2)]


heterogeneity_DBP_ON_CAD <- mr_heterogeneity(UVdat_DBP_ON_CAD)
str(heterogeneity_DBP_ON_CAD)
heterogeneity_DBP_ON_CAD <- heterogeneity_DBP_ON_CAD[,-c(1,2)]


singleSNP_DBP_ON_CAD <- mr_singlesnp(UVdat_DBP_ON_CAD)
funnel_plot_DBP_ON_CAD <- mr_funnel_plot(singleSNP_DBP_ON_CAD)
funnel_plot_DBP_ON_CAD[[1]]
ggsave(funnel_plot_DBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/funnel_plot_DBP_ON_CAD.pdf", width=7, height=7)


loo_DBP_ON_CAD <- mr_leaveoneout(UVdat_DBP_ON_CAD)
loo_plot_DBP_ON_CAD <- mr_leaveoneout_plot(loo_DBP_ON_CAD)
loo_plot_DBP_ON_CAD[[1]]
ggsave(loo_plot_DBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/loo_plot_DBP_ON_CAD.pdf", width=7, height=21)


singleSNP_plot_DBP_ON_CAD <- mr_forest_plot(singleSNP_DBP_ON_CAD)
singleSNP_plot_DBP_ON_CAD[[1]]
ggsave(singleSNP_plot_DBP_ON_CAD[[1]], file="Output/UVMR_Secondary Results/singleSNP_plot_DBP_ON_CAD.pdf", width=7, height=21)

PRESSO_DBP_ON_CAD <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                               OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = UVdat_DBP_ON_CAD, NbDistribution = 10000,  
                               SignifThreshold = 0.05)

result_DBP_ON_CAD <- generate_odds_ratios(result_DBP_ON_CAD)
result_DBP_ON_CAD <- result_DBP_ON_CAD[,-c(1,2)]


###### Main Results
SummaryResultsUVMR_Mediator_ON_CAD <- rbind(result_protein_ON_CAD,
                                            result_fat_ON_CAD,
                                            result_FreshFruit_ON_CAD,
                                            result_DriedFruit_ON_CAD,
                                            result_RawVegetables_ON_CAD,
                                            result_salt_ON_CAD,
                                            result_processed_ON_CAD,
                                            result_OilyFish_ON_CAD,
                                            result_wholegrain_ON_CAD,
                                            result_DeviceOverallActivity_ON_CAD,
                                            result_PA_ON_CAD,
                                            result_DeviceModerate_ON_CAD,
                                            result_DeviceSedentary_ON_CAD,
                                            result_LST_ON_CAD,
                                            result_television_ON_CAD,
                                            result_computer_ON_CAD,
                                            result_driving_ON_CAD,
                                            result_LifetimeSmoking_ON_CAD,
                                            result_SleepDuration_ON_CAD,
                                            result_DeviceSleepDuration_ON_CAD,
                                            result_Insomnia_ON_CAD,
                                            result_WBS_ON_CAD,
                                            result_MentalProblems_ON_CAD,
                                            result_depression_ON_CAD,
                                            result_BMI_ON_CAD,
                                            result_lipid_ON_CAD,
                                            result_glucose_ON_CAD,
                                            result_SBP_ON_CAD,
                                            result_DBP_ON_CAD)
write.xlsx(SummaryResultsUVMR_Mediator_ON_CAD, file = "Output/UVMR_Main Results/SummaryResultsUVMR_Mediator_ON_CAD.xlsx")


###### Secondary Results

Summary_UVPleiotropy_Mediator_ON_CAD <- rbind(pleiotropy_protein_ON_CAD,
                                              pleiotropy_fat_ON_CAD,
                                              pleiotropy_FreshFruit_ON_CAD,
                                              pleiotropy_DriedFruit_ON_CAD,
                                              pleiotropy_RawVegetables_ON_CAD,
                                              pleiotropy_salt_ON_CAD,
                                              pleiotropy_processed_ON_CAD,
                                              pleiotropy_OilyFish_ON_CAD,
                                              pleiotropy_wholegrain_ON_CAD,
                                              pleiotropy_DeviceOverallActivity_ON_CAD,
                                              pleiotropy_PA_ON_CAD,
                                              pleiotropy_DeviceSedentary_ON_CAD,
                                              pleiotropy_LST_ON_CAD,
                                              pleiotropy_television_ON_CAD,
                                              pleiotropy_computer_ON_CAD,
                                              pleiotropy_driving_ON_CAD,
                                              pleiotropy_LifetimeSmoking_ON_CAD,
                                              pleiotropy_SleepDuration_ON_CAD,
                                              pleiotropy_DeviceSleepDuration_ON_CAD,
                                              pleiotropy_Insomnia_ON_CAD,
                                              pleiotropy_WBS_ON_CAD,
                                              pleiotropy_MentalProblems_ON_CAD,
                                              pleiotropy_depression_ON_CAD,
                                              pleiotropy_BMI_ON_CAD,
                                              pleiotropy_lipid_ON_CAD,
                                              pleiotropy_glucose_ON_CAD,
                                              pleiotropy_SBP_ON_CAD,
                                              pleiotropy_DBP_ON_CAD)
write.xlsx(Summary_UVPleiotropy_Mediator_ON_CAD,"Output/UVMR_Secondary Results/Summary_UVPleiotropy_Mediator_ON_CAD.xlsx")


Summary_UVHeterogeneity_Mediator_ON_CAD <- rbind(heterogeneity_protein_ON_CAD,
                                                 heterogeneity_fat_ON_CAD,
                                                 heterogeneity_FreshFruit_ON_CAD,
                                                 heterogeneity_DriedFruit_ON_CAD,
                                                 heterogeneity_RawVegetables_ON_CAD,
                                                 heterogeneity_salt_ON_CAD,
                                                 heterogeneity_processed_ON_CAD,
                                                 heterogeneity_OilyFish_ON_CAD,
                                                 heterogeneity_wholegrain_ON_CAD,
                                                 heterogeneity_DeviceOverallActivity_ON_CAD,
                                                 heterogeneity_PA_ON_CAD,
                                                 heterogeneity_DeviceSedentary_ON_CAD,
                                                 heterogeneity_LST_ON_CAD,
                                                 heterogeneity_television_ON_CAD,
                                                 heterogeneity_computer_ON_CAD,
                                                 heterogeneity_driving_ON_CAD,
                                                 heterogeneity_LifetimeSmoking_ON_CAD,
                                                 heterogeneity_SleepDuration_ON_CAD,
                                                 heterogeneity_DeviceSleepDuration_ON_CAD,
                                                 heterogeneity_Insomnia_ON_CAD,
                                                 heterogeneity_WBS_ON_CAD,
                                                 heterogeneity_MentalProblems_ON_CAD,
                                                 heterogeneity_depression_ON_CAD,
                                                 heterogeneity_BMI_ON_CAD,
                                                 heterogeneity_lipid_ON_CAD,
                                                 heterogeneity_glucose_ON_CAD,
                                                 heterogeneity_SBP_ON_CAD,
                                                 heterogeneity_DBP_ON_CAD)
write.xlsx(Summary_UVHeterogeneity_Mediator_ON_CAD,"Output/UVMR_Secondary Results/Summary_UVHeterogeneity_Mediator_ON_CAD.xlsx")


## the mediator should have a direct causal effect on the outcome independently of education

#### ________M1.2 Physical activity_PA ####
str(IV_educationStringent)
str(IV_PA)

SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_PA <- IV_PA["SNP"]

SNP_IV_education.PA <- rbind(SNP_IV_educationStringent,SNP_IV_PA)# 364SNP
SNP_IV_education.PA <- unique(SNP_IV_education.PA) # 0SNP
SNP_IV_education.PA <- clump_data(SNP_IV_education.PA) # 348SNP
write.xlsx(SNP_IV_education.PA, file = "Output/MVMR_IV/SNP_IV_education.PA.xlsx")

exp_education_MVwith.PA <- format_data(
  GWAS_education,
  snps = SNP_IV_education.PA$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.PA$exposure <- "education"

exp_PA_MVwith.education <- format_data(
  GWAS_PA,
  snps = SNP_IV_education.PA$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_PA_MVwith.education$outcome <- "PA"

education.PA <- harmonise_data(exposure_dat = exp_education_MVwith.PA, outcome_dat = exp_PA_MVwith.education, action = 1)

out_CAD_MV.education.PA <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.PA$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.PA$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.PA, outcome_dat = out_CAD_MV.education.PA, action = 1)

education.PA.dat <- education.PA[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.PA.dat) <- c("SNP", "beta.education", "beta.PA", "se.education", "se.PA")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.PA_ON_CAD <- merge(education.PA.dat, education.CAD.dat, by = "SNP") # 339SNP
write.xlsx(MVdat_education.PA_ON_CAD, file = "Output/MVMR_IV/MVdat_education.PA_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.PA_ON_CAD[,c("beta.education","beta.PA")])
bxse = as.matrix(MVdat_education.PA_ON_CAD[,c("se.education","se.PA")])
by = as.vector(MVdat_education.PA_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.PA_ON_CAD$se.CAD)

MVdatForm_education.PA_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","PA"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.PA_ON_CAD <- mr_mvivw(MVdatForm_education.PA_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.PA_ON_CAD <- mr_mvegger(MVdatForm_education.PA_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.PA_ON_CAD <- mr_mvlasso(MVdatForm_education.PA_ON_CAD)

# MV Median
result_MV.Median_education.PA_ON_CAD <- mr_mvmedian(MVdatForm_education.PA_ON_CAD)



F_education.PA <- strength_mvmr(r_input = MVdatForm_education.PA_ON_CAD, gencov = 0)


mv_hete_education.PA <- pleiotropy_mvmr(r_input = MVdatForm_education.PA_ON_CAD, gencov = 0)


PRESSO_education.PA_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                        BetaExposure = c("beta.education", "beta.PA"), 
                                        SdOutcome = "se.CAD", 
                                        SdExposure = c("se.education", "se.PA"),
                                        OUTLIERtest = TRUE, 
                                        DISTORTIONtest = TRUE, 
                                        data = MVdat_education.PA_ON_CAD,
                                        NbDistribution = 1000, 
                                        SignifThreshold = 0.05)


#### ________M1.2 Physical activity_LST ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_LST <- IV_LST["SNP"]

SNP_IV_education.LST <- rbind(SNP_IV_educationStringent,SNP_IV_LST)# 463SNP
SNP_IV_education.LST <- unique(SNP_IV_education.LST) # 0SNP
SNP_IV_education.LST <- clump_data(SNP_IV_education.LST) # 355SNP
write.xlsx(SNP_IV_education.LST, file = "Output/MVMR_IV/SNP_IV_education.LST.xlsx")

exp_education_MVwith.LST <- format_data(
  GWAS_education,
  snps = SNP_IV_education.LST$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.LST$exposure <- "education"

exp_LST_MVwith.education <- format_data(
  GWAS_LST,
  snps = SNP_IV_education.LST$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_LST_MVwith.education$outcome <- "LST"

education.LST <- harmonise_data(exposure_dat = exp_education_MVwith.LST, outcome_dat = exp_LST_MVwith.education, action = 1)

out_CAD_MV.education.LST <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.LST$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.LST$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.LST, outcome_dat = out_CAD_MV.education.LST, action = 1)

education.LST.dat <- education.LST[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.LST.dat) <- c("SNP", "beta.education", "beta.LST", "se.education", "se.LST")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.LST_ON_CAD <- merge(education.LST.dat, education.CAD.dat, by = "SNP") # 342SNP
write.xlsx(MVdat_education.LST_ON_CAD, file = "Output/MVMR_IV/MVdat_education.LST_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.LST_ON_CAD[,c("beta.education","beta.LST")])
bxse = as.matrix(MVdat_education.LST_ON_CAD[,c("se.education","se.LST")])
by = as.vector(MVdat_education.LST_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.LST_ON_CAD$se.CAD)

MVdatForm_education.LST_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","LST"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.LST_ON_CAD <- mr_mvivw(MVdatForm_education.LST_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.LST_ON_CAD <- mr_mvegger(MVdatForm_education.LST_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.LST_ON_CAD <- mr_mvlasso(MVdatForm_education.LST_ON_CAD)

# MV Median
result_MV.Median_education.LST_ON_CAD <- mr_mvmedian(MVdatForm_education.LST_ON_CAD)



F_education.LST <- strength_mvmr(r_input = MVdatForm_education.LST_ON_CAD, gencov = 0)


mv_hete_education.LST <- pleiotropy_mvmr(r_input = MVdatForm_education.LST_ON_CAD, gencov = 0)


PRESSO_education.LST_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                        BetaExposure = c("beta.education", "beta.LST"), 
                                        SdOutcome = "se.CAD", 
                                        SdExposure = c("se.education", "se.LST"),
                                        OUTLIERtest = TRUE, 
                                        DISTORTIONtest = TRUE, 
                                        data = MVdat_education.LST_ON_CAD,
                                        NbDistribution = 1000, 
                                        SignifThreshold = 0.05)


#### ________M1.2 Physical activity_television ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_television <- IV_television["SNP"]

SNP_IV_education.television <- rbind(SNP_IV_educationStringent,SNP_IV_television)# 509SNP
SNP_IV_education.television <- unique(SNP_IV_education.television) # 1SNP
tt <- clump_data(SNP_IV_education.television)
SNP_IV_education.television <- clump_data(SNP_IV_education.television) # 376SNP
write.xlsx(SNP_IV_education.television, file = "Output/MVMR_IV/SNP_IV_education.television.xlsx")

exp_education_MVwith.television <- format_data(
  GWAS_education,
  snps = SNP_IV_education.television$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.television$exposure <- "education"

exp_television_MVwith.education <- format_data(
  GWAS_television,
  snps = SNP_IV_education.television$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_television_MVwith.education$outcome <- "television"

education.television <- harmonise_data(exposure_dat = exp_education_MVwith.television, outcome_dat = exp_television_MVwith.education, action = 1)

out_CAD_MV.education.television <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.television$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.television$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.television, outcome_dat = out_CAD_MV.education.television, action = 1)

education.television.dat <- education.television[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.television.dat) <- c("SNP", "beta.education", "beta.television", "se.education", "se.television")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.television_ON_CAD <- merge(education.television.dat, education.CAD.dat, by = "SNP") # 361SNP
write.xlsx(MVdat_education.television_ON_CAD, file = "Output/MVMR_IV/MVdat_education.television_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.television_ON_CAD[,c("beta.education","beta.television")])
bxse = as.matrix(MVdat_education.television_ON_CAD[,c("se.education","se.television")])
by = as.vector(MVdat_education.television_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.television_ON_CAD$se.CAD)

MVdatForm_education.television_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","television"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.television_ON_CAD <- mr_mvivw(MVdatForm_education.television_ON_CAD)
tt <- mvmr(MVdatForm_education.television_ON_CAD)#F-statistic: 42.29 on 2 and 359 DF

# MV MR-Egger
result_MV.Egger_education.television_ON_CAD <- mr_mvegger(MVdatForm_education.television_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.television_ON_CAD <- mr_mvlasso(MVdatForm_education.television_ON_CAD)

# MV Median
result_MV.Median_education.television_ON_CAD <- mr_mvmedian(MVdatForm_education.television_ON_CAD)



F_education.television <- strength_mvmr(r_input = MVdatForm_education.television_ON_CAD, gencov = 0)


mv_hete_education.television <- pleiotropy_mvmr(r_input = MVdatForm_education.television_ON_CAD, gencov = 0)


PRESSO_education.television_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.education", "beta.television"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.education", "se.television"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_education.television_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)

#### ________M1.2 Physical activity_computer ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_computer <- IV_computer["SNP"]

SNP_IV_education.computer <- rbind(SNP_IV_educationStringent,SNP_IV_computer)# 397SNP
SNP_IV_education.computer <- unique(SNP_IV_education.computer) # 2SNP
SNP_IV_education.computer <- clump_data(SNP_IV_education.computer) # 347SNP
write.xlsx(SNP_IV_education.computer, file = "Output/MVMR_IV/SNP_IV_education.computer.xlsx")

exp_education_MVwith.computer <- format_data(
  GWAS_education,
  snps = SNP_IV_education.computer$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.computer$exposure <- "education"

exp_computer_MVwith.education <- format_data(
  GWAS_computer,
  snps = SNP_IV_education.computer$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF") 
exp_computer_MVwith.education$outcome <- "computer"

education.computer <- harmonise_data(exposure_dat = exp_education_MVwith.computer, outcome_dat = exp_computer_MVwith.education, action = 1)

out_CAD_MV.education.computer <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.computer$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.computer$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.computer, outcome_dat = out_CAD_MV.education.computer, action = 1)

education.computer.dat <- education.computer[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.computer.dat) <- c("SNP", "beta.education", "beta.computer", "se.education", "se.computer")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.computer_ON_CAD <- merge(education.computer.dat, education.CAD.dat, by = "SNP") # 345SNP
write.xlsx(MVdat_education.computer_ON_CAD, file = "Output/MVMR_IV/MVdat_education.computer_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.computer_ON_CAD[,c("beta.education","beta.computer")])
bxse = as.matrix(MVdat_education.computer_ON_CAD[,c("se.education","se.computer")])
by = as.vector(MVdat_education.computer_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.computer_ON_CAD$se.CAD)

MVdatForm_education.computer_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","computer"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.computer_ON_CAD <- mr_mvivw(MVdatForm_education.computer_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.computer_ON_CAD <- mr_mvegger(MVdatForm_education.computer_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.computer_ON_CAD <- mr_mvlasso(MVdatForm_education.computer_ON_CAD)

# MV Median
result_MV.Median_education.computer_ON_CAD <- mr_mvmedian(MVdatForm_education.computer_ON_CAD)



F_education.computer <- strength_mvmr(r_input = MVdatForm_education.computer_ON_CAD, gencov = 0)


mv_hete_education.computer <- pleiotropy_mvmr(r_input = MVdatForm_education.computer_ON_CAD, gencov = 0)


PRESSO_education.computer_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                             BetaExposure = c("beta.education", "beta.computer"), 
                                             SdOutcome = "se.CAD", 
                                             SdExposure = c("se.education", "se.computer"),
                                             OUTLIERtest = TRUE, 
                                             DISTORTIONtest = TRUE, 
                                             data = MVdat_education.computer_ON_CAD,
                                             NbDistribution = 1000, 
                                             SignifThreshold = 0.05)



#### ________M1.2 Physical activity_driving ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_driving <- IV_driving["SNP"]

SNP_IV_education.driving <- rbind(SNP_IV_educationStringent,SNP_IV_driving)# 353SNP
SNP_IV_education.driving <- unique(SNP_IV_education.driving) # 0SNP
SNP_IV_education.driving <- clump_data(SNP_IV_education.driving) # 346SNP
write.xlsx(SNP_IV_education.driving, file = "Output/MVMR_IV/SNP_IV_education.driving.xlsx")

exp_education_MVwith.driving <- format_data(
  GWAS_education,
  snps = SNP_IV_education.driving$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.driving$exposure <- "education"

exp_driving_MVwith.education <- format_data(
  GWAS_driving,
  snps = SNP_IV_education.driving$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF") 
exp_driving_MVwith.education$outcome <- "driving"

education.driving <- harmonise_data(exposure_dat = exp_education_MVwith.driving, outcome_dat = exp_driving_MVwith.education, action = 1)

out_CAD_MV.education.driving <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.driving$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.driving$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.driving, outcome_dat = out_CAD_MV.education.driving, action = 1)

education.driving.dat <- education.driving[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.driving.dat) <- c("SNP", "beta.education", "beta.driving", "se.education", "se.driving")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.driving_ON_CAD <- merge(education.driving.dat, education.CAD.dat, by = "SNP") # 344SNP
write.xlsx(MVdat_education.driving_ON_CAD, file = "Output/MVMR_IV/MVdat_education.driving_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.driving_ON_CAD[,c("beta.education","beta.driving")])
bxse = as.matrix(MVdat_education.driving_ON_CAD[,c("se.education","se.driving")])
by = as.vector(MVdat_education.driving_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.driving_ON_CAD$se.CAD)

MVdatForm_education.driving_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","driving"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.driving_ON_CAD <- mr_mvivw(MVdatForm_education.driving_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.driving_ON_CAD <- mr_mvegger(MVdatForm_education.driving_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.driving_ON_CAD <- mr_mvlasso(MVdatForm_education.driving_ON_CAD)

# MV Median
result_MV.Median_education.driving_ON_CAD <- mr_mvmedian(MVdatForm_education.driving_ON_CAD)



F_education.driving <- strength_mvmr(r_input = MVdatForm_education.driving_ON_CAD, gencov = 0)


mv_hete_education.driving <- pleiotropy_mvmr(r_input = MVdatForm_education.driving_ON_CAD, gencov = 0)


PRESSO_education.driving_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.education", "beta.driving"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.education", "se.driving"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_education.driving_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)



#### ________M1.3 Smoking_LifetimeSmoking ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_LifetimeSmoking <- IV_LifetimeSmoking["SNP"]

SNP_IV_education.LifetimeSmoking <- rbind(SNP_IV_educationStringent,SNP_IV_LifetimeSmoking)# 474SNP
SNP_IV_education.LifetimeSmoking <- unique(SNP_IV_education.LifetimeSmoking) # 0SNP
SNP_IV_education.LifetimeSmoking <- clump_data(SNP_IV_education.LifetimeSmoking) # 366SNP
write.xlsx(SNP_IV_education.LifetimeSmoking, file = "Output/MVMR_IV/SNP_IV_education.LifetimeSmoking.xlsx")

exp_education_MVwith.LifetimeSmoking <- format_data(
  GWAS_education,
  snps = SNP_IV_education.LifetimeSmoking$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.LifetimeSmoking$exposure <- "education"

exp_LifetimeSmoking_MVwith.education <- format_data(
  GWAS_LifetimeSmoking,
  snps = SNP_IV_education.LifetimeSmoking$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "P") 
exp_LifetimeSmoking_MVwith.education$outcome <- "LifetimeSmoking"

education.LifetimeSmoking <- harmonise_data(exposure_dat = exp_education_MVwith.LifetimeSmoking, outcome_dat = exp_LifetimeSmoking_MVwith.education, action = 1)

out_CAD_MV.education.LifetimeSmoking <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.LifetimeSmoking$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.LifetimeSmoking$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.LifetimeSmoking, outcome_dat = out_CAD_MV.education.LifetimeSmoking, action = 1)

education.LifetimeSmoking.dat <- education.LifetimeSmoking[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.LifetimeSmoking.dat) <- c("SNP", "beta.education", "beta.LifetimeSmoking", "se.education", "se.LifetimeSmoking")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.LifetimeSmoking_ON_CAD <- merge(education.LifetimeSmoking.dat, education.CAD.dat, by = "SNP") # 365SNP
write.xlsx(MVdat_education.LifetimeSmoking_ON_CAD, file = "Output/MVMR_IV/MVdat_education.LifetimeSmoking_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.LifetimeSmoking_ON_CAD[,c("beta.education","beta.LifetimeSmoking")])
bxse = as.matrix(MVdat_education.LifetimeSmoking_ON_CAD[,c("se.education","se.LifetimeSmoking")])
by = as.vector(MVdat_education.LifetimeSmoking_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.LifetimeSmoking_ON_CAD$se.CAD)

MVdatForm_education.LifetimeSmoking_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","LifetimeSmoking"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.LifetimeSmoking_ON_CAD <- mr_mvivw(MVdatForm_education.LifetimeSmoking_ON_CAD)
tt <- mvmr(MVdatForm_education.LifetimeSmoking_ON_CAD)#F-statistic: 48.97 on 2 and 363 DF


# MV MR-Egger
result_MV.Egger_education.LifetimeSmoking_ON_CAD <- mr_mvegger(MVdatForm_education.LifetimeSmoking_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.LifetimeSmoking_ON_CAD <- mr_mvlasso(MVdatForm_education.LifetimeSmoking_ON_CAD)

# MV Median
result_MV.Median_education.LifetimeSmoking_ON_CAD <- mr_mvmedian(MVdatForm_education.LifetimeSmoking_ON_CAD)



F_education.LifetimeSmoking <- strength_mvmr(r_input = MVdatForm_education.LifetimeSmoking_ON_CAD, gencov = 0)


mv_hete_education.LifetimeSmoking <- pleiotropy_mvmr(r_input = MVdatForm_education.LifetimeSmoking_ON_CAD, gencov = 0)


PRESSO_education.LifetimeSmoking_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                         BetaExposure = c("beta.education", "beta.LifetimeSmoking"), 
                                         SdOutcome = "se.CAD", 
                                         SdExposure = c("se.education", "se.LifetimeSmoking"),
                                         OUTLIERtest = TRUE, 
                                         DISTORTIONtest = TRUE, 
                                         data = MVdat_education.LifetimeSmoking_ON_CAD,
                                         NbDistribution = 1000, 
                                         SignifThreshold = 0.05)


#### ________M1.4 Sleep health_SleepDuration ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_SleepDuration <- IV_SleepDuration["SNP"]

SNP_IV_education.SleepDuration <- rbind(SNP_IV_educationStringent,SNP_IV_SleepDuration) # 395SNP
SNP_IV_education.SleepDuration <- unique(SNP_IV_education.SleepDuration) # 0SNP
SNP_IV_education.SleepDuration <- clump_data(SNP_IV_education.SleepDuration) # 349SNP
write.xlsx(SNP_IV_education.SleepDuration, file = "Output/MVMR_IV/SNP_IV_education.SleepDuration.xlsx")

exp_education_MVwith.SleepDuration <- format_data(
  GWAS_education,
  snps = SNP_IV_education.SleepDuration$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.SleepDuration$exposure <- "education"

exp_SleepDuration_MVwith.education <- format_data(
  GWAS_SleepDuration,
  snps = SNP_IV_education.SleepDuration$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_SleepDuration_MVwith.education$outcome <- "SleepDuration"

education.SleepDuration <- harmonise_data(exposure_dat = exp_education_MVwith.SleepDuration, outcome_dat = exp_SleepDuration_MVwith.education, action = 1)

out_CAD_MV.education.SleepDuration <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.SleepDuration$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.SleepDuration$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.SleepDuration, outcome_dat = out_CAD_MV.education.SleepDuration, action = 1)

education.SleepDuration.dat <- education.SleepDuration[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.SleepDuration.dat) <- c("SNP", "beta.education", "beta.SleepDuration", "se.education", "se.SleepDuration")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.SleepDuration_ON_CAD <- merge(education.SleepDuration.dat, education.CAD.dat, by = "SNP") # 332SNP
write.xlsx(MVdat_education.SleepDuration_ON_CAD, file = "Output/MVMR_IV/MVdat_education.SleepDuration_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.SleepDuration_ON_CAD[,c("beta.education","beta.SleepDuration")])
bxse = as.matrix(MVdat_education.SleepDuration_ON_CAD[,c("se.education","se.SleepDuration")])
by = as.vector(MVdat_education.SleepDuration_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.SleepDuration_ON_CAD$se.CAD)

MVdatForm_education.SleepDuration_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","SleepDuration"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.SleepDuration_ON_CAD <- mr_mvivw(MVdatForm_education.SleepDuration_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.SleepDuration_ON_CAD <- mr_mvegger(MVdatForm_education.SleepDuration_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.SleepDuration_ON_CAD <- mr_mvlasso(MVdatForm_education.SleepDuration_ON_CAD)

# MV Median
result_MV.Median_education.SleepDuration_ON_CAD <- mr_mvmedian(MVdatForm_education.SleepDuration_ON_CAD)



F_education.SleepDuration <- strength_mvmr(r_input = MVdatForm_education.SleepDuration_ON_CAD, gencov = 0)


mv_hete_education.SleepDuration <- pleiotropy_mvmr(r_input = MVdatForm_education.SleepDuration_ON_CAD, gencov = 0)


PRESSO_education.SleepDuration_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                           BetaExposure = c("beta.education", "beta.SleepDuration"), 
                                           SdOutcome = "se.CAD", 
                                           SdExposure = c("se.education", "se.SleepDuration"),
                                           OUTLIERtest = TRUE, 
                                           DISTORTIONtest = TRUE, 
                                           data = MVdat_education.SleepDuration_ON_CAD,
                                           NbDistribution = 1000, 
                                           SignifThreshold = 0.05)


#### ________M1.4 Sleep health_Insomnia ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_Insomnia <- IV_Insomnia["SNP"]

SNP_IV_education.Insomnia <- rbind(SNP_IV_educationStringent,SNP_IV_Insomnia)# 361SNP
SNP_IV_education.Insomnia <- unique(SNP_IV_education.Insomnia) # 0SNP
SNP_IV_education.Insomnia <- clump_data(SNP_IV_education.Insomnia) # 347SNP
write.xlsx(SNP_IV_education.Insomnia, file = "Output/MVMR_IV/SNP_IV_education.Insomnia.xlsx")

exp_education_MVwith.Insomnia <- format_data(
  GWAS_education,
  snps = SNP_IV_education.Insomnia$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.Insomnia$exposure <- "education"

exp_Insomnia_MVwith.education <- format_data(
  GWAS_Insomnia,
  snps = SNP_IV_education.Insomnia$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "MAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_Insomnia_MVwith.education$outcome <- "Insomnia"

education.Insomnia <- harmonise_data(exposure_dat = exp_education_MVwith.Insomnia, outcome_dat = exp_Insomnia_MVwith.education, action = 1)

out_CAD_MV.education.Insomnia <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.Insomnia$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.Insomnia$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.Insomnia, outcome_dat = out_CAD_MV.education.Insomnia, action = 1)

education.Insomnia.dat <- education.Insomnia[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.Insomnia.dat) <- c("SNP", "beta.education", "beta.Insomnia", "se.education", "se.Insomnia")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.Insomnia_ON_CAD <- merge(education.Insomnia.dat, education.CAD.dat, by = "SNP") # 330SNP
write.xlsx(MVdat_education.Insomnia_ON_CAD, file = "Output/MVMR_IV/MVdat_education.Insomnia_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.Insomnia_ON_CAD[,c("beta.education","beta.Insomnia")])
bxse = as.matrix(MVdat_education.Insomnia_ON_CAD[,c("se.education","se.Insomnia")])
by = as.vector(MVdat_education.Insomnia_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.Insomnia_ON_CAD$se.CAD)

MVdatForm_education.Insomnia_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","Insomnia"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.Insomnia_ON_CAD <- mr_mvivw(MVdatForm_education.Insomnia_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.Insomnia_ON_CAD <- mr_mvegger(MVdatForm_education.Insomnia_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.Insomnia_ON_CAD <- mr_mvlasso(MVdatForm_education.Insomnia_ON_CAD)

# MV Median
result_MV.Median_education.Insomnia_ON_CAD <- mr_mvmedian(MVdatForm_education.Insomnia_ON_CAD)



F_education.Insomnia <- strength_mvmr(r_input = MVdatForm_education.Insomnia_ON_CAD, gencov = 0)


mv_hete_education.Insomnia <- pleiotropy_mvmr(r_input = MVdatForm_education.Insomnia_ON_CAD, gencov = 0)


PRESSO_education.Insomnia_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                     BetaExposure = c("beta.education", "beta.Insomnia"), 
                                                     SdOutcome = "se.CAD", 
                                                     SdExposure = c("se.education", "se.Insomnia"),
                                                     OUTLIERtest = TRUE, 
                                                     DISTORTIONtest = TRUE, 
                                                     data = MVdat_education.Insomnia_ON_CAD,
                                                     NbDistribution = 1000, 
                                                     SignifThreshold = 0.05)


#### ________M2.1 Well-being spectrum(WBS) ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_WBS <- IV_WBS["SNP"]
SNP_IV_education.WBS <- rbind(SNP_IV_educationStringent,SNP_IV_WBS)# 478SNP
SNP_IV_education.WBS <- unique(SNP_IV_education.WBS) # 1SNP
SNP_IV_education.WBS <- clump_data(SNP_IV_education.WBS) # 363SNP
write.xlsx(SNP_IV_education.WBS, file = "Output/MVMR_IV/SNP_IV_education.WBS.xlsx")

exp_education_MVwith.WBS <- format_data(
  GWAS_education,
  snps = SNP_IV_education.WBS$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.WBS$exposure <- "education"

exp_WBS_MVwith.education <- format_data(
  GWAS_WBS,
  snps = SNP_IV_education.WBS$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "PVAL") 
exp_WBS_MVwith.education$outcome <- "WBS"

education.WBS <- harmonise_data(exposure_dat = exp_education_MVwith.WBS, outcome_dat = exp_WBS_MVwith.education, action = 1)

out_CAD_MV.education.WBS <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.WBS$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.WBS$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.WBS, outcome_dat = out_CAD_MV.education.WBS, action = 1)

education.WBS.dat <- education.WBS[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.WBS.dat) <- c("SNP", "beta.education", "beta.WBS", "se.education", "se.WBS")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.WBS_ON_CAD <- merge(education.WBS.dat, education.CAD.dat, by = "SNP") # 310SNP
write.xlsx(MVdat_education.WBS_ON_CAD, file = "Output/MVMR_IV/MVdat_education.WBS_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.WBS_ON_CAD[,c("beta.education","beta.WBS")])
bxse = as.matrix(MVdat_education.WBS_ON_CAD[,c("se.education","se.WBS")])
by = as.vector(MVdat_education.WBS_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.WBS_ON_CAD$se.CAD)

MVdatForm_education.WBS_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","WBS"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.WBS_ON_CAD <- mr_mvivw(MVdatForm_education.WBS_ON_CAD)
tt <- mvmr(MVdatForm_education.WBS_ON_CAD)#F-statistic: 31.87 on 2 and 308 DF

# MV MR-Egger
result_MV.Egger_education.WBS_ON_CAD <- mr_mvegger(MVdatForm_education.WBS_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.WBS_ON_CAD <- mr_mvlasso(MVdatForm_education.WBS_ON_CAD)

# MV Median
result_MV.Median_education.WBS_ON_CAD <- mr_mvmedian(MVdatForm_education.WBS_ON_CAD)



F_education.WBS <- strength_mvmr(r_input = MVdatForm_education.WBS_ON_CAD, gencov = 0)


mv_hete_education.WBS <- pleiotropy_mvmr(r_input = MVdatForm_education.WBS_ON_CAD, gencov = 0)


PRESSO_education.WBS_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.education", "beta.WBS"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.education", "se.WBS"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_education.WBS_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)


#### ________M2.2 MentalProblems ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_MentalProblems <- IV_MentalProblems["SNP"]

SNP_IV_education.MentalProblems <- rbind(SNP_IV_educationStringent,SNP_IV_MentalProblems) # 362SNP
SNP_IV_education.MentalProblems <- unique(SNP_IV_education.MentalProblems) # 0SNP
SNP_IV_education.MentalProblems <- clump_data(SNP_IV_education.MentalProblems) # 348SNP
write.xlsx(SNP_IV_education.MentalProblems, file = "Output/MVMR_IV/SNP_IV_education.MentalProblems.xlsx")

exp_education_MVwith.MentalProblems <- format_data(
  GWAS_education,
  snps = SNP_IV_education.MentalProblems$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.MentalProblems$exposure <- "education"

exp_MentalProblems_MVwith.education <- format_data(
  GWAS_MentalProblems,
  snps = SNP_IV_education.MentalProblems$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_MentalProblems_MVwith.education$outcome <- "MentalProblems"

education.MentalProblems <- harmonise_data(exposure_dat = exp_education_MVwith.MentalProblems, outcome_dat = exp_MentalProblems_MVwith.education, action = 1)

out_CAD_MV.education.MentalProblems <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.MentalProblems$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.MentalProblems$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.MentalProblems, outcome_dat = out_CAD_MV.education.MentalProblems, action = 1)

education.MentalProblems.dat <- education.MentalProblems[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.MentalProblems.dat) <- c("SNP", "beta.education", "beta.MentalProblems", "se.education", "se.MentalProblems")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.MentalProblems_ON_CAD <- merge(education.MentalProblems.dat, education.CAD.dat, by = "SNP") # 345SNP
write.xlsx(MVdat_education.MentalProblems_ON_CAD, file = "Output/MVMR_IV/MVdat_education.MentalProblems_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.MentalProblems_ON_CAD[,c("beta.education","beta.MentalProblems")])
bxse = as.matrix(MVdat_education.MentalProblems_ON_CAD[,c("se.education","se.MentalProblems")])
by = as.vector(MVdat_education.MentalProblems_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.MentalProblems_ON_CAD$se.CAD)

MVdatForm_education.MentalProblems_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","MentalProblems"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.MentalProblems_ON_CAD <- mr_mvivw(MVdatForm_education.MentalProblems_ON_CAD)

# MV MR-Egger
result_MV.Egger_education.MentalProblems_ON_CAD <- mr_mvegger(MVdatForm_education.MentalProblems_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.MentalProblems_ON_CAD <- mr_mvlasso(MVdatForm_education.MentalProblems_ON_CAD)

# MV Median
result_MV.Median_education.MentalProblems_ON_CAD <- mr_mvmedian(MVdatForm_education.MentalProblems_ON_CAD)



F_education.MentalProblems <- strength_mvmr(r_input = MVdatForm_education.MentalProblems_ON_CAD, gencov = 0)


mv_hete_education.MentalProblems <- pleiotropy_mvmr(r_input = MVdatForm_education.MentalProblems_ON_CAD, gencov = 0)


PRESSO_education.MentalProblems_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.education", "beta.MentalProblems"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.education", "se.MentalProblems"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_education.MentalProblems_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)



#### ________M2.3 Depression ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_depression <- IV_depression["SNP"]

SNP_IV_education.depression <- rbind(SNP_IV_educationStringent,SNP_IV_depression) # 450SNP
SNP_IV_education.depression <- unique(SNP_IV_education.depression) # 1SNP
SNP_IV_education.depression <- clump_data(SNP_IV_education.depression) # 360SNP
write.xlsx(SNP_IV_education.depression, file = "Output/MVMR_IV/SNP_IV_education.depression.xlsx")

exp_education_MVwith.depression <- format_data(
  GWAS_education,
  snps = SNP_IV_education.depression$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.depression$exposure <- "education"

exp_depression_MVwith.education <- format_data(
  GWAS_depression,
  snps = SNP_IV_education.depression$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "LogOR",
  se_col = "StdErrLogOR",
  eaf_col = "Freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_depression_MVwith.education$outcome <- "depression"

education.depression <- harmonise_data(exposure_dat = exp_education_MVwith.depression, outcome_dat = exp_depression_MVwith.education, action = 1)

out_CAD_MV.education.depression <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.depression$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.depression$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.depression, outcome_dat = out_CAD_MV.education.depression, action = 1)

education.depression.dat <- education.depression[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.depression.dat) <- c("SNP", "beta.education", "beta.depression", "se.education", "se.depression")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.depression_ON_CAD <- merge(education.depression.dat, education.CAD.dat, by = "SNP") # 352SNP
write.xlsx(MVdat_education.depression_ON_CAD, file = "Output/MVMR_IV/MVdat_education.depression_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.depression_ON_CAD[,c("beta.education","beta.depression")])
bxse = as.matrix(MVdat_education.depression_ON_CAD[,c("se.education","se.depression")])
by = as.vector(MVdat_education.depression_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.depression_ON_CAD$se.CAD)

MVdatForm_education.depression_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","depression"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.depression_ON_CAD <- mr_mvivw(MVdatForm_education.depression_ON_CAD)
tt <- mvmr(MVdatForm_education.depression_ON_CAD)#F-statistic: 34.43 on 2 and 350 DF

# MV MR-Egger
result_MV.Egger_education.depression_ON_CAD <- mr_mvegger(MVdatForm_education.depression_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.depression_ON_CAD <- mr_mvlasso(MVdatForm_education.depression_ON_CAD)

# MV Median
result_MV.Median_education.depression_ON_CAD <- mr_mvmedian(MVdatForm_education.depression_ON_CAD)



F_education.depression <- strength_mvmr(r_input = MVdatForm_education.depression_ON_CAD, gencov = 0)


mv_hete_education.depression <- pleiotropy_mvmr(r_input = MVdatForm_education.depression_ON_CAD, gencov = 0)


PRESSO_education.depression_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                                BetaExposure = c("beta.education", "beta.depression"), 
                                                SdOutcome = "se.CAD", 
                                                SdExposure = c("se.education", "se.depression"),
                                                OUTLIERtest = TRUE, 
                                                DISTORTIONtest = TRUE, 
                                                data = MVdat_education.depression_ON_CAD,
                                                NbDistribution = 1000, 
                                                SignifThreshold = 0.05)


#### ________M3.1 BMI ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_BMI <- IV_BMI["SNP"]

SNP_IV_education.BMI <- rbind(SNP_IV_educationStringent,SNP_IV_BMI) # 858SNP
SNP_IV_education.BMI <- unique(SNP_IV_education.BMI) # 0SNP
SNP_IV_education.BMI <- clump_data(SNP_IV_education.BMI) # 536SNP
write.xlsx(SNP_IV_education.BMI, file = "Output/MVMR_IV/SNP_IV_education.BMI.xlsx")

exp_education_MVwith.BMI <- format_data(
  GWAS_education,
  snps = SNP_IV_education.BMI$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.BMI$exposure <- "education"

exp_BMI_MVwith.education <- format_data(
  GWAS_BMI,
  snps = SNP_IV_education.BMI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_BMI_MVwith.education$outcome <- "BMI"

education.BMI <- harmonise_data(exposure_dat = exp_education_MVwith.BMI, outcome_dat = exp_BMI_MVwith.education, action = 1)

out_CAD_MV.education.BMI <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.BMI$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.BMI$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.BMI, outcome_dat = out_CAD_MV.education.BMI, action = 1)

education.BMI.dat <- education.BMI[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.BMI.dat) <- c("SNP", "beta.education", "beta.BMI", "se.education", "se.BMI")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.BMI_ON_CAD <- merge(education.BMI.dat, education.CAD.dat, by = "SNP") # 441SNP
write.xlsx(MVdat_education.BMI_ON_CAD, file = "Output/MVMR_IV/MVdat_education.BMI_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.BMI_ON_CAD[,c("beta.education","beta.BMI")])
bxse = as.matrix(MVdat_education.BMI_ON_CAD[,c("se.education","se.BMI")])
by = as.vector(MVdat_education.BMI_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.BMI_ON_CAD$se.CAD)

MVdatForm_education.BMI_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","BMI"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.BMI_ON_CAD <- mr_mvivw(MVdatForm_education.BMI_ON_CAD)
tt <- mvmr(MVdatForm_education.BMI_ON_CAD)#F-statistic: 53.46 on 2 and 439 DF

# MV MR-Egger
result_MV.Egger_education.BMI_ON_CAD <- mr_mvegger(MVdatForm_education.BMI_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.BMI_ON_CAD <- mr_mvlasso(MVdatForm_education.BMI_ON_CAD)

# MV Median
result_MV.Median_education.BMI_ON_CAD <- mr_mvmedian(MVdatForm_education.BMI_ON_CAD)



F_education.BMI <- strength_mvmr(r_input = MVdatForm_education.BMI_ON_CAD, gencov = 0)


mv_hete_education.BMI <- pleiotropy_mvmr(r_input = MVdatForm_education.BMI_ON_CAD, gencov = 0)


PRESSO_education.BMI_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                         BetaExposure = c("beta.education", "beta.BMI"), 
                                         SdOutcome = "se.CAD", 
                                         SdExposure = c("se.education", "se.BMI"),
                                         OUTLIERtest = TRUE, 
                                         DISTORTIONtest = TRUE, 
                                         data = MVdat_education.BMI_ON_CAD,
                                         NbDistribution = 1000, 
                                         SignifThreshold = 0.05)


#### ________M3.2 blood lipids ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_lipid <- IV_lipid["SNP"]

SNP_IV_education.lipid <- rbind(SNP_IV_educationStringent,SNP_IV_lipid) # 733SNP
SNP_IV_education.lipid <- unique(SNP_IV_education.lipid) # 0SNP
SNP_IV_education.lipid <- clump_data(SNP_IV_education.lipid) # 492SNP
write.xlsx(SNP_IV_education.lipid, file = "Output/MVMR_IV/SNP_IV_education.lipid.xlsx")

exp_education_MVwith.lipid <- format_data(
  GWAS_education,
  snps = SNP_IV_education.lipid$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.lipid$exposure <- "education"

exp_lipid_MVwith.education <- format_data(
  GWAS_lipid,
  snps = SNP_IV_education.lipid$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_lipid_MVwith.education$outcome <- "lipid"

education.lipid <- harmonise_data(exposure_dat = exp_education_MVwith.lipid, outcome_dat = exp_lipid_MVwith.education, action = 1)

out_CAD_MV.education.lipid <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.lipid$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.lipid$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.lipid, outcome_dat = out_CAD_MV.education.lipid, action = 1)

education.lipid.dat <- education.lipid[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.lipid.dat) <- c("SNP", "beta.education", "beta.lipid", "se.education", "se.lipid")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.lipid_ON_CAD <- merge(education.lipid.dat, education.CAD.dat, by = "SNP") # 476SNP
write.xlsx(MVdat_education.lipid_ON_CAD, file = "Output/MVMR_IV/MVdat_education.lipid_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.lipid_ON_CAD[,c("beta.education","beta.lipid")])
bxse = as.matrix(MVdat_education.lipid_ON_CAD[,c("se.education","se.lipid")])
by = as.vector(MVdat_education.lipid_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.lipid_ON_CAD$se.CAD)

MVdatForm_education.lipid_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","lipid"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.lipid_ON_CAD <- mr_mvivw(MVdatForm_education.lipid_ON_CAD)
tt <- mvmr(MVdatForm_education.lipid_ON_CAD)#F-statistic: 133.91 on 2 and 474 DF

# MV MR-Egger
result_MV.Egger_education.lipid_ON_CAD <- mr_mvegger(MVdatForm_education.lipid_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.lipid_ON_CAD <- mr_mvlasso(MVdatForm_education.lipid_ON_CAD)

# MV Median
result_MV.Median_education.lipid_ON_CAD <- mr_mvmedian(MVdatForm_education.lipid_ON_CAD)



F_education.lipid <- strength_mvmr(r_input = MVdatForm_education.lipid_ON_CAD, gencov = 0)


mv_hete_education.lipid <- pleiotropy_mvmr(r_input = MVdatForm_education.lipid_ON_CAD, gencov = 0)


PRESSO_education.lipid_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                           BetaExposure = c("beta.education", "beta.lipid"),
                                           SdOutcome = "se.CAD",
                                           SdExposure = c("se.education", "se.lipid"),
                                           OUTLIERtest = TRUE,
                                           DISTORTIONtest = TRUE,
                                           data = MVdat_education.lipid_ON_CAD,
                                           NbDistribution = 1000,
                                           SignifThreshold = 0.05)


#### ________M3.3 blood glucose ####
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_glucose <- IV_glucose["SNP"]

SNP_IV_education.glucose <- rbind(SNP_IV_educationStringent,SNP_IV_glucose) # 423SNP
SNP_IV_education.glucose <- unique(SNP_IV_education.glucose) # 0SNP
SNP_IV_education.glucose <- clump_data(SNP_IV_education.glucose) # 371SNP
write.xlsx(SNP_IV_education.glucose, file = "Output/MVMR_IV/SNP_IV_education.glucose.xlsx")

exp_education_MVwith.glucose <- format_data(
  GWAS_education,
  snps = SNP_IV_education.glucose$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.glucose$exposure <- "education"

exp_glucose_MVwith.education <- format_data(
  GWAS_glucose,
  snps = SNP_IV_education.glucose$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_glucose_MVwith.education$outcome <- "glucose"

education.glucose <- harmonise_data(exposure_dat = exp_education_MVwith.glucose, outcome_dat = exp_glucose_MVwith.education, action = 1)

out_CAD_MV.education.glucose <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.glucose$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.glucose$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.glucose, outcome_dat = out_CAD_MV.education.glucose, action = 1)

education.glucose.dat <- education.glucose[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.glucose.dat) <- c("SNP", "beta.education", "beta.glucose", "se.education", "se.glucose")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.glucose_ON_CAD <- merge(education.glucose.dat, education.CAD.dat, by = "SNP") # 363SNP
write.xlsx(MVdat_education.glucose_ON_CAD, file = "Output/MVMR_IV/MVdat_education.glucose_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.glucose_ON_CAD[,c("beta.education","beta.glucose")])
bxse = as.matrix(MVdat_education.glucose_ON_CAD[,c("se.education","se.glucose")])
by = as.vector(MVdat_education.glucose_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.glucose_ON_CAD$se.CAD)

MVdatForm_education.glucose_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","glucose"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.glucose_ON_CAD <- mr_mvivw(MVdatForm_education.glucose_ON_CAD)
tt <- mvmr(MVdatForm_education.glucose_ON_CAD)#F-statistic: 40.23 on 2 and 361 DF

# MV MR-Egger
result_MV.Egger_education.glucose_ON_CAD <- mr_mvegger(MVdatForm_education.glucose_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.glucose_ON_CAD <- mr_mvlasso(MVdatForm_education.glucose_ON_CAD)

# MV Median
result_MV.Median_education.glucose_ON_CAD <- mr_mvmedian(MVdatForm_education.glucose_ON_CAD)



F_education.glucose <- strength_mvmr(r_input = MVdatForm_education.glucose_ON_CAD, gencov = 0)


mv_hete_education.glucose <- pleiotropy_mvmr(r_input = MVdatForm_education.glucose_ON_CAD, gencov = 0)


PRESSO_education.glucose_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                             BetaExposure = c("beta.education", "beta.glucose"), 
                                             SdOutcome = "se.CAD", 
                                             SdExposure = c("se.education", "se.glucose"),
                                             OUTLIERtest = TRUE, 
                                             DISTORTIONtest = TRUE, 
                                             data = MVdat_education.glucose_ON_CAD,
                                             NbDistribution = 1000, 
                                             SignifThreshold = 0.05)


#### ________M3.4 blood pressure(BP) ####
## SBP
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_SBP <- IV_SBP["SNP"]

SNP_IV_education.SBP <- rbind(SNP_IV_educationStringent,SNP_IV_SBP) # 804SNP
SNP_IV_education.SBP <- unique(SNP_IV_education.SBP) # 0SNP
SNP_IV_education.SBP <- clump_data(SNP_IV_education.SBP) # 527SNP
write.xlsx(SNP_IV_education.SBP, file = "Output/MVMR_IV/SNP_IV_education.SBP.xlsx")

exp_education_MVwith.SBP <- format_data(
  GWAS_education,
  snps = SNP_IV_education.SBP$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.SBP$exposure <- "education"

exp_SBP_MVwith.education <- format_data(
  GWAS_SBP,
  snps = SNP_IV_education.SBP$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_SBP_MVwith.education$outcome <- "SBP"

education.SBP <- harmonise_data(exposure_dat = exp_education_MVwith.SBP, outcome_dat = exp_SBP_MVwith.education, action = 1)

out_CAD_MV.education.SBP <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.SBP$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.SBP$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.SBP, outcome_dat = out_CAD_MV.education.SBP, action = 1)

education.SBP.dat <- education.SBP[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.SBP.dat) <- c("SNP", "beta.education", "beta.SBP", "se.education", "se.SBP")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.SBP_ON_CAD <- merge(education.SBP.dat, education.CAD.dat, by = "SNP") # 512SNP
write.xlsx(MVdat_education.SBP_ON_CAD, file = "Output/MVMR_IV/MVdat_education.SBP_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.SBP_ON_CAD[,c("beta.education","beta.SBP")])
bxse = as.matrix(MVdat_education.SBP_ON_CAD[,c("se.education","se.SBP")])
by = as.vector(MVdat_education.SBP_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.SBP_ON_CAD$se.CAD)

MVdatForm_education.SBP_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","SBP"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.SBP_ON_CAD <- mr_mvivw(MVdatForm_education.SBP_ON_CAD)
tt <- mvmr(MVdatForm_education.SBP_ON_CAD)#F-statistic: 67.96 on 2 and 510 DF

# MV MR-Egger
result_MV.Egger_education.SBP_ON_CAD <- mr_mvegger(MVdatForm_education.SBP_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.SBP_ON_CAD <- mr_mvlasso(MVdatForm_education.SBP_ON_CAD)

# MV Median
result_MV.Median_education.SBP_ON_CAD <- mr_mvmedian(MVdatForm_education.SBP_ON_CAD)



F_education.SBP <- strength_mvmr(r_input = MVdatForm_education.SBP_ON_CAD, gencov = 0)


mv_hete_education.SBP <- pleiotropy_mvmr(r_input = MVdatForm_education.SBP_ON_CAD, gencov = 0)


PRESSO_education.SBP_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                         BetaExposure = c("beta.education", "beta.SBP"), 
                                         SdOutcome = "se.CAD", 
                                         SdExposure = c("se.education", "se.SBP"),
                                         OUTLIERtest = TRUE, 
                                         DISTORTIONtest = TRUE, 
                                         data = MVdat_education.SBP_ON_CAD,
                                         NbDistribution = 1000, 
                                         SignifThreshold = 0.05)


## DBP
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_DBP <- IV_DBP["SNP"]

SNP_IV_education.DBP <- rbind(SNP_IV_educationStringent,SNP_IV_DBP) # 803SNP
SNP_IV_education.DBP <- unique(SNP_IV_education.DBP) # 0SNP
SNP_IV_education.DBP <- clump_data(SNP_IV_education.DBP) # 532SNP
write.xlsx(SNP_IV_education.DBP, file = "Output/MVMR_IV/SNP_IV_education.DBP.xlsx")

exp_education_MVwith.DBP <- format_data(
  GWAS_education,
  snps = SNP_IV_education.DBP$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.DBP$exposure <- "education"

exp_DBP_MVwith.education <- format_data(
  GWAS_DBP,
  snps = SNP_IV_education.DBP$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_DBP_MVwith.education$outcome <- "DBP"

education.DBP <- harmonise_data(exposure_dat = exp_education_MVwith.DBP, outcome_dat = exp_DBP_MVwith.education, action = 1)

out_CAD_MV.education.DBP <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.DBP$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.DBP$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.DBP, outcome_dat = out_CAD_MV.education.DBP, action = 1)

education.DBP.dat <- education.DBP[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.DBP.dat) <- c("SNP", "beta.education", "beta.DBP", "se.education", "se.DBP")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

MVdat_education.DBP_ON_CAD <- merge(education.DBP.dat, education.CAD.dat, by = "SNP") # 363SNP
write.xlsx(MVdat_education.DBP_ON_CAD, file = "Output/MVMR_IV/MVdat_education.DBP_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.DBP_ON_CAD[,c("beta.education","beta.DBP")])
bxse = as.matrix(MVdat_education.DBP_ON_CAD[,c("se.education","se.DBP")])
by = as.vector(MVdat_education.DBP_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.DBP_ON_CAD$se.CAD)

MVdatForm_education.DBP_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education","DBP"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.DBP_ON_CAD <- mr_mvivw(MVdatForm_education.DBP_ON_CAD)
tt <- mvmr(MVdatForm_education.DBP_ON_CAD)#F-statistic: 91.29 on 2 and 509 DF

# MV MR-Egger
result_MV.Egger_education.DBP_ON_CAD <- mr_mvegger(MVdatForm_education.DBP_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.DBP_ON_CAD <- mr_mvlasso(MVdatForm_education.DBP_ON_CAD)

# MV Median
result_MV.Median_education.DBP_ON_CAD <- mr_mvmedian(MVdatForm_education.DBP_ON_CAD)



F_education.DBP <- strength_mvmr(r_input = MVdatForm_education.DBP_ON_CAD, gencov = 0)


mv_hete_education.DBP <- pleiotropy_mvmr(r_input = MVdatForm_education.DBP_ON_CAD, gencov = 0)


PRESSO_education.DBP_ON_CAD <- mr_presso(BetaOutcome = "beta.CAD",
                                         BetaExposure = c("beta.education", "beta.DBP"), 
                                         SdOutcome = "se.CAD", 
                                         SdExposure = c("se.education", "se.DBP"),
                                         OUTLIERtest = TRUE, 
                                         DISTORTIONtest = TRUE, 
                                         data = MVdat_education.DBP_ON_CAD,
                                         NbDistribution = 1000, 
                                         SignifThreshold = 0.05)



###### Main method: product of coefficients method 
# Individual Mediation effect: β1×β2, 95%CI
# Individual Mediation proportion: β1×β2/β, 95%CI
# divide the result of the indirect effect (β1×β2) by the total causal effect of X on Y estimated by UVMR (β)

med_cal <- function(beta1,se1,beta2,se2,beta_total,n=1000000){
  mediate_effect <- beta1 * beta2
  mediate_se <- sqrt((beta1^2 * se2^2) + (beta2^2 * se1^2))
  m_low <- mediate_effect - 1.96* mediate_se
  m_upper <- mediate_effect + 1.96* mediate_se
  Z <- mediate_effect/mediate_se
  p <- 2*pnorm(q=abs(Z), lower.tail=FALSE)
  
  total_effect <- beta_total
  
  direct_effect <- total_effect - mediate_effect
  
  mediate_percentage <- (mediate_effect / total_effect) * 100
  
  mediate_percentages <- numeric(n)
  
  set.seed(123456789)  
  for (i in 1:n) {
    beta1_boot <- rnorm(1, beta1, se1)
    beta2_boot <- rnorm(1, beta2, se2)
    mediate_effect_boot <- beta1_boot * beta2_boot
    direct_effect_boot <- total_effect - mediate_effect_boot
    mediate_percentages[i] <- (mediate_effect_boot / total_effect) * 100
  }
  
  pro <- (beta1*beta2 / total_effect) * 100
  se <- sd(mediate_percentages)
  pro_low <- pro - 1.96*se 
  pro_upper <- pro + 1.96*se 
  res <- data.frame(mediate_effect,m_low,m_upper,p,pro,pro_low,pro_upper)
  return(res)
}



#### ________M1.2 Physical activity_PA ####
result_education_ON_PA
result_MV.IVW_education.PA_ON_CAD
result_education_ON_CAD

med_cal_PA_IN_edu.CAD <- med_cal(beta1 = 0.4941715, se1 = 0.01468776, beta2 = -0.165, se2 = 0.105, beta_total = -0.4790515)


#### ________M1.2 Physical activity_LST ####
result_education_ON_LST
result_MV.IVW_education.LST_ON_CAD
result_education_ON_CAD

med_cal_LST_IN_edu.CAD <- med_cal(beta1 = -0.6576202, se1 = 0.01831928, beta2 = 0.074, se2 = 0.074, beta_total = -0.4790515)


#### ________M1.2 Physical activity_television ####
result_education_ON_television
result_MV.IVW_education.television_ON_CAD
result_education_ON_CAD

med_cal_television_IN_edu.CAD <- med_cal(beta1 = -0.5582677, se1 = 0.01027837, beta2 = 0.251, se2 = 0.106, beta_total = -0.4790515)


#### ________M1.2 Physical activity_computer ####
result_education_ON_computer
result_MV.IVW_education.computer_ON_CAD
result_education_ON_CAD

med_cal_computer_IN_edu.CAD <- med_cal(beta1 = 0.3112223, se1 = 0.01002613, beta2 = 0.208, se2 = 0.148, beta_total = -0.4790515)


#### ________M1.2 Physical activity_driving ####
result_education_ON_driving
result_MV.IVW_education.driving_ON_CAD
result_education_ON_CAD

med_cal_driving_IN_edu.CAD <- med_cal(beta1 = -0.10453785, se1 = 0.008651701, beta2 = 0.162, se2 = 0.195, beta_total = -0.4790515)


#### ________M1.3 Smoking_LifetimeSmoking ####
result_education_ON_LifetimeSmoking
result_MV.IVW_education.LifetimeSmoking_ON_CAD
result_education_ON_CAD

med_cal_LifetimeSmoking_IN_edu.CAD <- med_cal(beta1 = -0.2322174, se1 = 0.006734792, beta2 = 0.692, se2 = 0.159, beta_total = -0.4790515) 


#### ________M1.4 Sleep health_SleepDuration ####
result_education_ON_SleepDuration
result_MV.IVW_education.SleepDuration_ON_CAD
result_education_ON_CAD

med_cal_SleepDuration_IN_edu.CAD <- med_cal(beta1 = 0.06460235, se1 = 0.01157275, beta2 = -0.105, se2 = 0.120, beta_total = -0.4790515) 


#### ________M1.4 Sleep health_Insomnia ####
result_education_ON_Insomnia
result_MV.IVW_education.Insomnia_ON_CAD
result_education_ON_CAD

med_cal_Insomnia_IN_edu.CAD <- med_cal(beta1 = -0.3298880, se1 = 0.02094483, beta2 = 0.162, se2 = 0.084, beta_total = -0.4790515) 


#### ________M2.1 Well-being spectrum(WBS) ####
result_education_ON_WBS
result_MV.IVW_education.WBS_ON_CAD
result_education_ON_CAD

med_cal_WBS_IN_edu.CAD <- med_cal(beta1 = 0.06175960, se1 = 0.005480236, beta2 = -0.603, se2 = 0.216, beta_total = -0.4790515) 


#### ________M2.2 MentalProblems ####
result_education_ON_MentalProblems
result_MV.IVW_education.MentalProblems_ON_CAD
result_education_ON_CAD

med_cal_MentalProblems_IN_edu.CAD <- med_cal(beta1 = -0.05703687, se1 = 0.005132001, beta2 = 0.411, se2 = 0.314, beta_total = -0.4790515) 


#### ________M2.3 Depression ####
result_education_ON_depression
result_MV.IVW_education.depression_ON_CAD
result_education_ON_CAD

med_cal_depression_IN_edu.CAD <- med_cal(beta1 = -0.2515700, se1 = 0.02090266, beta2 = 0.122, se2 = 0.061, beta_total = -0.4790515) 


#### ________M3.1 BMI ####
result_education_ON_BMI
result_MV.IVW_education.BMI_ON_CAD
result_education_ON_CAD

med_cal_BMI_IN_edu.CAD <- med_cal(beta1 = -0.1916290, se1 = 0.01458166, beta2 = 0.336, se2 = 0.045, beta_total = -0.4790515) 


#### ________M3.2 blood lipids ####
result_education_ON_lipid
result_MV.IVW_education.lipid_ON_CAD
result_education_ON_CAD

med_cal_lipid_IN_edu.CAD <- med_cal(beta1 = -0.11548920, se1 = 0.007464035, beta2 = 0.598, se2 = 0.039, beta_total = -0.4790515)  


#### ________M3.3 blood glucose ####
result_education_ON_glucose
result_MV.IVW_education.glucose_ON_CAD
result_education_ON_CAD

med_cal_glucose_IN_edu.CAD <- med_cal(beta1 = -0.0269212163, se1 = 0.005049031, beta2 = 0.495, se2 = 0.140, beta_total = -0.4790515) 


#### ________M3.4 blood pressure(BP) ####
result_education_ON_SBP
result_MV.IVW_education.SBP_ON_CAD
med_cal_SBP_IN_edu.CAD <- med_cal(beta1 = -1.7494468, se1 = 0.1826313, beta2 = 0.032, se2 = 0.003, beta_total = -0.4790515) 

result_education_ON_DBP
result_MV.IVW_education.DBP_ON_CAD
med_cal_DBP_IN_edu.CAD <- med_cal(beta1 = -0.7875675, se1 = 0.09982863, beta2 = 0.058, se2 = 0.005, beta_total = -0.4790515) 


###### ________Individual mediation######
SummaryResults_IndividualMediation_education <- rbind(med_cal_PA_IN_edu.CAD,
                                                      med_cal_LST_IN_edu.CAD,
                                                      med_cal_television_IN_edu.CAD,
                                                      med_cal_computer_IN_edu.CAD,
                                                      med_cal_driving_IN_edu.CAD,
                                                      med_cal_LifetimeSmoking_IN_edu.CAD,
                                                      med_cal_SleepDuration_IN_edu.CAD,
                                                      med_cal_Insomnia_IN_edu.CAD,
                                                      med_cal_WBS_IN_edu.CAD,
                                                      med_cal_MentalProblems_IN_edu.CAD,
                                                      med_cal_depression_IN_edu.CAD,
                                                      med_cal_BMI_IN_edu.CAD,
                                                      med_cal_lipid_IN_edu.CAD,
                                                      med_cal_glucose_IN_edu.CAD,
                                                      med_cal_SBP_IN_edu.CAD,
                                                      med_cal_DBP_IN_edu.CAD)
write.xlsx(SummaryResults_IndividualMediation_education,"Output/SummaryResults_IndividualMediation_education.xlsx")


###### ___3.2）Combined proportion mediated ######
# Combined proportions mediated by multiple mediators
# We examined the proportion mediated of different combinations of mediating variables

###### ________9：television + LifetimeSmoking + WBS + depression + BMI + lipid + glucose + SBP + DBP ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_television <- IV_television["SNP"]
SNP_IV_LifetimeSmoking <- IV_LifetimeSmoking["SNP"]
SNP_IV_WBS <- IV_WBS["SNP"]
SNP_IV_depression <- IV_depression["SNP"]
SNP_IV_BMI <- IV_BMI["SNP"]
SNP_IV_lipid <- IV_lipid["SNP"]
SNP_IV_glucose <- IV_glucose["SNP"]
SNP_IV_SBP <- IV_SBP["SNP"]
SNP_IV_DBP <- IV_DBP["SNP"]

SNP_IV_education.9mediators <- rbind(SNP_IV_educationStringent,
                                     SNP_IV_television,
                                     SNP_IV_LifetimeSmoking,
                                     SNP_IV_WBS,
                                     SNP_IV_depression,
                                     SNP_IV_BMI,
                                     SNP_IV_lipid,
                                     SNP_IV_glucose,
                                     SNP_IV_SBP,
                                     SNP_IV_DBP) # 2748SNP
# ttt <- SNP_IV_education.9mediators
# tttt <- ttt [duplicated(ttt),]
# ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
# SNP_IV_education.9mediators <- ttttt

SNP_IV_education.9mediators <- unique(SNP_IV_education.9mediators)
tt <- clump_data(SNP_IV_education.9mediators)
SNP_IV_education.9mediators <- clump_data(SNP_IV_education.9mediators) # 695SNP
write.xlsx(SNP_IV_education.9mediators, file = "Output/MVMR_IV/SNP_IV_education.9mediators.xlsx")

exp_education_MVwith.9mediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.9mediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.9mediators$exposure <- "education"

exp_television_MVwith.education.9mediators <- format_data(
  GWAS_television,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_television_MVwith.education.9mediators$outcome <- "television"

exp_LifetimeSmoking_MVwith.education.9mediators <- format_data(
  GWAS_LifetimeSmoking,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "P") 
exp_LifetimeSmoking_MVwith.education.9mediators$outcome <- "LifetimeSmoking"

exp_WBS_MVwith.education.9mediators <- format_data(
  GWAS_WBS,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "PVAL") 
exp_WBS_MVwith.education.9mediators$outcome <- "WBS"

exp_depression_MVwith.education.9mediators <- format_data(
  GWAS_depression,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "LogOR",
  se_col = "StdErrLogOR",
  eaf_col = "Freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_depression_MVwith.education.9mediators$outcome <- "depression"

exp_BMI_MVwith.education.9mediators <- format_data(
  GWAS_BMI,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_BMI_MVwith.education.9mediators$outcome <- "BMI"

exp_lipid_MVwith.education.9mediators <- format_data(
  GWAS_lipid,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_lipid_MVwith.education.9mediators$outcome <- "lipid"

exp_glucose_MVwith.education.9mediators <- format_data(
  GWAS_glucose,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_glucose_MVwith.education.9mediators$outcome <- "glucose"

exp_SBP_MVwith.education.9mediators <- format_data(
  GWAS_SBP,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_SBP_MVwith.education.9mediators$outcome <- "SBP"

exp_DBP_MVwith.education.9mediators <- format_data(
  GWAS_DBP,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_DBP_MVwith.education.9mediators$outcome <- "DBP"

education.television <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_television_MVwith.education.9mediators, action = 1)
education.LifetimeSmoking <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_LifetimeSmoking_MVwith.education.9mediators, action = 1)
education.WBS <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_WBS_MVwith.education.9mediators, action = 1)
education.depression <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_depression_MVwith.education.9mediators, action = 1)
education.BMI <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_BMI_MVwith.education.9mediators, action = 1)
education.lipid <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_lipid_MVwith.education.9mediators, action = 1)
education.glucose <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_glucose_MVwith.education.9mediators, action = 1)
education.SBP <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_SBP_MVwith.education.9mediators, action = 1)
education.DBP <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = exp_DBP_MVwith.education.9mediators, action = 1)

out_CAD_MV.education.9mediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.9mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.9mediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.9mediators, outcome_dat = out_CAD_MV.education.9mediators, action = 1)

education.television.dat <- education.television[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.television.dat) <- c("SNP", "beta.education", "beta.television", "se.education", "se.television")

education.LifetimeSmoking.dat <- education.LifetimeSmoking[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.LifetimeSmoking.dat) <- c("SNP", "beta.LifetimeSmoking", "se.LifetimeSmoking")

education.WBS.dat <- education.WBS[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.WBS.dat) <- c("SNP", "beta.WBS", "se.WBS")

education.depression.dat <- education.depression[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.depression.dat) <- c("SNP", "beta.depression", "se.depression")

education.BMI.dat <- education.BMI[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.BMI.dat) <- c("SNP", "beta.BMI", "se.BMI")

education.lipid.dat <- education.lipid[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.lipid.dat) <- c("SNP", "beta.lipid", "se.lipid")

education.glucose.dat <- education.glucose[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.glucose.dat) <- c("SNP", "beta.glucose", "se.glucose")

education.SBP.dat <- education.SBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.SBP.dat) <- c("SNP", "beta.SBP", "se.SBP")

education.DBP.dat <- education.DBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.DBP.dat) <- c("SNP", "beta.DBP", "se.DBP")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.television.dat, education.LifetimeSmoking.dat, by = "SNP")
mvmr.dat2 <- merge(mvmr.dat1, education.WBS.dat, by = "SNP")
mvmr.dat3 <- merge(mvmr.dat2, education.depression.dat, by = "SNP")
mvmr.dat4 <- merge(mvmr.dat3, education.BMI.dat, by = "SNP")
mvmr.dat5 <- merge(mvmr.dat4, education.lipid.dat, by = "SNP")
mvmr.dat6 <- merge(mvmr.dat5, education.glucose.dat, by = "SNP")
mvmr.dat7 <- merge(mvmr.dat6, education.SBP.dat, by = "SNP")
mvmr.dat8 <- merge(mvmr.dat7, education.DBP.dat, by = "SNP")
MVdat_education.9mediators_ON_CAD <- merge(mvmr.dat8, education.CAD.dat, by = "SNP")# 314SNP
write.xlsx(MVdat_education.9mediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.9mediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.9mediators_ON_CAD[,c("beta.education", "beta.television", "beta.LifetimeSmoking", "beta.WBS","beta.depression", "beta.BMI", "beta.lipid", "beta.glucose", "beta.SBP", "beta.DBP")])
bxse = as.matrix(MVdat_education.9mediators_ON_CAD[,c("se.education", "se.television", "se.LifetimeSmoking", "se.WBS","se.depression", "se.BMI", "se.lipid", "se.glucose", "se.SBP", "se.DBP")])
by = as.vector(MVdat_education.9mediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.9mediators_ON_CAD$se.CAD)

MVdatForm_education.9mediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "television", "LifetimeSmoking", "WBS","depression", "BMI", "lipid", "glucose", "SBP", "DBP"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.9mediators_ON_CAD <- mr_mvivw(MVdatForm_education.9mediators_ON_CAD)
tt <- mvmr(MVdatForm_education.9mediators_ON_CAD)#F-statistic: 8.37 on 10 and 304 DF

# MV MR-Egger
result_MV.Egger_education.9mediators_ON_CAD <- mr_mvegger(MVdatForm_education.9mediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.9mediators_ON_CAD <- mr_mvlasso(MVdatForm_education.9mediators_ON_CAD)

# MV Median
result_MV.Median_education.9mediators_ON_CAD <- mr_mvmedian(MVdatForm_education.9mediators_ON_CAD)



F_education.9mediators <- strength_mvmr(r_input = MVdatForm_education.9mediators_ON_CAD, gencov = 0)


mv_hete_education.9mediators <- pleiotropy_mvmr(r_input = MVdatForm_education.9mediators_ON_CAD, gencov = 0)


CombineMed_cal <- function(beta1,se1,beta2,se2,n=1000000){
  total_effect <- beta1
  direct_effect <- beta2
  direct_percentage <- beta2/beta1
  CombineMed_percentage <- (1-direct_percentage) * 100
  
  
  CombineMed_percentages <- numeric(n)
  
  set.seed(123456789)
  for (i in 1:n) {
    beta1_boot <- rnorm(1, beta1, se1)
    beta2_boot <- rnorm(1, beta2, se2)
    direct_percentage_boot <- beta2_boot/beta1_boot
    CombineMed_percentages[i] <- (1-direct_percentage_boot)*100
  }
  
  
  CombineMed_percentage_se <- sd(CombineMed_percentages)
  CombineMed_percentage_low <- CombineMed_percentage - 1.96*CombineMed_percentage_se
  CombineMed_percentage_upper <- CombineMed_percentage + 1.96*CombineMed_percentage_se
  
  Z <- CombineMed_percentage/CombineMed_percentage_se 
  p <- 2*pnorm(q=abs(Z), lower.tail=FALSE) 
  
  res <- data.frame(CombineMed_percentage,CombineMed_percentage_low,CombineMed_percentage_upper,p)
  return(res)
}

CombineMed_CI <- function(beta1,se1,beta2,se2){
  total_effect <- beta1
  direct_effect <- beta2
  direct_percentage <- beta2/beta1
  CombineMed_percentage <- (1-direct_percentage) * 100

  
  direct_percentage_low <- beta2/beta1 - 1.96*sqrt((se2/beta2)^2+(se1/beta1)^2)
  direct_percentage_upper <- beta2/beta1 + 1.96*sqrt((se2/beta2)^2+(se1/beta1)^2)
  
  
  CombineMed_percentage_low <- (1-direct_percentage_upper) * 100
  CombineMed_percentage_upper <- (1-direct_percentage_low) * 100
  
  res <- data.frame(CombineMed_percentage,CombineMed_percentage_low,CombineMed_percentage_upper)
  return(res)
}


result_education_ON_CAD
result_MV.IVW_education.9mediators_ON_CAD
CombineMed_cal_9mediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.097, se2 = 0.215) 
1-0.097/0.4790515 # 79.75166%

###### ________7：LifetimeSmoking + depression + BMI + lipid + glucose + SBP + DBP ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_LifetimeSmoking <- IV_LifetimeSmoking["SNP"]
SNP_IV_depression <- IV_depression["SNP"]
SNP_IV_BMI <- IV_BMI["SNP"]
SNP_IV_lipid <- IV_lipid["SNP"]
SNP_IV_glucose <- IV_glucose["SNP"]
SNP_IV_SBP <- IV_SBP["SNP"]
SNP_IV_DBP <- IV_DBP["SNP"]

SNP_IV_education.7mediators <- rbind(SNP_IV_educationStringent,
                                     SNP_IV_LifetimeSmoking,
                                     SNP_IV_depression,
                                     SNP_IV_BMI,
                                     SNP_IV_lipid,
                                     SNP_IV_glucose,
                                     SNP_IV_SBP,
                                     SNP_IV_DBP) # 2457SNP
ttt <- SNP_IV_education.7mediators
tttt <- ttt [duplicated(ttt),]
ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
SNP_IV_education.7mediators <- ttttt

tt <- clump_data(SNP_IV_education.7mediators)
SNP_IV_education.7mediators <- clump_data(SNP_IV_education.7mediators) # 679SNP
write.xlsx(SNP_IV_education.7mediators, file = "Output/MVMR_IV/SNP_IV_education.7mediators.xlsx")

exp_education_MVwith.7mediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.7mediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.7mediators$exposure <- "education"

exp_LifetimeSmoking_MVwith.education.7mediators <- format_data(
  GWAS_LifetimeSmoking,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "P") 
exp_LifetimeSmoking_MVwith.education.7mediators$outcome <- "LifetimeSmoking"

exp_depression_MVwith.education.7mediators <- format_data(
  GWAS_depression,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "LogOR",
  se_col = "StdErrLogOR",
  eaf_col = "Freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_depression_MVwith.education.7mediators$outcome <- "depression"

exp_BMI_MVwith.education.7mediators <- format_data(
  GWAS_BMI,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_BMI_MVwith.education.7mediators$outcome <- "BMI"

exp_lipid_MVwith.education.7mediators <- format_data(
  GWAS_lipid,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_lipid_MVwith.education.7mediators$outcome <- "lipid"

exp_glucose_MVwith.education.7mediators <- format_data(
  GWAS_glucose,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_glucose_MVwith.education.7mediators$outcome <- "glucose"

exp_SBP_MVwith.education.7mediators <- format_data(
  GWAS_SBP,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_SBP_MVwith.education.7mediators$outcome <- "SBP"

exp_DBP_MVwith.education.7mediators <- format_data(
  GWAS_DBP,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_DBP_MVwith.education.7mediators$outcome <- "DBP"

education.LifetimeSmoking <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = exp_LifetimeSmoking_MVwith.education.7mediators, action = 1)
education.depression <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = exp_depression_MVwith.education.7mediators, action = 1)
education.BMI <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = exp_BMI_MVwith.education.7mediators, action = 1)
education.lipid <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = exp_lipid_MVwith.education.7mediators, action = 1)
education.glucose <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = exp_glucose_MVwith.education.7mediators, action = 1)
education.SBP <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = exp_SBP_MVwith.education.7mediators, action = 1)
education.DBP <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = exp_DBP_MVwith.education.7mediators, action = 1)

out_CAD_MV.education.7mediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.7mediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.7mediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.7mediators, outcome_dat = out_CAD_MV.education.7mediators, action = 1)

education.LifetimeSmoking.dat <- education.LifetimeSmoking[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.LifetimeSmoking.dat) <- c("SNP", "beta.education", "beta.LifetimeSmoking", "se.education", "se.LifetimeSmoking")

education.depression.dat <- education.depression[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.depression.dat) <- c("SNP", "beta.depression", "se.depression")

education.BMI.dat <- education.BMI[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.BMI.dat) <- c("SNP", "beta.BMI", "se.BMI")

education.lipid.dat <- education.lipid[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.lipid.dat) <- c("SNP", "beta.lipid", "se.lipid")

education.glucose.dat <- education.glucose[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.glucose.dat) <- c("SNP", "beta.glucose", "se.glucose")

education.SBP.dat <- education.SBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.SBP.dat) <- c("SNP", "beta.SBP", "se.SBP")

education.DBP.dat <- education.DBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.DBP.dat) <- c("SNP", "beta.DBP", "se.DBP")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.LifetimeSmoking.dat, education.depression.dat, by = "SNP")
mvmr.dat2 <- merge(mvmr.dat1, education.BMI.dat, by = "SNP")
mvmr.dat3 <- merge(mvmr.dat2, education.lipid.dat, by = "SNP")
mvmr.dat4 <- merge(mvmr.dat3, education.glucose.dat, by = "SNP")
mvmr.dat5 <- merge(mvmr.dat4, education.SBP.dat, by = "SNP")
mvmr.dat6 <- merge(mvmr.dat5, education.DBP.dat, by = "SNP")
MVdat_education.7mediators_ON_CAD <- merge(mvmr.dat6, education.CAD.dat, by = "SNP")# 374SNP
write.xlsx(MVdat_education.7mediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.7mediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.7mediators_ON_CAD[,c("beta.education", "beta.LifetimeSmoking","beta.depression", "beta.BMI", "beta.lipid", "beta.glucose", "beta.SBP", "beta.DBP")])
bxse = as.matrix(MVdat_education.7mediators_ON_CAD[,c("se.education", "se.LifetimeSmoking","se.depression", "se.BMI", "se.lipid", "se.glucose", "se.SBP", "se.DBP")])
by = as.vector(MVdat_education.7mediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.7mediators_ON_CAD$se.CAD)

MVdatForm_education.7mediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "LifetimeSmoking","depression", "BMI", "lipid", "glucose", "SBP", "DBP"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.7mediators_ON_CAD <- mr_mvivw(MVdatForm_education.7mediators_ON_CAD)
tt <- mvmr(MVdatForm_education.7mediators_ON_CAD)#F-statistic: 15.43 on 8 and 366 DF

# MV MR-Egger
result_MV.Egger_education.7mediators_ON_CAD <- mr_mvegger(MVdatForm_education.7mediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.7mediators_ON_CAD <- mr_mvlasso(MVdatForm_education.7mediators_ON_CAD)

# MV Median
result_MV.Median_education.7mediators_ON_CAD <- mr_mvmedian(MVdatForm_education.7mediators_ON_CAD)



F_education.7mediators <- strength_mvmr(r_input = MVdatForm_education.7mediators_ON_CAD, gencov = 0)


mv_hete_education.7mediators <- pleiotropy_mvmr(r_input = MVdatForm_education.7mediators_ON_CAD, gencov = 0)


result_education_ON_CAD
result_MV.IVW_education.7mediators_ON_CAD
CombineMed_cal_7mediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.174, se2 = 0.196) 


###### ___________________________________________________________________________________________________________________________________ ######
###### ________Health behaviors2：television + LifetimeSmoking ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_television <- IV_television["SNP"]
SNP_IV_LifetimeSmoking <- IV_LifetimeSmoking["SNP"]

SNP_IV_education.BehaviorMediators <- rbind(SNP_IV_educationStringent,
                                            SNP_IV_television,
                                            SNP_IV_LifetimeSmoking) # 635SNP
ttt <- SNP_IV_education.BehaviorMediators
tttt <- ttt [duplicated(ttt),]
ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
SNP_IV_education.BehaviorMediators <- ttttt

SNP_IV_education.BehaviorMediators <- clump_data(SNP_IV_education.BehaviorMediators) # 390SNP
write.xlsx(SNP_IV_education.BehaviorMediators, file = "Output/MVMR_IV/SNP_IV_education.BehaviorMediators.xlsx")

exp_education_MVwith.BehaviorMediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.BehaviorMediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.BehaviorMediators$exposure <- "education"

exp_television_MVwith.education.BehaviorMediators <- format_data(
  GWAS_television,
  snps = SNP_IV_education.BehaviorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_television_MVwith.education.BehaviorMediators$outcome <- "television"

exp_LifetimeSmoking_MVwith.education.BehaviorMediators <- format_data(
  GWAS_LifetimeSmoking,
  snps = SNP_IV_education.BehaviorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "P") 
exp_LifetimeSmoking_MVwith.education.BehaviorMediators$outcome <- "LifetimeSmoking"

education.television <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorMediators, outcome_dat = exp_television_MVwith.education.BehaviorMediators, action = 1)
education.LifetimeSmoking <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorMediators, outcome_dat = exp_LifetimeSmoking_MVwith.education.BehaviorMediators, action = 1)

out_CAD_MV.education.BehaviorMediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.BehaviorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.BehaviorMediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorMediators, outcome_dat = out_CAD_MV.education.BehaviorMediators, action = 1)

education.television.dat <- education.television[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.television.dat) <- c("SNP", "beta.education", "beta.television", "se.education", "se.television")

education.LifetimeSmoking.dat <- education.LifetimeSmoking[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.LifetimeSmoking.dat) <- c("SNP", "beta.LifetimeSmoking", "se.LifetimeSmoking")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.television.dat, education.LifetimeSmoking.dat, by = "SNP")
MVdat_education.BehaviorMediators_ON_CAD <- merge(mvmr.dat1, education.CAD.dat, by = "SNP")# 373SNP
write.xlsx(MVdat_education.BehaviorMediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.BehaviorMediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.BehaviorMediators_ON_CAD[,c("beta.education", "beta.television","beta.LifetimeSmoking")])
bxse = as.matrix(MVdat_education.BehaviorMediators_ON_CAD[,c("se.education", "se.television","se.LifetimeSmoking")])
by = as.vector(MVdat_education.BehaviorMediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.BehaviorMediators_ON_CAD$se.CAD)

MVdatForm_education.BehaviorMediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "television","LifetimeSmoking"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.BehaviorMediators_ON_CAD <- mr_mvivw(MVdatForm_education.BehaviorMediators_ON_CAD)
tt <- mvmr(MVdatForm_education.BehaviorMediators_ON_CAD)#F-statistic: 31.63 on 3 and 370 DF

# MV MR-Egger
result_MV.Egger_education.BehaviorMediators_ON_CAD <- mr_mvegger(MVdatForm_education.BehaviorMediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.BehaviorMediators_ON_CAD <- mr_mvlasso(MVdatForm_education.BehaviorMediators_ON_CAD)

# MV Median
result_MV.Median_education.BehaviorMediators_ON_CAD <- mr_mvmedian(MVdatForm_education.BehaviorMediators_ON_CAD)




F_education.BehaviorMediators <- strength_mvmr(r_input = MVdatForm_education.BehaviorMediators_ON_CAD, gencov = 0)


mv_hete_education.BehaviorMediators <- pleiotropy_mvmr(r_input = MVdatForm_education.BehaviorMediators_ON_CAD, gencov = 0)

result_education_ON_CAD
result_MV.IVW_education.BehaviorMediators_ON_CAD
CombineMed_cal_BehaviorMediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.251, se2 = 0.101) 


###### ________Psychological Health2：WBS + depression ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_WBS <- IV_WBS["SNP"]
SNP_IV_depression <- IV_depression["SNP"]

SNP_IV_education.PsychologicalMediators <- rbind(SNP_IV_educationStringent,
                                                 SNP_IV_WBS,
                                                 SNP_IV_depression) # 580SNP
ttt <- SNP_IV_education.PsychologicalMediators
tttt <- ttt [duplicated(ttt),]
ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
SNP_IV_education.PsychologicalMediators <- ttttt

SNP_IV_education.PsychologicalMediators <- clump_data(SNP_IV_education.PsychologicalMediators) # 370SNP
write.xlsx(SNP_IV_education.PsychologicalMediators, file = "Output/MVMR_IV/SNP_IV_education.PsychologicalMediators.xlsx")

exp_education_MVwith.PsychologicalMediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.PsychologicalMediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.PsychologicalMediators$exposure <- "education"

exp_WBS_MVwith.education.PsychologicalMediators <- format_data(
  GWAS_WBS,
  snps = SNP_IV_education.PsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "PVAL") 
exp_WBS_MVwith.education.PsychologicalMediators$outcome <- "WBS"

exp_depression_MVwith.education.PsychologicalMediators <- format_data(
  GWAS_depression,
  snps = SNP_IV_education.PsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "LogOR",
  se_col = "StdErrLogOR",
  eaf_col = "Freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_depression_MVwith.education.PsychologicalMediators$outcome <- "depression"

education.WBS <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalMediators, outcome_dat = exp_WBS_MVwith.education.PsychologicalMediators, action = 1)
education.depression <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalMediators, outcome_dat = exp_depression_MVwith.education.PsychologicalMediators, action = 1)

out_CAD_MV.education.PsychologicalMediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.PsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.PsychologicalMediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalMediators, outcome_dat = out_CAD_MV.education.PsychologicalMediators, action = 1)

education.WBS.dat <- education.WBS[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.WBS.dat) <- c("SNP", "beta.education", "beta.WBS", "se.education", "se.WBS")

education.depression.dat <- education.depression[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.depression.dat) <- c("SNP", "beta.depression", "se.depression")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.WBS.dat, education.depression.dat, by = "SNP")
MVdat_education.PsychologicalMediators_ON_CAD <- merge(mvmr.dat1, education.CAD.dat, by = "SNP")# 313SNP
write.xlsx(MVdat_education.PsychologicalMediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.PsychologicalMediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.PsychologicalMediators_ON_CAD[,c("beta.education", "beta.WBS","beta.depression")])
bxse = as.matrix(MVdat_education.PsychologicalMediators_ON_CAD[,c("se.education", "se.WBS","se.depression")])
by = as.vector(MVdat_education.PsychologicalMediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.PsychologicalMediators_ON_CAD$se.CAD)

MVdatForm_education.PsychologicalMediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "WBS","depression"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.PsychologicalMediators_ON_CAD <- mr_mvivw(MVdatForm_education.PsychologicalMediators_ON_CAD)
tt <- mvmr(MVdatForm_education.PsychologicalMediators_ON_CAD)#F-statistic: 19.75 on 3 and 310 DF

# MV MR-Egger
result_MV.Egger_education.PsychologicalMediators_ON_CAD <- mr_mvegger(MVdatForm_education.PsychologicalMediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.PsychologicalMediators_ON_CAD <- mr_mvlasso(MVdatForm_education.PsychologicalMediators_ON_CAD)

# MV Median
result_MV.Median_education.PsychologicalMediators_ON_CAD <- mr_mvmedian(MVdatForm_education.PsychologicalMediators_ON_CAD)




F_education.PsychologicalMediators <- strength_mvmr(r_input = MVdatForm_education.PsychologicalMediators_ON_CAD, gencov = 0)


mv_hete_education.PsychologicalMediators <- pleiotropy_mvmr(r_input = MVdatForm_education.PsychologicalMediators_ON_CAD, gencov = 0)

result_education_ON_CAD
result_MV.IVW_education.PsychologicalMediators_ON_CAD
CombineMed_cal_PsychologicalMediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.423, se2 = 0.072) 


###### ________Health factors5：BMI + lipid + glucose + SBP + DBP ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_BMI <- IV_BMI["SNP"]
SNP_IV_lipid <- IV_lipid["SNP"]
SNP_IV_glucose <- IV_glucose["SNP"]
SNP_IV_SBP <- IV_SBP["SNP"]
SNP_IV_DBP <- IV_DBP["SNP"]

SNP_IV_education.FactorMediators <- rbind(SNP_IV_educationStringent,
                                          SNP_IV_BMI,
                                          SNP_IV_lipid,
                                          SNP_IV_glucose,
                                          SNP_IV_SBP,
                                          SNP_IV_DBP) # 2229SNP
ttt <- SNP_IV_education.FactorMediators
tttt <- ttt [duplicated(ttt),]
ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
SNP_IV_education.FactorMediators <- ttttt

SNP_IV_education.FactorMediators <- clump_data(SNP_IV_education.FactorMediators) # 666SNP
write.xlsx(SNP_IV_education.FactorMediators, file = "Output/MVMR_IV/SNP_IV_education.FactorMediators.xlsx")

exp_education_MVwith.FactorMediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.FactorMediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.FactorMediators$exposure <- "education"

exp_BMI_MVwith.education.FactorMediators <- format_data(
  GWAS_BMI,
  snps = SNP_IV_education.FactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_BMI_MVwith.education.FactorMediators$outcome <- "BMI"

exp_lipid_MVwith.education.FactorMediators <- format_data(
  GWAS_lipid,
  snps = SNP_IV_education.FactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_lipid_MVwith.education.FactorMediators$outcome <- "lipid"

exp_glucose_MVwith.education.FactorMediators <- format_data(
  GWAS_glucose,
  snps = SNP_IV_education.FactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_glucose_MVwith.education.FactorMediators$outcome <- "glucose"

exp_SBP_MVwith.education.FactorMediators <- format_data(
  GWAS_SBP,
  snps = SNP_IV_education.FactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_SBP_MVwith.education.FactorMediators$outcome <- "SBP"

exp_DBP_MVwith.education.FactorMediators <- format_data(
  GWAS_DBP,
  snps = SNP_IV_education.FactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_DBP_MVwith.education.FactorMediators$outcome <- "DBP"

education.BMI <- harmonise_data(exposure_dat = exp_education_MVwith.FactorMediators, outcome_dat = exp_BMI_MVwith.education.FactorMediators, action = 1)
education.lipid <- harmonise_data(exposure_dat = exp_education_MVwith.FactorMediators, outcome_dat = exp_lipid_MVwith.education.FactorMediators, action = 1)
education.glucose <- harmonise_data(exposure_dat = exp_education_MVwith.FactorMediators, outcome_dat = exp_glucose_MVwith.education.FactorMediators, action = 1)
education.SBP <- harmonise_data(exposure_dat = exp_education_MVwith.FactorMediators, outcome_dat = exp_SBP_MVwith.education.FactorMediators, action = 1)
education.DBP <- harmonise_data(exposure_dat = exp_education_MVwith.FactorMediators, outcome_dat = exp_DBP_MVwith.education.FactorMediators, action = 1)

out_CAD_MV.education.FactorMediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.FactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.FactorMediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.FactorMediators, outcome_dat = out_CAD_MV.education.FactorMediators, action = 1)

education.BMI.dat <- education.BMI[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.BMI.dat) <- c("SNP", "beta.education", "beta.BMI", "se.education", "se.BMI")

education.lipid.dat <- education.lipid[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.lipid.dat) <- c("SNP", "beta.lipid", "se.lipid")

education.glucose.dat <- education.glucose[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.glucose.dat) <- c("SNP", "beta.glucose", "se.glucose")

education.SBP.dat <- education.SBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.SBP.dat) <- c("SNP", "beta.SBP", "se.SBP")

education.DBP.dat <- education.DBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.DBP.dat) <- c("SNP", "beta.DBP", "se.DBP")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.BMI.dat, education.lipid.dat, by = "SNP")
mvmr.dat2 <- merge(mvmr.dat1, education.glucose.dat, by = "SNP")
mvmr.dat3 <- merge(mvmr.dat2, education.SBP.dat, by = "SNP")
mvmr.dat4 <- merge(mvmr.dat3, education.DBP.dat, by = "SNP")
MVdat_education.FactorMediators_ON_CAD <- merge(mvmr.dat4, education.CAD.dat, by = "SNP")# 376SNP
write.xlsx(MVdat_education.FactorMediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.FactorMediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.FactorMediators_ON_CAD[,c("beta.education", "beta.BMI","beta.lipid","beta.glucose","beta.SBP","beta.DBP")])
bxse = as.matrix(MVdat_education.FactorMediators_ON_CAD[,c("se.education", "se.BMI","se.lipid","se.glucose","se.SBP","se.DBP")])
by = as.vector(MVdat_education.FactorMediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.FactorMediators_ON_CAD$se.CAD)

MVdatForm_education.FactorMediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "BMI","lipid","glucose","SBP","DBP"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.FactorMediators_ON_CAD <- mr_mvivw(MVdatForm_education.FactorMediators_ON_CAD)
tt <- mvmr(MVdatForm_education.FactorMediators_ON_CAD)#F-statistic: 19.62 on 6 and 370 DF

# MV MR-Egger
result_MV.Egger_education.FactorMediators_ON_CAD <- mr_mvegger(MVdatForm_education.FactorMediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.FactorMediators_ON_CAD <- mr_mvlasso(MVdatForm_education.FactorMediators_ON_CAD)

# MV Median
result_MV.Median_education.FactorMediators_ON_CAD <- mr_mvmedian(MVdatForm_education.FactorMediators_ON_CAD)




F_education.FactorMediators <- strength_mvmr(r_input = MVdatForm_education.FactorMediators_ON_CAD, gencov = 0)


mv_hete_education.FactorMediators <- pleiotropy_mvmr(r_input = MVdatForm_education.FactorMediators_ON_CAD, gencov = 0)

result_education_ON_CAD
result_MV.IVW_education.FactorMediators_ON_CAD
CombineMed_cal_FactorMediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.290, se2 = 0.136) 


###### ________Health behaviors2+Psychological Health2：television + LifetimeSmoking + WBS + depression ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_television <- IV_television["SNP"]
SNP_IV_LifetimeSmoking <- IV_LifetimeSmoking["SNP"]
SNP_IV_WBS <- IV_WBS["SNP"]
SNP_IV_depression <- IV_depression["SNP"]

SNP_IV_education.BehaviorPsychologicalMediators <- rbind(SNP_IV_educationStringent,
                                                         SNP_IV_television,
                                                         SNP_IV_LifetimeSmoking,
                                                         SNP_IV_WBS,
                                                         SNP_IV_depression) # 867SNP
ttt <- SNP_IV_education.BehaviorPsychologicalMediators
tttt <- ttt [duplicated(ttt),]
ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
SNP_IV_education.BehaviorPsychologicalMediators <- ttttt

SNP_IV_education.BehaviorPsychologicalMediators <- clump_data(SNP_IV_education.BehaviorPsychologicalMediators) # 415SNP
write.xlsx(SNP_IV_education.BehaviorPsychologicalMediators, file = "Output/MVMR_IV/SNP_IV_education.BehaviorPsychologicalMediators.xlsx")

exp_education_MVwith.BehaviorPsychologicalMediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.BehaviorPsychologicalMediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.BehaviorPsychologicalMediators$exposure <- "education"

exp_television_MVwith.education.BehaviorPsychologicalMediators <- format_data(
  GWAS_television,
  snps = SNP_IV_education.BehaviorPsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_television_MVwith.education.BehaviorPsychologicalMediators$outcome <- "television"

exp_LifetimeSmoking_MVwith.education.BehaviorPsychologicalMediators <- format_data(
  GWAS_LifetimeSmoking,
  snps = SNP_IV_education.BehaviorPsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "P") 
exp_LifetimeSmoking_MVwith.education.BehaviorPsychologicalMediators$outcome <- "LifetimeSmoking"

exp_WBS_MVwith.education.BehaviorPsychologicalMediators <- format_data(
  GWAS_WBS,
  snps = SNP_IV_education.BehaviorPsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "PVAL") 
exp_WBS_MVwith.education.BehaviorPsychologicalMediators$outcome <- "WBS"

exp_depression_MVwith.education.BehaviorPsychologicalMediators <- format_data(
  GWAS_depression,
  snps = SNP_IV_education.BehaviorPsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "LogOR",
  se_col = "StdErrLogOR",
  eaf_col = "Freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_depression_MVwith.education.BehaviorPsychologicalMediators$outcome <- "depression"

education.television <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorPsychologicalMediators, outcome_dat = exp_television_MVwith.education.BehaviorPsychologicalMediators, action = 1)
education.LifetimeSmoking <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorPsychologicalMediators, outcome_dat = exp_LifetimeSmoking_MVwith.education.BehaviorPsychologicalMediators, action = 1)
education.WBS <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorPsychologicalMediators, outcome_dat = exp_WBS_MVwith.education.BehaviorPsychologicalMediators, action = 1)
education.depression <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorPsychologicalMediators, outcome_dat = exp_depression_MVwith.education.BehaviorPsychologicalMediators, action = 1)

out_CAD_MV.education.BehaviorPsychologicalMediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.BehaviorPsychologicalMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.BehaviorPsychologicalMediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorPsychologicalMediators, outcome_dat = out_CAD_MV.education.BehaviorPsychologicalMediators, action = 1)

education.television.dat <- education.television[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.television.dat) <- c("SNP", "beta.education", "beta.television", "se.education", "se.television")

education.LifetimeSmoking.dat <- education.LifetimeSmoking[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.LifetimeSmoking.dat) <- c("SNP", "beta.LifetimeSmoking", "se.LifetimeSmoking")

education.WBS.dat <- education.WBS[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.WBS.dat) <- c("SNP", "beta.WBS", "se.WBS")

education.depression.dat <- education.depression[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.depression.dat) <- c("SNP", "beta.depression", "se.depression")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.television.dat, education.LifetimeSmoking.dat, by = "SNP")
mvmr.dat2 <- merge(mvmr.dat1, education.WBS.dat, by = "SNP")
mvmr.dat3 <- merge(mvmr.dat2, education.depression.dat, by = "SNP")
MVdat_education.BehaviorPsychologicalMediators_ON_CAD <- merge(mvmr.dat3, education.CAD.dat, by = "SNP")# 255SNP
write.xlsx(MVdat_education.BehaviorPsychologicalMediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.BehaviorPsychologicalMediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.BehaviorPsychologicalMediators_ON_CAD[,c("beta.education", "beta.television","beta.LifetimeSmoking","beta.WBS","beta.depression")])
bxse = as.matrix(MVdat_education.BehaviorPsychologicalMediators_ON_CAD[,c("se.education", "se.television","se.LifetimeSmoking","se.WBS","se.depression")])
by = as.vector(MVdat_education.BehaviorPsychologicalMediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.BehaviorPsychologicalMediators_ON_CAD$se.CAD)

MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "television","LifetimeSmoking","WBS","depression"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.BehaviorPsychologicalMediators_ON_CAD <- mr_mvivw(MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD)
tt <- mvmr(MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD)#F-statistic: 10.37 on 5 and 250 DF

# MV MR-Egger
result_MV.Egger_education.BehaviorPsychologicalMediators_ON_CAD <- mr_mvegger(MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.BehaviorPsychologicalMediators_ON_CAD <- mr_mvlasso(MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD)

# MV Median
result_MV.Median_education.BehaviorPsychologicalMediators_ON_CAD <- mr_mvmedian(MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD)



F_education.BehaviorPsychologicalMediators <- strength_mvmr(r_input = MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD, gencov = 0)


mv_hete_education.BehaviorPsychologicalMediators <- pleiotropy_mvmr(r_input = MVdatForm_education.BehaviorPsychologicalMediators_ON_CAD, gencov = 0)

result_education_ON_CAD
result_MV.IVW_education.BehaviorPsychologicalMediators_ON_CAD
CombineMed_cal_BehaviorPsychologicalMediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.206, se2 = 0.130) 


###### ________Health behaviors2+Health factors5：television + LifetimeSmoking + BMI + lipid + glucose + SBP + DBP ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_television <- IV_television["SNP"]
SNP_IV_LifetimeSmoking <- IV_LifetimeSmoking["SNP"]
SNP_IV_BMI <- IV_BMI["SNP"]
SNP_IV_lipid <- IV_lipid["SNP"]
SNP_IV_glucose <- IV_glucose["SNP"]
SNP_IV_SBP <- IV_SBP["SNP"]
SNP_IV_DBP <- IV_DBP["SNP"]

SNP_IV_education.BehaviorFactorMediators <- rbind(SNP_IV_educationStringent,
                                                  SNP_IV_television,
                                                  SNP_IV_LifetimeSmoking,
                                                  SNP_IV_BMI,
                                                  SNP_IV_lipid,
                                                  SNP_IV_glucose,
                                                  SNP_IV_SBP,
                                                  SNP_IV_DBP) # 2516SNP
ttt <- SNP_IV_education.BehaviorFactorMediators
tttt <- ttt [duplicated(ttt),]
ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
SNP_IV_education.BehaviorFactorMediators <- ttttt

SNP_IV_education.BehaviorFactorMediators <- clump_data(SNP_IV_education.BehaviorFactorMediators) # 671SNP
write.xlsx(SNP_IV_education.BehaviorFactorMediators, file = "Output/MVMR_IV/SNP_IV_education.BehaviorFactorMediators.xlsx")

exp_education_MVwith.BehaviorFactorMediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.BehaviorFactorMediators$exposure <- "education"

exp_television_MVwith.education.BehaviorFactorMediators <- format_data(
  GWAS_television,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_television_MVwith.education.BehaviorFactorMediators$outcome <- "television"

exp_LifetimeSmoking_MVwith.education.BehaviorFactorMediators <- format_data(
  GWAS_LifetimeSmoking,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "P") 
exp_LifetimeSmoking_MVwith.education.BehaviorFactorMediators$outcome <- "LifetimeSmoking"

exp_BMI_MVwith.education.BehaviorFactorMediators <- format_data(
  GWAS_BMI,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_BMI_MVwith.education.BehaviorFactorMediators$outcome <- "BMI"

exp_lipid_MVwith.education.BehaviorFactorMediators <- format_data(
  GWAS_lipid,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_lipid_MVwith.education.BehaviorFactorMediators$outcome <- "lipid"

exp_glucose_MVwith.education.BehaviorFactorMediators <- format_data(
  GWAS_glucose,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_glucose_MVwith.education.BehaviorFactorMediators$outcome <- "glucose"

exp_SBP_MVwith.education.BehaviorFactorMediators <- format_data(
  GWAS_SBP,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_SBP_MVwith.education.BehaviorFactorMediators$outcome <- "SBP"

exp_DBP_MVwith.education.BehaviorFactorMediators <- format_data(
  GWAS_DBP,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_DBP_MVwith.education.BehaviorFactorMediators$outcome <- "DBP"

education.television <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = exp_television_MVwith.education.BehaviorFactorMediators, action = 1)
education.LifetimeSmoking <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = exp_LifetimeSmoking_MVwith.education.BehaviorFactorMediators, action = 1)
education.BMI <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = exp_BMI_MVwith.education.BehaviorFactorMediators, action = 1)
education.lipid <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = exp_lipid_MVwith.education.BehaviorFactorMediators, action = 1)
education.glucose <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = exp_glucose_MVwith.education.BehaviorFactorMediators, action = 1)
education.SBP <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = exp_SBP_MVwith.education.BehaviorFactorMediators, action = 1)
education.DBP <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = exp_DBP_MVwith.education.BehaviorFactorMediators, action = 1)

out_CAD_MV.education.BehaviorFactorMediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.BehaviorFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.BehaviorFactorMediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.BehaviorFactorMediators, outcome_dat = out_CAD_MV.education.BehaviorFactorMediators, action = 1)

education.television.dat <- education.television[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.television.dat) <- c("SNP", "beta.education", "beta.television", "se.education", "se.television")

education.LifetimeSmoking.dat <- education.LifetimeSmoking[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.LifetimeSmoking.dat) <- c("SNP", "beta.LifetimeSmoking", "se.LifetimeSmoking")

education.BMI.dat <- education.BMI[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.BMI.dat) <- c("SNP", "beta.BMI", "se.BMI")

education.lipid.dat <- education.lipid[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.lipid.dat) <- c("SNP", "beta.lipid", "se.lipid")

education.glucose.dat <- education.glucose[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.glucose.dat) <- c("SNP", "beta.glucose", "se.glucose")

education.SBP.dat <- education.SBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.SBP.dat) <- c("SNP", "beta.SBP", "se.SBP")

education.DBP.dat <- education.DBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.DBP.dat) <- c("SNP", "beta.DBP", "se.DBP")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.television.dat, education.LifetimeSmoking.dat, by = "SNP")
mvmr.dat2 <- merge(mvmr.dat1, education.BMI.dat, by = "SNP")
mvmr.dat3 <- merge(mvmr.dat2, education.lipid.dat, by = "SNP")
mvmr.dat4 <- merge(mvmr.dat3, education.glucose.dat, by = "SNP")
mvmr.dat5 <- merge(mvmr.dat4, education.SBP.dat, by = "SNP")
mvmr.dat6 <- merge(mvmr.dat5, education.DBP.dat, by = "SNP")
MVdat_education.BehaviorFactorMediators_ON_CAD <- merge(mvmr.dat6, education.CAD.dat, by = "SNP")# 362SNP
write.xlsx(MVdat_education.BehaviorFactorMediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.BehaviorFactorMediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.BehaviorFactorMediators_ON_CAD[,c("beta.education", "beta.television","beta.LifetimeSmoking","beta.BMI","beta.lipid","beta.glucose","beta.SBP","beta.DBP")])
bxse = as.matrix(MVdat_education.BehaviorFactorMediators_ON_CAD[,c("se.education", "se.television","se.LifetimeSmoking","se.BMI","se.lipid","se.glucose","se.SBP","se.DBP")])
by = as.vector(MVdat_education.BehaviorFactorMediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.BehaviorFactorMediators_ON_CAD$se.CAD)

MVdatForm_education.BehaviorFactorMediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "television","LifetimeSmoking","BMI","lipid","glucose","SBP","DBP"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.BehaviorFactorMediators_ON_CAD <- mr_mvivw(MVdatForm_education.BehaviorFactorMediators_ON_CAD)
tt <- mvmr(MVdatForm_education.BehaviorFactorMediators_ON_CAD)#F-statistic: 21.25 on 8 and 354 DF

# MV MR-Egger
result_MV.Egger_education.BehaviorFactorMediators_ON_CAD <- mr_mvegger(MVdatForm_education.BehaviorFactorMediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.BehaviorFactorMediators_ON_CAD <- mr_mvlasso(MVdatForm_education.BehaviorFactorMediators_ON_CAD)

# MV Median
result_MV.Median_education.BehaviorFactorMediators_ON_CAD <- mr_mvmedian(MVdatForm_education.BehaviorFactorMediators_ON_CAD)



F_education.BehaviorFactorMediators <- strength_mvmr(r_input = MVdatForm_education.BehaviorFactorMediators_ON_CAD, gencov = 0)


mv_hete_education.BehaviorFactorMediators <- pleiotropy_mvmr(r_input = MVdatForm_education.BehaviorFactorMediators_ON_CAD, gencov = 0)

result_education_ON_CAD
result_MV.IVW_education.BehaviorFactorMediators_ON_CAD
CombineMed_cal_BehaviorFactorMediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.115, se2 = 0.179) 


###### ________Psychological Health2+Health factors5：WBS + depression + BMI + lipid + glucose + SBP + DBP ######
SNP_IV_educationStringent <- IV_educationStringent["SNP"]
SNP_IV_WBS <- IV_WBS["SNP"]
SNP_IV_depression <- IV_depression["SNP"]
SNP_IV_BMI <- IV_BMI["SNP"]
SNP_IV_lipid <- IV_lipid["SNP"]
SNP_IV_glucose <- IV_glucose["SNP"]
SNP_IV_SBP <- IV_SBP["SNP"]
SNP_IV_DBP <- IV_DBP["SNP"]

SNP_IV_education.PsychologicalFactorMediators <- rbind(SNP_IV_educationStringent,
                                                       SNP_IV_WBS,
                                                       SNP_IV_depression,
                                                       SNP_IV_BMI,
                                                       SNP_IV_lipid,
                                                       SNP_IV_glucose,
                                                       SNP_IV_SBP,
                                                       SNP_IV_DBP) 
ttt <- SNP_IV_education.PsychologicalFactorMediators
tttt <- ttt [duplicated(ttt),]
ttttt <- ttt [!ttt$SNP %in% tttt, ,drop = F]
SNP_IV_education.PsychologicalFactorMediators <- ttttt

SNP_IV_education.PsychologicalFactorMediators <- clump_data(SNP_IV_education.PsychologicalFactorMediators) # 674SNP
write.xlsx(SNP_IV_education.PsychologicalFactorMediators, file = "Output/MVMR_IV/SNP_IV_education.PsychologicalFactorMediators.xlsx")

exp_education_MVwith.PsychologicalFactorMediators <- format_data(
  GWAS_education,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_education_MVwith.PsychologicalFactorMediators$exposure <- "education"

exp_WBS_MVwith.education.PsychologicalFactorMediators <- format_data(
  GWAS_WBS,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "PVAL") 
exp_WBS_MVwith.education.PsychologicalFactorMediators$outcome <- "WBS"

exp_depression_MVwith.education.PsychologicalFactorMediators <- format_data(
  GWAS_depression,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "LogOR",
  se_col = "StdErrLogOR",
  eaf_col = "Freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P") 
exp_depression_MVwith.education.PsychologicalFactorMediators$outcome <- "depression"

exp_BMI_MVwith.education.PsychologicalFactorMediators <- format_data(
  GWAS_BMI,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_BMI_MVwith.education.PsychologicalFactorMediators$outcome <- "BMI"

exp_lipid_MVwith.education.PsychologicalFactorMediators <- format_data(
  GWAS_lipid,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_lipid_MVwith.education.PsychologicalFactorMediators$outcome <- "lipid"

exp_glucose_MVwith.education.PsychologicalFactorMediators <- format_data(
  GWAS_glucose,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "P") 
exp_glucose_MVwith.education.PsychologicalFactorMediators$outcome <- "glucose"

exp_SBP_MVwith.education.PsychologicalFactorMediators <- format_data(
  GWAS_SBP,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_SBP_MVwith.education.PsychologicalFactorMediators$outcome <- "SBP"

exp_DBP_MVwith.education.PsychologicalFactorMediators <- format_data(
  GWAS_DBP,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P") 
exp_DBP_MVwith.education.PsychologicalFactorMediators$outcome <- "DBP"


education.WBS <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = exp_WBS_MVwith.education.PsychologicalFactorMediators, action = 1)
education.depression <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = exp_depression_MVwith.education.PsychologicalFactorMediators, action = 1)
education.BMI <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = exp_BMI_MVwith.education.PsychologicalFactorMediators, action = 1)
education.lipid <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = exp_lipid_MVwith.education.PsychologicalFactorMediators, action = 1)
education.glucose <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = exp_glucose_MVwith.education.PsychologicalFactorMediators, action = 1)
education.SBP <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = exp_SBP_MVwith.education.PsychologicalFactorMediators, action = 1)
education.DBP <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = exp_DBP_MVwith.education.PsychologicalFactorMediators, action = 1)


out_CAD_MV.education.PsychologicalFactorMediators <- format_data(
  GWAS_CAD,
  snps = SNP_IV_education.PsychologicalFactorMediators$SNP,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
out_CAD_MV.education.PsychologicalFactorMediators$outcome <- "CAD"

education.CAD <- harmonise_data(exposure_dat = exp_education_MVwith.PsychologicalFactorMediators, outcome_dat = out_CAD_MV.education.PsychologicalFactorMediators, action = 1)

education.WBS.dat <- education.WBS[c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")]
colnames(education.WBS.dat) <- c("SNP", "beta.education", "beta.WBS", "se.education", "se.WBS")

education.depression.dat <- education.depression[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.depression.dat) <- c("SNP", "beta.depression", "se.depression")

education.BMI.dat <- education.BMI[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.BMI.dat) <- c("SNP", "beta.BMI", "se.BMI")

education.lipid.dat <- education.lipid[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.lipid.dat) <- c("SNP", "beta.lipid", "se.lipid")

education.glucose.dat <- education.glucose[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.glucose.dat) <- c("SNP", "beta.glucose", "se.glucose")

education.SBP.dat <- education.SBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.SBP.dat) <- c("SNP", "beta.SBP", "se.SBP")

education.DBP.dat <- education.DBP[,c("SNP", "beta.outcome", "se.outcome")]
colnames(education.DBP.dat) <- c("SNP", "beta.DBP", "se.DBP")

education.CAD.dat <- education.CAD[c("SNP", "beta.outcome", "se.outcome")]
colnames(education.CAD.dat) <- c("SNP", "beta.CAD", "se.CAD")

mvmr.dat1 <- merge(education.WBS.dat, education.depression.dat, by = "SNP")
mvmr.dat2 <- merge(mvmr.dat1, education.BMI.dat, by = "SNP")
mvmr.dat3 <- merge(mvmr.dat2, education.lipid.dat, by = "SNP")
mvmr.dat4 <- merge(mvmr.dat3, education.glucose.dat, by = "SNP")
mvmr.dat5 <- merge(mvmr.dat4, education.SBP.dat, by = "SNP")
mvmr.dat6 <- merge(mvmr.dat5, education.DBP.dat, by = "SNP")
MVdat_education.PsychologicalFactorMediators_ON_CAD <- merge(mvmr.dat6, education.CAD.dat, by = "SNP")# 327SNP
write.xlsx(MVdat_education.PsychologicalFactorMediators_ON_CAD, file = "Output/MVMR_IV/MVdat_education.PsychologicalFactorMediators_ON_CAD.xlsx")

bx = as.matrix(MVdat_education.PsychologicalFactorMediators_ON_CAD[,c("beta.education", "beta.WBS","beta.depression","beta.BMI","beta.lipid","beta.glucose","beta.SBP","beta.DBP")])
bxse = as.matrix(MVdat_education.PsychologicalFactorMediators_ON_CAD[,c("se.education", "se.WBS","se.depression","se.BMI","se.lipid","se.glucose","se.SBP","se.DBP")])
by = as.vector(MVdat_education.PsychologicalFactorMediators_ON_CAD$beta.CAD)
byse = as.vector(MVdat_education.PsychologicalFactorMediators_ON_CAD$se.CAD)

MVdatForm_education.PsychologicalFactorMediators_ON_CAD <- mr_mvinput(
  bx = bx,
  bxse = bxse,
  by = by,
  byse = byse,
  correlation =matrix(),
  exposure = c("education", "WBS","depression","BMI","lipid","glucose","SBP","DBP"),outcome = "CAD")

# MV-IVW
result_MV.IVW_education.PsychologicalFactorMediators_ON_CAD <- mr_mvivw(MVdatForm_education.PsychologicalFactorMediators_ON_CAD)
tt <- mvmr(MVdatForm_education.PsychologicalFactorMediators_ON_CAD)#F-statistic: 12.98 on 8 and 319 DF

# MV MR-Egger
result_MV.Egger_education.PsychologicalFactorMediators_ON_CAD <- mr_mvegger(MVdatForm_education.PsychologicalFactorMediators_ON_CAD)

# MV MR-Lasso
result_MV.Lasso_education.PsychologicalFactorMediators_ON_CAD <- mr_mvlasso(MVdatForm_education.PsychologicalFactorMediators_ON_CAD)

# MV Median
result_MV.Median_education.PsychologicalFactorMediators_ON_CAD <- mr_mvmedian(MVdatForm_education.PsychologicalFactorMediators_ON_CAD)



F_education.PsychologicalFactorMediators <- strength_mvmr(r_input = MVdatForm_education.PsychologicalFactorMediators_ON_CAD, gencov = 0)


mv_hete_education.PsychologicalFactorMediators <- pleiotropy_mvmr(r_input = MVdatForm_education.PsychologicalFactorMediators_ON_CAD, gencov = 0)

result_education_ON_CAD
result_MV.IVW_education.PsychologicalFactorMediators_ON_CAD
CombineMed_cal_PsychologicalFactorMediators_IN_edu.CAD <- CombineMed_cal(beta1 = -0.4790515, se1 = 0.03212376, beta2 = -0.180, se2 = 0.149) 


