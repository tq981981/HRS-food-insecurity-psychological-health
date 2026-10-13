# 数据处理
library(tidyverse)
# 描述性统计
library(tableone)
# 缺失值分析
library(naniar)
# 回归模型
library(broom)
# 多重插补（后续使用）
library(mice)
# SEM（后续使用）
library(lavaan)
# 可视化
library(ggplot2)
library(psych)
library(lme4)
library(lmerTest)
library(broom)
library(broom.mixed)
library(mediation)
library(performance)
library(car)
library(gtsummary)
library(skimr)
library(flextable)
library(officer)
library(haven)

####导入数据集并进行基本的数据处理####
# 必须放在 read_csv() 之前
Sys.setenv(VROOM_CONNECTION_SIZE = "268435456")  # 256 MB

library(readr)

raw_file <- paste0(
  "/dataset/HRS/hcns2013_raw_merged_hrs2012_2020/",
  "hrs_hcns2013_raw_merged_2012_2020.csv.gz"
)

hrs_raw <- read_csv(
  raw_file,
  col_types = cols(.default = col_character()),
  na = character(),
  name_repair = "minimal",
  progress = TRUE,
  show_col_types = FALSE
)

# 核验
dim(hrs_raw)                         # 8073 × 10244
n_distinct(hrs_raw$HHIDPN)           # 8073
sum(duplicated(hrs_raw$HHIDPN))      # 0

# 基本核验
dim(hrs_raw)
# 应为：8073 × 10244

n_distinct(hrs_raw$HHIDPN)
# 应为：8073

sum(duplicated(hrs_raw$HHIDPN))
# 应为：0


library(haven)
library(dplyr)
library(stringr)

# 1. 从RAND文件读取2012年SES变量
rand_zip <- paste0(
  "/dataset/HRS/RAND HRS Longitudinal File 2022/",
  "randhrs1992_2022v1_SAS.zip"
)

rand_member <- "randhrs1992_2022v1.sas7bdat"
rand_dir <- file.path(tempdir(), "rand_hrs_2022")
dir.create(rand_dir, showWarnings = FALSE)

unzip(
  rand_zip,
  files = rand_member,
  exdir = rand_dir,
  overwrite = TRUE
)

rand_file <- file.path(rand_dir, rand_member)

rand_ses <- read_sas(
  rand_file,
  col_select = c(
    HHIDPN,
    H11ATOTB,  # 2012家庭总财富/总资产
    H11ITOT,   # 2012家庭总收入
    RAEDUC     # 统一教育分类
  )
) |>
  transmute(
    HHIDPN = str_pad(
      as.character(HHIDPN),
      width = 9,
      side = "left",
      pad = "0"
    ),
    wealth_raw_2012 = as.numeric(H11ATOTB),
    income_raw_2012 = as.numeric(H11ITOT),
    education_raw = as.integer(RAEDUC)
  )

# RAND文件中每个人必须唯一
stopifnot(sum(duplicated(rand_ses$HHIDPN)) == 0)


# 保证主数据ID仍为9位字符型
hrs_raw <- hrs_raw |>
  mutate(
    HHIDPN = str_pad(
      as.character(HHIDPN),
      width = 9,
      side = "left",
      pad = "0"
    )
  )

n_before <- nrow(hrs_raw)

hrs_ses <- hrs_raw |>
  left_join(rand_ses, by = "HHIDPN")

stopifnot(nrow(hrs_ses) == n_before)
stopifnot(sum(duplicated(hrs_ses$HHIDPN)) == 0) 


#将财富和收入分成三等分，教育分为三个等级
# 返回0、1、2分；缺失值保持缺失
make_tertile_score <- function(x) {
  q <- quantile(
    x,
    probs = c(1/3, 2/3),
    na.rm = TRUE,
    type = 2
  )
  
  if (q[1] >= q[2]) {
    stop("三等分切点重复，需检查变量分布")
  }
  
  case_when(
    is.na(x)  ~ NA_integer_,
    x <= q[1] ~ 0L,
    x <= q[2] ~ 1L,
    TRUE      ~ 2L
  )
}

hrs_ses <- hrs_ses |>
  mutate(
    # 财富：低、中、高
    wealth_tertile_score = make_tertile_score(wealth_raw_2012),
    
    # 收入：低、中、高
    income_tertile_score = make_tertile_score(income_raw_2012),
    
    # 教育：低于高中、高中/GED、部分大学及以上
    education_3cat = case_when(
      education_raw == 1L        ~ "Low",
      education_raw %in% 2:3     ~ "Middle",
      education_raw %in% 4:5     ~ "High",
      TRUE                       ~ NA_character_
    ),
    
    education_score = case_when(
      education_raw == 1L        ~ 0L,
      education_raw %in% 2:3     ~ 1L,
      education_raw %in% 4:5     ~ 2L,
      TRUE                       ~ NA_integer_
    ),
    
    education_3cat = factor(
      education_3cat,
      levels = c("Low", "Middle", "High"),
      ordered = TRUE
    )
  )

#构建三变量SES总分
hrs_ses <- hrs_ses |>
  mutate(
    ses_3var_complete =
      !is.na(wealth_tertile_score) &
      !is.na(income_tertile_score) &
      !is.na(education_score),
    
    # 仅三个组成变量均非缺失时计算，范围0–6
    ses_3var_score_0_6 = if_else(
      ses_3var_complete,
      wealth_tertile_score +
        income_tertile_score +
        education_score,
      NA_integer_
    )
  )


table(hrs_ses$wealth_tertile_score, useNA = "always")
table(hrs_ses$income_tertile_score, useNA = "always")
table(hrs_ses$education_3cat, useNA = "always")
table(hrs_ses$ses_3var_score_0_6, useNA = "always")

summary(hrs_ses$wealth_raw_2012)
summary(hrs_ses$income_raw_2012)

sum(hrs_ses$ses_3var_complete)

####处理基线SES变量缺失值####
hrs_ses = hrs_ses %>%
  filter(wealth_tertile_score != "")%>%   #189缺失
  filter(income_tertile_score != "")%>%   #189缺失
  filter(education_3cat != "") #17缺失
#📌 去掉基线SES三个子评分缺失的参与者，剩7882人

hrs_ses$ses_3cat <- case_when(
  hrs_ses$ses_3var_score_0_6 %in% c("0", "1", "2") ~ "Low SES",
  hrs_ses$ses_3var_score_0_6 %in% c("3", "4") ~ "Middle SES",
  hrs_ses$ses_3var_score_0_6 %in% c("5", "6") ~ "High SES",
  TRUE ~ NA_character_
)

####食品不安全FI--查看B1-B5原始回答和缺失情况####
fi_vars <- paste0(
  "hcns2013_questionnaire__hnb",
  1:5,
  "_13"
)

stopifnot(all(fi_vars %in% names(hrs_ses)))

fi_raw <- hrs_ses[fi_vars]

#查看每个变量的原始频数
for (v in fi_vars) {
  cat("\n", v, "\n")
  print(table(hrs_ses[[v]], useNA = "always"))
}

#汇总特殊缺失编码99、空字符和NA
data.frame(
  variable = fi_vars,
  n_99 = sapply(fi_raw, function(x) sum(x == "99", na.rm = TRUE)),
  n_blank = sapply(fi_raw, function(x) sum(x == "", na.rm = TRUE)),
  n_NA = sapply(fi_raw, function(x) sum(is.na(x)))
)

#查看五项原始回答组合的完整性
fi_missing <- sapply(
  fi_raw,
  function(x) is.na(x) | x == "" | x == "99"
)

table(
  missing_items = rowSums(fi_missing),
  useNA = "always"
)

#重点确认B3完整频率回答
table(
  hrs_ses$hcns2013_questionnaire__hnb3_13,
  useNA = "always"
)

#处理无效回答"99"
fi <- data.frame(
  lapply(
    hrs_ses[paste0("hcns2013_questionnaire__hnb", 1:5, "_13")],
    as.integer
  )
)

names(fi) <- paste0("fi_b", 1:5)

# 99为无效回答，转为缺失
fi[fi == 99] <- NA

# 加入工作数据，不覆盖原始变量
hrs_ses <- cbind(hrs_ses, fi)

lapply(fi, table, useNA = "always")


#统计一下所有参与者这五个问题的回答从全部缺失到全部都有的分类情况
fi_vars <- paste0(
  "hcns2013_questionnaire__hnb",
  1:5,
  "_13"
)

fi_raw <- hrs_ses[fi_vars]

fi_missing <- sapply(
  hrs_ses,
  function(x) is.na(x) | x == "" | x == "99"
)

n_answered <- 5 - rowSums(fi_missing)

tab <- table(factor(n_answered, levels = 0:5))

data.frame(
  answered_items = 0:5,
  missing_items = 5:0,
  n = as.integer(tab),
  percent = round(100 * as.integer(tab) / nrow(hrs_raw), 2)
)

#***********************************去掉FI_score的缺失值*********************************#
fi_vars <- paste0("fi_b", 1:5)

hrs_ses <- hrs_ses[
  complete.cases(hrs_ses[fi_vars]),
]
#还剩7388人

hrs_ses <- hrs_ses %>%    #计算总分
  mutate(
    fi_score_0_6 =
      ifelse(fi_b1 %in% c(1, 2), 1, 0) +
      ifelse(fi_b2 %in% c(1, 2), 1, 0) +
      ifelse(fi_b3 %in% c(1, 2), 2,
             ifelse(fi_b3 == 3, 1, 0)) +
      ifelse(fi_b4 == 1, 1, 0) +
      ifelse(fi_b5 == 1, 1, 0)
  )
table(hrs_ses$fi_score_0_6)

#***********************************按照USDA 6-item short form对FI_score进行分类***********************************#
hrs_ses <- hrs_ses %>%
  mutate(
    fi_USDA6_3cat = case_when(
      fi_score_0_6 %in% 0:1 ~ "Food secure/marginal",
      fi_score_0_6 %in% 2:4 ~ "Low food security",
      fi_score_0_6 %in% 5:6 ~ "Very low food security",
      TRUE ~ NA_character_
    ),
    
    fi_USDA6_binary = case_when(
      fi_score_0_6 %in% 0:1 ~ 0,
      fi_score_0_6 %in% 2:6 ~ 1,
      TRUE ~ NA_real_
    )
  )

#设置变量分类顺序
hrs_ses$fi_USDA6_3cat <- factor(
  hrs_ses$fi_USDA6_3cat,
  levels = c(
    "Food secure/marginal",
    "Low food security",
    "Very low food security"
  )
)



#检查分类情况
table(hrs_ses$fi_score_0_6, useNA = "always")
table(hrs_ses$fi_USDA6_3cat, useNA = "always")
table(hrs_ses$fi_USDA6_binary, useNA = "always")

####心理变量查看####
Sys.setenv(VROOM_CONNECTION_SIZE = "536870912")

cesd_raw <- read_csv(
  "/dataset/HRS/hcns2013_raw_merged_hrs2012_2020/hrs_hcns2013_raw_merged_2012_2020_v2.csv.gz",
  col_select = c(HHIDPN, matches("d_r__.*d11[0-7]$")),
  col_types = cols(.default = col_character())
)

hrs_ses <- hrs_ses %>%
  left_join(cesd_raw, by = "HHIDPN")

#原始编码中1=是、5=否，其他值均设为缺失
cesd_vars <- grep("d_r__.*d11[0-7]$", names(hrs_ses), value = TRUE)

for (v in cesd_vars) {
  hrs_ses[[v]] <- as.numeric(hrs_ses[[v]])
  hrs_ses[[v]][!hrs_ses[[v]] %in% c(1, 5)] <- NA
}

#计算CES-D的简单函数
cesd_score <- function(data, negative, positive) {
  ifelse(
    complete.cases(data[c(negative, positive)]),
    rowSums(data[negative] == 1) +
      rowSums(data[positive] == 5),
    NA
  )
}

#计算五个波次总分
for (x in c("n", "o", "p", "q", "r")) {
  
  year <- c(n = 2012, o = 2014, p = 2016, q = 2018, r = 2020)[x]
  prefix <- paste0("h", substr(year, 3, 4), "d_r__", x, "d")
  
  negative <- paste0(prefix, c(110, 111, 112, 114, 116, 117))
  positive <- paste0(prefix, c(113, 115))
  
  hrs_ses[[paste0("cesd_score_0_8_", year)]] <-
    cesd_score(hrs_ses, negative, positive)
}

# 计算2012年基线CES-D的均值和标准差
cesd_mean_2012 <- mean(
  hrs_ses$cesd_score_0_8_2012,
  na.rm = TRUE
)

cesd_sd_2012 <- sd(
  hrs_ses$cesd_score_0_8_2012,
  na.rm = TRUE
)
# 所有波次均使用2012年的均值和标准差进行标准化，以保证结果的可比性
hrs_ses$cesd_good_z_2012 <- 
  -(hrs_ses$cesd_score_0_8_2012 - cesd_mean_2012) / cesd_sd_2012

hrs_ses$cesd_good_z_2014 <-
  -(hrs_ses$cesd_score_0_8_2014 - cesd_mean_2012) / cesd_sd_2012

hrs_ses$cesd_good_z_2016 <-
  -(hrs_ses$cesd_score_0_8_2016 - cesd_mean_2012) / cesd_sd_2012

hrs_ses$cesd_good_z_2018 <-
  -(hrs_ses$cesd_score_0_8_2018 - cesd_mean_2012) / cesd_sd_2012

hrs_ses$cesd_good_z_2020 <-
  -(hrs_ses$cesd_score_0_8_2020 - cesd_mean_2012) / cesd_sd_2012

#***********************************去掉基线时期CESD缺失以及只有一次测量的参与者***********************************#
hrs_ses = hrs_ses %>%
  filter(!is.na(cesd_good_z_2012))
#还剩7106
cesd_score_vars <- c(
  "cesd_score_0_8_2012",
  "cesd_score_0_8_2014",
  "cesd_score_0_8_2016",
  "cesd_score_0_8_2018",
  "cesd_score_0_8_2020"
)

# 每人的有效CES-D测量次数
hrs_ses$cesd_measurement_n <- rowSums(
  !is.na(hrs_ses[cesd_score_vars])
)

table(hrs_ses$cesd_measurement_n)

#排除只有一次结果的人
hrs_ses <- hrs_ses %>%
  filter(
    !is.na(cesd_score_0_8_2012),
    cesd_measurement_n >= 2
  )
#在7106基础上再减去295

#***********************************主要变量处理后的数据情况，从这里导入***********************************#
save.image("/dataset/HRS/dietaryindex_standard_input/投稿——frontiers in nutrition/返修/code_revised_frontier.RData")
load("/dataset/HRS/dietaryindex_standard_input/投稿——frontiers in nutrition/返修/code_revised_frontier.RData")

####协变量####
#协变量使用2012年基线，包括年龄、性别、种族/族裔、婚姻/伴侣、吸烟、饮酒、BMI、自评健康、慢性病和保险。
#从RAND文件读取
raw_covar <- c(
  "R11AGEY_E", "RAGENDER", "RARACEM", "RAHISPAN",
  "R11MSTAT", "R11SMOKEN", "R11DRINK", "R11BMI", "R11SHLT",
  "R11HIBPE", "R11DIABE", "R11CANCRE", "R11LUNGE",
  "R11HEARTE", "R11STROKE", "R11ARTHRE", "R11PSYCHE",
  "R11HIGOV", "R11HIOTHP"
)
covar <- read_sas(
  rand_file,
  col_select = c(HHIDPN, all_of(raw_covar))
)
#统一ID
hrs_ses$HHIDPN <- sprintf("%09.0f", as.numeric(hrs_ses$HHIDPN))
covar$HHIDPN <- sprintf("%09.0f", as.numeric(covar$HHIDPN))

hrs_ses <- hrs_ses %>%
  left_join(covar, by = "HHIDPN")

#慢性病计数
disease_vars <- c(
  "R11HIBPE", #高血压
  "R11DIABE", #糖尿病
  "R11CANCRE",#癌症或恶性肿瘤
  "R11LUNGE", #慢性肺部疾病，如慢性支气管炎或肺气肿
  "R11HEARTE",#心脏病，包括冠心病、心绞痛、心力衰竭等
  "R11STROKE",#卒中或短暂性脑缺血发作
  "R11ARTHRE" #关节炎或风湿病
)

hrs_ses$chronic_condition_count_2012 <-
  rowSums(hrs_ses[disease_vars], na.rm = FALSE)


#R11PSYCHE单独保留。主模型的慢性病计数暂不纳入精神/心理疾病，因为结局就是CES-D，纳入可能造成过度调整；
#***********************************psych_history_2012用于敏感性分析。***********************************#
hrs_ses$psych_history_2012 <- hrs_ses$R11PSYCHE

#对协变量重新编码赋值处理
hrs_ses <- hrs_ses %>%
  mutate(
    age_2012 = R11AGEY_E,
    
    sex_1male_0female = case_when(
      RAGENDER == 1 ~ 1,
      RAGENDER == 2 ~ 0,
      TRUE ~ NA_real_ #对于缺失者，变为NA
    ),
    
    race_ethnicity_4cat = case_when(
      RAHISPAN == 1 ~ "Hispanic",
      RAHISPAN == 0 & RARACEM == 1 ~ "Non-Hispanic White",
      RAHISPAN == 0 & RARACEM == 2 ~ "Non-Hispanic Black",
      RAHISPAN == 0 & RARACEM == 3 ~ "Other",
      TRUE ~ NA_character_
    ),
    
    partnered_2012 = case_when(
      R11MSTAT %in% c(1, 2, 3) ~ 1,
      R11MSTAT %in% c(4, 5, 6, 7, 8) ~ 0,
      TRUE ~ NA_real_
    ),
    
    smoking_now_2012 = R11SMOKEN,
    drinking_any_2012 = R11DRINK,
    bmi_2012 = R11BMI,
    
    self_rated_health_2012 = case_when(
      R11SHLT %in% 1:5 ~ 6 - R11SHLT,
      TRUE ~ NA_real_
    ),
    
    insurance_3cat = case_when(
      R11HIOTHP == 1 ~ "Private",
      R11HIOTHP == 0 & R11HIGOV == 1 ~ "Public only",
      R11HIOTHP == 0 & R11HIGOV == 0 ~ "None",
      TRUE ~ NA_character_
    )
  )

#*******************************统计各个协变量的缺失情况******************************
summary(hrs_ses[c(
  "age_2012",
  "sex_1male_0female",
  "race_ethnicity_4cat",
  "partnered_2012",
  "smoking_now_2012",
  "drinking_any_2012",
  "bmi_2012",
#  "self_rated_health_2012", 暂时不处理，放在敏感性分析中
#  "chronic_condition_count_2012",暂时不处理，放在敏感性分析中
#  "psych_history_2012", 暂时不处理，放在敏感性分析中
  "insurance_3cat"
)])

sapply(hrs_ses[c(
  "age_2012",   #缺失177人
  "sex_1male_0female", #12
  "race_ethnicity_4cat", #24
  "partnered_2012", #181
  "smoking_now_2012", #213
  "drinking_any_2012", #178
  "bmi_2012", #258
#  "self_rated_health_2012",#183, 暂时不处理，放在敏感性分析中
#  "chronic_condition_count_2012",#177, 暂时不处理，放在敏感性分析中
#  "psych_history_2012",#177, 暂时不处理，放在敏感性分析中
  "insurance_3cat"#243
)], function(x) sum(is.na(x)))


#*******************************去掉协变量的缺失值******************************
covariates <- c(
  "age_2012",   #缺失177人
  "sex_1male_0female", #12
  "race_ethnicity_4cat", #24
  "partnered_2012", #181
  "smoking_now_2012", #213
  "drinking_any_2012", #178
  "bmi_2012" #258
  #  "self_rated_health_2012",#183, 暂时不处理，放在敏感性分析中
# "chronic_condition_count_2012",#177, 暂时不处理，放在敏感性分析中
  #  "psych_history_2012",#177, 暂时不处理，放在敏感性分析中
#  "insurance_3cat"#243, 暂时不处理，放在敏感性分析中
)

# 1=协变量完整，0=至少缺失一个协变量
hrs_ses$covariate_complete_2012 <-
  ifelse(complete.cases(hrs_ses[covariates]), 1, 0)

table(hrs_ses$covariate_complete_2012)

#删除协变量缺失者
hrs_ses <- hrs_ses %>%
  filter(covariate_complete_2012 == 1)
#还剩6990人（再减去只有一次心理测量结果的295人）

####绘制基线表####
#指定需要汇总的变量
dput(names(hrs_ses))
myVars <- c(
  # 协变量
  "age_2012", "sex_1male_0female","race_ethnicity_4cat","partnered_2012",
  "smoking_now_2012", "drinking_any_2012",  "bmi_2012", "insurance_3cat",
  "self_rated_health_2012", "chronic_condition_count_2012",
  # SES
  "wealth_tertile_score", "income_tertile_score", "education_3cat", "ses_3cat", "ses_3var_score_0_6",
  # CES-D
  "cesd_score_0_8_2012","cesd_good_z_2012",
  # 食品不安全
  "fi_score_0_6", "fi_USDA6_binary","fi_USDA6_3cat")

#指定分类变量
catVars <- c(
  "sex_1male_0female",
  "race_ethnicity_4cat",
  "partnered_2012",
  "smoking_now_2012",
  "drinking_any_2012",
  "wealth_tertile_score",
  "income_tertile_score",
  "education_3cat",
  "fi_USDA6_binary",
  "fi_USDA6_3cat",
  "insurance_3cat",
  "self_rated_health_2012"
)

### 优化基线特征表
tab2 <- CreateTableOne(
  vars = myVars,
  strata = "ses_3cat",
  data = hrs_ses,
  factorVars = catVars,
  addOverall = TRUE,
  test = FALSE
)

tab2

tab4Mat <- print(
  tab2,
  smd = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  showAllLevels = TRUE,
  contDigits = 1,
  catDigits = 1,
  smdDigits = 3
)

# 转为data.frame
tab4Mat <- as.data.frame(tab4Mat)

# 确认是data.frame
tab4Mat <- as.data.frame(tab4Mat)

# SMD所在列
smd_col <- ncol(tab4Mat)

# SMD保留3位小数
smd_value <- suppressWarnings(
  as.numeric(tab4Mat[[smd_col]])
)

tab4Mat[[smd_col]] <- ifelse(
  is.na(smd_value),
  "",
  sprintf("%.3f", smd_value)
)

# 除SMD外，其他列保留1位小数
other_cols <- setdiff(
  seq_len(ncol(tab4Mat)),
  smd_col
)

tab4Mat[other_cols] <- lapply(
  tab4Mat[other_cols],
  function(x) {
    gsub(
      "([0-9]+\\.[0-9])[0-9]+",
      "\\1",
      x
    )
  }
)
# 再转一次data.frame
tab4Mat <- as.data.frame(tab4Mat)
write.csv(tab4Mat, file = "/dataset/HRS/dietaryindex_standard_input/投稿——frontiers in nutrition/返修/baseline_table_FI.csv")


####无应答者分析####





####Table S2 SES和FI的关联性分析####
#***********************************连续SES评分与二分类FI***********************************#
table(hrs_ses$fi_USDA6_binary, useNA = "always")
table(hrs_ses$ses_3var_score_0_6, useNA = "always")

m_score_1 <- glm(
  fi_USDA6_binary ~
    ses_3var_score_0_6 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)
#计算HC0稳健标准误
robust_vcov_1 <- sandwich::vcovHC(
  m_score_1,
  type = "HC0"
)

robust_test_1 <- lmtest::coeftest(
  m_score_1,
  vcov. = robust_vcov_1
)

robust_test_1

#提取连续SES评分的结果
beta_1 <- robust_test_1["ses_3var_score_0_6", 1]
se_1 <- robust_test_1["ses_3var_score_0_6", 2]
p_1 <- robust_test_1["ses_3var_score_0_6", 4]

PR_1 <- exp(beta_1)
lower_1 <- exp(beta_1 - 1.96 * se_1)
upper_1 <- exp(beta_1 + 1.96 * se_1)

#整理结果
result_score_binary <- data.frame(
  SES_measure = "SES score, per 1-point increase",
  PR = round(PR_1, 3),
  lower_95CI = round(lower_1, 3),
  upper_95CI = round(upper_1, 3),
  p_value = ifelse(
    p_1 < 0.001,
    "<0.001",
    sprintf("%.3f", p_1)
  )
)

result_score_binary

#***********************************分类SES等级与二分类FI***********************************#
table(hrs_ses$ses_3cat, useNA = "always")
#设置高SES为参照组
hrs_ses$ses_3cat <- factor(
  hrs_ses$ses_3cat,
  levels = c("High SES", "Middle SES", "Low SES")
)

#检查参照顺序
levels(hrs_ses$ses_3cat)

m_category_1 <- glm(
  fi_USDA6_binary ~
    ses_3cat +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

#计算稳健标准误
robust_vcov_2 <- sandwich::vcovHC(
  m_category_1,
  type = "HC0"
)

robust_test_2 <- lmtest::coeftest(
  m_category_1,
  vcov. = robust_vcov_2
)

robust_test_2

#提取中SES结果
beta_middle <- robust_test_2["ses_3catMiddle SES", 1]
se_middle <- robust_test_2["ses_3catMiddle SES", 2]
p_middle <- robust_test_2["ses_3catMiddle SES", 4]

PR_middle <- exp(beta_middle)
lower_middle <- exp(beta_middle - 1.96 * se_middle)
upper_middle <- exp(beta_middle + 1.96 * se_middle)

#提取低SES结果
beta_low <- robust_test_2["ses_3catLow SES", 1]
se_low <- robust_test_2["ses_3catLow SES", 2]
p_low <- robust_test_2["ses_3catLow SES", 4]

PR_low <- exp(beta_low)
lower_low <- exp(beta_low - 1.96 * se_low)
upper_low <- exp(beta_low + 1.96 * se_low)

#整理结果
result_category_binary <- data.frame(
  SES_category = c("High SES", "Middle SES", "Low SES"),
  
  PR = c(
    1,
    PR_middle,
    PR_low
  ),
  
  lower_95CI = c(
    NA,
    lower_middle,
    lower_low
  ),
  
  upper_95CI = c(
    NA,
    upper_middle,
    upper_low
  ),
  
  p_value = c(
    "Reference",
    ifelse(p_middle < 0.001, "<0.001", sprintf("%.3f", p_middle)),
    ifelse(p_low < 0.001, "<0.001", sprintf("%.3f", p_low))
  )
)

result_category_binary$PR <- round(
  result_category_binary$PR,
  3
)

result_category_binary$lower_95CI <- round(
  result_category_binary$lower_95CI,
  3
)

result_category_binary$upper_95CI <- round(
  result_category_binary$upper_95CI,
  3
)

result_category_binary

#***********************************连续SES评分与5个食品不安全单项指标***********************************#
#先生成5个二分类结局
hrs_ses$fi_food_not_last <- ifelse(
  hrs_ses$fi_b1 %in% c(1, 2), 1,
  ifelse(hrs_ses$fi_b1 == 3, 0, NA)
)

hrs_ses$fi_cannot_afford_balanced <- ifelse(
  hrs_ses$fi_b2 %in% c(1, 2), 1,
  ifelse(hrs_ses$fi_b2 == 3, 0, NA)
)

hrs_ses$fi_cut_skip_meals <- ifelse(
  hrs_ses$fi_b3 %in% c(1, 2, 3), 1,
  ifelse(hrs_ses$fi_b3 == 5, 0, NA)
)

hrs_ses$fi_ate_less <- ifelse(
  hrs_ses$fi_b4 == 1, 1,
  ifelse(hrs_ses$fi_b4 == 5, 0, NA)
)

hrs_ses$fi_hungry <- ifelse(
  hrs_ses$fi_b5 == 1, 1,
  ifelse(hrs_ses$fi_b5 == 5, 0, NA)
)

#检查频数
table(hrs_ses$fi_cannot_afford_balanced)
table(hrs_ses$fi_cut_skip_meals)
table(hrs_ses$fi_ate_less)
table(hrs_ses$fi_hungry)
table(hrs_ses$fi_food_not_last)

#建立回归模型
m_score_balanced <- glm(
  fi_cannot_afford_balanced ~
    ses_3var_score_0_6 + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_score_cut_skip <- glm(
  fi_cut_skip_meals ~
    ses_3var_score_0_6 + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_score_ate_less <- glm(
  fi_ate_less ~
    ses_3var_score_0_6 + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_score_hungry <- glm(
  fi_hungry ~
    ses_3var_score_0_6 + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_score_food_last <- glm(
  fi_food_not_last ~
    ses_3var_score_0_6 + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

#计算稳健检验结果
test_balanced <- lmtest::coeftest(
  m_score_balanced,
  vcov. = sandwich::vcovHC(m_score_balanced, type = "HC0")
)

test_cut_skip <- lmtest::coeftest(
  m_score_cut_skip,
  vcov. = sandwich::vcovHC(m_score_cut_skip, type = "HC0")
)

test_ate_less <- lmtest::coeftest(
  m_score_ate_less,
  vcov. = sandwich::vcovHC(m_score_ate_less, type = "HC0")
)

test_hungry <- lmtest::coeftest(
  m_score_hungry,
  vcov. = sandwich::vcovHC(m_score_hungry, type = "HC0")
)

test_food_last <- lmtest::coeftest(
  m_score_food_last,
  vcov. = sandwich::vcovHC(m_score_food_last, type = "HC0")
)

#汇总连续SES系数
result_score_items <- data.frame(
  outcome = c(
    "Could not afford balanced meals",
    "Cut or skipped meals",
    "Ate less due to lack of money",
    "Hungry due to lack of money",
    "Food did not last"
  ),
  
  beta = c(
    test_balanced["ses_3var_score_0_6", 1],
    test_cut_skip["ses_3var_score_0_6", 1],
    test_ate_less["ses_3var_score_0_6", 1],
    test_hungry["ses_3var_score_0_6", 1],
    test_food_last["ses_3var_score_0_6", 1]
  ),
  
  robust_se = c(
    test_balanced["ses_3var_score_0_6", 2],
    test_cut_skip["ses_3var_score_0_6", 2],
    test_ate_less["ses_3var_score_0_6", 2],
    test_hungry["ses_3var_score_0_6", 2],
    test_food_last["ses_3var_score_0_6", 2]
  ),
  
  p_value = c(
    test_balanced["ses_3var_score_0_6", 4],
    test_cut_skip["ses_3var_score_0_6", 4],
    test_ate_less["ses_3var_score_0_6", 4],
    test_hungry["ses_3var_score_0_6", 4],
    test_food_last["ses_3var_score_0_6", 4]
  )
)

#计算PR及95%C
result_score_items$PR <- exp(result_score_items$beta)

result_score_items$lower_95CI <- exp(
  result_score_items$beta -
    1.96 * result_score_items$robust_se
)

result_score_items$upper_95CI <- exp(
  result_score_items$beta +
    1.96 * result_score_items$robust_se
)

result_score_items$PR_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  result_score_items$PR,
  result_score_items$lower_95CI,
  result_score_items$upper_95CI
)

result_score_items$p_report <- ifelse(
  result_score_items$p_value < 0.001,
  "<0.001",
  sprintf("%.3f", result_score_items$p_value)
)

result_score_items[
  c("outcome", "PR_95CI", "p_report")
]


#***********************************三分类SES与5个食品不安全单项指标，高SES为参照组***********************************#
#确认参照组
hrs_ses$ses_3cat <- factor(
  hrs_ses$ses_3cat,
  levels = c("High SES", "Middle SES", "Low SES")
)

levels(hrs_ses$ses_3cat)

#建立回归模型
m_cat_balanced <- glm(
  fi_cannot_afford_balanced ~
    ses_3cat + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_cat_cut_skip <- glm(
  fi_cut_skip_meals ~
    ses_3cat + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_cat_ate_less <- glm(
  fi_ate_less ~
    ses_3cat + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_cat_hungry <- glm(
  fi_hungry ~
    ses_3cat + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

m_cat_food_last <- glm(
  fi_food_not_last ~
    ses_3cat + age_2012 + sex_1male_0female +
    race_ethnicity_4cat + partnered_2012 +
    smoking_now_2012 + drinking_any_2012 + bmi_2012,
  family = poisson(link = "log"),
  data = hrs_ses
)

#计算稳健检验结果
test_cat_balanced <- lmtest::coeftest(
  m_cat_balanced,
  vcov. = sandwich::vcovHC(m_cat_balanced, type = "HC0")
)

test_cat_cut_skip <- lmtest::coeftest(
  m_cat_cut_skip,
  vcov. = sandwich::vcovHC(m_cat_cut_skip, type = "HC0")
)

test_cat_ate_less <- lmtest::coeftest(
  m_cat_ate_less,
  vcov. = sandwich::vcovHC(m_cat_ate_less, type = "HC0")
)

test_cat_hungry <- lmtest::coeftest(
  m_cat_hungry,
  vcov. = sandwich::vcovHC(m_cat_hungry, type = "HC0")
)

test_cat_food_last <- lmtest::coeftest(
  m_cat_food_last,
  vcov. = sandwich::vcovHC(m_cat_food_last, type = "HC0")
)

#整理中SES和低SES结果
ses_terms <- c(
  "ses_3catMiddle SES",
  "ses_3catLow SES"
)

result_cat_items <- data.frame(
  outcome = rep(
    c(
      "Could not afford balanced meals",
      "Cut or skipped meals",
      "Ate less due to lack of money",
      "Hungry due to lack of money",
      "Food did not last"
    ),
    each = 2
  ),
  
  SES_category = rep(
    c("Middle SES", "Low SES"),
    times = 5
  ),
  
  beta = c(
    test_cat_balanced[ses_terms, 1],
    test_cat_cut_skip[ses_terms, 1],
    test_cat_ate_less[ses_terms, 1],
    test_cat_hungry[ses_terms, 1],
    test_cat_food_last[ses_terms, 1]
  ),
  
  robust_se = c(
    test_cat_balanced[ses_terms, 2],
    test_cat_cut_skip[ses_terms, 2],
    test_cat_ate_less[ses_terms, 2],
    test_cat_hungry[ses_terms, 2],
    test_cat_food_last[ses_terms, 2]
  ),
  
  p_value = c(
    test_cat_balanced[ses_terms, 4],
    test_cat_cut_skip[ses_terms, 4],
    test_cat_ate_less[ses_terms, 4],
    test_cat_hungry[ses_terms, 4],
    test_cat_food_last[ses_terms, 4]
  )
)

#计算PR和95%CI
result_cat_items$PR <- exp(result_cat_items$beta)

result_cat_items$lower_95CI <- exp(
  result_cat_items$beta -
    1.96 * result_cat_items$robust_se
)

result_cat_items$upper_95CI <- exp(
  result_cat_items$beta +
    1.96 * result_cat_items$robust_se
)

result_cat_items$PR_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  result_cat_items$PR,
  result_cat_items$lower_95CI,
  result_cat_items$upper_95CI
)

result_cat_items$p_report <- ifelse(
  result_cat_items$p_value < 0.001,
  "<0.001",
  sprintf("%.3f", result_cat_items$p_value)
)

result_cat_items[
  c(
    "outcome",
    "SES_category",
    "PR_95CI",
    "p_report"
  )
]


#***********************************敏感性分析：连续SES评分与连续FI***********************************#
m_FIscore_SES_continuous <- lm(
  fi_score_0_6 ~
    ses_3var_score_0_6 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012,
  data = hrs_ses
)

vcov_FIscore_continuous <- sandwich::vcovHC(
  m_FIscore_SES_continuous,
  type = "HC0"
)

test_FIscore_continuous <- lmtest::coeftest(
  m_FIscore_SES_continuous,
  vcov. = vcov_FIscore_continuous
)

beta <- test_FIscore_continuous[
  "ses_3var_score_0_6",
  "Estimate"
]

se <- test_FIscore_continuous[
  "ses_3var_score_0_6",
  "Std. Error"
]

p <- test_FIscore_continuous[
  "ses_3var_score_0_6",
  "Pr(>|t|)"
]

continuous_result <- data.frame(
  SES_measure = "SES score, per 1-point increase",
  Beta_95CI = sprintf(
    "%.3f (%.3f, %.3f)",
    beta,
    beta - 1.96 * se,
    beta + 1.96 * se
  ),
  P = ifelse(
    p < 0.001,
    "<0.001",
    sprintf("%.3f", p)
  )
)

continuous_result


#***********************************敏感性分析：三分类SES与连续FI***********************************#
m_FIscore_SES_category <- lm(
  fi_score_0_6 ~
    ses_3cat +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012,
  data = hrs_ses
)

vcov_FIscore_category <- sandwich::vcovHC(
  m_FIscore_SES_category,
  type = "HC0"
)

test_FIscore_category <- lmtest::coeftest(
  m_FIscore_SES_category,
  vcov. = vcov_FIscore_category
)

test_FIscore_category

# Middle SES
beta_middle <- test_FIscore_category[
  "ses_3catMiddle SES",
  "Estimate"
]

se_middle <- test_FIscore_category[
  "ses_3catMiddle SES",
  "Std. Error"
]

p_middle <- test_FIscore_category[
  "ses_3catMiddle SES",
  "Pr(>|t|)"
]

# Low SES
beta_low <- test_FIscore_category[
  "ses_3catLow SES",
  "Estimate"
]

se_low <- test_FIscore_category[
  "ses_3catLow SES",
  "Std. Error"
]

p_low <- test_FIscore_category[
  "ses_3catLow SES",
  "Pr(>|t|)"
]


categorical_result <- data.frame(
  SES_measure = c(
    "High SES",
    "Middle SES",
    "Low SES"
  ),
  Beta_95CI = c(
    "Reference",
    sprintf(
      "%.3f (%.3f, %.3f)",
      beta_middle,
      beta_middle - 1.96 * se_middle,
      beta_middle + 1.96 * se_middle
    ),
    sprintf(
      "%.3f (%.3f, %.3f)",
      beta_low,
      beta_low - 1.96 * se_low,
      beta_low + 1.96 * se_low
    )
  ),
  P = c(
    "-",
    ifelse(p_middle < 0.001, "<0.001", sprintf("%.3f", p_middle)),
    ifelse(p_low < 0.001, "<0.001", sprintf("%.3f", p_low))
  )
)

categorical_result


#***********************************敏感性分析：SES三个组成成分与二分类FI***********************************#
# Model 1：财富
model1_component_FI <- glm(
  fi_USDA6_binary ~
    wealth_component +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat,
  family = poisson(link = "log"),
  data = hrs_ses
)

vcov1_component_FI <- sandwich::vcovHC(
  model1_component_FI,
  type = "HC0"
)

test1_component_FI <- lmtest::coeftest(
  model1_component_FI,
  vcov. = vcov1_component_FI
)

# Model 2：收入
model2_component_FI <- glm(
  fi_USDA6_binary ~
    income_component +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat,
  family = poisson(link = "log"),
  data = hrs_ses
)

vcov2_component_FI <- sandwich::vcovHC(
  model2_component_FI,
  type = "HC0"
)

test2_component_FI <- lmtest::coeftest(
  model2_component_FI,
  vcov. = vcov2_component_FI
)

# Model 3：教育
model3_component_FI <- glm(
  fi_USDA6_binary ~
    education_component +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat,
  family = poisson(link = "log"),
  data = hrs_ses
)

vcov3_component_FI <- sandwich::vcovHC(
  model3_component_FI,
  type = "HC0"
)

test3_component_FI <- lmtest::coeftest(
  model3_component_FI,
  vcov. = vcov3_component_FI
)

# Model 4：三个SES成分共同校正
model4_component_FI <- glm(
  fi_USDA6_binary ~
    wealth_component +
    income_component +
    education_component +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat,
  family = poisson(link = "log"),
  data = hrs_ses
)

vcov4_component_FI <- sandwich::vcovHC(
  model4_component_FI,
  type = "HC0"
)

test4_component_FI <- lmtest::coeftest(
  model4_component_FI,
  vcov. = vcov4_component_FI
)


# 分别提取三个单独模型
wealth_single <- test1_component_FI[
  grep("^wealth_component", rownames(test1_component_FI)),
  ,
  drop = FALSE
]

income_single <- test2_component_FI[
  grep("^income_component", rownames(test2_component_FI)),
  ,
  drop = FALSE
]

education_single <- test3_component_FI[
  grep("^education_component", rownames(test3_component_FI)),
  ,
  drop = FALSE
]

# 从共同校正模型中提取三个SES成分
wealth_joint <- test4_component_FI[
  grep("^wealth_component", rownames(test4_component_FI)),
  ,
  drop = FALSE
]

income_joint <- test4_component_FI[
  grep("^income_component", rownames(test4_component_FI)),
  ,
  drop = FALSE
]

education_joint <- test4_component_FI[
  grep("^education_component", rownames(test4_component_FI)),
  ,
  drop = FALSE
]

#生成最终表格
single_beta <- c(
  wealth_single[, 1],
  income_single[, 1],
  education_single[, 1]
)

single_se <- c(
  wealth_single[, 2],
  income_single[, 2],
  education_single[, 2]
)

single_p <- c(
  wealth_single[, 4],
  income_single[, 4],
  education_single[, 4]
)

joint_beta <- c(
  wealth_joint[, 1],
  income_joint[, 1],
  education_joint[, 1]
)

joint_se <- c(
  wealth_joint[, 2],
  income_joint[, 2],
  education_joint[, 2]
)

joint_p <- c(
  wealth_joint[, 4],
  income_joint[, 4],
  education_joint[, 4]
)

component_FI_table <- data.frame(
  SES_component = rep(
    c("Wealth", "Income", "Education"),
    each = 3
  ),
  Level = rep(
    c("High", "Middle", "Low"),
    times = 3
  ),
  Separate_model_aPR = c(
    "Reference",
    sprintf("%.3f (%.3f, %.3f)",
            exp(single_beta[1]),
            exp(single_beta[1] - 1.96 * single_se[1]),
            exp(single_beta[1] + 1.96 * single_se[1])),
    sprintf("%.3f (%.3f, %.3f)",
            exp(single_beta[2]),
            exp(single_beta[2] - 1.96 * single_se[2]),
            exp(single_beta[2] + 1.96 * single_se[2])),
    "Reference",
    sprintf("%.3f (%.3f, %.3f)",
            exp(single_beta[3]),
            exp(single_beta[3] - 1.96 * single_se[3]),
            exp(single_beta[3] + 1.96 * single_se[3])),
    sprintf("%.3f (%.3f, %.3f)",
            exp(single_beta[4]),
            exp(single_beta[4] - 1.96 * single_se[4]),
            exp(single_beta[4] + 1.96 * single_se[4])),
    "Reference",
    sprintf("%.3f (%.3f, %.3f)",
            exp(single_beta[5]),
            exp(single_beta[5] - 1.96 * single_se[5]),
            exp(single_beta[5] + 1.96 * single_se[5])),
    sprintf("%.3f (%.3f, %.3f)",
            exp(single_beta[6]),
            exp(single_beta[6] - 1.96 * single_se[6]),
            exp(single_beta[6] + 1.96 * single_se[6]))
  )
)

#将共同校正的结果合并进来
component_FI_table$Joint_model_aPR <- c(
  "Reference",
  sprintf(
    "%.3f (%.3f, %.3f)",
    exp(joint_beta[1]),
    exp(joint_beta[1] - 1.96 * joint_se[1]),
    exp(joint_beta[1] + 1.96 * joint_se[1])
  ),
  sprintf(
    "%.3f (%.3f, %.3f)",
    exp(joint_beta[2]),
    exp(joint_beta[2] - 1.96 * joint_se[2]),
    exp(joint_beta[2] + 1.96 * joint_se[2])
  ),
  
  "Reference",
  sprintf(
    "%.3f (%.3f, %.3f)",
    exp(joint_beta[3]),
    exp(joint_beta[3] - 1.96 * joint_se[3]),
    exp(joint_beta[3] + 1.96 * joint_se[3])
  ),
  sprintf(
    "%.3f (%.3f, %.3f)",
    exp(joint_beta[4]),
    exp(joint_beta[4] - 1.96 * joint_se[4]),
    exp(joint_beta[4] + 1.96 * joint_se[4])
  ),
  
  "Reference",
  sprintf(
    "%.3f (%.3f, %.3f)",
    exp(joint_beta[5]),
    exp(joint_beta[5] - 1.96 * joint_se[5]),
    exp(joint_beta[5] + 1.96 * joint_se[5])
  ),
  sprintf(
    "%.3f (%.3f, %.3f)",
    exp(joint_beta[6]),
    exp(joint_beta[6] - 1.96 * joint_se[6]),
    exp(joint_beta[6] + 1.96 * joint_se[6])
  )
)

#加入模型P值
component_FI_table$Separate_model_P <- c(
  "—",
  ifelse(single_p[1] < 0.001, "<0.001", sprintf("%.3f", single_p[1])),
  ifelse(single_p[2] < 0.001, "<0.001", sprintf("%.3f", single_p[2])),
  
  "—",
  ifelse(single_p[3] < 0.001, "<0.001", sprintf("%.3f", single_p[3])),
  ifelse(single_p[4] < 0.001, "<0.001", sprintf("%.3f", single_p[4])),
  
  "—",
  ifelse(single_p[5] < 0.001, "<0.001", sprintf("%.3f", single_p[5])),
  ifelse(single_p[6] < 0.001, "<0.001", sprintf("%.3f", single_p[6]))
)

component_FI_table$Joint_model_P <- c(
  "—",
  ifelse(joint_p[1] < 0.001, "<0.001", sprintf("%.3f", joint_p[1])),
  ifelse(joint_p[2] < 0.001, "<0.001", sprintf("%.3f", joint_p[2])),
  
  "—",
  ifelse(joint_p[3] < 0.001, "<0.001", sprintf("%.3f", joint_p[3])),
  ifelse(joint_p[4] < 0.001, "<0.001", sprintf("%.3f", joint_p[4])),
  
  "—",
  ifelse(joint_p[5] < 0.001, "<0.001", sprintf("%.3f", joint_p[5])),
  ifelse(joint_p[6] < 0.001, "<0.001", sprintf("%.3f", joint_p[6]))
)

component_FI_table <- component_FI_table[
  ,
  c(
    "SES_component",
    "Level",
    "Separate_model_aPR",
    "Separate_model_P",
    "Joint_model_aPR",
    "Joint_model_P"
  )
]

component_FI_table


####生成长格式数据，做LMM####
#由于CESD2012属于饮食问卷测量前的，这里不当作随访而当作协变量控制
library(dplyr)
library(tidyr)
library(lme4)

cesd_long <- hrs_ses %>%
  select(
    HHIDPN,
    
    # SES
    wealth_tertile_score,
    income_tertile_score,
    education_3cat,
    ses_3var_score_0_6,
    ses_3cat,
    
    # 食品不安全
    fi_score_0_6,
    fi_USDA6_binary,
    fi_USDA6_3cat,
    fi_food_not_last,
    fi_cannot_afford_balanced,
    fi_cut_skip_meals,
    fi_ate_less,
    fi_hungry,
    
    # 2012年基线CES-D
    cesd_score_0_8_2012,
    cesd_good_z_2012,
    
    # 2014—2020年随访CES-D
    cesd_good_z_2014,
    cesd_good_z_2016,
    cesd_good_z_2018,
    cesd_good_z_2020,
    
    # 基线协变量
    age_2012,
    sex_1male_0female,
    race_ethnicity_4cat,
    partnered_2012,
    smoking_now_2012,
    drinking_any_2012,
    bmi_2012,
    self_rated_health_2012,
    chronic_condition_count_2012,
    insurance_3cat
  ) %>%
  
  pivot_longer(
    cols = c(
      cesd_good_z_2014,
      cesd_good_z_2016,
      cesd_good_z_2018,
      cesd_good_z_2020
    ),
    names_to = "cesd_wave",
    values_to = "cesd_good_z"
  ) %>%
  
  mutate(
    time = case_when(
      cesd_wave == "cesd_good_z_2014" ~ 0,
      cesd_wave == "cesd_good_z_2016" ~ 1,
      cesd_wave == "cesd_good_z_2018" ~ 2,
      cesd_wave == "cesd_good_z_2020" ~ 3
    ),
    
    year = case_when(
      time == 0 ~ 2014,
      time == 1 ~ 2016,
      time == 2 ~ 2018,
      time == 3 ~ 2020
    ),
    
    HHIDPN = factor(HHIDPN),
    cesd_wave = factor(
      cesd_wave,
      levels = c(
        "cesd_good_z_2014",
        "cesd_good_z_2016",
        "cesd_good_z_2018",
        "cesd_good_z_2020"
      )
    )
  )

#检查长数据
nrow(cesd_long)

dplyr::n_distinct(cesd_long$HHIDPN)

table(cesd_long$year)

table(
  table(cesd_long$HHIDPN)
)

#建立空模型
null_model <- lmer(
  cesd_good_z ~ 1 + (1 | HHIDPN),
  data = cesd_long,
  REML = TRUE
)

summary(null_model)
VarCorr(null_model)

#直接计算ICC
person_variance <- as.numeric(
  VarCorr(null_model)$HHIDPN[1, 1]
)

residual_variance <- sigma(null_model)^2

ICC <- person_variance /
  (person_variance + residual_variance)

ICC


####图2：按照SES分组的CESD轨迹####
library(ggplot2)
library(patchwork)
trajectory_long <- hrs_ses %>%
  select(
    HHIDPN,
    ses_3cat,
    fi_USDA6_binary,
    cesd_good_z_2012,
    cesd_good_z_2014,
    cesd_good_z_2016,
    cesd_good_z_2018,
    cesd_good_z_2020
  ) %>%
  pivot_longer(
    cols = starts_with("cesd_good_z_"),
    names_to = "year",
    values_to = "cesd_good_z"
  ) %>%
  mutate(
    year = as.numeric(
      sub("cesd_good_z_", "", year)
    ),
    
    ses_3cat = factor(
      ses_3cat,
      levels = c(
        "Low SES",
        "Middle SES",
        "High SES"
      )
    ),
    
    fi_group = ifelse(
      as.character(fi_USDA6_binary) == "1",
      "Food insecure",
      "Food secure"
    ),
    
    fi_group = factor(
      fi_group,
      levels = c(
        "Food secure",
        "Food insecure"
      )
    )
  )

#计算均值与置信区间
ses_trajectory <- trajectory_long %>%
  group_by(ses_3cat, year) %>%
  summarise(
    n = sum(!is.na(cesd_good_z)),
    mean = mean(cesd_good_z, na.rm = TRUE),
    se = sd(cesd_good_z, na.rm = TRUE) / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  )

fi_trajectory <- trajectory_long %>%
  group_by(fi_group, year) %>%
  summarise(
    n = sum(!is.na(cesd_good_z)),
    mean = mean(cesd_good_z, na.rm = TRUE),
    se = sd(cesd_good_z, na.rm = TRUE) / sqrt(n),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  )

#panel A
plot_A <- ggplot(
  ses_trajectory,
  aes(
    x = year,
    y = mean,
    color = ses_3cat,
    fill = ses_3cat,
    group = ses_3cat
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.18,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c(
      "Low SES" = "#294D5B",
      "Middle SES" = "#78949D",
      "High SES" = "#BCCDD2"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Low SES" = "#294D5B",
      "Middle SES" = "#78949D",
      "High SES" = "#BCCDD2"
    )
  ) +
  scale_x_continuous(
    breaks = c(2012, 2014, 2016, 2018, 2020)
  ) +
  scale_y_continuous(
    breaks = seq(-0.75, 0.50, 0.25)
  ) +
  coord_cartesian(ylim = c(-0.75, 0.50)) +
  labs(
    x = NULL,
    y = "Standardized CES-D score",
    color = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.width = unit(0.7, "cm"),
    legend.spacing.x = unit(0.10, "cm"),
    panel.grid.major.y = element_line(
      color = "#E8E8E8",
      linewidth = 0.35
    ),
    axis.text = element_text(color = "#404040"),
    axis.title.y = element_text(size = 9),
    plot.margin = margin(4, 8, 4, 4)
  )


plot_B <- ggplot(
  fi_trajectory,
  aes(
    x = year,
    y = mean,
    color = fi_group,
    fill = fi_group,
    group = fi_group
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.18,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c(
      "Food secure" = "#294D5B",
      "Food insecure" = "#B65C5C"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Food secure" = "#294D5B",
      "Food insecure" = "#B65C5C"
    )
  ) +
  scale_x_continuous(
    breaks = c(2012, 2014, 2016, 2018, 2020)
  ) +
  scale_y_continuous(
    breaks = seq(-0.75, 0.50, 0.25)
  ) +
  coord_cartesian(ylim = c(-0.75, 0.50)) +
  labs(
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.width = unit(0.7, "cm"),
    legend.spacing.x = unit(0.10, "cm"),
    panel.grid.major.y = element_line(
      color = "#E8E8E8",
      linewidth = 0.35
    ),
    axis.text = element_text(color = "#404040"),
    plot.margin = margin(4, 4, 4, 8)
  )

figure_trajectory <- plot_A + plot_B +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      size = 12,
      face = "plain"
    ),
    plot.tag.position = c(0, 1)
  )

figure_trajectory
####SES和CESD####
#***********************************Model 0：基础模型***********************************#
model_0 <- lmer(
  cesd_good_z ~
    time +
    cesd_good_z_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

summary(model_0)

#检查模型
isSingular(
  model_0,
  tol = 1e-4
)

VarCorr(model_0)

#提取固定效应及95%CI
model_0_result <- coef(
  summary(model_0)
)

model_0_CI <- confint(
  model_0,
  parm = "beta_",
  method = "Wald"
)

model_0_result
model_0_CI

#***********************************Model 1：SES基础关联***********************************#
#Model 1检验：在控制随访时间和2012年基线心理健康后，连续SES评分是否与后续心理健康水平相关。
model1_SES_CESD <- lmer(
  cesd_good_z ~
    ses_3var_score_0_6 +
    time +
    cesd_good_z_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

#查看模型结果
summary(model1_SES_CESD)

# Model 1结果格式化
model1_result <- tidy(
  model1_SES_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model1_result <- model1_result[
  model1_result$term %in% c(
    "ses_3var_score_0_6",
    "time"
  ),
]

model1_result$Variable <- c(
  "SES score, per 1-point increase",
  "Time, per 2-year interval"
)

model1_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model1_result$estimate,
  model1_result$conf.low,
  model1_result$conf.high
)

model1_result$P <- ifelse(
  model1_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model1_result$p.value)
)

model1_result <- model1_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model1_result

#***********************************Model 2：主要混杂因素调整***********************************#
#Model 2在Model 1基础上加入核心混杂因素，然后直接格式化SES和time的结果。
# Model 2：加入核心混杂因素
model2_SES_CESD <- lmer(
  cesd_good_z ~
    ses_3var_score_0_6 +
    time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Model 2结果
model2_result <- tidy(
  model2_SES_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model2_result <- model2_result[
  model2_result$term %in% c(
    "ses_3var_score_0_6",
    "time"
  ),
]

model2_result$Variable <- c(
  "SES score, per 1-point increase",
  "Time, per 2-year interval"
)

model2_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model2_result$estimate,
  model2_result$conf.low,
  model2_result$conf.high
)

model2_result$P <- ifelse(
  model2_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model2_result$p.value)
)

model2_result <- model2_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model2_result


#***********************************Model 3：SES轨迹差异***********************************#
# Model 3：在Model 2基础上加入SES与时间的交互项
model3_SES_CESD <- lmer(
  cesd_good_z ~
    ses_3var_score_0_6 * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)
# 格式化Model 3结果
model3_result <- tidy(
  model3_SES_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model3_result <- model3_result[
  model3_result$term %in% c(
    "ses_3var_score_0_6",
    "time",
    "ses_3var_score_0_6:time"
  ),
]

# 固定行顺序
model3_result <- model3_result[
  match(
    c(
      "ses_3var_score_0_6",
      "time",
      "ses_3var_score_0_6:time"
    ),
    model3_result$term
  ),
]

model3_result$Variable <- c(
  "SES score, per 1-point increase",
  "Time, per 2-year interval",
  "SES score × time"
)

model3_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model3_result$estimate,
  model3_result$conf.low,
  model3_result$conf.high
)

model3_result$P <- ifelse(
  model3_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model3_result$p.value)
)

model3_result <- model3_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model3_result


#***********************************Model 4：在Model 2基础上加入扩展调整变量***********************************#
# Model 4：在Model 2基础上加入扩展调整变量
model4_SES_CESD <- lmer(
  cesd_good_z ~
    ses_3var_score_0_6 +
    time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Model 4结果
model4_result <- tidy(
  model4_SES_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model4_result <- model4_result[
  model4_result$term %in% c(
    "ses_3var_score_0_6",
    "time"
  ),
]

model4_result <- model4_result[
  match(
    c(
      "ses_3var_score_0_6",
      "time"
    ),
    model4_result$term
  ),
]

model4_result$Variable <- c(
  "SES score, per 1-point increase",
  "Time, per 2-year interval"
)

model4_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model4_result$estimate,
  model4_result$conf.low,
  model4_result$conf.high
)

model4_result$P <- ifelse(
  model4_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model4_result$p.value)
)

model4_result <- model4_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model4_result

#**********************************Model 5：在Model 4基础上加入SES与时间的交互项***********************************#
# Model 5：在Model 4基础上加入SES与时间的交互项
model5_SES_CESD <- lmer(
  cesd_good_z ~
    ses_3var_score_0_6 * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Model 5结果
model5_result <- tidy(
  model5_SES_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model5_result <- model5_result[
  model5_result$term %in% c(
    "ses_3var_score_0_6",
    "time",
    "ses_3var_score_0_6:time"
  ),
]

model5_result <- model5_result[
  match(
    c(
      "ses_3var_score_0_6",
      "time",
      "ses_3var_score_0_6:time"
    ),
    model5_result$term
  ),
]

model5_result$Variable <- c(
  "SES score, per 1-point increase",
  "Time, per 2-year interval",
  "SES score × time"
)

model5_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model5_result$estimate,
  model5_result$conf.low,
  model5_result$conf.high
)

model5_result$P <- ifelse(
  model5_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model5_result$p.value)
)

model5_result <- model5_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model5_result


#***********************************补充分析：单个SES和CESD的纵向分析***********************************#
id_match <- match(
  as.character(cesd_long$HHIDPN),
  as.character(hrs_ses$HHIDPN)
)

sum(is.na(id_match))

cesd_long$wealth_component <-
  hrs_ses$wealth_component[id_match]

cesd_long$income_component <-
  hrs_ses$income_component[id_match]

cesd_long$education_component <-
  hrs_ses$education_component[id_match]

table(cesd_long$wealth_component, useNA = "always")
table(cesd_long$income_component, useNA = "always")
table(cesd_long$education_component, useNA = "always")

levels(cesd_long$wealth_component)
levels(cesd_long$income_component)
levels(cesd_long$education_component)

# Table S7, Model 1：财富
model1_component_CESD <- lmer(
  cesd_good_z ~
    wealth_component * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

model1_component_result <- tidy(
  model1_component_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model1_component_result <- model1_component_result[
  model1_component_result$term %in% c(
    "wealth_componentMiddle",
    "wealth_componentLow",
    "time",
    "wealth_componentMiddltime",
    "wealth_componentLow:time"
  ),
]

model1_component_result <- model1_component_result[
  match(
    c(
      "wealth_componentMiddle",
      "wealth_componentLow",
      "time",
      "wealth_componentMiddltime",
      "wealth_componentLow:time"
    ),
    model1_component_result$term
  ),
]

model1_component_result$Variable <- c(
  "Middle wealth",
  "Low wealth",
  "Time, per 2-year interval",
  "Middle wealth × time",
  "Low wealth × time"
)

model1_component_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model1_component_result$estimate,
  model1_component_result$conf.low,
  model1_component_result$conf.high
)

model1_component_result$P <- ifelse(
  model1_component_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model1_component_result$p.value)
)

model1_component_result <- model1_component_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model1_component_result

# Table S7, Model 2：收入
model2_component_CESD <- lmer(
  cesd_good_z ~
    income_component * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Model 2结果
model2_component_result <- tidy(
  model2_component_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model2_component_result <- model2_component_result[
  model2_component_result$term %in% c(
    "income_componentMiddle",
    "income_componentLow",
    "time",
    "income_componentMiddltime",
    "income_componentLow:time"
  ),
]

model2_component_result <- model2_component_result[
  match(
    c(
      "income_componentMiddle",
      "income_componentLow",
      "time",
      "income_componentMiddltime",
      "income_componentLow:time"
    ),
    model2_component_result$term
  ),
]

model2_component_result$Variable <- c(
  "Middle income",
  "Low income",
  "Time, per 2-year interval",
  "Middle income × time",
  "Low income × time"
)

model2_component_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model2_component_result$estimate,
  model2_component_result$conf.low,
  model2_component_result$conf.high
)

model2_component_result$P <- ifelse(
  model2_component_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model2_component_result$p.value)
)

model2_component_result <- model2_component_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model2_component_result


# Table S7, Model 3：教育
model3_component_CESD <- lmer(
  cesd_good_z ~
    education_component * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 提取并格式化Model 3结果
model3_component_result <- tidy(
  model3_component_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model3_component_result <- model3_component_result[
  model3_component_result$term %in% c(
    "education_componentMiddle",
    "education_componentLow",
    "time",
    "education_componentMiddltime",
    "education_componentLow:time"
  ),
]


model3_component_result <- model3_component_result[
  match(
    c(
      "education_componentMiddle",
      "education_componentLow",
      "time",
      "education_componentMiddltime",
      "education_componentLow:time"
    ),
    model3_component_result$term
  ),
]

model3_component_result$Variable <- c(
  "Middle education",
  "Low education",
  "Time, per 2-year interval",
  "Middle education × time",
  "Low education × time"
)

model3_component_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model3_component_result$estimate,
  model3_component_result$conf.low,
  model3_component_result$conf.high
)

model3_component_result$P <- ifelse(
  model3_component_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model3_component_result$p.value)
)

model3_component_result <- model3_component_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model3_component_result



#Table S7 三者共同校正
# 明确分类顺序，最高组作为参考组
cesd_long$wealth_component <- factor(
  cesd_long$wealth_tertile_score,
  levels = c(2, 1, 0),
  labels = c(
    "High wealth",
    "Middle wealth",
    "Low wealth"
  )
)

cesd_long$income_component <- factor(
  cesd_long$income_tertile_score,
  levels = c(2, 1, 0),
  labels = c(
    "High income",
    "Middle income",
    "Low income"
  )
)

cesd_long$education_component <- factor(
  as.character(cesd_long$education_3cat),
  levels = c("High", "Middle", "Low"),
  labels = c(
    "High education",
    "Middle education",
    "Low education"
  )
)

# 检查参考组和频数
table(cesd_long$wealth_component, useNA = "always")
table(cesd_long$income_component, useNA = "always")
table(cesd_long$education_component, useNA = "always")

# 三个组成成分共同纳入
model4_components_joint <- lmer(
  cesd_good_z ~
    wealth_component * time +
    income_component * time +
    education_component * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

summary(model4_components_joint)

model4_components_tidy <- broom.mixed::tidy(
  model4_components_joint,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model4_components_tidy


# 需要保留的模型项
model4_keep_terms <- c(
  "wealth_componentMiddle wealth",
  "wealth_componentLow wealth",
  "income_componentMiddle income",
  "income_componentLow income",
  "education_componentMiddle education",
  "education_componentLow education",
  "time",
  "wealth_componentMiddle wealth:time",
  "wealth_componentLow wealth:time",
  "income_componentMiddle incomtime",
  "income_componentLow incomtime",
  "education_componentMiddle education:time",
  "education_componentLow education:time"
)

model4_component_result <- model4_components_tidy[
  model4_components_tidy$term %in% model4_keep_terms,
]

# 按目标顺序排列
model4_component_result <- model4_component_result[
  match(
    model4_keep_terms,
    model4_component_result$term
  ),
]

# 显示名称
model4_labels <- c(
  "Middle wealth",
  "Low wealth",
  "Middle income",
  "Low income",
  "Middle education",
  "Low education",
  "Time, per 2-year interval",
  "Middle wealth × time",
  "Low wealth × time",
  "Middle income × time",
  "Low income × time",
  "Middle education × time",
  "Low education × time"
)

model4_component_result$Variable <- model4_labels

model4_component_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model4_component_result$estimate,
  model4_component_result$conf.low,
  model4_component_result$conf.high
)

model4_component_result$P <- ifelse(
  model4_component_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model4_component_result$p.value)
)

model4_component_result <- model4_component_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model4_component_result


####FI和CESD####
#Panel B采用二分类食品不安全作为暴露。先确保0为食品安全/边缘食品安全参考组，1为食品不安全。
cesd_long$fi_USDA6_binary <- factor(
  cesd_long$fi_USDA6_binary,
  levels = c(0, 1)
)

#***********************************Model 1***********************************#
# Panel B, Model 1
model1_FI_CESD <- lmer(
  cesd_good_z ~
    fi_USDA6_binary +
    time +
    cesd_good_z_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Panel B Model 1结果
model1_FI_result <- tidy(
  model1_FI_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model1_FI_result <- model1_FI_result[
  model1_FI_result$term %in% c(
    "fi_USDA6_binary1",
    "time"
  ),
]

model1_FI_result <- model1_FI_result[
  match(
    c(
      "fi_USDA6_binary1",
      "time"
    ),
    model1_FI_result$term
  ),
]

model1_FI_result$Variable <- c(
  "Food insecurity",
  "Time, per 2-year interval"
)

model1_FI_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model1_FI_result$estimate,
  model1_FI_result$conf.low,
  model1_FI_result$conf.high
)

model1_FI_result$P <- ifelse(
  model1_FI_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model1_FI_result$p.value)
)

model1_FI_result <- model1_FI_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model1_FI_result

#***********************************Model 2***********************************#
# Panel B, Model 2：加入核心混杂因素
model2_FI_CESD <- lmer(
  cesd_good_z ~
    fi_USDA6_binary +
    time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Panel B Model 2结果
model2_FI_result <- tidy(
  model2_FI_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model2_FI_result <- model2_FI_result[
  model2_FI_result$term %in% c(
    "fi_USDA6_binary1",
    "time"
  ),
]

model2_FI_result <- model2_FI_result[
  match(
    c(
      "fi_USDA6_binary1",
      "time"
    ),
    model2_FI_result$term
  ),
]

model2_FI_result$Variable <- c(
  "Food insecurity",
  "Time, per 2-year interval"
)

model2_FI_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model2_FI_result$estimate,
  model2_FI_result$conf.low,
  model2_FI_result$conf.high
)

model2_FI_result$P <- ifelse(
  model2_FI_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model2_FI_result$p.value)
)

model2_FI_result <- model2_FI_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model2_FI_result


#***********************************Model 3***********************************#
# Panel B, Model 3：加入食品不安全与时间的交互项
model3_FI_CESD <- lmer(
  cesd_good_z ~
    fi_USDA6_binary * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Panel B Model 3结果
model3_FI_result <- tidy(
  model3_FI_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model3_FI_result <- model3_FI_result[
  model3_FI_result$term %in% c(
    "fi_USDA6_binary1",
    "time",
    "fi_USDA6_binary1:time"
  ),
]

model3_FI_result <- model3_FI_result[
  match(
    c(
      "fi_USDA6_binary1",
      "time",
      "fi_USDA6_binary1:time"
    ),
    model3_FI_result$term
  ),
]

model3_FI_result$Variable <- c(
  "Food insecurity",
  "Time, per 2-year interval",
  "Food insecurity × time"
)

model3_FI_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model3_FI_result$estimate,
  model3_FI_result$conf.low,
  model3_FI_result$conf.high
)

model3_FI_result$P <- ifelse(
  model3_FI_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model3_FI_result$p.value)
)

model3_FI_result <- model3_FI_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model3_FI_result

anova(model2_FI_CESD, model3_FI_CESD)


#***********************************Model 4***********************************#
# Panel B, Model 4：加入扩展调整变量
model4_FI_CESD <- lmer(
  cesd_good_z ~
    fi_USDA6_binary +
    time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Panel B Model 4结果
model4_FI_result <- tidy(
  model4_FI_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model4_FI_result <- model4_FI_result[
  model4_FI_result$term %in% c(
    "fi_USDA6_binary1",
    "time"
  ),
]

model4_FI_result <- model4_FI_result[
  match(
    c(
      "fi_USDA6_binary1",
      "time"
    ),
    model4_FI_result$term
  ),
]

model4_FI_result$Variable <- c(
  "Food insecurity",
  "Time, per 2-year interval"
)

model4_FI_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model4_FI_result$estimate,
  model4_FI_result$conf.low,
  model4_FI_result$conf.high
)

model4_FI_result$P <- ifelse(
  model4_FI_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model4_FI_result$p.value)
)

model4_FI_result <- model4_FI_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model4_FI_result


#***********************************Model 5***********************************#
# Panel B, Model 5：在扩展调整模型中加入食品不安全与时间的交互项
model5_FI_CESD <- lmer(
  cesd_good_z ~
    fi_USDA6_binary * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化Panel B Model 5结果
model5_FI_result <- tidy(
  model5_FI_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model5_FI_result <- model5_FI_result[
  model5_FI_result$term %in% c(
    "fi_USDA6_binary1",
    "time",
    "fi_USDA6_binary1:time"
  ),
]

model5_FI_result <- model5_FI_result[
  match(
    c(
      "fi_USDA6_binary1",
      "time",
      "fi_USDA6_binary1:time"
    ),
    model5_FI_result$term
  ),
]

model5_FI_result$Variable <- c(
  "Food insecurity",
  "Time, per 2-year interval",
  "Food insecurity × time"
)

model5_FI_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model5_FI_result$estimate,
  model5_FI_result$conf.low,
  model5_FI_result$conf.high
)

model5_FI_result$P <- ifelse(
  model5_FI_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model5_FI_result$p.value)
)

model5_FI_result <- model5_FI_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model5_FI_result


#***********************************补充分析***********************************#
# 补充分析：FI严重程度评分，Model 1
model1_FIscore_CESD <- lmer(
  cesd_good_z ~
    fi_score_0_6 +
    time +
    cesd_good_z_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)
# 格式化结果
model1_FIscore_result <- tidy(
  model1_FIscore_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model1_FIscore_result <- model1_FIscore_result[
  model1_FIscore_result$term %in% c(
    "fi_score_0_6",
    "time"
  ),
]

model1_FIscore_result <- model1_FIscore_result[
  match(
    c("fi_score_0_6", "time"),
    model1_FIscore_result$term
  ),
]

model1_FIscore_result$Variable <- c(
  "Food insecurity score, per 1-point increase",
  "Time, per 2-year interval"
)

model1_FIscore_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model1_FIscore_result$estimate,
  model1_FIscore_result$conf.low,
  model1_FIscore_result$conf.high
)

model1_FIscore_result$P <- ifelse(
  model1_FIscore_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model1_FIscore_result$p.value)
)

model1_FIscore_result <- model1_FIscore_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model1_FIscore_result


# 补充分析：FI严重程度评分，Model 2
model2_FIscore_CESD <- lmer(
  cesd_good_z ~
    fi_score_0_6 +
    time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化结果
model2_FIscore_result <- tidy(
  model2_FIscore_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model2_FIscore_result <- model2_FIscore_result[
  model2_FIscore_result$term %in% c(
    "fi_score_0_6",
    "time"
  ),
]

model2_FIscore_result <- model2_FIscore_result[
  match(
    c("fi_score_0_6", "time"),
    model2_FIscore_result$term
  ),
]

model2_FIscore_result$Variable <- c(
  "Food insecurity score, per 1-point increase",
  "Time, per 2-year interval"
)

model2_FIscore_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model2_FIscore_result$estimate,
  model2_FIscore_result$conf.low,
  model2_FIscore_result$conf.high
)

model2_FIscore_result$P <- ifelse(
  model2_FIscore_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model2_FIscore_result$p.value)
)

model2_FIscore_result <- model2_FIscore_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model2_FIscore_result

# 补充分析：FI严重程度评分，Model 3
model3_FIscore_CESD <- lmer(
  cesd_good_z ~
    fi_score_0_6 * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化结果
model3_FIscore_result <- tidy(
  model3_FIscore_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model3_FIscore_result <- model3_FIscore_result[
  model3_FIscore_result$term %in% c(
    "fi_score_0_6",
    "time",
    "fi_score_0_6:time"
  ),
]

model3_FIscore_result <- model3_FIscore_result[
  match(
    c(
      "fi_score_0_6",
      "time",
      "fi_score_0_6:time"
    ),
    model3_FIscore_result$term
  ),
]

model3_FIscore_result$Variable <- c(
  "Food insecurity score, per 1-point increase",
  "Time, per 2-year interval",
  "Food insecurity score × time"
)

model3_FIscore_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model3_FIscore_result$estimate,
  model3_FIscore_result$conf.low,
  model3_FIscore_result$conf.high
)

model3_FIscore_result$P <- ifelse(
  model3_FIscore_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model3_FIscore_result$p.value)
)

model3_FIscore_result <- model3_FIscore_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model3_FIscore_result


# 补充分析：FI严重程度评分，Model 4
model4_FIscore_CESD <- lmer(
  cesd_good_z ~
    fi_score_0_6 +
    time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化结果
model4_FIscore_result <- tidy(
  model4_FIscore_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model4_FIscore_result <- model4_FIscore_result[
  model4_FIscore_result$term %in% c(
    "fi_score_0_6",
    "time"
  ),
]

model4_FIscore_result <- model4_FIscore_result[
  match(
    c("fi_score_0_6", "time"),
    model4_FIscore_result$term
  ),
]

model4_FIscore_result$Variable <- c(
  "Food insecurity score, per 1-point increase",
  "Time, per 2-year interval"
)

model4_FIscore_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model4_FIscore_result$estimate,
  model4_FIscore_result$conf.low,
  model4_FIscore_result$conf.high
)

model4_FIscore_result$P <- ifelse(
  model4_FIscore_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model4_FIscore_result$p.value)
)

model4_FIscore_result <- model4_FIscore_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model4_FIscore_result


# 补充分析：FI严重程度评分，Model 5
model5_FIscore_CESD <- lmer(
  cesd_good_z ~
    fi_score_0_6 * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

# 格式化结果
model5_FIscore_result <- tidy(
  model5_FIscore_CESD,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
)

model5_FIscore_result <- model5_FIscore_result[
  model5_FIscore_result$term %in% c(
    "fi_score_0_6",
    "time",
    "fi_score_0_6:time"
  ),
]

model5_FIscore_result <- model5_FIscore_result[
  match(
    c(
      "fi_score_0_6",
      "time",
      "fi_score_0_6:time"
    ),
    model5_FIscore_result$term
  ),
]

model5_FIscore_result$Variable <- c(
  "Food insecurity score, per 1-point increase",
  "Time, per 2-year interval",
  "Food insecurity score × time"
)

model5_FIscore_result$Beta_95CI <- sprintf(
  "%.3f (%.3f, %.3f)",
  model5_FIscore_result$estimate,
  model5_FIscore_result$conf.low,
  model5_FIscore_result$conf.high
)

model5_FIscore_result$P <- ifelse(
  model5_FIscore_result$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", model5_FIscore_result$p.value)
)

model5_FIscore_result <- model5_FIscore_result[
  ,
  c("Variable", "Beta_95CI", "P")
]

model5_FIscore_result


####绘制SES和时间交互的图####
library(emmeans)
library(ggplot2)
library(patchwork)
model_SEScat_time <- lmer(
  cesd_good_z ~
    ses_3cat * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)
summary(model_SEScat_time)
#检查交互是否显著
# 删除SES × time交互项，保留SES和time主效应
model_SEScat_no_interaction <- update(
  model_SEScat_time,
  formula = . ~ . - ses_3cat:time,
  REML = FALSE
)

# 确认两个模型使用相同观测数
c(
  reduced_model = nobs(model_SEScat_no_interaction),
  full_model = nobs(model_SEScat_time)
)

# 似然比检验
panelA_interaction_test <- anova(
  model_SEScat_no_interaction,
  model_SEScat_time
)

panelA_interaction_test

#检查模型真正使用的参考组，便于解释结果
levels(
  model.frame(model_SEScat_time)$ses_3cat
)

#提取两个单独交互项：低VS高和中VS高
panelA_coefficients <- broom.mixed::tidy(
  model_SEScat_time,
  effects = "fixed",
  conf.int = TRUE
)

panelA_coefficients[
  grepl(
    "ses_3cat.*:time",
    panelA_coefficients$term
  ),
  c(
    "term",
    "estimate",
    "conf.low",
    "conf.high",
    "p.value"
  )
]


#提取调整后预测值，为后面画图准备
# Panel A：SES × time
pred_A <- as.data.frame(
  emmeans(
    model_SEScat_time,
    ~ ses_3cat * time,
    at = list(time = 0:3)
  )
)

pred_A$year <- 2014 + 2 * pred_A$time

pred_A$ses_3cat <- factor(
  pred_A$ses_3cat,
  levels = c("Low SES", "Middle SES", "High SES")
)

pred_A$lower <- pred_A$emmean -
  qt(0.975, pred_A$df) * pred_A$SE

pred_A$upper <- pred_A$emmean +
  qt(0.975, pred_A$df) * pred_A$SE



plot_A_interaction <- ggplot(
  pred_A,
  aes(
    x = year,
    y = emmean,
    color = ses_3cat,
    fill = ses_3cat,
    group = ses_3cat
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.14,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c(
      "Low SES" = "#2B83BA",
      "Middle SES" = "#C75D52",
      "High SES" = "#E6A14A"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Low SES" = "#2B83BA",
      "Middle SES" = "#C75D52",
      "High SES" = "#E6A14A"
    )
  ) +
  scale_x_continuous(
    breaks = c(2014, 2016, 2018, 2020)
  ) +
  labs(
    x = NULL,
    y = "Adjusted psychological health score"
  ) +
  interaction_theme



# 1. 从hrs_ses重新加入正确的FI二分类变量
id_match <- match(
  as.character(cesd_long$HHIDPN),
  as.character(hrs_ses$HHIDPN)
)

stopifnot(sum(is.na(id_match)) == 0)

fi_original <- as.character(
  hrs_ses$fi_USDA6_binary[id_match]
)

# 原始变量必须只有0和1
table(fi_original, useNA = "always")
stopifnot(all(fi_original %in% c("0", "1")))

cesd_long$fi_group <- factor(
  fi_original,
  levels = c("0", "1"),
  labels = c("Food secure", "Food insecure")
)

# 高SES作为模型参考组
cesd_long$ses_3cat <- factor(
  as.character(cesd_long$ses_3cat),
  levels = c(
    "High SES",
    "Middle SES",
    "Low SES"
  )
)

table(cesd_long$fi_group, useNA = "always")
table(cesd_long$ses_3cat, useNA = "always")


# 2. Panel B模型：FI × time
model_FI_time <- lmer(
  cesd_good_z ~
    fi_group * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

#查看交互项的结果
anova(
  model_FI_time,
  type = 3,
  ddf = "Satterthwaite"
)["fi_group:time", ]

result_B <- broom.mixed::tidy(
  model_FI_time,
  effects = "fixed",
  conf.int = TRUE
)

result_B[
  grepl("fi_group.*:time", result_B$term),
  c("term", "estimate", "conf.low", "conf.high", "p.value")
]

#检查一下参考组
levels(model.frame(model_FI_time)$fi_group)


# 3. Panel C模型：SES × FI × time
model_SES_FI_time <- lmer(
  cesd_good_z ~
    ses_3cat * fi_group * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    self_rated_health_2012 +
    chronic_condition_count_2012 +
    insurance_3cat +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

nobs(model_FI_time)
nobs(model_SES_FI_time)


#检查panelC的交互效应
anova(
  model_SES_FI_time,
  type = 3,
  ddf = "Satterthwaite"
)["ses_3cat:fi_group:time", ]

result_C <- broom.mixed::tidy(
  model_SES_FI_time,
  effects = "fixed",
  conf.int = TRUE
)

result_C[
  grepl(
    "ses_3cat.*:fi_group.*:time",
    result_C$term
  ),
  c(
    "term",
    "estimate",
    "conf.low",
    "conf.high",
    "p.value"
  )
]

#确认参考组
levels(model.frame(model_SES_FI_time)$ses_3cat)
levels(model.frame(model_SES_FI_time)$fi_group)

# 4. 生成Panel B预测值
pred_B <- as.data.frame(
  emmeans(
    model_FI_time,
    ~ fi_group * time,
    at = list(time = 0:3)
  )
)

pred_B$year <- 2014 + 2 * pred_B$time

pred_B$lower <- pred_B$emmean -
  qt(0.975, pred_B$df) * pred_B$SE

pred_B$upper <- pred_B$emmean +
  qt(0.975, pred_B$df) * pred_B$SE

table(pred_B$fi_group, useNA = "always")


# 5. 生成Panel C预测值
pred_C <- as.data.frame(
  emmeans(
    model_SES_FI_time,
    ~ ses_3cat * fi_group * time,
    at = list(time = 0:3)
  )
)

pred_C$year <- 2014 + 2 * pred_C$time

pred_C$lower <- pred_C$emmean -
  qt(0.975, pred_C$df) * pred_C$SE

pred_C$upper <- pred_C$emmean +
  qt(0.975, pred_C$df) * pred_C$SE

# 调整绘图中的SES顺序
pred_C$ses_3cat <- factor(
  pred_C$ses_3cat,
  levels = c(
    "Low SES",
    "Middle SES",
    "High SES"
  )
)

table(pred_C$fi_group, useNA = "always")
table(pred_C$ses_3cat, useNA = "always")


# 6. 统一绘图主题
interaction_theme <- theme_classic(base_size = 9) +
  theme(
    panel.grid.major.y = element_line(
      color = "#E8E8E8",
      linewidth = 0.35
    ),
    axis.text = element_text(
      color = "#404040"
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.width = grid::unit(0.7, "cm"),
    plot.margin = margin(5, 7, 5, 5)
  )


# 7. 绘制Panel B
plot_B_interaction <- ggplot(
  pred_B,
  aes(
    x = year,
    y = emmean,
    color = fi_group,
    fill = fi_group,
    group = fi_group
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.14,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Food secure" = "#315D78",
      "Food insecure" = "#B65C5C"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Food secure" = "#315D78",
      "Food insecure" = "#B65C5C"
    )
  ) +
  scale_x_continuous(
    breaks = c(2014, 2016, 2018, 2020)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  interaction_theme


# 8. 绘制Panel C
plot_C_interaction <- ggplot(
  pred_C,
  aes(
    x = year,
    y = emmean,
    color = fi_group,
    fill = fi_group,
    linetype = fi_group,
    shape = fi_group,
    group = fi_group
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.13,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 0.75
  ) +
  geom_point(
    size = 1.7
  ) +
  facet_wrap(
    ~ ses_3cat,
    nrow = 1
  ) +
  scale_color_manual(
    values = c(
      "Food secure" = "#315D78",
      "Food insecure" = "#B65C5C"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Food secure" = "#315D78",
      "Food insecure" = "#B65C5C"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Food secure" = "solid",
      "Food insecure" = "dashed"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Food secure" = 16,
      "Food insecure" = 17
    )
  ) +
  scale_x_continuous(
    breaks = c(2014, 2016, 2018, 2020)
  ) +
  labs(
    x = "Year",
    y = "Adjusted psychological health score"
  ) +
  interaction_theme +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 9,
      face = "plain"
    ),
    legend.position = "top"
  )


# 9. 与已经生成的Panel A合并
figure3 <- (
  plot_A_interaction |
    plot_B_interaction
) /
  plot_C_interaction +
  plot_layout(
    heights = c(1, 1.25)
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      size = 12,
      face = "plain"
    )
  )

figure3
####反事实推断####
#第1步：拟合最终的GEE反事实模型
library(geepack)

# SES参考组：High SES
cesd_gee$ses_3cat <- factor(
  cesd_gee$ses_3cat,
  levels = c(
    "High SES",
    "Middle SES",
    "Low SES"
  )
)

# FI参考组：Food secure
cesd_gee$fi_group <- factor(
  cesd_gee$fi_group,
  levels = c(
    "Food secure",
    "Food insecure"
  )
)

# 随访年份作为分类变量
cesd_gee$year_factor <- factor(
  2014 + 2 * cesd_gee$time,
  levels = c(2014, 2016, 2018, 2020)
)

# 按参与者和时间排序
cesd_gee <- cesd_gee[
  order(cesd_gee$HHIDPN, cesd_gee$time),
]

# 检查分类
table(cesd_gee$ses_3cat, useNA = "always")
table(cesd_gee$fi_group, useNA = "always")
table(cesd_gee$year_factor, useNA = "always")

# 最终GEE模型
model_cf_GEE <- geeglm(
  cesd_good_z ~
    ses_3cat * fi_group * year_factor +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012,
  id = HHIDPN,
  waves = time,
  data = cesd_gee,
  family = gaussian(link = "identity"),
  corstr = "exchangeable",
  std.err = "san.se"
)

summary(model_cf_GEE)

# 工作相关系数
model_cf_GEE$geese$alpha

# 确认样本量
nrow(cesd_gee)
length(unique(cesd_gee$HHIDPN))

#第2步：建立三种反事实情景并生成预测值
# 每名参与者保留一行基线信息
baseline_cf_GEE <- cesd_gee[
  !duplicated(cesd_gee$HHIDPN),
  c(
    "HHIDPN",
    "ses_3cat",
    "fi_group",
    "cesd_good_z_2012",
    "age_2012",
    "sex_1male_0female",
    "race_ethnicity_4cat",
    "partnered_2012",
    "smoking_now_2012",
    "drinking_any_2012",
    "bmi_2012"
  )
]

nrow(baseline_cf_GEE)
# 应为6695


# 扩展为每人4个随访年份
cf_observed_GEE <- baseline_cf_GEE[
  rep(
    seq_len(nrow(baseline_cf_GEE)),
    each = 4
  ),
]

cf_observed_GEE$year <- rep(
  c(2014, 2016, 2018, 2020),
  times = nrow(baseline_cf_GEE)
)

cf_observed_GEE$time <- rep(
  0:3,
  times = nrow(baseline_cf_GEE)
)

cf_observed_GEE$year_factor <- factor(
  cf_observed_GEE$year,
  levels = c(2014, 2016, 2018, 2020)
)


# 情景1：保持每名参与者实际观察到的FI状态
cf_observed_GEE$predicted_observed <- predict(
  model_cf_GEE,
  newdata = cf_observed_GEE,
  type = "response"
)


# 情景2：所有参与者均设为食品安全
cf_noFI_GEE <- cf_observed_GEE

cf_noFI_GEE$fi_group <- factor(
  rep(
    "Food secure",
    nrow(cf_noFI_GEE)
  ),
  levels = c(
    "Food secure",
    "Food insecure"
  )
)

cf_noFI_GEE$predicted_noFI <- predict(
  model_cf_GEE,
  newdata = cf_noFI_GEE,
  type = "response"
)


# 情景3：所有参与者均设为食品不安全
cf_allFI_GEE <- cf_observed_GEE

cf_allFI_GEE$fi_group <- factor(
  rep(
    "Food insecure",
    nrow(cf_allFI_GEE)
  ),
  levels = c(
    "Food secure",
    "Food insecure"
  )
)

cf_allFI_GEE$predicted_allFI <- predict(
  model_cf_GEE,
  newdata = cf_allFI_GEE,
  type = "response"
)

#检查数据结果
# 应为6695 × 4 = 26780
nrow(cf_observed_GEE)
nrow(cf_noFI_GEE)
nrow(cf_allFI_GEE)

# 预测值不应缺失
sum(is.na(cf_observed_GEE$predicted_observed))
sum(is.na(cf_noFI_GEE$predicted_noFI))
sum(is.na(cf_allFI_GEE$predicted_allFI))

# 每个年份应有6695行
table(cf_observed_GEE$year)

# 三种情景中的FI分布
table(cf_observed_GEE$fi_group)
table(cf_noFI_GEE$fi_group)
table(cf_allFI_GEE$fi_group)

# 简单查看预测值范围
summary(cf_observed_GEE$predicted_observed)
summary(cf_noFI_GEE$predicted_noFI)
summary(cf_allFI_GEE$predicted_allFI)

#第3步先计算图中所有指标的点估计
# 合并三个反事实情景的预测结果
cf_results_GEE <- cf_observed_GEE

cf_results_GEE$predicted_noFI <-
  cf_noFI_GEE$predicted_noFI

cf_results_GEE$predicted_allFI <-
  cf_allFI_GEE$predicted_allFI


# PIB：全部食品安全减去观察情景
cf_results_GEE$PIB <-
  cf_results_GEE$predicted_noFI -
  cf_results_GEE$predicted_observed


# ATE/CATE：全部食品安全减去全部食品不安全
cf_results_GEE$secure_vs_insecure <-
  cf_results_GEE$predicted_noFI -
  cf_results_GEE$predicted_allFI

#计算总体PIB
PIB_overall_GEE <- aggregate(
  PIB ~ year,
  data = cf_results_GEE,
  FUN = mean
)

PIB_overall_GEE

#计算各SES组、各年份PIB，对应图中Panel A
PIB_by_SES_year_GEE <- aggregate(
  PIB ~ ses_3cat + year,
  data = cf_results_GEE,
  FUN = mean
)

PIB_by_SES_year_GEE

#计算各SES组跨4个年份的平均PIB，对应Panel C
PIB_average_by_SES_GEE <- aggregate(
  PIB ~ ses_3cat,
  data = cf_results_GEE,
  FUN = mean
)

PIB_average_by_SES_GEE

#计算总体ATE
ATE_overall_GEE <- aggregate(
  secure_vs_insecure ~ year,
  data = cf_results_GEE,
  FUN = mean
)

names(ATE_overall_GEE)[2] <- "ATE"

ATE_overall_GEE

#计算各SES组、各年份CATE，对应Panel B
CATE_by_SES_year_GEE <- aggregate(
  secure_vs_insecure ~ ses_3cat + year,
  data = cf_results_GEE,
  FUN = mean
)

names(CATE_by_SES_year_GEE)[3] <- "CATE"

CATE_by_SES_year_GEE

#计算低SES与高SES之间的CATE差异，对应Panel D
CATE_high_GEE <- CATE_by_SES_year_GEE[
  CATE_by_SES_year_GEE$ses_3cat == "High SES",
  c("year", "CATE")
]

names(CATE_high_GEE)[2] <- "CATE_high"

CATE_low_GEE <- CATE_by_SES_year_GEE[
  CATE_by_SES_year_GEE$ses_3cat == "Low SES",
  c("year", "CATE")
]

names(CATE_low_GEE)[2] <- "CATE_low"

CATE_difference_GEE <- merge(
  CATE_low_GEE,
  CATE_high_GEE,
  by = "year"
)

# 正值表示FI对低SES组的不利影响更强
CATE_difference_GEE$CATE_difference <-
  CATE_difference_GEE$CATE_low -
  CATE_difference_GEE$CATE_high

CATE_difference_GEE

#最后检查所有结果
PIB_overall_GEE
PIB_by_SES_year_GEE
PIB_average_by_SES_GEE
ATE_overall_GEE
CATE_by_SES_year_GEE
CATE_difference_GEE


#第4步统一计算图A–D所需的95% CI和P值
# GEE系数和稳健协方差矩阵
beta_GEE <- coef(model_cf_GEE)
vcov_GEE <- vcov(model_cf_GEE)

coef_names <- names(beta_GEE)

# 建立三个情景的设计矩阵
gee_terms <- delete.response(
  terms(model_cf_GEE)
)

X_observed <- model.matrix(
  gee_terms,
  data = cf_observed_GEE
)

X_noFI <- model.matrix(
  gee_terms,
  data = cf_noFI_GEE
)

X_allFI <- model.matrix(
  gee_terms,
  data = cf_allFI_GEE
)

# 与模型系数顺序保持一致
X_observed <- X_observed[
  , coef_names, drop = FALSE
]

X_noFI <- X_noFI[
  , coef_names, drop = FALSE
]

X_allFI <- X_allFI[
  , coef_names, drop = FALSE
]

vcov_GEE <- vcov_GEE[
  coef_names,
  coef_names,
  drop = FALSE
]

#定义delta-method函数
delta_result <- function(L) {
  
  estimate <- as.numeric(
    L %*% beta_GEE
  )
  
  se <- sqrt(
    as.numeric(
      L %*% vcov_GEE %*% L
    )
  )
  
  p_value <- 2 * pnorm(
    abs(estimate / se),
    lower.tail = FALSE
  )
  
  c(
    estimate = estimate,
    SE = se,
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    P = p_value
  )
}

#Panel A：各SES组PIB
panelA_GEE <- data.frame()

for (current_year in c(2014, 2016, 2018, 2020)) {
  
  for (current_SES in c(
    "Low SES",
    "Middle SES",
    "High SES"
  )) {
    
    index <-
      cf_observed_GEE$year == current_year &
      cf_observed_GEE$ses_3cat == current_SES
    
    L <- colMeans(
      X_noFI[index, , drop = FALSE] -
        X_observed[index, , drop = FALSE]
    )
    
    result <- delta_result(L)
    
    panelA_GEE <- rbind(
      panelA_GEE,
      data.frame(
        year = current_year,
        ses_3cat = current_SES,
        estimate = unname(result["estimate"]),
        SE = unname(result["SE"]),
        lower = unname(result["lower"]),
        upper = unname(result["upper"]),
        P = unname(result["P"])
      )
    )
  }
}

panelA_GEE

#Panel B：各SES组CATE
panelB_GEE <- data.frame()

for (current_year in c(2014, 2016, 2018, 2020)) {
  
  for (current_SES in c(
    "Low SES",
    "Middle SES",
    "High SES"
  )) {
    
    index <-
      cf_observed_GEE$year == current_year &
      cf_observed_GEE$ses_3cat == current_SES
    
    L <- colMeans(
      X_noFI[index, , drop = FALSE] -
        X_allFI[index, , drop = FALSE]
    )
    
    result <- delta_result(L)
    
    panelB_GEE <- rbind(
      panelB_GEE,
      data.frame(
        year = current_year,
        ses_3cat = current_SES,
        estimate = unname(result["estimate"]),
        SE = unname(result["SE"]),
        lower = unname(result["lower"]),
        upper = unname(result["upper"]),
        P = unname(result["P"])
      )
    )
  }
}

panelB_GEE


#Panel C：各SES组跨4年的平均PIB
panelC_GEE <- data.frame()

for (current_SES in c(
  "Low SES",
  "Middle SES",
  "High SES"
)) {
  
  index <-
    cf_observed_GEE$ses_3cat == current_SES
  
  L <- colMeans(
    X_noFI[index, , drop = FALSE] -
      X_observed[index, , drop = FALSE]
  )
  
  result <- delta_result(L)
  
  panelC_GEE <- rbind(
    panelC_GEE,
    data.frame(
      ses_3cat = current_SES,
      estimate = unname(result["estimate"]),
      SE = unname(result["SE"]),
      lower = unname(result["lower"]),
      upper = unname(result["upper"]),
      P = unname(result["P"])
    )
  )
}

panelC_GEE

#Panel D：低SES与高SES的CATE差异
panelD_GEE <- data.frame()

# 保存4个年份的对比向量，用于总体Wald检验
R_CATE_difference <- matrix(
  NA_real_,
  nrow = 4,
  ncol = length(beta_GEE)
)

colnames(R_CATE_difference) <- names(beta_GEE)
rownames(R_CATE_difference) <- c(
  "2014", "2016", "2018", "2020"
)

years <- c(2014, 2016, 2018, 2020)

for (j in seq_along(years)) {
  
  current_year <- years[j]
  
  index_low <-
    cf_observed_GEE$year == current_year &
    cf_observed_GEE$ses_3cat == "Low SES"
  
  index_high <-
    cf_observed_GEE$year == current_year &
    cf_observed_GEE$ses_3cat == "High SES"
  
  L_low <- colMeans(
    X_noFI[index_low, , drop = FALSE] -
      X_allFI[index_low, , drop = FALSE]
  )
  
  L_high <- colMeans(
    X_noFI[index_high, , drop = FALSE] -
      X_allFI[index_high, , drop = FALSE]
  )
  
  # 正值：FI对低SES组的不利影响更强
  L_difference <- L_low - L_high
  
  R_CATE_difference[j, ] <- L_difference
  
  result <- delta_result(L_difference)
  
  panelD_GEE <- rbind(
    panelD_GEE,
    data.frame(
      year = current_year,
      estimate = unname(result["estimate"]),
      SE = unname(result["SE"]),
      lower = unname(result["lower"]),
      upper = unname(result["upper"]),
      P = unname(result["P"])
    )
  )
}

panelD_GEE


#计算Panel D的4自由度总体Wald检验
wald_estimate <- as.numeric(
  R_CATE_difference %*% beta_GEE
)

wald_covariance <-
  R_CATE_difference %*%
  vcov_GEE %*%
  t(R_CATE_difference)

wald_statistic <- as.numeric(
  t(wald_estimate) %*%
    solve(wald_covariance) %*%
    wald_estimate
)

wald_df <- qr(wald_covariance)$rank

wald_P <- pchisq(
  wald_statistic,
  df = wald_df,
  lower.tail = FALSE
)

wald_heterogeneity_GEE <- data.frame(
  Wald_chisq = wald_statistic,
  df = wald_df,
  P = wald_P
)

wald_heterogeneity_GEE

panelA_GEE
panelB_GEE
panelC_GEE
panelD_GEE
wald_heterogeneity_GEE

#第5步绘制四Panel图
library(ggplot2)
library(patchwork)

# SES顺序和颜色
ses_levels <- c(
  "Low SES",
  "Middle SES",
  "High SES"
)

ses_colors <- c(
  "Low SES" = "#C95F59",
  "Middle SES" = "#2F8C8C",
  "High SES" = "#6D6596"
)

panelA_GEE$ses_3cat <- factor(
  panelA_GEE$ses_3cat,
  levels = ses_levels
)

panelB_GEE$ses_3cat <- factor(
  panelB_GEE$ses_3cat,
  levels = ses_levels
)

panelC_GEE$ses_3cat <- factor(
  panelC_GEE$ses_3cat,
  levels = rev(ses_levels)
)

panelD_GEE$year_factor <- factor(
  panelD_GEE$year,
  levels = c(2020, 2018, 2016, 2014)
)

# 数值标签
panelA_GEE$label <- sprintf(
  "%.3f",
  panelA_GEE$estimate
)

panelB_GEE$label <- sprintf(
  "%.3f",
  panelB_GEE$estimate
)

panelC_GEE$label <- sprintf(
  "%.3f\n(%.3f, %.3f)",
  panelC_GEE$estimate,
  panelC_GEE$lower,
  panelC_GEE$upper
)

format_P <- function(x) {
  ifelse(
    x < 0.001,
    "P<0.001",
    sprintf("P=%.3f", x)
  )
}

panelD_GEE$label <- paste0(
  sprintf(
    "%.3f (%.3f, %.3f)",
    panelD_GEE$estimate,
    panelD_GEE$lower,
    panelD_GEE$upper
  ),
  "\n",
  format_P(panelD_GEE$P)
)

wald_label <- sprintf(
  "Overall heterogeneity across 4 years:\nWald \u03C7\u00B2(%d)=%.2f, P=%.3f",
  wald_heterogeneity_GEE$df,
  wald_heterogeneity_GEE$Wald_chisq,
  wald_heterogeneity_GEE$P
)

#设置绘图的主题
cf_theme <- theme_classic(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9.5),
    
    plot.tag = element_text(
      size = 14,
      face = "plain"
    ),
    
    plot.margin = margin(10, 18, 10, 10)
  )
ses_levels <- c(
  "Low SES",
  "Middle SES",
  "High SES"
)

panelA_GEE$ses_3cat <- factor(
  panelA_GEE$ses_3cat,
  levels = ses_levels
)

panelB_GEE$ses_3cat <- factor(
  panelB_GEE$ses_3cat,
  levels = ses_levels
)

panelA_GEE$label_vjust <- ifelse(
  panelA_GEE$ses_3cat == "High SES",
  1.5,
  -0.8
)

panelB_GEE$label_vjust <- ifelse(
  panelB_GEE$ses_3cat == "Middle SES",
  1.5,
  -0.8
)

# 三组在年份两侧错开的宽度
ses_dodge <- position_dodge(width = 0.48)

#Panel A
plot_A <- ggplot(
  panelA_GEE,
  aes(
    x = year,
    y = estimate,
    color = ses_3cat,
    group = ses_3cat
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey65"
  ) +
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.10,
    linewidth = 0.45,
    position = ses_dodge
  ) +
  geom_line(
    linewidth = 0.55,
    position = ses_dodge
  ) +
  geom_point(
    size = 1.8,
    position = ses_dodge
  ) +
  geom_text(
    aes(
      label = label,
      vjust = label_vjust
    ),
    position = ses_dodge,
    size = 4.0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = ses_colors,
    limits = ses_levels,
    breaks = ses_levels
  ) +
  scale_x_continuous(
    breaks = c(2014, 2016, 2018, 2020),
    limits = c(2013.35, 2020.65)
  ) +
  labs(
    x = "Year",
    y = "PIB (SD)",
    tag = "A"
  ) +
  cf_theme

#Panel B
plot_B <- ggplot(
  panelB_GEE,
  aes(
    x = year,
    y = estimate,
    color = ses_3cat,
    group = ses_3cat
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey65"
  ) +
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.10,
    linewidth = 0.45,
    position = ses_dodge
  ) +
  geom_line(
    linewidth = 0.55,
    position = ses_dodge
  ) +
  geom_point(
    size = 1.8,
    position = ses_dodge
  ) +
  geom_text(
    aes(
      label = label,
      vjust = label_vjust
    ),
    position = ses_dodge,
    size = 4.0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = ses_colors,
    limits = ses_levels,
    breaks = ses_levels
  ) +
  scale_x_continuous(
    breaks = c(2014, 2016, 2018, 2020),
    limits = c(2013.35, 2020.65)
  ) +
  labs(
    x = "Year",
    y = "CATE (SD)",
    tag = "B"
  ) +
  cf_theme

#Panel C
panelC_text_x <- max(panelC_GEE$upper) + 0.025

plot_C <- ggplot(
  panelC_GEE,
  aes(
    x = estimate,
    y = ses_3cat,
    color = ses_3cat
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey65"
  ) +
  geom_errorbar(
    aes(
      xmin = lower,
      xmax = upper
    ),
    orientation = "y",
    width = 0.12,
    linewidth = 0.55,
    color = "#C95F59"
  ) +
  geom_point(size = 2.2) +
  geom_text(
    aes(
      x = panelC_text_x,
      label = label
    ),
    hjust = 0,
    size = 3.6,
    show.legend = FALSE
  ) +
  scale_color_manual(values = ses_colors) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.05, 0.65)
    )
  ) +
  labs(
    x = "Average PIB (SD)",
    y = NULL,
    tag = "C"
  ) +
  cf_theme +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")


#Panel D
panelD_text_x <- max(panelD_GEE$upper) + 0.06

plot_D <- ggplot(
  panelD_GEE,
  aes(
    x = estimate,
    y = year_factor
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey65"
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.12,
    linewidth = 0.55,
    color = "#C95F59"
  ) +
  geom_point(
    size = 2.2,
    color = "#C95F59"
  ) +
  geom_text(
    aes(
      x = panelD_text_x,
      label = label
    ),
    hjust = 0,
    size = 3.6
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.05, 0.75)
    )
  ) +
  labs(
    x = "Difference in CATE (SD)",
    y = "Year",
    tag = "D",
    caption = wald_label
  ) +
  cf_theme +
  theme(
    legend.position = "none",
    plot.caption = element_text(
      hjust = 0.5,
      size = 10
    )
  ) +
  coord_cartesian(clip = "off")

#将四组panel组合
figure4_GEE <- (
  plot_A + plot_B
) / (
  plot_C + plot_D
) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

figure4_GEE


####检验三阶交互####
#建立包含全部两阶交互、但不包含三阶交互的模型
model_no_threeway <- lmer(
  cesd_good_z ~
    ses_3cat * fi_group * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)
summary(model_no_threeway)
anova(
  model_no_threeway,
  model_cf_primary
)

#检验总体SES异质性
#建立不包含SES × FI及三阶交互的模型
model_no_SES_heterogeneity <- lmer(
  cesd_good_z ~
    ses_3cat * time +
    fi_group * time +
    cesd_good_z_2012 +
    age_2012 +
    sex_1male_0female +
    race_ethnicity_4cat +
    partnered_2012 +
    smoking_now_2012 +
    drinking_any_2012 +
    bmi_2012 +
    (1 | HHIDPN),
  data = cesd_long,
  REML = FALSE,
  na.action = na.exclude
)

anova(
  model_no_SES_heterogeneity,
  model_cf_primary
)

