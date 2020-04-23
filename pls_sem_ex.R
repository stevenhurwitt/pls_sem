library(seminr)
library(plspm)
library(ggplot2)
library(reshape)


mobi = mobi
write.table(mobi, file = 'mobi.csv', sep = ",", row.names = F)

#### seminr ####

#create measurement model matrix
mobi_mm <- constructs(
  composite("Image",         multi_items("IMAG", 1:5), weights = mode_A),
  composite("Expectation",   multi_items("CUEX", 1:3), weights = mode_A),
  composite("Quality",       multi_items("PERQ", 1:7), weights = mode_A),
  composite("Value",         multi_items("PERV", 1:2), weights = mode_A),
  composite("Satisfaction", multi_items("CUSA", 1:3), weights = mode_A),
  composite("Complaints",   single_item("CUSCO"), weights = mode_A),
  higher_composite("HOC", c("Value", "Satisfaction"), orthogonal, mode_A),
  interaction_term(iv = "Image", moderator = "Expectation", method =  orthogonal, weights = mode_A),
  composite("Loyalty",      multi_items("CUSL", 1:3), weights = mode_A)
)

#create structural model
mobi_sm <- relationships(
  paths(from = "Image",        to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation",  to = c("Quality", "Value", "Satisfaction")),
  paths(from = "Quality",      to = c("Value", "Satisfaction")),
  paths(from = "Value",        to = c("Satisfaction")),
  paths(from = "Satisfaction", to = c("Complaints", "Loyalty")),
  paths(from = "Complaints",   to = "Loyalty")
)

#model estimation
mobi_pls <- estimate_pls(data = mobi,
                         measurement_model = mobi_mm,
                         structural_model = mobi_sm,
                         inner_weights = path_weighting)

summary(mobi_pls)

mobi_pls$path_coef
mobi_pls$outer_loadings
mobi_pls$outer_weights
mobi_pls$rSquared

#bootstrap for SE
boot_mobi_pls <- bootstrap_model(seminr_model = mobi_pls,
                                 nboot = 1000,
                                 cores = 2)

summary(boot_mobi_pls)


#### plspm ####

data(russett)

# path matrix (inner model realtionships)
AGRIN = c(0, 0, 0)
INDEV = c(0, 0, 0)
POLINS = c(1, 1, 0)
rus_path = rbind(AGRIN, INDEV, POLINS)

innerplot(rus_path)

# list indicating what variables are associated with what latent variables
rus_blocks = list(1:3, 4:5, 6:11)

# list indicating what variables are associated with what latent variables
rus_blocks = list(
  c("gini", "farm", "rent"),
  c("gnpr", "labo"),
  c("inst", "ecks", "death", "demostab", "demoinst", "dictator"))

rus_modes = rep("A", 3)

rus_pls = plspm(russett, rus_path, rus_blocks, modes = rus_modes)
summary(rus_pls)

rus_pls$inner_model
plot(rus_pls)

plot(rus_pls, what = "loadings", arr.width = 0.1)
plot(rus_pls, what = "weights", arr.width = 0.1)

xloads = melt(rus_pls$crossloadings, id.vars = c("name", "block"),variable_name = "LV")

ggplot(data = xloads,
       aes(x = name, y = value, fill = block)) +
  geom_hline(yintercept = 0, color = "gray75") +
  geom_hline(yintercept = c(-0.5, 0.5), color = "gray70", linetype = 2) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(block ~ LV) +
  theme(axis.text.x = element_text(angle = 90),
        line = element_blank()) +
  ggtitle("Crossloadings")


### education
getwd()
setwd('C:/Users/steve/Documents/gerson')
education = read.csv('education.csv')

cor(education[, 1:4])

# rows of path matrix
Support = c(0, 0, 0, 0, 0, 0)
Advising = c(0, 0, 0, 0, 0, 0)
Tutoring = c(0, 0, 0, 0, 0, 0)
Value = c(1, 1, 1, 0, 0, 0)
Satisfaction = c(1, 1, 1, 1, 0, 0)
Loyalty = c(0, 0, 0, 0, 1, 0)
# matrix (by row binding)
edu_path = rbind(Support, Advising, Tutoring, Value, Satisfaction, Loyalty)

#add col names
colnames(edu_path) = rownames(edu_path)

innerplot(edu_path, box.size = 0.1)

# outer model
edu_blocks = list(1:4, 5:8, 9:12, 13:16, 17:19, 20:23)
# modes (reflective blocks)
edu_modes = rep("A", 6)

# apply plspm
edu_pls1 = plspm(education, edu_path, edu_blocks, modes = edu_modes)
summary(edu_pls1)

edu_pls1$unidim

plot(edu_pls1, what = "loadings")

edu_blocks2 = list(c(1, 27, 3, 4), 5:8, 9:12, 13:16, 17:19, c(20, 21, 28,23))

edu_pls2 = plspm(education, edu_path, edu_blocks2, modes = edu_modes)

edu_pls2$outer_model
