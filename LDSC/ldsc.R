library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(xtable)
setwd("~/Dropbox/McGill/Julien & Sahir/projects/MDD_UKBB/LDSC")

# Read results for Probable MDD without irritability
res_case1 <- fread("saige2_case1.log", skip = 1230, nrows = 46, header = TRUE) %>%
  select(p1, p2, rg, se, z, p) %>%
  mutate(p1 = "Probable MDD without irritability") %>%
  mutate(p2 = ifelse(str_detect(p2, "adhd"), "ADHD",
              ifelse(str_detect(p2, "case2"), "Probable MDD with irritability",
              ifelse(str_detect(p2, "antisocial"), "Antisocial behaviour",
              ifelse(str_detect(p2, "anxiety"), "Anxiety",
              ifelse(str_detect(p2, "asd"), "ASD",
              ifelse(str_detect(p2, "bip"), "Bipolar disorder",
              ifelse(str_detect(p2, "insomnia"), "Insomnia",
              ifelse(str_detect(p2, "panic"), "Panic",
              ifelse(str_detect(p2, "scz"), "Schyzophrenia",
              ifelse(str_detect(p2, "alzheimer"), "Alzheimer",
              ifelse(str_detect(p2, "intelligence"), "Intelligence",
              ifelse(str_detect(p2, "canabis"), "Canabis",
              ifelse(str_detect(p2, "cigarettes"), "Cigarettes",
              ifelse(str_detect(p2, "drinks"), "Drinks",
              ifelse(str_detect(p2, "smoking"), "Smoking initiation",
              ifelse(str_detect(p2, "education"), "Education",
              ifelse(str_detect(p2, "maltreatment"), "Childhood maltreatment",
              ifelse(str_detect(p2, "suicide_attempt"), "Suicide attempt",
              ifelse(str_detect(p2, "suicide"), "Suicide death", 
              ifelse(str_detect(p2, "ocd"), "OCD", 
              ifelse(str_detect(p2, "ptsd"), "PTSD", 
              ifelse(str_detect(p2, "swb"), "Subjective well-being",
              ifelse(str_detect(p2, "executive"), "Executive functions",
              ifelse(str_detect(p2, "substance_use"), "Substance use",
              ifelse(str_detect(p2, "social_depriv"), "Social deprivation",
              ifelse(str_detect(p2, "children"), "Number of children",
              ifelse(str_detect(p2, "agefirstbirth"), "Age at first birth",
              ifelse(str_detect(p2, "income"), "Income",
              ifelse(str_detect(p2, "loneliness"), "Loneliness",
              ifelse(str_detect(p2, "social_isolation"), "Social Isolation",
              ifelse(str_detect(p2, "risk"), "Risk tolerance",
              ifelse(str_detect(p2, "externalizing"), "Externalizing behaviors",
              ifelse(str_detect(p2, "agreableness"), "Agreableness",
              ifelse(str_detect(p2, "conscientiousness"), "Conscientiousness",
              ifelse(str_detect(p2, "extraversion"), "Extraversion",
              ifelse(str_detect(p2, "neuroticism"), "Neuroticism",
              ifelse(str_detect(p2, "openneness"), "Openneness",
              ifelse(str_detect(p2, "sensitivity_env"), "Stress adversity",
              ifelse(str_detect(p2, "birth_weight"), "Birth weight",
              ifelse(str_detect(p2, "BMI"), "BMI",
              ifelse(str_detect(p2, "height"), "Height",
              ifelse(str_detect(p2, "arthritis"), "Arthritis",
              ifelse(str_detect(p2, "chd"), "Chronic Hearth Disease",
              ifelse(str_detect(p2, "diabetes"), "Diabetes",
              ifelse(str_detect(p2, "metabolic"), "Metabolic syndrome",
              ifelse(str_detect(p2, "ed"), "Eating disorder", p2
                     ))))))))))))))))))))))))))))))))))))))))))))))
         )
  
# Read results for Probable MDD with irritability
res_case2 <- fread("saige2_case2.log", skip = 1204, nrows = 45, header = TRUE) %>%
  select(p1, p2, rg, se, z, p) %>%
  mutate(p1 = "Probable MDD with irritability") %>%
  mutate(p2 = ifelse(str_detect(p2, "adhd"), "ADHD",
              ifelse(str_detect(p2, "case2"), "Probable MDD with irritability",
              ifelse(str_detect(p2, "antisocial"), "Antisocial behaviour",
              ifelse(str_detect(p2, "anxiety"), "Anxiety",
              ifelse(str_detect(p2, "asd"), "ASD",
              ifelse(str_detect(p2, "bip"), "Bipolar disorder",
              ifelse(str_detect(p2, "insomnia"), "Insomnia",
              ifelse(str_detect(p2, "panic"), "Panic",
              ifelse(str_detect(p2, "scz"), "Schyzophrenia",
              ifelse(str_detect(p2, "alzheimer"), "Alzheimer",
              ifelse(str_detect(p2, "intelligence"), "Intelligence",
              ifelse(str_detect(p2, "canabis"), "Canabis",
              ifelse(str_detect(p2, "cigarettes"), "Cigarettes",
              ifelse(str_detect(p2, "drinks"), "Drinks",
              ifelse(str_detect(p2, "smoking"), "Smoking initiation",
              ifelse(str_detect(p2, "education"), "Education",
              ifelse(str_detect(p2, "maltreatment"), "Childhood maltreatment",
              ifelse(str_detect(p2, "suicide_attempt"), "Suicide attempt",
              ifelse(str_detect(p2, "suicide"), "Suicide death", 
              ifelse(str_detect(p2, "ocd"), "OCD", 
              ifelse(str_detect(p2, "ptsd"), "PTSD", 
              ifelse(str_detect(p2, "swb"), "Subjective well-being",
              ifelse(str_detect(p2, "executive"), "Executive functions",
              ifelse(str_detect(p2, "substance_use"), "Substance use",
              ifelse(str_detect(p2, "social_depriv"), "Social deprivation",
              ifelse(str_detect(p2, "children"), "Number of children",
              ifelse(str_detect(p2, "agefirstbirth"), "Age at first birth",
              ifelse(str_detect(p2, "income"), "Income",
              ifelse(str_detect(p2, "loneliness"), "Loneliness",
              ifelse(str_detect(p2, "social_isolation"), "Social Isolation",
              ifelse(str_detect(p2, "risk"), "Risk tolerance",
              ifelse(str_detect(p2, "externalizing"), "Externalizing behaviors",
              ifelse(str_detect(p2, "agreableness"), "Agreableness",
              ifelse(str_detect(p2, "conscientiousness"), "Conscientiousness",
              ifelse(str_detect(p2, "extraversion"), "Extraversion",
              ifelse(str_detect(p2, "neuroticism"), "Neuroticism",
              ifelse(str_detect(p2, "openneness"), "Openneness",
              ifelse(str_detect(p2, "sensitivity_env"), "Stress adversity",
              ifelse(str_detect(p2, "birth_weight"), "Birth weight",
              ifelse(str_detect(p2, "BMI"), "BMI",
              ifelse(str_detect(p2, "height"), "Height",
              ifelse(str_detect(p2, "arthritis"), "Arthritis",
              ifelse(str_detect(p2, "chd"), "Chronic Hearth Disease",
              ifelse(str_detect(p2, "diabetes"), "Diabetes",
              ifelse(str_detect(p2, "metabolic"), "Metabolic syndrome",
              ifelse(str_detect(p2, "ed"), "Eating disorder", p2
                     ))))))))))))))))))))))))))))))))))))))))))))))
         )

# Combine results for figures
res_df <- rbind(res_case1, res_case2) %>%
  filter(p2 != "Agreableness") %>% 
  mutate(lower = ifelse(rg - 1.96*se < -1, -1, rg - 1.96*se), upper = ifelse(rg + 1.96*se > 1, 1, rg + 1.96*se))

ggplot(res_df, aes(x=rg, y=p2, xmin=lower, xmax=upper, col = p1)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbarh(height=.2, position = position_dodge(0.5)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_x_continuous(breaks = seq(-1, 1, 0.2))+
  scale_y_discrete(limits=rev)+
  labs(y="", x=expression(r[g]))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())

# Table with test of r1 = r2
tab_df <- merge(
            rename(select(res_case1, -p1), r1 = rg, se1 = se, z1 = z, trait = p2, p1 = p),
            rename(select(res_case2, -p1), r2 = rg, se2 = se, z2 = z, trait = p2, p2 = p),
            by = "trait"
          ) %>%
          select(trait, r1, se1, p1, r2, se2, p2) %>%
          mutate(Z = (r2-r1)/sqrt(se1^2+se2^2)) %>%
          mutate(p.Z = 2*pnorm(q=abs(Z), lower.tail=FALSE)) %>%
          mutate(r1 = paste0(round(r1, 2), " (", round(r1-1.96*se1, 2), ", ", round(r1+1.96*se1, 2), ")")) %>%
          mutate(r2 = paste0(round(r2, 2), " (", round(r2-1.96*se2, 2), ", ", round(r2+1.96*se2, 2), ")")) %>%
          select(trait, r1, p1, r2, p2, p.Z)

# Print to latex
print(
  xtable(tab_df,
         digits = -2
  ),
  include.rownames = FALSE
)
