library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(miamiplot)

setwd("/Users/julienst-pierre/Library/CloudStorage/Dropbox/McGill/Julien & Sahir/projects/MDD_UKBB/LDSC")

#------------------------------
# 1. Miami plot
#------------------------------
# Read GWAS results file
res_case1 <- fread("saige2_case1_full.gz", header = TRUE) %>% mutate(study = "case1")
res_case2 <- fread("saige2_case2_full.gz", header = TRUE) %>% mutate(study = "case2")

# Combine into same file
res_gwas <- bind_rows(res_case1, res_case2) %>%
  select(MarkerID, CHR, POS, BETA, SE, p.value, study) %>%
  mutate(Chi2 = ifelse(study == "case1", (BETA/SE)^2/1.0373, (BETA/SE)^2/1.027)) %>%
  mutate(p.value.adj = pchisq(Chi2, df=1, lower.tail = FALSE))

# Miami plot
my_upper_colors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]
my_lower_colors <- RColorBrewer::brewer.pal(4, "Paired")[3:4]

ggmiami(data = filter(res_gwas, p.value.adj <= 5e-2), 
        split_by = "study", 
        split_at = "case1", 
        p = "p.value.adj",
        chr = "CHR" ,
        pos = "POS" ,
        upper_ylab = "Probable MDD without irritability",
        lower_ylab = "Probable MDD with irritability",
        chr_colors = NULL,
        upper_chr_colors = my_upper_colors, 
        lower_chr_colors = my_lower_colors
)

#-------------------------------
# 2. Comparison of effect sizes
#-------------------------------
# Read result files from FUMA
res_case1_ind <- fread("FUMA_case1/GenomicRiskLoci.txt", header = TRUE) %>% 
                 mutate(study = "case1")

res_case2_ind <- fread("FUMA_case2/GenomicRiskLoci.txt", header = TRUE) %>% 
  mutate(study = "case2")

# Combine and merge with main file
res_gwas_ind <- bind_rows(res_case1_ind, res_case2_ind) %>%
                rename(MarkerID = rsID) %>%
                select(MarkerID, chr, pos, study) %>%
                merge(res_gwas, by = c("MarkerID", "study")) %>%
                mutate(BETA1 = ifelse(study == "case1", BETA, NA)) %>%
                mutate(BETA2 = ifelse(study == "case2", BETA, NA)) %>%
                mutate(SE1 = ifelse(study == "case1", SE, NA)) %>%
                mutate(SE2 = ifelse(study == "case2", SE, NA)) %>%
                mutate(p1 = ifelse(study == "case1", p.value, NA)) %>%
                mutate(p2 = ifelse(study == "case2", p.value, NA)) %>%
                mutate(Chi2 = ifelse(study == "case1", (BETA/SE)^2/1.0373, (BETA/SE)^2/1.027)) %>%
                mutate(p.adj1 = ifelse(study == "case1", pchisq(Chi2, df=1, lower.tail = FALSE), NA)) %>%
                mutate(p.adj2 = ifelse(study == "case2", pchisq(Chi2, df=1, lower.tail = FALSE), NA)) %>%
                mutate(study = ifelse(study == "case1", "case2", "case1")) %>%
                select(MarkerID, chr, pos, study, BETA1, BETA2, SE1, SE2, p1, p2, p.adj1, p.adj2) %>%
                merge(res_gwas, by = c("MarkerID", "study")) %>%
                mutate(BETA1 = ifelse(is.na(BETA1), BETA, BETA1)) %>%
                mutate(BETA2 = ifelse(is.na(BETA2), BETA, BETA2)) %>%
                mutate(SE1 = ifelse(is.na(SE1), SE, SE1)) %>%
                mutate(SE2 = ifelse(is.na(SE2), SE, SE2)) %>%
                mutate(p1 = ifelse(is.na(p1), p.value, p1)) %>%
                mutate(p2 = ifelse(is.na(p2), p.value, p2)) %>%
                mutate(p.adj1 = ifelse(is.na(p.adj1), pchisq((BETA1/SE1)^2/1.0373, df=1, lower.tail = FALSE), p.adj1)) %>%
                mutate(p.adj2 = ifelse(is.na(p.adj2), pchisq((BETA2/SE2)^2/1.027, df=1, lower.tail = FALSE), p.adj2)) %>%
                select(MarkerID, chr, pos, BETA1, BETA2, SE1, SE2, p1, p2, p.adj1, p.adj2) %>%
                mutate(rho = sqrt(35857 * 23613 / (303869 * 291625))) %>%
                mutate(Z = (BETA1 - BETA2) / sqrt(SE1^2*1.0373 + SE2^2*1.027 - 2 * rho * SE1 * SE2 * sqrt(1.0373 * 1.027))) %>%
                mutate(p.Z = 2*pnorm(q=abs(Z), lower.tail=FALSE)) %>%
                select(-p1,-p2) %>%
                filter(p.adj1 <= 5e-6 | p.adj2 <= 5e-6) %>%
                mutate(SE1 = SE1 * sqrt(1.0373), SE2 = SE2 * sqrt(1.027))
                
# Create plot
ggplot(res_gwas_ind, aes(x=BETA1, y = BETA2, col = -log10(p.Z))) +
  geom_point(size=1) +
  geom_abline(lty=2)+
  geom_errorbar(aes(ymin=BETA2-1.96*SE2, ymax=BETA2+1.96*SE2, col = -log10(p.Z)), linewidth=.25)+
  geom_errorbar(aes(xmin=BETA1-1.96*SE1, xmax=BETA1+1.96*SE1, col = -log10(p.Z)), linewidth=.25)+
  scale_color_gradient2(name=bquote(atop(-log[10](P-value), H[0]:beta[1]==beta[2])), low = "yellow", mid = "red", high = "darkblue", midpoint = mean(-log10(res_gwas_ind$p.Z))) + 
  scale_y_continuous(limits = c(-0.3, 0.4)) +
  scale_x_continuous(limits = c(-0.3, 0.3)) +
  labs(y=expression(MDD ~ with ~ irritability ~ effect ~ size~(beta[2])), x=expression(MDD ~ without ~ irritability ~ effect ~ size~(beta[1])))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        plot.title = element_text(hjust = 0.5))

# Print table to latex
print(
  xtable::xtable(
    arrange(res_gwas_ind, chr, pos)
    , digits = c(0, rep(1, 3), rep(3, 4), rep(-2, 2), 2, -2)
  )
  , include.rownames = FALSE
)

#-------------------------------
# 3. Table with lead SNPs
#-------------------------------
# Read result files from FUMA
res_case1_ind <- fread("FUMA_case1/GenomicRiskLoci.txt", header = TRUE) %>%
  merge(select(fread("FUMA_case1/snps.txt", header = TRUE), uniqID, non_effect_allele,effect_allele, MAF, beta, se, nearestGene, dist, func), by = c("uniqID")) %>%
  mutate(OR = exp(beta), A1 = non_effect_allele, A2 = effect_allele) %>%
  mutate(Chi2 =(beta/se)^2/1.0373) %>%
  mutate(p.adj = pchisq(Chi2, df=1, lower.tail = FALSE)) %>%
  select(GenomicLocus, chr, pos, rsID, nearestGene, dist, func, A1, A2, MAF, OR, p.adj) %>%
  filter(p.adj <= 5e-6) %>%
  arrange(chr, pos) %>%
  mutate(GenomicLocus = row_number()) %>%
  print()

print(
  xtable::xtable(
    res_case1_ind
    , digits = c(0, rep(2, 9), 3, 3, -2)
  )
  , include.rownames = FALSE
)


res_case2_ind <- fread("FUMA_case2/GenomicRiskLoci.txt", header = TRUE) %>%
  merge(select(fread("FUMA_case2/snps.txt", header = TRUE), uniqID, non_effect_allele,effect_allele, MAF, beta, se, nearestGene, dist, func), by = c("uniqID")) %>%
  mutate(OR = exp(beta), A1 = non_effect_allele, A2 = effect_allele) %>%
  mutate(Chi2 = (beta/se)^2/1.027) %>%
  mutate(p.adj = pchisq(Chi2, df=1, lower.tail = FALSE)) %>%
  select(GenomicLocus, chr, pos, rsID, nearestGene, dist, func, A1, A2, MAF, OR, p.adj) %>%
  filter(p.adj <= 5e-6) %>%
  arrange(chr, pos) %>%
  mutate(GenomicLocus = row_number()) %>%
  print()

print(
  xtable::xtable(
    res_case2_ind
    , digits = c(0, rep(2, 9), 3, 3, -2)
  )
  , include.rownames = FALSE
)

