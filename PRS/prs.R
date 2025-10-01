library(pROC)
library(dplyr)
library(DescTools)

#---------------------------------------------------------------
# Function to calculate Nagelkerke pseudo R2 between two models
#---------------------------------------------------------------
PseudoR2 <- function(res_model, full_model, method="Nagelkerke"){
  stopifnot(length(res_model$y) == length(full_model$y))
  N <- length(res_model$y)
  as.numeric((1-exp(2/N * (logLik(res_model) - logLik(full_model))))/(1-exp(2/N * logLik(res_model))))  
}

# PseudoR2(model_null, model_prs)

#---------------------------------------------------
# Results for PRS
#---------------------------------------------------
for (case in c("case1","case2")){
  
  # Initialize empty list
  auc_all <- array(NA, dim = c(5, 4, 50))

  for (j in 1:50){
    dat  <- read.table(file = paste0("Results/", case, "_set", j, "_prscs.profile"), header = TRUE) %>%
      mutate(score = scale(SCORESUM, center = FALSE)) %>%
      group_by(set) %>%
      mutate(decile = ntile(score, 10)) %>%
      ungroup() %>%
      mutate(score_strata = ifelse(decile==10, 1, ifelse(decile==1, 0, NA)))
    
    # Null model
    model_null <- glm(eval(parse(text=case)) ~ sex+age+batch,
                      family = "binomial", 
                      data = filter(dat, set == "valid"))
    
    # Null model + 10 PCs                                  
    model_null_pcs <- glm(eval(parse(text=case))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
                          family = "binomial", 
                          data = filter(dat, set == "valid"))
    
    # Null model + 10 PCs + PRS-CS
    model_prscs <- glm(eval(parse(text=case))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+score, 
                     family = "binomial", 
                     data = filter(dat, set == "valid"))
    
    model_prscs_decile <- glm(eval(parse(text=case))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+decile, 
                       family = "binomial", 
                       data = filter(dat, set == "valid"))
    
    model_prscs_strata <- glm(eval(parse(text=case))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+score_strata, 
                            family = "binomial", 
                            data = filter(dat, set == "valid"))
    
    for (set in c("valid","test")){
      
      # Null model
      pred_null <- predict(model_null, 
                           newdata = filter(dat, set == !!set),
                           type = "response")
      
      # Null model + 10 PCs                                  
      pred_null_pcs <- predict(model_null_pcs, 
                               newdata = filter(dat, set == !!set),
                               type = "response")
      
      # Null model + 10 PCs + PRS-CS
      pred_prscs <- predict(model_prscs, 
                          newdata = filter(dat, set == !!set),
                          type = "response")
      
      pred_prscs_decile <- predict(model_prscs_decile, 
                            newdata = filter(dat, set == !!set),
                            type = "response")
      
      pred_prscs_strata <- predict(model_prscs_strata, 
                                 newdata = filter(dat, set == !!set),
                                 type = "response")
      
      if (set == "valid"){
        auc_ <- cbind(sapply(list(pred_null, pred_null_pcs,  pred_prscs, pred_prscs_decile, pred_prscs_strata), 
               function(x) auc(unlist(dat[dat$set == set, case]), x)
               ))
      } else{
        auc_ <- cbind(auc_, sapply(list(pred_null, pred_null_pcs, pred_prscs, pred_prscs_decile, pred_prscs_strata), 
                             function(x) auc(unlist(dat[dat$set == set, case]), x)
        ))
        auc_ <- cbind(auc_,
                      cbind(DescTools::PseudoR2(model_null, "Nagelkerke"),
                            DescTools::PseudoR2(model_null_pcs, "Nagelkerke"),
                            DescTools::PseudoR2(model_prscs, "Nagelkerke"),
                            DescTools::PseudoR2(model_prscs_decile, "Nagelkerke"),
                            DescTools::PseudoR2(model_prscs_strata, "Nagelkerke"))
                      %>% t()
        )
        auc_ <- cbind(auc_,
                      cbind(NA,
                            NA,
                            summary(model_prscs)$coeff["score", "Estimate"],
                            summary(model_prscs_decile)$coeff["decile", "Estimate"],
                            summary(model_prscs_strata)$coeff["score_strata", "Estimate"]
                            )
                      %>% t()
        )
        
      }
    }
    
    auc_[, 4] = exp(auc_[, 4])
    auc_all[, , j] <- auc_
  }
  
  auc_mean <- round(apply(auc_all, c(1,2), mean, na.rm = TRUE), 3)
  lcl <- apply(auc_all, c(1,2), quantile,  na.rm = TRUE, .025)
  ucl <- apply(auc_all, c(1,2), quantile,  na.rm = TRUE, .975)
  for (i in 1:length(auc_mean)){
    auc_[i] <- paste0(auc_mean[i], " (", round(lcl[i], 4), ", ", round(ucl[i], 4), ")")
  }
  colnames(auc_) <- c("AUC_val", "AUC_test", "Nagelkerke R2", "OR")
  rownames(auc_) <- c("(1) Age+Sex+Batch", "(2) Model1 + 10PCs", "(3) Model2 + PRS-CS", "(4) Model2 + PRS-CS by decile", "(5) Model2 + PRS-CS low vs high")
  
  if (case == "case1"){print("PRS for probable depression without irritability"); print(auc_)}
  else if (case == "case2"){print("PRS for probable depression with irritability"); print(auc_)}
}

#---------------------------------------------------
# # Results for trans-phenotype PRS
#---------------------------------------------------
for (case in c("case1","case2")){
  
  # Define other phenotype
  case_ <- ifelse(case == "case1", "case2", "case1")
  
  # Initialize empty list
  auc_all <- array(NA, dim = c(3, 1, 50))
  
  for (j in 1:50){
    dat  <- read.table(file = paste0("Results/", case, "_set", j, "_prscs.profile"), header = TRUE) %>%
      mutate(score = scale(SCORESUM, center = FALSE)) %>%
      group_by(set) %>%
      mutate(decile = ntile(score, 10)) %>%
      ungroup() %>%
      mutate(score_strata = ifelse(decile==10, 1, ifelse(decile==1, 0, NA)))
    
    dat_ <- read.table(file = paste0("Results/", case_, "_set", j, "_prscs.profile"), header = TRUE) %>%
              mutate(score = scale(SCORESUM, center = FALSE)) %>%
              group_by(set) %>%
              mutate(decile = ntile(score, 10)) %>%
              ungroup() %>%
              mutate(score_strata = ifelse(decile==10, 1, ifelse(decile==1, 0, NA))) %>%
              filter(!(IID %in% filter(dat, set %in% c("train", "valid"))$IID)) %>%
              filter(set != "train")
  
    # Train the model on the validation set
    model_prs <- glm(eval(parse(text=case))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+score, 
                     family = "binomial", 
                     data = filter(dat, set == "valid"))
    
    model_prs_decile <- glm(eval(parse(text=case))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+decile, 
                              family = "binomial", 
                              data = filter(dat, set == "valid"))
    
    model_prs_strata <- glm(eval(parse(text=case))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+score_strata, 
                            family = "binomial", 
                            data = filter(dat, set == "valid"))
      
    # Make predictions
    pred_prs <- predict(model_prs, 
                        newdata = dat_,
                        type = "response")
    
    pred_prs_decile <- predict(model_prs_decile, 
                               newdata = dat_,
                               type = "response")
    
    pred_prs_strata <- predict(model_prs_strata, 
                               newdata = dat_,
                               type = "response")
  
    
    auc_all[, , j] <- cbind(sapply(list(pred_prs, pred_prs_decile,  pred_prs_strata), 
                           function(x) auc(unlist(dat_[, case_]), x)
      ))

  }
  
  auc_mean <- round(apply(auc_all, c(1,2), mean, na.rm = TRUE), 3)
  lcl <- apply(auc_all, c(1,2), quantile,  na.rm = TRUE, .025)
  ucl <- apply(auc_all, c(1,2), quantile,  na.rm = TRUE, .975)
  for (i in 1:length(auc_mean)){
    auc_[i] <- paste0(auc_mean[i], " (", round(lcl[i], 4), ", ", round(ucl[i], 4), ")")
  }
  
  colnames(auc_) <- c("AUC")
  rownames(auc_) <- c("PRS", "PRS decile", "PRS high vs low")
  
  if (case == "case1"){print("PRS for probable depression with irritability"); print(auc_)}
  else if (case == "case2"){print("PRS for probable depression without irritability"); print(auc_)}
}

#---------------------------------------------------
# Boxplots
#---------------------------------------------------
library(ggplot2)
library(tidyr)
library(ggpubr)

for (case in 1:2){
  for (j in 1:50){
    
    dat <- read.table(file = paste0("Results/case", case, "_set", j, "_prscs.profile"), header = TRUE) %>%
      rename(score = SCORESUM) %>%
      group_by(set) %>%
      mutate(decile = ntile(score, 10)) %>%
      ungroup()
    
    # Calculate OR for each decile
    OR_decile <- sapply(1:10, function(i){
      m <- glm(
        eval(parse(text=paste0("case", case))) ~ sex + age + batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
          PC9 + PC10 + score_strata,
        family = "binomial",
        data = mutate(filter(dat, set == "valid"), score_strata = ifelse(decile ==
                                                      !!i, 1, ifelse(decile == 1, 0, NA)))
      )
      exp(c(coef(m)["score_strata"]
            #, confint(m, parm = "score_strata")
          ))
    }
    )
    
    # Train the predictive model using the validation set
    model_prs <- glm(eval(parse(text=paste0("case", case)))~sex+age+batch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+score, 
                     family = "binomial", 
                     data = filter(dat, set == "valid"))
    
    # Make predictions using the test set
    pred_prs <- predict(model_prs, 
                        newdata = filter(dat, set == "test"),
                        type = "response")
    
    if(case==1 && j==1){
      dat_ <- data.frame(dat[dat$set == "test", paste0("case", case)], 
                         pred_prs, 
                         pheno = "Probable depression without irritability",
                         j = j
      ) %>% rename(case = paste0("case", case))
      
      OR_decile_dat <- cbind(Decile = 1:10, OR_decile, 
                             case = "Probable depression without irritability"
      )
      OR_decile_dat[1, 2] <- 1
      colnames(OR_decile_dat)[2] <- "OR"
      
    } else if (case==1){
      dat_ <- rbind(dat_,
                    data.frame(dat[dat$set == "test", paste0("case", case)], 
                               pred_prs, 
                               pheno = "Probable depression without irritability",
                               j = j
                    ) %>% rename(case = paste0("case", case))
      )
      
      OR_decile_dat <- rbind(cbind(Decile = 1:10, OR_decile, 
                                   case = "Probable depression without irritability"
                                  )
                             , OR_decile_dat
      )
      OR_decile_dat[1, 2] <- 1
      colnames(OR_decile_dat)[2] <- "OR"
      
    } else if (case==2){
      dat_ <- rbind(dat_,
                    data.frame(dat[dat$set == "test", paste0("case", case)], 
                               pred_prs, 
                               pheno = "Probable depression with irritability",
                               j = j
                    ) %>% rename(case = paste0("case", case))
      )
      
      OR_decile_dat <- rbind(cbind(Decile = 1:10, OR_decile, 
                                   case = "Probable depression with irritability"
                                  )
                            , OR_decile_dat
      )
      OR_decile_dat[1, 2] <- 1
      colnames(OR_decile_dat)[2] <- "OR"
    }
    
  }
  
  p <- paste0("p", case)
  
  assign(p, 
      data.frame(dat_) %>%
        mutate(case = factor(case),
               pred = 100*ntile(pred_prs, length(pred_prs))/length(pred_prs)) %>%
      ggplot(aes(y=pred, group=case, x =case, fill = case)) +
      geom_boxplot() +
      theme_bw()+
      scale_fill_discrete(name = "Group", labels = c("controls", "cases"))+
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      #theme(panel.border= element_blank())+
      theme(axis.line.x = element_line(color="black", size = 0.5),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.y = element_line(color="black", size = 0.5),
            plot.title = element_text(hjust = 0.5))
    )
  
  if (case == 2){
    p3 <- dat_ %>%
              group_by(pheno, j) %>%
              mutate(decile = ntile(pred_prs, 10)) %>%
              group_by(pheno, decile, j) %>%
              summarise(prev = 100*mean(as.numeric(case))) %>%
              group_by(pheno, decile) %>%
              summarise(LCL = quantile(prev, .025), UCL = quantile(prev, .975), prev = mean(prev)) %>%
              ggplot(aes(y=prev, group=decile, x = decile, colour = pheno)) +
              geom_point() +
              geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2, position = position_dodge(0.2)) +
              geom_line(aes(group=pheno)) +
              labs(y = "Prevalence (%) (95% CI)") +
              scale_x_continuous(name = "PRS decile", breaks = 1:10, labels = seq(1, 10, 1)) +
              theme_bw()+
              theme(plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())+
              #theme(panel.border= element_blank())+
              theme(axis.line.x = element_line(color="black", size = 0.5),
                    axis.line.y = element_line(color="black", size = 0.5),
                    plot.title = element_text(hjust = 0.5)) +
              theme(legend.title=element_blank())
  
    p4 <- data.frame(OR_decile_dat) %>%
             group_by(case, Decile) %>%
             summarise(LCL = quantile(as.numeric(OR), .025), UCL = quantile(as.numeric(OR), .975), OR = mean(as.numeric(OR))) %>%
             ggplot(aes(y=as.numeric(OR), x=as.numeric(Decile), colour=case)) +
             geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.2, position = position_dodge(0.2)) +
             geom_line(position = position_dodge(0.1)) +
             geom_point(position = position_dodge(0.1)) +
             labs(y = "OR (95% CI)") +
             scale_x_continuous(name = "PRS decile", breaks = 1:10, labels = seq(1, 10, 1)) +
             theme_bw()+
             theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
             #theme(panel.border= element_blank())+
             theme(axis.line.x = element_line(color="black", size = 0.5),
                   axis.line.y = element_line(color="black", size = 0.5),
                   plot.title = element_text(hjust = 0.5)) +
             theme(legend.title=element_blank())
    
  }
  
  if(case == 1) {
    assign(p, eval(parse(text=p)) + labs(y = "Prediction percentile", title = "PRS for probable depression without irritability"))
  } else {
    assign(p, eval(parse(text=p)) + labs(y = "Prediction percentile", title = "PRS for probable depression with irritability"))
  }
  
}

ggarrange(p1, p2, common.legend = TRUE, legend="bottom")
ggarrange(p3, p4, common.legend = TRUE, legend="bottom")

#---------------------------------------------------
# Compare betas between PLINK and SAIGE
#---------------------------------------------------

betas <- readRDS("betas.rds")
mean(betas$A1 == betas$Allele1)
mean(betas$A1 == betas$Allele2)

library(ggplot2)
library(dplyr)

betas %>%
  mutate(LEGEND = ifelse(A1 == Allele1, 1, 0)) %>%
  ggplot(aes(x=BETA_SAIGE, y=BETA_PLINK, col = as.factor(LEGEND))) +
  scale_color_manual(labels = c("A1 = Allele2", "A1 = Allele1"), values = c("blue", "red")) +
  geom_point() +
  ylim(-1, 1) +
  labs(x = "BETA_SAIGE") +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())

betas %>%
  mutate(BETA_SAIGE = ifelse(AF_Allele2 > 0.5, -BETA_SAIGE, BETA_SAIGE), 
         LEGEND = ifelse(A1 == Allele1, 1, 0)) %>%
  ggplot(aes(x=BETA_SAIGE, y=BETA_PLINK, col = as.factor(LEGEND))) +
  scale_color_manual(labels = c("A1 = Allele2", "A1 = Allele1"), values = c("blue", "red")) +
  geom_point() +
  ylim(-1, 1) +
  labs(x = "BETA_SAIGE") +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())

