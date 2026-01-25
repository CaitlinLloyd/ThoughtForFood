#recode age
dat <- subset(dat,!is.na(dat$degree))
dat <- dat %>%  dplyr::mutate(age2 = as.numeric(as.factor(dat$age))) 
dat$age2 <- dat$age2-1
Obj_info <-read.csv("nutrient_info_updated.csv")
Obj_info$stimulus <- gsub(".jpg", "", Obj_info$Food_image)
t <- merge(dat, Obj_info, by="stimulus")


# code whether a SC trial at all - also count out of trial with neutral response (note distinction between sc and sc_fac 
# fac is for the healthiness/tastiness factors, not the ratings)
ref <- subset(t, is.na(t$choice.rating))
ref <- ref %>% dplyr::select("subjectId", "stimulus",starts_with("rating"),"Tastiness","Healthiness","Umami", ends_with("cal"),ends_with("_g"))
colnames(ref)[2:31] <- paste0("neut.",colnames(ref)[2:31])

dat <- merge(dat, ref,by="subjectId",all.x=TRUE)

dat$sc_bin<-NULL
dat$sc_bin[(dat$rating.tasty > dat$neut.rating.tasty & dat$rating.healthy < dat$neut.rating.healthy) | (dat$rating.tasty < dat$neut.rating.tasty & dat$rating.healthy >  dat$neut.rating.healthy)] <- 1
dat$sc_bin[(dat$rating.tasty > dat$neut.rating.tasty & dat$rating.healthy > dat$neut.rating.healthy) | (dat$rating.tasty < dat$neut.rating.tasty & dat$rating.healthy < dat$neut.rating.healthy)] <- 0


# use of self-control (1 vs 0) - NOT expanded bins 
dat$sc<-NULL
dat$sc[(dat$rating.tasty > dat$neut.rating.tasty & dat$rating.healthy < dat$neut.rating.healthy & dat$choice.rating <5)  | (dat$rating.tasty < dat$neut.rating.tasty & dat$rating.healthy > dat$neut.rating.healthy & dat$choice.rating >5)] <- 1
dat$sc[(dat$rating.tasty > dat$neut.rating.tasty & dat$rating.healthy < dat$neut.rating.healthy & dat$choice.rating >5)  | (dat$rating.tasty < dat$neut.rating.tasty & dat$rating.healthy > dat$neut.rating.healthy & dat$choice.rating <5)] <- 0
dat$sc[dat$choice.rating == 5] <- NA
dat$sc[dat$sc_bin == 0] <- NA

dat$sc_bin_fac<-NULL
dat$sc_bin_fac[(dat$Tastiness > dat$neut.Tastiness & dat$Healthiness < dat$neut.Healthiness) | (dat$Tastiness < dat$neut.Tastiness & dat$Healthiness >  dat$neut.Healthiness)] <- 1
dat$sc_bin_fac[(dat$Tastiness > dat$neut.Tastiness & dat$Healthiness > dat$neut.Healthiness) | (dat$Tastiness < dat$neut.Tastiness & dat$Healthiness < dat$neut.Healthiness)] <- 0

# use of self-control (1 vs 0) - NOT expanded bins
dat$sc_fac<-NULL
dat$sc_fac[(dat$Tastiness > dat$neut.Tastiness & dat$Healthiness < dat$neut.Healthiness & dat$choice.rating <5)  | (dat$Tastiness < dat$neut.Tastiness & dat$Healthiness > dat$neut.Healthiness & dat$choice.rating >5)] <- 1
dat$sc_fac[(dat$Tastiness > dat$neut.Tastiness &  dat$Healthiness < dat$neut.Healthiness & dat$choice.rating >5)  | (dat$Tastiness < dat$neut.Tastiness & dat$Healthiness > dat$neut.Healthiness & dat$choice.rating <5)] <- 0
dat$sc_fac[dat$choice.rating == 5] <- NA
dat$sc_fac[dat$sc_bin_fac == 0] <- NA



# use of self-control (1 vs 0) - NOT expanded bins

# this calculates n self-control opps and usage per participant
sc <- dat %>% subset(dat$sc_bin==1) %>% dplyr::select(subjectId, sc_bin, sc) %>% group_by(subjectId,sc) %>% tally() %>% spread(sc, n) %>% dplyr::select(-4)
sc$sc.rat <- (rowSums(sc[,c("0", "1")], na.rm=TRUE))
sc$sc.rat <- sc$`1`/sc$sc.rat
sc <- as.data.frame(sc[,c(1,4)])
colnames(sc) <- c("subjectId", "sc.rat")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)



sc <- dat %>% subset(dat$sc_bin_fac==1) %>% dplyr::select(subjectId, sc_bin_fac, sc_fac) %>% group_by(subjectId,sc_fac) %>% tally () %>% spread(sc_fac, n) %>% dplyr::select(-4)
sc$sc.rat <- (rowSums(sc[,c("0", "1")], na.rm=TRUE))
sc$sc.rat <- sc$`1`/sc$sc.rat
sc <- as.data.frame(sc[,c(1,4)])
colnames(sc) <- c("subjectId", "sc.rat.factor")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)


sc <- dat %>% subset(dat$sc_bin_fac==1) %>% dplyr::select(subjectId, sc_bin_fac, sc_fac) %>% group_by(subjectId,sc_fac) %>% tally () %>% spread(sc_fac, n) %>% dplyr::select(-4)
sc$sc.rat <- (rowSums(sc[,c("0", "1")]))

sc <- as.data.frame(sc[,c(1,3)])
colnames(sc) <- c("subjectId", "sc.factor.success.n")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)

# types of self-control (1 vs 0) - NOT expanded bins
dat$sc_fac_approach_health_trial <- 0
dat$sc_fac_approach_health_trial[(dat$Tastiness < dat$neut.Tastiness & dat$Healthiness > dat$neut.Healthiness)] <- 1

dat$sc_fac_avoid_taste_trial <- 0
dat$sc_fac_avoid_taste_trial[(dat$Tastiness > dat$neut.Tastiness & dat$Healthiness < dat$neut.Healthiness)] <- 1

dat$sc_fac_approach_health <- NA
dat$sc_fac_approach_health[(dat$sc_fac_approach_health_trial==1 & dat$choice.rating > 5)] <- 1
dat$sc_fac_approach_health[(dat$sc_fac_approach_health_trial==1 & dat$choice.rating < 5)] <- 0

dat$sc_fac_avoid_taste <- NA
dat$sc_fac_avoid_taste[(dat$sc_fac_avoid_taste_trial==1 & dat$choice.rating < 5)] <- 1
dat$sc_fac_avoid_taste[(dat$sc_fac_avoid_taste_trial==1 & dat$choice.rating > 5)] <- 0


sc <- dat %>% subset(dat$sc_fac_approach_health_trial==1) %>% dplyr::select(subjectId, sc_fac_approach_health_trial, sc_fac_approach_health) %>% group_by(subjectId,sc_fac_approach_health) %>% tally () %>% spread(sc_fac_approach_health, n) %>% dplyr::select(-4)
sc$`1`[is.na(sc$`1`)] <- 0
sc$sc.total <- (rowSums(sc[,c("0", "1")]))
sc$sc.rat <- sc$`1`/sc$sc.total
sc <- as.data.frame(sc[,c(1,4,5)])
colnames(sc) <- c("subjectId", "sc.approach.n","sc.approach.prop")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)

sc <- dat %>% subset(dat$sc_fac_avoid_taste_trial==1) %>% dplyr::select(subjectId, sc_fac_avoid_taste_trial, sc_fac_avoid_taste) %>% group_by(subjectId,sc_fac_avoid_taste) %>% tally () %>% spread(sc_fac_avoid_taste, n) %>% dplyr::select(-4)
sc$`1`[is.na(sc$`1`)] <- 0
sc$sc.total <- (rowSums(sc[,c("0", "1")]))
sc$sc.rat <- sc$`1`/sc$sc.total
sc <- as.data.frame(sc[,c(1,4,5)])
colnames(sc) <- c("subjectId", "sc.avoid.n","sc.avoid.prop")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)

dat$choice.rt <- ifelse(dat$choice.rt < 200,NA, dat$choice.rt)
dat$choice.rt <- ifelse(dat$choice.rt > 5000,NA, dat$choice.rt)

dat$choice.rt.unscale <-  dat$choice.rt
dat <- dat %>% group_by(subjectId) %>%
  mutate(choice.rt = scale(choice.rt)) 
dat$choice.rt <- (dat$choice.rt)[,1]

sc <- dat %>% subset(dat$sc_bin==1) %>% dplyr::select(subjectId, sc_bin) %>% group_by(subjectId,sc_bin) %>% tally() 
sc <- sc[,c(1,3)]
colnames(sc) <- c("subjectId", "n.sctrial")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)

sc <- dat %>% subset(dat$sc_bin_fac==1) %>% dplyr::select(subjectId, sc_bin_fac) %>% group_by(subjectId,sc_bin_fac) %>% tally() 
sc <- sc[,c(1,3)]
colnames(sc) <- c("subjectId", "n.sctrial_fac")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)

sc <- dat %>% group_by(subjectId,sc_fac) %>% summarize_at("choice.rt.unscale",c(mean),na.rm=TRUE) %>% pivot_wider(1,names_from=2,values_from=3)
colnames(sc) <- c("subjectId", "sc.fail.rt","no.sc.rt","sc.success.rt")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)

sc <- dat %>% group_by(subjectId) %>% summarize_at("choice.rt.unscale",c(mean),na.rm=TRUE) 
colnames(sc) <- c("subjectId", "avg.rt")
dat <- merge(dat,sc, by.x="subjectId", by.y="subjectId", all.x = TRUE)


# change the demo codings
dat$ethnicity <- ifelse(dat$ethnicity=="White","White","Non-white")

dat$wt.cat <- cut(as.numeric(dat$bmi), breaks=c(0,18.5,25,30,50,80), labels=c("underweight", "health", "overweight", "obese", "super obese"), order=TRUE)
dat$rest.cat <- cut(as.numeric(dat$s.cognitive.restraint), breaks=c(0,33.3,44.4,61.11,100), labels=c("low restraint", "moderate restraint", "high restraint", "very high restraint"), order=TRUE)
dat$age.cat <- cut(as.numeric(dat$age2), breaks=c(0,2,5,8,11), labels=c("under 30", "30 to 45", "45 to 60", "60+"), order=TRUE)
dat$social.class <- dplyr::recode(dat$social.class, "Below the poverty level"=0, "Lower middle class"=1, "Middle class"=2, "Upper middle class"=3, "Upper class"=4)
dat$income <- dplyr::recode(dat$income, "Less than $20,000"=0, "$20,000 to $34,999"=1, "$35,000 to $49,999"=2, "$50,000 to $74,999"=3, "$75,000 to $99,999"=4, "Over $100,000"=5)

dat$income <- as.numeric(dat$income)
dat$social.class <- as.numeric(dat$social.class)

dat$education <- NA
dat$education <- ifelse(dat$degree == "Some college" | dat$degree == "High School", 0,dat$education)
dat$education <- ifelse(dat$degree == "Associates", 1,dat$education)
dat$education <- ifelse(dat$degree == "Bachelors", 2,dat$education)
dat$education <- ifelse(dat$degree == "Masters"| dat$degree == "Professional" | dat$degree == "Doctorate", 3,dat$education)
dat$education <- as.numeric(dat$education, order=TRUE)

dat$region <- as.character(dat$region)
dat$region <- ifelse(dat$region=="Fareast","Other",dat$region)
dat$region <- ifelse(dat$region=="Farwest","Other",dat$region)
dat$region <- ifelse(dat$region=="Hawaii/Alaska/US territories","Other",dat$region)
dat$region <- ifelse(dat$region=="Outside of the US","Other",dat$region)

dat <- distinct(dat, stimulus, subjectId, .keep_all= TRUE)
dat$X <- NULL

#
single_f <- dat[!duplicated(dat$subjectId),]

single_f$gender_num <- ifelse(single_f$gender=="Female",1,0)
for(var in c("gender_num","bmi","s.cognitive.restraint","age2","education","income")){
  name <- paste0(var,"_scale")
  single_f[[name]] <- scale(single_f[[var]],center=TRUE)[,1]
}

selection <- single_f %>% dplyr::select("subjectId","gender_num_scale","bmi_scale","s.cognitive.restraint_scale","age2_scale","education_scale","income_scale")
dat <- merge(selection,dat)

Obj_info <- read.csv("~/Dropbox/mTurk/nutrient_info_updated.csv")
Obj_info$stimulus <- gsub(".jpg", "", Obj_info$Food_image)
t <- merge(dat, Obj_info, by="stimulus")
t$yes <- ifelse(t$choice.rating > 5,1,0)
t$rt.bin <- ifelse(t$choice.rt <0,0,1)

t.1 <- t %>% dplyr::group_by(subjectId,HI_LO_fat) %>% tally(yes)
te <- pivot_wider(data = t.1,
                  id_cols = subjectId,
                  names_from = c(HI_LO_fat),
                  values_from = n)
colnames(te) <- c("subjectId", "lofat.yes", "hifat.yes")
t <- merge(t,te, by="subjectId", all.x = TRUE)




t.1 <- t %>% dplyr::group_by(subjectId,HI_LO_fat) %>% dplyr::count()
te <- pivot_wider(data = t.1,
                  id_cols = subjectId,
                  names_from = c(HI_LO_fat),
                  values_from = n)

colnames(te) <- c("subjectId", "lofat", "hifat")
t <- merge(t,te, by="subjectId", all.x = TRUE)

t$prop.h.f <- t$hifat.yes/t$hifat
t$prop.l.f <- t$lofat.yes/t$lofat


t$choice_bin <- NULL

# now unhealthy choice is opposite of healthy
t$choice_bin <- ifelse(t$choice.rating > 5,1,0)
t$healthy_choice <- ifelse(t$Healthiness > t$neut.Healthiness & t$choice_bin==1 | t$Healthiness < t$neut.Healthiness & t$choice_bin==0,1,0)


t.1 <- t %>% dplyr::group_by(subjectId,healthy_choice) %>% tally()
te <- pivot_wider(data = t.1,
                  id_cols = subjectId,
                  names_from = c(healthy_choice),
                  values_from = n)
te$total <- te$`0` + te$`1` 
te$choice_healthy <- te$`1`/te$total


te <- te[,c(1,5)]
colnames(te) <- c("subjectId", "healthy.choice.prop")

t <- merge(t,te, by="subjectId", all.x = TRUE)


t.1 <- t %>% dplyr::group_by(subjectId,choice_bin) %>% tally()
te <- pivot_wider(data = t.1,
                  id_cols = subjectId,
                  names_from = c(choice_bin),
                  values_from = n)
te <- te[,c(1,3)]
te$`1`[is.na(te$`1`)] <- 0
colnames(te) <- c("subjectId", "choice.alternate")

t <- merge(t,te, by="subjectId", all.x = TRUE)

## summarize RT by trial type
t.1 <- t %>% dplyr::group_by(subjectId,sc) %>% summarize_at("choice.rt.unscale",c(mean),na.rm=TRUE)
te <- pivot_wider(data = t.1,
                  id_cols = subjectId,
                  names_from = c(sc),
                  values_from = choice.rt.unscale)

colnames(te) <- c("subjectId", "rt_no_sc","rt_no_sc_trial","rt_sc")

t <- merge(t,te, by="subjectId", all.x = TRUE)

t.1 <- t %>% dplyr::group_by(subjectId) %>% dplyr::summarize(correlation_fat_pct=cor(Fat_pctKcal,rating.fat), correlation_kcal=cor(Total.kcal,rating.calories),
                                                      correlation_pro_pct=cor(PRO_pctKcal,rating.protein),correlation_carb_pct=cor(CHO_pctKcal,rating.carbohydrates),
                                                      correlation_fat_g=cor(Fat_g,rating.fat), correlation_pro_g=cor(PRO_g,rating.protein),correlation_carb_g=cor(CHO_g,rating.carbohydrates),correlation_sugar_g=cor(sugar_g,rating.sugar),correlation_sugar_pct=cor(sugar_pct_kcal,rating.sugar),
                                                      correlation_fat_Healthiness_g=cor(Fat_g,Healthiness), correlation_fat_Healthiness=cor(Fat_pctKcal,Healthiness),correlation_kcal_Healthiness=cor(Total.kcal,Healthiness),correlation_sugar_Healthiness=cor(sugar_pct_kcal,Healthiness))

t <-merge(t,t.1, by="subjectId", all = TRUE) 

health_select <- t

## info about ref item

health_select <- health_select[,c("subjectId", "healthy.choice.prop","prop.h.f","prop.l.f", "choice.alternate","rt_no_sc","rt_no_sc_trial","rt_sc","correlation_fat_pct", "correlation_kcal",
                                  "correlation_pro_pct", "correlation_carb_pct","correlation_sugar_pct", "correlation_fat_Healthiness","correlation_sugar_Healthiness", "correlation_kcal_Healthiness","correlation_fat_g","correlation_pro_g","correlation_carb_g","correlation_fat_Healthiness_g","correlation_sugar_g")]

dedup <- health_select[!duplicated(health_select$subjectId),]
single <- merge(single_f,dedup,all.x = TRUE)

plot_template <- function(dataset,pred,out,col_out,xlabel,ylims,ybreaks,title,yname){
  ylabel=ybreaks
if(class(dataset$pred)=="numeric"){
  xbreaks <- range(round(dataset$pred,0),na.rm=TRUE)

  labs <- xbreaks
  if(pred=="age2"){xbreaks=c(1,10);labs=c("18-25","71-75")}
  p <- (ggplot(dataset,aes_string(x=pred,y=out)) + geom_smooth(method = "lm",size=3,color=col_out,fill=col_out,alpha=0.6) + geom_point(size=2,color=col_out,fill=col_out,alpha=0.6) + scale_alpha_identity()  + coord_cartesian(ylim=ylims) +xlab(xlabel)+ ylab(yname) + scale_x_continuous(breaks = xbreaks, labels = labs)) + scale_y_continuous(breaks = ybreaks, labels = ylabel)
  p <- p + theme_blank() + theme(axis.line.x=element_line(color="black",size=1),axis.line.y=element_line(color="black",size=1),axis.text = element_text(size=24),axis.title = element_text(size=24),plot.title = element_text(color=col_out, size=14,face="bold",hjust = .5))
} else {
  dataset$pred <- as.factor(dataset$pred)
  p <- ggplot(dataset,aes(1,group=pred,y=out,alpha=fct_rev(pred)),fill=col_out) + ggrain::geom_rain(fill=col_out,boxplot.args.pos= list(width=0.1,position = ggpp::position_dodgenudge(0.2)),point.args=list(alpha=0)) +   coord_cartesian(ylim=ylims,clip="on") + ylab(yname) + xlab("Gender") + scale_alpha_manual(values=c(0.7,0.2)) + scale_y_continuous(breaks = ybreaks, labels = ylabel) + scale_x_continuous(breaks = c(0.8,1.4), labels = c("",""))
  p <- p + theme_blank() + theme(axis.line.x=element_line(color="black",size=1),axis.line.y=element_line(color="black",size=1),axis.text = element_text(size=24),axis.title = element_text(size=24),plot.title = element_text(color=col_out, size=14,face="bold",hjust = .5),legend.title = element_blank(),legend.position=c(.9,0.1),legend.key.size = unit(1.5,"line"))}
return(p)
}

determine_var_type <- function(x) {
  if (is.factor(x) || is.character(x)) {
    unique_vals <- length(unique(x))
    if (unique_vals == 2) return("binary")
    else if (unique_vals <= 10 && all(sort(unique(x)) == 1:unique_vals)) return("ordinal")
    else return("categorical")
  } else if (is.numeric(x)) {
    unique_vals <- length(unique(x))
    if (unique_vals == 2) return("binary")
    else if (unique_vals <= 10) return("ordinal")
    else return("continuous")
  }
  return("unknown")
}

get_association <- function(x, y) {
  complete_cases <- complete.cases(x, y)
  x_clean <- x[complete_cases]
  y_clean <- y[complete_cases]
  
  type_x <- determine_var_type(x_clean)
  type_y <- determine_var_type(y_clean)
  
  if (type_x == "continuous" && type_y == "continuous") {
    result <- cor.test(x_clean, y_clean, method = "pearson",na.action="na.omit")
    return(list(value = result$estimate, pval = result$p.value, test = "Pearson"))
  }
    else if (type_x == "ordinal" && type_y == "continuous" | type_y == "ordinal" && type_x == "continuous")  {
      result <- cor.test(x_clean, y_clean, method = "spearman",na.action="na.omit")
      return(list(value = result$estimate, pval = result$p.value, test = "Spearman"))
  } else if (type_x == "ordinal" && type_y == "ordinal") {
    result <- cor.test(x_clean, y_clean, method = "spearman",na.action="na.omit")
    return(list(value = result$estimate, pval = result$p.value, test = "Spearman"))
  } else if ((type_x == "binary" && type_y == "ordinal") || (type_x == "ordinal" && type_y == "binary")) {
    if (type_x == "binary") {
      groups <- as.factor(x_clean)
      values <- y_clean
    } else {
      groups <- as.factor(y_clean)
      values <- x_clean
    }
    result <- wilcox.test(values ~ groups)
    return(list(value = result$statistic, pval = result$p.value, test = "Wilcoxon"))
  } else if ((type_x == "binary" && type_y == "continuous") || (type_x == "continuous" && type_y == "binary")) {
    if (type_x == "binary") {
      groups <- as.factor(x_clean)
      values <- y_clean
    } else {
      groups <- as.factor(y_clean)
      values <- x_clean
    }
    result <- t.test(values ~ groups)
    return(list(value = result$statistic, pval = result$p.value, test = "t-test"))
  }
  return(list(value = NA, pval = NA, test = "None"))
}

create_correlation_table <- function(data, vars = NULL) {
  if (is.null(vars)) vars <- names(data)
  n <- length(vars)
  
  cor_matrix <- matrix(NA, n, n, dimnames = list(vars, vars))
  test_matrix <- matrix("", n, n, dimnames = list(vars, vars))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        cor_matrix[i, j] <- 1.0
        test_matrix[i, j] <- "-"
      } else {
        result <- get_association(data[[vars[i]]], data[[vars[j]]])
        cor_matrix[i, j] <- round(result$value, 3)
        test_matrix[i, j] <- result$test
      }
    }
  }
  
  # Create final display matrix
  display_matrix <- matrix("", n, n, dimnames = list(vars, vars))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        display_matrix[i, j] <- "-"
      } else if (i < j) {
        # Upper triangle: test names
        display_matrix[i, j] <- test_matrix[i, j]
      } else {
        # Lower triangle: correlation values with p-values
        pval <- get_association(data[[vars[i]]], data[[vars[j]]])$pval
        display_matrix[i, j] <- paste0(cor_matrix[i, j], " (p=", round(pval, 3), ")")
      }
    }
  }
  
  return(display_matrix)
}
