library(openxlsx)
library(dplyr)
library(vegan)
library(compositions)
library(ggplot2)
library(mediation)
library(miMediation)
library(tidyverse)
library(lsr)
library(Maaslin2)
library(tidyverse)
library(microeco)
library(magrittr)
library(patchwork)
library(glmnet)
library(tableone)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggsci)
library(psych)
library(ggcorrplot)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(car)



# Distribution of pesticide and their correlation; Figure 2A, Table S7 -------------------------

name_all <- colnames(data_part1)
anemia <- c(24:25)
pesticide <- c(80:108)
anemia_parameter_24 <- c(26:27)
anemia_parameter_32 <- c(28:29)
anemia_parameter <- c(anemia_parameter_24,anemia_parameter_32)

geo_mean <- function(data) {
  log_data=log(data)
  gm=round(exp(mean(log_data[is.finite(log_data)])),2)
  return(gm)
}
quar <- function(data){
  quartile=round(quantile(data,probs = c(0.025,0.25,0.5,0.75,0.975)),2)
  return(quartile)
}
range_diff <- function(data){
  low=round(range(data)[1],2)
  high=round(range(data)[2],2)
  return(paste0(low,"–",high))
}

result <- as.data.frame(matrix(NA,ncol = 9))
names(result) <- c("pesticide","mean","GM","p2.5","p25","p50","p75","p97.5","range")

for (i in 1:length(pesticide)){
  result[i,1] <- names(data_part1)[pesticide[i]]
  result[i,2] <- round(mean(data_part1[,pesticide[i]]),2)
  result[i,3] <- round(geo_mean(data_part1[,pesticide[i]]),2)
  result[i,4] <- quar(data_part1[,pesticide[i]])[1]
  result[i,5] <- quar(data_part1[,pesticide[i]])[2]
  result[i,6] <- quar(data_part1[,pesticide[i]])[3]
  result[i,7] <- quar(data_part1[,pesticide[i]])[4]
  result[i,8] <- quar(data_part1[,pesticide[i]])[5]
  result[i,9] <- range_diff(data_part1[,pesticide[i]])
}
result <- result %>% arrange(pesticide)
write.csv(result,"pesticide_distribution.csv",quote = F,row.names = F) # 手动修改，保留两位小数

# pesticide abundance 

data_part1_pesticide <- data_part1[,pesticide]
data_part1_pesticide <- log10(data_part1_pesticide) #log10转换
data_part1_pesticide_melt <- melt(data_part1_pesticide,variable.name = "pesticide",value.name = "concentration")
data_part1_pesticide_melt$pesticide <- factor(data_part1_pesticide_melt$pesticide,levels = sort(colnames(data_part1_pesticide)))

opar <- circos.par()
pdf("pesticide_abundance_circos.pdf",height = 10,width = 10)
circos.par("track.height" = 0.3,cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(data_part1_pesticide_melt$pesticide, x = data_part1_pesticide_melt$concentration,sector.width = 0.15,xlim = c(-2,2))
circos.trackHist(data_part1_pesticide_melt$pesticide, x = data_part1_pesticide_melt$concentration,
                 bg.col = "#e6e5e5",border = NA,
                 col=c(brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"))[1:length(pesticide)] #Display the color range of the connection according to the correlation size
)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index,
              facing = "clockwise", niceFacing = T, adj = c(-0.5, 0.5), cex = 0.8)
  circos.axis(h = "top", labels.cex = 0.5,minor.ticks = 1,major.at = c(-2,-1,0,1,2))
}, bg.border = NA)
dev.off()
circos.clear()

# Spearman correlation of pesticides
data_part1_pesticide <- as.matrix(data_part1[,pesticide])
data_part1_pesticide <- log(data_part1_pesticide)
pesticide_cor <- corr.test(data_part1_pesticide, method = "spearman",adjust = "BH")
pesticide_corr <- data.frame(pesticide_cor$r);diag(pesticide_corr) <-0
pesticide_p <- data.frame(pesticide_cor$p)
write.csv(pesticide_corr,"pesticide_corr.csv",quote = F,row.names = T)
write.csv(pesticide_p,"pesticide_p.csv",quote = F,row.names = T)

for (i in 1:length(pesticide)){
  for (j in 1:i){
    pesticide_p[j,i] = pesticide_p[i,j]
  }
}

pesticide_corr$rowname <- rownames(pesticide_corr);pesticide_corr_melt <- melt(pesticide_corr)
pesticide_p$rowname <- rownames(pesticide_p);pesticide_p_melt <- melt(pesticide_p)
pesticide_corr_melt <- pesticide_corr_melt[which(pesticide_p_melt$value<0.05),]
pesticide_corr_melt <- subset(pesticide_corr_melt,value!=0)

pesticide_corr_melt$rowname <- factor(pesticide_corr_melt$rowname,levels = sort(colnames(data_part1_pesticide)))
pesticide_corr_melt$variable <- factor(pesticide_corr_melt$variable,levels = sort(colnames(data_part1_pesticide)))

pdf("pesticide_chordgram.pdf",height = 10,width = 10)
circos.par(gap.after = rep(8, nrow(pesticide_corr))) # gap数和变量数保持一致
chordDiagram(pesticide_corr_melt, 
             directional = 0, # 0 represents no direction, 1 represents forward direction, -1 represents reverse direction, 2 represents two-way
             grid.col = c(brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"))[1:length(pesticide)], # ring color
             order = sort(colnames(data_part1_pesticide)),
             annotationTrack = "grid", # Draw the outer arc area and display the name and scale axis
             col = colorRamp2(c(-1, 0, 1), c('#238dce', 'white', '#d12918'), transparency = 0), #Display the color range of the connection according to the correlation size
             annotationTrackHeight = c(0.04, 0.1) #The distance of the name from the arc, and the width of the arc
)
circos.track(track.index = 1, panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
  #             facing = "clockwise", niceFacing = T, adj = c(-0.5, 0.5), cex = 0.8)
  circos.axis(h = "top", labels.cex = 0.4,minor.ticks = 4,major.at = c(0,10,20,30))
}, bg.border = NA)
dev.off()
circos.clear()

bar_plot <- data.frame(X = 1:length(seq(-1,1,by=2/500)),Y = seq(-1,1,by=2/500))
bar_plot$X <- as.factor(bar_plot$X)
ggplot(bar_plot,aes(x=X, y=1, fill=Y)) + 
  geom_tile(alpha=1) + theme_minimal() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colours=colorRampPalette(c('#238dce', 'white', '#d12918'))(500), limits=c(-1, 1), guide=FALSE)
ggsave("pesticide_class_legend.pdf",height = 1,width = 6)

# selection of gut microbiota  --------------------------------------------

microbiota_T2_part1 <- read.xlsx("data_part1.xlsx",sheet = 2,rowNames = T)
shannon_T2_part1 <- diversity(microbiota_T2_part1,index = "shannon") %>% data.frame()
# raw counts and relative abundance are allowed
colnames(shannon_T2_part1) <- "shannon"

result <- data.frame(matrix(NA,ncol = 3))
colnames(result) <- c("microbiota","prevalence","relative_abundance")
for (i in 1:ncol(microbiota_T2_part1)){
  result[i,1]=colnames(microbiota_T2_part1)[i]
  result[i,2]=sum(microbiota_T2_part1[,i]!=0)/nrow(microbiota_T2_part1)
  result[i,3]=sum(microbiota_T2_part1[,i])/sum(microbiota_T2_part1)
}
genus_select <- result[result$prevalence>0.1 & result$relative_abundance>0.00005,]$microbiota
microbiota_T2_part1_select <- microbiota_T2_part1[,genus_select]
microbiota_T2_part1_select <- microbiota_T2_part1_select/rowSums(microbiota_T2_part1_select)
write.csv(microbiota_T2_part1_select,"microbiota_T2_part1_select.csv",quote = F,row.names = T)

# distribution of gut microbiota


# correlations of bacteria; Figure S2 ------------------------------------------------

bacteria_cor <- corr.test(as.matrix(clr(microbiota_T2_part1_select)), method = "spearman",adjust = "BH")
bacteria_corr <- data.frame(bacteria_cor$r);diag(bacteria_corr) <- NA
bacteria_p <- data.frame(bacteria_cor$p);diag(bacteria_p) <- NA
write.csv(bacteria_corr,"bacteria_corr.csv",quote = F,row.names = T)
write.csv(bacteria_p,"bacteria_p.csv",quote = F,row.names = T)

for (i in 1:length(bacteria_p)){
  for (j in 1:i){
    bacteria_p[j,i] = bacteria_p[i,j]
  }
}
bacteria_p <- matrix(ifelse(bacteria_p<0.05,"*",""),nrow(bacteria_p))
bacteria_p[is.na(bacteria_p)] <- ""

bacteria_corr$rowname <- rownames(bacteria_corr)
bacteria_corr$final_name <- gsub(";g__uncultured","",rownames(bacteria_corr))
bacteria_corr$final_name <- gsub(";__","",bacteria_corr$final_name) # 确定final_name,即细菌的最终名称
bacteria_corr$final_name <- gsub(".*;","",bacteria_corr$final_name)
bacteria_corr$final_name[grep("f__Lachnospiraceae;g__uncultured",bacteria_corr$rowname)] <- "f__Lachnospiraceae;g__uncultured"
bacteria_corr$final_name[grep("f__Lachnospiraceae;__",bacteria_corr$rowname)] <- "f__Lachnospiraceae;__"
bacteria_corr$final_name[grep("f__Oscillospiraceae;g__uncultured",bacteria_corr$rowname)] <- "f__Oscillospiraceae;g__uncultured"
bacteria_corr$final_name[grep("f__Oscillospiraceae;__",bacteria_corr$rowname)] <- "f__Oscillospiraceae;__"
bacteria_corr$final_name[grep("f__Ruminococcaceae;g__uncultured",bacteria_corr$rowname)] <- "f__Ruminococcaceae;g__uncultured"
bacteria_corr$final_name[grep("f__Ruminococcaceae;__",bacteria_corr$rowname)] <- "f__Ruminococcaceae;__"

rownames(bacteria_corr) <- bacteria_corr$final_name
colnames(bacteria_corr)[1:c(dim(bacteria_corr)[2]-2)] <- bacteria_corr$final_name
bacteria_corr <- bacteria_corr[,1:c(dim(bacteria_corr)[2]-2)]

pdf("bacteria_correlation.pdf",width = 18,height = 18)
pheatmap(bacteria_corr,cluster_cols = F,cluster_rows = F,scale = "none",
         fontsize = 7,border="black",
         cellwidth = 9,cellheight = 9,
         display_numbers = bacteria_p,
         fontsize_number = 8,number_color = "black",
         color = c(colorRampPalette(colors = c("#238dce","white"))(378),colorRampPalette(colors = c("white","#d12918"))(479))
         )
dev.off()

#
bacteria_corr <- read.csv("bacteria_corr.csv",row.names = 1)
bacteria_p <- read.csv("bacteria_p.csv",row.names = 1)
bacteria_corr$rowname <- rownames(bacteria_corr);bacteria_corr_melt <- melt(bacteria_corr)
bacteria_p$rowname <- rownames(bacteria_p);bacteria_p_melt <- melt(bacteria_p)
bacteria_corr_melt <- bacteria_corr_melt[which(bacteria_p_melt$value<0.05),]
bacteria_corr_melt <- subset(bacteria_corr_melt,value!=0)
range(bacteria_corr_melt$value)

# RCS for assessing non-linearity pesticide and anemia; Figure S4, Table S4 -----------------------------------------

library(rms)
library(ggrcs)
library(nnet)
# set dummy variable
dd <- data_part1
dd[,name_all[pesticide]] <- log(dd[,name_all[pesticide]])
a <- class.ind(dd$education);colnames(a) <- paste0("education",colnames(a))
dd <- cbind(dd,a)

# Furalaxyl anemia_24 (glm sig)

dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- lrm(anemia_24 ~ rcs(Furalaxyl,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_24
           ,data = dd)
AIC(fit)
fit <- lrm(anemia_24 ~ rcs(Furalaxyl,4)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_24
           ,data = dd)
AIC(fit)
fit <- lrm(anemia_24 ~ rcs(Furalaxyl,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_24
           ,data = dd)
AIC(fit)

dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- lrm(anemia_24 ~ rcs(Furalaxyl,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_24
           ,data = dd)
dd_dist$limits$Furalaxyl[2] <- median(dd$Furalaxyl) # set as the median
fit <- update(fit)
OR<-Predict(fit,Furalaxyl,fun=exp,ref.zero = TRUE)

ggplot(data=OR)+
  geom_line(aes(Furalaxyl,yhat),linetype="solid",size=1,colour="#e3191b")+
  geom_ribbon(aes(Furalaxyl,ymin = lower, ymax = upper),alpha = 0.3,fill="#e3191b")+
  geom_hline(yintercept = 1,linetype=2,size=1)+
  labs(x="ln(Furalaxyl)",y="OR (95%CI)")+
  ggtitle(paste0("P for overall association: ",round(anova(fit)[1,3],3),
                 "\nP for non-linearity: ",round(anova(fit)[2,3],3)))+
  theme_classic()+
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        title = element_text(color="black",size=10))
ggsave("Furalaxyl_anemia_24_rcs.pdf",width = 4,height = 4)

# Atrazine anemia_32  (glm sig)
dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- lrm(anemia_32 ~ rcs(Atrazine,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- lrm(anemia_32 ~ rcs(Atrazine,4)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- lrm(anemia_32 ~ rcs(Atrazine,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)

dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- lrm(anemia_32 ~ rcs(Atrazine,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
dd_dist$limits$Atrazine[2] <- median(dd$Atrazine) # set as the median
fit <- update(fit)
OR<-Predict(fit,Atrazine,fun=exp,ref.zero = TRUE)

ggplot(data=OR)+
  geom_line(aes(Atrazine,yhat),linetype="solid",size=1,colour="#e3191b")+
  geom_ribbon(aes(Atrazine,ymin = lower, ymax = upper),alpha = 0.3,fill="#e3191b")+
  geom_hline(yintercept = 1,linetype=2,size=1)+
  labs(x="ln(Atrazine)",y="OR (95%CI)")+
  ggtitle(paste0("P for overall association: ",round(anova(fit)[1,3],3),
                 "\nP for non-linearity: ",round(anova(fit)[2,3],3)))+
  theme_classic()+ 
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        title = element_text(color="black",size=10))
ggsave("Atrazine_anemia_32_rcs.pdf",width = 4,height = 4)

# Pyrimethanil Hb_32 (glm sig)
dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(Hb_32 ~ rcs(Pyrimethanil,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(Hb_32 ~ rcs(Pyrimethanil,4)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(Hb_32 ~ rcs(Pyrimethanil,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)


dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(Hb_32 ~ rcs(Pyrimethanil,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
dd_dist$limits$Pyrimethanil[2] <- median(dd$Pyrimethanil) # set as the median
fit <- update(fit)
beta<-Predict(fit,Pyrimethanil,ref.zero = TRUE)

ggplot(data=beta)+
  geom_line(aes(Pyrimethanil,yhat),linetype="solid",size=1,colour="#6a3d99")+
  geom_ribbon(aes(Pyrimethanil,ymin = lower, ymax = upper),alpha = 0.3,fill="#6a3d99")+
  geom_hline(yintercept = 0,linetype=2,size=1)+
  labs(x="ln(Pyrimethanil)",y="beta (95%CI)")+
  ggtitle(paste0("P for overall association: ",round(anova(fit)[1,3],3),
                 "\nP for non-linearity: ",round(anova(fit)[2,3],3)))+
  theme_classic()+ 
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        title = element_text(color="black",size=10))
ggsave("Pyrimethanil_Hb_32_rcs.pdf",width = 4,height = 4)

# Clomazone Hb_32 (glm sig)
dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(Hb_32 ~ rcs(Clomazone,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(Hb_32 ~ rcs(Clomazone,4)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(Hb_32 ~ rcs(Clomazone,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)

dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(Hb_32 ~ rcs(Clomazone,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
dd_dist$limits$Clomazone[2] <- median(dd$Clomazone) # set as the median
fit <- update(fit)
beta<-Predict(fit,Clomazone,ref.zero = TRUE)

ggplot(data=beta)+
  geom_line(aes(Clomazone,yhat),linetype="solid",size=1,colour="#6a3d99")+
  geom_ribbon(aes(Clomazone,ymin = lower, ymax = upper),alpha = 0.3,fill="#6a3d99")+
  geom_hline(yintercept = 0,linetype=2,size=1)+
  labs(x="ln(Clomazone)",y="beta (95%CI)")+
  ggtitle(paste0("P for overall association: ",round(anova(fit)[1,3],3),
                 "\nP for non-linearity: ",round(anova(fit)[2,3],3)))+
  theme_classic()+ 
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        title = element_text(color="black",size=10))
ggsave("Clomazone_Hb_32_rcs.pdf",width = 4,height = 4)

# Clomazone RBC_32 (glm sig)
dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(RBC_32 ~ rcs(Clomazone,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(RBC_32 ~ rcs(Clomazone,4)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(RBC_32 ~ rcs(Clomazone,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)

dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(RBC_32 ~ rcs(Clomazone,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
dd_dist$limits$Clomazone[2] <- median(dd$Clomazone) # set as the median
fit <- update(fit)
beta<-Predict(fit,Clomazone,ref.zero = TRUE)

ggplot(data=beta)+
  geom_line(aes(Clomazone,yhat),linetype="solid",size=1,colour="#339f2c")+
  geom_ribbon(aes(Clomazone,ymin = lower, ymax = upper),alpha = 0.3,fill="#339f2c")+
  geom_hline(yintercept = 0,linetype=2,size=1)+
  labs(x="ln(Clomazone)",y="beta (95%CI)")+
  ggtitle(paste0("P for overall association: ",round(anova(fit)[1,3],3),
                 "\nP for non-linearity: ",round(anova(fit)[2,3],3)))+
  theme_classic()+ 
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        title = element_text(color="black",size=10))
ggsave("Clomazone_RBC_32_rcs.pdf",width = 4,height = 4)

# Atrazine RBC_32 (glm sig)
dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(RBC_32 ~ rcs(Atrazine,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(RBC_32 ~ rcs(Atrazine,4)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(RBC_32 ~ rcs(Atrazine,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)


dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(RBC_32 ~ rcs(Atrazine,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
dd_dist$limits$Atrazine[2] <- median(dd$Atrazine) # set as the median
fit <- update(fit)
beta<-Predict(fit,Atrazine,ref.zero = TRUE)

ggplot(data=beta)+
  geom_line(aes(Atrazine,yhat),linetype="solid",size=1,colour="#339f2c")+
  geom_ribbon(aes(Atrazine,ymin = lower, ymax = upper),alpha = 0.3,fill="#339f2c")+
  geom_hline(yintercept = 0,linetype=2,size=1)+
  labs(x="ln(Atrazine)",y="beta (95%CI)")+
  ggtitle(paste0("P for overall association: ",round(anova(fit)[1,3],3),
                 "\nP for non-linearity: ",round(anova(fit)[2,3],3)))+
  theme_classic()+ 
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        title = element_text(color="black",size=10))
ggsave("Atrazine_RBC_32_rcs.pdf",width = 4,height = 4)

# Metribuzin RBC_32 (glm sig)
dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(RBC_32 ~ rcs(Metribuzin,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(RBC_32 ~ rcs(Metribuzin,4)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)
fit <- Glm(RBC_32 ~ rcs(Metribuzin,5)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
AIC(fit)

dd_dist <- datadist(dd)
options(datadist='dd_dist')
fit <- Glm(RBC_32 ~ rcs(Metribuzin,3)+age+BMI+parity+income+education2+education3+smoking_history+folic_acid+iron+alcohol+coffee+GW_32
           ,data = dd)
dd_dist$limits$Metribuzin[2] <- median(dd$Metribuzin) # set as the median
fit <- update(fit)
beta<-Predict(fit,Metribuzin,ref.zero = TRUE)

ggplot(data=beta)+
  geom_line(aes(Metribuzin,yhat),linetype="solid",size=1,colour="#339f2c")+
  geom_ribbon(aes(Metribuzin,ymin = lower, ymax = upper),alpha = 0.3,fill="#339f2c")+
  geom_hline(yintercept = 0,linetype=2,size=1)+
  labs(x="ln(Metribuzin)",y="beta (95%CI)")+
  ggtitle(paste0("P for overall association: ",round(anova(fit)[1,3],3),
                 "\nP for non-linearity: ",round(anova(fit)[2,3],3)))+
  theme_classic()+ 
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        title = element_text(color="black",size=10))
ggsave("Metribuzin_RBC_32_rcs.pdf",width = 4,height = 4)

# generalized estimation  equation; pesticide and anemia; Figure 4 -------------------------------------------------
# https://zhuanlan.zhihu.com/p/161497000
# https://zhuanlan.zhihu.com/p/434644253
# lnstruction of seelcting suitable corstr
# https://library.virginia.edu/data/articles/getting-started-with-generalized-estimating-equations

library(gee)
library(geepack)
library(reshape2)

covariates <- c("age","BMI","parity","income","education","smoking_history","folic_acid","iron","alcohol","coffee","GW_24","GW_32")
input_data <- part1_all[,c("id",name_all[c(pesticide,anemia,anemia_parameter)],covariates)]
input_data[,name_all[pesticide]] <- log(input_data[,name_all[pesticide]])

# firstly melt the data frame containing other variables and only one type of variables that need to be melted, GW
input_data_gee <- melt(input_data[,c("id",name_all[c(pesticide)],covariates)],measure.vars = c("GW_24","GW_32"),variable.name = "GW_variable",value.name = "GW")
# secondly melt the data frame containing only one type of variables that need to be melted, GW
anemia_melt <- melt(input_data[,c("id","anemia_24","anemia_32")],measure.vars = c("anemia_24","anemia_32"),variable.name = "anemia_variable",value.name = "anemia")
RBC_melt <- melt(input_data[,c("id","RBC_24","RBC_32")],measure.vars = c("RBC_24","RBC_32"),variable.name = "RBC_variable",value.name = "RBC")
Hb_melt <- melt(input_data[,c("id","Hb_24","Hb_32")],measure.vars = c("Hb_24","Hb_32"),variable.name = "Hb_variable",value.name = "Hb")

input_data_gee <- cbind(input_data_gee,anemia_melt[2:3],RBC_melt[2:3],Hb_melt[2:3])
# write.csv(input_data_gee,"input_data_gee.csv",quote = F,row.names = F)

# anemia 
anemia_gee <- data.frame(matrix(NA,ncol = 6))
names(anemia_gee) <- c("exposure","outcome","OR","l95","u95","p")
for (j in 1:length(pesticide)) {
  m = gee(anemia ~ input_data_gee[,name_all[pesticide[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = binomial,corstr = "exchangeable")
  anemia_gee[j,1] = name_all[pesticide[j]]
  anemia_gee[j,2] = "anemia"
  anemia_gee[j,3] = exp(coef(summary(m))[2,1])
  anemia_gee[j,4] = exp(coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4])
  anemia_gee[j,5] = exp(coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4])
  anemia_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
anemia_gee <- anemia_gee %>% arrange(p)
write.csv(anemia_gee,"pesticide_anemia_gee.csv",quote = F,row.names = F)

# anemia_gee <- read.csv("pesticide_anemia_gee.csv",header = T)
# anemia_gee <- anemia_gee %>% arrange(OR)
# anemia_gee$exposure <- factor(anemia_gee$exposure,levels = anemia_gee$exposure)
# ggplot(anemia_gee,aes(x=exposure,y=log(OR)))+
#   geom_point(size=3,fill="#e3191b",color="#e3191b")+
#   geom_errorbar(aes(ymax=log(u95),ymin=log(l95)),color="#e3191b",size=1,width=0.5)+
#   geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="grey")+
#   labs(x="",y="ln(OR)")+
#   theme_bw()+
#   theme(axis.title=element_text(color="black",size=15),
#         axis.text.x=element_text(color='black',size=15),
#         axis.text.y=element_text(color='black',size=15),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())+
#   coord_flip()
# ggsave("pesticide_anemia_gee.pdf",width = 6,height = 10)

anemia_gee <- read.csv("pesticide_anemia_gee.csv",header = T)
anemia_gee <- anemia_gee %>% arrange(exposure)
anemia_gee$exposure <- factor(anemia_gee$exposure,levels = anemia_gee$exposure)
ggplot(anemia_gee,aes(x=exposure,y=log(OR)))+
  geom_point(size=3,fill="#e3191b",color="#e3191b")+
  geom_errorbar(aes(ymax=log(u95),ymin=log(l95)),color="#e3191b",size=1,width=0.5)+
  geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="black")+
  labs(x="",y="ln(OR)")+
  theme_bw()+
  theme(axis.title=element_text(color="black",size=15),
        axis.text.x=element_text(color='black',size=15,angle = 45,hjust = 1,vjust = 1),
        axis.text.y=element_text(color='black',size=15),
        axis.ticks.x = element_blank())
ggsave("pesticide_anemia_gee.pdf",width = 12,height = 4)

# Hb 
Hb_gee <- data.frame(matrix(NA,ncol = 6))
names(Hb_gee) <- c("exposure","outcome","beta","l95","u95","p")
for (j in 1:length(pesticide)) {
  m = gee(Hb ~ input_data_gee[,name_all[pesticide[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = gaussian,corstr = "exchangeable")
  Hb_gee[j,1] = name_all[pesticide[j]]
  Hb_gee[j,2] = "Hb"
  Hb_gee[j,3] = coef(summary(m))[2,1]
  Hb_gee[j,4] = coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4]
  Hb_gee[j,5] = coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4]
  Hb_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
Hb_gee <- Hb_gee %>% arrange(p)
write.csv(Hb_gee,"pesticide_Hb_gee.csv",quote = F,row.names = F)

# Hb_gee <- read.csv("pesticide_Hb_gee.csv",header = T)
# Hb_gee <- Hb_gee %>% arrange(beta)
# Hb_gee$exposure <- factor(Hb_gee$exposure,levels = Hb_gee$exposure)
# ggplot(Hb_gee,aes(x=exposure,y=beta))+
#   geom_point(size=3,fill="#6a3d99",color="#6a3d99")+
#   geom_errorbar(aes(ymax=u95,ymin=l95),color="#6a3d99",size=1,width=0.5)+
#   geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="grey")+
#   labs(x="",y="beta")+
#   theme_bw()+
#   theme(axis.title=element_text(color="black",size=15),
#         axis.text.x=element_text(color='black',size=15),
#         axis.text.y=element_text(color='black',size=15),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())+
#   coord_flip()
# ggsave("pesticide_Hb_gee.pdf",width = 6,height = 10)

Hb_gee <- read.csv("pesticide_Hb_gee.csv",header = T)
Hb_gee <- Hb_gee %>% arrange(exposure)
Hb_gee$exposure <- factor(Hb_gee$exposure,levels = Hb_gee$exposure)
ggplot(Hb_gee,aes(x=exposure,y=beta))+
  geom_point(size=3,fill="#6a3d99",color="#6a3d99")+
  geom_errorbar(aes(ymax=u95,ymin=l95),color="#6a3d99",size=1,width=0.5)+
  geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="black")+
  labs(x="",y="beta")+
  theme_bw()+
  theme(axis.title=element_text(color="black",size=15),
        axis.text.x=element_text(color='black',size=15,angle = 45,hjust = 1,vjust = 1),
        axis.text.y=element_text(color='black',size=15),
        axis.ticks.x = element_blank())
ggsave("pesticide_Hb_gee.pdf",width = 12,height = 4)

# RBC 
RBC_gee <- data.frame(matrix(NA,ncol = 6))
names(RBC_gee) <- c("exposure","outcome","beta","l95","u95","p")
for (j in 1:length(pesticide)) {
  m = gee(RBC ~ input_data_gee[,name_all[pesticide[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = gaussian,corstr = "exchangeable")
  RBC_gee[j,1] = name_all[pesticide[j]]
  RBC_gee[j,2] = "RBC"
  RBC_gee[j,3] = coef(summary(m))[2,1]
  RBC_gee[j,4] = coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4]
  RBC_gee[j,5] = coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4]
  RBC_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
RBC_gee <- RBC_gee %>% arrange(p)
write.csv(RBC_gee,"pesticide_RBC_gee.csv",quote = F,row.names = F)

# RBC_gee <- read.csv("pesticide_RBC_gee.csv",header = T)
# RBC_gee <- RBC_gee %>% arrange(beta)
# RBC_gee$exposure <- factor(RBC_gee$exposure,levels = RBC_gee$exposure)
# ggplot(RBC_gee,aes(x=exposure,y=beta))+
#   geom_point(size=3,fill="#339f2c",color="#339f2c")+
#   geom_errorbar(aes(ymax=u95,ymin=l95),color="#339f2c",size=1,width=0.5)+
#   geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="grey")+
#   labs(x="",y="beta")+
#   theme_bw()+
#   theme(axis.title=element_text(color="black",size=15),
#         axis.text.x=element_text(color='black',size=15),
#         axis.text.y=element_text(color='black',size=15),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())+
#   coord_flip()
# ggsave("pesticide_RBC_gee.pdf",width = 6,height = 10)

RBC_gee <- read.csv("pesticide_RBC_gee.csv",header = T)
RBC_gee <- RBC_gee %>% arrange(exposure)
RBC_gee$exposure <- factor(RBC_gee$exposure,levels = RBC_gee$exposure)
ggplot(RBC_gee,aes(x=exposure,y=beta))+
  geom_point(size=3,fill="#339f2c",color="#339f2c")+
  geom_errorbar(aes(ymax=u95,ymin=l95),color="#339f2c",size=1,width=0.5)+
  geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="black")+
  labs(x="",y="beta")+
  theme_bw()+
  theme(axis.title=element_text(color="black",size=15),
        axis.text.x=element_text(color='black',size=15,angle = 45,hjust = 1,vjust = 1),
        axis.text.y=element_text(color='black',size=15))
ggsave("pesticide_RBC_gee.pdf",width = 12,height = 4)

# temporal difference analysis; Figure S5 -------------------------------------------------------------

sankey_data <- part1_all[,name_all[anemia]]
rownames(sankey_data) <- part1_all$id
sankey_data$type <- ifelse(sankey_data$anemia_24==0 & sankey_data$anemia_32==0,"1",ifelse(
  sankey_data$anemia_24==0 & sankey_data$anemia_32==1,"2",ifelse(
    sankey_data$anemia_24==1 & sankey_data$anemia_32==0,"3" ,"4")))
table(sankey_data$anemia_24,sankey_data$anemia_32)
prop.table(table(sankey_data$anemia_24,sankey_data$anemia_32),margin = 1)
chisq.test(table(sankey_data$anemia_24,sankey_data$anemia_32))

# https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram
library(networkD3)
sankey_data <- part1_all[,name_all[anemia]]
rownames(sankey_data) <- part1_all$id
links <- data.frame(table(sankey_data))
links$anemia_24 <- ifelse(links$anemia_24==0,"Non-GA (T2)","GA (T2)")
links$anemia_32 <- ifelse(links$anemia_32==0,"Non-GA (T3)","GA (T3)")
nodes <- data.frame(name=c(links$anemia_24,links$anemia_32)%>%unique())
links$source <- match(links$anemia_24,nodes$name)-1
links$target <- match(links$anemia_32,nodes$name)-1
nodes$name <- as.factor(nodes$name)
mycolor <- 'd3.scaleOrdinal().domain(["Non-GA (T2)","GA (T2)","Non-GA (T3)","GA (T3)"]).range(["#fc6883","#007cc4","#ef7d0a","#00a03d"])'

# links$link_color <- c("#feb6bc","#81bedc","#f7bf84","#81d09a")
# nodes$nodes_color <- c("#fc6883","#007cc4","#ef7d0a","#00a03d")
p <- sankeyNetwork(Links = links,Nodes = nodes,Source = "source",Target = "target",Value = "Freq",NodeID = "name",sinksRight = F,
              nodeWidth = 10,fontSize = 20,fontFamily = "Arial",colourScale = mycolor,height = "500px",width = "500px")
p
library(htmlwidgets)
saveWidget(p, "sankey.html")

# cluster1 vs cluster2, pesticide
input_data <- part1_all
input_data$type <- as.numeric(sankey_data$type)

input_data <- input_data[input_data$type==1 |input_data$type==2, ]
input_data$type <- input_data$type-1

result <- data.frame(matrix(NA,ncol = 6))
names(result) <- c("exposure","outcome","OR","l95","u95","p")
for (j in 1:length(pesticide)) {
  m = glm(input_data$type ~ log(input_data[,name_all[pesticide[j]]])+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24+GW_32
          ,data = input_data, family = binomial())
  result[j,1] = name_all[pesticide[j]]
  result[j,2] = "type"
  result[j,3] = exp(summary(m)$coefficients[2,1])
  result[j,4] = exp(confint(m)[2,1])
  result[j,5] = exp(confint(m)[2,2])
  result[j,6] = summary(m)$coefficients[2,4]
}
result <- result %>% arrange(p)
write.csv(result,"pesticide_anemia_cluster1_2.csv",quote = F,row.names = F)

# result <- read.csv("pesticide_anemia_cluster1_2.csv",header = T)
# result <- result %>% arrange(OR)
# result$exposure <- factor(result$exposure,levels = result$exposure)
# ggplot(result,aes(x=exposure,y=log(OR)))+
#   geom_point(size=3,fill="#e3191b",color="#e3191b")+
#   geom_errorbar(aes(ymax=log(u95),ymin=log(l95)),color="#e3191b",size=1,width=0.5)+
#   geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="grey")+
#   labs(x="",y="ln(OR)")+
#   theme_bw()+
#   theme(axis.title=element_text(color="black",size=15),
#         axis.text.x=element_text(color='black',size=15),
#         axis.text.y=element_text(color='black',size=15),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())+
#   coord_flip()
# ggsave("pesticide_anemia_cluster1_2.pdf",width = 6,height = 10)

result <- read.csv("pesticide_anemia_cluster1_2.csv",header = T)
result <- result %>% arrange(exposure)
result$exposure <- factor(result$exposure,levels = result$exposure)
ggplot(result,aes(x=exposure,y=log(OR)))+
  geom_point(size=3,fill="#e3191b",color="#e3191b")+
  geom_errorbar(aes(ymax=log(u95),ymin=log(l95)),color="#e3191b",size=1,width=0.5)+
  geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="black")+
  labs(x="",y="ln(OR)")+
  theme_bw()+
  theme(axis.title=element_text(color="black",size=15),
        axis.text.x=element_text(color='black',size=15,angle = 45,hjust = 1,vjust = 1),
        axis.text.y=element_text(color='black',size=15),
        axis.ticks.x = element_blank())
ggsave("pesticide_anemia_cluster1_2.pdf",width = 12,height = 4)

# cluster3 vs cluster4, pesticide
input_data <- part1_all
input_data$type <- as.numeric(sankey_data$type)

input_data <- input_data[input_data$type==3 |input_data$type==4, ]
input_data$type[input_data$type==3] <- input_data$type[input_data$type==3]-2
input_data$type[input_data$type==4] <- input_data$type[input_data$type==4]-4

result <- data.frame(matrix(NA,ncol = 6))
names(result) <- c("exposure","outcome","OR","l95","u95","p")
for (j in 1:length(pesticide)) {
  m = glm(input_data$type ~ log(input_data[,name_all[pesticide[j]]])+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24+GW_32
          ,data = input_data, family = binomial())
  result[j,1] = name_all[pesticide[j]]
  result[j,2] = "type"
  result[j,3] = exp(summary(m)$coefficients[2,1])
  result[j,4] = exp(confint(m)[2,1])
  result[j,5] = exp(confint(m)[2,2])
  result[j,6] = summary(m)$coefficients[2,4]
}
result <- result %>% arrange(p)
write.csv(result,"pesticide_anemia_cluster3_4.csv",quote = F,row.names = F)

# result <- read.csv("pesticide_anemia_cluster3_4.csv",header = T)
# result <- result %>% arrange(OR)
# result$exposure <- factor(result$exposure,levels = result$exposure)
# ggplot(result,aes(x=exposure,y=log(OR)))+
#   geom_point(size=3,fill="#e3191b",color="#e3191b")+
#   geom_errorbar(aes(ymax=log(u95),ymin=log(l95)),color="#e3191b",size=1,width=0.5)+
#   geom_hline(aes(yintercept=0),lty=2,linewidth=1,color="black")+
#   labs(x="",y="ln(OR)")+
#   theme_bw()+
#   theme(axis.title=element_text(color="black",size=15),
#         axis.text.x=element_text(color='black',size=15),
#         axis.text.y=element_text(color='black',size=15),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())+
#   coord_flip()
# ggsave("pesticide_anemia_cluster3_4.pdf",width = 6,height = 10)


# generalized estimation  equation; gut microbiota and anemia; Figure 6A-C -------------------------------------------------
# https://zhuanlan.zhihu.com/p/161497000
# https://zhuanlan.zhihu.com/p/434644253
# Instruction of selecting suitable corstr
# https://library.virginia.edu/data/articles/getting-started-with-generalized-estimating-equations

library(gee)
library(geepack)
library(reshape2)

covariates <- c("age","BMI","parity","income","education","smoking_history","folic_acid","iron","alcohol","coffee","GW_24","GW_32")
input_data <- part1_all[,c("id",name_all[c(microbiota,microbiota_parameter,anemia,anemia_parameter)],covariates)]
input_data[,name_all[c(microbiota)]] <- clr(input_data[,name_all[c(microbiota)]])

# firstly melt the data frame containing other variables and only one type of variables that need to be melted (GW)
input_data_gee <- melt(input_data[,c("id",name_all[c(microbiota,microbiota_parameter)],covariates)],measure.vars = c("GW_24","GW_32"),variable.name = "GW_variable",value.name = "GW")
# secondly melt the data frame containing only one type of variables that need to be melted (anemia, RBC, Hb)
anemia_melt <- melt(input_data[,c("id","anemia_24","anemia_32")],measure.vars = c("anemia_24","anemia_32"),variable.name = "anemia_variable",value.name = "anemia")
RBC_melt <- melt(input_data[,c("id","RBC_24","RBC_32")],measure.vars = c("RBC_24","RBC_32"),variable.name = "RBC_variable",value.name = "RBC")
Hb_melt <- melt(input_data[,c("id","Hb_24","Hb_32")],measure.vars = c("Hb_24","Hb_32"),variable.name = "Hb_variable",value.name = "Hb")

input_data_gee <- cbind(input_data_gee,anemia_melt[2:3],RBC_melt[2:3],Hb_melt[2:3])
# write.csv(input_data_gee,"input_data_gutmicrobiota_gee.csv",quote = F,row.names = F)

# anemia microbiota
anemia_gee <- data.frame(matrix(NA,ncol = 6))
names(anemia_gee) <- c("exposure","outcome","OR","l95","u95","p")
for (j in 1:length(microbiota)) {
  m = gee(anemia ~ input_data_gee[,name_all[microbiota[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = binomial,corstr = "exchangeable")
  anemia_gee[j,1] = name_all[microbiota[j]]
  anemia_gee[j,2] = "anemia"
  anemia_gee[j,3] = exp(coef(summary(m))[2,1])
  anemia_gee[j,4] = exp(coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4])
  anemia_gee[j,5] = exp(coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4])
  anemia_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
anemia_gee$padj <- p.adjust(anemia_gee$p)
anemia_gee <- anemia_gee %>% arrange(padj)
write.csv(anemia_gee,"microbiota_anemia_gee.csv",quote = F,row.names = F)

# anemia microbiota_parameter
anemia_gee <- data.frame(matrix(NA,ncol = 6))
names(anemia_gee) <- c("exposure","outcome","OR","l95","u95","p")
for (j in 1:length(microbiota_parameter)) {
  m = gee(anemia ~ input_data_gee[,name_all[microbiota_parameter[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = binomial,corstr = "exchangeable")
  anemia_gee[j,1] = name_all[microbiota_parameter[j]]
  anemia_gee[j,2] = "anemia"
  anemia_gee[j,3] = exp(coef(summary(m))[2,1])
  anemia_gee[j,4] = exp(coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4])
  anemia_gee[j,5] = exp(coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4])
  anemia_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
anemia_gee <- anemia_gee %>% arrange(p)
write.csv(anemia_gee,"microbiota_parameter_anemia_gee.csv",quote = F,row.names = F)

# Hb 
Hb_gee <- data.frame(matrix(NA,ncol = 6))
names(Hb_gee) <- c("exposure","outcome","beta","l95","u95","p")
for (j in 1:length(microbiota)) {
  m = gee(Hb ~ input_data_gee[,name_all[microbiota[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = gaussian,corstr = "exchangeable")
  Hb_gee[j,1] = name_all[microbiota[j]]
  Hb_gee[j,2] = "Hb"
  Hb_gee[j,3] = coef(summary(m))[2,1]
  Hb_gee[j,4] = coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4]
  Hb_gee[j,5] = coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4]
  Hb_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
Hb_gee$padj <- p.adjust(Hb_gee$p)
Hb_gee <- Hb_gee %>% arrange(padj)
write.csv(Hb_gee,"microbiota_Hb_gee.csv",quote = F,row.names = F)

Hb_gee <- data.frame(matrix(NA,ncol = 6))
names(Hb_gee) <- c("exposure","outcome","beta","l95","u95","p")
for (j in 1:length(microbiota_parameter)) {
  m = gee(Hb ~ input_data_gee[,name_all[microbiota_parameter[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = gaussian,corstr = "exchangeable")
  Hb_gee[j,1] = name_all[microbiota_parameter[j]]
  Hb_gee[j,2] = "Hb"
  Hb_gee[j,3] = coef(summary(m))[2,1]
  Hb_gee[j,4] = coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4]
  Hb_gee[j,5] = coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4]
  Hb_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
Hb_gee <- Hb_gee %>% arrange(p)
write.csv(Hb_gee,"microbiota_parameter_Hb_gee.csv",quote = F,row.names = F)

# RBC 
RBC_gee <- data.frame(matrix(NA,ncol = 6))
names(RBC_gee) <- c("exposure","outcome","beta","l95","u95","p")
for (j in 1:length(microbiota)) {
  m = gee(RBC ~ input_data_gee[,name_all[microbiota[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = gaussian,corstr = "exchangeable")
  RBC_gee[j,1] = name_all[microbiota[j]]
  RBC_gee[j,2] = "RBC"
  RBC_gee[j,3] = coef(summary(m))[2,1]
  RBC_gee[j,4] = coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4]
  RBC_gee[j,5] = coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4]
  RBC_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
RBC_gee$padj <- p.adjust(RBC_gee$p)
RBC_gee <- RBC_gee %>% arrange(padj)
write.csv(RBC_gee,"microbiota_RBC_gee.csv",quote = F,row.names = F)

RBC_gee <- data.frame(matrix(NA,ncol = 6))
names(RBC_gee) <- c("exposure","outcome","beta","l95","u95","p")
for (j in 1:length(microbiota_parameter)) {
  m = gee(RBC ~ input_data_gee[,name_all[microbiota_parameter[j]]]+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW
          ,id=id,data = input_data_gee, family = gaussian,corstr = "exchangeable")
  RBC_gee[j,1] = name_all[microbiota_parameter[j]]
  RBC_gee[j,2] = "RBC"
  RBC_gee[j,3] = coef(summary(m))[2,1]
  RBC_gee[j,4] = coef(summary(m))[2,1]-1.96*coef(summary(m))[2,4]
  RBC_gee[j,5] = coef(summary(m))[2,1]+1.96*coef(summary(m))[2,4]
  RBC_gee[j,6] = 2*pnorm(abs(coef(summary(m))[2,5]),lower.tail = F)
}
RBC_gee <- RBC_gee %>% arrange(p)
write.csv(RBC_gee,"microbiota_parameter_RBC_gee.csv",quote = F,row.names = F)

# Manhattan plot; only include microbiota; also include longitudinal anemia parameter Figure 6D
result1 <- read.csv("microbiota_anemia_parameter_T2_part1.csv",header = T)
result2 <- read.csv("microbiota_anemia_T2_part1.csv",header = T);result2$outcome <- "Anemia_24"
result3 <- read.csv("microbiota_anemia_parameter_T3_part1.csv",header = T)
result4 <- read.csv("microbiota_anemia_T3_part1.csv",header = T);result4$outcome <- "Anemia_32"
result5 <- read.csv("microbiota_Hb_gee.csv",header = T);result5$outcome <- "ZH_Hb" # ZH最后，根据字母排序，让这些选项出现在最后
result6 <- read.csv("microbiota_RBC_gee.csv",header = T);result6$outcome <- "ZH_RBC"
result7 <- read.csv("microbiota_anemia_gee.csv",header = T);result7$outcome <- "ZH_Anemia"

result1$sig <- ifelse(result1$padj>=0.05,"non-sig",ifelse(result1$padj<0.05 & result1$beta>0,"positive","negative"))
result2$sig <- ifelse(result2$padj>=0.05,"non-sig",ifelse(result2$padj<0.05 & result2$OR>0,"positive","negative"))
result3$sig <- ifelse(result3$padj>=0.05,"non-sig",ifelse(result3$padj<0.05 & result3$beta>0,"positive","negative"))
result4$sig <- ifelse(result4$padj>=0.05,"non-sig",ifelse(result4$padj<0.05 & result4$OR>0,"positive","negative"))
result5$sig <- ifelse(result5$padj>=0.05,"non-sig",ifelse(result5$padj<0.05 & result5$beta>0,"positive","negative"))
result6$sig <- ifelse(result6$padj>=0.05,"non-sig",ifelse(result6$padj<0.05 & result6$beta>0,"positive","negative"))
result7$sig <- ifelse(result7$padj>=0.05,"non-sig",ifelse(result7$padj<0.05 & result7$OR>0,"positive","negative"))

result <- rbind(result1[,-3],result2[,-3],result3[,-3],result4[,-3],result5[,-3],result6[,-3],result7[,-3])
result <- result %>% arrange(outcome,exposure)
for (i in 1:length(c(anemia,anemia_parameter,"anemia_gee","Hb_gee","RBC_gee"))){
  for (j in 1:length(c(microbiota)))
    result$position[100*(i-1)+j]=130*(i-1)+j
}
position_break <- aggregate(result$position, by = list(result$outcome), FUN = median)
result$rank <- rep(rep(c(1,2,3),each=length(c(microbiota))),3) %>% factor()
result[result$padj<0.05,]$exposure <- sub(".*;","",result[result$padj<0.05,]$exposure)

p <- ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "black",linetype=2,size = 0.35) +
  geom_point(data=result[result$sig=="non-sig",],aes(position, -log(padj, 10),fill=rank),size=2,alpha=0.6,show.legend = FALSE,shape=21,stroke=0.1) + # stroke控制描边
  geom_point(data=result[result$sig!="non-sig",],aes(position, -log(padj, 10),color = sig), size=2,alpha=0.8,shape = 19,stroke=0.1) +
  scale_fill_manual(values = c("#857bbb","#8b9071","#a96976")) +
  scale_color_manual(values = c("#238dce","#d12918")) +
  # scale_color_manual(values = rep(c("#85bce5","#fad76f","#857bbb","#8b9071","#e85fa2","#ef9b6a","#e4aed0","#b0e3c9","#a96976"),times=4)) +
  scale_x_continuous(breaks = position_break$x, labels = position_break$Group.1, expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0,1.0,-log10(0.05),2.0,3.0), labels = c("0.0","1.0","","2.0","3.0")) +
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black'), 
        panel.background = element_rect(fill = 'transparent'),
        axis.text.x = element_text(color = "black",size=6,angle=45,hjust = 1,vjust = 1),
        axis.text.y = element_text(color = "black",size=6),
        axis.title = element_text(color = "black",size=6)) +
  labs(x = '', y = expression(''~-log[10]~'(FDR-corrected P)'))+
  geom_text_repel(data=result[result$sig!="non-sig",],aes(position, -log(padj, 10),label=exposure), size = 2,box.padding = unit(0.3, 'lines'), segment.color = "black", show.legend = FALSE, color = "black")
p
ggsave("bacteria_anemia.pdf",width = 5,height = 3)


# LEfSe; Figure 6E-F  ------------------------------------------------------------------

count_part1 <- read.xlsx("data_part1.xlsx",sheet = 2,rowNames = T)
mediators <- count_part1[,name_all[microbiota]]
# 需要为count数据，而非relative abundance数据
rownames(mediators) <- paste0("id",rownames(mediators))
tree <- data.frame(name_all[microbiota]);colnames(tree) <- c("bacteria")
tree <- tree %>% separate(col = "bacteria",into = c("Kingdom","Phylum","Class","Order","Family","Genus"),sep = ";")
tree$Family[which(tree$Family=="__")] <- paste0("f__Unclassified",c(1:length(which(tree$Family=="__"))))
tree$Genus[which(tree$Genus=="__")] <- paste0("g__Unclassified",c(1:length(which(tree$Genus=="__"))))
tree$Genus[which(tree$Genus=="g__uncultured")] <- paste0("g__Uncultured",c(1:length(which(tree$Genus=="g__uncultured"))))
# 同一层级下物种名称需要不一致，同时需要保留各层次的前缀
tree$Kingdom <- gsub("^d__","k__",tree$Kingdom) 
rownames(tree) <- tree$Genus
colnames(mediators) <- tree$Genus
mediators <- data.frame(t(mediators))
sample <- part1_all[,c("id",name_all[anemia])]
rownames(sample) <- paste0("id",part1_all$id)
#创建microeco包可识别的对象
dataset <- microtable$new(otu_table = mediators,sample_table = sample,tax_table = tree)

# anemia_24
#开始LEfse分析
lefse <- trans_diff$new(dataset = dataset,method = "lefse",group = "anemia_24",alpha = 0.05,p_adjust_method="none",lefse_subgroup = NULL)
#查看组间差异结果
View(lefse$res_diff)
write.csv(lefse$res_diff,"lefse_diff_anemia24_part1.csv",quote = F,row.names = F)
#后续只挑选属水平细菌进行分析
# 绘制前10个具有最高LDA（log10）的分类单元的差异特征柱状图
lefse$plot_diff_bar(color_values = c("#aad9c9","#f7beaa"),width=0.6,group_order = c("0","1"))+
  geom_hline(aes(yintercept = 0),lty=1,linewidth=1,color="black")+
  geom_hline(aes(yintercept = -2),lty=2,linewidth=1,color="grey")+
  geom_hline(aes(yintercept = 2),lty=2,linewidth=1,color="grey")+
  theme_classic()+
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        legend.title = element_blank())
ggsave("diff_bar_anemia24_part1.pdf",width = 5,height = 4)

# 展示前200个分类单元和前50个特征的分支进化图
lefse$plot_diff_cladogram(color = c("#aad9c9","#f7beaa"),group_order = c("0","1"),clade_label_level = 6)
ggsave("diff_cladogram_anemia24_part1.pdf",width = 6,height = 4)

# anemia_32
#开始LEfse分析
lefse <- trans_diff$new(dataset = dataset,method = "lefse",group = "anemia_32",alpha = 0.05,p_adjust_method="none",lefse_subgroup = NULL)
#查看组间差异结果
View(lefse$res_diff)
write.csv(lefse$res_diff,"lefse_diff_anemia32_part1.csv",quote = F,row.names = F)
# 绘制前10个具有最高LDA（log10）的分类单元的差异特征柱状图
lefse$plot_diff_bar(color_values = c("#aad9c9","#f7beaa"),width=0.6,group_order = c("0","1"))+
  geom_hline(aes(yintercept = 0),lty=1,linewidth=1,color="black")+
  geom_hline(aes(yintercept = -2),lty=2,linewidth=1,color="grey")+
  geom_hline(aes(yintercept = 2),lty=2,linewidth=1,color="grey")+
  theme_classic()+
  theme(axis.title = element_text(color="black",size=10),
        axis.text = element_text(color="black",size=10),
        legend.title = element_blank())
ggsave("diff_bar_anemia32_part1.pdf",width = 5,height = 4)

# 展示前200个分类单元和前50个特征的分支进化图
lefse$plot_diff_cladogram(color = c("#aad9c9","#f7beaa"),group_order = c("0","1"),clade_label_level = 6)
ggsave("diff_cladogram_anemia32_part1.pdf",width = 6,height = 4)

# mediation analysis; Figure 7, Table S8 ------------------------------------------------------

# T2 shannon

# Diphenamid-shannon-RBC_24
m1 <- glm(shannon ~ log(Diphenamid)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
m2 <- glm(RBC_24 ~ log(Diphenamid)+shannon+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Diphenamid)",mediator = "shannon")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME            1.17e-03     3.89e-05         0.00    0.04 *
#   ADE             8.55e-04    -1.21e-02         0.01    0.88  
# Total Effect    2.03e-03    -1.09e-02         0.01    0.73  
# Prop. Mediated  9.05e-02    -2.75e+00         2.73    0.74  

m3 <- glm(RBC_24 ~ shannon+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
summary(m1)
summary(m2)
summary(m3)

# Metribuzin-shannon-RBC_24
m1 <- glm(shannon ~ log(Metribuzin)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
m2 <- glm(RBC_24 ~ log(Metribuzin)+shannon+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Metribuzin)",mediator = "shannon")
summary(model)

# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME           -0.00124     -0.00329         0.00   0.046 *
#   ADE            -0.00327     -0.01826         0.01   0.678  
# Total Effect   -0.00450     -0.01950         0.01   0.550  
# Prop. Mediated  0.09056     -2.46066         2.64   0.576  

m3 <- glm(RBC_24 ~ shannon+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
summary(m1)
summary(m2)
summary(m3)


# Dimethipin-shannon-RBC_24
m1 <- glm(shannon ~ log(Dimethipin)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
m2 <- glm(RBC_24 ~ log(Dimethipin)+shannon+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Dimethipin)",mediator = "shannon")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME            2.57e-03     3.43e-05         0.01   0.048 *
#   ADE             1.35e-03    -2.85e-02         0.03   0.922  
# Total Effect    3.93e-03    -2.61e-02         0.03   0.794  
# Prop. Mediated  6.86e-02    -3.00e+00         1.92   0.798 

m3 <- glm(RBC_24 ~ shannon+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_24,data = part1_all, family = gaussian())
summary(m1)
summary(m2)
summary(m3)

# T3 Roseburia

# Dicrotophos-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Dicrotophos)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Dicrotophos)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Dicrotophos)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.00933     -0.01898         0.00   0.016 *
#   ACME (treated)           -0.00943     -0.01983         0.00   0.016 *
#   ADE (control)             0.01443     -0.03637         0.08   0.682  
# ADE (treated)             0.01432     -0.03516         0.08   0.682  
# Total Effect              0.00499     -0.04319         0.07   0.922  
# Prop. Mediated (control) -0.11565     -6.62255         4.65   0.922  
# Prop. Mediated (treated) -0.12112     -6.78296         4.74   0.922  
# ACME (average)           -0.00938     -0.01940         0.00   0.016 *
#   ADE (average)             0.01437     -0.03573         0.08   0.682  
# Prop. Mediated (average) -0.11839     -6.71163         4.69   0.922  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Paclobutrazol-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Paclobutrazol)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Paclobutrazol)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Paclobutrazol)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)            -0.0172      -0.0370         0.00   0.018 *
#   ACME (treated)            -0.0173      -0.0394         0.00   0.018 *
#   ADE (control)              0.0410      -0.0595         0.17   0.602  
# ADE (treated)              0.0409      -0.0565         0.18   0.602  
# Total Effect               0.0237      -0.0682         0.16   0.796  
# Prop. Mediated (control)  -0.1357      -5.2592         4.19   0.798  
# Prop. Mediated (treated)  -0.1363      -5.4641         4.33   0.798  
# ACME (average)            -0.0173      -0.0377         0.00   0.018 *
#   ADE (average)              0.0410      -0.0576         0.18   0.602  
# Prop. Mediated (average)  -0.1360      -5.3964         4.26   0.798  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Parathion.methyl-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Parathion.methyl)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Parathion.methyl)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Parathion.methyl)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.00392     -0.00839         0.00   0.016 *
#   ACME (treated)           -0.00401     -0.00855         0.00   0.016 *
#   ADE (control)             0.01323     -0.01706         0.04   0.410  
# ADE (treated)             0.01314     -0.01684         0.04   0.410  
# Total Effect              0.00921     -0.02017         0.04   0.534  
# Prop. Mediated (control) -0.16253     -4.18340         2.64   0.546  
# Prop. Mediated (treated) -0.16917     -4.23760         2.65   0.546  
# ACME (average)           -0.00397     -0.00847         0.00   0.016 *
#   ADE (average)             0.01318     -0.01699         0.04   0.410  
# Prop. Mediated (average) -0.16585     -4.21505         2.65   0.546  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Clomazone-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Clomazone)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Clomazone)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Clomazone)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)

# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)            -0.0133      -0.0290         0.00   0.016 *
#   ACME (treated)            -0.0135      -0.0311         0.00   0.016 *
#   ADE (control)              0.0535      -0.0451         0.17   0.396  
# ADE (treated)              0.0533      -0.0435         0.17   0.396  
# Total Effect               0.0400      -0.0538         0.16   0.524  
# Prop. Mediated (control)  -0.1459      -3.1575         3.02   0.536  
# Prop. Mediated (treated)  -0.1520      -3.2971         3.10   0.536  
# ACME (average)            -0.0134      -0.0297         0.00   0.016 *
#   ADE (average)              0.0534      -0.0442         0.17   0.396  
# Prop. Mediated (average)  -0.1490      -3.2274         3.06   0.536  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Pyrimethanil-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Pyrimethanil)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Pyrimethanil)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Pyrimethanil)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)

# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.00888     -0.01883         0.00   0.012 *
#   ACME (treated)           -0.00932     -0.01993         0.00   0.012 *
#   ADE (control)             0.07474      0.00313         0.15   0.038 *
#   ADE (treated)             0.07430      0.00308         0.15   0.038 *
#   Total Effect              0.06542     -0.00473         0.14   0.068 .
# Prop. Mediated (control) -0.12689     -1.12789         0.45   0.080 .
# Prop. Mediated (treated) -0.13374     -1.15868         0.44   0.080 .
# ACME (average)           -0.00910     -0.01932         0.00   0.012 *
#   ADE (average)             0.07452      0.00311         0.15   0.038 *
#   Prop. Mediated (average) -0.13031     -1.14329         0.44   0.080 .


m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Diphenamid-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Diphenamid)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Diphenamid)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Diphenamid)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.00248     -0.00537         0.00   0.014 *
#   ACME (treated)           -0.00253     -0.00551         0.00   0.014 *
#   ADE (control)             0.01268     -0.00680         0.03   0.232  
# ADE (treated)             0.01263     -0.00680         0.03   0.232  
# Total Effect              0.01015     -0.00964         0.03   0.350  
# Prop. Mediated (control) -0.15273     -1.63710         2.81   0.360  
# Prop. Mediated (treated) -0.15645     -1.64852         2.82   0.360  
# ACME (average)           -0.00250     -0.00544         0.00   0.014 *
#   ADE (average)             0.01266     -0.00680         0.03   0.232  
# Prop. Mediated (average) -0.15459     -1.64281         2.81   0.360  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Monocrotophos-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Monocrotophos)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Monocrotophos)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Monocrotophos)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.00913     -0.02038         0.00   0.016 *
#   ACME (treated)           -0.00928     -0.02102         0.00   0.016 *
#   ADE (control)             0.05650     -0.02748         0.15   0.256  
# ADE (treated)             0.05635     -0.02697         0.15   0.256  
# Total Effect              0.04722     -0.03372         0.14   0.352  
# Prop. Mediated (control) -0.12330     -1.35811         1.73   0.364  
# Prop. Mediated (treated) -0.12384     -1.37610         1.75   0.364  
# ACME (average)           -0.00920     -0.02081         0.00   0.016 *
#   ADE (average)             0.05642     -0.02723         0.15   0.256  
# Prop. Mediated (average) -0.12357     -1.36711         1.74   0.364  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Atrazine-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Atrazine)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Atrazine)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Atrazine)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.002188    -0.004907         0.00   0.016 *
#   ACME (treated)           -0.002271    -0.005061         0.00   0.016 *
#   ADE (control)             0.021831     0.002491         0.04   0.024 *
#   ADE (treated)             0.021748     0.002483         0.04   0.024 *
#   Total Effect              0.019561     0.000511         0.04   0.044 *
#   Prop. Mediated (control) -0.104454    -0.931901         0.02   0.060 .
# Prop. Mediated (treated) -0.108341    -0.937944         0.02   0.060 .
# ACME (average)           -0.002229    -0.004983         0.00   0.016 *
#   ADE (average)             0.021790     0.002487         0.04   0.024 *
#   Prop. Mediated (average) -0.106398    -0.934923         0.02   0.060 .

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2);exp(summary(m2)$coefficients[2,1])
summary(m3);exp(summary(m3)$coefficients[2,1])

# Alachlor-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Alachlor)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Alachlor)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Alachlor)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.00429     -0.00986         0.00   0.022 *
#   ACME (treated)           -0.00427     -0.00977         0.00   0.022 *
#   ADE (control)            -0.00356     -0.04295         0.04   0.870  
# ADE (treated)            -0.00355     -0.04262         0.04   0.870  
# Total Effect             -0.00783     -0.04801         0.03   0.700  
# Prop. Mediated (control)  0.11471     -3.46733         3.37   0.710  
# Prop. Mediated (treated)  0.11014     -3.52465         3.40   0.710  
# ACME (average)           -0.00428     -0.00986         0.00   0.022 *
#   ADE (average)            -0.00355     -0.04278         0.04   0.870  
# Prop. Mediated (average)  0.11242     -3.49599         3.39   0.710  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# Monolinuron-Roseburia-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`)) ~ log(Monolinuron)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Monolinuron)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Monolinuron)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.00359     -0.00859         0.00   0.028 *
#   ACME (treated)           -0.00346     -0.00818         0.00   0.028 *
#   ADE (control)            -0.01842     -0.05169         0.02   0.298  
# ADE (treated)            -0.01829     -0.05108         0.02   0.298  
# Total Effect             -0.02188     -0.05453         0.01   0.214  
# Prop. Mediated (control)  0.12195     -1.24942         1.53   0.238  
# Prop. Mediated (treated)  0.11500     -1.26679         1.53   0.238  
# ACME (average)           -0.00353     -0.00839         0.00   0.028 *
#   ADE (average)            -0.01836     -0.05138         0.02   0.298  
# Prop. Mediated (average)  0.11847     -1.25810         1.53   0.238  

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# T3 Oscillibacter

# Nuarimol-Oscillibacter-anemia_32
m1 <- glm(as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Oscillibacter`)) ~ log(Nuarimol)+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = gaussian())
m2 <- glm(anemia_32 ~ log(Nuarimol)+as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Oscillibacter`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
set.seed(133)
model <- mediation::mediate(m1,m2,sims = 1000,treat = "log(Nuarimol)",mediator = "as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Oscillibacter`))")
summary(model)
# Estimate 95% CI Lower 95% CI Upper p-value
# ACME (control)            0.00212     -0.00623         0.01    0.62
# ACME (treated)            0.00182     -0.00567         0.01    0.62
# ADE (control)            -0.03832     -0.07251         0.02    0.14
# ADE (treated)            -0.03862     -0.07327         0.02    0.14
# Total Effect             -0.03650     -0.07058         0.02    0.15
# Prop. Mediated (control) -0.04179     -0.72531         0.43    0.68
# Prop. Mediated (treated) -0.03447     -0.69355         0.43    0.68
# ACME (average)            0.00197     -0.00604         0.01    0.62
# ADE (average)            -0.03847     -0.07279         0.02    0.14
# Prop. Mediated (average) -0.03813     -0.70649         0.43    0.68

m3 <- glm(anemia_32 ~ as.numeric(clr(`d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Oscillibacter`))+age+BMI+as.factor(parity)+income+as.factor(education)+as.factor(smoking_history)+as.factor(folic_acid)+as.factor(iron)+as.factor(alcohol)+as.factor(coffee)+GW_32,data = part1_all, family = binomial())
summary(m1)
summary(m2)
summary(m3)

# 绘图

links <- read.csv("mediation_linkages.csv")
nodes <- data.frame(name=c(links$source,links$target)%>%unique())
links$IDsource <- match(links$source,nodes$name)-1
links$IDtarget <- match(links$target,nodes$name)-1
p <- sankeyNetwork(Links = links, # 输入数据1
                   Nodes = nodes, # 输入数据2
                   Source = "IDsource", # 来源变量
                   Target = "IDtarget", # 接受变量
                   Value = "weight", # 关系权重
                   NodeID = "name", #节点名称
                   LinkGroup = 'source', # 颜色分组
                   sinksRight = FALSE, # 设置最后一层标签位置在左/右
                   nodeWidth = 10, #节点格子宽度
                   fontSize = 15, #文本标签字体的大小
                   nodePadding = 4,
                   fontFamily = "Arial",
                   height = "400px",width = "500px") #节点格子间空隙宽度
p
# 保存   
saveWidget(p,"sankey_mediation.html")






