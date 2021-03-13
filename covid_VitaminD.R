#### Cargamos datos ####
library(dummies);library(readr); library(tidyverse); library(survival); library(mediation); library(ggpubr); library(rms); library(CAinterprTools)
library(survminer); library(haven); library(rsq); library(ResourceSelection); library(ggsci);library(timereg); library(coxme); library(FactoMineR)
library(pROC);library(sf); library(rgdal); library(ggsci); library(ggmap); library(scales); library(jtools); library(cowplot); library(factoextra)
library(ggstance); library(flextable); library(simPH); library(ggthemes); library(lme4); library(lmerTest); library(prismatic); library(viridis)
library(simPH); library(dplyr); library(MatchIt)
setwd("C:/Users/Usuario Lab. Datos/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID- Roopa")
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID- Roopa/Vitamina D")
setwd("/Users/nefoantonio/UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/OMAR YAXMEHEN BELLO CHAVOLLA - COVID- Roopa")


data <-read_sav("COVID 28 SEPTIEMBRE 2020.sav")
data$DIAS_INTERNAMIENTO[is.na(data$DIAS_INTERNAMIENTO)]<-0
data$INICIODESINTOMAS[is.na(data$INICIODESINTOMAS)]<-0
data$INTUBACION[is.na(data$INTUBACION)]<-0
data$EDO_ALTA[is.na(data$EDO_ALTA)]<-1
data$defuncion[data$EDO_ALTA==0]<-1;data$defuncion[data$EDO_ALTA==1]<-0
data$critical<-data$defuncion+data$INTUBACION
data$critical<-factor(data$critical, labels = c("Non-critical", "Critical", "Lethal"))
data$fu_time<-data$INICIODESINTOMAS+data$DIAS_INTERNAMIENTO
data$pericardial<-data$PROMEDIO_AVD+data$PROMEDIO_AVI+data$PROMEDIO_SIV
data$defuncion1<-factor(data$defuncion, labels = c("Non-lethal", "Lethal"))
data$intubacion1<-factor(data$INTUBACION, labels = c("No-IVM", "IVM"))
data$GENERO<-factor(data$GENERO, labels = c("Female", "Male"))
data$critical1[(data$defuncion+data$INTUBACION)>0]<-1;data$critical1[is.na(data$critical1)]<-0
data$vitd20[data$VitD<20]<-1;data$vitd20[data$VitD>=20]<-0
data$vitd12[data$VitD<12]<-1;data$vitd12[data$VitD>=12]<-0
data$imc_cat2[data$IMC>30]<-1;data$imc_cat2[data$IMC<25]<-0;data$imc_cat2[data$IMC<30 & data$IMC>=25]<-0
data$imc_cat2<-factor(data$imc_cat2, labels = c("Normal weight", "Obesity"))
data<-data %>% filter(VitD>0)
data$epi_ob[data$GROSOR_PROMEDIO>12.3]<-1;data$epi_ob[data$GROSOR_PROMEDIO<=12.3]<-0
data$num_sintomas<-data$FIEBRE+data$TOS_SECA+data$DISNEA+data$MIALGIAS+data$CEFALEA+data$ANOSMIA+data$DIARREA+data$EDO_MENTAL
data$num_sintomas<-as.numeric(data$num_sintomas)
data$comorb<-as.numeric(data$DMT2+data$OBESIDAD+data$ECV+data$EPOC+data$CANCER+data$ERC+data$INMUNOCOMPROMETIDOS+data$CIRROSIS)
data$comorb<-na.tools::na.replace(data$comorb, 0)
data$edad_4C[data$EDAD<50]<-0;data$edad_4C[data$EDAD>=50 & data$EDAD<60]<-2
data$edad_4C[data$EDAD>=60 & data$EDAD<70]<-4;data$edad_4C[data$EDAD>=70 & data$EDAD<80]<-6
data$edad_4C[data$EDAD>=80]<-7
data$RR_4C[data$FR<20]<-0;data$RR_4C[data$FR>=20 & data$FR<30]<-1;data$RR_4C[data$FR>=30]<-2
data$spo2_4C[data$SAT_O2<92]<-2;data$spo2_4C[data$SAT_O2>=92]<-0
data$BUN_4C[data$BUN_INICIAL*0.357<7]<-0;data$BUN_4C[data$BUN_INICIAL*0.357>=7 & data$BUN_INICIAL*0.357<15]<-1;data$BUN_4C[data$BUN_INICIAL*0.357>14]<-3
data$PCR_4C[data$PCR_INICIAL<50]<-0;data$PCR_4C[data$PCR_INICIAL>=50 & data$PCR_INICIAL<100]<-1;data$PCR_4C[data$PCR_INICIAL>=100]<-2
data$comorb_4C[data$comorb==1]<-1;data$comorb_4C[data$comorb==0]<-0
data$comorb_4C[data$comorb>1]<-2

data$score_4C<-data$edad_4C+ifelse(data$GENERO=="Male",1,0)+data$comorb_4C+data$RR_4C+data$spo2_4C+(2*data$EDO_MENTAL)+data$BUN_4C+data$PCR_4C

#### Vitamin D ####
f1<-data %>% dplyr::select(intubacion1, VitD, GENERO) %>% drop_na() %>%
  ggplot(aes(x=intubacion1, y=VitD, fill=intubacion1))+
  geom_boxplot()+stat_compare_means()+theme_classic()+scale_y_log10()+
  labs(fill="IVM")+ylab("Vitamin D")+facet_wrap(~GENERO)+xlab("")+scale_fill_jama()
f2<-data %>% dplyr::select(defuncion1, VitD, GENERO) %>% drop_na() %>%
  ggplot(aes(x=defuncion1, y=VitD, fill=defuncion1))+
  geom_boxplot()+stat_compare_means()+theme_classic()+scale_y_log10()+
  labs(fill="Lethality")+ylab("Vitamin D")+facet_wrap(~GENERO)+xlab("")+scale_fill_jama()
f3<-data %>% dplyr::select(intubacion1, VitD, imc_cat2) %>% drop_na() %>%
  ggplot(aes(x=intubacion1, y=VitD, fill=intubacion1))+
  geom_boxplot()+stat_compare_means()+theme_classic()+scale_y_log10()+
  labs(fill="IVM")+ylab("Vitamin D")+facet_wrap(~imc_cat2)+xlab("")+scale_fill_jama()
f4<-data %>% dplyr::select(defuncion1, VitD, imc_cat2) %>% drop_na() %>%
  ggplot(aes(x=defuncion1, y=VitD, fill=defuncion1))+
  geom_boxplot()+stat_compare_means()+theme_classic()+scale_y_log10()+
  labs(fill="Lethality")+ylab("Vitamin D")+facet_wrap(~imc_cat2)+xlab("")+scale_fill_jama()

fig1<-ggarrange(f1, f2, f3, f4, labels = LETTERS[1:4])

ggsave(file = "Figure1.jpg", 
       fig1,
       width = 15, 
       height = 10,
       units=c("in"),
       dpi = 300,
       limitsize = FALSE)


#### Comparaciones ####

tapply(data$EDAD, data$vitd20, mean,na.rm=T)
tapply(data$EDAD, data$vitd20, sd, na.rm=T)
t.test(data$EDAD~data$vitd20)

table(data$vitd20, data$GENERO)
table(data$vitd20, data$GENERO)%>%prop.table(1)
chisq.test(data$vitd20, data$GENERO)

table(data$vitd20, data$NIVEL_SOCIOECO==1)
table(data$vitd20, data$NIVEL_SOCIOECO==1)%>%prop.table(1)*100
chisq.test(data$vitd20, data$NIVEL_SOCIOECO==1)

table(data$vitd20, data$critical1)
table(data$vitd20, data$critical1)%>%prop.table(1)*100
chisq.test(data$vitd20, data$critical1)

table(data$vitd20, data$intubacion1)
table(data$vitd20, data$intubacion1)%>%prop.table(1)*100
chisq.test(data$vitd20, data$intubacion1)

table(data$vitd20, data$defuncion1)
table(data$vitd20, data$defuncion1)%>%prop.table(1)*100
chisq.test(data$vitd20, data$defuncion1)

table(data$vitd20, data$HIPERTENSION_ARTERIAL==1)
table(data$vitd20, data$HIPERTENSION_ARTERIAL==1)%>%prop.table(1)*100
chisq.test(data$vitd20, data$HIPERTENSION_ARTERIAL==1)

table(data$vitd20, data$DMT2)
table(data$vitd20, data$DMT2)%>%prop.table(1)
chisq.test(data$vitd20, data$DMT2)

tapply(data$TIEMPO_DX, data$vitd20, mean,na.rm=T)
tapply(data$TIEMPO_DX, data$vitd20, sd, na.rm=T)
t.test(data$TIEMPO_DX~data$vitd20)


table(data$vitd20, data$OBESIDAD)
table(data$vitd20, data$OBESIDAD)%>%prop.table(1)*100
chisq.test(data$vitd20, data$OBESIDAD)

table(data$vitd20, data$TABAQUISMO_ACTUAL>=2)
table(data$vitd20, data$TABAQUISMO_ACTUAL>=2)%>%prop.table(1)*100
chisq.test(data$vitd20, data$TABAQUISMO_ACTUAL>=1)

table(data$vitd20, data$ERC)
table(data$vitd20, data$ERC)%>%prop.table(1)*100
chisq.test(data$vitd20, data$ERC)

table(data$vitd20, data$ECV)
table(data$vitd20, data$ECV)%>%prop.table(1)*100
chisq.test(data$vitd20, data$ECV)

table(data$vitd20, data$CIRROSIS)
table(data$vitd20, data$CIRROSIS)%>%prop.table(1)*100
chisq.test(data$vitd20, data$CIRROSIS)

table(data$vitd20, data$NEFROPATIA)
table(data$vitd20, data$NEFROPATIA)%>%prop.table(1)*100
chisq.test(data$vitd20, data$NEFROPATIA)

tapply(data$IMC, data$vitd20, mean,na.rm=T)
tapply(data$IMC, data$vitd20, sd, na.rm=T)
t.test(data$IMC~data$vitd20)

tapply(data$SAT_O2, data$vitd20, mean,na.rm=T)
tapply(data$SAT_O2, data$vitd20, sd, na.rm=T)
t.test(data$SAT_O2~data$vitd20)

tapply(data$FR, data$vitd20, mean,na.rm=T)
tapply(data$FR, data$vitd20, sd, na.rm=T)
t.test(data$FR~data$vitd20)

tapply(data$FC, data$vitd20, mean,na.rm=T)
tapply(data$FC, data$vitd20, sd, na.rm=T)
t.test(data$FC~data$vitd20)

quantile(data[data$vitd20==1,]$TAS, na.rm = T)
quantile(data[data$vitd20==0,]$TAS, na.rm = T)
wilcox.test(data$TAS~data$vitd20)

quantile(data[data$vitd20==1,]$TAD, na.rm = T)
quantile(data[data$vitd20==0,]$TAD, na.rm = T)
wilcox.test(data$TAD~data$vitd20)

quantile(data[data$vitd20==1,]$PCR_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$PCR_INICIAL, na.rm = T)
wilcox.test(data$PCR_INICIAL~data$vitd20)

tapply(data$GLUCOSA_INICIAL, data$vitd20, mean,na.rm=T)
tapply(data$GLUCOSA_INICIAL, data$vitd20, sd, na.rm=T)
t.test(data$GLUCOSA_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$HbA1c, na.rm = T)
quantile(data[data$vitd20==0,]$HbA1c, na.rm = T)
wilcox.test(data$HbA1c~data$vitd20)

quantile(data[data$vitd20==1,]$TRIGLICERIDOS, na.rm = T)
quantile(data[data$vitd20==0,]$TRIGLICERIDOS, na.rm = T)
wilcox.test(data$TRIGLICERIDOS~data$vitd20)

tapply(data$C_HDL, data$vitd20, mean,na.rm=T)
tapply(data$C_HDL, data$vitd20, sd, na.rm=T)
t.test(data$C_HDL~data$vitd20)

tapply(data$C_LDL, data$vitd20, mean,na.rm=T)
tapply(data$C_LDL, data$vitd20, sd, na.rm=T)
t.test(data$C_LDL~data$vitd20)

tapply(data$HbA1c, data$vitd20, mean,na.rm=T)
tapply(data$HbA1c, data$vitd20, sd, na.rm=T)
t.test(data$HbA1c~data$vitd20)

tapply(data$COLESTEROL_T, data$vitd20, mean,na.rm=T)
tapply(data$COLESTEROL_T, data$vitd20, sd, na.rm=T)
t.test(data$COLESTEROL_T~data$vitd20)

tapply(data$LEUCOCITOS_TOTAL, data$vitd20, mean,na.rm=T)
tapply(data$LEUCOCITOS_TOTAL, data$vitd20, sd, na.rm=T)
t.test(data$LEUCOCITOS_TOTAL~data$vitd20)

tapply(data$NEUTRO_TOTAL, data$vitd20, mean,na.rm=T)
tapply(data$NEUTRO_TOTAL, data$vitd20, sd, na.rm=T)
t.test(data$NEUTRO_TOTAL~data$vitd20)

tapply(data$HEMOGLOBINA_INICIAL, data$vitd20, mean,na.rm=T)
tapply(data$HEMOGLOBINA_INICIAL, data$vitd20, sd, na.rm=T)
t.test(data$HEMOGLOBINA_INICIAL~data$vitd20)

tapply(data$PLAQUETAS_INICIAL, data$vitd20, mean,na.rm=T)
tapply(data$PLAQUETAS_INICIAL, data$vitd20, sd, na.rm=T)
t.test(data$PLAQUETAS_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$CREATININA, na.rm = T)
quantile(data[data$vitd20==0,]$CREATININA, na.rm = T)
wilcox.test(data$CREATININA~data$vitd20)

quantile(data[data$vitd20==1,]$FERRITINA, na.rm = T)
quantile(data[data$vitd20==0,]$FERRITINA, na.rm = T)
wilcox.test(data$FERRITINA~data$vitd20)

quantile(data[data$vitd20==1,]$D_DIMER, na.rm = T)
quantile(data[data$vitd20==0,]$D_DIMER, na.rm = T)
wilcox.test(data$D_DIMER~data$vitd20)

quantile(data[data$vitd20==1,]$PROTROMBINA_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$PROTROMBINA_INICIAL, na.rm = T)
wilcox.test(data$PROTROMBINA_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$FIBRINOGENO, na.rm = T)
quantile(data[data$vitd20==0,]$FIBRINOGENO, na.rm = T)
wilcox.test(data$FIBRINOGENO~data$vitd20)

tapply(data$BUN_INICIAL, data$vitd20, mean,na.rm=T)
tapply(data$BUN_INICIAL, data$vitd20, sd, na.rm=T)
t.test(data$BUN_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$AST_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$AST_INICIAL, na.rm = T)
wilcox.test(data$AST_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$ALT_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$ALT_INICIAL, na.rm = T)
wilcox.test(data$ALT_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$ALBUMINA_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$ALBUMINA_INICIAL, na.rm = T)
wilcox.test(data$ALBUMINA_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$DHL_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$DHL_INICIAL, na.rm = T)
wilcox.test(data$DHL_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$CK_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$CK_INICIAL, na.rm = T)
wilcox.test(data$CK_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$BD_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$BD_INICIAL, na.rm = T)
wilcox.test(data$BD_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$BI_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$BI_INICIAL, na.rm = T)
wilcox.test(data$BI_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$PROCALCITONINA_INICIAL, na.rm = T)
quantile(data[data$vitd20==0,]$PROCALCITONINA_INICIAL, na.rm = T)
wilcox.test(data$PROCALCITONINA_INICIAL~data$vitd20)

quantile(data[data$vitd20==1,]$num_sintomas, na.rm = T)
quantile(data[data$vitd20==0,]$num_sintomas, na.rm = T)
wilcox.test(data$num_sintomas~data$vitd20)

quantile(data[data$vitd20==1,]$TOTALCOMORBILIDAD, na.rm = T)
quantile(data[data$vitd20==0,]$TOTALCOMORBILIDAD, na.rm = T)
wilcox.test(data$TOTALCOMORBILIDAD~data$vitd20)

quantile(data[data$vitd20==1,]$DIAS_INTERNAMIENTO, na.rm = T)
quantile(data[data$vitd20==0,]$DIAS_INTERNAMIENTO, na.rm = T)
wilcox.test(data$DIAS_INTERNAMIENTO~data$vitd20)

quantile(data[data$vitd20==1,]$GROSOR_PROMEDIO, na.rm = T)
quantile(data[data$vitd20==0,]$GROSOR_PROMEDIO, na.rm = T)
wilcox.test(data$GROSOR_PROMEDIO~data$vitd20)

table(data$vitd20, data$GRASA_EXTRAPERICARDICA)
table(data$vitd20, data$GRASA_EXTRAPERICARDICA)%>%prop.table(1)*100
chisq.test(data$vitd20, data$GRASA_EXTRAPERICARDICA)

quantile(data[data$vitd20==1,]$GRASA_SUBCUTTORACICA, na.rm = T)
quantile(data[data$vitd20==0,]$GRASA_SUBCUTTORACICA, na.rm = T)
wilcox.test(data$GRASA_SUBCUTTORACICA~data$vitd20)

table(data$vitd20, data$TC_CONSOLIDATION, useNA="always")
table(data$vitd20, data$TC_CONSOLIDATION,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd20, data$TC_CONSOLIDATION)

table(data$vitd20, data$TC_GGO, useNA="always")
table(data$vitd20, data$TC_GGO,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd20, data$TC_GGO)

table(data$vitd20, data$TC_GGO.CONS, useNA="always")
table(data$vitd20, data$TC_GGO.CONS,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd20, data$TC_GGO.CONS)

table(data$vitd20, data$TC_GRADOAFECCION, useNA="always")
table(data$vitd20, data$TC_GRADOAFECCION,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd20, data$TC_GRADOAFECCION)

table(data$vitd20, data$TC_ESTEATOSIS, useNA="always")
table(data$vitd20, data$TC_ESTEATOSIS,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd20, data$TC_ESTEATOSIS)

quantile(data[data$vitd20==1,]$NEWS_PUNTOS, na.rm = T)
quantile(data[data$vitd20==0,]$NEWS_PUNTOS, na.rm = T)
wilcox.test(data$NEWS_PUNTOS~data$vitd20)

quantile(data[data$vitd20==1,]$QSOFA, na.rm = T)
quantile(data[data$vitd20==0,]$QSOFA, na.rm = T)
wilcox.test(data$QSOFA~data$vitd20)

quantile(data[data$vitd20==1,]$CURB65, na.rm = T)
quantile(data[data$vitd20==0,]$CURB65, na.rm = T)
wilcox.test(data$CURB65~data$vitd20)




####Comparaciones con VITD 12#####

tapply(data$EDAD, data$vitd12, mean,na.rm=T)
tapply(data$EDAD, data$vitd12, sd, na.rm=T)
t.test(data$EDAD~data$vitd12)

table(data$vitd12, data$GENERO)
table(data$vitd12, data$GENERO)%>%prop.table(1)
chisq.test(data$vitd12, data$GENERO)

table(data$vitd12, data$critical1)
table(data$vitd12, data$critical1)%>%prop.table(1)*100
chisq.test(data$vitd12, data$critical1)

table(data$vitd12, data$intubacion1)
table(data$vitd12, data$intubacion1)%>%prop.table(1)*100
chisq.test(data$vitd12, data$intubacion1)

table(data$vitd12, data$defuncion1)
table(data$vitd12, data$defuncion1)%>%prop.table(1)*100
chisq.test(data$vitd12, data$defuncion1)

table(data$vitd12, data$GENERO)
table(data$vitd12, data$GENERO)%>%prop.table(1)
chisq.test(data$vitd12, data$GENERO)

table(data$vitd12, data$NIVEL_SOCIOECO==1)
table(data$vitd12, data$NIVEL_SOCIOECO==1)%>%prop.table(1)*100
chisq.test(data$vitd12, data$NIVEL_SOCIOECO==1)

table(data$vitd12, data$HIPERTENSION_ARTERIAL==1)
table(data$vitd12, data$HIPERTENSION_ARTERIAL==1)%>%prop.table(1)*100
chisq.test(data$vitd12, data$HIPERTENSION_ARTERIAL==1)

table(data$vitd12, data$DMT2)
table(data$vitd12, data$DMT2)%>%prop.table(1)
chisq.test(data$vitd12, data$DMT2)

tapply(data$TIEMPO_DX, data$vitd12, mean,na.rm=T)
tapply(data$TIEMPO_DX, data$vitd12, sd, na.rm=T)
t.test(data$TIEMPO_DX~data$vitd12)

table(data$vitd12, data$OBESIDAD)
table(data$vitd12, data$OBESIDAD)%>%prop.table(1)*100
chisq.test(data$vitd12, data$OBESIDAD)

table(data$vitd12, data$TABAQUISMO_ACTUAL>=2)
table(data$vitd12, data$TABAQUISMO_ACTUAL>=2)%>%prop.table(1)*100
chisq.test(data$vitd12, data$TABAQUISMO_ACTUAL>=1)

table(data$vitd12, data$ERC)
table(data$vitd12, data$ERC)%>%prop.table(1)*100
chisq.test(data$vitd12, data$ERC)

table(data$vitd12, data$ECV)
table(data$vitd12, data$ECV)%>%prop.table(1)*100
chisq.test(data$vitd12, data$ECV)


table(data$vitd12, data$CIRROSIS)
table(data$vitd12, data$CIRROSIS)%>%prop.table(1)*100
chisq.test(data$vitd12, data$CIRROSIS)

table(data$vitd12, data$NEFROPATIA)
table(data$vitd12, data$NEFROPATIA)%>%prop.table(1)*100
chisq.test(data$vitd12, data$NEFROPATIA)

tapply(data$IMC, data$vitd12, mean,na.rm=T)
tapply(data$IMC, data$vitd12, sd, na.rm=T)
t.test(data$IMC~data$vitd12)

tapply(data$SAT_O2, data$vitd12, mean,na.rm=T)
tapply(data$SAT_O2, data$vitd12, sd, na.rm=T)
t.test(data$SAT_O2~data$vitd12)

tapply(data$FR, data$vitd12, mean,na.rm=T)
tapply(data$FR, data$vitd12, sd, na.rm=T)
t.test(data$FR~data$vitd12)

tapply(data$FC, data$vitd12, mean,na.rm=T)
tapply(data$FC, data$vitd12, sd, na.rm=T)
t.test(data$FC~data$vitd12)

quantile(data[data$vitd12==1,]$TAS, na.rm = T)
quantile(data[data$vitd12==0,]$TAS, na.rm = T)
wilcox.test(data$TAS~data$vitd12)

quantile(data[data$vitd12==1,]$TAD, na.rm = T)
quantile(data[data$vitd12==0,]$TAD, na.rm = T)
wilcox.test(data$TAD~data$vitd12)

quantile(data[data$vitd12==1,]$PCR_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$PCR_INICIAL, na.rm = T)
wilcox.test(data$PCR_INICIAL~data$vitd12)

tapply(data$GLUCOSA_INICIAL, data$vitd12, mean,na.rm=T)
tapply(data$GLUCOSA_INICIAL, data$vitd12, sd, na.rm=T)
t.test(data$GLUCOSA_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$HbA1c, na.rm = T)
quantile(data[data$vitd12==0,]$HbA1c, na.rm = T)
wilcox.test(data$HbA1c~data$vitd12)

quantile(data[data$vitd12==1,]$TRIGLICERIDOS, na.rm = T)
quantile(data[data$vitd12==0,]$TRIGLICERIDOS, na.rm = T)
wilcox.test(data$TRIGLICERIDOS~data$vitd12)

tapply(data$C_HDL, data$vitd12, mean,na.rm=T)
tapply(data$C_HDL, data$vitd12, sd, na.rm=T)
t.test(data$C_HDL~data$vitd12)

tapply(data$C_LDL, data$vitd12, mean,na.rm=T)
tapply(data$C_LDL, data$vitd12, sd, na.rm=T)
t.test(data$C_LDL~data$vitd12)

tapply(data$HbA1c, data$vitd12, mean,na.rm=T)
tapply(data$HbA1c, data$vitd12, sd, na.rm=T)
t.test(data$HbA1c~data$vitd12)

tapply(data$COLESTEROL_T, data$vitd12, mean,na.rm=T)
tapply(data$COLESTEROL_T, data$vitd12, sd, na.rm=T)
t.test(data$COLESTEROL_T~data$vitd12)

tapply(data$LEUCOCITOS_TOTAL, data$vitd12, mean,na.rm=T)
tapply(data$LEUCOCITOS_TOTAL, data$vitd12, sd, na.rm=T)
t.test(data$LEUCOCITOS_TOTAL~data$vitd12)

tapply(data$NEUTRO_TOTAL, data$vitd12, mean,na.rm=T)
tapply(data$NEUTRO_TOTAL, data$vitd12, sd, na.rm=T)
t.test(data$NEUTRO_TOTAL~data$vitd12)

tapply(data$HEMOGLOBINA_INICIAL, data$vitd12, mean,na.rm=T)
tapply(data$HEMOGLOBINA_INICIAL, data$vitd12, sd, na.rm=T)
t.test(data$HEMOGLOBINA_INICIAL~data$vitd12)

tapply(data$PLAQUETAS_INICIAL, data$vitd12, mean,na.rm=T)
tapply(data$PLAQUETAS_INICIAL, data$vitd12, sd, na.rm=T)
t.test(data$PLAQUETAS_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$CREATININA, na.rm = T)
quantile(data[data$vitd12==0,]$CREATININA, na.rm = T)
wilcox.test(data$CREATININA~data$vitd12)

quantile(data[data$vitd12==1,]$FERRITINA, na.rm = T)
quantile(data[data$vitd12==0,]$FERRITINA, na.rm = T)
wilcox.test(data$FERRITINA~data$vitd12)

quantile(data[data$vitd12==1,]$D_DIMER, na.rm = T)
quantile(data[data$vitd12==0,]$D_DIMER, na.rm = T)
wilcox.test(data$D_DIMER~data$vitd12)

quantile(data[data$vitd12==1,]$PROTROMBINA_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$PROTROMBINA_INICIAL, na.rm = T)
wilcox.test(data$PROTROMBINA_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$FIBRINOGENO, na.rm = T)
quantile(data[data$vitd12==0,]$FIBRINOGENO, na.rm = T)
wilcox.test(data$FIBRINOGENO~data$vitd12)

tapply(data$BUN_INICIAL, data$vitd12, mean,na.rm=T)
tapply(data$BUN_INICIAL, data$vitd12, sd, na.rm=T)
t.test(data$BUN_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$AST_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$AST_INICIAL, na.rm = T)
wilcox.test(data$AST_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$ALT_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$ALT_INICIAL, na.rm = T)
wilcox.test(data$ALT_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$ALBUMINA_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$ALBUMINA_INICIAL, na.rm = T)
wilcox.test(data$ALBUMINA_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$DHL_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$DHL_INICIAL, na.rm = T)
wilcox.test(data$DHL_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$CK_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$CK_INICIAL, na.rm = T)
wilcox.test(data$CK_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$PROCALCITONINA_INICIAL, na.rm = T)
quantile(data[data$vitd12==0,]$PROCALCITONINA_INICIAL, na.rm = T)
wilcox.test(data$PROCALCITONINA_INICIAL~data$vitd12)

quantile(data[data$vitd12==1,]$num_sintomas, na.rm = T)
quantile(data[data$vitd12==0,]$num_sintomas, na.rm = T)
wilcox.test(data$num_sintomas~data$vitd12)

quantile(data[data$vitd12==1,]$TOTALCOMORBILIDAD, na.rm = T)
quantile(data[data$vitd12==0,]$TOTALCOMORBILIDAD, na.rm = T)
wilcox.test(data$TOTALCOMORBILIDAD~data$vitd12)

quantile(data[data$vitd12==1,]$DIAS_INTERNAMIENTO, na.rm = T)
quantile(data[data$vitd12==0,]$DIAS_INTERNAMIENTO, na.rm = T)
wilcox.test(data$DIAS_INTERNAMIENTO~data$vitd12)

quantile(data[data$vitd12==1,]$GROSOR_PROMEDIO, na.rm = T)
quantile(data[data$vitd12==0,]$GROSOR_PROMEDIO, na.rm = T)
wilcox.test(data$GROSOR_PROMEDIO~data$vitd12)

table(data$vitd12, data$GRASA_EXTRAPERICARDICA)
table(data$vitd12, data$GRASA_EXTRAPERICARDICA)%>%prop.table(1)*100
chisq.test(data$vitd12, data$GRASA_EXTRAPERICARDICA)

quantile(data[data$vitd12==1,]$GRASA_SUBCUTTORACICA, na.rm = T)
quantile(data[data$vitd12==0,]$GRASA_SUBCUTTORACICA, na.rm = T)
wilcox.test(data$GRASA_SUBCUTTORACICA~data$vitd12)

table(data$vitd12, data$TC_CONSOLIDATION, useNA="always")
table(data$vitd12, data$TC_CONSOLIDATION,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd12, data$TC_CONSOLIDATION)

table(data$vitd12, data$TC_GGO, useNA="always")
table(data$vitd12, data$TC_GGO,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd12, data$TC_GGO)

table(data$vitd12, data$TC_GGO.CONS, useNA="always")
table(data$vitd12, data$TC_GGO.CONS,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd12, data$TC_GGO.CONS)

table(data$vitd12, data$TC_GRADOAFECCION, useNA="always")
table(data$vitd12, data$TC_GRADOAFECCION,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd12, data$TC_GRADOAFECCION)

table(data$vitd12, data$TC_ESTEATOSIS, useNA="always")
table(data$vitd12, data$TC_ESTEATOSIS,useNA="always")%>%prop.table(1)*100
chisq.test(data$vitd12, data$TC_ESTEATOSIS)

quantile(data[data$vitd12==1,]$NEWS_PUNTOS, na.rm = T)
quantile(data[data$vitd12==0,]$NEWS_PUNTOS, na.rm = T)
wilcox.test(data$NEWS_PUNTOS~data$vitd12)

quantile(data[data$vitd12==1,]$QSOFA, na.rm = T)
quantile(data[data$vitd12==0,]$QSOFA, na.rm = T)
wilcox.test(data$QSOFA~data$vitd12)

quantile(data[data$vitd12==1,]$CURB65, na.rm = T)
quantile(data[data$vitd12==0,]$CURB65, na.rm = T)
wilcox.test(data$CURB65~data$vitd12)

#### Models ####

m1<-coxph(Surv(fu_time, defuncion)~VitD, data=data)
summary(m1)

m2<-coxph(Surv(fu_time, defuncion)~VitD+EDAD+GENERO+IMC+PCR_INICIAL, data=data)
summary(m2)

data$log_epi<-log(data$GROSOR_PROMEDIO+1)
m3<-coxph(Surv(fu_time, defuncion)~VitD+IMC+EDAD+GENERO+log(PCR_INICIAL)+ERC+DMT2+log_epi+log(D_DIMER)+EDAD+SAT_O2+log(data$TROPONINA_ULTRASENSIBLE), data=data)
summary(m3)

m3.i<-coxph(Surv(fu_time, defuncion)~VitD*IMC+EDAD+GENERO+log(PCR_INICIAL)+ERC+DMT2+log_epi+log(D_DIMER)+EDAD+SAT_O2+log(data$TROPONINA_ULTRASENSIBLE), data=data)
summary(m3.i)


m4<-coxph(Surv(fu_time, defuncion)~(VitD<20), data=data)
summary(m4)
m4<-coxph(Surv(fu_time, defuncion)~(VitD<20)*imc_cat2+EDAD+GENERO+IMC+PCR_INICIAL+ERC+GROSOR_PROMEDIO+DMT2+log_epi, data=data)
summary(m4)

m4<-coxph(Surv(fu_time, defuncion)~(VitD<=12), data=data)
summary(m4)
m5<-coxph(Surv(fu_time, defuncion)~(VitD<=12)+EDAD+GENERO+log(PCR_INICIAL)+ERC+DMT2+log_epi+log(data$D_DIMER)+data$SAT_O2, data=data)
summary(m5)


l1<-glm(INTUBACION~VitD+EDAD+GENERO+IMC+DMT2+ERC+PCR_INICIAL, data=data, family="binomial")
summary(l1)
exp(cbind(coef(l1), confint(l1)))

l1<-glm(data$critical~VitD+EDAD+GENERO+IMC+DMT2+PCR_INICIAL+SAT_O2+D_DIMER+GROSOR_PROMEDIO, data=data, family="binomial")
summary(l1)
exp(cbind(coef(l1), confint(l1)))

l1<-glm(data$critical~vitd20+EDAD+GENERO+IMC+DMT2+PCR_INICIAL+SAT_O2, data=data, family="binomial")
summary(l1)
exp(cbind(coef(l1), confint(l1)))


### Vitamin D determinants ###
data$log_epi<-log(data$GROSOR_PROMEDIO+1)
l1<-lm(log(VitD)~GENERO+log(data$TROPONINA_ULTRASENSIBLE)+log(data$D_DIMER)+log(PCR_INICIAL)+log_epi+factor(data$DMT2),data=data)
car::vif(l1)
cbind(coef(l1),confint(l1), summary(l1)$coefficients[,4])
summary(l1)
hist(l1$residuals)
l1<-glm((data$VitD<20)~GENERO+log(data$PCR_INICIAL)+log_epi+log(data$GLUCOSA_INICIAL)+log(data$D_DIMER)+EDAD+data$SAT_O2, data=data, family="binomial")
summary(l1)
exp(cbind(coef(l1), confint(l1)))

#### Propensity score matching ####

data1<- data %>% dplyr::select(defuncion, critical,intubacion1,fu_time, VitD, EDAD, GENERO, DMT2, ERC, INICIODESINTOMAS, IMC, PCR_INICIAL, D_DIMER, SAT_O2, GROSOR_PROMEDIO) %>% drop_na()
set.seed(123)
match<-matchit(defuncion~EDAD+GENERO+IMC+DMT2+ERC+INICIODESINTOMAS+PCR_INICIAL+D_DIMER+SAT_O2+GROSOR_PROMEDIO, data=data1, method = "nearest", caliper=0.25)
X2<-c(match$match.matrix); X2<-as.numeric(X2)
X1<-rownames(match$match.matrix); X1<-as.numeric(X1)
pairs<-data.frame(X1, X2)
pairs<-na.omit(pairs)
id<-as.numeric(c(pairs$X1, pairs$X2))
pairs2<-data.frame(id, rep(1:nrow(pairs),2))
names(pairs2)<-c("id", "match")
pairs2$id<-as.numeric(pairs2$id)
pairs2<- pairs2 %>%arrange(id)
base2<-match.data(match)
base2<-base2[pairs2$id %in% rownames(base2),]
base2$id<-pairs2$id;base2$match<-pairs2$match
table(base2$defuncion)

m1<-coxphw::coxphw(Surv(fu_time, defuncion)~VitD+frailty(match), data=base2)
summary(m1)

m1<-coxphw::coxphw(Surv(fu_time, defuncion)~(VitD<=12)+frailty(match), data=data)
summary(m1)

m1<-glmer(critical~VitD+(1|match), data=base2, family="binomial")
summary(m1)

#### Medicación ####
library(mediation)
data$log_VitD<-log(data$VitD)
data$log_DD<-log(data$D_DIMER)
data$log_tropo<-log(data$TROPONINA_ULTRASENSIBLE)
data$log_pcr<-log(data$PCR_INICIAL)
data$log_epi<-log(data$GROSOR_PROMEDIO+1)
data$log_sat<-log(data$SAT_O2)

##D-dimer

m1<-lm(log_DD~log_VitD+log_epi+IMC+EDAD+GENERO, data=data)
m2<-glm(defuncion~log_DD+log_epi+IMC+log_VitD+EDAD+GENERO, data=data, family="binomial")

set.seed(123)
med<-mediate(model.m = m1, model.y = m2, treat = "log_VitD", mediator = "log_DD", boot = T, boot.ci.type = "bca",sims = 1000)
summary(med)


#Troponins
m1<-lm(log_tropo~log_VitD+log_epi+IMC+EDAD+GENERO, data=data)
m2<-glm(defuncion~log_tropo+log_epi+IMC+log_VitD+EDAD+GENERO, data=data, family="binomial")
summary(m1)


set.seed(123)
med<-mediate(model.m = m1, model.y = m2, treat = "log_VitD", mediator = "log_tropo", boot = T, boot.ci.type = "bca",sims = 1000)
summary(med)

#Oxygen saturation
m1<-lm(log_sat~log_VitD+EDAD+GENERO+IMC+log_epi, data=data)
m2<-glm(defuncion~log_sat+log_epi+IMC+log_VitD+EDAD+GENERO, data=data, family="binomial")

set.seed(123)
med<-mediate(model.m = m1, model.y = m2, treat = "log_VitD", mediator = "log_sat", boot = T, boot.ci.type = "bca",sims = 1000)
summary(med)


#### Simulación ####
Sim1<-coxsimLinear(m3, b = "VitD", Xj = seq(0, 80),  Xl = seq(20,20), qi = "Hazard Ratio",spin = FALSE, ci = 0.95)
f2a<-simGG(Sim1)+xlab("Vitamin D levels (ng/mL)")
Sim2<-coxsimLinear(m3, b = "VitD", Xj = seq(0, 80),  Xl = seq(12,12), qi = "Hazard Ratio",spin = FALSE, ci = 0.95)
f2b<-simGG(Sim2)+xlab("Vitamin D levels (ng/mL)")

f2<-ggarrange(f2a, f2b, labels=LETTERS[1:2])

ggsave(file = "Figure2.jpg", 
       f2,
       width = 30, 
       height = 15,
       units=c("cm"),
       dpi = 600,
       limitsize = FALSE)

