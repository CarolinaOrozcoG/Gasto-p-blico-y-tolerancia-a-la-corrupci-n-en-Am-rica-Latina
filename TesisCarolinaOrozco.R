####Paquetería####
library(tidyverse) 
library(extrafont) 
library(readxl)  
library(stargazer) 
library(lme4)  
library(MuMIn) 
library(MASS)
library(margins)
library(ggdag) 
library(wesanderson) 
library(sensemakr) 
library(sem) 
library(dagitty)  
library(parameters) 
library(dplyr)
library(lavaan)
library(lavaanPlot)
library(piecewiseSEM)
library(pals)
library(ordinal) 
library(car)
library(estimatr)
library(texreg)
library(patchwork)
library(broom.mixed)

font_import()

loadfonts(device ="win", quiet =TRUE)
rm(list=ls())

setwd("C:\\Users\\carol\\OneDrive - correounivalle.edu.co\\MAESTRIA GAP\\Seminario\\Datos\\TESIS")


####Cargar la base ####

BaseCompleta <- read_xlsx("Basecompleta.xlsx")

#### Se excluyen de la muestra los países no democráticos

BaseFiltrada <- BaseCompleta %>% 
  filter(Pais!="VEN" & Pais!="NIC")

#### Se crea la variable año como factor para su inclusión en los modelos

BaseFiltrada <- BaseFiltrada %>% 
  mutate(año = as.factor(year))

#### Recodificando la variable de resultado

BaseFiltrada <- BaseFiltrada %>% 
  mutate(Toledecod= (ifelse(tolerancia ==1,4,
                            ifelse(tolerancia==2,3,
                                   ifelse(tolerancia==3,2,
                                          ifelse(tolerancia==4,1,0))))))
BaseFiltrada <- BaseFiltrada %>% 
  mutate(Toledecod = as.factor(Toledecod)) %>% 
  mutate(Toledecod= as.ordered(Toledecod))

#### Generación de DAGs ####

## Crear objeto con las relaciones entre las variables

simple_dag <- dagify(
  TC ~ G+ED+NE+I+P+SE+PIB,       
  SE ~ CE ,          
  G ~ PIB+CE+I,
  CE ~ ED,
  NE~ L,
  PIB~ L,
  I ~ M,
  P~ M,
latent = c("L","M"),
exposure = "G",          
  outcome = "TC"            
)

dag <- ggdag_status(simple_dag, layout= "nicely")

#Cambiar títulos de las variables

dag$data$status <-  case_when(is.na(dag$data$status) ~ "Covariable",
                              dag$data$status == "latent" ~ "Latente",
                              dag$data$status == "exposure" ~ "Tratamiento",
                              dag$data$status == "outcome" ~ "Resultado")

#Cambiar colores del gráfico

dag <- dag +
  theme_dag()+
  labs(title= "Gráfico Acíclico Direccional",
       caption = "Fuente: elaboración propia")+
  scale_color_manual(values= c("cadetblue","mediumpurple","seagreen","deeppink2")
  )
dag
dag

#Guardar la imagen

ggsave(filename = "dag.png", dpi = 500)

# Obtener trayectorias de puerta trasera

paths(simple_dag)

# Combinaciones que bloquean las trayectorias de puerta de trasera

adjustmentSets(simple_dag)

# Visualizar las alternativas de bloqueo

ggdag_adjustment_set(simple_dag, shadow = TRUE)+ # zoom y pantalla completa
  theme_dag()+
  labs(title= "Alternativas de Bloqueo de Trayectorias de Puerta Trasera",
       caption = "Fuente: elaboración propia")+
  scale_color_manual(values= c("cadetblue","deeppink2")) +
  theme(plot.title=element_text(size=15, face='bold', color='black', hjust = 0.5))

# Exportar imagen con alternativas 

ggsave(filename = "spt.png", dpi = 500)

#### Modelo Logístico Ordinal ####

ordered <- clmm(Toledecod ~ gasto2+PosPol+Cuadrado+año+Rule_BM+SEcon+
                  log(valor.pib)+(1|Countryear), data = BaseFiltrada)

summary(ordered)

ordered1 <- clmm(Toledecod ~ gasto2+PosPol+Cuadrado+año+
                   crec.pib+log(valor.pib)+(1 |Countryear), data = BaseFiltrada)

summary(ordered1)

#Intervalos de confianza

ci <- confint(ordered)
ci1 <- confint(ordered1)

#Calculo de probabilidades

coef <- ordered$coefficients
coef1 <- ordered1$coefficients

prob <- (exp(coef)-1)*100
prob1 <- (exp(coef1)-1)*100

####Modelos lineales ####

lineal1 <- lmer(toledummy ~ gasto2+Rule_BM+PosPol+Cuadrado+año+
                  LogPib+SEcon+(1|Countryear), data = BaseFiltrada)
summary(lineal1)

lineal <- lmer(toledummy ~ gasto2+PosPol+Cuadrado+año+
                 LogPib+crec.pib+(1 |Countryear), data = BaseFiltrada)
summary(lineal)

####Con errores robustos####

robusto1 <- glmmPQL(toledummy ~ gasto2+PosPol+Cuadrado+
                      as.factor(year)+LogPib+crec.pib
                    , data = BaseFiltrada,random = ~ 1|Countryear,
                    family= gaussian )

summary(robusto1)

robusto2 <- glmmPQL(toledummy ~ gasto2+Rule_BM+PosPol+Cuadrado+
                      año+LogPib+
                      SEcon, data = BaseFiltrada,random = ~ 1|Countryear,
                    family= gaussian )

summary(robusto2)

#### Análisis de sensibilidad ####

prueba<- lm(toledummy ~ gasto2+Rule_BM+PosPol+Cuadrado+año+
              LogPib+SEcon, data = BaseFiltrada)
summary(prueba)

sensibilidad <- sensemakr(model= prueba, treatment= "gasto2", 
                          benchmark_covariates= "SEcon", kd=1:3,
                          ky=1:3, q=1, alpha=0.05, reduce= TRUE)

summary(sensibilidad)

#Para tabla en Latex

ovb_minimal_reporting(sensibilidad, format = "latex")

plot(sensibilidad)
plot(sensibilidad, sensitivity.of="t-value")

prueba1<- lm(toledummy ~ gasto2+PosPol+Cuadrado+año+
               LogPib+crec.pib, data = BaseFiltrada)

summary(prueba1)

sensibilidad1 <- sensemakr(model= prueba1, treatment= "gasto2", 
                           benchmark_covariates= "PosPol", kd=1:3,
                           ky=1:3, q=1, alpha=0.05, reduce= TRUE)

ovb_minimal_reporting(sensibilidad1, format = "latex")

plot(sensibilidad1)
plot(sensibilidad1, sensitivity.of="t-value")


#### Gráficos de coeficientes ####

coeficientes1 <- jtools::plot_coefs(lineal,lineal1,
                   coefs = c("Gasto Público" = "gasto2"),
                   ci_level = .95, inner_ci_level = .90,
                   model.names = c("Modelo 1 (I,CE,PIB)",
                                   "Modelo 2 (ED,I,PIB,SE)"),
                   colors = wes_palette("BottleRocket1"),
                   legend.title = "Modelo") +
  labs(x = "Estimación", y = "", 
       title= "Efecto del gasto público sobre la tolerancia a la corrupción")

coeficientes2 <- jtools::plot_coefs(ordered,ordered1,
                                    coefs = c("Gasto Público" = "gasto2"),
                                    ci_level = .99, inner_ci_level = .95,
                                    model.names = c("Modelo 1 (I,CE,PIB)",
                                                    "Modelo 2 (ED,I,PIB,SE)"),
                                    colors = wes_palette("BottleRocket1"),
                                    legend.title = "Modelo") +
  labs(x = "Estimación", y = "", 
       title= "Efecto del gasto público sobre la tolerancia a la corrupción")

#ggsave(filename="modelos.png",dpi=300)

#### Modelo de Ecuaciones estructurales ####

Basesita <- BaseFiltrada %>% 
  select_("toledummy","Ajustado", "gasto2", "recaudo","LogPib","crec.pib")

Basesita <- Basesita %>% 
  na.omit(toledummy)
  
Model <- 'toledummy ~ gasto2 + recaudo + Ajustado
          recaudo ~ gasto2
          Ajustado ~ gasto2'

model1.fit <- sem(Model, data = Basesita) 

summary(model1.fit, rsq = TRUE, fit.measures = TRUE, standardized = TRUE) 

plotLA <- lavaanPlot(name = "model1", model1.fit, coefs = TRUE) 

#### Debido al uso de los decimales que por defecto hace la función,
# se añade de manera manual en la siguiente linea de código el valor de
# la relación entre las variables respetando los resultados de la línea 245.

plotLA$x[[1]] <- " digraph model1 { \n graph [ overlap = true, fontsize = 10 ]
\n node [ shape = box ] \n node [shape = box] \n GP; RI; PC; TC \n node 
[shape = oval] \n  \n \n edge [ color = black ] \n GP->TC [label = \"-0.001\"] 
RI->TC [label = \"-0.004***\"] PC->TC [label = \"0.06***\"] GP->RI [label = \"0.54***\"]
GP->PC [label = \"-0.03***\"]  \n}"

plotLA

#### Test de submuestras ####

#### Se crea una función que aplique el modelo a cada submuestra

f <- function(base){
  p <- clmm(Toledecod ~ gasto2+PosPol+Cuadrado+año+Rule_BM+SEcon+
              log(valor.pib)+(1|Countryear), data = base)
  q <- p$coefficients[4]
  return(q)
}

##### Se crea una función que extraiga la estimación del gasto y su p-valor.

f1 <- function(base){
  p <- clmm(Toledecod ~ gasto2+PosPol+Cuadrado+año+Rule_BM+SEcon+
              log(valor.pib)+(1|Countryear), data = base)
  q <- summary(p)
  q <- q$coefficients[4,4]
  return(q)
}

#### Se crean 5000 submuestras de 1000 observaciones 

muestras <- list()
n <- 5000
set.seed(1)
seed <- sample(1:1000000,size = n, replace= F)

for(i in 1:n){
  set.seed(seed[i])
  x <- sample(x=1:nrow(BaseFiltrada),size=1000, replace= F)
  muestras[[i]] <- BaseFiltrada[x,]
}

#### Se aplican las funciones a cada muestra

coeficientes <- c()
pvalor <- c()

for (i in 1:n) {
  coeficientes[i] <- f(muestras[[i]])
  pvalor[i] <- f1(muestras[[i]])
}

## Igual para el Modelo 2

f2 <- function(base){
  p <- clmm(Toledecod ~ gasto2+PosPol+Cuadrado+año+
              crec.pib+log(valor.pib)+(1 |Countryear), data = base)
  q <- p$coefficients[4]
  return(q)
}

f3 <- function(base){
  p <- clmm(Toledecod ~ gasto2+PosPol+Cuadrado+año+
              crec.pib+log(valor.pib)+(1 |Countryear), data = base)
  q <- summary(p)
  q <- q$coefficients[4,4]
  return(q)
}

coeficientes1 <- c()
pvalor1 <- c()

for (i in 1:n) {
  coeficientes1[i] <- f2(muestras[[i]])
  pvalor1[i] <- f3(muestras[[i]])
}

#### Se crea tabla con los resultados

df <- tibble(pvaloresModel1 = pvalor,
             pvaloresModel2 = pvalor1,
             coeficientesModel1 = coeficientes,
             coeficientesModel2 = coeficientes1)

df <- df %>% 
  mutate(negativos1= ifelse(coeficientesModel1<0,1,0)) %>% 
  mutate(negativos2= ifelse(coeficientesModel2<0,1,0))

df1 <- df %>% 
  filter(!is.na(pvaloresModel1))

df2 <- df %>% 
  filter(!is.na(pvaloresModel2))

df1 <- df1 %>% 
  mutate(sign1= ifelse(pvaloresModel1<=0.05,"<0.05",">0.05"))

df2 <- df2 %>% 
  mutate(sign2= ifelse(pvaloresModel2<=0.05,"<0.05",">0.05"))

#### Se grafican los resultados

Model1 <- df1 %>% ggplot(mapping = aes(x= coeficientesModel1, fill= as.factor(sign1)))+
  geom_density(alpha = 0.5)+
  scale_y_continuous(limits = c(0,75))+
  scale_x_continuous(breaks = c(-0.06,-0.04,-0.02,0))+
  geom_vline(aes(xintercept=0),
             color="skyblue3", linetype="dashed", linewidth=1)+
  labs(title="Curva de densidad submuestras modelo 1",x="Coeficientes", 
       y = "Densidad", fill= "P-valor")+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(size= 11, hjust = 0.5),
        axis.title = element_text(size = 9),
        legend.title = element_text(size=9)) 

Model2 <- df2 %>% ggplot(mapping = aes(x= coeficientesModel2, fill = as.factor(sign2)))+
  geom_density(alpha = 0.5)+
  scale_y_continuous(limits = c(0,75))+
  geom_vline(aes(xintercept=0),
             color="skyblue3", linetype="dashed", linewidth=1)+
  labs(title="Curva de densidad submuestras modelo 2",x="Coeficientes", y =NULL,
       fill= "P-valor")+
  theme_bw() + 
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(size= 11, hjust = 0.5),
        axis.title = element_text(size=9),
        legend.title = element_text(size=9)) 

Model1+Model2

ggsave(filename="placebo.png",
       width = 7,
       height = 3,
       dpi=300)

writexl::write_xlsx(df, path = "placebo.xlsx")
