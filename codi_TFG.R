library(readxl)
library(dplyr)
library(ggplot2)
library(car)
library(minpack.lm)
library(nls2)
library(segmented)
library(drc)
library(devtools)
library(aomisc)
library(openxlsx)
library(Matrix)

datos <- read_excel("C:/Users/gerar/OneDrive/Documentos/data_lactate.xlsx")
datos

lactato_reposo <- datos$lactate[1]
lactato_reposo

datos_sin_reposo <- datos[-1,]
datos_sin_reposo
datos$log_power <-log(datos$power)
datos$log_lactate <-log(datos$lactate)

ggplot(datos, aes(x = power, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto")

LT <- segmented(lm(lactate ~ power, data=datos), npsi = 1)
LT

ggplot(datos, aes(x = power, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = 292.4, color = "red") +  
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT")


LT12 <- segmented(lm(lactate ~ power, data=datos), npsi = 2)
LT12


ggplot(datos, aes(x = power, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = 183.1 , color = "red") + 
  geom_vline(xintercept = 319.4 , color = "green") + 
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT1 y el LT2")

mod_lin <- lm(power ~ lactate, data = datos_sin_reposo)

ggplot(datos_sin_reposo, aes(x = lactate, y = power)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_abline(intercept = coef(mod_lin)[1], slope = coef(mod_lin)[2], col = "red") +
  labs(x = "Lactato", y = "Potencia", title = "Regresión lineal ajustada a la potencia sobre el lactato")


modelo_mm <- nls(power ~ SSmicmen(lactate, a, b), 
                 start = list(a = max(datos_sin_reposo$power), b = 2),
                 data = datos_sin_reposo)
modelo_mm
coef(modelo_mm)

confint(modelo_mm)

Power_mm <- function(lactate) {
  return((452.90683 * lactate) / (1.148256 + lactate))
}

plot(datos_sin_reposo$lactate, datos_sin_reposo$power, pch = 16, col = "blue", xlab = "lactate", ylab = "power")
lines(datos_sin_reposo$lactate, predict(modelo_mm), col = "red", lwd = 2)
legend("bottomright", legend = c("Datos", "Ajuste del modelo Michaelis-Menten"), col = c("blue", "red"), lty = 1, lwd = 2)


lactate_mm <- function(Power) {
  return((1.148256 * Power) / (452.90683 - Power))
}

power_values <- seq(100, 400, by = 0.05)


lactate_values_mm <- sapply(power_values, lactate_mm)

valores_mm <- data.frame(power = power_values, lactate = lactate_values_mm)


ggplot(datos_sin_reposo, aes(x = power, y = lactate),col = "blue") +
  geom_point() + 
  geom_line(data = valores_mm, aes(x = power, y = lactate), color = "red") +  
  labs(x = "Potencia", y = "Lactato")+theme_minimal()


modelo_asintotico <- drm(datos_sin_reposo$power ~ datos_sin_reposo$lactate, fct = DRC.asymReg())
coef(modelo_asintotico)

# Definir función
Power_asintotico <- function(x) {
  367.0111684  - (353.6983618) * exp(-0.8132733 * x)
}

confint(modelo_asintotico)


plot(datos_sin_reposo$lactate, datos_sin_reposo$power, pch = 16, col = "blue", xlab = "lactate", ylab = "power")
lines(datos_sin_reposo$lactate, predict(modelo_asintotico), col = "green", lwd = 2)
legend("bottomright", legend = c("Datos", "Ajuste del modelo asintótico"), col = c("blue", "green"), lty = 1, lwd = 2)


lactate_as <- function(Power) {
  return(-log((367.0111684-Power)/353.6983618)/0.8132733)
}

power_values <- seq(100, 400, by = 0.05)

lactate_values_as <- sapply(power_values, lactate_as)

valores_as <- data.frame(power = power_values, lactate = lactate_values_as)


ggplot(datos_sin_reposo, aes(x = power, y = lactate),col = "blue") +
  geom_point() + 
  geom_line(data = valores_as, aes(x = power, y = lactate), color = "green") +  
  labs(x = "Potencia", y = "Lactato")+theme_minimal()


newfit_model <- function(x, rest, p1, p2) {
  ifelse(x <= p1, rest, rest * exp((x - p1) / p2))
}

nls_formula <- lactate ~ newfit_model(power, rest, p1, p2)

modelo_segmentado<- nls(nls_formula, data = data.frame(power=datos_sin_reposo$power,lactate=datos_sin_reposo$lactate), start = list(rest = lactato_reposo, p1 = 150, p2 = 100))
params <- coef(modelo_segmentado)
params

lactato_reposo_opt <- params["rest"]
p1_opt <- params["p1"]
p2_opt <- params["p2"]

print(paste("rest:",lactato_reposo_opt))
print(paste("p1:", p1_opt))
print(paste("p2:", p2_opt))


confint(modelo_segmentado)


lactate_segmentado <- function(Power) {
  ifelse(Power <= p1_opt, lactato_reposo_opt, lactato_reposo_opt * exp((Power - p1_opt) / p2_opt))
}

power_values <- seq(100, 400, by = 0.05)

lactate_values_segmentado <- sapply(power_values, lactate_segmentado)

valores_segmentado <- data.frame(power = power_values, lactate = lactate_values_segmentado)


ggplot(datos_sin_reposo, aes(x = power, y = lactate),col = "blue") +
  geom_point() + 
  geom_line(data = valores_segmentado, aes(x = power, y = lactate), color = "tan2") +  
  labs(x = "Potencia", y = "Lactato")+theme_minimal()


ggplot() +
  geom_line(data = valores_mm, aes(x = power, y = lactate, color = "Ajuste del modelo Michaelis-Menten"), size = 1) +  
  geom_line(data = valores_as, aes(x = power, y = lactate, color = "Ajuste del modelo asintótico"), size = 1) +  
  geom_line(data = valores_segmentado, aes(x = power, y = lactate, color = "Ajuste del modelo segmentado"), size = 1) +  
  geom_point(data = datos_sin_reposo, aes(x = power, y = lactate, color = "Datos"), size = 2) + 
  labs(x = "Potencia", y = "Lactato", title = "Comparación de Modelos de predicción de Lactato") +
  scale_color_manual(
    name = "Modelos",
    values = c("Ajuste del modelo Michaelis-Menten" = "red3", 
               "Ajuste del modelo asintótico" = "green", 
               "Ajuste del modelo segmentado" = "tan2",
               "Datos" = "blue"),
    labels = c("Ajuste del modelo asintótico ","Ajuste del modelo Michaelis-Menten", "Ajuste del modelo segmentado" ,"Datos")
  ) +
  theme_minimal() +
  theme(legend.position = c(0.3, 0.75)) + # Coordenadas para la esquina superior izquierda
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1)))

power_values_2 <- seq(100, 375, by = 25)
power_values_2
power_values_as <- seq(100, 350, by = 25)
power_values_as

predicciones_mm <- lactate_mm(power_values_2)
predicciones_asintotico <- lactate_as(power_values_as)
predicciones_segmentado <- lactate_segmentado(power_values_2)

error_cuadratico_mm <- mean((datos_sin_reposo$lactate - predicciones_mm)^2)
error_cuadratico_mm
error_cuadratico_asintotico <- mean(( datos_sin_reposo$lactate- predicciones_asintotico)^2)
error_cuadratico_asintotico
error_cuadratico_segmentado <- mean((datos_sin_reposo$lactate - predicciones_segmentado)^2)
error_cuadratico_segmentado



AIC_mm <- AIC(modelo_mm)
AIC_mm
BIC_mm <- BIC(modelo_mm)
BIC_mm

AIC_asintotico <- AIC(modelo_asintotico)
AIC_asintotico
BIC_asintotico <- BIC(modelo_asintotico)
BIC_asintotico

AIC_segmentado <- AIC(modelo_segmentado)
AIC_segmentado
BIC_segmentado <- BIC(modelo_segmentado)
BIC_segmentado


resultados <- data.frame(Modelo = c("modelo_mm", "modelo_asintotico"),
                         ECM = c(error_cuadratico_mm, error_cuadratico_asintotico),
                         AIC = c(AIC_mm, AIC_asintotico),
                         BIC = c(BIC_mm, BIC_asintotico))
print(resultados)

calcular_VO2_max <- function(potencia,masa){
  vo2_max <- (10.8 * potencia/masa) + 7
  return(vo2_max)
}

calcular_VO2_max_2 <- function(potencia_entre_masa){
  vo2_max <- (10.8 * potencia_entre_masa) + 7
  return(vo2_max)
}


P_max <- 452.91
masa <- 75

VO_2_max <- calcular_VO2_max(P_max,masa)
VO_2_max


calcular_VLa_max <- function(max_lactato,lactato_reposo, tiempo_ejercicio, tiempo_atp){
  vla_max <- (max_lactato-lactato_reposo)/(tiempo_ejercicio- tiempo_atp)
  return(vla_max)
}


max_lactato <- 6.20
lactato_reposo <- 0.6
tiempo_ejercicio <- 26
tiempo_atp <- 3

VLa_max <- calcular_VLa_max(max_lactato,lactato_reposo, tiempo_ejercicio, tiempo_atp)
VLa_max


calcular_C <- function(METs, masa){
  C <- 0.0175*METs*masa
  return(C)
}

C <- calcular_C(18,masa)
C

calcular_endurance_perfomance <- function(LT2, VO2_max, C){
  return(LT2*(VO2_max/C))
}


LT2_porcentaje <- 70.52
VO2_max <- 74.43
endurance_performance <- calcular_endurance_perfomance(LT2_porcentaje, VO2_max, C)
endurance_performance


datos2 <- read_excel("C:/Users/gerar/OneDrive/Documentos/DATA_PMEA.xlsx")
datos2


datos_individuo_1 <- subset(datos2, id == 1)
datos_individuo_1
datos_individuo_1_sin_reposo <- datos_individuo_1[-1,]
datos_individuo_1_sin_reposo
datos_individuo_2 <- subset(datos2, id == 2)
datos_individuo_2_sin_reposo <- datos_individuo_2[-1,]
datos_individuo_3 <- subset(datos2, id == 3)
datos_individuo_3_sin_reposo <- datos_individuo_3[-1,]
datos_individuo_4 <- subset(datos2, id == 4)
datos_individuo_4_sin_reposo <- datos_individuo_4[-1,]
datos_individuo_5 <- subset(datos2, id == 7)
datos_individuo_5_sin_reposo <- datos_individuo_5[-1,]
datos_individuo_6 <- subset(datos2, id == 9)
datos_individuo_6_sin_reposo <- datos_individuo_6[-1,]


ggplot() +
  geom_point(data = datos_individuo_1, aes(x = watts_kg, y = lactate, color = "Individuo 1")) +
  geom_line(data = datos_individuo_1, aes(x = watts_kg, y = lactate, color = "Individuo 1")) +
  geom_point(data = datos_individuo_2, aes(x = watts_kg, y = lactate, color = "Individuo 2")) +
  geom_line(data = datos_individuo_2, aes(x = watts_kg, y = lactate, color = "Individuo 2")) +
  geom_point(data = datos_individuo_3, aes(x = watts_kg, y = lactate, color = "Individuo 3")) +
  geom_line(data = datos_individuo_3, aes(x = watts_kg, y = lactate, color = "Individuo 3")) +
  geom_point(data = datos_individuo_4, aes(x = watts_kg, y = lactate, color = "Individuo 4")) +
  geom_line(data = datos_individuo_4, aes(x = watts_kg, y = lactate, color = "Individuo 4")) +
  geom_point(data = datos_individuo_5, aes(x = watts_kg, y = lactate, color = "Individuo 5")) +
  geom_line(data = datos_individuo_5, aes(x = watts_kg, y = lactate, color = "Individuo 5")) +
  geom_point(data = datos_individuo_6, aes(x = watts_kg, y = lactate, color = "Individuo 6")) +
  geom_line(data = datos_individuo_6, aes(x = watts_kg, y = lactate, color = "Individuo 6")) +
  labs(x = "Watts/kg", y = "Lactato", color = "Individuo") + # Etiquetas de los ejes y la leyenda
  scale_color_manual(values = c("Individuo 1" = "blue", "Individuo 2" = "red", 
                                "Individuo 3" = "green", "Individuo 4" = "orange", 
                                "Individuo 5" = "purple", "Individuo 6" = "skyblue")) +
  scale_x_continuous(breaks = seq(0, 5.5, by = 0.5))  

LT12_id1 <- segmented(lm(lactate ~ watts_kg, data=datos_individuo_1), npsi = 2)
LT12_id1
LT1_id1 <- 2.463
LT2_id1 <- 4.344

LT12_id2 <- segmented(lm(lactate ~ watts_kg, data=datos_individuo_2), npsi = 2)
LT12_id2
LT1_id2 <- 1.916
LT2_id2 <- 3.156

LT12_id3 <- segmented(lm(lactate ~ watts_kg, data=datos_individuo_3), npsi = 2)
LT12_id3
LT1_id3 <- 3.299
LT2_id3 <- 4.424

LT12_id4 <- segmented(lm(lactate ~ watts_kg, data=datos_individuo_4), npsi = 2)
LT12_id4
LT1_id4 <- 2.085
LT2_id4 <- 3.250

LT12_id5 <- segmented(lm(lactate ~ watts_kg, data=datos_individuo_5), npsi = 2)
LT12_id5
LT1_id5 <- 2.229
LT2_id5 <- 3.949

LT12_id6 <- segmented(lm(lactate ~ watts_kg, data=datos_individuo_6), npsi = 2)
LT12_id6
LT1_id6 <- 2.483
LT2_id6 <- 3.462

tabla_lt <- data.frame(
  Individuo = c(1, 2, 3, 4, 5, 6),
  LT1 = c(LT1_id1, LT1_id2, LT1_id3, LT1_id4, LT1_id5, LT1_id6),
  LT2 = c(LT2_id1, LT2_id2, LT2_id3, LT2_id4, LT2_id5, LT2_id6)
)

print(tabla_lt)


ggplot(datos_individuo_1, aes(x = watts_kg, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = LT1_id1 , color = "red") + 
  geom_vline(xintercept = LT2_id1 , color = "green") + 
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT1 y el LT2")

ggplot(datos_individuo_2, aes(x = watts_kg, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = LT1_id2 , color = "red") + 
  geom_vline(xintercept = LT2_id2 , color = "green") + 
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT1 y el LT2")


ggplot(datos_individuo_3, aes(x = watts_kg, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = LT1_id3 , color = "red") + 
  geom_vline(xintercept = LT2_id3 , color = "green") + 
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT1 y el LT2")

ggplot(datos_individuo_4, aes(x = watts_kg, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = LT1_id4 , color = "red") + 
  geom_vline(xintercept = LT2_id4 , color = "green") + 
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT1 y el LT2")

ggplot(datos_individuo_5, aes(x = watts_kg, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = LT1_id5 , color = "red") + 
  geom_vline(xintercept = LT2_id5 , color = "green") + 
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT1 y el LT2")

ggplot(datos_individuo_6, aes(x = watts_kg, y = lactate)) +
  geom_point() +
  geom_line(col = "blue") +
  geom_vline(xintercept = LT1_id6 , color = "red") + 
  geom_vline(xintercept = LT2_id6 , color = "green") + 
  labs(x = "Potencia", y = "Lactacto", title = "Gráfico de Dispersión Potencia vs. Lactacto con el LT1 y el LT2")

modelo_mm_id1 <- nls(watts_kg ~ SSmicmen(lactate, a, b), 
                     start = list(a = max(datos_individuo_1_sin_reposo$watts_kg), b = 2),
                     data = datos_individuo_1_sin_reposo)
modelo_mm_id1
power_id1 <- function(lactate) {
  return((coef(modelo_mm_id1)[1] * lactate) / (coef(modelo_mm_id1)[2] + lactate))
}


modelo_mm_id2 <- nls(watts_kg ~ SSmicmen(lactate, a, b), 
                     start = list(a = max(datos_individuo_2_sin_reposo$watts_kg), b = 2),
                     data = datos_individuo_2_sin_reposo)
modelo_mm_id2
coef(modelo_mm_id2)
power_id2 <- function(lactate) {
  return((coef(modelo_mm_id2)[1] * lactate) / (coef(modelo_mm_id2)[2] + lactate))
}

modelo_mm_id3 <- nls(watts_kg ~ SSmicmen(lactate, a, b), 
                     start = list(a = max(datos_individuo_3_sin_reposo$watts_kg), b = 2),
                     data = datos_individuo_3_sin_reposo)
modelo_mm_id3
coef(modelo_mm_id3)
power_id3 <- function(lactate) {
  return((coef(modelo_mm_id3)[1] * lactate) / (coef(modelo_mm_id3)[2] + lactate))
}

modelo_mm_id4 <- nls(watts_kg ~ SSmicmen(lactate, a, b), 
                     start = list(a = max(datos_individuo_4_sin_reposo$watts_kg), b = 2),
                     data = datos_individuo_4_sin_reposo)
modelo_mm_id4
coef(modelo_mm_id4)
power_id4 <- function(lactate) {
  return((coef(modelo_mm_id4)[1] * lactate) / (coef(modelo_mm_id4)[2] + lactate))
}

modelo_mm_id5 <- nls(watts_kg ~ SSmicmen(lactate, a, b), 
                     start = list(a = max(datos_individuo_5_sin_reposo$watts_kg), b = 2),
                     data = datos_individuo_5_sin_reposo)
modelo_mm_id5
coef(modelo_mm_id5)
power_id5 <- function(lactate) {
  return((coef(modelo_mm_id5)[1] * lactate) / (coef(modelo_mm_id5)[2] + lactate))
}


modelo_mm_id6 <- nls(watts_kg ~ SSmicmen(lactate, a, b), 
                     start = list(a = max(datos_individuo_6_sin_reposo$watts_kg), b = 2),
                     data = datos_individuo_6_sin_reposo)
modelo_mm_id6
coef(modelo_mm_id6)
power_id6 <- function(lactate) {
  return((coef(modelo_mm_id6)[1]* lactate) / (coef(modelo_mm_id6)[2] + lactate))
}


tabla_parametros_mm <- data.frame(
  Individuo = c(1, 2, 3, 4, 5, 6),
  V_max = c(coef(modelo_mm_id1)[1],coef(modelo_mm_id2)[1],coef(modelo_mm_id3)[1],coef(modelo_mm_id4)[1],coef(modelo_mm_id5)[1],coef(modelo_mm_id6)[1]),
  K_M = c(coef(modelo_mm_id1)[2],coef(modelo_mm_id2)[2],coef(modelo_mm_id3)[2],coef(modelo_mm_id4)[2],coef(modelo_mm_id5)[2],coef(modelo_mm_id6)[2])
)

print(tabla_parametros_mm )



predicciones_mm_id1 <- predict(modelo_mm_id1)

error_cuadratico_mm_id1 <- mean((datos_individuo_1_sin_reposo$watts_kg - predicciones_mm_id1)^2)
error_cuadratico_mm_id1


predicciones_mm_id2 <- predict(modelo_mm_id2)

error_cuadratico_mm_id2 <- mean((datos_individuo_2_sin_reposo$watts_kg - predicciones_mm_id2)^2)
error_cuadratico_mm_id2


predicciones_mm_id3 <- predict(modelo_mm_id3)

error_cuadratico_mm_id3 <- mean((datos_individuo_3_sin_reposo$watts_kg - predicciones_mm_id3)^2)
error_cuadratico_mm_id3


predicciones_mm_id4 <- predict(modelo_mm_id4)

error_cuadratico_mm_id4 <- mean((datos_individuo_4_sin_reposo$watts_kg - predicciones_mm_id4)^2)
error_cuadratico_mm_id4


predicciones_mm_id5 <- predict(modelo_mm_id5)

error_cuadratico_mm_id5 <- mean((datos_individuo_5_sin_reposo$watts_kg - predicciones_mm_id5)^2)
error_cuadratico_mm_id5



predicciones_mm_id6 <- predict(modelo_mm_id6)

error_cuadratico_mm_id6 <- mean((datos_individuo_6_sin_reposo$watts_kg - predicciones_mm_id6)^2)
error_cuadratico_mm_id6


AIC_mm_id1 <- AIC(modelo_mm_id1)
AIC_mm_id1
BIC_mm_id1 <- BIC(modelo_mm_id1)
BIC_mm_id1

AIC_mm_id2 <- AIC(modelo_mm_id2)
AIC_mm_id2
BIC_mm_id2 <- BIC(modelo_mm_id2)
BIC_mm_id2

AIC_mm_id3 <- AIC(modelo_mm_id3)
AIC_mm_id3
BIC_mm_id3 <- BIC(modelo_mm_id3)
BIC_mm_id3

AIC_mm_id4 <- AIC(modelo_mm_id4)
AIC_mm_id4
BIC_mm_id4 <- BIC(modelo_mm_id4)
BIC_mm_id4

AIC_mm_id5 <- AIC(modelo_mm_id5)
AIC_mm_id5
BIC_mm_id5 <- BIC(modelo_mm_id5)
BIC_mm_id5

AIC_mm_id6 <- AIC(modelo_mm_id6)
AIC_mm_id6
BIC_mm_id6 <- BIC(modelo_mm_id6)
BIC_mm_id6

tabla_metricas_mm <- data.frame(
  Individuo = c(1, 2, 3, 4, 5, 6),
  EMC = c(error_cuadratico_mm_id1,error_cuadratico_mm_id2,error_cuadratico_mm_id3,error_cuadratico_mm_id4,error_cuadratico_mm_id5,error_cuadratico_mm_id6),
  AIC = c(AIC_mm_id1,AIC_mm_id2,AIC_mm_id3,AIC_mm_id4,AIC_mm_id5,AIC_mm_id6),
  BIC = c(BIC_mm_id1,BIC_mm_id2,BIC_mm_id3,BIC_mm_id4,BIC_mm_id5,BIC_mm_id6)
)

print(tabla_metricas_mm)


P_max_entre_masa_id1 <- coef(modelo_mm_id1)[1]
P_max_entre_masa_id2 <- coef(modelo_mm_id2)[1]
P_max_entre_masa_id3 <- coef(modelo_mm_id3)[1]
P_max_entre_masa_id4 <- coef(modelo_mm_id4)[1]
P_max_entre_masa_id5 <- coef(modelo_mm_id5)[1]
P_max_entre_masa_id6 <- coef(modelo_mm_id6)[1]

VO2_max_id1 <- calcular_VO2_max_2(P_max_entre_masa_id1)
VO2_max_id2 <- calcular_VO2_max_2(P_max_entre_masa_id2)
VO2_max_id3 <- calcular_VO2_max_2(P_max_entre_masa_id3)
VO2_max_id4 <- calcular_VO2_max_2(P_max_entre_masa_id4)
VO2_max_id5 <- calcular_VO2_max_2(P_max_entre_masa_id5)
VO2_max_id6 <- calcular_VO2_max_2(P_max_entre_masa_id6)

tabla_VO2_max <- data.frame(
  Individuo = c(1, 2, 3, 4, 5, 6),
  VO_2_max = c(VO2_max_id1,VO2_max_id2,VO2_max_id3,VO2_max_id4,VO2_max_id5,VO2_max_id6)
)

print(tabla_VO2_max)


lactate_id1 <- function(Pm) {
  return((3.627108 * Pm) / (7.257428 - Pm))
}

lactate_id2 <- function(Pm) {
  return((2.314049 * Pm) / (4.878828 - Pm))
}

lactate_id3 <- function(Pm) {
  return((3.13428 * Pm) / (6.582378 - Pm))
}

lactate_id4 <- function(Pm) {
  return((3.405558 * Pm) / (5.547301 - Pm))
}

lactate_id5 <- function(Pm) {
  return((2.567460 * Pm) / (5.289942 - Pm))
}


lactate_id6 <- function(Pm) {
  return((3.23895 * Pm) / (6.112095 - Pm))
}

  

P_m_ranges <- list(
  id1 = seq(1, 7.25, by = 0.02),
  id2 = seq(1, 4.87, by = 0.02),
  id3 = seq(1, 6.58, by = 0.02),
  id4 = seq(1, 5.54, by = 0.02),
  id5 = seq(1, 5.28, by = 0.02),
  id6 = seq(1, 6.11, by = 0.02)
)


# Generación de predicciones de lactato para cada individuo
predictions <- list(
  id1 = lactate_id1(P_m_ranges$id1),
  id2 = lactate_id2(P_m_ranges$id2),
  id3 = lactate_id3(P_m_ranges$id3),
  id4 = lactate_id4(P_m_ranges$id4),
  id5 = lactate_id5(P_m_ranges$id5),
  id6 = lactate_id6(P_m_ranges$id6)
)


results <- data.frame(
  Pm = c(P_m_ranges$id1, P_m_ranges$id2, P_m_ranges$id3, P_m_ranges$id4, P_m_ranges$id5, P_m_ranges$id6),
  Lactate = c(predictions$id1, predictions$id2, predictions$id3, predictions$id4, predictions$id5, predictions$id6),
  Individuo = rep(c("Id 1", "Id 2", "Id 3", "Id 4", "Id 5", "Id 6"),
                  times = sapply(P_m_ranges, length))
)

ggplot(results, aes(x = Pm, y = Lactate, color = Individuo, group = Individuo)) +
  geom_line(size = 1) +
  labs(title = "Predicción del nivel de lactato según la potencia entre masa",
       x = "Potencia entre masa (P/m)",
       y = "Nivel de lactato (mmol/l)",
       color = "Individuo") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(0, 30)+xlim(1,7)+scale_x_continuous(breaks = seq(1, 6.6 , by = 1))

calcula_diferencia_LT2 <- function(datos_individuo_1, datos_individuo_2) {
  # Ajustar modelos segmentados para los dos individuos
  LT12_id1 <- segmented(lm(lactate ~ watts_kg, data = datos_individuo_1), npsi = 2)
  LT12_id2 <- segmented(lm(lactate ~ watts_kg, data = datos_individuo_2), npsi = 2)
  
  # Obtener los umbrales LT2
  LT2_id1 <- LT12_id1$psi[2, "Est."]
  LT2_id2 <- LT12_id2$psi[2, "Est."]
  
  # Obtener las matrices de varianza-covarianza de los modelos segmentados
  vcov_id1 <- vcov(LT12_id1)
  vcov_id2 <- vcov(LT12_id2)
  
  # Calcular la diferencia de LT2 entre los dos individuos
  diferencia_LT2 <- LT2_id1 - LT2_id2
  
  # Crear la matriz de varianza-covarianza combinada
  vcov_combined <- bdiag(vcov_id1, vcov_id2)
  
  # Extraer las varianzas relevantes
  var_LT2_id1 <- vcov_id1[2, 2]
  var_LT2_id2 <- vcov_id2[2, 2]
  
  # Calcular el error estándar de la diferencia
  se_diferencia <- sqrt(var_LT2_id1 + var_LT2_id2)
  
  # Nivel de confianza
  z <- qnorm(0.975)
  
  # Calcular el intervalo de confianza
  ci_lower <- diferencia_LT2 - z * se_diferencia
  ci_upper <- diferencia_LT2 + z * se_diferencia
  
  # Calcular el p valor
  z <- diferencia_LT2 / se_diferencia
  p_value <- 2 * (1 - pnorm(abs(z)))
  
  # Crear una lista con los resultados
  resultado <- list(
    Estimate = diferencia_LT2,
    SE = se_diferencia,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    p_valor=p_value
  )
  
  return(resultado)
}


resultado1_4 <- calcula_diferencia_LT2(datos_individuo_1, datos_individuo_4)
print(resultado1_4)


funcion_comparacion_potencias <- function(funcion_power_id1, funcion_power_id2,lactate){
  return(funcion_power_id1(lactate)-funcion_power_id2(lactate))
}
power_id1(4)
power_id4(4)
funcion_comparacion_potencias(power_id1,power_id4, 4)

