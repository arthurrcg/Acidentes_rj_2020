# Acidentes_rj_2020
# Análise de Correspondência Simples de variáveis de acidentes de trânsito na cidade do Rio de Janeiro em 2020 em R


# Pacotes a serem instalados e carregados ---------------------------------

library(plotly)
library(tidyverse)
library(ggrepel)
library(sjPlot)
library(FactoMineR)
library(amap)
library(ade4)
library(fastDummies)
library(knitr)
library(kableExtra)

# Análise de Correspondência Simples (ANACOR) - Exemplo 1 -----------------

#Carregando a base de dados datatran2020


#Observado os dados carregados
datatran2020 %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                font_size = 12)

#Tabelas de frequências
summary(datatran2020)

#Criando uma tabela de contingências
tab <- table(datatran2020$classificacao_acidente, 
             datatran2020$tipo_acidente)

tab

#Exemplo de uma tabela de contingências mais elegante
sjt.xtab(var.row = datatran2020$classificacao_acidente,
         var.col = datatran2020$tipo_acidente)

#Exemplo de uma tabela de contingências mais elegante
sjt.xtab(var.row = datatran2020$classificacao_acidente,
         var.col = datatran2020$tipo_acidente,
         show.exp = TRUE)

#Teste Qui-Quadrado
chi2 <- chisq.test(tab)
chi2

chi2$statistic
chi2$parameter
chi2$p.value
chi2$method
chi2$data.name
chi2$observed
chi2$expected
chi2$residuals #Resíduos PADRONIZADOS
chi2$stdres #Resíduos PADRONIZADOS AJUSTADOS

#Mapa de calor dos resíduos padronizados ajustados
data.frame(chi2$stdres) %>%
  rename(perfil = 1,
         aplicacao = 2) %>% 
  ggplot(aes(x = fct_rev(perfil), y = aplicacao, fill = Freq, label = round(Freq,3))) +
  geom_tile() +
  geom_text(size = 3) +
  scale_fill_gradient2(low = "darkorchid", 
                       mid = "white", 
                       high = "orange",
                       midpoint = 0) +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.text.x = element_text())

#Decomposição da inércia principal total /// IT = inercia total
It <- chi2$statistic/nrow(datatran2020)
It

#Construindo a matriz P
P <- 1/nrow(datatran2020) * tab
P

#Column profile
column_profile <- apply(tab, MARGIN = 1, FUN = sum) / nrow(datatran2020)
column_profile

#Row profile
row_profile <- apply(tab, MARGIN = 2, FUN = sum) / nrow(datatran2020)
row_profile

#Matriz Dl
Dl <- diag(column_profile)
Dl

#Matriz Dc
Dc <- diag(row_profile)
Dc

#Matriz lc'
lc <- column_profile %o% row_profile
lc

#Matriz A
A <- diag(diag(Dl) ^ (-1/2)) %*% (P - lc) %*% diag(diag(Dc) ^ (-1/2))
A


#Curiosidade:
A_matriz <- chi2$residuals / sqrt(nrow(datatran2020))
A_matriz

#Matriz W
W_matriz <- t(A_matriz) %*% A_matriz
W_matriz

#Extraindo os eigenvalues da matriz W
eigenvalues <- eigen(W_matriz)
eigenvalues

sum(eigenvalues$values) #It
It

#Dimensionalidade dos dados
dimensoes <- min(nrow(A_matriz) - 1, ncol(A_matriz) - 1)
dimensoes

#Variabilidade dos dados explicada
var_explicada <- eigenvalues$values[1:2] / It
var_explicada

#Cálculo das coordenadas do mapa perceptual

#Decomposição do valor singular da matriz A
decomp <- svd(x = A_matriz,
              nu = dimensoes,
              nv = dimensoes)

decomp

#Variável em linha - coordenada no eixo das abcissas
Xl_perfil <- diag((decomp$d[1]) * diag(diag(Dl)^(-1/2)) * decomp$u[,1])
Xl_perfil

#Variável em linha - coordenada no eixo das ordenadas
Yl_perfil <- diag((decomp$d[2]) * diag(diag(Dl)^(-1/2)) * decomp$u[,2])
Yl_perfil

#Variável em coluna - coordenada no eixo das abcissas
Xc_aplicacao <- diag((decomp$d[1]) * diag(diag(Dc)^(-1/2)) * decomp$v[,1])
Xc_aplicacao

#Variável em coluna - coordenada no eixo das ordenadas
Yc_aplicacao <- diag((decomp$d[2]) * diag(diag(Dc)^(-1/2)) * decomp$v[,2])
Yc_aplicacao

#Elaborando a ANACOR diretamente
anacor <- CA(tab)
anacor$row$coord
anacor$col$coord

#Plotando o mapa perceptual de maneira mais elegante:

#Capturando todas as coordenadas num só objeto
ca_coordenadas <- rbind(anacor$row$coord, anacor$col$coord)
ca_coordenadas

#Capturando a quantidade de categorias por variável
id_var <- apply(datatran2020[,4:3],
                MARGIN=2,
                FUN = function(x)nlevels(as.factor(x)))

id_var

#Juntando as coordenadas e as categorias capturadas anteriormente
ca_coordenadas_final <- data.frame(ca_coordenadas, 
                                   Variable = rep(names(id_var), id_var))

ca_coordenadas_final

#Mapa perceptual bidimensional
ca_coordenadas_final %>% 
  rownames_to_column() %>% 
  rename(Category = 1) %>% 
  ggplot(aes(x = Dim.1, y = Dim.2, label = Category, color = Variable)) +
  geom_point() +
  geom_label_repel() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = paste("Dimension 1:", paste0(round(anacor$eig[1,2], digits = 2), "%")),
       y = paste("Dimension 2:", paste0(round(anacor$eig[2,2], digits = 2), "%"))) +
  scale_color_manual("Variable:",
                     values = c("darkorchid", "orange")) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect("NA"),
        panel.grid = element_line("gray95"),
        legend.position = "none")

#Repetindo o mapa de calor dos resíduos padronizados ajustados
data.frame(chi2$stdres) %>%
  rename(perfil = 1,
         aplicacao = 2) %>% 
  ggplot(aes(x = fct_rev(perfil), y = aplicacao, fill = Freq, label = round(Freq,3))) +
  geom_tile() +
  geom_text(size = 3) +
  scale_fill_gradient2(low = "darkorchid", 
                       mid = "white", 
                       high = "orange",
                       midpoint = 0) +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.text.x = element_text())

