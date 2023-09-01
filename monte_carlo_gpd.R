
# FUNÇÕES -----------------------------------------------------------------

# # Gerador de números aleatórios da GPD, pelo método da transformada inversa
rgpd1 <- function(n, sigma, ksi) {
  
  x <- c()
  
  for(i in 1:n) {
    u <- runif(1, min = 0.001, max = 1)
    if(ksi != 0) {
      x[i] <-  (sigma * ( (u)^(-ksi) - 1 ) ) / (ksi)
    } else {
      x[i] <- -sigma * log(u)
    }
  }
  return(x)
}


# Negativo do log da função de máxima verossimilhança
gpdloglik <- function(theta, x){
  
  sigma <- theta[1]
  ksi <- theta[2]
  n <- length(x)
  
  if(ksi != 0) {
    loglik <- (-n)*log(sigma) - ( (1/ksi + 1)*(sum( log( 1 + (ksi*x/sigma) ) ) ) )
  } else {
    loglik <- (-n)*log(sigma) - (1/sigma)*sum(x)
  }
  
  
  
  return(-loglik)
  
}

bias <- function(par, est) mean(est) - par 

eqm <- function(par, est) var(est) - bias(par, est)^2

# TESTAR FUNÇÕES ----------------------------------------------------------

n <- 1000
sigma <- 1
ksi <- -0.1

# Amostra de GPD 
Z <- rgpd1(n, sigma, ksi)

gpd_par <- c(sigma, ksi)
optim(par = gpd_par, fn = gpdloglik, x = Z, method = "BFGS", hessian = F)


# SIMULAR -----------------------------------------------------------------


faz_MC_gpd <- function(n, R, B, sigma, ksi) {
  
  # Criar vetores para receber os resultados das estimativas de MV
  sigma_hat_MV <- numeric(R)
  ksi_hat_MV <- numeric(R)
  
  # Criar vetores para receber os resultados dos parâmetros corrigidos por bootstrap
  sigma_hat_corrigido <- numeric(R)
  ksi_hat_corrigido <- numeric(R)
  
  # Criar vetores para receber os valores do bootstrap
  boot1 <- numeric(B)
  boot2 <- numeric(B)
  
  tentativas_MV <- 0
  prop_BS_missing <- 0
  
  message(paste0("Simulando valores da GPD, com sigma = ", sigma, ", ksi = ", ksi, " e amostra de tamanho = ", n))
  
  # MONTE CARLO (com R repetições)
  for(k in 1:R) {
    
    # Gerar amostra da dsitribuição GPD
    # Z <- evmix::rgpd(n, sigmau = sigma, xi = ksi)
    Z <- rgpd1(n, sigma = sigma, ksi = ksi)
    
    # Obter estimativa de MV por otimização
    esti_MV <- tryCatch(optim(par = c(sigma, ksi), fn = gpdloglik, x = Z, method = "BFGS"), error = function(e) NA)
    
    while (any(is.na(esti_MV))) {
      
      tentativas_MV <- tentativas_MV + 1
      # Z <- evmix::rgpd(n, sigmau = sigma, xi = ksi)
      Z <- rgpd1(n, sigma = sigma, ksi = ksi)
      esti_MV <- tryCatch(optim(par = c(sigma, ksi), fn = gpdloglik, x = Z, method = "BFGS"), error = function(e) NA)
    
    }
    
    sigma_hat_MV[k] <- esti_MV$par[1]
    ksi_hat_MV[k] <- esti_MV$par[2]
    
    # Fazer bootstrap para cada amostra da MC
    for(i in 1:B) {
      
      # Obter amostra usando as estimativas de MV
      # Z_boot <- evmix::rgpd(n, sigmau = sigma_hat_MV[k], xi = ksi_hat_MV[k])
      Z_boot <- rgpd1(n, sigma = sigma_hat_MV[k], ksi = ksi_hat_MV[k])
      est_boot <- tryCatch(optim(par = c(sigma_hat_MV[k], ksi_hat_MV[k]), fn = gpdloglik, x = Z_boot, method = "BFGS"), error = function(e) NA)
      
      if (!is.na(est_boot[1])) {
        boot1[i] <- est_boot$par[1]
        boot2[i] <- est_boot$par[2]
      } else {
        boot1[i] <- NA
        boot2[i] <- NA
      }

    }
    
    # Obter estimativas dos parâmetros corrigidos por bootstrap
    sigma_hat_corrigido[k] <- 2*sigma_hat_MV[k] - mean(boot1, na.rm = T)
    ksi_hat_corrigido[k] <- 2*ksi_hat_MV[k] - mean(boot2, na.rm = T)
    
    # Proporção de amostras BS que falharam na estimativa, para a k-ésima MC
    prop_BS_missing[k] <- sum(is.na(boot1))/B
    
  } 
  
  # Juntar os resultados em vetores e combiná-los em um data.frame
  media <- c(mean(sigma_hat_MV), mean(sigma_hat_corrigido), mean(ksi_hat_MV), mean(ksi_hat_corrigido))
  vies <- c(bias(sigma, sigma_hat_MV), bias(sigma, sigma_hat_corrigido), bias(ksi, ksi_hat_MV), bias(ksi, ksi_hat_corrigido))
  variancia <- c(var(sigma_hat_MV), var(sigma_hat_corrigido), var(ksi_hat_MV), var(ksi_hat_corrigido))
  EQM <- c(eqm(sigma, sigma_hat_MV), eqm(sigma, sigma_hat_corrigido), eqm(ksi, ksi_hat_MV), eqm(ksi, ksi_hat_corrigido))
  
  # Nomear as linhas do data.frame
  rn <- c(
    paste("sigma_hat_MV", ksi, n, sep = "_"), 
    paste("sigma_hat_corrigido", ksi, n, sep = "_"), 
    paste("ksi_hat_MV", ksi, n, sep = "_"),
    paste("ksi_hat_corrigido", ksi, n, sep = "_")
  )
  
  res <- data.frame(media, vies, variancia, EQM, row.names = rn)
  
  # Tentativas de BS perdidas (proporção média em relação ao total de BS)
  mean_prop_BS_missing <- round(mean(prop_BS_missing), 4)

  
  # Retornar dataframe final
  return(list(est = res, tentativas_MV = tentativas_MV, BS_missing = mean_prop_BS_missing))  
  
}

# OBTER RESULTADOS --------------------------------------------------------


# Testes
n <- 100
sigma <- 1
ksi <- -0.1
R <- 100
B <- 100

faz_MC_gpd(n, R, B, sigma, ksi)

n <- 100
sigma <- 1
ksi <- -0.1
R <- 1000
B <- 1000

faz_MC_gpd(n, R, B, sigma, ksi)

####

n_values <- c(50, 100, 200, 400)
ksi_values <- c(-0.2,-0.1, 0, 0.1, 0.2)
R <- 5000
B <- 500


##  Rodar simulação em paralelo -----
library(doParallel)

# Setar número de CPUs a serem utilizados
ncpu <- detectCores() - 1
cluster <- parallel::makeCluster(ncpu, type = "PSOCK")
registerDoParallel(cl = cluster)

# Checar se os workers foram devidamente registrados
foreach::getDoParRegistered() 
foreach::getDoParWorkers()

# Todas as combinações possíveis entre tamanho das amostras e valores de ksi
cb <- expand.grid(n_values, ksi_values)

# Obter os resultados
system.time(
  
  res2 <- foreach(x = 1:nrow(cb), .combine = "rbind") %dopar% {
    
    faz_MC_gpd(n = cb[x, 1], R = R, B = B, sigma = 1, ksi = cb[x, 2])
    
  }
  
)

stopCluster(cluster)

# Tempo decorrido (aproximadamente 100 min)

# usuário   sistema decorrido 
# 18.282    21.533  5990.957 

save(res2, file = "/media/iaradsouza/DATA1/Estatística/8_periodo/Métodos Computacionais/Projeto/Ex_RB_1000_02.rda")
