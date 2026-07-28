################################################################################
# ESTUDO DE SIMULACAO ROBUSTO PARA AVALIAR MIBLUP VS GBLUP
# Autor: Ricardo Antonio Ruiz Cardozo
# Objetivo:
#   - Comparar GBLUP, MIBLUP, RF e MLP em varios cenarios simulados
#   - Quantificar desempenho medio, variabilidade e IC95%
#   - Verificar se SNPs com maior MI estao enriquecidos em QTLs causais
################################################################################

rm(list = ls())

################################################################################
# 0) PACOTES
################################################################################
required_pkgs <- c("rrBLUP", "randomForest", "dplyr", "purrr", "tibble", "tidyr")

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(rrBLUP)
library(randomForest)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

# Deep Learning opcional
use_mlp <- FALSE
if (use_mlp) {
  if (!requireNamespace("keras", quietly = TRUE)) install.packages("keras")
  library(keras)
}

################################################################################
# 1) FUNCOES AUXILIARES
################################################################################

#------------------------------------------------------------------------------
# 1.1 Simular frequencias alélicas realistas
#------------------------------------------------------------------------------
simulate_maf <- function(m, min_maf = 0.05, max_maf = 0.50) {
  runif(m, min = min_maf, max = max_maf)
}

#------------------------------------------------------------------------------
# 1.2 Simular matriz de SNPs sob HWE
#     Genotipos: 0, 1, 2
#------------------------------------------------------------------------------
simulate_genotypes <- function(n, m, maf = NULL) {
  if (is.null(maf)) maf <- simulate_maf(m)
  
  Z <- sapply(maf, function(p) {
    rbinom(n = n, size = 2, prob = p)
  })
  
  Z <- matrix(Z, nrow = n, ncol = m)
  storage.mode(Z) <- "numeric"
  Z
}

#------------------------------------------------------------------------------
# 1.3 Padronizacao simples
#------------------------------------------------------------------------------
scale_vec <- function(x) {
  s <- sd(x)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x)) / s)
}

#------------------------------------------------------------------------------
# 1.4 Cálculo da Informação Mútua entre SNP discreto e fenotipo discretizado
#------------------------------------------------------------------------------
calc_mi_discrete <- function(x, y_disc) {
  px  <- prop.table(table(x))
  py  <- prop.table(table(y_disc))
  pxy <- prop.table(table(x, y_disc))
  
  px  <- px[px > 0]
  py  <- py[py > 0]
  pxy <- pxy[pxy > 0]
  
  Hx  <- -sum(px * log2(px))
  Hy  <- -sum(py * log2(py))
  Hxy <- -sum(pxy * log2(pxy))
  
  mi <- Hx + Hy - Hxy
  max(0, mi)
}

#------------------------------------------------------------------------------
# 1.5 Calcular vetor de MI para todos os SNPs
#------------------------------------------------------------------------------
compute_mi_weights <- function(Z_train, y_train, n_bins = 10, eps = 1e-8) {
  y_disc <- cut(y_train, breaks = n_bins, labels = FALSE, include.lowest = TRUE)
  
  mi <- apply(Z_train, 2, function(snp) calc_mi_discrete(snp, y_disc))
  
  # Evitar divisao por zero em situacoes degeneradas
  if (mean(mi) <= eps) {
    w <- rep(1, length(mi))
  } else {
    w <- mi / mean(mi)
  }
  
  # Evita pesos exatamente zero
  w <- pmax(w, eps)
  
  list(mi = mi, weights = w)
}

#------------------------------------------------------------------------------
# 1.6 Centralizacao estilo VanRaden
#------------------------------------------------------------------------------
center_genotypes <- function(Z_train, Z_test = NULL) {
  p_freq <- colMeans(Z_train) / 2
  P_vec  <- 2 * p_freq
  
  Z_train_c <- sweep(Z_train, 2, P_vec, FUN = "-")
  
  if (!is.null(Z_test)) {
    Z_test_c <- sweep(Z_test, 2, P_vec, FUN = "-")
    return(list(Z_train = Z_train_c, Z_test = Z_test_c, p_freq = p_freq))
  }
  
  list(Z_train = Z_train_c, p_freq = p_freq)
}

#------------------------------------------------------------------------------
# 1.7 Gerar fenotipo sob diferentes arquiteturas genéticas
#------------------------------------------------------------------------------
simulate_phenotype <- function(Z,
                               architecture = c("additive", "add_epi", "sparse"),
                               n_qtl = 20,
                               epi_strength = 0,
                               h2 = 0.50,
                               effect_sd = 1.0,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  architecture <- match.arg(architecture)
  n <- nrow(Z)
  m <- ncol(Z)
  
  qtl_idx <- sort(sample(seq_len(m), size = n_qtl, replace = FALSE))
  
  # efeitos aditivos
  beta <- rep(0, m)
  
  if (architecture == "sparse") {
    # poucos efeitos fortes
    beta_qtl <- rnorm(n_qtl, mean = 0, sd = effect_sd * 2.5)
  } else {
    beta_qtl <- rnorm(n_qtl, mean = 0, sd = effect_sd)
  }
  
  beta[qtl_idx] <- beta_qtl
  g_add <- as.numeric(Z[, qtl_idx, drop = FALSE] %*% beta_qtl)
  
  # componente epistatico
  g_epi <- rep(0, n)
  epi_pairs <- NULL
  
  if (architecture == "add_epi" && epi_strength > 0 && n_qtl >= 4) {
    # cria pares epistaticos entre QTLs
    qtl_perm <- sample(qtl_idx)
    pair_mat <- matrix(qtl_perm[1:(2 * floor(n_qtl / 2))], ncol = 2, byrow = TRUE)
    
    epi_pairs <- pair_mat
    epi_terms <- apply(pair_mat, 1, function(idx) Z[, idx[1]] * Z[, idx[2]])
    epi_terms <- as.matrix(epi_terms)
    
    gamma <- rnorm(n = ncol(epi_terms), mean = 0, sd = epi_strength)
    g_epi <- as.numeric(epi_terms %*% gamma)
  }
  
  g_true_raw <- g_add + g_epi
  vg <- var(g_true_raw)
  
  if (vg <= 0) stop("Variância genética não positiva.")
  
  # Ajustar erro para atingir herdabilidade aproximada
  ve <- vg * (1 - h2) / h2
  e  <- rnorm(n, mean = 0, sd = sqrt(ve))
  
  y <- g_true_raw + e
  
  list(
    y = y,
    g_true = g_true_raw,
    qtl_idx = qtl_idx,
    beta = beta,
    epi_pairs = epi_pairs
  )
}

#------------------------------------------------------------------------------
# 1.8 Métricas
#------------------------------------------------------------------------------
safe_cor <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y)
}

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

#------------------------------------------------------------------------------
# 1.9 Intervalo de confianca normal aproximado
#------------------------------------------------------------------------------
mean_ci95 <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(tibble(mean = NA_real_, sd = NA_real_, se = NA_real_,
                  lwr = NA_real_, upr = NA_real_, n = 0))
  }
  m  <- mean(x)
  s  <- sd(x)
  se <- s / sqrt(n)
  z  <- 1.96
  tibble(mean = m, sd = s, se = se, lwr = m - z * se, upr = m + z * se, n = n)
}

#------------------------------------------------------------------------------
# 1.10 Enriquecimento dos top SNPs por MI em QTLs verdadeiros
#      Mede quantos QTLs estão entre os top-k SNPs por MI
#------------------------------------------------------------------------------
qtl_enrichment <- function(mi_scores, qtl_idx, top_prop = 0.05) {
  m <- length(mi_scores)
  k <- max(1, ceiling(m * top_prop))
  top_idx <- order(mi_scores, decreasing = TRUE)[1:k]
  
  hits <- sum(top_idx %in% qtl_idx)
  expected <- k * (length(qtl_idx) / m)
  fold_enrichment <- ifelse(expected > 0, hits / expected, NA_real_)
  
  tibble(
    top_k = k,
    qtl_hits = hits,
    expected_hits = expected,
    fold_enrichment = fold_enrichment
  )
}

################################################################################
# 2) AJUSTE DOS MODELOS
################################################################################

fit_gblup <- function(Z_train, y_train, Z_test) {
  t0 <- proc.time()[3]
  fit <- rrBLUP::mixed.solve(y = y_train, Z = Z_train)
  pred <- as.numeric(fit$beta) + as.numeric(Z_test %*% fit$u)
  t1 <- proc.time()[3]
  
  list(pred = pred, time = t1 - t0)
}

fit_miblup <- function(Z_train, y_train, Z_test, n_bins = 10) {
  t0 <- proc.time()[3]
  
  mi_obj <- compute_mi_weights(Z_train, y_train, n_bins = n_bins)
  w <- mi_obj$weights
  mi_scores <- mi_obj$mi
  
  # ponderacao diagonal via sqrt(w)
  Zw_train <- sweep(Z_train, 2, sqrt(w), `*`)
  Zw_test  <- sweep(Z_test,  2, sqrt(w), `*`)
  
  fit <- rrBLUP::mixed.solve(y = y_train, Z = Zw_train)
  pred <- as.numeric(fit$beta) + as.numeric(Zw_test %*% fit$u)
  
  t1 <- proc.time()[3]
  
  list(pred = pred, time = t1 - t0, mi_scores = mi_scores, weights = w)
}

fit_rf <- function(Z_train, y_train, Z_test, ntree = 300, nodesize = 5) {
  t0 <- proc.time()[3]
  fit <- randomForest::randomForest(
    x = Z_train,
    y = y_train,
    ntree = ntree,
    nodesize = nodesize
  )
  pred <- predict(fit, Z_test)
  t1 <- proc.time()[3]
  
  list(pred = as.numeric(pred), time = t1 - t0)
}

fit_mlp <- function(Z_train, y_train, Z_test,
                    epochs = 80, batch_size = 32,
                    lr = 0.001, verbose = 0) {
  if (!use_mlp) stop("use_mlp = FALSE. Ative para usar MLP.")
  
  t0 <- proc.time()[3]
  
  keras::k_clear_session()
  
  Xtr <- scale(Z_train)
  Xte <- scale(Z_test,
               center = attr(Xtr, "scaled:center"),
               scale  = attr(Xtr, "scaled:scale"))
  
  model <- keras_model_sequential() |>
    layer_dense(units = 128, activation = "relu", input_shape = ncol(Xtr)) |>
    layer_dropout(rate = 0.30) |>
    layer_dense(units = 64, activation = "relu") |>
    layer_dropout(rate = 0.20) |>
    layer_dense(units = 1)
  
  model |>
    compile(
      optimizer = optimizer_adam(learning_rate = lr),
      loss = "mse"
    )
  
  model |>
    fit(
      x = Xtr,
      y = y_train,
      epochs = epochs,
      batch_size = batch_size,
      validation_split = 0.10,
      verbose = verbose
    )
  
  pred <- model |> predict(Xte, verbose = 0)
  pred <- as.numeric(pred)
  t1 <- proc.time()[3]
  
  list(pred = pred, time = t1 - t0)
}

################################################################################
# 3) UMA REPETICAO COMPLETA DE UM CENARIO
################################################################################

run_one_replication <- function(n = 800,
                                m = 2000,
                                architecture = "add_epi",
                                n_qtl = 20,
                                epi_strength = 1.0,
                                h2 = 0.50,
                                train_prop = 0.80,
                                seed = 1,
                                rf_ntree = 300,
                                mi_bins = 10) {
  set.seed(seed)
  
  # 1. Simular genotipos
  maf <- simulate_maf(m)
  Z <- simulate_genotypes(n = n, m = m, maf = maf)
  
  # 2. Simular fenotipo
  pheno <- simulate_phenotype(
    Z = Z,
    architecture = architecture,
    n_qtl = n_qtl,
    epi_strength = epi_strength,
    h2 = h2,
    seed = seed + 1000
  )
  
  y <- pheno$y
  g_true <- pheno$g_true
  qtl_idx <- pheno$qtl_idx
  
  # 3. Split treino/teste
  idx_train <- sample(seq_len(n), size = floor(train_prop * n), replace = FALSE)
  idx_test  <- setdiff(seq_len(n), idx_train)
  
  Z_train <- Z[idx_train, , drop = FALSE]
  Z_test  <- Z[idx_test,  , drop = FALSE]
  
  y_train <- y[idx_train]
  y_test  <- y[idx_test]
  
  g_test_true <- g_true[idx_test]
  
  # 4. Centralizacao
  cen <- center_genotypes(Z_train, Z_test)
  Z_train_c <- cen$Z_train
  Z_test_c  <- cen$Z_test
  
  # 5. GBLUP
  gblup <- fit_gblup(Z_train = Z_train_c, y_train = y_train, Z_test = Z_test_c)
  
  # 6. MIBLUP
  miblup <- fit_miblup(Z_train = Z_train_c, y_train = y_train, Z_test = Z_test_c, n_bins = mi_bins)
  
  # 7. RF
  rf <- fit_rf(Z_train = Z_train, y_train = y_train, Z_test = Z_test, ntree = rf_ntree)
  
  # 8. MLP opcional
  mlp <- NULL
  if (use_mlp) {
    mlp <- fit_mlp(Z_train = Z_train, y_train = y_train, Z_test = Z_test)
  }
  
  # 9. Enriquecimento MI
  enrich <- qtl_enrichment(miblup$mi_scores, qtl_idx = qtl_idx, top_prop = 0.05)
  
  # 10. Organizar resultados
  res <- tibble(
    model = c("GBLUP", "MIBLUP", "RF", if (use_mlp) "MLP" else NULL),
    cor_y = c(
      safe_cor(y_test, gblup$pred),
      safe_cor(y_test, miblup$pred),
      safe_cor(y_test, rf$pred),
      if (use_mlp) safe_cor(y_test, mlp$pred) else NULL
    ),
    cor_g = c(
      safe_cor(g_test_true, gblup$pred),
      safe_cor(g_test_true, miblup$pred),
      safe_cor(g_test_true, rf$pred),
      if (use_mlp) safe_cor(g_test_true, mlp$pred) else NULL
    ),
    rmse = c(
      rmse(y_test, gblup$pred),
      rmse(y_test, miblup$pred),
      rmse(y_test, rf$pred),
      if (use_mlp) rmse(y_test, mlp$pred) else NULL
    ),
    time_sec = c(
      gblup$time,
      miblup$time,
      rf$time,
      if (use_mlp) mlp$time else NULL
    )
  )
  
  list(
    results = res,
    enrichment = enrich,
    qtl_idx = qtl_idx
  )
}

################################################################################
# 4) GRID DE CENARIOS
################################################################################

# scenario_grid <- tidyr::crossing(
#   scenario_id = 1:6,
#   n = c(600, 800, 1000),
#   m = c(1500, 2000),
#   architecture = c("additive", "add_epi", "sparse"),
#   h2 = c(0.30, 0.60),
#   n_qtl = c(10, 30),
#   epi_strength = c(0.0, 1.2)
# ) |>
#   # filtrar combinacoes incoerentes
#   dplyr::filter(
#     !(architecture == "additive" & epi_strength > 0),
#     !(architecture == "sparse"   & epi_strength > 0)
#   ) |>
#   # reduzir volume para algo computacionalmente mais viavel
#   dplyr::slice_sample(n = 12)

scenario_grid <- tidyr::crossing(
  n = c(600, 800, 1000),
  m = c(1500, 2000),
  architecture = c("additive", "add_epi", "sparse"),
  h2 = c(0.30, 0.60),
  n_qtl = c(10, 30),
  epi_strength = c(0.0, 1.2)
) |>
  dplyr::filter(
    !(architecture == "additive" & epi_strength > 0),
    !(architecture == "sparse"   & epi_strength > 0)
  ) |>
  dplyr::slice_sample(n = 12) |>
  dplyr::mutate(scenario_row = dplyr::row_number())

scenario_grid
print(scenario_grid)

################################################################################
# 5) RODAR O ESTUDO
################################################################################

n_reps <- 100
base_seed <- 20260414

all_results <- list()
all_enrichment <- list()

counter <- 1

for (i in seq_len(nrow(scenario_grid))) {
  sc <- scenario_grid[i, ]
  
  cat("\n====================================================\n")
  cat("Cenário", i, "de", nrow(scenario_grid), "\n")
  print(sc)
  cat("====================================================\n")
  
  for (r in seq_len(n_reps)) {
    seed_now <- base_seed + i * 10000 + r
    
    out <- run_one_replication(
      n = sc$n,
      m = sc$m,
      architecture = sc$architecture,
      n_qtl = sc$n_qtl,
      epi_strength = sc$epi_strength,
      h2 = sc$h2,
      seed = seed_now,
      rf_ntree = 300,
      mi_bins = 10
    )
    
    tmp_res <- out$results |>
      mutate(
        scenario_row = i,
        replication = r,
        n = sc$n,
        m = sc$m,
        architecture = sc$architecture,
        h2 = sc$h2,
        n_qtl = sc$n_qtl,
        epi_strength = sc$epi_strength
      )
    
    tmp_enr <- out$enrichment |>
      mutate(
        scenario_row = i,
        replication = r,
        n = sc$n,
        m = sc$m,
        architecture = sc$architecture,
        h2 = sc$h2,
        n_qtl = sc$n_qtl,
        epi_strength = sc$epi_strength
      )
    
    all_results[[counter]] <- tmp_res
    all_enrichment[[counter]] <- tmp_enr
    counter <- counter + 1
    
    if (r %% 5 == 0) {
      cat("  repetição", r, "de", n_reps, "concluída\n")
    }
  }
}

results_df <- bind_rows(all_results)
enrichment_df <- bind_rows(all_enrichment)

################################################################################
# 6) RESUMOS POR CENARIO E MODELO
################################################################################

summary_perf <- results_df |>
  group_by(n, m, architecture, h2, n_qtl, epi_strength, model) |>
  summarise(
    cor_y_mean = mean(cor_y, na.rm = TRUE),
    cor_y_sd   = sd(cor_y, na.rm = TRUE),
    cor_y_lwr  = cor_y_mean - 1.96 * cor_y_sd / sqrt(sum(is.finite(cor_y))),
    cor_y_upr  = cor_y_mean + 1.96 * cor_y_sd / sqrt(sum(is.finite(cor_y))),
    
    cor_g_mean = mean(cor_g, na.rm = TRUE),
    cor_g_sd   = sd(cor_g, na.rm = TRUE),
    cor_g_lwr  = cor_g_mean - 1.96 * cor_g_sd / sqrt(sum(is.finite(cor_g))),
    cor_g_upr  = cor_g_mean + 1.96 * cor_g_sd / sqrt(sum(is.finite(cor_g))),
    
    rmse_mean  = mean(rmse, na.rm = TRUE),
    rmse_sd    = sd(rmse, na.rm = TRUE),
    
    time_mean  = mean(time_sec, na.rm = TRUE),
    time_sd    = sd(time_sec, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_perf)

################################################################################
# 7) GANHO RELATIVO DO MIBLUP SOBRE GBLUP
################################################################################

wide_perf <- results_df |>
  select(n, m, architecture, h2, n_qtl, epi_strength, replication, model, cor_y, cor_g, rmse, time_sec) |>
  tidyr::pivot_wider(
    names_from = model,
    values_from = c(cor_y, cor_g, rmse, time_sec)
  )

# wide_perf <- results_df |>
#   select(n, m, architecture, h2, n_qtl, epi_strength, replication, model, cor_y, cor_g, rmse, time_sec) |>
#   tidyr::pivot_wider(
#     names_from = model,
#     values_from = c(cor_y, cor_g, rmse, time_sec)
#   )

# wide_perf <- results_df |>
#   select(
#     scenario_row, n, m, architecture, h2, n_qtl, epi_strength,
#     replication, model, cor_y, cor_g, rmse, time_sec
#   ) |>
#   tidyr::pivot_wider(
#     names_from = model,
#     values_from = c(cor_y, cor_g, rmse, time_sec)
#   )

gain_summary <- wide_perf |>
  mutate(
    gain_cor_y = cor_y_MIBLUP - cor_y_GBLUP,
    gain_cor_g = cor_g_MIBLUP - cor_g_GBLUP,
    gain_rmse  = rmse_GBLUP - rmse_MIBLUP,   # positivo = MIBLUP melhor
    extra_time = time_sec_MIBLUP - time_sec_GBLUP
  ) |>
  group_by(n, m, architecture, h2, n_qtl, epi_strength) |>
  summarise(
    gain_cor_y_mean = mean(gain_cor_y, na.rm = TRUE),
    gain_cor_y_sd   = sd(gain_cor_y, na.rm = TRUE),
    gain_cor_g_mean = mean(gain_cor_g, na.rm = TRUE),
    gain_rmse_mean  = mean(gain_rmse, na.rm = TRUE),
    extra_time_mean = mean(extra_time, na.rm = TRUE),
    prop_miblup_better_y = mean(gain_cor_y > 0, na.rm = TRUE),
    prop_miblup_better_g = mean(gain_cor_g > 0, na.rm = TRUE),
    .groups = "drop"
  )

gain_summary <- wide_perf |>
  mutate(
    gain_cor_y = cor_y_MIBLUP - cor_y_GBLUP,
    gain_cor_g = cor_g_MIBLUP - cor_g_GBLUP,
    gain_rmse  = rmse_GBLUP - rmse_MIBLUP,
    extra_time = time_sec_MIBLUP - time_sec_GBLUP
  ) |>
  group_by(n, m, architecture, h2, n_qtl, epi_strength) |>
  summarise(
    gain_cor_y_mean = mean(gain_cor_y, na.rm = TRUE),
    gain_cor_y_sd   = sd(gain_cor_y, na.rm = TRUE),
    gain_cor_g_mean = mean(gain_cor_g, na.rm = TRUE),
    gain_cor_g_sd   = sd(gain_cor_g, na.rm = TRUE),
    gain_rmse_mean  = mean(gain_rmse, na.rm = TRUE),
    gain_rmse_sd    = sd(gain_rmse, na.rm = TRUE),
    extra_time_mean = mean(extra_time, na.rm = TRUE),
    extra_time_sd   = sd(extra_time, na.rm = TRUE),
    prop_miblup_better_y = mean(gain_cor_y > 0, na.rm = TRUE),
    prop_miblup_better_g = mean(gain_cor_g > 0, na.rm = TRUE),
    .groups = "drop"
  )

print(gain_summary)

################################################################################
# 8) RESUMO DO ENRIQUECIMENTO DE QTLs ENTRE OS TOP SNPs POR MI
################################################################################

summary_enrichment <- enrichment_df |>
  group_by(n, m, architecture, h2, n_qtl, epi_strength) |>
  summarise(
    qtl_hits_mean = mean(qtl_hits, na.rm = TRUE),
    fold_enrichment_mean = mean(fold_enrichment, na.rm = TRUE),
    fold_enrichment_sd   = sd(fold_enrichment, na.rm = TRUE),
    .groups = "drop"
  )


print(summary_enrichment)

################################################################################
# 9) TESTES PAREADOS SIMPLES MIBLUP VS GBLUP POR CENARIO
#    (mostrar consistencia)
################################################################################

paired_tests <- wide_perf |>
  group_by(n, m, architecture, h2, n_qtl, epi_strength) |>
  group_modify(~{
    x <- .x$cor_y_MIBLUP
    y <- .x$cor_y_GBLUP
    
    if (length(x) >= 3 && all(is.finite(x)) && all(is.finite(y))) {
      tt <- t.test(x, y, paired = TRUE)
      tibble(
        mean_diff = mean(x - y),
        p_value = tt$p.value,
        ci_lwr = tt$conf.int[1],
        ci_upr = tt$conf.int[2]
      )
    } else {
      tibble(
        mean_diff = NA_real_,
        p_value = NA_real_,
        ci_lwr = NA_real_,
        ci_upr = NA_real_
      )
    }
  }) |>
  ungroup()

print(paired_tests)

################################################################################
# 10) TABELA GERAL CONSOLIDADA
################################################################################

final_summary <- summary_perf |>
  left_join(
    gain_summary,
    by = c("n", "m", "architecture", "h2", "n_qtl", "epi_strength")
  ) |>
  left_join(
    summary_enrichment,
    by = c("n", "m", "architecture", "h2", "n_qtl", "epi_strength")
  )

print(final_summary)

################################################################################
# 11) EXPORTAR RESULTADOS
################################################################################

write.csv(results_df, "results_replications.csv", row.names = FALSE)
write.csv(summary_perf, "summary_performance_by_scenario.csv", row.names = FALSE)
write.csv(gain_summary, "gain_miblup_vs_gblup.csv", row.names = FALSE)
write.csv(summary_enrichment, "summary_qtl_enrichment.csv", row.names = FALSE)
write.csv(paired_tests, "paired_tests_miblup_vs_gblup.csv", row.names = FALSE)
write.csv(final_summary, "final_summary_all.csv", row.names = FALSE)

cat("\nArquivos exportados com sucesso.\n")


library(ggplot2)
library(dplyr)

# Carregando seus dados originais
dados <- read.csv("results_replications.csv")

# Filtrando apenas o cenário mais difícil (Epistasia com baixa herdabilidade)
dados_epi <- dados %>% 
  filter(architecture == "add_epi", h2 == 0.3)

# 1. Gráfico de Acurácia (cor_g)
g1 <- ggplot(dados_epi, aes(x = model, y = cor_g, fill = model)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = c("GBLUP" = "#E74C3C", "MIBLUP" = "#2ECC71", "RF" = "#3498DB")) +
  labs(title = "Acurácia de Predição (Valor Genético Real)",
       subtitle = "Cenário com Interação Epistática e Baixa Herdabilidade",
       x = "Modelo", y = "Acurácia (cor_g)") +
  theme(legend.position = "none", text = element_text(size = 14))

# 2. Gráfico de Tempo Computacional
g2 <- ggplot(dados_epi, aes(x = model, y = time_sec, fill = model)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = c("GBLUP" = "#E74C3C", "MIBLUP" = "#2ECC71", "RF" = "#3498DB")) +
  labs(title = "Custo Computacional",
       subtitle = "Tempo médio de execução (segundos)",
       x = "Modelo", y = "Tempo (s)") +
  theme(legend.position = "none", text = element_text(size = 14))

print(g1)
print(g2)



## Gráficos Ganhos
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# Ler os resultados agregados
gain_summary <- read_csv("gain_miblup_vs_gblup.csv")

# Criar rótulo do cenário
gain_plot <- gain_summary %>%
  mutate(
    scenario = paste0(
      "n=", n, ", m=", m,
      ", ", architecture,
      ", h²=", h2,
      ", QTL=", n_qtl,
      ", epi=", epi_strength
    ),
    # IC aproximado usando 100 repetições
    se = gain_cor_y_sd / sqrt(100),
    lwr = gain_cor_y_mean - 1.96 * se,
    upr = gain_cor_y_mean + 1.96 * se
  ) %>%
  arrange(gain_cor_y_mean) %>%
  mutate(scenario = factor(scenario, levels = scenario))

gain_plot$cenarios <- 

# Gráfico
ggplot(gain_plot, aes(x = gain_cor_y_mean, y = scenario)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), width = 0.2) +
  geom_point(size = 2.8) +
  labs(
    x = expression(Delta*"cor"[y]*" (MIBLUP - GBLUP)"),
    y = "Cenário",
    title = "Ganhos do MIBLUP sobre o GBLUP por cenário"
  ) +
  theme_minimal(base_size = 11)


gain_plot_g <- gain_summary %>%
  mutate(
    scenario = paste0(
      "n=", n, ", m=", m,
      ", ", architecture,
      ", h²=", h2,
      ", QTL=", n_qtl,
      ", epi=", epi_strength
    ),
    se = gain_cor_g_sd / sqrt(100),
    lwr = gain_cor_g_mean - 1.96 * se,
    upr = gain_cor_g_mean + 1.96 * se
  ) %>%
  arrange(gain_cor_g_mean) %>%
  mutate(scenario = factor(scenario, levels = scenario))

ggplot(gain_plot_g, aes(x = gain_cor_g_mean, y = scenario)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), width = 0.2) +
  geom_point(size = 2.8) +
  labs(
    x = expression(Delta*"cor"[g]*" (MIBLUP - GBLUP)"),
    y = "Cenário",
    title = "Ganhos do MIBLUP sobre o GBLUP para o valor genético verdadeiro"
  ) +
  theme_minimal(base_size = 11)

str(gain_plot)


## =========================================================
## Gráficos de ganhos do MIBLUP vs GBLUP (cenários C1–C12)
## =========================================================

library(readr)
library(dplyr)
library(ggplot2)

# ---------------------------------------------------------
# 1. Ler dados
# ---------------------------------------------------------
gain_summary <- read_csv("gain_miblup_vs_gblup.csv")

# Padronizar (evita erro de join)
gain_summary <- gain_summary %>%
  mutate(
    h2 = round(h2, 1),
    epi_strength = round(epi_strength, 1)
  )

# ---------------------------------------------------------
# 2. Mapa dos cenários (igual à tabela do LaTeX)
# ---------------------------------------------------------
mapa_cenarios <- tibble::tribble(
  ~cenario, ~n, ~m, ~architecture, ~h2, ~n_qtl, ~epi_strength,
  "C1",   600, 1500, "add_epi",  0.6, 10, 1.2,
  "C2",   600, 2000, "additive", 0.3, 30, 0.0,
  "C3",   600, 2000, "add_epi",  0.3, 10, 0.0,
  "C4",   600, 2000, "add_epi",  0.6, 10, 1.2,
  "C5",   600, 2000, "add_epi",  0.6, 30, 0.0,
  "C6",   600, 2000, "sparse",   0.6, 10, 0.0,
  "C7",   800, 1500, "sparse",   0.3, 30, 0.0,
  "C8",   800, 2000, "add_epi",  0.3, 10, 0.0,
  "C9",   800, 2000, "add_epi",  0.6, 10, 0.0,
  "C10", 1000, 2000, "add_epi",  0.3, 10, 0.0,
  "C11", 1000, 2000, "add_epi",  0.6, 10, 0.0,
  "C12", 1000, 2000, "sparse",   0.6, 30, 0.0
)

# ---------------------------------------------------------
# 3. Juntar cenários aos dados
# ---------------------------------------------------------
gain_summary_cenarios <- gain_summary %>%
  left_join(
    mapa_cenarios,
    by = c("n", "m", "architecture", "h2", "n_qtl", "epi_strength")
  )

# Checagem obrigatória
if (any(is.na(gain_summary_cenarios$cenario))) {
  stop("Erro: existem cenários sem correspondência.")
}

# ---------------------------------------------------------
# 4. Número de repetições
# ---------------------------------------------------------
n_rep <- 100

# ---------------------------------------------------------
# 5. Gráfico - cor_y (ordenado por ganho)
# ---------------------------------------------------------
gain_plot_y <- gain_summary_cenarios %>%
  mutate(
    se  = gain_cor_y_sd / sqrt(n_rep),
    lwr = gain_cor_y_mean - 1.96 * se,
    upr = gain_cor_y_mean + 1.96 * se
  ) %>%
  arrange(gain_cor_y_mean) %>%
  mutate(cenario = factor(cenario, levels = cenario))

grafico_y <- ggplot(gain_plot_y,
                    aes(x = gain_cor_y_mean, y = cenario)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), width = 0.2) +
  geom_point(size = 2.8) +
  labs(
    x = expression(Delta*"cor"[y]*" (MIBLUP - GBLUP)"),
    y = "Cenário",
    title = ""
  ) +
  theme_minimal(base_size = 11)

grafico_y

# ---------------------------------------------------------
# 6. Gráfico - cor_g (ordenado por ganho)
# ---------------------------------------------------------
gain_plot_g <- gain_summary_cenarios %>%
  mutate(
    se  = gain_cor_g_sd / sqrt(n_rep),
    lwr = gain_cor_g_mean - 1.96 * se,
    upr = gain_cor_g_mean + 1.96 * se
  ) %>%
  arrange(gain_cor_g_mean) %>%
  mutate(cenario = factor(cenario, levels = cenario))

grafico_g <- ggplot(gain_plot_g,
                    aes(x = gain_cor_g_mean, y = cenario)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), width = 0.2) +
  geom_point(size = 2.8) +
  labs(
    x = expression(Delta*"cor"[g]*" (MIBLUP - GBLUP)"),
    y = "Cenário",
    title = ""
  ) +
  theme_minimal(base_size = 11)

grafico_g

# ---------------------------------------------------------
# 7. Salvar figuras
# ---------------------------------------------------------
ggsave("Figuras/Ganho_cor_y_cenarios.png",
       grafico_y, width = 8, height = 4.5, dpi = 300)

ggsave("Figuras/Ganho_cor_g_cenarios.png",
       grafico_g, width = 8, height = 4.5, dpi = 300)


