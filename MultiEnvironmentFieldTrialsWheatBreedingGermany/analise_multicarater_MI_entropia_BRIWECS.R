# ==============================================================================
# BRIWECS
# ÍNDICE MULTICARÁTER GUIADO POR INFORMAÇÃO MÚTUA
# E ESTABILIDADE MULTIVARIADA VIA ENTROPIA DE SHANNON
#
# Este script é uma extensão de "analise_MI_BRIWECS.R".
#
# Ele NÃO utiliza abordagem bayesiana. Os valores de entrada são os BLUEs
# frequencistas já estimados para cada:
#
#   Genótipo x Condição de cultivo x Característica
#
# Saídas principais:
#   1) pesos informacionais das características;
#   2) índice multicaráter local em cada ambiente;
#   3) mérito multicaráter médio de cada genótipo;
#   4) entropia de Shannon normalizada como estabilidade informacional;
#   5) identificação do quadrante "superior e estável";
#   6) fronteira de Pareto mérito x estabilidade;
#   7) gráfico semelhante ao usado no repositório anterior.
#
# IMPORTANTE:
#   - A Informação Mútua mede associação, não direção de melhoramento.
#   - A direção desejável de cada característica deve ser definida pelo usuário.
#   - Variáveis derivadas diretamente de produtividade, como Grain e Biomass,
#     foram excluídas para evitar circularidade matemática.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURAÇÕES
# ==============================================================================

if (!exists("PASTA_DADOS")) {
  PASTA_DADOS <- "."
}

if (!exists("PASTA_SAIDAS")) {
  PASTA_SAIDAS <- file.path(
    PASTA_DADOS,
    "resultados_MI_BRIWECS"
  )
}

dir.create(
  PASTA_SAIDAS,
  showWarnings = FALSE,
  recursive = TRUE
)

SEMENTE_MULTITRAIT <- 20260716

# Características incluídas no índice.
#
# Todas são observadas diretamente ou resultam apenas da harmonização de
# medições equivalentes. Foram excluídas:
#   Grain          = função de Seedyield e TGW;
#   Biomass        = função de Seedyield e Harvest_Index;
#   Straw          = função de Biomass e Harvest_Index;
#   Harvest_Index  = parcialmente reconstruído a partir de yield/biomass;
#   Grain_per_spike = função de yield, TGW e número de espigas.
TRAITS_INDEX <- c(
  "Seedyield",
  "Crude_protein",
  "TGW",
  "Spike_number"
)

# Direção de seleção:
#    1 = valores maiores são desejáveis
#   -1 = valores menores são desejáveis
#
# A MI não define essas direções.
DIRECOES <- c(
  Seedyield = 1,
  Crude_protein = 1,
  TGW = 1,
  Spike_number = 1
)

# Peso explícito da produtividade no índice final.
# Os demais pesos compartilham o restante.
PESO_PRODUTIVIDADE <- 0.50

# "MI"   = pesos secundários proporcionais somente à relevância com yield.
# "mRMR" = relevância com yield penalizada pela redundância entre secundárias.
# A opção mRMR é recomendada para evitar premiar características repetitivas.
MODO_PESOS <- "mRMR"

# Discretização para estimar Informação Mútua.
N_BINS_PESOS <- 6

# Permutações para correção de viés da NMI.
N_PERM_PESOS <- 499

# Estados de desempenho usados na entropia.
# Com 4 estados, cada genótipo fica classificado em quartis dentro do ambiente.
N_ESTADOS_ENTROPIA <- 4

# Um índice local é aceito quando pelo menos esta proporção do peso total
# está representada pelas características disponíveis.
COBERTURA_MINIMA_PESO <- 0.75

# Um genótipo precisa estar presente em pelo menos esta fração das condições
# para entrar na seleção final.
FRACAO_MINIMA_AMBIENTES <- 0.80

# Número máximo de genótipos do quadrante superior/estável a rotular.
N_ROTULOS_QUADRANTE <- 20


# ==============================================================================
# 1. PACOTES
# ==============================================================================

PACOTES_MULTITRAIT <- c(
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "tibble",
  "ggplot2",
  "ggrepel",
  "stringr"
)

pacotes_ausentes <- setdiff(
  PACOTES_MULTITRAIT,
  rownames(
    installed.packages()
  )
)

if (length(pacotes_ausentes) > 0) {
  install.packages(
    pacotes_ausentes
  )
}

invisible(
  lapply(
    PACOTES_MULTITRAIT,
    library,
    character.only = TRUE
  )
)

options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

set.seed(
  SEMENTE_MULTITRAIT
)


# ==============================================================================
# 2. FUNÇÕES AUXILIARES
# ==============================================================================

titulo_console <- function(texto) {
  cat(
    "\n",
    paste(
      rep("=", 80),
      collapse = ""
    ),
    "\n",
    texto,
    "\n",
    paste(
      rep("=", 80),
      collapse = ""
    ),
    "\n",
    sep = ""
  )
}


discretizar_por_postos_multi <- function(
  x,
  nbins = 4
) {

  out <- rep(
    NA_integer_,
    length(x)
  )

  ok <- is.finite(x)

  if (
    sum(ok) < nbins ||
      length(unique(x[ok])) < 2
  ) {
    return(out)
  }

  postos <- rank(
    x[ok],
    ties.method = "average"
  )

  classes <- ceiling(
    nbins *
      postos /
      length(postos)
  )

  classes <- pmin(
    nbins,
    pmax(
      1,
      classes
    )
  )

  out[ok] <- as.integer(
    classes
  )

  out
}


entropia_discreta_multi <- function(x) {

  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0) {
    return(NA_real_)
  }

  p <- prop.table(
    table(x)
  )

  -sum(
    p * log(p)
  )
}


mi_discreta_multi <- function(
  x,
  y
) {

  ok <- complete.cases(
    x,
    y
  )

  x <- x[ok]
  y <- y[ok]

  if (
    length(x) == 0 ||
      length(unique(x)) < 2 ||
      length(unique(y)) < 2
  ) {
    return(NA_real_)
  }

  tab <- table(
    x,
    y
  )

  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)

  idx <- which(
    pxy > 0,
    arr.ind = TRUE
  )

  sum(
    pxy[idx] *
      log(
        pxy[idx] /
          (
            px[idx[, 1]] *
              py[idx[, 2]]
          )
      )
  )
}


nmi_discreta_multi <- function(
  x,
  y
) {

  mi <- mi_discreta_multi(
    x,
    y
  )

  hx <- entropia_discreta_multi(x)
  hy <- entropia_discreta_multi(y)

  if (
    is.na(mi) ||
      is.na(hx) ||
      is.na(hy) ||
      (hx + hy) == 0
  ) {
    return(NA_real_)
  }

  2 * mi / (hx + hy)
}


calcular_nmi_corrigida <- function(
  x,
  y,
  nbins = 6,
  n_perm = 499
) {

  ok <- complete.cases(
    x,
    y
  ) &
    is.finite(x) &
    is.finite(y)

  x <- as.numeric(
    x[ok]
  )

  y <- as.numeric(
    y[ok]
  )

  n <- length(x)

  if (
    n < 30 ||
      length(unique(x)) < 2 ||
      length(unique(y)) < 2
  ) {
    return(
      tibble::tibble(
        n = n,
        NMI_observada = NA_real_,
        NMI_nula_media = NA_real_,
        NMI_corrigida = NA_real_,
        p_permutacao = NA_real_
      )
    )
  }

  xd <- discretizar_por_postos_multi(
    x,
    nbins = nbins
  )

  yd <- discretizar_por_postos_multi(
    y,
    nbins = nbins
  )

  nmi_obs <- nmi_discreta_multi(
    xd,
    yd
  )

  nmi_perm <- replicate(
    n_perm,
    nmi_discreta_multi(
      xd,
      sample(
        yd,
        replace = FALSE
      )
    )
  )

  media_nula <- mean(
    nmi_perm,
    na.rm = TRUE
  )

  p_perm <- (
    1 +
      sum(
        nmi_perm >= nmi_obs,
        na.rm = TRUE
      )
  ) / (
    n_perm + 1
  )

  tibble::tibble(
    n = n,
    NMI_observada = nmi_obs,
    NMI_nula_media = media_nula,
    NMI_corrigida = max(
      0,
      nmi_obs - media_nula
    ),
    p_permutacao = p_perm
  )
}


normalizar_vetor <- function(x) {

  x <- as.numeric(x)

  if (
    all(!is.finite(x)) ||
      sum(
        x,
        na.rm = TRUE
      ) <= 0
  ) {
    return(
      rep(
        1 / length(x),
        length(x)
      )
    )
  }

  x[
    !is.finite(x)
  ] <- 0

  x / sum(x)
}


calcular_pareto <- function(
  merito,
  entropia
) {

  n <- length(merito)

  sapply(
    seq_len(n),
    function(i) {

      dominado <- (
        merito >= merito[i] &
          entropia <= entropia[i] &
          (
            merito > merito[i] |
              entropia < entropia[i]
          )
      )

      !any(
        dominado,
        na.rm = TRUE
      )
    }
  )
}


# ==============================================================================
# 3. LEITURA DOS BLUEs
# ==============================================================================

titulo_console(
  "1. LEITURA DOS BLUEs"
)

if (!exists("blues")) {

  arquivo_blues <- file.path(
    PASTA_SAIDAS,
    "03_BLUEs_45_condicoes_9_traits.csv"
  )

  if (!file.exists(arquivo_blues)) {

    stop(
      "O objeto 'blues' não está na memória e o arquivo:\n",
      arquivo_blues,
      "\n não foi encontrado.\n",
      "Execute primeiro o script analise_MI_BRIWECS.R."
    )
  }

  blues <- readr::read_csv(
    arquivo_blues,
    show_col_types = FALSE
  )
}

cat(
  "Dimensão da tabela de BLUEs:\n"
)

print(
  dim(blues)
)

traits_ausentes <- setdiff(
  TRAITS_INDEX,
  unique(
    blues$Trait
  )
)

if (length(traits_ausentes) > 0) {
  stop(
    "As seguintes características não estão disponíveis em 'blues': ",
    paste(
      traits_ausentes,
      collapse = ", "
    )
  )
}

if ("status_modelo" %in% names(blues)) {

  cat(
    "\nStatus dos modelos usados na análise multicaráter:\n"
  )

  print(
    blues |>
      dplyr::filter(
        Trait %in% TRAITS_INDEX
      ) |>
      dplyr::distinct(
        Condition,
        Trait,
        status_modelo
      ) |>
      dplyr::count(
        status_modelo,
        name = "n_modelos"
      )
  )
}


# ==============================================================================
# 4. PREPARAÇÃO DOS BLUEs MULTICARÁTER
# ==============================================================================

titulo_console(
  "2. PREPARAÇÃO DAS CARACTERÍSTICAS"
)

tabela_direcoes <- tibble::tibble(
  Trait = names(DIRECOES),
  Direcao = as.numeric(
    DIRECOES
  )
)

blues_index <- blues |>
  dplyr::filter(
    Trait %in% TRAITS_INDEX,
    is.finite(BLUE)
  ) |>
  dplyr::inner_join(
    tabela_direcoes,
    by = "Trait"
  )

disponibilidade <- blues_index |>
  dplyr::group_by(
    Trait
  ) |>
  dplyr::summarise(
    n_BLUEs = dplyr::n(),
    n_genotipos = dplyr::n_distinct(
      BRISONr
    ),
    n_condicoes = dplyr::n_distinct(
      Condition
    ),
    .groups = "drop"
  )

cat(
  "Disponibilidade das características:\n"
)

print(
  disponibilidade,
  n = nrow(disponibilidade)
)


# ==============================================================================
# 5. MÉDIAS GLOBAIS DOS GENÓTIPOS
# ==============================================================================

titulo_console(
  "3. MÉDIAS GLOBAIS DOS GENÓTIPOS"
)

medias_globais <- blues_index |>
  dplyr::group_by(
    BRISONr,
    Trait
  ) |>
  dplyr::summarise(
    BLUE_medio = mean(
      BLUE,
      na.rm = TRUE
    ),
    n_ambientes_trait = dplyr::n_distinct(
      Condition
    ),
    .groups = "drop"
  )

medias_wide <- medias_globais |>
  dplyr::select(
    BRISONr,
    Trait,
    BLUE_medio
  ) |>
  tidyr::pivot_wider(
    names_from = Trait,
    values_from = BLUE_medio
  )

dados_pesos <- medias_wide |>
  dplyr::filter(
    complete.cases(
      dplyr::across(
        dplyr::all_of(
          TRAITS_INDEX
        )
      )
    )
  )

cat(
  "Genótipos completos para estimar os pesos:",
  nrow(dados_pesos),
  "\n"
)


# ==============================================================================
# 6. PESOS INFORMACIONAIS
# ==============================================================================

titulo_console(
  "4. PESOS INFORMACIONAIS DAS CARACTERÍSTICAS"
)

traits_secundarias <- setdiff(
  TRAITS_INDEX,
  "Seedyield"
)

set.seed(
  SEMENTE_MULTITRAIT
)

relevancia_lista <- purrr::map(
  traits_secundarias,
  function(trait_i) {

    calcular_nmi_corrigida(
      x = dados_pesos$Seedyield,
      y = dados_pesos[[trait_i]],
      nbins = N_BINS_PESOS,
      n_perm = N_PERM_PESOS
    ) |>
      dplyr::mutate(
        Trait = trait_i
      )
  }
)

relevancia <- dplyr::bind_rows(
  relevancia_lista
) |>
  dplyr::select(
    Trait,
    dplyr::everything()
  )


# Redundância média de cada característica secundária em relação às outras.
if (
  MODO_PESOS == "mRMR" &&
    length(traits_secundarias) > 1
) {

  pares_secundarias <- t(
    utils::combn(
      traits_secundarias,
      2
    )
  )

  redundancia_pares <- purrr::map_dfr(
    seq_len(
      nrow(pares_secundarias)
    ),
    function(i) {

      t1 <- pares_secundarias[i, 1]
      t2 <- pares_secundarias[i, 2]

      calcular_nmi_corrigida(
        x = dados_pesos[[t1]],
        y = dados_pesos[[t2]],
        nbins = N_BINS_PESOS,
        n_perm = N_PERM_PESOS
      ) |>
        dplyr::transmute(
          Trait_1 = t1,
          Trait_2 = t2,
          NMI_redundancia = NMI_corrigida
        )
    }
  )

  redundancia <- dplyr::bind_rows(
    redundancia_pares |>
      dplyr::transmute(
        Trait = Trait_1,
        NMI_redundancia
      ),

    redundancia_pares |>
      dplyr::transmute(
        Trait = Trait_2,
        NMI_redundancia
      )
  ) |>
    dplyr::group_by(
      Trait
    ) |>
    dplyr::summarise(
      Redundancia_media = mean(
        NMI_redundancia,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

} else {

  redundancia <- tibble::tibble(
    Trait = traits_secundarias,
    Redundancia_media = 0
  )
}


pesos_secundarios <- relevancia |>
  dplyr::left_join(
    redundancia,
    by = "Trait"
  ) |>
  dplyr::mutate(
    Escore_informacional = dplyr::case_when(
      MODO_PESOS == "mRMR" ~
        NMI_corrigida /
        (
          1 +
            Redundancia_media
        ),

      TRUE ~ NMI_corrigida
    ),

    Peso_relativo_secundario = normalizar_vetor(
      Escore_informacional
    ),

    Peso_final = (
      1 -
        PESO_PRODUTIVIDADE
    ) *
      Peso_relativo_secundario
  )


pesos_finais <- dplyr::bind_rows(
  tibble::tibble(
    Trait = "Seedyield",
    n = nrow(dados_pesos),
    NMI_observada = NA_real_,
    NMI_nula_media = NA_real_,
    NMI_corrigida = NA_real_,
    p_permutacao = NA_real_,
    Redundancia_media = NA_real_,
    Escore_informacional = NA_real_,
    Peso_relativo_secundario = NA_real_,
    Peso_final = PESO_PRODUTIVIDADE
  ),

  pesos_secundarios
) |>
  dplyr::left_join(
    tabela_direcoes,
    by = "Trait"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      Peso_final
    )
  )


cat(
  "Modo de ponderação:",
  MODO_PESOS,
  "\n\n"
)

cat(
  "Pesos finais do índice:\n"
)

print(
  pesos_finais |>
    dplyr::select(
      Trait,
      Direcao,
      NMI_corrigida,
      Redundancia_media,
      Peso_final
    ),
  n = nrow(pesos_finais)
)

cat(
  "\nSoma dos pesos:",
  sum(
    pesos_finais$Peso_final
  ),
  "\n"
)

readr::write_csv(
  pesos_finais,
  file.path(
    PASTA_SAIDAS,
    "19_pesos_indice_multicarater_MI.csv"
  )
)


# ==============================================================================
# 7. ÍNDICE MULTICARÁTER LOCAL EM CADA CONDIÇÃO
# ==============================================================================

titulo_console(
  "5. ÍNDICE MULTICARÁTER LOCAL"
)

# Padronização dentro de cada ambiente e característica.
# Assim, o índice representa a posição relativa do genótipo naquele ambiente.
blues_padronizados <- blues_index |>
  dplyr::group_by(
    Condition,
    Trait
  ) |>
  dplyr::mutate(
    media_condicao_trait = mean(
      BLUE,
      na.rm = TRUE
    ),

    dp_condicao_trait = stats::sd(
      BLUE,
      na.rm = TRUE
    ),

    Z_trait = dplyr::if_else(
      is.finite(dp_condicao_trait) &
        dp_condicao_trait > 0,

      (
        BLUE -
          media_condicao_trait
      ) /
        dp_condicao_trait,

      NA_real_
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    pesos_finais |>
      dplyr::select(
        Trait,
        Peso_final,
        Direcao
      ),
    by = c(
      "Trait",
      "Direcao"
    )
  ) |>
  dplyr::mutate(
    Z_desejavel = Z_trait *
      Direcao
  )


indice_local <- blues_padronizados |>
  dplyr::group_by(
    BRISONr,
    Condition
  ) |>
  dplyr::summarise(
    tem_yield = any(
      Trait == "Seedyield" &
        is.finite(
          Z_desejavel
        )
    ),

    peso_disponivel = sum(
      Peso_final[
        is.finite(
          Z_desejavel
        )
      ],
      na.rm = TRUE
    ),

    n_traits_disponiveis = sum(
      is.finite(
        Z_desejavel
      )
    ),

    MI_Index_Local = dplyr::if_else(
      tem_yield &
        peso_disponivel >= COBERTURA_MINIMA_PESO,

      sum(
        Z_desejavel *
          Peso_final,
        na.rm = TRUE
      ) /
        peso_disponivel,

      NA_real_
    ),

    .groups = "drop"
  ) |>
  dplyr::filter(
    is.finite(
      MI_Index_Local
    )
  )


cat(
  "Número de observações Genótipo × Condição com índice válido:",
  nrow(indice_local),
  "\n"
)

cat(
  "Número de genótipos:",
  dplyr::n_distinct(
    indice_local$BRISONr
  ),
  "\n"
)

cat(
  "Número de condições:",
  dplyr::n_distinct(
    indice_local$Condition
  ),
  "\n"
)

readr::write_csv(
  indice_local,
  file.path(
    PASTA_SAIDAS,
    "20_indice_multicarater_local.csv"
  )
)


# ==============================================================================
# 8. ESTADOS DE DESEMPENHO E ENTROPIA DE SHANNON
# ==============================================================================

titulo_console(
  "6. ESTABILIDADE MULTIVARIADA VIA ENTROPIA"
)

# Cada genótipo é classificado em um quartil do índice multicaráter dentro
# de cada condição.
indice_estados <- indice_local |>
  dplyr::group_by(
    Condition
  ) |>
  dplyr::mutate(
    Estado_multicarater = discretizar_por_postos_multi(
      MI_Index_Local,
      nbins = N_ESTADOS_ENTROPIA
    )
  ) |>
  dplyr::ungroup()


n_condicoes_total <- dplyr::n_distinct(
  indice_estados$Condition
)

resultado_estabilidade <- indice_estados |>
  dplyr::group_by(
    BRISONr
  ) |>
  dplyr::summarise(
    Merito_multicarater = mean(
      MI_Index_Local,
      na.rm = TRUE
    ),

    DP_indice_local = stats::sd(
      MI_Index_Local,
      na.rm = TRUE
    ),

    n_ambientes = dplyr::n_distinct(
      Condition
    ),

    Fracao_ambientes = n_ambientes /
      n_condicoes_total,

    Entropia_Shannon = entropia_discreta_multi(
      Estado_multicarater
    ),

    Entropia_normalizada = Entropia_Shannon /
      log(
        N_ESTADOS_ENTROPIA
      ),

    Fracao_quartil_superior = mean(
      Estado_multicarater ==
        N_ESTADOS_ENTROPIA,
      na.rm = TRUE
    ),

    Fracao_quartil_inferior = mean(
      Estado_multicarater == 1,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) |>
  dplyr::filter(
    Fracao_ambientes >=
      FRACAO_MINIMA_AMBIENTES,
    is.finite(
      Merito_multicarater
    ),
    is.finite(
      Entropia_normalizada
    )
  )


# ==============================================================================
# 9. METADADOS DOS CULTIVARES
# ==============================================================================

arquivo_cultivares <- file.path(
  PASTA_DADOS,
  "BRIWECs_cultivar_info.csv"
)

if (
  exists("cultivares") &&
    all(
      c(
        "BRISONr",
        "genotype"
      ) %in% names(cultivares)
    )
) {

  info_cultivares <- cultivares

} else if (
  file.exists(
    arquivo_cultivares
  )
) {

  info_cultivares <- readr::read_csv(
    arquivo_cultivares,
    locale = readr::locale(
      encoding = "ISO-8859-1"
    ),
    show_col_types = FALSE
  )

} else {

  info_cultivares <- tibble::tibble(
    BRISONr = resultado_estabilidade$BRISONr,
    genotype = resultado_estabilidade$BRISONr,
    RYear = NA_real_
  )
}


resultado_final_multitraits <- resultado_estabilidade |>
  dplyr::left_join(
    info_cultivares |>
      dplyr::select(
        dplyr::any_of(
          c(
            "BRISONr",
            "genotype",
            "RYear",
            "baking_qulaity",
            "country",
            "breeder"
          )
        )
      ),
    by = "BRISONr"
  ) |>
  dplyr::mutate(
    Genotipo = dplyr::coalesce(
      genotype,
      BRISONr
    )
  )


# ==============================================================================
# 10. QUADRANTES, PARETO E RANKING
# ==============================================================================

titulo_console(
  "7. CLASSIFICAÇÃO DOS GENÓTIPOS"
)

limite_merito <- median(
  resultado_final_multitraits$Merito_multicarater,
  na.rm = TRUE
)

limite_entropia <- median(
  resultado_final_multitraits$Entropia_normalizada,
  na.rm = TRUE
)


resultado_final_multitraits <- resultado_final_multitraits |>
  dplyr::mutate(
    Quadrante = dplyr::case_when(
      Merito_multicarater >= limite_merito &
        Entropia_normalizada <= limite_entropia ~
        "Superior e estável",

      Merito_multicarater >= limite_merito &
        Entropia_normalizada > limite_entropia ~
        "Superior e instável",

      Merito_multicarater < limite_merito &
        Entropia_normalizada <= limite_entropia ~
        "Inferior e estável",

      TRUE ~
        "Inferior e instável"
    )
  )


resultado_final_multitraits$Pareto <- calcular_pareto(
  merito = resultado_final_multitraits$Merito_multicarater,
  entropia = resultado_final_multitraits$Entropia_normalizada
)


resultado_final_multitraits <- resultado_final_multitraits |>
  dplyr::mutate(
    Z_merito = as.numeric(
      scale(
        Merito_multicarater
      )
    ),

    Z_estabilidade = -as.numeric(
      scale(
        Entropia_normalizada
      )
    ),

    Escore_integrado_exploratorio =
      Z_merito +
      Z_estabilidade
  ) |>
  dplyr::arrange(
    dplyr::desc(
      Escore_integrado_exploratorio
    )
  ) |>
  dplyr::mutate(
    Rank_exploratorio = dplyr::row_number()
  )


cat(
  "Limite de mérito utilizado:",
  round(
    limite_merito,
    4
  ),
  "\n"
)

cat(
  "Limite de entropia utilizado:",
  round(
    limite_entropia,
    4
  ),
  "\n\n"
)


cat(
  "Número de genótipos por quadrante:\n"
)

print(
  resultado_final_multitraits |>
    dplyr::count(
      Quadrante,
      name = "n_genotipos"
    ) |>
    dplyr::arrange(
      dplyr::desc(
        n_genotipos
      )
    ),
  n = 4
)


cat(
  "\n20 genótipos do quadrante superior e estável com maior mérito:\n"
)

top_superiores_estaveis <- resultado_final_multitraits |>
  dplyr::filter(
    Quadrante ==
      "Superior e estável"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      Merito_multicarater
    ),
    Entropia_normalizada
  ) |>
  dplyr::select(
    Genotipo,
    BRISONr,
    Merito_multicarater,
    Entropia_normalizada,
    Fracao_quartil_superior,
    n_ambientes,
    Pareto,
    RYear
  ) |>
  dplyr::slice_head(
    n = 20
  )

print(
  top_superiores_estaveis,
  n = 20
)


cat(
  "\nGenótipos na fronteira de Pareto:\n"
)

tabela_pareto <- resultado_final_multitraits |>
  dplyr::filter(
    Pareto
  ) |>
  dplyr::arrange(
    Entropia_normalizada,
    dplyr::desc(
      Merito_multicarater
    )
  ) |>
  dplyr::select(
    Genotipo,
    BRISONr,
    Merito_multicarater,
    Entropia_normalizada,
    Fracao_quartil_superior,
    n_ambientes,
    RYear
  )

print(
  tabela_pareto,
  n = nrow(
    tabela_pareto
  )
)


cat(
  "\n20 maiores escores integrados exploratórios:\n"
)

print(
  resultado_final_multitraits |>
    dplyr::select(
      Genotipo,
      BRISONr,
      Merito_multicarater,
      Entropia_normalizada,
      Fracao_quartil_superior,
      Quadrante,
      Pareto,
      Escore_integrado_exploratorio,
      Rank_exploratorio
    ) |>
    dplyr::slice_head(
      n = 20
    ),
  n = 20
)


readr::write_csv(
  resultado_final_multitraits,
  file.path(
    PASTA_SAIDAS,
    "21_selecao_multicarater_entropia_genotipos.csv"
  )
)

readr::write_csv(
  top_superiores_estaveis,
  file.path(
    PASTA_SAIDAS,
    "22_top_superiores_estaveis.csv"
  )
)

readr::write_csv(
  tabela_pareto,
  file.path(
    PASTA_SAIDAS,
    "23_fronteira_pareto_merito_estabilidade.csv"
  )
)


# ==============================================================================
# 11. DEFINIÇÃO DOS RÓTULOS DO GRÁFICO
# ==============================================================================

rotulos_quadrante <- resultado_final_multitraits |>
  dplyr::filter(
    Quadrante ==
      "Superior e estável"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      Merito_multicarater
    ),
    Entropia_normalizada
  ) |>
  dplyr::slice_head(
    n = N_ROTULOS_QUADRANTE
  ) |>
  dplyr::pull(
    BRISONr
  )


resultado_plot <- resultado_final_multitraits |>
  dplyr::mutate(
    Rotular = Pareto |
      BRISONr %in%
      rotulos_quadrante,

    Destaque = dplyr::case_when(
      Pareto ~
        "Fronteira de Pareto",

      Quadrante ==
        "Superior e estável" ~
        "Superior e estável",

      TRUE ~
        "Demais genótipos"
    )
  )


# ==============================================================================
# 12. GRÁFICO PRINCIPAL
# ==============================================================================

titulo_console(
  "8. GRÁFICO MÉRITO MULTICARÁTER × ESTABILIDADE"
)

grafico_multitraits_entropia <- ggplot(
  resultado_plot,
  aes(
    x = Entropia_normalizada,
    y = Merito_multicarater
  )
) +
  geom_point(
    aes(
      color = Destaque
    ),
    size = 2.6,
    alpha = 0.75
  ) +
  geom_vline(
    xintercept = limite_entropia,
    linetype = "dashed",
    color = "red",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = limite_merito,
    linetype = "dashed",
    color = "red",
    linewidth = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = resultado_plot |>
      dplyr::filter(
        Rotular
      ),
    aes(
      label = Genotipo
    ),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      "Fronteira de Pareto" = "#006400",
      "Superior e estável" = "#5AAE61",
      "Demais genótipos" = "grey65"
    )
  ) +
  coord_cartesian(
    xlim = c(
      0,
      1
    )
  ) +
  labs(
    title = "Mérito multicaráter × estabilidade informacional",
    subtitle = paste0(
      "BLUEs frequencistas; pesos ",
      MODO_PESOS,
      "; linhas tracejadas representam as medianas"
    ),
    x = paste0(
      "Entropia de Shannon normalizada\n",
      "(mais próximo de zero = maior estabilidade)"
    ),
    y = paste0(
      "Índice multicaráter informacional\n",
      "(maior = melhor desempenho global)"
    ),
    color = NULL,
    caption = paste0(
      "Quadrante superior esquerdo: genótipos superiores e estáveis. ",
      "Pontos verde-escuros: fronteira de Pareto."
    )
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


print(
  grafico_multitraits_entropia
)


ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig08_merito_multicarater_vs_entropia.png"
  ),
  plot = grafico_multitraits_entropia,
  width = 13,
  height = 9,
  dpi = 300
)


# ==============================================================================
# 13. GRÁFICO COMPLEMENTAR:
#     FREQUÊNCIA NO QUARTIL SUPERIOR × ENTROPIA
# ==============================================================================

grafico_quartil_superior <- ggplot(
  resultado_plot,
  aes(
    x = Entropia_normalizada,
    y = Fracao_quartil_superior
  )
) +
  geom_point(
    aes(
      color = Destaque
    ),
    size = 2.6,
    alpha = 0.75
  ) +
  ggrepel::geom_text_repel(
    data = resultado_plot |>
      dplyr::filter(
        Rotular
      ),
    aes(
      label = Genotipo
    ),
    size = 3.1,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      "Fronteira de Pareto" = "#006400",
      "Superior e estável" = "#5AAE61",
      "Demais genótipos" = "grey65"
    )
  ) +
  coord_cartesian(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +
  labs(
    title = "Estabilidade e frequência de desempenho superior",
    subtitle = "A frequência no quartil superior evita confundir estabilidade boa com estabilidade ruim",
    x = "Entropia de Shannon normalizada",
    y = "Proporção de ambientes no quartil superior",
    color = NULL
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


print(
  grafico_quartil_superior
)


ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig09_entropia_vs_frequencia_quartil_superior.png"
  ),
  plot = grafico_quartil_superior,
  width = 13,
  height = 9,
  dpi = 300
)


# ==============================================================================
# 14. RESUMO FINAL
# ==============================================================================

titulo_console(
  "9. RESUMO FINAL"
)

cat(
  "Análise concluída.\n\n"
)

cat(
  "Interpretação do gráfico principal:\n",
  "  - esquerda: menor entropia, portanto maior estabilidade de classificação;\n",
  "  - direita: maior instabilidade entre os ambientes;\n",
  "  - acima: maior mérito multicaráter médio;\n",
  "  - abaixo: menor mérito multicaráter médio;\n",
  "  - superior esquerdo: região desejável;\n",
  "  - verde-escuro: genótipos não dominados na fronteira de Pareto.\n",
  sep = ""
)

cat(
  "\nSaídas principais:\n",
  "  19_pesos_indice_multicarater_MI.csv\n",
  "  20_indice_multicarater_local.csv\n",
  "  21_selecao_multicarater_entropia_genotipos.csv\n",
  "  22_top_superiores_estaveis.csv\n",
  "  23_fronteira_pareto_merito_estabilidade.csv\n",
  "  fig08_merito_multicarater_vs_entropia.png\n",
  "  fig09_entropia_vs_frequencia_quartil_superior.png\n",
  sep = ""
)

cat(
  "\nCuidado de interpretação:\n",
  "Baixa entropia significa comportamento previsível, não necessariamente bom.\n",
  "Por isso, a entropia deve ser interpretada em conjunto com o mérito e com a\n",
  "frequência de presença no quartil superior.\n",
  sep = ""
)

cat(
  "\nFim do script.\n"
)
