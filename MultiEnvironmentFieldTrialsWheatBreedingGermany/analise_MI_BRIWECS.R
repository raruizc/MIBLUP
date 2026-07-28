# ==============================================================================
# BRIWECS - Informação Mútua em Ensaios Multiambiente de Trigo
#
# Objetivo principal:
#   Comparar a consistência do desempenho genotípico entre condições de cultivo
#   (Ano x Local x Manejo) usando:
#     - R² associado à regressão SMA (equivalente a r de Pearson ao quadrado
#       no caso bivariado)
#     - correlação de Spearman
#     - Informação Mútua Normalizada (NMI)
#     - correção do viés da NMI por permutação
#
# Objetivos complementares:
#   1) Identificar pares de ambientes com dependência informacional relativamente
#      alta, mas associação linear relativamente baixa.
#   2) Comparar a consistência de nove características agronômicas.
#   3) Avaliar, dentro de cada ambiente, a associação de outras características
#      com a produtividade.
#   4) Auditar arquivos de clima, solo, fertilização e proteção de plantas.
#
# IMPORTANTE:
#   - A análise principal segue o subconjunto balanceado usado na validação
#     técnica da publicação: 2015-2017, GGE/HAN/KAL/KIE/QLB e os manejos
#     HN_WF_RF, HN_NF_RF e LN_NF_RF.
#   - Os arquivos de correção GG2019.xlsx e RHH2016.xlsx do repositório oficial
#     não foram fornecidos. Por isso, os anos/locais afetados não entram na
#     análise inferencial principal.
#   - O script imprime os principais resultados no console E também salva
#     tabelas e figuras.
#
# Como usar:
#   1) Coloque este script na pasta que contém os arquivos enviados.
#   2) Ajuste PASTA_DADOS, se necessário.
#   3) Execute o script inteiro.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURAÇÕES
# ==============================================================================

PASTA_DADOS  <- "."   # altere se os arquivos estiverem em outra pasta
PASTA_SAIDAS <- file.path(PASTA_DADOS, "resultados_MI_BRIWECS")

dir.create(PASTA_SAIDAS, showWarnings = FALSE, recursive = TRUE)

SEMENTE <- 20260713

# Número de permutações para a análise principal de produtividade.
# Para uma análise final mais rigorosa, pode-se aumentar para 999 ou 1999.
N_PERM_MI <- 999

# Permutações para relações característica x produtividade.
N_PERM_TRAIT_YIELD <- 199

# Número principal de classes para discretização por postos.
# A sensibilidade também será verificada com 4, 6 e 8 classes.
N_BINS_MI <- 6
BINS_SENSIBILIDADE <- c(4, 6, 8)

# Número mínimo de cultivares em comum para calcular associação.
MIN_N_PAR <- 30

# Características usadas na extensão da análise da publicação.
TRAITS_ANALISE <- c(
  "Seedyield",
  "Harvest_Index",
  "Straw",
  "Biomass",
  "Grain",
  "Crude_protein",
  "Grain_per_spike",
  "Spike_number",
  "TGW"
)

ABREV_TRAITS <- c(
  Seedyield       = "GY",
  Harvest_Index   = "HI",
  Straw           = "Straw",
  Biomass         = "SDM",
  Grain           = "GN",
  Crude_protein   = "GP",
  Grain_per_spike = "GpS",
  Spike_number    = "SN",
  TGW             = "TGW"
)


# ==============================================================================
# 1. PACOTES
# ==============================================================================

PACOTES <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "readr",
  "readxl",
  "stringr",
  "ggplot2",
  "lme4",
  "tibble"
)

faltantes <- setdiff(PACOTES, rownames(installed.packages()))

if (length(faltantes) > 0) {
  cat("\nInstalando pacotes ausentes:\n")
  print(faltantes)
  install.packages(faltantes)
}

invisible(
  lapply(
    PACOTES,
    library,
    character.only = TRUE
  )
)

options(
  dplyr.summarise.inform = FALSE,
  scipen = 999
)

set.seed(SEMENTE)


# ==============================================================================
# 2. FUNÇÕES AUXILIARES GERAIS
# ==============================================================================

titulo_console <- function(texto) {
  cat(
    "\n",
    paste(rep("=", 80), collapse = ""),
    "\n",
    texto,
    "\n",
    paste(rep("=", 80), collapse = ""),
    "\n",
    sep = ""
  )
}


numero_seguro <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "-", ".")] <- NA_character_

  # Os arquivos misturam separadores decimais em alguns pontos.
  x <- gsub(",", ".", x, fixed = TRUE)

  suppressWarnings(as.numeric(x))
}


media_sem_nan <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}


mediana_sem_nan <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(x, na.rm = TRUE)
}


detectar_delimitador <- function(arquivo) {
  linha <- readLines(
    arquivo,
    n = 1,
    warn = FALSE,
    encoding = "latin1"
  )

  n_pv <- stringr::str_count(linha, fixed(";"))
  n_vg <- stringr::str_count(linha, fixed(","))

  ifelse(n_pv > n_vg, ";", ",")
}


ler_csv_flexivel <- function(arquivo) {
  delim <- detectar_delimitador(arquivo)

  readr::read_delim(
    file = arquivo,
    delim = delim,
    col_types = cols(.default = col_character()),
    trim_ws = TRUE,
    show_col_types = FALSE,
    name_repair = "unique",
    locale = locale(encoding = "ISO-8859-1")
  )
}


padronizar_local <- function(x) {
  dplyr::case_when(
    x == "DKI" ~ "KIE",
    x == "KAD" ~ "KAL",
    TRUE ~ x
  )
}


# Regra presente no script público de pré-processamento:
#   - DKI representa o ensaio rain-out shelter em Kiel.
#   - os sufixos _D e _DD são retirados;
#   - condições sem IR ou RO recebem RF;
#   - DKI é posteriormente renomeado para KIE.
#
# Essa regra reproduz 100 condições Ano x Local x Manejo nos arquivos enviados
# e gera exatamente as 45 condições do subconjunto balanceado da publicação.
padronizar_tratamento <- function(Treatment, Location) {
  out <- Treatment

  out <- ifelse(
    Location == "DKI",
    paste0(out, "_RO"),
    out
  )

  out <- stringr::str_replace(out, "_D{1,2}$", "")

  out <- ifelse(
    !stringr::str_detect(out, "(IR|RO)$"),
    paste0(out, "_RF"),
    out
  )

  out
}


limpar_sufixo_repetido <- function(x) {
  x |>
    stringr::str_replace("_IR_IR$", "_IR") |>
    stringr::str_replace("_RF_RF$", "_RF") |>
    stringr::str_replace("_RO_RO$", "_RO")
}


# ==============================================================================
# 3. LOCALIZAÇÃO DOS ARQUIVOS
# ==============================================================================

titulo_console("1. LOCALIZAÇÃO E AUDITORIA DOS ARQUIVOS")

arquivos_talhao <- list.files(
  PASTA_DADOS,
  pattern = "^(DKI|GGE|HAN|KAL|KIE|QLB|RHH)_20[0-9]{2}.*\\.csv$",
  full.names = TRUE
)

arquivos_clima <- list.files(
  PASTA_DADOS,
  pattern = "^DWD_Intpol_.*\\.csv$",
  full.names = TRUE
)

arquivos_manejo_brutos <- list.files(
  PASTA_DADOS,
  pattern = "^Management_information_.*\\.xlsx$",
  full.names = TRUE
)

cat("Arquivos fenotípicos encontrados:", length(arquivos_talhao), "\n")
cat("Arquivos meteorológicos encontrados:", length(arquivos_clima), "\n")
cat("Arquivos brutos de manejo encontrados:", length(arquivos_manejo_brutos), "\n\n")

if (length(arquivos_talhao) == 0) {
  stop(
    "Nenhum arquivo fenotípico foi encontrado. ",
    "Verifique o caminho definido em PASTA_DADOS."
  )
}

cat("Primeiros arquivos fenotípicos:\n")
print(head(basename(arquivos_talhao), 10))


# ==============================================================================
# 4. IMPORTAÇÃO E HARMONIZAÇÃO DOS DADOS FENOTÍPICOS
# ==============================================================================

titulo_console("2. IMPORTAÇÃO E HARMONIZAÇÃO DOS DADOS FENOTÍPICOS")

lista_fenotipos <- purrr::map(
  arquivos_talhao,
  function(arq) {

    dat <- ler_csv_flexivel(arq)

    # Harmoniza Sedimentation/Sedi.
    if ("Sedi" %in% names(dat)) {

      if (!"Sedimentation" %in% names(dat)) {
        dat <- dat |>
          dplyr::rename(Sedimentation = Sedi)
      } else {
        dat <- dat |>
          dplyr::mutate(
            Sedimentation = dplyr::coalesce(
              Sedimentation,
              Sedi
            )
          ) |>
          dplyr::select(-Sedi)
      }
    }

    dat |>
      dplyr::mutate(
        source_file = basename(arq)
      )
  }
)

dados_brutos <- dplyr::bind_rows(lista_fenotipos)

cat("Dimensão após empilhamento dos arquivos:\n")
print(dim(dados_brutos))

cat("\nNúmero de colunas encontradas:", ncol(dados_brutos), "\n")
cat("Número de linhas encontradas:", nrow(dados_brutos), "\n")


# ==============================================================================
# 5. PADRONIZAÇÃO DE NOMES, TIPOS E CORREÇÕES RELEVANTES
# ==============================================================================

titulo_console("3. PADRONIZAÇÃO DOS DADOS")

# Garante que todas as colunas esperadas existam.
colunas_esperadas <- c(
  "BRISONr", "Treatment", "Block", "Row", "Column",
  "Year", "Location",
  "Sowing_date", "Emergence_date", "BBCH59", "BBCH87",
  "Plantheight", "Seedyield", "Seedyield_bio", "Biomass_bio",
  "Harvest_Index", "TKW_plot", "TKW_bio", "Spike_number",
  "Stripe_rust", "Powdery_mildew", "Leaf_rust",
  "Septoria", "DTR", "Fusarium",
  "Falling_number", "Crude_protein", "Sedimentation"
)

faltantes_colunas <- setdiff(colunas_esperadas, names(dados_brutos))

if (length(faltantes_colunas) > 0) {
  for (nm in faltantes_colunas) {
    dados_brutos[[nm]] <- NA_character_
  }
}


# Correção usada no repositório para KAL 2018-2020:
# BRISONr_214:220 -> BRISONr_222:228.
mapa_kal <- tibble::tibble(
  BRISONr_antigo = paste0("BRISONr_", 214:220),
  BRISONr_novo   = paste0("BRISONr_", 222:228)
)

dados <- dados_brutos |>
  dplyr::mutate(
    Year = numero_seguro(Year),
    Location_original = Location,
    Treatment_original = Treatment,

    BRISONr = stringr::str_replace(
      BRISONr,
      "BRISONR",
      "BRISONr"
    )
  ) |>
  dplyr::left_join(
    mapa_kal,
    by = c("BRISONr" = "BRISONr_antigo")
  ) |>
  dplyr::mutate(
    BRISONr = dplyr::case_when(
      Location_original == "KAD" &
        Year %in% 2018:2020 &
        !is.na(BRISONr_novo) ~ BRISONr_novo,

      TRUE ~ BRISONr
    )
  ) |>
  dplyr::select(-BRISONr_novo) |>
  dplyr::mutate(
    Treatment = padronizar_tratamento(
      Treatment_original,
      Location_original
    ),

    Location = padronizar_local(Location_original),

    Block = dplyr::case_when(
      Block == "1" ~ "B1",
      Block == "2" ~ "B2",
      TRUE ~ Block
    )
  )


# Conversão das variáveis quantitativas.
colunas_numericas <- c(
  "Row", "Column", "Year",
  "Sowing_date", "Emergence_date", "BBCH59", "BBCH87",
  "Plantheight", "Seedyield", "Seedyield_bio", "Biomass_bio",
  "Harvest_Index", "TKW_plot", "TKW_bio", "Spike_number",
  "Stripe_rust", "Powdery_mildew", "Leaf_rust",
  "Septoria", "DTR", "Fusarium",
  "Falling_number", "Crude_protein", "Sedimentation"
)

dados <- dados |>
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(colunas_numericas),
      numero_seguro
    )
  )


# Filtros de faixas irreais usados no pré-processamento da publicação.
dados <- dados |>
  dplyr::mutate(
    Plantheight = ifelse(
      Plantheight > 150,
      NA_real_,
      Plantheight
    ),

    Seedyield = ifelse(
      Seedyield < 0 | Seedyield > 3000,
      NA_real_,
      Seedyield
    ),

    Seedyield_bio = ifelse(
      Seedyield_bio > 2000,
      NA_real_,
      Seedyield_bio
    ),

    Biomass_bio = ifelse(
      Biomass_bio > 3500,
      NA_real_,
      Biomass_bio
    ),

    TKW_plot = ifelse(
      TKW_plot > 80,
      NA_real_,
      TKW_plot
    ),

    Spike_number = ifelse(
      Spike_number > 1500,
      NA_real_,
      Spike_number
    ),

    Leaf_rust = ifelse(
      Leaf_rust < 0,
      NA_real_,
      Leaf_rust
    ),

    Crude_protein = ifelse(
      Crude_protein < 5,
      NA_real_,
      Crude_protein
    )
  )


# Completa o índice de colheita quando possível.
dados <- dados |>
  dplyr::mutate(
    Harvest_Index = dplyr::if_else(
      is.na(Harvest_Index) &
        !is.na(Seedyield_bio) &
        !is.na(Biomass_bio) &
        Biomass_bio != 0,

      Seedyield_bio / Biomass_bio,
      Harvest_Index
    ),

    Harvest_Index = ifelse(
      Harvest_Index < 0.1 | Harvest_Index > 0.8,
      NA_real_,
      Harvest_Index
    )
  )


# Tratamento das doenças:
# se uma doença foi avaliada dentro da condição, NA e valores negativos são 0;
# se toda a condição for NA, permanece NA.
doencas <- c(
  "Stripe_rust",
  "Powdery_mildew",
  "Leaf_rust",
  "Septoria",
  "DTR",
  "Fusarium"
)

dados <- dados |>
  dplyr::group_by(
    Treatment,
    Location,
    Year
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(doencas),
      ~ {
        if (all(is.na(.x))) {
          .x
        } else {
          ifelse(is.na(.x) | .x < 0, 0, .x)
        }
      }
    )
  ) |>
  dplyr::ungroup()


# Variáveis derivadas do script público.
dados <- dados |>
  dplyr::mutate(
    Grain_per_spike = dplyr::case_when(
      !is.na(TKW_bio) &
        !is.na(Seedyield_bio) &
        !is.na(Spike_number) &
        TKW_bio > 0 &
        Spike_number > 0 ~
        1000 * Seedyield_bio / (TKW_bio * Spike_number),

      is.na(TKW_bio) &
        !is.na(TKW_plot) &
        !is.na(Seedyield_bio) &
        !is.na(Spike_number) &
        TKW_plot > 0 &
        Spike_number > 0 ~
        1000 * Seedyield_bio / (TKW_plot * Spike_number),

      TRUE ~ NA_real_
    ),

    TGW = dplyr::coalesce(
      TKW_plot,
      TKW_bio
    ),

    Grain = dplyr::if_else(
      !is.na(Seedyield) &
        !is.na(TGW) &
        TGW > 0,

      Seedyield * 1000 / TGW,
      NA_real_
    ),

    Protein_yield = Seedyield * Crude_protein / 100
  )


# Outliers de Harvest_Index e Spike_number:
# média +/- 4 desvios-padrão por condição e cultivar.
dados <- dados |>
  dplyr::group_by(
    Year,
    Location,
    Treatment,
    BRISONr
  ) |>
  dplyr::mutate(
    HI_media = mean(Harvest_Index, na.rm = TRUE),
    HI_dp = sd(Harvest_Index, na.rm = TRUE),

    SN_media = mean(Spike_number, na.rm = TRUE),
    SN_dp = sd(Spike_number, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    HI_media = ifelse(is.nan(HI_media), NA_real_, HI_media),
    SN_media = ifelse(is.nan(SN_media), NA_real_, SN_media),

    Harvest_Index = ifelse(
      !is.na(HI_dp) &
        !is.na(Harvest_Index) &
        (
          Harvest_Index < HI_media - 4 * HI_dp |
            Harvest_Index > HI_media + 4 * HI_dp
        ),

      NA_real_,
      Harvest_Index
    ),

    Spike_number = ifelse(
      !is.na(SN_dp) &
        !is.na(Spike_number) &
        (
          Spike_number < SN_media - 4 * SN_dp |
            Spike_number > SN_media + 4 * SN_dp
        ),

      NA_real_,
      Spike_number
    )
  ) |>
  dplyr::select(
    -HI_media,
    -HI_dp,
    -SN_media,
    -SN_dp
  )


# Derivação após filtro do HI.
dados <- dados |>
  dplyr::mutate(
    Biomass = dplyr::if_else(
      !is.na(Seedyield) &
        !is.na(Harvest_Index) &
        Harvest_Index > 0,

      Seedyield / Harvest_Index,
      NA_real_
    ),

    Straw = dplyr::if_else(
      !is.na(Biomass) &
        !is.na(Harvest_Index),

      Biomass * (1 - Harvest_Index),
      NA_real_
    )
  ) |>
  dplyr::filter(
    Treatment != "LLN_WF_RF",
    !is.na(BRISONr),
    BRISONr != "BRISONr_NA"
  ) |>
  dplyr::mutate(
    BRISONr = ifelse(
      BRISONr == "BRISONr_?",
      "BRISONr_229",
      BRISONr
    ),

    Condition = paste(
      Year,
      Location,
      Treatment,
      sep = "__"
    ),

    Nitrogen = stringr::word(
      Treatment,
      1,
      sep = "_"
    ),

    Fungicide = stringr::word(
      Treatment,
      2,
      sep = "_"
    ),

    Water = stringr::word(
      Treatment,
      3,
      sep = "_"
    )
  )


cat("Dimensão do banco harmonizado:\n")
print(dim(dados))

cat("\nNúmero total de condições Ano x Local x Manejo:\n")
print(dplyr::n_distinct(dados$Condition))

cat("\nManejos encontrados após padronização:\n")
print(sort(unique(dados$Treatment)))

cat("\nNúmero de cultivares únicos:\n")
print(dplyr::n_distinct(dados$BRISONr))


# ==============================================================================
# 6. AUDITORIA DAS 100 CONDIÇÕES
# ==============================================================================

titulo_console("4. AUDITORIA DAS CONDIÇÕES DE CULTIVO")

auditoria_condicoes <- dados |>
  dplyr::group_by(
    Year,
    Location,
    Treatment,
    Condition,
    Nitrogen,
    Fungicide,
    Water
  ) |>
  dplyr::summarise(
    n_parcelas = dplyr::n(),
    n_cultivares = dplyr::n_distinct(BRISONr),
    n_yield = sum(!is.na(Seedyield)),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    Year,
    Location,
    Treatment
  )

print(auditoria_condicoes, n = nrow(auditoria_condicoes))

readr::write_csv(
  auditoria_condicoes,
  file.path(
    PASTA_SAIDAS,
    "01_auditoria_condicoes.csv"
  )
)


# ==============================================================================
# 7. SUBCONJUNTO BALANCEADO DA ANÁLISE PRINCIPAL
# ==============================================================================

titulo_console("5. CONSTRUÇÃO DO SUBCONJUNTO BALANCEADO")

anos_principais <- 2015:2017

locais_principais <- c(
  "GGE",
  "HAN",
  "KAL",
  "KIE",
  "QLB"
)

manejos_principais <- c(
  "HN_WF_RF",
  "HN_NF_RF",
  "LN_NF_RF"
)

dados_principal <- dados |>
  dplyr::filter(
    Year %in% anos_principais,
    Location %in% locais_principais,
    Treatment %in% manejos_principais
  )

auditoria_principal <- dados_principal |>
  dplyr::group_by(
    Year,
    Location,
    Treatment,
    Condition
  ) |>
  dplyr::summarise(
    n_parcelas = dplyr::n(),
    n_cultivares = dplyr::n_distinct(BRISONr),
    n_yield = sum(!is.na(Seedyield)),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    Year,
    Location,
    Treatment
  )

cat("Número de condições no subconjunto principal:\n")
print(nrow(auditoria_principal))

cat("\nNúmero mínimo e máximo de cultivares por condição:\n")
print(
  range(
    auditoria_principal$n_cultivares,
    na.rm = TRUE
  )
)

cat("\nAuditoria das 45 condições:\n")
print(auditoria_principal, n = nrow(auditoria_principal))

if (nrow(auditoria_principal) != 45) {
  warning(
    "O subconjunto principal não possui 45 condições. ",
    "Verifique a importação/padronização dos arquivos."
  )
}

readr::write_csv(
  auditoria_principal,
  file.path(
    PASTA_SAIDAS,
    "02_auditoria_subconjunto_45_condicoes.csv"
  )
)


# ==============================================================================
# 8. FUNÇÃO PARA ESTIMAR BLUEs
# ==============================================================================

ajustar_blue <- function(dat, trait) {

  dd <- dat |>
    dplyr::filter(
      !is.na(.data[[trait]]),
      !is.na(BRISONr),
      !is.na(Row),
      !is.na(Column)
    ) |>
    dplyr::mutate(
      BRISONr = factor(BRISONr),
      Row_f = factor(Row),
      Column_f = factor(Column)
    )

  n_obs <- nrow(dd)
  n_gen <- dplyr::n_distinct(dd$BRISONr)

  if (
    n_obs < 10 ||
      n_gen < MIN_N_PAR ||
      dplyr::n_distinct(dd$Row_f) < 2 ||
      dplyr::n_distinct(dd$Column_f) < 2
  ) {
    return(
      tibble::tibble(
        BRISONr = character(),
        BLUE = numeric(),
        Trait = trait,
        n_obs_modelo = n_obs,
        n_genotipos_modelo = n_gen,
        singular = NA,
        status_modelo = "dados_insuficientes"
      )
    )
  }

  formula_modelo <- stats::as.formula(
    paste0(
      trait,
      " ~ 0 + BRISONr + (1 | Row_f) + (1 | Column_f)"
    )
  )

  ajuste <- tryCatch(
    lme4::lmer(
      formula_modelo,
      data = dd,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(
          maxfun = 200000
        )
      )
    ),
    error = function(e) e
  )

  if (inherits(ajuste, "error")) {

    # Fallback explícito: médias por cultivar.
    # O status é preservado para que o usuário saiba que não é BLUE espacial.
    return(
      dd |>
        dplyr::group_by(BRISONr) |>
        dplyr::summarise(
          BLUE = mean(
            .data[[trait]],
            na.rm = TRUE
          ),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          BRISONr = as.character(BRISONr),
          Trait = trait,
          n_obs_modelo = n_obs,
          n_genotipos_modelo = n_gen,
          singular = NA,
          status_modelo = paste0(
            "fallback_media: ",
            ajuste$message
          )
        )
    )
  }

  coefs <- lme4::fixef(ajuste)

  tibble::tibble(
    termo = names(coefs),
    BLUE = as.numeric(coefs)
  ) |>
    dplyr::mutate(
      BRISONr = stringr::str_remove(
        termo,
        "^BRISONr"
      ),

      Trait = trait,
      n_obs_modelo = n_obs,
      n_genotipos_modelo = n_gen,
      singular = lme4::isSingular(
        ajuste,
        tol = 1e-5
      ),

      status_modelo = "lmer_ok"
    ) |>
    dplyr::select(
      BRISONr,
      BLUE,
      Trait,
      n_obs_modelo,
      n_genotipos_modelo,
      singular,
      status_modelo
    )
}


# ==============================================================================
# 9. ESTIMAÇÃO DOS BLUEs NAS 45 CONDIÇÕES
# ==============================================================================

titulo_console("6. ESTIMAÇÃO DOS BLUEs")

metadados_condicoes <- dados_principal |>
  dplyr::distinct(
    Condition,
    Year,
    Location,
    Treatment,
    Nitrogen,
    Fungicide,
    Water
  ) |>
  dplyr::arrange(
    Year,
    Location,
    Treatment
  )

resultados_blues <- vector(
  "list",
  nrow(metadados_condicoes) *
    length(TRAITS_ANALISE)
)

contador <- 1L

for (i in seq_len(nrow(metadados_condicoes))) {

  cond_i <- metadados_condicoes$Condition[i]

  cat(
    sprintf(
      "[%02d/%02d] Ajustando condição: %s\n",
      i,
      nrow(metadados_condicoes),
      cond_i
    )
  )

  dat_i <- dados_principal |>
    dplyr::filter(
      Condition == cond_i
    )

  for (trait_i in TRAITS_ANALISE) {

    res_i <- ajustar_blue(
      dat = dat_i,
      trait = trait_i
    ) |>
      dplyr::mutate(
        Condition = cond_i
      )

    resultados_blues[[contador]] <- res_i
    contador <- contador + 1L
  }
}

blues <- dplyr::bind_rows(resultados_blues) |>
  dplyr::left_join(
    metadados_condicoes,
    by = "Condition"
  ) |>
  dplyr::mutate(
    Trait_abbrev = unname(
      ABREV_TRAITS[Trait]
    )
  )

cat("\nResumo dos ajustes dos BLUEs:\n")
print(
  blues |>
    dplyr::distinct(
      Condition,
      Trait,
      status_modelo,
      singular
    ) |>
    dplyr::count(
      status_modelo,
      singular,
      name = "n_modelos"
    )
)

cat("\nNúmero de BLUEs estimados por característica:\n")
print(
  blues |>
    dplyr::count(
      Trait_abbrev,
      name = "n_BLUEs"
    ) |>
    dplyr::arrange(
      Trait_abbrev
    )
)

readr::write_csv(
  blues,
  file.path(
    PASTA_SAIDAS,
    "03_BLUEs_45_condicoes_9_traits.csv"
  )
)


# ==============================================================================
# 10. FUNÇÕES DE INFORMAÇÃO MÚTUA
# ==============================================================================

# Discretização por postos:
#   - produz classes aproximadamente equiprováveis;
#   - mantém empates com o mesmo posto médio;
#   - evita depender de limites absolutos de escala entre ambientes.
discretizar_por_postos <- function(x, nbins = 6) {

  if (all(is.na(x))) {
    return(
      rep(
        NA_integer_,
        length(x)
      )
    )
  }

  r <- rank(
    x,
    ties.method = "average",
    na.last = "keep"
  )

  n_validos <- sum(!is.na(r))

  out <- ceiling(
    nbins * r / n_validos
  )

  out <- pmin(
    nbins,
    pmax(
      1,
      out
    )
  )

  as.integer(out)
}


entropia_discreta <- function(x) {

  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  p <- prop.table(table(x))

  -sum(
    p * log(p)
  )
}


mi_discreta <- function(x, y) {

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


nmi_discreta <- function(x, y) {

  mi <- mi_discreta(
    x,
    y
  )

  hx <- entropia_discreta(x)
  hy <- entropia_discreta(y)

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


# avaliar_associacao <- function(
#   x,
#   y,
#   nbins = N_BINS_MI,
#   n_perm = 0
# ) {
# 
#   ok <- complete.cases(
#     x,
#     y
#   )
# 
#   x <- x[ok]
#   y <- y[ok]
# 
#   n <- length(x)
# 
#   if (
#     n < MIN_N_PAR ||
#       sd(x) == 0 ||
#       sd(y) == 0
#   ) {
#     return(
#       tibble::tibble(
#         n = n,
#         R2_SMA = NA_real_,
#         Spearman = NA_real_,
#         NMI_obs = NA_real_,
#         NMI_perm_media = NA_real_,
#         NMI_corrigida = NA_real_,
#         Z_MI = NA_real_,
#         p_MI = NA_real_
#       )
#     )
#   }
# 
#   r_pearson <- suppressWarnings(
#     cor(
#       x,
#       y,
#       method = "pearson"
#     )
#   )
# 
#   rho <- suppressWarnings(
#     cor(
#       x,
#       y,
#       method = "spearman"
#     )
#   )
# 
#   xd <- discretizar_por_postos(
#     x,
#     nbins = nbins
#   )
# 
#   yd <- discretizar_por_postos(
#     y,
#     nbins = nbins
#   )
# 
#   nmi_obs <- nmi_discreta(
#     xd,
#     yd
#   )
# 
#   if (n_perm <= 0) {
# 
#     return(
#       tibble::tibble(
#         n = n,
#         R2_SMA = r_pearson^2,
#         Spearman = rho,
#         NMI_obs = nmi_obs,
#         NMI_perm_media = NA_real_,
#         NMI_corrigida = NA_real_,
#         Z_MI = NA_real_,
#         p_MI = NA_real_
#       )
#     )
#   }
# 
#   nmi_perm <- replicate(
#     n_perm,
#     nmi_discreta(
#       xd,
#       sample(
#         yd,
#         replace = FALSE
#       )
#     )
#   )
# 
#   media_perm <- mean(
#     nmi_perm,
#     na.rm = TRUE
#   )
# 
#   dp_perm <- sd(
#     nmi_perm,
#     na.rm = TRUE
#   )
# 
#   p_perm <- (
#     1 +
#       sum(
#         nmi_perm >= nmi_obs,
#         na.rm = TRUE
#       )
#   ) / (
#     n_perm + 1
#   )
# 
#   tibble::tibble(
#     n = n,
#     R2_SMA = r_pearson^2,
#     Spearman = rho,
#     NMI_obs = nmi_obs,
#     NMI_perm_media = media_perm,
#     NMI_corrigida = max(
#       0,
#       nmi_obs - media_perm
#     ),
#     Z_MI = ifelse(
#       dp_perm > 0,
#       (nmi_obs - media_perm) / dp_perm,
#       NA_real_
#     ),
#     p_MI = p_perm
#   )
# }

avaliar_associacao <- function(
    x,
    y,
    nbins = N_BINS_MI,
    n_perm = 0
) {
  
  # Caso algum ambiente não tenha a característica
  if (
    is.null(x) ||
    is.null(y)
  ) {
    
    return(
      tibble::tibble(
        n = 0,
        R2_SMA = NA_real_,
        Spearman = NA_real_,
        NMI_obs = NA_real_,
        NMI_perm_media = NA_real_,
        NMI_corrigida = NA_real_,
        Z_MI = NA_real_,
        p_MI = NA_real_
      )
    )
  }
  
  # Garante formato numérico
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Mantém apenas pares completos e finitos
  ok <- complete.cases(
    x,
    y
  ) &
    is.finite(x) &
    is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  # Verifica tamanho amostral e variabilidade
  if (
    n < MIN_N_PAR ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    
    return(
      tibble::tibble(
        n = n,
        R2_SMA = NA_real_,
        Spearman = NA_real_,
        NMI_obs = NA_real_,
        NMI_perm_media = NA_real_,
        NMI_corrigida = NA_real_,
        Z_MI = NA_real_,
        p_MI = NA_real_
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Associação linear
  # ---------------------------------------------------------------------------
  
  r_pearson <- suppressWarnings(
    cor(
      x,
      y,
      method = "pearson"
    )
  )
  
  rho <- suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )
  
  # ---------------------------------------------------------------------------
  # Discretização para MI
  # ---------------------------------------------------------------------------
  
  xd <- discretizar_por_postos(
    x,
    nbins = nbins
  )
  
  yd <- discretizar_por_postos(
    y,
    nbins = nbins
  )
  
  nmi_obs <- nmi_discreta(
    xd,
    yd
  )
  
  # ---------------------------------------------------------------------------
  # Sem permutação
  # ---------------------------------------------------------------------------
  
  if (n_perm <= 0) {
    
    return(
      tibble::tibble(
        n = n,
        R2_SMA = r_pearson^2,
        Spearman = rho,
        NMI_obs = nmi_obs,
        NMI_perm_media = NA_real_,
        NMI_corrigida = NA_real_,
        Z_MI = NA_real_,
        p_MI = NA_real_
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Permutação
  # ---------------------------------------------------------------------------
  
  nmi_perm <- replicate(
    n_perm,
    nmi_discreta(
      xd,
      sample(
        yd,
        replace = FALSE
      )
    )
  )
  
  media_perm <- mean(
    nmi_perm,
    na.rm = TRUE
  )
  
  dp_perm <- sd(
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
    R2_SMA = r_pearson^2,
    Spearman = rho,
    NMI_obs = nmi_obs,
    NMI_perm_media = media_perm,
    NMI_corrigida = max(
      0,
      nmi_obs - media_perm
    ),
    Z_MI = ifelse(
      is.finite(dp_perm) &&
        dp_perm > 0,
      (nmi_obs - media_perm) / dp_perm,
      NA_real_
    ),
    p_MI = p_perm
  )
}


# ==============================================================================
# 11. ANÁLISE PRINCIPAL:
#     CONSISTÊNCIA DA PRODUTIVIDADE ENTRE AS 45 CONDIÇÕES
# ==============================================================================

titulo_console(
  "7. ANÁLISE PRINCIPAL: PRODUTIVIDADE ENTRE CONDIÇÕES"
)

blues_yield <- blues |>
  dplyr::filter(
    Trait == "Seedyield"
  ) |>
  dplyr::select(
    BRISONr,
    Condition,
    BLUE
  )

yield_wide <- blues_yield |>
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = BLUE
  )

condicoes <- metadados_condicoes$Condition

pares_condicoes <- t(
  utils::combn(
    condicoes,
    2
  )
)

resultados_pares <- vector(
  "list",
  nrow(pares_condicoes)
)

set.seed(SEMENTE)

for (i in seq_len(nrow(pares_condicoes))) {

  if (
    i == 1 ||
      i %% 100 == 0 ||
      i == nrow(pares_condicoes)
  ) {
    cat(
      sprintf(
        "Par %d de %d\n",
        i,
        nrow(pares_condicoes)
      )
    )
  }

  e1 <- pares_condicoes[i, 1]
  e2 <- pares_condicoes[i, 2]

  x <- yield_wide[[e1]]
  y <- yield_wide[[e2]]

  met <- avaliar_associacao(
    x,
    y,
    nbins = N_BINS_MI,
    n_perm = N_PERM_MI
  )

  # Sensibilidade ao número de classes.
  sens <- purrr::map_dbl(
    BINS_SENSIBILIDADE,
    function(b) {

      ok <- complete.cases(
        x,
        y
      )

      xd <- discretizar_por_postos(
        x[ok],
        nbins = b
      )

      yd <- discretizar_por_postos(
        y[ok],
        nbins = b
      )

      nmi_discreta(
        xd,
        yd
      )
    }
  )

  resultados_pares[[i]] <- dplyr::bind_cols(
    tibble::tibble(
      Environment_1 = e1,
      Environment_2 = e2
    ),
    met,
    tibble::as_tibble_row(
      stats::setNames(
        as.list(sens),
        paste0(
          "NMI_bins_",
          BINS_SENSIBILIDADE
        )
      )
    )
  )
}

pares_yield <- dplyr::bind_rows(
  resultados_pares
) |>
  dplyr::mutate(
    p_MI_FDR = p.adjust(
      p_MI,
      method = "BH"
    )
  )


# Junta os metadados dos dois ambientes.
meta_1 <- metadados_condicoes |>
  dplyr::rename_with(
    ~ paste0(.x, "_1"),
    -Condition
  ) |>
  dplyr::rename(
    Environment_1 = Condition
  )

meta_2 <- metadados_condicoes |>
  dplyr::rename_with(
    ~ paste0(.x, "_2"),
    -Condition
  ) |>
  dplyr::rename(
    Environment_2 = Condition
  )

pares_yield <- pares_yield |>
  dplyr::left_join(
    meta_1,
    by = "Environment_1"
  ) |>
  dplyr::left_join(
    meta_2,
    by = "Environment_2"
  ) |>
  dplyr::mutate(
    same_year = Year_1 == Year_2,
    same_location = Location_1 == Location_2,
    same_management = Treatment_1 == Treatment_2,
    same_nitrogen = Nitrogen_1 == Nitrogen_2,
    same_fungicide = Fungicide_1 == Fungicide_2,
    same_water = Water_1 == Water_2
  )


# Índice exploratório:
# pares com NMI relativamente alta em comparação ao seu R².
pares_yield <- pares_yield |>
  dplyr::mutate(
    z_R2 = as.numeric(
      scale(R2_SMA)
    ),

    z_NMI = as.numeric(
      scale(NMI_corrigida)
    ),

    MI_advantage = z_NMI - z_R2
  ) |>
  dplyr::arrange(
    dplyr::desc(MI_advantage)
  )


cat("\nResumo das 990 comparações de produtividade:\n")
print(
  pares_yield |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      n_medio_cultivares = mean(n),
      R2_medio = mean(R2_SMA, na.rm = TRUE),
      R2_mediano = median(R2_SMA, na.rm = TRUE),
      NMI_media = mean(NMI_obs, na.rm = TRUE),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      NMI_corrigida_mediana = median(
        NMI_corrigida,
        na.rm = TRUE
      ),
      pares_MI_FDR_5pct = sum(
        p_MI_FDR < 0.05,
        na.rm = TRUE
      )
    )
)

cat(
  "\n10 pares com maior Informação Mútua corrigida:\n"
)

print(
  pares_yield |>
    dplyr::arrange(
      dplyr::desc(NMI_corrigida)
    ) |>
    dplyr::select(
      Environment_1,
      Environment_2,
      n,
      R2_SMA,
      Spearman,
      NMI_obs,
      NMI_corrigida,
      Z_MI,
      p_MI_FDR
    ) |>
    dplyr::slice_head(
      n = 10
    ),
  n = 10
)


cat(
  "\n10 pares com maior MI relativa ao R² linear:\n"
)

print(
  pares_yield |>
    dplyr::select(
      Environment_1,
      Environment_2,
      n,
      R2_SMA,
      Spearman,
      NMI_corrigida,
      MI_advantage,
      p_MI_FDR
    ) |>
    dplyr::slice_head(
      n = 10
    ),
  n = 10
)


cat(
  "\n10 pares com menor Informação Mútua corrigida:\n"
)

print(
  pares_yield |>
    dplyr::arrange(
      NMI_corrigida
    ) |>
    dplyr::select(
      Environment_1,
      Environment_2,
      n,
      R2_SMA,
      Spearman,
      NMI_corrigida,
      p_MI_FDR
    ) |>
    dplyr::slice_head(
      n = 10
    ),
  n = 10
)


readr::write_csv(
  pares_yield,
  file.path(
    PASTA_SAIDAS,
    "04_pares_produtividade_R2_MI.csv"
  )
)


# ==============================================================================
# 12. SENSIBILIDADE À DISCRETIZAÇÃO
# ==============================================================================

titulo_console("8. SENSIBILIDADE DA MI AO NÚMERO DE CLASSES")

sensibilidade_mi <- tibble::tibble(
  comparacao = c(
    "4 vs 6 classes",
    "4 vs 8 classes",
    "6 vs 8 classes"
  ),

  correlacao_spearman = c(
    cor(
      pares_yield$NMI_bins_4,
      pares_yield$NMI_bins_6,
      method = "spearman",
      use = "complete.obs"
    ),

    cor(
      pares_yield$NMI_bins_4,
      pares_yield$NMI_bins_8,
      method = "spearman",
      use = "complete.obs"
    ),

    cor(
      pares_yield$NMI_bins_6,
      pares_yield$NMI_bins_8,
      method = "spearman",
      use = "complete.obs"
    )
  )
)

print(sensibilidade_mi)

readr::write_csv(
  sensibilidade_mi,
  file.path(
    PASTA_SAIDAS,
    "05_sensibilidade_numero_classes_MI.csv"
  )
)


# ==============================================================================
# 13. RESUMOS POR ANO, LOCAL E MANEJO
# ==============================================================================

titulo_console(
  "9. RESUMOS DA CONSISTÊNCIA POR AGRUPAMENTOS"
)

resumo_mesmo_diferente <- dplyr::bind_rows(

  pares_yield |>
    dplyr::group_by(
      grupo = "Ano",
      mesmo = same_year
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  pares_yield |>
    dplyr::group_by(
      grupo = "Local",
      mesmo = same_location
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  pares_yield |>
    dplyr::group_by(
      grupo = "Manejo",
      mesmo = same_management
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  pares_yield |>
    dplyr::group_by(
      grupo = "Nitrogênio",
      mesmo = same_nitrogen
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  pares_yield |>
    dplyr::group_by(
      grupo = "Fungicida",
      mesmo = same_fungicide
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)

print(
  resumo_mesmo_diferente,
  n = nrow(resumo_mesmo_diferente)
)

readr::write_csv(
  resumo_mesmo_diferente,
  file.path(
    PASTA_SAIDAS,
    "06_resumo_mesmo_diferente.csv"
  )
)


# Resumos dentro do mesmo nível:
# por exemplo, comparações entre ambientes de Hannover.
resumo_por_nivel <- dplyr::bind_rows(

  pares_yield |>
    dplyr::filter(
      same_location
    ) |>
    dplyr::group_by(
      tipo = "Location",
      nivel = Location_1
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  pares_yield |>
    dplyr::filter(
      same_year
    ) |>
    dplyr::group_by(
      tipo = "Year",
      nivel = as.character(Year_1)
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  pares_yield |>
    dplyr::filter(
      same_management
    ) |>
    dplyr::group_by(
      tipo = "Treatment",
      nivel = Treatment_1
    ) |>
    dplyr::summarise(
      n_pares = dplyr::n(),
      R2_medio = mean(
        R2_SMA,
        na.rm = TRUE
      ),
      NMI_corrigida_media = mean(
        NMI_corrigida,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)

cat("\nResumo dentro de cada nível:\n")
print(
  resumo_por_nivel,
  n = nrow(resumo_por_nivel)
)

readr::write_csv(
  resumo_por_nivel,
  file.path(
    PASTA_SAIDAS,
    "07_resumo_por_nivel.csv"
  )
)


# ==============================================================================
# 14. TESTE POR PERMUTAÇÃO DE RÓTULOS DE AMBIENTES
#
# Evita tratar os 990 pares como observações completamente independentes.
# A estatística é:
#   média da métrica para pares do mesmo grupo
#   menos
#   média da métrica para pares de grupos diferentes.
# ==============================================================================

teste_rotulos_ambientes <- function(
  pares,
  metadados,
  coluna_grupo,
  coluna_metrica,
  B = 999,
  seed = SEMENTE
) {

  mapa <- metadados |>
    dplyr::select(
      Condition,
      grupo = dplyr::all_of(coluna_grupo)
    )

  if (
    dplyr::n_distinct(mapa$grupo) < 2
  ) {
    return(
      tibble::tibble(
        grupo = coluna_grupo,
        metrica = coluna_metrica,
        diferenca_observada = NA_real_,
        p_permutacao = NA_real_,
        B = B
      )
    )
  }

  pdat <- pares |>
    dplyr::select(
      Environment_1,
      Environment_2,
      valor = dplyr::all_of(coluna_metrica)
    ) |>
    dplyr::left_join(
      mapa |>
        dplyr::rename(
          Environment_1 = Condition,
          grupo_1 = grupo
        ),
      by = "Environment_1"
    ) |>
    dplyr::left_join(
      mapa |>
        dplyr::rename(
          Environment_2 = Condition,
          grupo_2 = grupo
        ),
      by = "Environment_2"
    )

  stat_fun <- function(g1, g2, valor) {

    mesmo <- g1 == g2

    if (
      all(mesmo) ||
        all(!mesmo)
    ) {
      return(NA_real_)
    }

    mean(
      valor[mesmo],
      na.rm = TRUE
    ) -
      mean(
        valor[!mesmo],
        na.rm = TRUE
      )
  }

  obs <- stat_fun(
    pdat$grupo_1,
    pdat$grupo_2,
    pdat$valor
  )

  set.seed(seed)

  labels_orig <- mapa$grupo
  nomes_env <- mapa$Condition

  estatisticas_perm <- replicate(
    B,
    {

      labels_perm <- sample(
        labels_orig,
        replace = FALSE
      )

      lookup <- stats::setNames(
        labels_perm,
        nomes_env
      )

      g1 <- unname(
        lookup[pdat$Environment_1]
      )

      g2 <- unname(
        lookup[pdat$Environment_2]
      )

      stat_fun(
        g1,
        g2,
        pdat$valor
      )
    }
  )

  p <- (
    1 +
      sum(
        abs(estatisticas_perm) >= abs(obs),
        na.rm = TRUE
      )
  ) / (
    B + 1
  )

  tibble::tibble(
    grupo = coluna_grupo,
    metrica = coluna_metrica,
    diferenca_observada = obs,
    p_permutacao = p,
    B = B
  )
}


titulo_console(
  "10. TESTES DE PERMUTAÇÃO DOS RÓTULOS DOS AMBIENTES"
)

grupos_testar <- c(
  "Year",
  "Location",
  "Treatment",
  "Nitrogen",
  "Fungicide"
)

metricas_testar <- c(
  "R2_SMA",
  "NMI_corrigida"
)

testes_rotulos <- tidyr::crossing(
  grupo = grupos_testar,
  metrica = metricas_testar
) |>
  purrr::pmap_dfr(
    function(grupo, metrica) {

      teste_rotulos_ambientes(
        pares = pares_yield,
        metadados = metadados_condicoes,
        coluna_grupo = grupo,
        coluna_metrica = metrica,
        B = 999,
        seed = SEMENTE
      )
    }
  )

print(
  testes_rotulos,
  n = nrow(testes_rotulos)
)

readr::write_csv(
  testes_rotulos,
  file.path(
    PASTA_SAIDAS,
    "08_testes_permutacao_rotulos_ambientes.csv"
  )
)


# ==============================================================================
# 15. GRÁFICOS DA ANÁLISE PRINCIPAL
# ==============================================================================

titulo_console("11. GRÁFICOS DA ANÁLISE PRINCIPAL")


# 15.1 NMI corrigida x R²
p_mi_r2 <- ggplot(
  pares_yield,
  aes(
    x = R2_SMA,
    y = NMI_corrigida
  )
) +
  geom_point(
    alpha = 0.55,
    size = 2
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    title = "Informação Mútua e consistência linear entre condições de cultivo",
    subtitle = "Cada ponto representa um par de condições Ano × Local × Manejo",
    x = expression(R[SMA]^2),
    y = "Informação Mútua Normalizada corrigida por permutação"
  ) +
  theme_minimal(
    base_size = 12
  )

print(p_mi_r2)

ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig01_NMI_corrigida_vs_R2.png"
  ),
  plot = p_mi_r2,
  width = 9,
  height = 7,
  dpi = 300
)


# Função para transformar pares em matriz longa simétrica.
matriz_pares_longa <- function(
  pares,
  coluna_valor,
  diagonal = 1
) {

  a <- pares |>
    dplyr::transmute(
      Environment_1,
      Environment_2,
      valor = .data[[coluna_valor]]
    )

  b <- a |>
    dplyr::transmute(
      Environment_1 = Environment_2,
      Environment_2 = Environment_1,
      valor = valor
    )

  diag_df <- tibble::tibble(
    Environment_1 = condicoes,
    Environment_2 = condicoes,
    valor = diagonal
  )

  dplyr::bind_rows(
    a,
    b,
    diag_df
  ) |>
    dplyr::mutate(
      Environment_1 = factor(
        Environment_1,
        levels = condicoes
      ),
      Environment_2 = factor(
        Environment_2,
        levels = rev(condicoes)
      )
    )
}


# 15.2 Heatmap NMI
heat_nmi <- matriz_pares_longa(
  pares_yield,
  "NMI_corrigida",
  diagonal = 1
)

p_heat_nmi <- ggplot(
  heat_nmi,
  aes(
    x = Environment_1,
    y = Environment_2,
    fill = valor
  )
) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "NMI\ncorrigida",
    option = "C"
  ) +
  labs(
    title = "Informação compartilhada entre condições de cultivo",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid = element_blank()
  )

print(p_heat_nmi)

ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig02_heatmap_NMI_corrigida.png"
  ),
  plot = p_heat_nmi,
  width = 16,
  height = 14,
  dpi = 300
)


# 15.3 Heatmap R²
heat_r2 <- matriz_pares_longa(
  pares_yield,
  "R2_SMA",
  diagonal = 1
)

p_heat_r2 <- ggplot(
  heat_r2,
  aes(
    x = Environment_1,
    y = Environment_2,
    fill = valor
  )
) +
  geom_tile() +
  scale_fill_viridis_c(
    name = expression(R[SMA]^2),
    option = "C",
    limits = c(0, 1)
  ) +
  labs(
    title = "Consistência linear entre condições de cultivo",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid = element_blank()
  )

print(p_heat_r2)

ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig03_heatmap_R2_SMA.png"
  ),
  plot = p_heat_r2,
  width = 16,
  height = 14,
  dpi = 300
)


# 15.4 Pares com maior vantagem informacional relativa ao R².
top_discrepantes <- pares_yield |>
  dplyr::slice_head(
    n = 6
  )

dados_discrepantes <- purrr::map2_dfr(
  top_discrepantes$Environment_1,
  top_discrepantes$Environment_2,
  function(e1, e2) {

    blues_yield |>
      dplyr::filter(
        Condition %in% c(
          e1,
          e2
        )
      ) |>
      dplyr::select(
        BRISONr,
        Condition,
        BLUE
      ) |>
      tidyr::pivot_wider(
        names_from = Condition,
        values_from = BLUE
      ) |>
      dplyr::transmute(
        BRISONr,
        x = .data[[e1]],
        y = .data[[e2]],
        par = paste(
          e1,
          "×",
          e2
        )
      ) |>
      dplyr::filter(
        complete.cases(
          x,
          y
        )
      )
  }
)

p_discrepantes <- ggplot(
  dados_discrepantes,
  aes(
    x = x,
    y = y
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.6
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linetype = 2
  ) +
  facet_wrap(
    ~ par,
    scales = "free",
    ncol = 2
  ) +
  labs(
    title = "Pares com maior informação relativa à associação linear",
    subtitle = "Linha contínua: ajuste linear; linha tracejada: tendência LOESS",
    x = "BLUE de produtividade no ambiente 1",
    y = "BLUE de produtividade no ambiente 2"
  ) +
  theme_minimal(
    base_size = 10
  )

print(p_discrepantes)

ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig04_pares_MI_alta_relativa_ao_R2.png"
  ),
  plot = p_discrepantes,
  width = 13,
  height = 12,
  dpi = 300
)


# ==============================================================================
# 16. CONSISTÊNCIA DAS NOVE CARACTERÍSTICAS ENTRE AMBIENTES
# ==============================================================================

titulo_console(
  "12. CONSISTÊNCIA DE NOVE CARACTERÍSTICAS ENTRE AMBIENTES"
)

resultados_traits_ambientes <- list()
contador_trait_par <- 1L

for (trait_i in TRAITS_ANALISE) {

  cat(
    "\nCaracterística:",
    ABREV_TRAITS[[trait_i]],
    "\n"
  )

  wide_i <- blues |>
    dplyr::filter(
      Trait == trait_i
    ) |>
    dplyr::select(
      BRISONr,
      Condition,
      BLUE
    ) |>
    tidyr::pivot_wider(
      names_from = Condition,
      values_from = BLUE
    )

  pares_i <- t(
    utils::combn(
      condicoes,
      2
    )
  )

  lista_i <- vector(
    "list",
    nrow(pares_i)
  )

  for (j in seq_len(nrow(pares_i))) {

    e1 <- pares_i[j, 1]
    e2 <- pares_i[j, 2]

    # x <- wide_i[[e1]]
    # y <- wide_i[[e2]]
    x <- if (e1 %in% names(wide_i)) {
      wide_i[[e1]]
    } else {
      NULL
    }
    
    y <- if (e2 %in% names(wide_i)) {
      wide_i[[e2]]
    } else {
      NULL
    }

    met <- avaliar_associacao(
      x,
      y,
      nbins = N_BINS_MI,
      n_perm = 0
    )

    lista_i[[j]] <- dplyr::bind_cols(
      tibble::tibble(
        Trait = trait_i,
        Trait_abbrev = ABREV_TRAITS[[trait_i]],
        Environment_1 = e1,
        Environment_2 = e2
      ),
      met
    )
  }

  resultados_traits_ambientes[[contador_trait_par]] <- dplyr::bind_rows(
    lista_i
  )

  contador_trait_par <- contador_trait_par + 1L
}

traits_consistencia <- dplyr::bind_rows(
  resultados_traits_ambientes
)

resumo_traits_consistencia <- traits_consistencia |>
  dplyr::group_by(
    Trait,
    Trait_abbrev
  ) |>
  dplyr::summarise(
    n_pares_validos = sum(
      !is.na(R2_SMA)
    ),
    R2_medio = mean(
      R2_SMA,
      na.rm = TRUE
    ),
    R2_mediano = median(
      R2_SMA,
      na.rm = TRUE
    ),
    NMI_media = mean(
      NMI_obs,
      na.rm = TRUE
    ),
    NMI_mediana = median(
      NMI_obs,
      na.rm = TRUE
    ),
    correlacao_R2_NMI = cor(
      R2_SMA,
      NMI_obs,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(NMI_media)
  )

print(
  resumo_traits_consistencia,
  n = nrow(resumo_traits_consistencia)
)

readr::write_csv(
  traits_consistencia,
  file.path(
    PASTA_SAIDAS,
    "09_consistencia_9_traits_por_pares.csv"
  )
)

readr::write_csv(
  resumo_traits_consistencia,
  file.path(
    PASTA_SAIDAS,
    "10_resumo_consistencia_9_traits.csv"
  )
)


p_box_traits_mi <- ggplot(
  traits_consistencia,
  aes(
    x = reorder(
      Trait_abbrev,
      NMI_obs,
      FUN = median,
      na.rm = TRUE
    ),
    y = NMI_obs
  )
) +
  geom_boxplot(
    outlier.alpha = 0.25
  ) +
  labs(
    title = "Consistência informacional das características entre ambientes",
    x = "Característica",
    y = "Informação Mútua Normalizada"
  ) +
  theme_minimal(
    base_size = 12
  )

print(p_box_traits_mi)

ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig05_boxplot_NMI_9_traits.png"
  ),
  plot = p_box_traits_mi,
  width = 9,
  height = 7,
  dpi = 300
)


p_box_traits_r2 <- ggplot(
  traits_consistencia,
  aes(
    x = reorder(
      Trait_abbrev,
      R2_SMA,
      FUN = median,
      na.rm = TRUE
    ),
    y = R2_SMA
  )
) +
  geom_boxplot(
    outlier.alpha = 0.25
  ) +
  labs(
    title = "Consistência linear das características entre ambientes",
    x = "Característica",
    y = expression(R[SMA]^2)
  ) +
  theme_minimal(
    base_size = 12
  )

print(p_box_traits_r2)

ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig06_boxplot_R2_9_traits.png"
  ),
  plot = p_box_traits_r2,
  width = 9,
  height = 7,
  dpi = 300
)


# ==============================================================================
# 17. ASSOCIAÇÃO DAS CARACTERÍSTICAS COM PRODUTIVIDADE DENTRO DE CADA AMBIENTE
# ==============================================================================

titulo_console(
  "13. ASSOCIAÇÃO DAS CARACTERÍSTICAS COM PRODUTIVIDADE"
)

traits_preditoras <- setdiff(
  TRAITS_ANALISE,
  "Seedyield"
)

resultados_trait_yield <- list()
contador_ty <- 1L

set.seed(
  SEMENTE + 1000
)

for (i in seq_len(nrow(metadados_condicoes))) {

  cond_i <- metadados_condicoes$Condition[i]

  gy_i <- blues |>
    dplyr::filter(
      Condition == cond_i,
      Trait == "Seedyield"
    ) |>
    dplyr::select(
      BRISONr,
      GY = BLUE
    )

  for (trait_i in traits_preditoras) {

    x_i <- blues |>
      dplyr::filter(
        Condition == cond_i,
        Trait == trait_i
      ) |>
      dplyr::select(
        BRISONr,
        X = BLUE
      )

    xy <- dplyr::inner_join(
      gy_i,
      x_i,
      by = "BRISONr"
    )

    met <- avaliar_associacao(
      x = xy$X,
      y = xy$GY,
      nbins = N_BINS_MI,
      n_perm = N_PERM_TRAIT_YIELD
    )

    resultados_trait_yield[[contador_ty]] <- dplyr::bind_cols(
      tibble::tibble(
        Condition = cond_i,
        Trait = trait_i,
        Trait_abbrev = ABREV_TRAITS[[trait_i]]
      ),
      met
    )

    contador_ty <- contador_ty + 1L
  }
}

trait_yield <- dplyr::bind_rows(
  resultados_trait_yield
) |>
  dplyr::left_join(
    metadados_condicoes,
    by = "Condition"
  ) |>
  dplyr::group_by(
    Trait
  ) |>
  dplyr::mutate(
    p_MI_FDR_dentro_trait = p.adjust(
      p_MI,
      method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    z_R2 = as.numeric(
      ave(
        R2_SMA,
        Trait,
        FUN = function(z) {
          as.numeric(
            scale(z)
          )
        }
      )
    ),

    z_NMI = as.numeric(
      ave(
        NMI_corrigida,
        Trait,
        FUN = function(z) {
          as.numeric(
            scale(z)
          )
        }
      )
    ),

    MI_advantage = z_NMI - z_R2
  )


resumo_trait_yield <- trait_yield |>
  dplyr::group_by(
    Trait,
    Trait_abbrev
  ) |>
  dplyr::summarise(
    n_ambientes = dplyr::n(),
    R2_medio = mean(
      R2_SMA,
      na.rm = TRUE
    ),
    R2_mediano = median(
      R2_SMA,
      na.rm = TRUE
    ),
    NMI_corrigida_media = mean(
      NMI_corrigida,
      na.rm = TRUE
    ),
    NMI_corrigida_mediana = median(
      NMI_corrigida,
      na.rm = TRUE
    ),
    ambientes_MI_FDR_5pct = sum(
      p_MI_FDR_dentro_trait < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      NMI_corrigida_media
    )
  )

cat(
  "\nResumo da associação de cada característica com a produtividade:\n"
)

print(
  resumo_trait_yield,
  n = nrow(resumo_trait_yield)
)


cat(
  "\n15 relações característica × ambiente com maior MI relativa ao R²:\n"
)

print(
  trait_yield |>
    dplyr::arrange(
      dplyr::desc(MI_advantage)
    ) |>
    dplyr::select(
      Condition,
      Trait_abbrev,
      n,
      R2_SMA,
      Spearman,
      NMI_corrigida,
      MI_advantage,
      p_MI_FDR_dentro_trait
    ) |>
    dplyr::slice_head(
      n = 15
    ),
  n = 15
)


readr::write_csv(
  trait_yield,
  file.path(
    PASTA_SAIDAS,
    "11_trait_vs_yield_por_ambiente.csv"
  )
)

readr::write_csv(
  resumo_trait_yield,
  file.path(
    PASTA_SAIDAS,
    "12_resumo_trait_vs_yield.csv"
  )
)


p_trait_yield <- ggplot(
  trait_yield,
  aes(
    x = Trait_abbrev,
    y = NMI_corrigida
  )
) +
  geom_boxplot(
    outlier.alpha = 0.25
  ) +
  labs(
    title = "Informação das características sobre a produtividade",
    subtitle = "Associação calculada entre BLUEs dos cultivares dentro de cada ambiente",
    x = "Característica",
    y = "NMI corrigida por permutação"
  ) +
  theme_minimal(
    base_size = 12
  )

print(p_trait_yield)

ggsave(
  filename = file.path(
    PASTA_SAIDAS,
    "fig07_trait_vs_yield_NMI.png"
  ),
  plot = p_trait_yield,
  width = 9,
  height = 7,
  dpi = 300
)


# ==============================================================================
# 18. AUDITORIA DOS DADOS DE CONTEXTO:
#     CLIMA, SOLO, FERTILIZAÇÃO E PROTEÇÃO DE PLANTAS
#
# Estes dados são importados e resumidos, mas não são usados de forma ingênua
# como milhares de observações independentes na análise de MI.
# ==============================================================================

titulo_console(
  "14. AUDITORIA DE CLIMA, SOLO E MANEJO"
)


# ------------------------------
# 18.1 Clima
# ------------------------------

if (length(arquivos_clima) > 0) {

  clima <- purrr::map_dfr(
    arquivos_clima,
    ler_csv_flexivel
  ) |>
    dplyr::mutate(
      Date = as.Date(Date),
      Location = padronizar_local(Location),

      dplyr::across(
        -c(
          Date,
          Location
        ),
        numero_seguro
      )
    )

  cat("\nDimensão da base meteorológica:\n")
  print(dim(clima))

  cat("\nCobertura temporal por local:\n")
  print(
    clima |>
      dplyr::group_by(
        Location
      ) |>
      dplyr::summarise(
        data_inicial = min(
          Date,
          na.rm = TRUE
        ),
        data_final = max(
          Date,
          na.rm = TRUE
        ),
        n_dias = dplyr::n(),
        .groups = "drop"
      )
  )

  readr::write_csv(
    clima,
    file.path(
      PASTA_SAIDAS,
      "13_clima_harmonizado.csv"
    )
  )
}


# ------------------------------
# 18.2 Solo
# ------------------------------

arquivo_solo <- file.path(
  PASTA_DADOS,
  "soil.xlsx"
)

if (file.exists(arquivo_solo)) {

  solo <- readxl::read_excel(
    arquivo_solo
  ) |>
    dplyr::mutate(
      Location = padronizar_local(
        Location
      ),
      Year = as.integer(Year)
    )

  cat("\nResumo da base de solo:\n")
  print(
    solo |>
      dplyr::count(
        Location,
        name = "n_anos"
      )
  )

  readr::write_csv(
    solo,
    file.path(
      PASTA_SAIDAS,
      "14_solo_harmonizado.csv"
    )
  )
}


# ------------------------------
# 18.3 Fertilização
# ------------------------------

arquivo_fert <- file.path(
  PASTA_DADOS,
  "fertilizer.xlsx"
)

if (file.exists(arquivo_fert)) {

  fertilizante <- readxl::read_excel(
    arquivo_fert
  ) |>
    dplyr::mutate(
      Location = padronizar_local(
        Location
      ),
      Treatment = limpar_sufixo_repetido(
        Treatment
      ),
      Year = as.integer(Year),
      Amount = numero_seguro(Amount)
    )

  resumo_fertilizante <- fertilizante |>
    dplyr::group_by(
      Location,
      Year,
      Treatment
    ) |>
    dplyr::summarise(
      n_aplicacoes_N = dplyr::n(),
      N_total_registrado = sum(
        Amount,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  cat("\nPrimeiras linhas do resumo de fertilização:\n")
  print(
    head(
      resumo_fertilizante,
      20
    )
  )

  readr::write_csv(
    resumo_fertilizante,
    file.path(
      PASTA_SAIDAS,
      "15_resumo_fertilizacao.csv"
    )
  )
}


# ------------------------------
# 18.4 Proteção de plantas
# ------------------------------

arquivo_pp <- file.path(
  PASTA_DADOS,
  "plant_protection.xlsx"
)

if (file.exists(arquivo_pp)) {

  protecao <- readxl::read_excel(
    arquivo_pp
  ) |>
    dplyr::mutate(
      Location = padronizar_local(
        Location
      ),
      Treatment = limpar_sufixo_repetido(
        Treatment
      ),
      Year = as.integer(Year),

      tipo_aplicacao = dplyr::case_when(
        stringr::str_detect(
          stringr::str_to_lower(
            `Plant protection`
          ),
          "fungicide"
        ) ~ "fungicide",

        stringr::str_detect(
          stringr::str_to_lower(
            `Plant protection`
          ),
          "herbicide"
        ) ~ "herbicide",

        stringr::str_detect(
          stringr::str_to_lower(
            `Plant protection`
          ),
          "insecticide"
        ) ~ "insecticide",

        stringr::str_detect(
          stringr::str_to_lower(
            `Plant protection`
          ),
          "growth regulator"
        ) ~ "growth_regulator",

        TRUE ~ "other"
      )
    )

  resumo_protecao <- protecao |>
    dplyr::count(
      Location,
      Year,
      Treatment,
      tipo_aplicacao,
      name = "n_aplicacoes"
    )

  cat("\nPrimeiras linhas do resumo de proteção de plantas:\n")
  print(
    head(
      resumo_protecao,
      20
    )
  )

  readr::write_csv(
    resumo_protecao,
    file.path(
      PASTA_SAIDAS,
      "16_resumo_protecao_plantas.csv"
    )
  )
}


# ------------------------------
# 18.5 Arquivos brutos de manejo
# ------------------------------

if (length(arquivos_manejo_brutos) > 0) {

  catalogo_manejo <- tibble::tibble(
    arquivo = basename(
      arquivos_manejo_brutos
    )
  ) |>
    tidyr::extract(
      arquivo,
      into = c(
        "Location_original",
        "Year"
      ),
      regex = "^Management_information_([A-Z]+)_([0-9]{4})\\.xlsx$",
      remove = FALSE
    ) |>
    dplyr::mutate(
      Year = as.integer(Year),
      Location = padronizar_local(
        Location_original
      )
    )

  cat("\nCatálogo dos arquivos brutos de manejo:\n")
  print(
    catalogo_manejo,
    n = nrow(catalogo_manejo)
  )

  readr::write_csv(
    catalogo_manejo,
    file.path(
      PASTA_SAIDAS,
      "17_catalogo_arquivos_manejo.csv"
    )
  )
}


# ==============================================================================
# 19. INFORMAÇÕES DOS CULTIVARES
# ==============================================================================

titulo_console(
  "15. METADADOS DOS CULTIVARES"
)

arquivo_cultivares <- file.path(
  PASTA_DADOS,
  "BRIWECs_cultivar_info.csv"
)

if (file.exists(arquivo_cultivares)) {

  cultivares <- readr::read_csv(
    arquivo_cultivares,
    show_col_types = FALSE
  )

  cat("\nDimensão da tabela de cultivares:\n")
  print(dim(cultivares))

  cat("\nResumo do ano de liberação dos cultivares:\n")
  print(
    summary(
      cultivares$RYear
    )
  )

  cat("\nNúmero de cultivares por fase:\n")
  print(
    cultivares |>
      dplyr::summarise(
        fase_I = sum(
          phaseI %in% TRUE,
          na.rm = TRUE
        ),
        fase_II = sum(
          phaseII %in% TRUE,
          na.rm = TRUE
        )
      )
  )

  readr::write_csv(
    cultivares,
    file.path(
      PASTA_SAIDAS,
      "18_metadados_cultivares.csv"
    )
  )
}


# ==============================================================================
# 20. RESUMO FINAL NO CONSOLE
# ==============================================================================

titulo_console(
  "16. RESUMO FINAL"
)

cat(
  "Análise concluída.\n\n"
)

cat(
  "Principais objetos criados:\n",
  "  dados                  -> banco fenotípico harmonizado\n",
  "  auditoria_condicoes    -> catálogo das condições disponíveis\n",
  "  dados_principal        -> subconjunto balanceado de 45 condições\n",
  "  blues                  -> BLUEs das nove características\n",
  "  pares_yield            -> R², Spearman e MI entre ambientes para produtividade\n",
  "  traits_consistencia    -> consistência das nove características\n",
  "  trait_yield            -> associação das características com produtividade\n",
  sep = ""
)

cat(
  "\nPasta com as saídas:\n",
  normalizePath(
    PASTA_SAIDAS,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "\nArquivos mais importantes para interpretação:\n",
  "  04_pares_produtividade_R2_MI.csv\n",
  "  08_testes_permutacao_rotulos_ambientes.csv\n",
  "  10_resumo_consistencia_9_traits.csv\n",
  "  12_resumo_trait_vs_yield.csv\n",
  "  fig01_NMI_corrigida_vs_R2.png\n",
  "  fig02_heatmap_NMI_corrigida.png\n",
  "  fig03_heatmap_R2_SMA.png\n",
  "  fig04_pares_MI_alta_relativa_ao_R2.png\n",
  sep = ""
)

cat(
  "\nObservação metodológica:\n",
  "A MI é estimada após discretização por postos e sua análise principal usa\n",
  "correção de viés por permutação. A sensibilidade a 4, 6 e 8 classes é\n",
  "explicitamente verificada. Os testes de diferenças entre grupos usam\n",
  "permutação dos rótulos dos ambientes, evitando tratar os pares como\n",
  "réplicas independentes em um teste convencional.\n",
  sep = ""
)

cat(
  "\nFim do script.\n"
)
