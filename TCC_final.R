# ==============================================================================
# TCC — Os impactos da expansão das apostas esportivas no endividamento e
#       no consumo das famílias no Brasil
# Autor: Kaio Leonardo Teixeira de Carvalho — Ciências Econômicas / UCB
#
# Script único que reproduz todos os resultados do texto.
#
# Modelo: VARX(1) em primeiras diferenças, ordem de Cholesky GT -> PMC -> END,
# FIRs acumuladas a um choque de 1 desvio-padrão no interesse por apostas.
# Intervalos de 95% por bootstrap wild em desenho recursivo, com pesos de
# Rademacher comuns às equações e 999 replicações (Gonçalves; Kilian, 2004).
#
# Séries
#   BCB/SGS : 29038 endividamento das famílias exceto habitacional (% da renda 12m)
#             21127 inadimplência PF no cartão rotativo (%)
#             4189 Selic acumulada no mês anualizada; 433 IPCA (variação mensal)
#             24369 taxa de desocupação PNAD Contínua (sem ajuste; X-13 aqui)
#             29029 RNDBF acumulada em 12 meses; 29027 RNDBF mensal real com ajuste
#   IBGE    : PMC 8880 (varejo, volume, com ajuste sazonal)
#             PMC 8882 (segmentos do varejo) e PMS 8688 (serviços às famílias)
#   Google Trends: Bet365, Betano, Sportingbet, KTO e Superbet, Brasil, uma
#             consulta; série salva com a data da coleta e dessazonalizada (X-13)
#
# Organização
#   1. Dados (reutiliza os arquivos salvos; ATUALIZAR_DADOS = TRUE baixa de novo)
#   2. Tabela 1: estatísticas descritivas e raiz unitária
#   3. Modelo, Tabela 2 (diagnósticos, com o ARCH) e Tabela 3 (FEVD e Granger)
#   4. FIRs com bootstrap wild e Figura 1
#   5. Robustez: Tabela 4 e Figura 2
#   6. Segmentos do varejo e serviços às famílias: Tabela 5
#   7. Renda, decomposição, inadimplência e proxies (frases do texto)
#   8. Números citados no texto
# Saídas em "resultados/final/". Tempo de execução: cerca de 10 a 15 minutos.
# A dependência de abril e maio de 2020 é avaliada com dummies de impulso nesses meses.
# ==============================================================================

# 0. Configuração ------------------------------------------------------------------
rm(list = ls())
options(scipen = 999)
# Ajuste para a pasta onde o repositório foi baixado:
# setwd("caminho/para/a/pasta")

# dplyr e tidyr por último: o vars carrega o MASS, que esconderia o select().
pacotes <- c("GetBCBData", "sidrar", "gtrendsR", "seasonal", "lubridate", "vars",
             "tseries", "strucchange", "sandwich", "tidyr", "dplyr")
faltando <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(faltando) > 0) install.packages(faltando)
invisible(lapply(pacotes, library, character.only = TRUE))
select <- dplyr::select
filter <- dplyr::filter

ATUALIZAR_DADOS <- FALSE                  # FALSE = usa os dados salvos (resultados do texto)
DATA_INICIO <- as.Date("2020-01-01")
DATA_FIM    <- as.Date("2026-08-01")      # corta no último mês comum às séries
DATA_REG    <- as.Date("2025-01-01")      # início do mercado regulado
H        <- 12                            # horizonte das FIRs (meses)
R_WILD   <- 999                           # replicações do bootstrap wild
SEMENTE  <- 2026
OUT <- file.path("resultados", "final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

baixar_sgs <- function(codigos, inicio, fim) {
  x <- gbcbd_get_series(codigos, inicio, fim)
  names(x) <- gsub(".", "_", names(x), fixed = TRUE)
  x %>% transmute(data = floor_date(as.Date(ref_date), "month"), serie = series_name, value) %>%
    pivot_wider(names_from = serie, values_from = value) %>% arrange(data)
}
dessaz <- function(x, datas) {
  s <- ts(x, start = c(year(min(datas)), month(min(datas))), frequency = 12)
  as.numeric(final(seas(s)))
}

# 1. Dados ------------------------------------------------------------------------
# 1.1 Base principal
arq_base <- "resultados/base_tcc.csv"
if (!ATUALIZAR_DADOS && file.exists(arq_base)) {
  base <- read.csv(arq_base); base$data <- as.Date(base$data)
} else {
  bacen <- baixar_sgs(c(END = 29038, INAD = 21127, SELIC = 4189, IPCA = 433),
                      DATA_INICIO, DATA_FIM + 30)
  # Desocupação desde 2015, para o X-13 estimar os fatores sazonais
  desoc <- baixar_sgs(c(DESOC_NSA = 24369), as.Date("2015-01-01"), DATA_FIM + 30) %>%
    filter(!is.na(DESOC_NSA))
  desoc$DESOC <- dessaz(desoc$DESOC_NSA, desoc$data)
  pmc <- get_sidra(api = "/t/8880/n1/all/v/7170/p/all/c11046/56734") %>%
    transmute(data = as.Date(paste0(`Mês (Código)`, "01"), "%Y%m%d"), PMC = Valor)
  # Google Trends: reutiliza a coleta salva (o Google reamostra a cada consulta)
  arq_gt <- "resultados/google_trends_bruto.csv"
  if (file.exists(arq_gt)) {
    gt_bruto <- read.csv(arq_gt, stringsAsFactors = FALSE)
  } else {
    gt <- gtrends(keyword = c("Bet365", "Betano", "Sportingbet", "KTO", "Superbet"),
                  geo = "BR", time = paste(DATA_INICIO, DATA_FIM + 30), onlyInterest = TRUE)
    gt_bruto <- gt$interest_over_time %>%
      transmute(date = as.character(as.Date(date)), keyword,
                hits = as.numeric(gsub("<1", "0.5", hits)), data_coleta = as.character(Sys.Date()))
    write.csv(gt_bruto, arq_gt, row.names = FALSE)
  }
  gt_m <- gt_bruto %>% mutate(data = floor_date(as.Date(date), "month")) %>%
    group_by(data) %>% summarise(GT_NSA = sum(hits), .groups = "drop") %>% arrange(data)
  gt_m$GT <- dessaz(gt_m$GT_NSA, gt_m$data)
  base <- bacen %>% inner_join(select(desoc, data, DESOC), by = "data") %>%
    inner_join(pmc, by = "data") %>% inner_join(gt_m, by = "data") %>%
    filter(data >= DATA_INICIO, data <= DATA_FIM) %>% arrange(data) %>%
    mutate(ln_GT = log(GT + 1), ln_PMC = log(PMC),
           D_COVID = as.numeric(data >= as.Date("2020-03-01") & data <= as.Date("2020-12-01")),
           D_TETO  = as.numeric(data >= as.Date("2024-01-01")),
           D_REG   = as.numeric(data >= DATA_REG)) %>%
    drop_na()
  write.csv(base, arq_base, row.names = FALSE)
}

# 1.2 Renda disponível bruta das famílias (RNDBF)
arq_rnd <- file.path(OUT, "dados_rndbf.csv")
if (!ATUALIZAR_DADOS && file.exists(arq_rnd)) {
  rnd <- read.csv(arq_rnd); rnd$data <- as.Date(rnd$data)
} else {
  rnd <- baixar_sgs(c(RNDBF12 = 29029, RNDBF_SA = 29027), min(base$data), max(base$data) + 30)
  write.csv(rnd, arq_rnd, row.names = FALSE)
}

# 1.3 Segmentos do varejo (PMC 8882) e serviços às famílias (PMS 8688)
SEGMENTOS <- c(S1 = "Essencial: hipermercados e supermercados",
               S2 = "Discricionário: móveis e eletrodomésticos",
               S3 = "Discricionário: tecidos, vestuário e calçados",
               S4 = "Discricionário: outros artigos de uso pessoal",
               S5 = "Lazer: serviços prestados às famílias",
               S6 = "Lazer: alojamento e alimentação",
               S7 = "Lazer: outros serviços às famílias")
arq_seg <- "resultados/segmentos/series_segmentos.csv"
if (!ATUALIZAR_DADOS && file.exists(arq_seg)) {
  segs <- read.csv(arq_seg, stringsAsFactors = FALSE); segs$data <- as.Date(segs$data)
} else {
  dir.create(dirname(arq_seg), recursive = TRUE, showWarnings = FALSE)
  baixa_tab <- function(tab) do.call(rbind, lapply(   # o SIDRA limita cada consulta a 50 mil valores
    c("201601-201912", "202001-202212", "202301-202512", "202601-202712"),
    function(p) tryCatch(get_sidra(x = tab, period = p, geo = "Brazil"), error = function(e) NULL)))
  serie <- function(x, padrao, rotulo) {
    ca <- names(x)[grepl("tividade", names(x)) & !grepl("Código", names(x))]
    y <- x[grepl(padrao, x[[ca]]) & grepl("Número-índice", x[["Variável"]]), ]
    if ("Tipos de índice" %in% names(y)) y <- y[grepl("volume", y[["Tipos de índice"]], ignore.case = TRUE), ]
    sa <- y[grepl("com ajuste", y[["Variável"]]), ]
    usar_sa <- nrow(sa) > 0
    y <- if (usar_sa) sa else y[!grepl("ajuste", y[["Variável"]]), ]
    y <- y %>% transmute(data = as.Date(paste0(`Mês (Código)`, "01"), "%Y%m%d"),
                         valor = as.numeric(Valor)) %>% filter(!is.na(valor)) %>% arrange(data)
    if (!usar_sa) y$valor <- dessaz(y$valor, y$data)
    y$segmento <- rotulo
    y
  }
  pmc_s <- baixa_tab(8882); pms_s <- baixa_tab(8688)
  segs <- rbind(serie(pmc_s, "^Hipermercados, supermercados", SEGMENTOS[["S1"]]),
                serie(pmc_s, "^Móveis e eletrodomésticos",    SEGMENTOS[["S2"]]),
                serie(pmc_s, "^Tecidos, vestuário",           SEGMENTOS[["S3"]]),
                serie(pmc_s, "^Outros artigos de uso pessoal",SEGMENTOS[["S4"]]),
                serie(pms_s, "^1\\. Serviços prestados às famílias", SEGMENTOS[["S5"]]),
                serie(pms_s, "^1\\.1 Serviços de alojamento",        SEGMENTOS[["S6"]]),
                serie(pms_s, "^1\\.2 Outros serviços prestados",     SEGMENTOS[["S7"]]))
  write.csv(segs, arq_seg, row.names = FALSE)
}
seg_largo <- segs %>% mutate(cod = names(SEGMENTOS)[match(segmento, SEGMENTOS)]) %>%
  filter(!is.na(cod)) %>% select(data, cod, valor) %>%
  pivot_wider(names_from = cod, values_from = valor)

# 1.4 Base de estimação (variáveis em diferença)
dif <- function(x) c(NA, diff(x))
d <- base %>% left_join(rnd, by = "data") %>% left_join(seg_largo, by = "data") %>%
  arrange(data) %>%
  mutate(P       = cumprod(1 + IPCA / 100), P = P / last(P),   # índice do IPCA (deflator)
         ln_DIV  = log(END / 100 * RNDBF12 / P),               # dívida real
         ln_REND = log(RNDBF12 / P),                           # renda real acumulada em 12 meses
         ln_RSA  = log(RNDBF_SA))                              # renda mensal real dessazonalizada
df <- d %>%
  mutate(d_ln_GT = dif(ln_GT), d_ln_PMC = dif(ln_PMC), d_END = dif(END),
         d_SELIC = dif(SELIC), d_DESOC = dif(DESOC), d_INAD = dif(INAD),
         d_ln_DIV = dif(ln_DIV), d_ln_REND = dif(ln_REND), d_ln_RSA = dif(ln_RSA),
         across(all_of(names(SEGMENTOS)), ~ dif(log(.x)), .names = "d_ln_{.col}"),
         D_IMP2020   = as.numeric(data %in% as.Date(c("2020-04-01", "2020-05-01"))),
         D_COPA      = as.numeric(data %in% as.Date(c("2022-11-01", "2022-12-01",
                                                      "2026-06-01", "2026-07-01"))),
         D_DESENROLA = as.numeric(data >= as.Date("2023-07-01")),
         D_BLOQ      = as.numeric(data >= as.Date("2024-10-01"))) %>%
  slice(-1)

ENDOG <- c("d_ln_GT", "d_ln_PMC", "d_END")                    # ordem de Cholesky
EXOG  <- c("d_SELIC", "IPCA", "d_DESOC", "D_COVID", "D_TETO", "D_REG")

cat("Amostra:", format(min(base$data)), "a", format(max(base$data)), "| T =", nrow(base),
    "| observações efetivas =", nrow(df) - 1, "\n")
if (file.exists("resultados/google_trends_bruto.csv"))
  cat("Data da coleta do Google Trends:",
      unique(read.csv("resultados/google_trends_bruto.csv")$data_coleta), "\n")

# 2. Tabela 1: descritivas e raiz unitária -----------------------------------------
# ADF com constante e tendência, defasagens trunc((n-1)^(1/3)) = 4; KPSS em nível
# e em torno de tendência, defasagem "short" = 3. Em nível, a série é não
# estacionária se o ADF não rejeitar a 5% OU o KPSS rejeitar; em diferença,
# basta um dos testes indicar estacionariedade.
raiz <- function(x) c(
  ADF       = suppressWarnings(adf.test(x))$p.value,
  KPSS      = suppressWarnings(kpss.test(x, null = "Level"))$p.value,
  KPSS_tend = suppressWarnings(kpss.test(x, null = "Trend"))$p.value,
  ADF_d     = suppressWarnings(adf.test(diff(x)))$p.value,
  KPSS_d    = suppressWarnings(kpss.test(diff(x), null = "Level"))$p.value)
nivel <- c(GT = "ln_GT", PMC = "ln_PMC", END = "END", SELIC = "SELIC", IPCA = "IPCA", DESOC = "DESOC")
tab1 <- do.call(rbind, lapply(names(nivel), function(v) {
  x <- base[[v]]; r <- raiz(base[[nivel[[v]]]])
  data.frame(variavel = v, media = mean(x), dp = sd(x), minimo = min(x), maximo = max(x),
             t(round(r, 3)),
             tratamento = ifelse(r[["ADF"]] > 0.05 | r[["KPSS"]] < 0.05, "1ª diferença", "nível"))
}))
tab1[, 2:5] <- round(tab1[, 2:5], 2)
print(tab1, row.names = FALSE)
write.csv(tab1, file.path(OUT, "tab1_descritivas_raiz.csv"), row.names = FALSE)


# 3. Modelo e diagnósticos -----------------------------------------------------------
prepara <- function(dados, endog, exog) {
  dados <- dados[complete.cases(dados[, c(endog, exog)]), ]
  X <- as.matrix(dados[, exog, drop = FALSE])
  list(dados = dados, Y = as.matrix(dados[, endog]),
       X = X[, apply(X, 2, var) > 0, drop = FALSE])          # retira dummies constantes
}
estima <- function(Y, X, tipo = "const", p = 1)
  do.call(VAR, list(y = Y, p = p, type = tipo, exogen = X))

pp <- prepara(df, ENDOG, EXOG)
sel <- VARselect(pp$Y, lag.max = 6, type = "const", exogen = pp$X)
print(sel$selection)
p <- as.integer(sel$selection["HQ(n)"])
repeat {                                                     # Hannan-Quinn + Breusch-Godfrey
  m0 <- estima(pp$Y, pp$X, p = p)
  if (serial.test(m0, lags.bg = 2, type = "BG")$serial$p.value > 0.05 || p >= 6) break
  p <- p + 1
}
cat("Defasagem: p =", p, "\n")
if (p != 1) stop("O bootstrap wild abaixo foi escrito para p = 1.")

bg1  <- serial.test(m0, lags.bg = 1, type = "BG")$serial
bg2  <- serial.test(m0, lags.bg = 2, type = "BG")$serial
jb   <- normality.test(m0)$jb.mul$JB
arch <- arch.test(m0, lags.multi = 2, multivariate.only = TRUE)$arch.mul
est  <- stability(m0, type = "OLS-CUSUM")
tab2 <- rbind(
  data.frame(teste = c("Maior raiz inversa (módulo)", "Breusch-Godfrey LM (1 defasagem)",
                       "Breusch-Godfrey LM (2 defasagens)", "Jarque-Bera multivariado",
                       "ARCH multivariado (2 defasagens)"),
             estatistica = round(c(max(roots(m0)), bg1$statistic, bg2$statistic,
                                   jb$statistic, arch$statistic), 3),
             p_valor = c(NA, round(c(bg1$p.value, bg2$p.value, jb$p.value, arch$p.value), 3))),
  data.frame(teste = paste0("OLS-CUSUM – equação ", names(est$stability)),
             estatistica = sapply(est$stability, function(e) round(sctest(e)$statistic, 3)),
             p_valor = sapply(est$stability, function(e) round(sctest(e)$p.value, 3)),
             row.names = NULL))
print(tab2, row.names = FALSE)
write.csv(tab2, file.path(OUT, "tab2_diagnosticos.csv"), row.names = FALSE)

# Tabela 3: decomposição da variância e Granger
dv <- fevd(m0, n.ahead = H)
tab3 <- data.frame(variavel = ENDOG[2:3],
                   t(sapply(ENDOG[2:3], function(r) round(100 * dv[[r]][c(1, 3, 6, 12), "d_ln_GT"], 2))),
                   row.names = NULL)
names(tab3)[-1] <- paste0("FEVD_h", c(1, 3, 6, 12))
granger <- do.call(rbind, lapply(ENDOG[2:3], function(eq) {
  full <- m0$varresult[[eq]]; Xf <- model.matrix(full); y <- full$model[[1]]
  lg <- startsWith(colnames(Xf), "d_ln_GT.l")
  rss_f <- sum(residuals(full)^2); rss_r <- sum(lm.fit(Xf[, !lg, drop = FALSE], y)$residuals^2)
  gl <- nrow(Xf) - ncol(Xf); F <- ((rss_r - rss_f) / sum(lg)) / (rss_f / gl)
  data.frame(variavel = eq, Granger_F = round(F, 3),
             Granger_p = round(pf(F, sum(lg), gl, lower.tail = FALSE), 4))
}))
tab3 <- merge(tab3, granger, by = "variavel", sort = FALSE)
gr <- causality(m0, cause = "d_ln_GT")$Granger
cat("Granger conjunto: F =", round(gr$statistic, 3), "| p =", round(gr$p.value, 4), "\n")
print(tab3, row.names = FALSE)
write.csv(tab3, file.path(OUT, "tab3_fevd_granger.csv"), row.names = FALSE)

# 4. FIRs com bootstrap wild ------------------------------------------------------------
# FIR pontual; respostas em diferença são acumuladas (efeito sobre o nível).
fir_pontual <- function(m, imp, resp, acum) {
  f <- irf(m, impulse = imp, response = resp, n.ahead = H, ortho = TRUE,
           cumulative = FALSE, boot = FALSE)$irf[[1]][, resp, drop = FALSE]
  for (k in seq_along(resp)) if (acum[k]) f[, k] <- cumsum(f[, k])
  f
}
# Bootstrap wild recursivo: em cada replicação, um sinal (+1 ou -1) por mês, o
# mesmo nas três equações; as séries são reconstruídas mês a mês com os
# coeficientes estimados, as mesmas exógenas e os resíduos com o sinal sorteado;
# o VAR é reestimado e as FIRs recalculadas. IC 95% = percentis 2,5% e 97,5%.
wild <- function(Y, X, tipo, imp, resp, acum) {
  set.seed(SEMENTE)
  m <- estima(Y, X, tipo)
  det <- if (tipo == "both") cbind(const = 1, trend = seq_len(nrow(Y))) else cbind(const = rep(1, nrow(Y)))
  B <- Bcoef(m)[, c(paste0(colnames(Y), ".l1"), colnames(det), colnames(X))]
  U <- residuals(m)
  sims <- replicate(R_WILD, {
    eta <- sample(c(-1, 1), nrow(U), replace = TRUE)
    Ys <- Y
    for (t in 2:nrow(Y)) Ys[t, ] <- as.numeric(B %*% c(Ys[t - 1, ], det[t, ], X[t, ])) + eta[t - 1] * U[t - 1, ]
    fir_pontual(estima(Ys, X, tipo), imp, resp, acum)
  }, simplify = "array")                                     # (H+1) x respostas x R
  list(m = m, sims = sims)
}
roda <- function(rotulo, dados, endog, exog, resp, acum = rep(TRUE, length(resp)),
                 tipo = "const", imp = "d_ln_GT") {
  pr <- prepara(dados, endog, exog)
  w  <- wild(pr$Y, pr$X, tipo, imp, resp, acum)
  pt <- fir_pontual(w$m, imp, resp, acum)
  fir <- do.call(rbind, lapply(seq_along(resp), function(k) data.frame(
    especificacao = rotulo, resposta = resp[k], horizonte = 0:H, fir = pt[, k],
    li = apply(w$sims[, k, , drop = FALSE], 1, quantile, 0.025),
    ls = apply(w$sims[, k, , drop = FALSE], 1, quantile, 0.975),
    n_obs = nrow(pr$Y) - 1, row.names = NULL)))
  fir$significativo <- fir$li > 0 | fir$ls < 0
  cat("OK:", rotulo, "| T efetivo =", nrow(pr$Y) - 1, "\n")
  list(fir = fir, modelo = w$m)
}
# Resumo: valor em 12 meses (log x 100 = %; END e INAD em p.p.) e meses significantes
resumir <- function(fir) fir %>%
  mutate(especificacao = factor(especificacao, levels = unique(especificacao)),
         k = ifelse(grepl("ln_", resposta), 100, 1)) %>%
  group_by(especificacao, resposta) %>%
  summarise(h12 = round(first(k) * fir[horizonte == 12], 3),
            li  = round(first(k) * li[horizonte == 12], 3),
            ls  = round(first(k) * ls[horizonte == 12], 3),
            meses_signif = ifelse(any(significativo), paste(horizonte[significativo], collapse = ","), "nenhum"),
            .groups = "drop") %>% as.data.frame()

cat("\nEstimando com bootstrap wild (", R_WILD, " replicações)...\n", sep = "")
esp <- list()
esp$principal <- roda("Principal", df, ENDOG, EXOG, c("d_ln_PMC", "d_END"))

# Tamanho e persistência do choque (resposta do próprio GT)
g <- irf(esp$principal$modelo, impulse = "d_ln_GT", response = "d_ln_GT", n.ahead = H,
         ortho = TRUE, cumulative = FALSE, boot = FALSE)$irf[[1]][, 1]
choque_pct <- 100 * (exp(g[1]) - 1); persist_pct <- 100 * (exp(sum(g)) - 1)

# Figura 1
h <- 0:H
fp <- esp$principal$fir
rot <- c(d_ln_PMC = "Varejo (PMC)", d_END = "Endividamento (END)")
painel <- function(r, extra = NULL, k = ifelse(r == "d_ln_PMC", 100, 1)) {
  x <- fp[fp$resposta == r, ]
  plot(h, k * x$fir, type = "n", ylim = range(k * x$li, k * x$ls, extra, 0),
       xlab = "Meses após o choque",
       ylab = ifelse(r == "d_ln_PMC", "Resposta acumulada (%)", "Resposta acumulada (p.p.)"),
       main = rot[[r]])
  polygon(c(h, rev(h)), c(k * x$li, rev(k * x$ls)), col = "grey85", border = NA)
  abline(h = 0, lty = 2, col = "grey40")
  x
}
png(file.path(OUT, "fig1_fir.png"), width = 2000, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (r in c("d_ln_PMC", "d_END")) { x <- painel(r); lines(h, ifelse(r == "d_ln_PMC", 100, 1) * x$fir, lwd = 2) }
dev.off()

# 5. Robustez ----------------------------------------------------------------------------
esp$invertida <- roda("Ordem invertida", df, c("d_ln_PMC", "d_END", "d_ln_GT"), EXOG, c("d_ln_PMC", "d_END"))
# IRF generalizada: com o GT em 1º, é idêntica à Cholesky em cada replicação
# (conferência numérica mais abaixo), então repete os valores do modelo principal.
esp$girf <- list(fir = transform(esp$principal$fir, especificacao = "IRF generalizada"))
esp$pre2025   <- roda("Subamostra 2020-2024", filter(df, data < DATA_REG), ENDOG, EXOG, c("d_ln_PMC", "d_END"))
esp$impulso   <- roda("Impulso abr-mai/2020", df, ENDOG, c(EXOG, "D_IMP2020"), c("d_ln_PMC", "d_END"))
esp$nivel     <- roda("Varejo em nível (VAR com tendência)", df, c("d_ln_GT", "ln_PMC", "d_END"), EXOG,
                      c("ln_PMC", "d_END"), acum = c(FALSE, TRUE), tipo = "both")
esp$copa      <- roda("Dummies de Copa do Mundo", df, ENDOG, c(EXOG, "D_COPA"), c("d_ln_PMC", "d_END"))
esp$dummies   <- roda("Dummies separadas (Desenrola, teto, bloqueio)", df, ENDOG,
                      c("d_SELIC", "IPCA", "d_DESOC", "D_COVID", "D_DESENROLA", "D_TETO", "D_BLOQ", "D_REG"),
                      c("d_ln_PMC", "d_END"))

# IRF generalizada (Pesaran; Shin, 1998): com o GT em 1º, coincide com a Cholesky
P0 <- Psi(m0, nstep = H); Ph <- Phi(m0, nstep = H)
S  <- P0[, , 1] %*% t(P0[, , 1]); j <- 1
girf <- sapply(0:H, function(k) Ph[, , k + 1] %*% S[, j] / sqrt(S[j, j]))
chol <- sapply(0:H, function(k) P0[, j, k + 1])
cat("Diferença máxima IRF generalizada x Cholesky:", signif(max(abs(girf - chol)), 3), "\n")

# Projeções locais (Jordà, 2005), erros-padrão de Newey-West
choque <- residuals(m0)[, "d_ln_GT"]; choque <- choque / sd(choque)
dd <- pp$dados; n <- length(choque)
ctrl <- cbind(dd[1:n, ENDOG], dd[2:(n + 1), EXOG]); names(ctrl)[1:3] <- paste0(ENDOG, "_l1")
idx <- match(dd$data[2:(n + 1)], base$data)
lp <- do.call(rbind, lapply(c("ln_PMC", "END"), function(v) do.call(rbind, lapply(0:H, function(k) {
  x <- data.frame(y = base[[v]][idx + k] - base[[v]][idx - 1], choque = choque, ctrl)
  x <- x[complete.cases(x), ]
  x <- x[, c(TRUE, TRUE, apply(x[, -(1:2), drop = FALSE], 2, var) > 0)]
  fit <- lm(y ~ ., data = x)
  se <- sqrt(diag(NeweyWest(fit, lag = k + 1, prewhite = FALSE)))[["choque"]]
  b <- coef(fit)[["choque"]]
  data.frame(especificacao = "Projeções locais", resposta = ifelse(v == "END", "d_END", "d_ln_PMC"),
             horizonte = k, fir = b, li = b - 1.96 * se, ls = b + 1.96 * se, n_obs = nrow(x))
}))))
lp$significativo <- lp$li > 0 | lp$ls < 0

rob <- rbind(do.call(rbind, lapply(esp, `[[`, "fir")), lp)
rob$resposta[rob$resposta == "ln_PMC"] <- "d_ln_PMC"       # varejo em nível: mesma coluna da tabela
write.csv(rob, file.path(OUT, "fir_robustez.csv"), row.names = FALSE)
tab4 <- resumir(rob) %>%
  mutate(resposta = ifelse(resposta == "d_END", "Endividamento (p.p.)", "Varejo (%)"))
print(tab4, row.names = FALSE)
write.csv(tab4, file.path(OUT, "tab4_robustez.csv"), row.names = FALSE)

# Figura 2: principal (IC 95%) x ordem invertida x subamostra 2020-2024
png(file.path(OUT, "fig2_fir_robustez.png"), width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), oma = c(2.5, 0, 0, 0))
for (r in c("d_ln_PMC", "d_END")) {
  k <- ifelse(r == "d_ln_PMC", 100, 1)
  curvas <- k * sapply(c("principal", "invertida", "pre2025"), function(e) {
    x <- esp[[e]]$fir; x$fir[x$resposta == r] })
  painel(r, extra = curvas)
  matlines(h, curvas, lty = 1:3, lwd = 2, col = "black")
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", horiz = TRUE, bty = "n", lty = 1:3, lwd = 2, cex = 0.8,
       legend = c("Modelo principal (IC 95%)", "Ordem invertida", "Subamostra 2020-2024"))
dev.off()

# 6. Segmentos do varejo e serviços às famílias ----------------------------------------
# Cada segmento substitui o varejo total, na mesma posição da ordem de Cholesky
# (principal) e com o interesse por apostas por último (ordem invertida).
seg <- list()
for (s in names(SEGMENTOS)) {
  v <- paste0("d_ln_", s)
  seg[[length(seg) + 1]] <- roda(paste(SEGMENTOS[[s]], "| principal"), df,
                                 c("d_ln_GT", v, "d_END"), EXOG, v)$fir
  seg[[length(seg) + 1]] <- roda(paste(SEGMENTOS[[s]], "| ordem invertida"), df,
                                 c(v, "d_END", "d_ln_GT"), EXOG, v)$fir
}
tab5 <- resumir(do.call(rbind, seg))
print(tab5, row.names = FALSE)
write.csv(tab5, file.path(OUT, "tab5_segmentos.csv"), row.names = FALSE)

# Segmentos com queda na ordem invertida: pandemia e Copas (frase do texto)
seg_extra <- list()
for (s in c("S2", "S4")) {
  v <- paste0("d_ln_", s); inv <- c(v, "d_END", "d_ln_GT")
  seg_extra[[length(seg_extra) + 1]] <- roda(paste(SEGMENTOS[[s]], "| invertida + impulso 2020"), df,
                                             inv, c(EXOG, "D_IMP2020"), v)$fir
  seg_extra[[length(seg_extra) + 1]] <- roda(paste(SEGMENTOS[[s]], "| invertida + Copa"), df,
                                             inv, c(EXOG, "D_COPA"), v)$fir
}
seg_extra <- resumir(do.call(rbind, seg_extra))

# 7. Renda, decomposição, inadimplência e proxies ----------------------------------------
comp <- list(
  roda("END com controle de renda (RNDBF)", df, ENDOG, c(EXOG, "d_ln_RSA"), "d_END")$fir,
  roda("Dívida real", df, c("d_ln_GT", "d_ln_PMC", "d_ln_DIV"), EXOG, "d_ln_DIV")$fir,
  roda("Renda real (12 meses)", df, c("d_ln_GT", "d_ln_PMC", "d_ln_REND"), EXOG, "d_ln_REND")$fir,
  roda("Inadimplência do rotativo", df, c("d_ln_GT", "d_ln_PMC", "d_INAD"), EXOG, "d_INAD")$fir)

# Proxies alternativas (usa as coletas salvas, se existirem)
proxy <- function(rotulo, gt) {          # gt: data.frame com data e ln_GT da proxy
  x <- df %>% select(-any_of("d_ln_GTx")) %>%
    left_join(gt %>% arrange(data) %>% transmute(data, d_ln_GTx = c(NA, diff(ln_GT))), by = "data")
  roda(rotulo, x, c("d_ln_GTx", "d_ln_PMC", "d_END"), EXOG, "d_END", imp = "d_ln_GTx")$fir
}
arq_ap <- "resultados/variaveis_alternativas/gt_apostas.csv"
if (file.exists(arq_ap)) {
  ap <- read.csv(arq_ap) %>% mutate(data = floor_date(as.Date(date), "month")) %>%
    group_by(data) %>% summarise(GT = sum(hits), .groups = "drop") %>% arrange(data)
  ap$ln_GT <- log(dessaz(ap$GT, ap$data) + 1)
  comp[[length(comp) + 1]] <- proxy("Proxy: termo \"apostas\"", ap)
}
arq_med <- "resultados/coletas_gt/gt_medio.csv"
if (file.exists(arq_med)) {
  med <- read.csv(arq_med); med$data <- as.Date(med$data)
  comp[[length(comp) + 1]] <- proxy("Proxy: média das coletas do GT", med)
}
comp <- resumir(do.call(rbind, comp))
print(rbind(seg_extra, comp), row.names = FALSE)
write.csv(rbind(seg_extra, comp), file.path(OUT, "complementares.csv"), row.names = FALSE)

# 8. Números citados no texto ------------------------------------------------------------
rnd12 <- tail(na.omit(d$RNDBF12), 1)                        # R$ milhões, último mês disponível
e12 <- tab4[tab4$especificacao == "Principal" & tab4$resposta == "Endividamento (p.p.)", ]
fe  <- fp[fp$resposta == "d_END", ]
txt <- c(
  sprintf("Amostra: %s a %s; T = %d; observações efetivas = %d",
          format(min(base$data)), format(max(base$data)), nrow(base), nrow(df) - 1),
  sprintf("Choque de 1 dp: +%.1f%% no nível das buscas; %.1f%% persistem após 12 meses",
          choque_pct, persist_pct),
  sprintf("END principal em 12 meses: %.3f p.p. (IC 95%%: %.3f a %.3f); meses significantes: %s",
          e12$h12, e12$li, e12$ls, e12$meses_signif),
  sprintf("END até o mês 3: %.3f p.p.", fe$fir[fe$horizonte == 3]),
  sprintf("Magnitude: %.1f bi de reais de dívida (IC: %.1f a %.1f), com RNDBF 12m = %.0f milhões",
          e12$h12 / 100 * rnd12 / 1000, e12$li / 100 * rnd12 / 1000, e12$ls / 100 * rnd12 / 1000, rnd12),
  sprintf("ARCH multivariado: estatística = %.2f; p = %.4f", arch$statistic, arch$p.value),
  sprintf("Granger conjunto: F = %.3f; p = %.4f", gr$statistic, gr$p.value))
writeLines(txt); writeLines(txt, file.path(OUT, "numeros_texto.txt"))

cat("\nArquivos salvos em:", normalizePath(OUT), "\n")
