# Apostas esportivas, endividamento e consumo das famílias no Brasil

Código e dados do Trabalho de Conclusão de Curso "Os impactos da expansão das apostas esportivas no endividamento e no consumo das famílias no Brasil", de Kaio Leonardo Teixeira de Carvalho (Ciências Econômicas, Universidade Católica de Brasília, 2026).

## Conteúdo

| Arquivo | Descrição |
|---|---|
| `TCC_final.R` | Script que reproduz todas as tabelas, figuras e números do texto |
| `resultados/base_tcc.csv` | Base mensal usada no modelo (jan./2020 a jul./2026) |
| `resultados/google_trends_bruto.csv` | Coleta original do Google Trends (26/10/2026) |
| `resultados/segmentos/series_segmentos.csv` | Segmentos da PMC (tabela 8882) e da PMS (tabela 8688) |
| `resultados/final/dados_rndbf.csv` | Renda disponível bruta das famílias (SGS 29027 e 29029) |
| `resultados/variaveis_alternativas/gt_apostas.csv` | Google Trends do termo "apostas" (proxy alternativa) |
| `resultados/coletas_gt/gt_medio.csv` | Média de três coletas do Google Trends das marcas |

## Como reproduzir

1. Baixe o repositório (botão "Code" > "Download ZIP") e descompacte.
2. Abra o `TCC_final.R` no RStudio. O script define sozinho a pasta do repositório como diretório de trabalho. Fora do RStudio, rode-o a partir da própria pasta.
3. Rode `TCC_final.R`. Os pacotes necessários são instalados automaticamente. A execução leva de 10 a 15 minutos, por causa do bootstrap wild com 999 replicações.
4. As saídas são salvas em `resultados/final/`.

Por padrão, o script usa os dados salvos, que são os do texto. Com `ATUALIZAR_DADOS <- TRUE`, ele baixa novamente as séries do BCB e do IBGE. O Google Trends reamostra os dados a cada consulta, por isso a coleta original está incluída.

## Fontes

Banco Central do Brasil (SGS), IBGE (SIDRA) e Google Trends.
