# MIBLUP: Seleção Genómica Ponderada por Informação Mútua no Melhoramento de Plantas 🌱

![Status](https://img.shields.io/badge/Status-Concluído-success)
![Prémio](https://img.shields.io/badge/Prémio-1º_Lugar_ESALQ_Desafio_GDM-gold)
![Genética Quantitativa](https://img.shields.io/badge/Domínio-Genética_Quantitativa-blue)

## 📌 Visão Geral

O **MIBLUP** é um método inovador de seleção genómica que utiliza a **Teoria da Informação de Shannon** (Informação Mútua) para filtrar o ruído ambiental e isolar o verdadeiro potencial genético das plantas. 

Esta solução otimiza o melhoramento vegetal ao aumentar drasticamente a precisão na deteção de marcadores genéticos (SNPs que realmente controlam características agronómicas importantes, como a produtividade), sem incorrer nos elevados custos computacionais e na falta de transparência típicos da inteligência artificial complexa (Deep Learning).

Este projeto foi o **grande vencedor do 1º lugar (melhor trabalho da ESALQ)** no **Desafio GDM Universidades!** 🏆

## 🚀 Principais Diferenciais e Resultados

- **Filtro de Ruído (Gargalo de Informação):** O algoritmo atua de forma análoga ao princípio do *Information Bottleneck* das redes neurais, esmagando matematicamente o peso de SNPs inúteis e exponenciando a importância dos SNPs que sistematicamente acertam a classe de produtividade.
- **Fator de Enriquecimento Elevado:** O modelo demonstrou um fator de enriquecimento médio de **9,125**, atingindo picos de até **12,58**. Isto significa que o MIBLUP deteta os *loci* de influência (QTLs verdadeiros) a uma taxa cerca de 9 a 12,5 vezes maior do que ocorreria ao acaso. Um verdadeiro "garimpo biológico" preciso.
- **Robustez Estatística:** Validado através de simulação estocástica de Monte Carlo com 100 repetições independentes, incluindo testes de *stress* com forte epistasia multiplicativa.
- **Interpretabilidade Biológica (Fim da Caixa-Preta):** Ao contrário de modelos de IA opacos, a ponderação por Informação Mútua permite rastrear fisicamente as regiões mais importantes do genoma, oferecendo total transparência aos melhoristas de campo.

## 🔬 Metodologia

O MIBLUP ancora-se na Equação do Melhorista para maximizar o ganho genético ($\Delta G$). O fluxo de trabalho inclui:
1. **Discretização:** A variável alvo (ex: produtividade) é dividida em classes (ex: 10).
2. **Cálculo da Informação Mútua (IM):** Avalia-se a redução da incerteza sobre o fenótipo ao conhecer o genótipo de cada SNP.
3. **Ponderação Genómica:** SNPs com alta IM (QTLs causais/epistáticos) recebem pesos elevados na matriz genómica ($G_{MI}$), enquanto os ruídos recebem pesos próximos de zero.
4. **Previsão:** Resolução das equações de modelos mistos utilizando a matriz otimizada.

## 💻 Como Utilizar

```bash
# Exemplo de clonagem do repositório
git clone [https://github.com/seu-usuario/miblup-genomics.git](https://github.com/seu-usuario/MIBLUP.git)
cd MIBLUP

# Executar o script principal (Exemplo em R)
Rscript Simulacao_MI.R 
