# CONVENTIONS.md — as regras do Daedalus

Especificação normativa do laboratório de caminhadas quânticas de tempo
contínuo. **Este documento manda no código.** Se o código discordar daqui, o
código está errado.

Cada regra abaixo existe porque quebrá-la produz resultado *plausível e
errado*, não erro visível. Onde há um teste que segura a regra, ele está
nomeado.

---

## 1. Um único núcleo

Existe **um** núcleo, em `core/`. Ele vira três coisas:

1. **WebAssembly**, por `emcc`, para o navegador;
2. **binário nativo**, por `gcc`/`clang`, para testes e benchmarks;
3. o **`.cpp` exportado**, por concatenação do mesmo núcleo com um `main` de
   template.

`tools/amalgamate.mjs` lê `core/amalgam.list`, remove os `#include "dae_*.h"`
internos e emite um tradutor único. É esse texto — e não uma segunda
implementação — que vai para dentro do WASM e para dentro do `.cpp`. Por isso o
teste de aceitação 7 compara **o mesmo código sob dois compiladores**, e não
duas implementações que concordam hoje e divergem no mês que vem.

Consequências obrigatórias, verificadas por `make cxx-check`:

- o núcleo é **C99 que também é C++17 válido**: sem `_Complex`, sem VLA, sem
  inicializador designado, sem literal composto, todo `malloc` convertido
  explicitamente;
- **zero estado global mutável** — o template exportado roda realizações sob
  `#pragma omp parallel for`, e qualquer estado compartilhado seria corrida.
  Toda função do núcleo é reentrante e recebe seu contexto por argumento;
- sem `<stdio.h>` no núcleo. Impressão só em `native/` e nos templates;
- **`-ffast-math` é proibido**. Quebraria o NaN de `p_alvo` (parte 6) e a
  reprodutibilidade bit a bit.

## 2. O `spec.json` é analisado só em C

O TypeScript **monta** o texto JSON; quem o **interpreta** é `dae_spec.c`, do
lado C. O `.cpp` exportado embute o mesmo texto num *raw string literal* e usa
o mesmo parser.

A interface valida chamando o parser via WASM — é barato — e exibe
`dae_error.line/col`. **Não existe validação paralela em TypeScript.** Se
existisse, um valor-padrão divergiria entre os dois lados, e a divergência
apareceria como discrepância física meses depois.

O parser (`core/dae_spec.c`) é **estrito**: chave desconhecida é ERRO, com
linha e coluna. Um emissor que escrevesse `seam` onde o formato diz
`seam_shift` produziria, com um parser tolerante, um arquivo internamente
consistente resolvendo o problema **errado** — e o oráculo de verificação não
pegaria, porque ele valida o propagador, não o emissor. Campo novo exige subir
`format_version`.

O `spec.json` é a única entrada da ponte WASM. Parâmetros atravessando a ponte
num vetor de `double` com posições fixas são uma fresta por onde os padrões
entram zerados sem aviso (parte 10.2).

**Emenda à regra 3 do núcleo:** `dae_spec.c` e `dae_csv.c` incluem `<stdio.h>`,
só por `snprintf`. Formatar para um buffer não é entrada e saída — nenhum
`printf`, nenhum arquivo. A alternativa seria serializar fora do núcleo, e aí o
JSON canônico existiria duas vezes: uma no TypeScript e outra no `main` do
`.cpp`. Duas versões do canônico é exatamente a divergência que a amalgamação
existe para impedir.

## 3. As duas disciplinas metodológicas

Não são conforto de interface. Sem elas os resultados são artefato.

### 3.1 Religar, não acrescentar

Acrescentar arestas de longo alcance aumenta `|E|` e melhora o transporte
trivialmente — o resultado mede o número de arestas, não a topologia. O modo
padrão do controle de conectividade (`DAE_REWIRE`) mantém `|E|` **fixo** e
religa. `DAE_ADD` existe, e a interface avisa explicitamente quando está ligado.

### 3.2 Fixar a escala de energia

Ao comparar topologias diferentes é obrigatório normalizar `‖H‖`, senão "mais
coerência" pode ser apenas "hopping maior". A normalização vem **ligada** por
padrão (`DAE_NORM_SPECTRAL` ou `DAE_NORM_MEAN_DEGREE`); desligar é escolha
explícita e fica registrada no `spec.json`.

Semântica da escala: o fator é medido em `H_cru` com `gamma = 1`, e então
`H = gamma · H_cru / escala`. Assim `gamma` continua sendo o botão físico e a
normalização remove só a escala imposta pela topologia. O fator aplicado sai em
`scale_out`, aparece na interface e entra no cabeçalho do CSV.

## 4. Reprodutibilidade

- **O PRNG mora no núcleo** (`dae_rng`, xoshiro256++ semeado por splitmix64).
  Nunca `Math.random()` do JavaScript: o grafo visto na tela tem de ser o grafo
  que o cluster rodou.
- **Nada de numérica consome o fluxo compartilhado.** Todo Lanczos do projeto
  — `dae_csr_lanczos_bounds`, a normalização espectral em `dae_hamiltonian`, o
  `λ₂` em `dae_metrics_compute` — usa vetor inicial DETERMINÍSTICO, gerado
  internamente por splitmix64 sobre uma constante fixa
  (`dae_csr_lanczos_start`). Nenhuma delas recebe `dae_rng*`, e isso é imposto
  pela assinatura, não pela disciplina.

  > Duas razões. **Reprodutibilidade:** se o Lanczos puxasse do mesmo fluxo que
  > os geradores de grafo, a ordem das chamadas viraria parte do contrato de
  > forma invisível — bastaria o `.cpp` exportado montar o Hamiltoniano antes
  > (ou depois) da religação para o grafo mudar. **Física:** com
  > `DAE_NORM_SPECTRAL`, a escala de energia passaria a depender de uma
  > estimativa estocástica, e dois grafos comparados com sementes ou
  > `lanczos_steps` diferentes receberiam γ efetivo diferente — contaminando
  > exatamente a comparação a escala fixa que a parte 3.2 existe para
  > garantir.
- **A ordem das chamadas faz parte do contrato**, não só a semente. Todo
  sorteio varre índices em ordem crescente. Nenhum gerador pode sortear
  iterando sobre tabela de dispersão ou `Set` — a ordem mudaria entre alvos e a
  semente deixaria de reproduzir o grafo.
- `spec.json` **mais** `DAE_CORE_HASH` determinam o resultado bit a bit.
  `tools/hash_core.mjs` calcula o hash sobre todos os fontes de
  `core/amalgam.list` (menos o próprio arquivo gerado) e o grava em
  `core/dae_version.generated.h`.
- **Todo arquivo exportado e todo CSV de saída carrega, no cabeçalho, o
  `spec.json` canônico e o `core_hash`.** Um arquivo de resultado solto, achado
  num diretório seis meses depois, continua rastreável até o que o gerou. Cada
  figura de artigo tem de poder voltar ao JSON que a produziu.

Testes: `t90_determinism.c` (mesma semente reproduz a mesma CSR bit a bit;
sementes vizinhas divergem; a semente 0 não colapsa o estado do xoshiro).

## 5. O microtúbulo: costura e pontas

Índice `j = m·N⊥ + n`, com `m ∈ [0, N∥)` ao longo do eixo e `n ∈ [0, N⊥)` entre
protofilamentos, `N⊥ = 13` por padrão.

O acoplamento transversal liga `n → n+1`. Ao dar a volta, o sítio `(m, N⊥−1)`
liga com `(m + seam_shift, 0)`. `seam_shift = 0` é o cilindro periódico ideal;
`seam_shift = 3` é o deslocamento helicoidal realista, que quebra a simetria
translacional transversal e muda o espectro.

**As pontas longitudinais são abertas por padrão** (`longitudinal_closed = 0`),
e isso é física, não conveniência: o microtúbulo real tem extremidades abertas
e **inequivalentes** (plus end e minus end).

> **Por que as duas pontas não são iguais.** Com `seam_shift = s` e pontas
> abertas, as `s` arestas de costura com `m + s ≥ N∥` são descartadas. O
> resultado, medido sítio a sítio em `t92_microtubule.c`, é: **`s` sítios
> deficientes na coluna `0`, na ponta de `m` baixo, e `s` sítios deficientes na
> coluna `N⊥−1`, na ponta de `m` alto.** As duas extremidades perdem ligação,
> mas em **protofilamentos diferentes** — e é essa a inequivalência plus/minus
> end, não uma ponta intacta e outra defeituosa.
>
> A assimetria é consequência da costura helicoidal, não defeito da montagem.
> Quem quiser as pontas equivalentes usa `longitudinal_closed = 1`, e aí perde
> a fisiologia das extremidades.

## 6. Numérica

### 6.1 O limite espectral tem de ser garantido, não estimado

Este é o modo de falha clássico do propagador de Chebyshev, e ele é
**silencioso**. Se `a` ficar menor que a semi-largura verdadeira, `H̃` tem
autovalores fora de `[−1, 1]`, os `T_k` crescem exponencialmente e o resultado
vira lixo sem que nada avise.

O erro é **assimétrico**: subestimar `a` destrói o resultado; superestimar
custa alguns termos a mais em `K`. As regras seguem essa assimetria:

- **Gershgorin é cota rigorosa e é o teto.** O intervalo final nunca a excede.
- **Lanczos só aperta por dentro.** Os valores de Ritz são interiores ao
  espectro (`θ_max ≤ λ_max` sempre), e `β_m = ‖r‖` é a margem rigorosa em
  aritmética exata. Em ponto flutuante a ortogonalidade se perde, então há
  ainda **5% de inflação da semi-largura** e o recorte final contra Gershgorin.
- **Rede de segurança em tempo de execução.** `dae_cheb_step` compara a norma
  antes e depois e devolve `DAE_ERR_NORM` se ela variar mais que `1e-10`.
  Custa duas reduções `O(N)` contra `K` produtos matriz–vetor: ruído no
  orçamento. Vale para qualquer estado de entrada, normalizado ou não, porque
  compara a razão e não o valor absoluto.

Testes: `t04_complete.c` (em `K_50` Gershgorin dá `a = 49` e Lanczos aperta
para `25,6`, ainda contendo o espectro `{−49, +1}`); `t91_conventions.c`
(sabota `a` pela metade e exige `DAE_ERR_NORM`).

### 6.2 Quem decide a ordem é a cauda, não a fórmula

`K ≈ 1.2·(a·Δt) + 20` é **chute inicial**, não teto. A ordem que o propagador
usa cresce até `|2·J_K(α)| < 1e-16`, com `α = a·Δt`.

Usar a fórmula como teto é um defeito silencioso e **não-monotônico em α**, o
que o esconde de qualquer teste feito num único α. A margem que ela deixa acima
do ponto de retorno é `0.2α + 20`; a transição de Airy de `J_k(α)` em `k = α`
tem largura `~α^{1/3}`, de modo que a margem é suficiente para α pequeno,
**insuficiente para α ≈ 20–100**, e volta a sobrar para α grande. Medido, com o
teto valendo:

| α | K do teto | último `|2 J_K|` retido | deriva da norma por passo |
|---|---|---|---|
| 5 | 27 | 9e-16 | −5,5e-17 |
| 20 | 45 | **1,8e-12** | −3,6e-15 |
| 50 | 81 | **1,9e-11** | **+4,2e-13** |
| 120 | 165 | — | +3,9e-14 |
| 300 | 373 | — | ok |

Em `α = 50` isso dá 4,2e-9 em 10⁴ passos: violaria o critério de aceitação por
quatro mil vezes, com uma escolha de `Δt` perfeitamente legítima. Deixando a
cauda mandar, a mesma medida cai para 4,2e-16 por passo.

Com a cauda decidindo a ordem, a deriva obedece

```
|norma − 1| ≈ n_passos · K · 5e-18
```

com `deriva/(passo·K)` constante dentro de um fator 5 ao longo de α ∈ [1, 120]
— o piso de arredondamento do `f64`. Baixar o limiar da cauda abaixo de `1e-16`
não melhora mais nada, e somação compensada de Kahan no acumulador também não
(testada: sem efeito, 38% mais lenta). **Consequência a saber:** o critério
"< 1e-12 em 10⁴ passos" só é alcançável para `K` moderado; com `K ≈ 174` a
deriva em 10⁴ passos é ~2e-11, e não há o que fazer em precisão dupla.

`t01_norm.c` guarda essa lei como **envoltória** varrendo α: é o teste que pega
um truncamento cedo demais.

### 6.3 Bessel por recorrência regressiva

`J_0..J_K` saem por Miller com a normalização `J_0 + 2ΣJ_{2m} = 1`. A
recorrência progressiva é instável para `k > x`: a solução `Y_k`, que cresce,
contamina e destrói a cauda — que é exatamente onde o truncamento decide o
erro.

A rotina é validada contra tabela externa (`t00_bessel_table.c`), gerada no
Wolfram com 30 dígitos e conferida de forma independente contra `mpmath`; as
duas fontes concordam em 2,2e-29. **Isso não é zelo excessivo**: no teste 2,
`J_k(a·dt)` aparece nos coeficientes de Chebyshev *e* `J_j(2γt)` aparece na
solução analítica da linha — a mesma rotina dos dois lados da igualdade. Um
erro de fase ou de normalização se cancelaria parcialmente e passaria.

### 6.4 Grade uniforme é premissa do cache

Os `c_k` dependem só de `dt` e ficam em cache. Uma **grade logarítmica** anula
esse cache e recalcula Bessel a cada passo: funciona, mas custa. Quem for
adicionar "escala log no eixo do tempo" precisa saber disso. Teste:
`t91_conventions.c` verifica que alternar `dt` recalcula os coeficientes.

### 6.5 `λ₂` não é "rodou m passos e pronto"

`λ₂` é reportado contra `p` e vira resultado de manchete. O regime que
interessa — rede fortemente modular — tem `λ₂ → 0`, e é exatamente onde a
convergência é mais lenta. Um número **fixo** de passos devolve `λ₂`
**superestimado** ali, achatando a curva na região mais interessante; e um
platô numérico é indistinguível de um platô físico quando se olha só a curva.

Por isso `dae_metrics_compute` tem critério de convergência, teto de passos,
**resíduo devolvido** e bandeira. Se `lambda2_converged` for 0, o valor é um
**limite superior** (quociente de Rayleigh de um vetor em `1⊥`), não uma
medida, e a interface tem de dizer isso. `lambda2_residual = ‖Lx − λx‖` é cota
rigorosa: existe autovalor de `L` a menos disso.

Quanto isso importa, medido em SBM `n = 200`, `M = 4`, `p_in = 0.35`,
`p_out = 0.002`: com o critério de convergência, `λ₂ = 0,178453`, batendo com
Jacobi denso em `1e-8`; com **2 passos fixos**, `λ₂ = 10,36` — 58× maior.

A deflação do vetor constante aparece em quatro pontos de `dae_fiedler`, e três
deles são **mutuamente redundantes de propósito**. Está escrito no cabeçalho da
função, junto com a medida que motivou cada um; não os "simplifique" achando
que são a mesma coisa repetida.

## 7. Observáveis

### 7.1 Concurrence: a convenção, e o termo espúrio do bloco diagonal

Para estado puro, `C_ij = 2|ψ_i||ψ_j|`. Com `s_M = Σ_{i∈M}|ψ_i|` e
`q_M = Σ_{i∈M}|ψ_i|²`, somando sobre **pares não-ordenados** `{i,j}`, `i ≠ j`:

```
C_MN = 2 s_M s_N        (M ≠ N)
C_MM = s_M² − q_M       (o produto externo incluiria i == j, que não é par;
                         o termo espúrio é exatamente q_M)
```

A matriz é armazenada simétrica, e **a soma do triângulo superior *com* a
diagonal reproduz `C_ℓ1` = `(Σ|ψ_j|)² − 1`**. Essa identidade é o teste da
convenção (`t91_conventions.c`), e é ela que amarra a concurrence agregada à
coerência `ℓ1`.

Custo: `O(N + M²)`, não `O(N²)` — a versão agregada por módulo é um produto
externo. Só a matriz completa é cara, e ela é oferecida apenas para
`N ≤ DAE_FULL_CONC_MAX` (2000).

### 7.2 Sem alvo, `p_alvo` é NaN

`target = −1` produz `p_alvo = NaN`, **nunca zero**. Zero é um valor
fisicamente válido — o alvo pode simplesmente estar vazio — e acabaria plotado
como se fosse dado.

## 8. Montagem da CSR: duplicata é descartada

`dae_csr_from_edges` **descarta** aresta repetida; não soma os pesos. Com
`seam_shift` e `N⊥` pequeno dá para gerar duplicata legítima, e somar dobraria
`j_perp` em alguns sítios sem nenhum aviso. O número de entradas descartadas
sai em `n_dropped`, para a interface poder avisar.

Invariante: dentro de cada linha, `colind` é estritamente crescente.

## 9. Layout de memória

- Matriz esparsa em **CSR**; grau típico ~4, `nnz ≈ 4N`.
- `H` é **real simétrica** (adjacência ou laplaciana com pesos reais), então o
  SpMV aplica a mesma matriz a `re[]` e a `im[]` — metade do tráfego de
  memória.
- Complexos como **dois arrays reais separados** (`re[]`, `im[]`). Nunca
  `double complex` intercalado: o SIMD do WASM só rende no layout SoA.
- `f64` na propagação; `f32` apenas no buffer que alimenta a textura WebGL.
- Índices em `int32_t`.

### 9.1 A fronteira WASM: nunca guarde uma view

Todo buffer que o núcleo devolve ao JavaScript vive no heap do WASM. **Qualquer
`malloc` do lado C que faça o heap crescer troca o `ArrayBuffer` inteiro**, e
toda view que o JavaScript tenha guardado passa a apontar para memória
*detached* — leitura silenciosa de zeros, sem exceção, sem aviso.

O modo de falha é o pior possível: funciona em toda rede pequena, e quebra
quando o usuário aumenta `N`.

A regra, em `wasm/daedalus.mjs`: as views saem sempre de `wasmMemory.buffer`,
que é o buffer corrente por definição, e são construídas **na hora de cada
leitura**. Construir uma view é criar um objeto, não copiar dados — continua
sendo leitura sem cópia, que é o ponto do contrato. Existem exatamente três
fábricas de view no projeto (`f64_`, `f32_`, `i32_`); nada mais constrói uma.

O que atravessa `postMessage` do worker para a interface é **cópia** (`slice()`),
porque uma view do heap do WASM não sobrevive nem ao `postMessage` nem ao
crescimento. O caminho sem cópia é o de dentro do worker, do WASM para a
textura WebGL.

`wasm/teste_memoria.mjs` verifica as **duas** metades: que a view guardada de
fato morre (senão o contrato seria superstição) e que a refeita continua certa.

## 10. A interface: erro deixa de ser número errado

No núcleo, todo erro aparece como número errado, e há oráculo analítico para
comparar. Na interface, erro aparece como **imagem plausível**:
um heatmap com os eixos `m`/`n` trocados, ou com o colormap invertido, é
indistinguível de física interessante para quem está olhando — e vira figura de
artigo.

Três defesas, e nenhuma delas é "olhar e achar bonito".

### 10.1 O mapeamento e a escala de cor moram sozinhos, e são testados

A transposição sítio→texel está em `src/nucleo/indices.ts`, sozinha e
exportada. A escala de cor está em `src/nucleo/paleta.ts`, e o shader monta a
LUT a partir dela em vez de ter uma cópia. Os testes atacam essas duas funções
diretamente, não o pixel na tela: `paleta.test.ts` exige monotonicidade em
luminosidade (senão "mais claro = mais população" vira convenção, não verdade)
e `indices.test.ts` exige que um delta caia num texel específico — com
`N∥ ≠ N⊥`, uma transposição muda esse índice na hora.

### 10.1.5 Toda invariância precisa de uma companheira

Regra geral, e não um caso particular das fixtures: **toda asserção de que algo
se conserva precisa de uma companheira que prove que aquilo pode ser
quebrado.** A pergunta a fazer de cada teste novo é *existe um estado
degenerado em que isto passa sem testar nada?* — e quase sempre existe:

| asserção | passa trivialmente quando |
|---|---|
| o padrão é simétrico | o estado não evoluiu (um delta é simétrico) |
| a norma se conserva | nada evoluiu |
| `p_alvo` bate nos dois lados | é zero dos dois lados |
| os dois CSV são idênticos | os dois estão vazios ou truncados |
| a varredura cobre o parâmetro | todos os pontos caíram no mesmo valor |
| o leitor aceitou o arquivo | o leitor aceita qualquer coisa |

As companheiras estão no código: `seam_shift = 3` tem de QUEBRAR a simetria
transversal; o comparador de CSV exige norma em 1 **e** pacote espalhado; a
varredura de α falha se a faixa coberta não passar de duas ordens de grandeza;
cada teste de aceitação do reimportador tem um gêmeo de recusa.

### 10.1.6 Ausência não se afirma sem provar que a sonda enxerga

Regra de **método**, e não só de código: vale para busca em texto, para varredura
bibliográfica, para consulta a banco, e para qualquer lugar onde "não achei" vire
conclusão.

**Nenhuma afirmação de ausência sem a medida que prova que a sonda enxerga.**

Uma sonda quebrada e uma ausência real produzem exatamente o mesmo resultado — vazio — e
a diferença entre as duas é toda a diferença entre um achado e um engano. Um fundo de
canvas da mesma cor que `cor(0)` é um exemplo: o mapa pode estar correto, a rede estar lá,
e a medida não a enxergar.

Na prática, duas exigências:

1. **Consulta-canário.** Antes de confiar num "não encontrei", rode uma consulta cuja
   resposta você já conhece e exija que ela apareça. Se o canário não volta, a busca está
   quebrada e todo "não encontrei" dela é ruído.
2. **A sonda é conferida no próprio alvo.** Contar zero ocorrências num PDF só vale depois
   de conferir que o PDF produziu texto: 6 páginas e 23 939 caracteres extraídos é medida;
   `grep -c` devolvendo zero num arquivo que não extraiu nada é o mesmo zero, e não quer
   dizer nada.

Registre quantos resultados cada consulta devolveu. **Consulta com zero resultados é
suspeita, não conclusão** — e uma afirmação de "ninguém fez X" precisa listar as consultas
que procuraram X e falharam, senão a ausência de evidência não mostra onde procurou.

**Corolário sobre COBERTURA.** Instrumentação pedida por cobertura, e não por investigação,
revela o que ela não foi construída para procurar — cinco pontos em γ em vez de dois, `p`
alto com o grafo variando. O achado tende a aparecer no extremo da faixa: o colapso da
extrapolação de Richardson só aparece onde o viés é grande e o valor verdadeiro é pequeno.

A lição prática é mais estreita do que "meça bastante": **cobertura não serve só para validar
a faixa, serve para descobrir onde a faixa acaba.** Um portão medido num ponto só diz que ele
funciona naquele ponto; medido nos extremos, ele diz onde para de funcionar — e é essa
segunda informação que impede que a extrapolação seja feita por otimismo.

### 10.1.7 Toda comparação numérica declara a sua RESOLUÇÃO

Uma tolerância diz o que a comparação aceita. Ela não diz o que a comparação
**consegue ver** — e são perguntas diferentes, com respostas que podem estar em
lados opostos.

**Para toda comparação numérica, meça o menor desvio que ela detecta, e registre
esse número ao lado da tolerância.** Um teste cujo limiar de sensibilidade nunca
foi medido pode estar passando por não conseguir enxergar, e ninguém saberia: um
verde por concordância e um verde por cegueira são o mesmo verde.

Exemplo em `t97_modos_saida`: a sabotagem de 1 ULP na
maior amplitude **não** quebra a igualdade entre os dois modos de saída, e isso
não é defeito da comparação — é a resolução dela: numa média sobre 24
trajetórias, uma perturbação relativa de ~1e-16 numa contribuição vira ~4e-18 no
resultado, abaixo do ULP do próprio acumulado. O teste declara isso, e usa
sabotagens acima desse piso (precisão de f32, acervo truncado).

Corolário, e é a outra metade da mesma lição: **a sabotagem precisa atingir uma
quantidade que participa do resultado.** Perturbar a amostra `t = 0`, onde
ψ = δ e quase toda amplitude vale exatamente zero, não testa nada:
`nextafter(0)` dá um denormal de 5e-324 que desaparece na primeira
multiplicação. "Alterei um valor" não é "alterei um valor que importa", e um
teste que ataca o lugar errado é tão inerte quanto um que não ataca.

### 10.1.8 Antes de culpar o gerador, teste um estimador LINEAR

Quando o desvio de uma média sobre amostras não cai como `1/√n`, a conclusão
óbvia é correlação entre as amostras — e o suspeito óbvio é o gerador de números
aleatórios. A conclusão óbvia pode estar errada, e o custo de investigá-la é
alto: o PRNG é a peça mais frágil e mais cara de auditar do projeto.

**Rode o mesmo diagnóstico sobre um estimador LINEAR nas mesmas amostras.**

Um funcional linear da média amostral é não enviesado e obedece `1/√n` por
construção. Se ele obedecer e o outro não, as amostras estão bem e o problema é
a **não linearidade** do funcional. Se nem ele obedecer, aí sim há correlação
entre amostras, e nada mais importa até isso ser resolvido.

Exemplo: `C_inter = Σ|ρ_ij|` dá `sd ~ n^(-0,16)` contra os `n^(-0,5)`
esperados. O controle linear — população de um módulo, `Σ_{i∈M} ρ_ii` — dá
`n^(-0,56)` no mesmo conjunto de trajetórias. Diagnóstico em uma execução, sem
tocar no gerador.

O mecanismo é conhecido, com nome, em pelo menos três literaturas: `|·|` é
convexo, então por Jensen `E[Σ|ρ̂_ij|] ≥ Σ|E ρ̂_ij|` — **viés positivo
garantido**, que cai como `1/√n`, a mesma taxa do ruído.

### 10.2 Fixtures com resposta conhecida, e anti-vacuidade junto

Com `seam_shift = 0` a rede é o produto cartesiano `P_41 × C_13`: a dinâmica
fatoriza e o padrão é simétrico nos dois eixos, com `Var(m) = 2(γt)²` exato da
caminhada livre e centroide longitudinal igual em **todos** os anéis.

Isso é previsão analítica, não instantâneo medido — mas sozinho seria **vácuo**:
um estado parado também é simétrico. Por isso todo teste de simetria vem
acompanhado de (a) uma verificação explícita de que o pacote se espalhou e (b)
o caso `seam_shift = 3`, que tem de QUEBRAR a simetria transversal.

> O modo de falha é concreto. Um wrapper JavaScript que montasse o vetor de
> parâmetros com zeros, em vez de partir dos padrões do núcleo, zeraria
> `j_par` e `j_perp`: o grafo sairia com o número certo de arestas, todas de
> peso nulo, `H` seria a matriz nula e o estado não se moveria. **Na tela, isso
> não parece bug — parece localização.** Quem pega é o caso de anti-vacuidade.
>
> A defesa é a mesma regra da parte 2: quem sabe os padrões é o C
> (`dae_ws_params_default`), e o JavaScript parte deles.

Sobre a assinatura congelada (`src/nucleo/assinatura.ts`): ela precisa ser
**conjunta**, não marginal. Duas escolhas marginais são cegas — um único número
que seja o primeiro momento (ou seja, o centroide, cego a qualquer rearranjo
simétrico), e estatística circular no eixo `q`, correta e ainda assim cega, porque a costura ACOPLA `m` e `q` e marginalizar sobre `m` destrói
exatamente a informação que ela cria.

### 10.2.5 Layout: desenhar a rede, não a rotulagem

Grafo sem geometria própria — SBM, hipercubo, `K_N`, lista importada — não tem
coordenada nenhuma que o gerador saiba dar. Desenhar os sítios **em ordem de
índice**, dobrados em linhas, desenha a ORDEM DE ROTULAGEM, não a rede: para um SBM a numeração é arbitrária, e a figura
resultante não diz nada sobre topologia enquanto parece que diz.

O recurso é o **layout espectral**: autovetores 2 e 3 da laplaciana como
coordenadas `(x, y)`. Ele reaproveita o Lanczos deflacionado do `λ₂` — é
literalmente o vetor de Fiedler que já se calcula, mais o próximo — e, para
rede modular, separa os módulos no plano, que é o que se quer enxergar.

`dae_graph.geom` diz o que o gerador soube dar: rede desenrolada (microtúbulo),
geometria própria (linha, ciclo, grade) ou nenhuma. O padrão da vista sai daí,
e **o nome do layout em uso fica visível na interface** — a figura muda de
significado conforme ele, e adivinhar qual está ativo é o começo de ler errado.

A fixture (`t95_layout.c`) mede a razão entre a distância média ENTRE módulos e
DENTRO de um módulo. Com `p_in = 0.35`, `p_out = 0.004`: **20,5**. A companheira
usa `p_in = p_out`, onde não existe módulo — a partição continua rotulando os
vértices, mas não corresponde a estrutura nenhuma — e ali a razão tem de ser
**0,996**. Sem esse par, a asserção de separação estaria medindo a numeração
dos vértices, não a topologia.

Normalização para a caixa do canvas: margem fixa em pixels e escala
**uniforme** nos dois eixos. Esticar `x` e `y` de forma independente encheria a
caixa por completo, mas faria a distância no plano mentir — e o ponto do
embedding é justamente que distância ali significa alguma coisa.

### 10.3 A norma fica na tela, sempre

O diagnóstico de norma não mora num painel de depuração: fica no rodapé,
permanente, com `|1 − norma|` em notação científica e o rodapé inteiro em
vermelho quando passa de `1e-9`. É o único indicador que o usuário — você,
daqui a um ano, ajustando `dt` num regime novo — pode olhar para saber se o que
está vendo é confiável. Junto dele ficam `λ₂` com o aviso de não-convergência,
`Q`, `|E|`, componentes, arestas duplicadas descartadas e religações sem
destino livre.

### 10.4 Teste de fumaça mede pixel, não intenção

`npm run fumaca` abre o aplicativo num Chrome sem cabeça pelo protocolo
DevTools, espera a propagação terminar e lê os pixels de volta do canvas. Sem
isso, um mapa preto com séries corretas passa despercebido.

O alcance está escrito no próprio arquivo: ele garante que a figura final
existe, não a animação quadro a quadro. Quem garante o conteúdo são as
fixtures.

> Nota de ferramenta: `--dump-dom` e `--virtual-time-budget` do Chrome sem
> cabeça **não esperam trabalho assíncrono em Web Worker**. Eles devolvem
> silêncio, e silêncio é indistinguível de falha. `tools/navegador.mjs` fala o
> protocolo DevTools direto, com o WebSocket embutido do Node e zero
> dependências, e espera em tempo real.

## 11. Código

Comentários em **português**, identificadores em **inglês**. Prefixo `dae_` em
tudo que é público.

Comentário explica **por que**, não o que. As regras deste documento aparecem
citadas no ponto do código em que valem.

### 11.1 Exceção: o pacote Wolfram tem símbolos públicos em português

`templates/wolfram/Daedalus.wl` exporta `DaedalusRede`, `DaedalusPropagar`,
`DaedalusObservaveis`, `DaedalusVarredura`, `DaedalusVerificarReferencia` e o
resto em **português**, CamelCase com prefixo `Daedalus`. É a única exceção à
regra de identificadores em inglês, e a razão fica registrada porque uma
exceção sem razão escrita vira precedente.

O núcleo em C é lido por quem trabalha no motor, e ali o inglês é a convenção do
campo. Este pacote é lido e **modificado** por quem trabalha em Mathematica —
ele não é um oráculo que se roda e se confere, é código que se altera. Os
comentários que explicam a física estão em português; ter os símbolos em inglês
no meio deles produziria um texto que muda de língua a cada linha, e a fricção
recairia justamente sobre quem mais precisa entender o que está lendo.

O prefixo `Daedalus` não é decorativo: Wolfram tem um único espaço de nomes
global por contexto, e um símbolo curto colide com o caderno do usuário.

A regra do resto do projeto **não muda**. Esta exceção vale para
`templates/wolfram/` e só para ele.

### 11.1.1 Por que o pacote Wolfram não é uma tradução do núcleo

Entre `.cpp` exportado, WASM e binário nativo a concordância é **estrutural**:
mesmo texto, três compiladores, não há o que divergir (parte 1). Com o Wolfram
isso é impossível — Wolfram não inclui C, e uma biblioteca compilada por
LibraryLink seria caixa-preta, o oposto de inspeção.

Então a garantia muda de natureza e passa a ser **empírica**: dois métodos,
mesmo resultado. Para isso valer alguma coisa o método tem de ser
**deliberadamente diferente**. O núcleo propaga por expansão de Chebyshev; o
pacote propaga por decomposição espectral (`Eigensystem`) e por Krylov
(`MatrixExp[-I H t, psi]`) quando N cresce. **Não há Chebyshev no pacote, de
propósito**: duas implementações do mesmo algoritmo podem compartilhar o mesmo
erro conceitual, e a concordância entre elas não provaria nada.

O que TEM de bater exatamente, sem tolerância, é o que não é numérico: o fluxo
do PRNG, a impressão digital do grafo e o contador de religações que falharam.
Ver parte 12.1.

## 11.2 Reimportação: a única entrada que não nasce aqui

O CSV que volta do cluster é a única entrada do sistema que não passa pelo
parser estrito por construção — e por isso é a que mais precisa desconfiar. Um
arquivo de uma versão anterior do núcleo, ou com as colunas em outra ordem,
plotaria e plotaria **errado**, sem nada quebrar.

Três defesas, em `src/nucleo/reimportar.ts`:

1. **Procedência obrigatória.** Sem `#! spec` e `#! core_hash` o arquivo é
   recusado. Não dá para saber que simulação ele descreve, e um gráfico sem
   essa resposta é pior que nenhum gráfico.
2. **O spec passa pelo parser estrito em C**, via WASM (`dae_ws_valida`), o
   mesmo do resto do sistema. Se não passa, o CSV não entra — e a mensagem que
   aparece na tela é a do parser, com linha e coluna.
3. **Colunas por nome, nunca por posição.** Ordem diferente tem de funcionar;
   coluna ausente tem de aparecer como ausente, não como outra coluna.

`core_hash` diferente é **aviso, não bloqueio**: reler resultado antigo é
legítimo. Mas o aviso vai para a **tela**, não para o console, e o seletor de
modo vira bronze — o número deixou de ser calculado aqui.

**`#! implementacao` diz QUEM calculou**, que é outra pergunta e o `core_hash`
não responde: o mesmo core_hash pode ter produzido o número pelo núcleo em C
(`c`) ou pelo pacote Wolfram (`wolfram`), que usa outro método. A interface
exibe a origem no selo, ao lado de "Reimportado". Ausência do campo — CSV de
versão anterior — aparece como **"origem não declarada"**, nunca como "núcleo":
supor a origem seria atribuir ao arquivo uma procedência que ele não afirma.

## 11.3 Internacionalização

Quatro línguas: cada entrada do catálogo é uma
**função**, não um texto com marcadores, para que cada língua resolva a própria
gramática sem biblioteca de pluralização no meio. Chave ausente devolve a
própria chave — o buraco aparece na tela em vez de virar string vazia.

`src/i18n/catalog.test.ts` assere que as quatro línguas têm **exatamente** o
mesmo conjunto de chaves, que toda chave usada na interface existe no catálogo,
e que nenhuma tradução devolve string vazia. São dez linhas e pegam o modo de
falha inteiro na CI. O teste que varre o código em busca de `t('chave')` tem
sua própria anti-vacuidade: falha se encontrar menos de dez chaves, porque um
varredor quebrado acharia zero e passaria.

## 11.5 Identidade visual

Uma das regras de aparência é **funcional, não decorativa**:
**Azul Egeu** identifica o laboratório interativo — números calculados
localmente, no navegador, agora. **Bronze** identifica tudo que passou pelo
ciclo de exportação e reimportação. Nenhum valor reimportado aparece em azul, e
a distinção não se reaproveita para outra coisa.

O mapa de cor de probabilidade tem luminância monotônica para sobreviver à
impressão em escala de cinza; o de fase é cíclico, porque 0 e 2π são a mesma
coisa e um mapa sequencial inventaria uma descontinuidade que a física não tem.
Os dois vivem em `src/nucleo/paleta.ts`, e figuras feitas fora da interface
(matplotlib, Wolfram) usam os mesmos valores — a figura da tela e a figura do
artigo têm de ter a mesma leitura. Todo valor numérico da
interface vai em IBM Plex Mono.

Nenhum elemento de microtúbulo, de biologia ou de Orch-OR entra no logotipo ou
na identidade: o motor é geral para qualquer grafo, e a identidade acompanha
isso. O microtúbulo é um cartão da galeria de geradores, com o mesmo peso
visual dos demais.

## 12. Portões de qualidade

```
make -C native test        # todos os testes de aceitação
make -C native cxx-check   # o núcleo amalgamado compila limpo em C99 e C++17
make -C native asan        # ASan + UBSan
make -C native tsan        # ThreadSanitizer (setarch -R por causa do ASLR)
make -C native bench       # metas de desempenho
make -C native mutants     # arnês de mutação (manual, não entra na CI)

source ~/emsdk/emsdk_env.sh
make -C native wasm        # módulo + artefatos de linha de comando
make -C native test6       # ACEITAÇÃO 6: WASM contra nativo
npm run teste7             # ACEITAÇÃO 7: navegador × .cpp exportado, varrendo α
make -C native wasm-test   # a suíte INTEIRA sob WASM + fronteira de memória
make -C native wasm-bench  # desempenho no alvo real

npm test                   # fixtures da interface (mapeamento, paleta)
npm run fumaca             # o aplicativo abre, propaga e DESENHA
npm run build              # bundle estático
```

### 11.6 A ordem da soma é contrato, não consequência do laço

Vale para toda redução sobre um ensemble, como a média de trajetórias.

Soma de ponto flutuante não é associativa. Se a ordem em que as contribuições
entram depender do laço — ou pior, do escalonador — dois caminhos que deveriam
concordar bit a bit passam a concordar por sorte, e a divergência aparece quando
alguém muda `OMP_NUM_THREADS`.

`dae_rho_acc_somar` **recusa** índice fora de sequência com `DAE_ERR_PARAM`. Não
é comentário pedindo boa-fé: quem somar fora de ordem recebe erro em vez de
receber outro número.

**Por isso a paralelização é por CÉLULA da grade, não por trajetória dentro de
uma célula.** Acumular por thread e reduzir em ordem de thread parece resolver e
não resolve: o número de threads entra no resultado, e o mesmo `spec.json` com
outra máquina dá outros dígitos. Com uma célula por thread, cada acumulador é
privado, a soma interna é serial e ordenada, e a reprodutibilidade não depende
de quem terminou primeiro.

E a semente é **derivada do índice**, nunca puxada de um gerador compartilhado
(`dae_traj_semente`): um fluxo comum faria o consumo depender da ordem de
chegada das threads, que é a mesma doença por outra porta.

### 12.1 Verificação do pacote Wolfram: manual, e por quê

A CI **não tem Mathematica**. Este portão é manual, e por isso cada execução
deve ser registrada com **data, `core_hash` e resultado**. Verificação manual
sem registro é verificação que ninguém sabe se aconteceu.

```
wolframscript -e 'Get["templates/wolfram/Daedalus.wl"];
  DaedalusVerificarReferencia["specs/oraculo"]'
wolframscript -e 'Get["templates/wolfram/Daedalus.wl"];
  DaedalusAntiVacuidade["specs/oraculo"]'
```

**A ordem da conferência é o diagnóstico**, e é ela que separa causas que se
corrigem em lugares diferentes:

1. **PRNG contra tabela de referência** (`specs/oraculo/prng.json`, gerada por
   `native/tests/t96_prng_tabela.c`). Barata e a mais diagnóstica: se o
   xoshiro256++ divergir, redes estocásticas saem diferentes e a discordância
   não aparece como erro numérico, aparece como **física diferente**. Mesmo
   raciocínio do teste 0 contra a tabela de Bessel.
2. **Impressão digital do grafo**, exigida **idêntica**, sem tolerância nenhuma.
3. **Contador de religações que falharam** (`#! rewire_failed`). A digital não
   vê isto: uma tentativa que colide consome 100 sorteios e não muda aresta
   nenhuma. O caso `completo-religado` existe só para isso — em K₁₂ **toda**
   tentativa colide, o grafo sai idêntico ao original, e esse contador é a única
   testemunha de que os dois gastaram o mesmo fluxo. Um caso com
   `rewire_failed = 0` compararia zero com zero e não provaria nada.
4. **Observáveis**, com tolerância declarada.

**Três medidas, e nenhuma sozinha.** Erro relativo puro não serve: em t = 0 o
núcleo escreve `p_alvo = 0` exato (Chebyshev em α = 0 dá a identidade sem
arredondamento) e a soma espectral devolve 3e-32; o erro relativo é 1, o pior
possível, para uma concordância de trinta e duas casas. Há uma lei por trás: se
as amplitudes diferem de δ em valor absoluto, p = |ψ|² difere relativamente de
~2δ/√p, o que explica cinco ordens de grandeza de erro relativo com uma única
concordância de ~1e-15 na amplitude. Então:

| medida | pior observado | declarada | margem |
|---|---|---|---|
| absoluta (tabela de escalares) | 9.4e-12 — `linha`, `coh_l1` que é O(n) | 1e-10 | 10× |
| relativa (só acima do piso 1e-10) | 9.6e-11 — `microtubulo-seam0` | 1e-9 | 10× |
| amplitude (estado final) | 1.5e-14 — **grade 2D**, n = 240 | 1e-12 | 68× |

O piso **não afrouxa** o teste: abaixo dele o critério absoluto continua valendo
e é ele que reprova discordância real. E os dois critérios se dividem o trabalho
sem sobra — `norm`, `ipr` e `coh_l1` estão sempre muito acima do piso e caem no
relativo; `p_alvo` e as populações por módulo atravessam caudas de 1e-30 e caem
no absoluto.

A grade 2D é mesmo o pior em amplitude, como se espera de espectro muito
degenerado. O piso de ruído da soma espectral escala com n (≈ nε/4): quem rodar
N ≫ 300 deve esperar a tolerância de amplitude subir junto — e reajustar sabendo
por quê, não para fazer passar.

**A escala é imposta.** O núcleo normaliza por estimativa de Lanczos com margem,
o pacote pelo raio espectral exato, e os dois discordam por construção (em
`microtubulo-espectral`, 5.799 contra 3.975). A verificação impõe a escala do
núcleo para comparar o **propagador**, e reporta a diferença das duas escalas
como número — senão uma discordância de norma contaminaria todas as colunas e
esconderia qualquer discordância de propagador atrás dela.

**Anti-vacuidade obrigatória.** `DaedalusAntiVacuidade` sabota o grafo
(`seam_shift` 3→2), o fluxo de sorteios (`rewire_failed` 30→29 com a digital
intacta), os números (deslocamento de 1e-6 **e** de 1e-15, para mostrar que o
limiar tem dois lados), o arquivo (truncado, e sem o bloco de estado), e confere
a identidade Σ_{M≤N} C_MN = C_ℓ1, que é álgebra da convenção e não comparação
com o núcleo. Cada sabotagem tem de produzir **o veredito correspondente**: sem
isso, "os 12 casos passaram" e "a verificação está quebrada e sempre verde" são
indistinguíveis pelo resultado.

**`make -C native mutants`** aplica um defeito deliberado por vez e imprime
quem mordeu. Comentário não roda na CI; este roda. Onde o projeto tem defesa
redundante de propósito, o mutante fica lá marcado `esperado: CEGO` — a
cegueira aparece na tabela em vez de ficar escondida num cabeçalho.

**`make -C native wasm-test`** é o portão que pega `-ffast-math` herdado por
flag: `t91_conventions` exige que `p_alvo` seja NaN sem alvo, e `-ffast-math`
transformaria isso em zero silenciosamente. Rodar a suíte inteira no alvo real
é mais barato e mais confiável do que auditar flags.

Compilação com `-Wall -Wextra -Werror -pedantic`. O `tsan` vigia a regra de
"zero estado global mutável", da qual depende o laço OpenMP sobre realizações
do template exportado.

---

## 13. Setor magnético: conversão de helicidade e escape

Vale para todo código que propagar `ν_L ↔ ν_R` (Dirac) ou `ν ↔ ν̄` (Majorana) pelo momento magnético
do neutrino. As regras das partes 1–12 continuam valendo; esta parte só
acrescenta o que a física nova impõe.

### 13.1 O piso de factibilidade lê o meio pelo mesmo caminho que a propagação

`escala_em` (o raio espectral que fixa quantos subpassos uma fatia exige no
mínimo) e `passo_fatia` (que propaga) **têm de avaliar `ρ` e `Y_e` pela mesma
função**, `dae_nu_densidade` e `dae_nu_ye` no raio. Se o piso usasse o `Y_e`
escalar do perfil enquanto a propagação usa `ye_table`, com o escalar abaixo
da tabela o Chebyshev recusaria fatias que o método resolve
(`DAE_ERR_TOOBIG`). Teste: `t105_piso_ye_tabela`.

Regra geral, e é por ela que esta seção existe: toda estimativa usada para
dimensionar um passo (raio espectral, largura de janela, localização de
cruzamento) consulta o meio **pela mesma função** que a propagação, nunca por um
atalho escalar. Quando o meio ganhar campo vetorial, polarizações ou fundo de
neutrinos, a regra se estende a eles.

### 13.2 Base, unidades e referencial

- Base Dirac: `Ψ = (ν_eL, ν_μL, ν_τL, ν_eR, ν_μR, ν_τR)`. Base Majorana:
  `Ψ = (ν_e, ν_μ, ν_τ, ν̄_e, ν̄_μ, ν̄_τ)`. O núcleo implementa só o par
  `(ν_eL, ν_eR)`; Majorana e mais de um sabor devolvem `DAE_ERR_UNIMPLEMENTED`.
- No spec: `μ_ν` em `μ_B`, `B` em gauss, distâncias em km, como no resto do
  GHOST. `H` em eV; `t = L[km] · DAE_NU_T_POR_KM`.
- `B_⊥`, `φ_B` e o ângulo `ϑ` entre `k̂` e `B` são medidos no **referencial
  ortonormal local do observador estático**, não em coordenadas.
- **O comprimento `s`** (coluna `s_km` do bloco `trajectory`) é o comprimento
  espacial **local medido pelo observador estático** ao longo do raio nulo,
  `dl² = e^λ dr² + r² dΩ²`. Não é o intervalo próprio da geodésica — que é
  **zero** para um raio nulo —, nem o parâmetro afim, nem a coordenada `r`. É
  nesse comprimento que `H` (em eV, no referencial do observador estático) gera
  a evolução: `i dψ/dl = H ψ`, `t = l · DAE_NU_T_POR_KM`.
- **A base transversal é TRANSPORTADA ao longo da trajetória**, e `φ_B` é medido
  nela. Uma base reconstruída em cada ponto por uma regra local (`ẑ × k̂`, ou a
  base esférica local) gira em torno de `k̂` quando o raio se curva, e essa
  rotação entra na evolução exatamente como uma rotação do campo: um termo
  `dφ_B/ds` em `Δ_ef` que não existe — **ressonância artificial**. Numa métrica
  estática esférica o raio fica num plano e a base transportada é
  `e1 = normal à órbita` (constante), `e2 = k̂ × e1` (sem rotação de Faraday
  gravitacional; ela aparece em Kerr, e aí a rotação da base transportada
  também entra na evolução e tem de vir da geometria). Medido: a base
  `ẑ × k̂` gira até 0,072 rad/km perto de
  3,5 M, um `Δ_ef` falso de 1,4×10⁻¹¹ eV — desprezível diante de
  `Ω = 5,8×10⁻⁶ eV` (`μ = 10⁻¹² μ_B`, `B = 10¹⁵ G`), mas **7,7 vezes** o `Ω` do
  momento mínimo de Dirac. Por isso **o núcleo não constrói base**: recebe
  `b1`, `b2` já projetados, e o spec declara
  `"transverse_basis": "parallel_transported"` (obrigatório).
  `dae_numag_base_transversa` (base fixa) só vale para raio reto.
  **O transporte é verificado pela derivada covariante, não só pela simetria:**
  calcula-se `De/dl = de/dl + Γ(p, e)` ao longo do raio e exige-se `De1 → 0` e
  `ω = g(e1, De2) → 0` em 2ª ordem, fora e dentro da estrela. Uma métrica nova
  (por exemplo, TOV tabelada) entra nesse teste antes de produzir qualquer
  trajetória.
- Tabelas ao longo de `s` interpolam as **componentes** `b1`, `b2`, não `φ_B`:
  evita o salto de 2π. Nos nós, o valor tabelado vale exatamente.
- Acoplamento: `M_B = μ_ν B_⊥ e^{iφ_B}`, com `φ_B` em radianos, medido na base
  transversa fixa de `dae_numag_base_transversa` (`e1 = ẑ × k̂` normalizado, ou
  `x̂ × k̂` perto de `ẑ`; `e2 = k̂ × e1`). Quem monta um campo que gira usa a mesma
  base. Sempre que
  `φ_B` varia ao longo da trajetória, `H` é complexo e exige o mergulho
  `2n × 2n` (seção *How it works* do `README`).

### 13.3 O desajuste efetivo e as suas unidades

Com campo transversal que gira, a ressonância é o zero do desajuste
**efetivo**, não de `Δ`:

```
Δ_ef [eV] = Δ [eV] + (dφ_B/ds)[rad/km] / DAE_NU_T_POR_KM
```

O sinal **não é escolhido e não é universal**: é derivado da convenção de
`M_B` acima e **só vale com ela** — com `H_LR = Ω e^{−iφ_B}` seria
`Δ − (dφ_B/ds)/T`, e com `φ_B` medido na orientação oposta (`e1 × e2 = −k̂`)
também troca. Com
`H_LR = Ω e^{+iφ_B}`, a passagem ao referencial que gira com `B_⊥`,
`ψ = diag(e^{iφ/2}, e^{−iφ/2}) ψ̃`, soma `−iU†dU/ds = diag(φ'/2, −φ'/2)` a `H`.
**Confirmado** por `t107` (bloco 7): ressonância completa em
`Δ = −(dφ_B/ds)/DAE_NU_T_POR_KM`, onde `Δ ≠ 0`, e amplitude 0,2 no sinal oposto.
O mutante com `H_LR = Ω e^{−iφ_B}` falha exatamente essas quatro verificações e
nenhuma outra — as demais são cegas ao sinal, como diz 13.7. Quem mudar a
convenção de `M_B` muda o sinal aqui e tem de ver `t107` falhar antes.

Toda grandeza com gradiente usa unidades naturais explícitas:
`Δ'_nat [eV²] = Δ' [eV/km] / DAE_NU_T_POR_KM`, `γ = 2πΩ²/|Δ'_nat|`,
`ℓ_LZ [km] = |Δ'_nat|^{-1/2} / DAE_NU_T_POR_KM` e
`N_col = ∫ Γ(s) · DAE_NU_T_POR_KM ds` (Γ em eV, `s` em km).

### 13.3b O fundo de neutrinos no desajuste

Com `ν_e` aprisionados, o potencial do `ν_eL` ganha o termo de Nötzold e
Raffelt (1988):

```
Δ = √2 G_F n_b [(3Y_e − 1)/2 + 2Y_νe],   Y_νe = (n_νe − n_ν̄e)/n_b
```

e a ressonância (com `φ_B` constante) passa de `Y_e = 1/3` para
`Y_e + (4/3)Y_νe = 1/3`. Entra pelo spec como `trajectory.ynu_table`,
opcional; sem ele, `Δ` e o canônico são bit a bit os do caso sem fundo (`t112`).
Regras:

- `Y_νe` pode ser negativo (excesso de antineutrinos); `|Y_νe| ≤ 1`.
- Só `ν_e`: assimetrias de `ν_μ`, `ν_τ` entrariam com peso 1 e são tomadas
  nulas (μ = 0). Não há termo térmico (de ordem
  `G_F² T⁴/M_W²`, desprezado) nem polarização do meio.
- O integrador por janelas deriva `Δ` com `dY_νe/ds`. O teste é a
  **equivalência** com a tabela `Y_e' = Y_e + (4/3)Y_νe` sem fundo: estado, cota
  adiabática, `ε̃` fora da janela e posição do cruzamento iguais. O mutante que
  esquece `dY_νe/ds` passa na comparação com força bruta (erro 3×10⁻⁶ < 10⁻⁴) e
  só cai na equivalência; por isso ela é a verificação, e não a força bruta.

### 13.4 Traço sobre o espaço inteiro

No sistema ativo–estéril, só a identidade do **espaço inteiro** é fase global.
Remover o traço de cada bloco separadamente apaga a diferença física entre eles
e desloca a ressonância. A remoção do traço de `dae_nu_hamiltoniano` (que é por
sabor ativo) **não** é reaproveitada no módulo magnético.

Consequência para os testes: um Hamiltoniano novo e o antigo podem diferir por
um múltiplo da identidade. Comparações entre eles são feitas em
**probabilidades**, ou em operadores **depois de alinhar a fase global** — nunca
em amplitudes bit a bit. Pela mesma razão, transformações exatas (por exemplo
`H → aH`, `s → s/a`) são testadas dentro da tolerância declarada, não bit a bit:
caminhos equivalentes não garantem arredondamentos idênticos.

### 13.5 Validade física não é critério numérico

A janela `|Δ_ef| < K · 2|Ω|` do integrador por janelas é **parâmetro de
convergência**: dobrar `K` dobra a janela sem mudar a física. Ela não entra em
nenhum critério de validade. O domínio em que a conversão coerente vale é
decidido por grandezas físicas — `Γ/Ω_⋆`, `N_col` sobre a transição (`ℓ_res`,
`ℓ_LZ`) e `N_col` entre cruzamentos — e qualquer triagem só
vale depois de confirmada pelo teste amortecido, com erro absoluto em
`P` ≤ 10⁻³.

A defasagem de `dae_traj` é **aproximação diagnóstica**, não cota do termo de
colisão: colisões também mudam populações, direções e energia, e produzem `ν_R`
diretamente.

**Triagem de colisões** (mapa conjunto `γ_LZ × Γ/Ω`): o desvio `|ΔP|` NÃO
colapsa em `N_col`; no regime linear `|ΔP|/N_col` cresce com `γ_LZ` até 0,41
(perto de `γ_LZ = 5`). Vale `|ΔP| ≤ 0,41 N_col` na grade `γ_LZ ∈ [0,1; 10]`,
logo **`N_col < 2,4×10⁻³` basta para `|ΔP| < 10⁻³`** (condição suficiente
nessa grade, não necessária). A triagem ingênua (`Γ/Ω < 10⁻²` e
`N_col < 10⁻²`) é violada em 4 pontos da grade e **não deve ser usada**.

**Colisões: o operador e os fatores** (`dae_numag_colisoes.c`). Operador
físico `√Γ |L⟩⟨L|` (só o ativo interage): coerência amortecida a `Γ/2`.
Desdobramento usado: simétrico, `√γ P_j` nos dois estados, taxa total `γ`
independente do estado; mesmo Liouvilliano se **`γ = Γ/2`**, ou seja
`γ [1/km] = Γ [eV] · T / 2` (`dae_numag_salto_por_km`). Verificado em `t111`
contra `½ e^{−γs}` com `H = 0`. O par `{j, n+j}` do mergulho representa `P_j`;
não é o par `ν_L/ν_R`. Colisões só com `brute_force`: o modo incoerente das
janelas é hipótese ENTRE cruzamentos, não evolução amortecida DURANTE um.
Somas do ensemble por Welford, na ordem do índice.

### 13.5b O integrador por janelas: fora da janela não quer dizer adiabático

`dae_numag_propaga_janelas`. Regras que decidem o resultado:

- **O critério não é `|Δ|` grande.** Um trecho só é tratado como adiabático se,
  no referencial que gira com `B_⊥` (onde `H~` é real e a ressonância é o zero
  de `Δ_ef`), a cota de primeira ordem da amplitude de transição
  `|ε~(início)| + |ε~(fim)| + variação total de ε~ + Σ saltos de θ~ nos nós`
  ficar abaixo de `eps_tol` (ver o item seguinte: na implementação, um
  estimador amostrado), com `ε~ = θ~'/(T E~)` e
  `θ~' = (Δ_ef Ω' − Ω Δ_ef')/E~²`. A variação de `|B|` (Ω') e da orientação
  (φ'' dentro de `Δ_ef'`) entra aí. Uma modulação em ressonância com o gap
  acumula variação total e é pega mesmo com `ε~` local pequeno: `t110`, caso 3,
  com `|Δ| = 300 × 2Ω` fora da janela K — sem o critério, P_R sai errado; com
  ele, igual à força bruta.
- **O que se grava é um ESTIMADOR, não uma cota certificada.** A integração por partes dá a desigualdade com a variação total
  CONTÍNUA de `ε~`; o código usa a variação AMOSTRADA nos extremos e no meio de
  cada uma das `samples_per_interval` células — em geral um limite inferior
  dela. Nos perfis de um cruzamento, o estimador sobe 1 % de `sub` = 8 a 128 e
  a classificação não muda; isso não
  garante que uma estrutura mais estreita que uma célula seja vista. O
  candidato conhecido: `b` passando perto de zero DENTRO de uma célula com
  campo girante — `φ' = N/|b|²` tem um pico de largura `|b|_min/|b'|`, que pode
  esconder dois zeros de `Δ_ef` entre amostras. Tabelas de campo têm de pôr nós
  nesses mínimos de `|b|`. A desigualdade, além disso, limita a amplitude de
  1ª ordem, não o erro total do operador: quadratura de fase, força bruta nas
  janelas, ordens superiores e acúmulo entre trechos são outras parcelas.
  Apesar do nome, a coluna `cota_adiabatica` do CSV é esse estimador.
- **Zeros de `Δ_ef`** (`t113`): um zero exatamente num ponto de amostragem é
  localizado nele (sem bisseção) e pertence à célula à direita; a
  classificação pelos sinais dos dois lados separa cruzamento, contato
  (`contatos`) e intervalo identicamente nulo (`intervalos_nulos`,
  `comprimento_nulo_km`). Uma raiz ascendente no início da célula não pode ir
  para o fim dela, senão um perfil inteiramente ressonante vira um
  "cruzamento" por célula. **O lado esquerdo é o da trajetória, não o da
  chamada:** em `s0` interior à tabela, a classificação usa
  uma célula virtual na grade global de subcélulas antes de `s0`, e a
  contagem de intervalos nulos continua entre chamadas. Sem isso, o CLI — uma
  chamada por fatia de saída — transformaria um contato que caísse na
  fronteira de uma fatia em cruzamento, e um intervalo nulo em dois:
  diagnóstico dependente de `nt` (P_R não muda). Só no início da TABELA um
  zero seguido de valor não nulo conta como cruzamento sem lado esquerdo.
- **`K` é parâmetro de convergência, não de física** (13.5). A precisão é
  controlada por `eps_max` (o erro segue `ε~` nas bordas das janelas).
- **Nenhum passo de força bruta atravessa um nó da tabela.** O meio interpolado
  tem quina nos nós; dentro de um passo ela derruba a regra do ponto médio para
  1ª ordem e o controle de passo vai ao teto (medido: 2²¹ subpassos).
  `dae_numag_propaga_tabela` corta nos nós; o `ghost run` com `brute_force` usa
  ela.
- **Modo incoerente: nunca fechar uma janela no meio.** Dois pedaços de uma
  ressonância combinados por probabilidade descartam a coerência DENTRO dela
  (medido: P_R 0,839 contra 0,914). Entre fatias de saída a janela fica aberta
  (`dae_numag_janela_aberta`).
- **Comparar populações adiabáticas, não P_R, para isolar erro de integrador e
  de tabela.** P_R num extremo finito carrega o termo cruzado de borda
  (~1/K_ext), que depende da fase dinâmica (~10⁴–10⁵ rad) — e a fase muda O(1)
  rad com a tabulação. Partir do autoestado e medir a população adiabática no
  fim tira a borda da comparação (`t110`, casos 1 e 4).

### 13.6 Convergência

- Modo incoerente entre cruzamentos: populações dos autoestados de matéria.
- Modo coerente: populações **e** o operador de transferência de cada janela.
  Populações convergidas podem esconder interferência errada.

### 13.7 O que um teste de dois estados prova, e o que não prova

`P_{L→R}` é par em `Δ` e em `Ω`. Um teste que só compara `P` (como `t106`)
valida forma e escala de tempo, **nunca sinais**. Todo teste de sinal precisa
de um observável ímpar — no setor magnético, o campo girante.

### 13.8 Leitura dos resultados

Limiares fixados antes de qualquer mapa: um efeito é relevante quando
`L_R/L_total ≥ 1%` ou `|ΔL|/L_ref ≥ 1%` **e** supera 3× o erro numérico e
estatístico combinado; a robustez entre perfis é avaliada à parte. Nenhum efeito
acima dos limiares produz um **mapa de ausência de efeito acima dos limiares**,
nunca um "mapa de exclusão": excluir parâmetros exige confronto com dados.

### 13.9 Alvos de verificação

O setor magnético passa pela suíte nativa, `cxx-check` (amálgama em C99 e
C++17), ASan + UBSan, typecheck e vitest da interface. Um resultado magnético
só vale para o alvo do navegador e para o `.cpp` exportado depois de passar
também pelos portões de WebAssembly (`make -C native wasm`, `test6`,
`wasm-test`, `teste7`).
