(* ::Package:: *)

(*__DAE_BANNER__*)

(* ==========================================================================
   Daedalus — caminhadas quânticas de tempo contínuo sobre grafos

   Este pacote é AUTÔNOMO: gera a rede a partir dos parâmetros, propaga,
   calcula os observáveis, faz varredura e escreve o CSV no mesmo formato do
   núcleo em C. Não depende de biblioteca compilada, e é para ser LIDO.

   ---------------------------------------------------------------------------
   POR QUE ELE NÃO É UMA TRADUÇÃO DO NÚCLEO

   No resto do projeto existe um núcleo em C que compila para WebAssembly, para
   binário nativo e, por concatenação, para o .cpp exportado — mesmo código, três
   alvos, e a concordância entre eles é ESTRUTURAL: não há o que divergir.

   Aqui não é possível: Wolfram não inclui C, e uma biblioteca compilada seria
   caixa-preta, o oposto de inspeção. Então a garantia muda de natureza e passa a
   ser EMPÍRICA: dois métodos, mesmo resultado.

   Para que isso valha alguma coisa, o método tem de ser deliberadamente
   diferente. O núcleo propaga por expansão de Chebyshev; este pacote propaga por
   DECOMPOSIÇÃO ESPECTRAL (e por Krylov, via MatrixExp, quando N cresce).
   Reimplementar Chebyshev aqui seria pior que inútil: as duas versões poderiam
   compartilhar o mesmo erro conceitual, e a concordância não provaria nada.
   ---------------------------------------------------------------------------
*)

BeginPackage["Daedalus`"];

(* ---------------------------------------------------------- especificação *)

DaedalusEspecPadrao::usage =
  "DaedalusEspecPadrao[] devolve a Association com os padrões do núcleo. \
Normalização de \[LeftDoubleBracketingBar]H\[RightDoubleBracketingBar] LIGADA e \
religação mantendo |E| fixo: as duas são disciplina metodológica, não \
conveniência, e desligá-las é escolha explícita.";

DaedalusEspecLer::usage =
  "DaedalusEspecLer[\"arquivo.json\"] lê um spec.json e devolve a Association \
correspondente, completada com os padrões. Chave desconhecida é ERRO, como no \
parser em C: um campo escrito errado produziria um arquivo internamente \
consistente resolvendo o problema errado.";

DaedalusEspecCanonico::usage =
  "DaedalusEspecCanonico[espec] devolve a string JSON canônica, na mesma ordem \
de chaves do núcleo. É o par (spec, core_hash) que torna um resultado \
reproduzível \[LongDash] o spec sozinho não basta, porque o núcleo muda.";

(* ------------------------------------------------------------------- rede *)

DaedalusRede::usage =
  "DaedalusRede[espec] constrói a rede a partir dos PARÂMETROS, não de uma \
lista de arestas pronta \[LongDash] o gerador é para ser mexido. Devolve \
<|\"A\", \"N\", \"Modulos\", \"NModulos\", \"NPar\", \"NPerp\", \"Geometria\", \
\"ArestasDescartadas\", \"ReligacoesFalhas\", \"Digital\", \"CSR\"|>. \
Geradores: \"microtubule\", \"sbm\", \"path\", \"cycle\", \"grid2d\", \
\"hypercube\", \"complete\".";

DaedalusDigital::usage =
  "DaedalusDigital[rede] devolve a impressão digital FNV-1a de 64 bits da \
estrutura CSR, idêntica à de dae_graph_fingerprint em C. É a PRIMEIRA coisa que \
a verificação compara: \"o gerador discordou\" e \"o propagador discordou\" têm \
causas e correções diferentes, e sem a digital as duas chegam como \"os números \
não batem\".";

DaedalusFluxoAleatorio::usage =
  "DaedalusFluxoAleatorio[semente, quantidade] devolve os primeiros inteiros de \
64 bits do xoshiro256++ semeado por splitmix64. Exposto para inspeção: é a peça \
que precisa bater bit a bit com o núcleo para que redes estocásticas tenham a \
mesma digital, e uma divergência aqui não aparece como erro \[LongDash] aparece \
como grafo diferente.";

DaedalusMetricas::usage =
  "DaedalusMetricas[rede] devolve <|\"Lambda2\", \"Q\", \"GrauMedio\", \
\"Arestas\", \"Componentes\", \"CaminhoMedio\"|>. Aqui \[Lambda]2 vem de \
Eigenvalues da laplaciana e é EXATO, sem o problema de convergência do Lanczos \
que o núcleo enfrenta em rede fortemente modular. Para N moderado este é o valor \
de referência, e é o núcleo que se confere contra ele.";

(* ----------------------------------------------------------------- física *)

DaedalusHamiltoniano::usage =
  "DaedalusHamiltoniano[rede, espec] devolve <|\"H\", \"Escala\"|>. \
H = -\[Gamma]A (adjacência) ou \[Gamma]L (laplaciana), dividido pelo fator de \
normalização, que é medido em H cru com \[Gamma] = 1 \[LongDash] assim \[Gamma] \
continua sendo o botão físico e a normalização remove só a escala imposta pela \
topologia.";

DaedalusPropagar::usage =
  "DaedalusPropagar[ham, espec] evolui |\[Psi](t)\[RightAngleBracket] = \
exp(-iHt)|\[Psi](0)\[RightAngleBracket] na grade de tempo do espec.\n\n\
Method -> \"Espectral\" usa Eigensystem e soma espectral: exato, legível, e erra \
de forma DIFERENTE do Chebyshev do núcleo \[LongDash] é essa diferença de método \
que faz da concordância uma evidência.\n\
Method -> \"Krylov\" usa MatrixExp[-I H dt, \[Psi]], que aplica a exponencial ao \
vetor sem formá-la.\n\
Method -> Automatic escolhe pela ordem da matriz.\n\n\
Não há implementação de Chebyshev aqui, de propósito.";

DaedalusObservaveis::usage =
  "DaedalusObservaveis[propagacao, rede, espec] devolve <|\"Norma\", \"IPR\", \
\"CoerenciaL1\", \"PAlvo\", \"PopulacaoPorModulo\", \"ConcurrencePorModulo\", \
\"Populacao\"|>.\n\n\
Sem sítio-alvo, PAlvo é Missing[\"NotAvailable\"] \[LongDash] nunca zero, que é \
valor fisicamente válido e seria plotado como se fosse dado.";

(* -------------------------------------------------------- varredura, saída *)

DaedalusVarredura::usage =
  "DaedalusVarredura[espec, {pMin, pMax, passos}, realizacoes] varre a fração \
religada e devolve, para cada p, a média TEMPORAL de p_alvo com desvio e número \
de realizações. Média temporal, não valor final: o valor final oscila.";

DaedalusCSV::usage =
  "DaedalusCSV[espec, rede, metricas, propagacao, observaveis] devolve o CSV no \
mesmo formato do núcleo, com cabeçalho de procedência, para que o resultado \
possa ser reimportado na interface. A linha \"#! implementacao wolfram\" diz \
QUEM calculou; \"#! core_hash\" diz de qual versão do projeto o spec veio.";

DaedalusExportarCSV::usage =
  "DaedalusExportarCSV[\"arquivo.csv\", espec, rede, metricas, propagacao, \
observaveis] escreve o CSV.";

DaedalusRodar::usage =
  "DaedalusRodar[espec] faz tudo \[LongDash] rede, métricas, hamiltoniano, \
propagação e observáveis \[LongDash] e devolve uma Association com as partes. É \
a linha única para quem quer o resultado; as funções acima são para quem quer \
mexer.";

DaedalusVerificarReferencia::usage =
  "DaedalusVerificarReferencia[\"pasta\"] compara este pacote com os CSV de \
referência gerados pelo núcleo em C. Confere, nesta ordem: o fluxo do PRNG, a \
digital do grafo (exigida IDÊNTICA), e só então os observáveis, com tolerância \
declarada.\n\n\
Esta verificação é MANUAL: a integração contínua do projeto não tem Mathematica.";

DaedalusLindblad::usage =
  "DaedalusLindblad[espec] resolve a equação de Lindblad de defasagem pura \
EXATAMENTE, montando o superoperador vetorizado e exponenciando-o. Devolve \
<|\"Tempos\", \"Rho\", \"Escala\"|>.\n\n\
Este é o método que NÃO é trajetória nenhuma. A média de trajetórias do núcleo \
em C tem de convergir para ele — e sem essa companheira, \"as trajetórias \
convergiram\" e \"as trajetórias convergiram para a coisa errada\" são \
indistinguíveis: com 200 realizações a barra de erro parece respeitável nos dois \
casos.";

DaedalusVerificarTrajetorias::usage =
  "DaedalusVerificarTrajetorias[{arquivos}] compara os ρ médios que \
`daedalus traj` escreveu contra a solução exata de Lindblad.\n\n\
O critério NÃO é um limiar de desvio: é a LEI DE CONVERGÊNCIA. O erro de Monte \
Carlo cai como 1/√n_traj, então quadruplicar as trajetórias tem de reduzir o \
desvio pela metade. Um desdobramento enviesado converge — para outro lugar — e \
o desvio ESTABILIZA em vez de cair. É a diferença que um limiar sozinho não vê.";

DaedalusAntiVacuidade::usage =
  "DaedalusAntiVacuidade[\"pasta\"] mostra que DaedalusVerificarReferencia PODE \
reprovar, e reprovar pelo motivo certo. Sabota deliberadamente o grafo, o fluxo \
de sorteios, os números e o próprio arquivo, e exige que cada sabotagem produza \
o veredito correspondente.\n\n\
Sem isto, \"todos os casos passaram\" não é informação: uma verificação sempre \
verde e uma verificação correta são indistinguíveis pelo resultado. Rode as duas \
juntas, sempre.";

Begin["`Private`"];

(* ==========================================================================
   PRNG — xoshiro256++ semeado por splitmix64

   Reimplementado aqui em aritmética exata módulo 2^64. É a peça mais frágil do
   pacote: se ela divergir do núcleo, redes estocásticas (SBM, qualquer
   religação) saem DIFERENTES, e a discordância não se apresenta como erro
   numérico, se apresenta como física diferente.

   O estado é mutável de propósito, espelhando o C: a ORDEM dos sorteios faz
   parte do contrato de reprodutibilidade, e uma versão funcional esconderia
   essa ordem em vez de exibi-la.
   ========================================================================== *)

$mascara64 = 2^64 - 1;
$estadoRNG = {0, 0, 0, 0};

girarEsq[x_, k_] :=
  BitAnd[BitOr[BitShiftLeft[x, k], BitShiftRight[x, 64 - k]], $mascara64];

splitmix64[z_] := Module[{t},
  t = BitAnd[BitXor[z, BitShiftRight[z, 30]]*16^^BF58476D1CE4E5B9, $mascara64];
  t = BitAnd[BitXor[t, BitShiftRight[t, 27]]*16^^94D049BB133111EB, $mascara64];
  BitXor[t, BitShiftRight[t, 31]]];

semear[semente_Integer] := Module[{z = BitAnd[semente, $mascara64]},
  $estadoRNG = Table[
    z = BitAnd[z + 16^^9E3779B97F4A7C15, $mascara64]; splitmix64[z], {4}];
  (* estado todo-zero é ponto fixo do xoshiro; splitmix nunca produz isso, mas
     a guarda é barata e o custo de errar seria silencioso *)
  If[Total[$estadoRNG] == 0, $estadoRNG[[1]] = 16^^9E3779B97F4A7C15];
  $estadoRNG];

proximoU64[] := Module[{s = $estadoRNG, r, t},
  r = BitAnd[girarEsq[BitAnd[s[[1]] + s[[4]], $mascara64], 23] + s[[1]], $mascara64];
  t = BitAnd[BitShiftLeft[s[[2]], 17], $mascara64];
  s[[3]] = BitXor[s[[3]], s[[1]]];
  s[[4]] = BitXor[s[[4]], s[[2]]];
  s[[2]] = BitXor[s[[2]], s[[3]]];
  s[[1]] = BitXor[s[[1]], s[[4]]];
  s[[3]] = BitXor[s[[3]], t];
  s[[4]] = girarEsq[s[[4]], 45];
  $estadoRNG = s;
  r];

(* [0,1) com 53 bits, exatamente como o núcleo: o SBM decide cada aresta com
   isto, e meio ulp de diferença muda o grafo *)
proximoUniforme[] := N[BitShiftRight[proximoU64[], 11]/9007199254740992];

(* Inteiro em [0,n) sem viés, por rejeição — o mesmo laço do núcleo, inclusive
   no número de sorteios que ele consome quando rejeita *)
proximoAbaixo[n_Integer] := Module[{faixa = n, blocos, limite, r},
  If[n <= 1, Return[0]];
  blocos = Quotient[$mascara64, faixa];
  limite = blocos*faixa;
  While[True, r = proximoU64[]; If[r < limite, Return[Mod[r, faixa]]]]];

DaedalusFluxoAleatorio[semente_Integer, quantidade_Integer] := (
  semear[semente];
  Table[proximoU64[], {quantidade}]);

(* ==========================================================================
   CSR — a mesma construção do núcleo, incluindo a política de duplicatas

   Aresta repetida é DESCARTADA, nunca somada: com seam_shift e N_perp pequeno
   dá para gerar duplicata legítima, e somar dobraria j_perp em alguns sítios
   sem nenhum aviso. A ordem importa: dentro de cada linha, ordenação estável
   por coluna e a PRIMEIRA ocorrência é que fica.
   ========================================================================== *)

construirCSR[n_Integer, arestas_List] := Module[
  {dirigidas, porLinha, linhas, colunas, valores, contagem, descartadas = 0,
   colFinal = {}, valFinal = {}, ponteiro},

  (* cada aresta não-dirigida vira duas entradas, salvo o laço *)
  dirigidas = Flatten[
    Map[If[#[[1]] === #[[2]], {{#[[1]], #[[2]], #[[3]]}},
           {{#[[1]], #[[2]], #[[3]]}, {#[[2]], #[[1]], #[[3]]}}] &, arestas], 1];

  porLinha = ConstantArray[{}, n];
  If[Length[dirigidas] > 0,
    Module[{grupos = GatherBy[dirigidas, First]},
      Do[porLinha[[grupos[[k, 1, 1]] + 1]] = grupos[[k]], {k, Length[grupos]}]]];

  contagem = ConstantArray[0, n];
  Do[Module[{g = porLinha[[i]], ordem, cols, vals, mantidas},
      If[Length[g] == 0, Continue[]];
      ordem = Ordering[g[[All, 2]]];         (* estável: empate mantém a ordem *)
      cols = g[[ordem, 2]]; vals = g[[ordem, 3]];
      mantidas = Reap[
        Do[If[k == 1 || cols[[k]] =!= cols[[k - 1]], Sow[k], descartadas++],
           {k, Length[cols]}]][[2]];
      mantidas = If[mantidas === {}, {}, First[mantidas]];
      contagem[[i]] = Length[mantidas];
      colFinal = Join[colFinal, cols[[mantidas]]];
      valFinal = Join[valFinal, vals[[mantidas]]]],
    {i, n}];

  ponteiro = Prepend[Accumulate[contagem], 0];
  <|"N" -> n, "RowPtr" -> ponteiro, "ColInd" -> colFinal, "Val" -> valFinal,
    "NNZ" -> Length[colFinal], "Descartadas" -> descartadas|>];

matrizDeCSR[csr_Association] := Module[{n = csr["N"], pares},
  If[csr["NNZ"] == 0, Return[SparseArray[{}, {n, n}]]];
  pares = Flatten[Table[
    Table[{i, csr["ColInd"][[p]] + 1}, {p, csr["RowPtr"][[i]] + 1, csr["RowPtr"][[i + 1]]}],
    {i, n}], 1];
  SparseArray[pares -> csr["Val"], {n, n}]];

(* ==========================================================================
   Impressão digital — FNV-1a de 64 bits sobre rowptr, colind, os BITS dos
   pesos e a partição em módulos. Idêntica a dae_graph_fingerprint.

   Os pesos entram pelo padrão de bits IEEE-754, não pelo valor: j_par trocado
   por j_perp muda a digital mesmo com a estrutura idêntica.
   ========================================================================== *)

bitsDeDouble[x_] := First@ImportString[
  ExportString[{N[x]}, "Binary", "DataFormat" -> "Real64"],
  "Binary", "DataFormat" -> "UnsignedInteger64"];

fnvMistura[h_, v_] := BitAnd[BitXor[h, BitAnd[v, $mascara64]]*1099511628211, $mascara64];

(* A base é 1469598103934665603, e NÃO a constante canônica do FNV-1a
   (14695981039346656037): o núcleo em C usa este valor, e a digital só serve
   se os dois lados usarem a MESMA semente, canônica ou não. Copiar a constante
   "certa" da literatura em vez da que está no outro programa faria as duas
   digitais divergirem sempre, para todo grafo. *)
DaedalusDigital[rede_Association] := Module[{csr = rede["CSR"], h = 1469598103934665603},
  h = Fold[fnvMistura, h, csr["RowPtr"]];
  h = Fold[fnvMistura, h, csr["ColInd"]];
  h = Fold[fnvMistura, h, bitsDeDouble /@ csr["Val"]];
  h = Fold[fnvMistura, h, rede["Modulos"]];
  h];

(* ==========================================================================
   Geradores

   A ORDEM em que as arestas entram importa duas vezes: ela decide qual
   duplicata sobrevive na CSR, e nos geradores estocásticos ela decide a
   sequência de sorteios. Por isso os laços aqui seguem a mesma ordem do núcleo,
   e não uma reorganização "mais idiomática" que produziria outro grafo.
   ========================================================================== *)

geraMicrotubulo[p_Association] := Module[
  {np = p["n_par"], nq = p["n_perp"], desloc = p["seam_shift"],
   fechado = TrueQ[p["longitudinal_closed"]], arestas, ms, mm},
  arestas = Reap[
    Do[
      Module[{j = m*nq + q},
        If[m + 1 < np, Sow[{j, (m + 1)*nq + q, p["j_par"]}],
          If[fechado, Sow[{j, q, p["j_par"]}]]];
        If[q + 1 < nq, Sow[{j, m*nq + q + 1, p["j_perp"]}],
          (* COSTURA: (m, N\[Perpendicular]-1) liga com (m + desloc, 0).
             Com as pontas abertas, as arestas que caem fora do cilindro somem —
             e somem em protofilamentos DIFERENTES nas duas extremidades. É a
             assimetria plus/minus end, física, não defeito de montagem. *)
          If[nq > 2,
            ms = m + desloc;
            If[fechado,
              mm = Mod[ms, np]; Sow[{j, mm*nq, p["j_perp"]}],
              If[0 <= ms < np, Sow[{j, ms*nq, p["j_perp"]}]]]]]],
      {m, 0, np - 1}, {q, 0, nq - 1}]][[2]];
  If[arestas === {}, {}, First[arestas]]];

geraSBM[p_Association, modulos_List] := Module[{n = p["n"], arestas},
  arestas = Reap[
    Do[If[proximoUniforme[] < If[modulos[[i + 1]] === modulos[[j + 1]], p["p_in"], p["p_out"]],
          Sow[{i, j, 1.0}]],
      {i, 0, n - 1}, {j, i + 1, n - 1}]][[2]];
  If[arestas === {}, {}, First[arestas]]];

geraReferencia[tipo_String, p_Association] := Module[{n, arestas},
  Switch[tipo,
    "path" | "cycle",
      n = p["n"];
      arestas = Table[{i, i + 1, 1.0}, {i, 0, n - 2}];
      If[tipo === "cycle", AppendTo[arestas, {n - 1, 0, 1.0}]];
      arestas,
    "complete",
      n = p["n"];
      Flatten[Table[{i, j, 1.0}, {i, 0, n - 1}, {j, i + 1, n - 1}], 1],
    "hypercube",
      n = 2^p["dim"];
      Flatten[Table[
        Module[{k = BitXor[i, 2^b]}, If[k > i, {{i, k, 1.0}}, {}]],
        {i, 0, n - 1}, {b, 0, p["dim"] - 1}], 2],
    "grid2d",
      Flatten[Table[
        Module[{j = r*p["cols"] + c, saida = {}},
          If[c + 1 < p["cols"], AppendTo[saida, {j, j + 1, 1.0}]];
          If[r + 1 < p["rows"], AppendTo[saida, {j, j + p["cols"], 1.0}]];
          saida],
        {r, 0, p["rows"] - 1}, {c, 0, p["cols"] - 1}], 2],
    _, {}]];

(* Religação Watts-Strogatz.
   No modo "rewire" o número de arestas é INVARIANTE: cada aresta sorteada troca
   um extremo, e a troca só é aceita se não colidir com aresta existente. Uma
   colisão aceita seria descartada depois pela deduplicação da CSR, |E| cairia
   sem aviso, e transporte pior por menos arestas seria lido como efeito de
   topologia. *)
chaveAresta[a_, b_, n_] := Min[a, b]*n + Max[a, b] + 2;

religar[arestas0_List, n_Integer, p_, modo_String] := Module[
  {arestas = arestas0, base = Length[arestas0], conjunto, falhas = 0, chave, aceito},
  If[p <= 0 || base == 0, Return[{arestas, 0}]];
  conjunto = Association[
    Table[chaveAresta[arestas[[e, 1]], arestas[[e, 2]], n] -> True, {e, base}]];
  Do[
    If[proximoUniforme[] < p,
      aceito = False;
      Do[
        Module[{b = proximoAbaixo[n]},
          If[b =!= arestas[[e, 1]],
            chave = chaveAresta[arestas[[e, 1]], b, n];
            If[! KeyExistsQ[conjunto, chave],
              If[modo === "rewire",
                conjunto = KeyDrop[conjunto,
                  chaveAresta[arestas[[e, 1]], arestas[[e, 2]], n]];
                arestas[[e, 2]] = b,
                AppendTo[arestas, {arestas[[e, 1]], b, arestas[[e, 3]]}]];
              conjunto[chave] = True;
              aceito = True]]];
        If[aceito, Break[]],
        {100}];
      If[! aceito, falhas++]],
    {e, base}];
  {arestas, falhas}];

DaedalusRede[espec_Association] := Module[
  {g = espec["graph"], p, tipo, n, np, nq, nmod, arestas, modulos, falhas = 0, csr, rede},
  tipo = g["generator"];
  p = g["params"];
  semear[espec["seed"]];

  Switch[tipo,
    "microtubule",
      np = p["n_par"]; nq = p["n_perp"]; n = np*nq;
      nmod = Max[1, p["n_modules"]];
      modulos = Flatten[Table[Quotient[m*nmod, np], {m, 0, np - 1}, {q, nq}]];
      arestas = geraMicrotubulo[p],
    "sbm",
      n = p["n"]; np = 0; nq = 0; nmod = Max[1, p["n_modules"]];
      modulos = Table[Quotient[i*nmod, n], {i, 0, n - 1}];
      arestas = geraSBM[p, modulos],
    "hypercube",
      n = 2^p["dim"]; np = 0; nq = 0; nmod = 1;
      modulos = ConstantArray[0, n]; arestas = geraReferencia[tipo, p],
    "grid2d",
      n = p["rows"]*p["cols"]; np = 0; nq = 0; nmod = 1;
      modulos = ConstantArray[0, n]; arestas = geraReferencia[tipo, p],
    _,
      n = p["n"]; np = 0; nq = 0; nmod = 1;
      modulos = ConstantArray[0, n]; arestas = geraReferencia[tipo, p]];

  {arestas, falhas} = religar[arestas, n, p["ws_p"], p["conn_mode"]];
  csr = construirCSR[n, arestas];
  rede = <|"N" -> n, "CSR" -> csr, "A" -> matrizDeCSR[csr],
    "Modulos" -> modulos, "NModulos" -> nmod, "NPar" -> np, "NPerp" -> nq,
    "Geometria" -> Which[tipo === "microtubule", 2,
                         MemberQ[{"path", "cycle", "grid2d"}, tipo], 1, True, 0],
    "ArestasDescartadas" -> csr["Descartadas"], "ReligacoesFalhas" -> falhas|>;
  Append[rede, "Digital" -> DaedalusDigital[rede]]];

(* ==========================================================================
   Hamiltoniano e propagação
   ========================================================================== *)

Options[DaedalusHamiltoniano] = {"EscalaImposta" -> Automatic};

DaedalusHamiltoniano[rede_Association, espec_Association, OptionsPattern[]] := Module[
  {h = espec["hamiltonian"], A = rede["A"], n = rede["N"], Hcru, escala, graus,
   imposta = OptionValue["EscalaImposta"]},
  graus = Total[A, {2}];
  Hcru = If[h["kind"] === "adjacency", -A, DiagonalMatrix[SparseArray[graus]] - A];
  escala = Switch[h["normalization"],
    "none", 1.0,
    "mean_degree", Total[Abs[Flatten[Normal[A]]]]/n,
    (* EXATO aqui, ESTIMADO no núcleo. O núcleo usa poucos passos de Lanczos
       com margem, e o valor que ele obtém é maior que o raio verdadeiro. Os
       dois programas discordam nesta escala por construção, e a verificação
       reporta a diferença em vez de exigir igualdade. *)
    "spectral", Max[Abs[Eigenvalues[N[Normal[Hcru]]]]],
    _, 1.0];
  If[! (escala > 0), escala = 1.0];
  (* Com "EscalaImposta" o fator vem de fora — é assim que a verificação
     compara o PROPAGADOR contra o núcleo sem que a diferença entre a estimativa
     de Lanczos dele e o raio exato daqui contamine todas as colunas. A escala
     calculada continua sendo devolvida em "EscalaExata", para que a diferença
     apareça como número em vez de sumir. *)
  <|"H" -> SparseArray[(h["gamma"]/If[NumericQ[imposta] && imposta > 0, imposta, escala])*Hcru],
    "Escala" -> If[NumericQ[imposta] && imposta > 0, imposta, escala],
    "EscalaExata" -> escala|>];



Options[DaedalusPropagar] = {Method -> Automatic};

DaedalusPropagar[ham_Association, espec_Association, OptionsPattern[]] := Module[
  {H = ham["H"], n = Length[ham["H"]], nt = espec["time"]["nt"], t1 = espec["time"]["t1"],
   psi0, dt, metodo, tempos, estados, vals, vecs, coef},
  dt = t1/nt;
  tempos = Table[k*dt, {k, nt}];
  psi0 = ConstantArray[0.0 + 0.0 I, n];
  psi0[[espec["initial"]["site"] + 1]] = 1.0 + 0.0 I;

  metodo = OptionValue[Method] /. Automatic -> If[n <= 1500, "Espectral", "Krylov"];

  estados = Switch[metodo,
    "Espectral",
      (* Decomposição espectral: exata, e o erro dela não se parece com o do
         Chebyshev do núcleo. É essa diferença que faz a concordância significar
         alguma coisa. *)
      {vals, vecs} = Eigensystem[N[Normal[H]]];
      vecs = Normalize /@ vecs;
      coef = vecs . psi0;
      Table[Transpose[vecs] . (Exp[-I vals t] coef), {t, tempos}],
    _,
      (* Krylov: MatrixExp com o terceiro argumento aplica a exponencial ao
         vetor sem formá-la. Uma linha, em vez de um laço de recorrência. *)
      Rest[FoldList[MatrixExp[-I N[H] dt, #1] &, psi0, Range[nt]]]];

  <|"Tempos" -> tempos, "Estados" -> estados, "Metodo" -> metodo, "dt" -> dt|>];

(* ==========================================================================
   Observáveis

   CONVENÇÃO DA CONCURRENCE, que precisa fechar com a coerência ℓ1. Para estado
   puro, C_ij = 2|ψ_i||ψ_j|. Somando sobre PARES NÃO-ORDENADOS {i,j}, i ≠ j:

       C_MN = 2 s_M s_N        (M ≠ N)
       C_MM = s_M² − q_M       (o produto externo incluiria i == j, que não é
                                par; o termo espúrio é exatamente q_M)

   A soma do triângulo superior COM a diagonal reproduz C_ℓ1 = (Σ|ψ|)² − 1, e é
   essa identidade — não a fórmula copiada — que o teste verifica.
   ========================================================================== *)

observaveisDeUmEstado[psi_List, modulos_List, nmod_Integer, alvo_Integer] := Module[
  {p = Abs[psi]^2, amp = Abs[psi], sM, qM, conc},
  sM = Table[Total[Pick[amp, modulos, m]], {m, 0, nmod - 1}];
  qM = Table[Total[Pick[p, modulos, m]], {m, 0, nmod - 1}];
  conc = Table[If[a === b, sM[[a]]^2 - qM[[a]], 2 sM[[a]] sM[[b]]],
    {a, nmod}, {b, nmod}];
  <|"Norma" -> Total[p], "IPR" -> Total[p^2],
    "CoerenciaL1" -> Total[amp]^2 - 1,
    (* Sem alvo: Missing, o análogo do NaN. Zero é valor fisicamente válido e
       seria plotado como se fosse dado. *)
    "PAlvo" -> If[alvo >= 0, p[[alvo + 1]], Missing["NotAvailable"]],
    "PopulacaoPorModulo" -> qM, "ConcurrencePorModulo" -> conc,
    "Populacao" -> p|>];

DaedalusObservaveis[prop_Association, rede_Association, espec_Association] := Module[
  {alvo = espec["observables"]["target"], serie},
  If[alvo >= rede["N"], alvo = -1];
  serie = observaveisDeUmEstado[#, rede["Modulos"], rede["NModulos"], alvo] & /@ prop["Estados"];
  (* Os parênteses não são decorativos: `&` liga mais fraco que `->`, e sem
     eles `"Norma" -> #["Norma"] &` vira uma função cujo corpo é a regra. *)
  <|"Norma" -> (#["Norma"] & /@ serie),
    "IPR" -> (#["IPR"] & /@ serie),
    "CoerenciaL1" -> (#["CoerenciaL1"] & /@ serie),
    "PAlvo" -> (#["PAlvo"] & /@ serie),
    "PopulacaoPorModulo" -> (#["PopulacaoPorModulo"] & /@ serie),
    "ConcurrencePorModulo" -> (#["ConcurrencePorModulo"] & /@ serie),
    "Populacao" -> (#["Populacao"] & /@ serie)|>];

(* ==========================================================================
   Métricas de rede

   λ2 aqui é EXATO: segundo menor autovalor da laplaciana, por Eigenvalues. O
   núcleo o obtém por Lanczos deflacionado com critério de convergência, e em
   rede fortemente modular a convergência é lenta. Para N moderado, este
   valor é a referência.
   ========================================================================== *)

DaedalusMetricas[rede_Association] := Module[
  {A = rede["A"], n = rede["N"], nmod = rede["NModulos"], modulos = rede["Modulos"],
   L, autovals, graus, W, Q, grafo, comps, dists, finitas},
  graus = Total[A, {2}];
  L = DiagonalMatrix[SparseArray[graus]] - A;
  autovals = Sort[Re[Eigenvalues[N[Normal[L]]]]];
  W = Total[graus]/2;
  Q = If[W > 0,
    Total[Table[
      Module[{dentro = Total[Flatten[Normal[A[[#, #]]]]] &@ Flatten[Position[modulos, m]],
              kM = Total[Pick[graus, modulos, m]]},
        dentro/(2 W) - (kM/(2 W))^2], {m, 0, nmod - 1}]], 0];
  grafo = AdjacencyGraph[Unitize[Normal[A]]];
  comps = Length[ConnectedComponents[grafo]];
  dists = GraphDistanceMatrix[grafo];
  finitas = Select[Flatten[dists], # > 0 && # < Infinity &];
  <|"Lambda2" -> autovals[[2]], "Q" -> Q,
    "GrauMedio" -> N[rede["CSR"]["NNZ"]/n],
    "Arestas" -> Quotient[rede["CSR"]["NNZ"], 2],
    "Componentes" -> comps,
    "CaminhoMedio" -> If[Length[finitas] > 0, N[Mean[finitas]], 0.0]|>];

(* ==========================================================================
   Saída — o mesmo CSV do núcleo
   ========================================================================== *)

formatoG17[x_] := Module[{v, sinal, digitos, expo},
  If[MissingQ[x], Return["nan"]];
  v = N[x];
  If[! NumberQ[v], Return["nan"]];
  If[v == 0, Return["0"]];
  sinal = If[v < 0, "-", ""];
  {digitos, expo} = RealDigits[Abs[SetPrecision[v, 20]], 10, 17];
  While[Length[digitos] > 1 && Last[digitos] === 0, digitos = Most[digitos]];
  sinal <> Which[
    expo - 1 < -4 || expo - 1 >= 17,
      StringJoin[ToString /@ Take[digitos, 1]] <>
        If[Length[digitos] > 1, "." <> StringJoin[ToString /@ Rest[digitos]], ""] <>
        "e" <> If[expo - 1 >= 0, "+", "-"] <> IntegerString[Abs[expo - 1], 10, 2],
    expo <= 0,
      "0." <> StringJoin[Table["0", {-expo}]] <> StringJoin[ToString /@ digitos],
    expo >= Length[digitos],
      StringJoin[ToString /@ digitos] <> StringJoin[Table["0", {expo - Length[digitos]}]],
    True,
      StringJoin[ToString /@ Take[digitos, expo]] <> "." <>
        StringJoin[ToString /@ Drop[digitos, expo]]]];

DaedalusEspecCanonico[espec_Association] := Module[{g = espec["graph"], p = espec["graph"]["params"]},
  StringJoin[
    "{\"format_version\":", ToString[espec["format_version"]],
    ",\"core_hash\":\"", espec["core_hash"], "\"",
    ",\"seed\":", ToString[espec["seed"]],
    ",\"graph\":{\"generator\":\"", g["generator"], "\"",
    ",\"params\":{\"n_par\":", ToString[p["n_par"]],
    ",\"n_perp\":", ToString[p["n_perp"]],
    ",\"seam_shift\":", ToString[p["seam_shift"]],
    ",\"longitudinal_closed\":", If[TrueQ[p["longitudinal_closed"]], "true", "false"],
    ",\"j_par\":", formatoG17[p["j_par"]],
    ",\"j_perp\":", formatoG17[p["j_perp"]],
    ",\"n_modules\":", ToString[p["n_modules"]],
    ",\"ws_p\":", formatoG17[p["ws_p"]],
    ",\"conn_mode\":\"", p["conn_mode"], "\"",
    ",\"p_in\":", formatoG17[p["p_in"]],
    ",\"p_out\":", formatoG17[p["p_out"]],
    ",\"n\":", ToString[p["n"]],
    ",\"dim\":", ToString[p["dim"]],
    ",\"rows\":", ToString[p["rows"]],
    ",\"cols\":", ToString[p["cols"]], "}}",
    ",\"hamiltonian\":{\"kind\":\"", espec["hamiltonian"]["kind"], "\"",
    ",\"gamma\":", formatoG17[espec["hamiltonian"]["gamma"]],
    ",\"normalization\":\"", espec["hamiltonian"]["normalization"], "\"",
    ",\"lanczos_steps\":", ToString[espec["hamiltonian"]["lanczos_steps"]], "}",
    ",\"initial\":{\"site\":", ToString[espec["initial"]["site"]], "}",
    ",\"time\":{\"t1\":", formatoG17[espec["time"]["t1"]],
    ",\"nt\":", ToString[espec["time"]["nt"]], "}",
    ",\"observables\":{\"target\":", ToString[espec["observables"]["target"]],
    ",\"population\":", If[TrueQ[espec["observables"]["population"]], "true", "false"],
    ",\"module_concurrence\":", If[TrueQ[espec["observables"]["module_concurrence"]], "true", "false"],
    ",\"full_concurrence\":", If[TrueQ[espec["observables"]["full_concurrence"]], "true", "false"],
    ",\"pop_stride\":", ToString[espec["observables"]["pop_stride"]], "}",
    ",\"realizations\":", ToString[espec["realizations"]], "}"]];

DaedalusCSV[espec_Association, rede_Association, metricas_Association,
            prop_Association, obs_Association] := Module[
  {nmod = rede["NModulos"], nt = Length[prop["Tempos"]], linhas},
  linhas = {
    "#! daedalus " <> espec["core_versao"],
    "#! core_hash " <> espec["core_hash"],
    (* QUEM calculou. core_hash diz de qual versão do projeto o spec veio; esta
       linha diz qual programa produziu os números. Sem ela, este CSV se
       apresentaria na interface como saída do núcleo em C. *)
    "#! implementacao wolfram",
    "#! metodo " <> prop["Metodo"],
    "#! graph_fingerprint " <> ToString[rede["Digital"]],
    "#! n " <> ToString[rede["N"]],
    "#! nnz " <> ToString[rede["CSR"]["NNZ"]],
    "#! nmod " <> ToString[nmod],
    "#! scale " <> formatoG17[espec["escala"]],
    "#! dt " <> formatoG17[prop["dt"]],
    "#! lambda2 " <> formatoG17[metricas["Lambda2"]],
    "#! modularity_Q " <> formatoG17[metricas["Q"]],
    "#! mean_degree " <> formatoG17[metricas["GrauMedio"]],
    "#! mean_path_len " <> formatoG17[metricas["CaminhoMedio"]],
    "#! n_edges " <> ToString[metricas["Arestas"]],
    "#! n_components " <> ToString[metricas["Componentes"]],
    "#! spec " <> DaedalusEspecCanonico[espec],
    "t,norm,ipr,coh_l1,p_target" <> StringJoin[Table[",pmod" <> ToString[m], {m, 0, nmod - 1}]]};
  linhas = Join[linhas,
    Table[StringRiffle[Join[
      {formatoG17[prop["Tempos"][[k]]], formatoG17[obs["Norma"][[k]]],
       formatoG17[obs["IPR"][[k]]], formatoG17[obs["CoerenciaL1"][[k]]],
       formatoG17[obs["PAlvo"][[k]]]},
      formatoG17 /@ obs["PopulacaoPorModulo"][[k]]], ","], {k, nt}]];
  linhas = Join[linhas, {"", "# estado final", "j,re,im"},
    Table[ToString[j - 1] <> "," <> formatoG17[Re[Last[prop["Estados"]][[j]]]] <>
          "," <> formatoG17[Im[Last[prop["Estados"]][[j]]]], {j, rede["N"]}]];
  StringRiffle[linhas, "\n"] <> "\n"];

DaedalusExportarCSV[arquivo_String, espec_, rede_, metricas_, prop_, obs_] :=
  Export[arquivo, DaedalusCSV[espec, rede, metricas, prop, obs], "Text"];

(* ==========================================================================
   Conveniência e varredura
   ========================================================================== *)

Options[DaedalusRodar] = {"EscalaImposta" -> Automatic};

DaedalusRodar::semtrajetoria = "Este spec pede `1` trajetórias, e o pacote \
Wolfram é unitário: ele não implementa o desdobramento de Haken-Strobl. O \
resultado seria a dinâmica sem ruído, com aparência de resposta. Rode as \
trajetórias pelo núcleo em C.";

DaedalusRodar[espec_Association, OptionsPattern[]] := Module[
  {rede, met, ham, prop, obs, e = espec},
  If[Lookup[Lookup[espec, "trajectories", <||>], "n_traj", 0] > 0,
    Message[DaedalusRodar::semtrajetoria,
      espec["trajectories"]["n_traj"]]; Abort[]];
  rede = DaedalusRede[e];
  met = DaedalusMetricas[rede];
  ham = DaedalusHamiltoniano[rede, e, "EscalaImposta" -> OptionValue["EscalaImposta"]];
  e = Append[e, "escala" -> ham["Escala"]];
  prop = DaedalusPropagar[ham, e];
  obs = DaedalusObservaveis[prop, rede, e];
  <|"Espec" -> e, "Rede" -> rede, "Metricas" -> met, "Hamiltoniano" -> ham,
    "Propagacao" -> prop, "Observaveis" -> obs|>];

DaedalusVarredura[espec_Association, {pMin_, pMax_, passos_Integer}, realizacoes_Integer] :=
  Table[
    Module[{p = If[passos == 1, pMin, pMin + (pMax - pMin) (k - 1)/(passos - 1)], amostras},
      amostras = Table[
        Module[{e = espec, r},
          e["graph"]["params"]["ws_p"] = N[p];
          e["seed"] = espec["seed"] + rr - 1;
          r = DaedalusRodar[e];
          (* média TEMPORAL de p_alvo — o valor final oscila *)
          Mean[DeleteMissing[r["Observaveis"]["PAlvo"]]]],
        {rr, realizacoes}];
      <|"p" -> N[p], "Media" -> Mean[amostras],
        "Desvio" -> If[Length[amostras] > 1, StandardDeviation[amostras], 0.0],
        "N" -> Length[amostras]|>],
    {k, passos}];

(* ==========================================================================
   Especificação
   ========================================================================== *)

DaedalusEspecPadrao[] := <|
  "format_version" -> 2, "core_hash" -> "", "core_versao" -> "0.1.0", "seed" -> 12345,
  (* Trajetórias. O bloco existe no padrão só para que o parser estrito ACEITE um
     spec que o traga; este pacote não implementa trajetórias, e diz isso em vez
     de recusar com "chave desconhecida", que mandaria o leitor procurar erro de
     digitação onde há funcionalidade ausente. *)
  "trajectories" -> <|"n_traj" -> 0, "gamma_deph" -> 0., "rho_stride" -> 1,
                      "output_mode" -> "accumulate_rho"|>,
  "graph" -> <|"generator" -> "path",
    "params" -> <|"n_par" -> 200, "n_perp" -> 13, "seam_shift" -> 3,
      "longitudinal_closed" -> False, "j_par" -> 1.0, "j_perp" -> 1.0,
      "n_modules" -> 1, "ws_p" -> 0.0, "conn_mode" -> "rewire",
      "p_in" -> 0.0, "p_out" -> 0.0, "n" -> 64, "dim" -> 6,
      "rows" -> 8, "cols" -> 8|>|>,
  (* As duas disciplinas metodológicas, ligadas por padrão: normalizar ‖H‖ ao
     comparar topologias (senão "mais coerência" pode ser só "hopping maior") e
     RELIGAR em vez de acrescentar (senão o resultado mede |E|, não topologia). *)
  "hamiltonian" -> <|"kind" -> "adjacency", "gamma" -> 1.0,
    "normalization" -> "spectral", "lanczos_steps" -> 40|>,
  "initial" -> <|"site" -> 0|>,
  "time" -> <|"t1" -> 50.0, "nt" -> 500|>,
  "observables" -> <|"target" -> -1, "population" -> True,
    "module_concurrence" -> True, "full_concurrence" -> False, "pop_stride" -> 1|>,
  "realizations" -> 1, "escala" -> 1.0|>;

DaedalusEspecLer::chave = "Chave desconhecida em `1`: `2`. O parser do núcleo é \
estrito de propósito — um campo escrito errado produziria um arquivo \
internamente consistente resolvendo o problema errado.";

funde[padrao_Association, lido_Association, onde_String] := Module[{saida = padrao},
  Do[
    If[! KeyExistsQ[padrao, k],
      Message[DaedalusEspecLer::chave, onde, k]; Abort[]];
    saida[k] = If[AssociationQ[padrao[k]] && AssociationQ[lido[k]],
      funde[padrao[k], lido[k], onde <> "." <> k], lido[k]],
    {k, Keys[lido]}];
  saida];

(* Aceita caminho de arquivo ou o JSON em si. O segundo caso existe porque a
   interface mostra o spec canônico numa linha só, e copiar essa linha para
   dentro do caderno é o caminho mais curto entre "vi na tela" e "rodei aqui". *)
DaedalusEspecLer[entrada_String] := Module[{lido},
  lido = If[StringStartsQ[StringTrim[entrada], "{"],
    ImportString[entrada, "RawJSON"],
    Import[entrada, "RawJSON"]];
  funde[DaedalusEspecPadrao[], lido, "raiz"]];

(* ==========================================================================
   Verificação contra a referência do núcleo

   Ordem: PRNG, depois digital do grafo, depois observáveis. A digital não
   tolera nada — ela é inteira. E é ela que separa "o gerador discordou" de "o
   propagador discordou", que têm causas e correções completamente diferentes.
   ========================================================================== *)

leCabecalhoCSV[texto_String] := Association[
  Cases[StringSplit[texto, "\n"],
    linha_ /; StringStartsQ[linha, "#! "] :>
      Module[{corpo = StringDrop[linha, 3], i},
        i = StringPosition[corpo, " ", 1];
        If[i === {}, corpo -> "",
          StringTake[corpo, i[[1, 1]] - 1] -> StringDrop[corpo, i[[1, 1]]]]]]];

(* O CSV vem em notação C ("1.15e-05"); em Wolfram o expoente é "*^". Sem esta
   conversão, ToExpression lê o "e" como símbolo indefinido e devolve uma
   expressão simbólica que se propaga silenciosamente por toda a comparação. *)
paraNumero[s_String] := If[s === "nan" || s === "",
  Missing["NotAvailable"],
  ToExpression[StringReplace[s,
    {"e+" -> "*^", "e-" -> "*^-", "E+" -> "*^", "E-" -> "*^-",
     "e" -> "*^", "E" -> "*^"}]]];

leTabelaCSV[texto_String] := Module[{linhas, iCab, cab, dados},
  linhas = StringSplit[texto, "\n"];
  iCab = Module[{pos = Position[linhas, _String?(StringStartsQ[#, "t,norm,"] &), 1]},
    If[pos === {}, 0, pos[[1, 1]]]];
  If[iCab === 0, Return[<||>]];
  cab = StringSplit[linhas[[iCab]], ","];
  dados = TakeWhile[Drop[linhas, iCab], StringContainsQ[#, ","] &];
  dados = Select[dados, ! StringStartsQ[#, "#"] &];
  (* Sem linhas de dados devolve <||> em vez de uma Association de colunas
     vazias: a segunda é uma estrutura malformada que atravessa a comparação
     inteira antes de dar erro, e num arquivo truncado o que se quer é a
     ausência declarada, não um objeto com a forma certa e nada dentro. *)
  If[dados === {}, Return[<||>]];
  AssociationThread[cab ->
    Transpose[Map[Function[l, paraNumero /@ StringSplit[l, ","]], dados]]]];

(* Bloco "# estado final" do CSV: as AMPLITUDES, que é onde a comparação tem
   sensibilidade uniforme. Uma população de 1e-23 e uma de 1e-3 têm a mesma
   incerteza absoluta na amplitude; em probabilidade, não. *)
leEstadoCSV[texto_String] := Module[{linhas, i},
  linhas = StringSplit[texto, "\n"];
  i = Position[linhas, "j,re,im", 1];
  If[i === {}, Return[{}]];
  Map[
    Function[l, Module[{c = StringSplit[l, ","]},
      paraNumero[c[[2]]] + I paraNumero[c[[3]]]]],
    TakeWhile[Drop[linhas, i[[1, 1]]],
      StringContainsQ[#, ","] && ! StringStartsQ[#, "#"] &]]];

DaedalusVerificarReferencia::vazio = "O caso `1` não tem linhas para comparar: \
uma comparação entre dois conjuntos vazios passa em qualquer tolerância.";

DaedalusVerificarReferencia::semestado = "O caso `1` não traz o bloco de estado \
final: gere o CSV de referência com --estado, senão a comparação de amplitude \
não roda e o veredito sai bom por ausência de teste.";

(* ==========================================================================
   COMO A CONCORDÂNCIA É MEDIDA, E POR QUE NÃO BASTA UM ERRO RELATIVO

   A tentação é comparar coluna a coluna com erro relativo e declarar uma
   tolerância única. Isso não funciona aqui, e o modo como falha é instrutivo.

   Em t = 0 o núcleo escreve p_alvo = 0 EXATO: a expansão de Chebyshev em
   \[Alpha] = 0 dá J_0(0) = 1 e J_k(0) = 0, ou seja, a identidade, sem
   arredondamento. A soma espectral daqui devolve V V^T \[Psi]_0, que é \[Psi]_0
   a menos de \[Epsilon] de máquina, e a probabilidade sai 3e-32. O erro
   relativo entre 0 e 3e-32 é 1 — o pior possível — para uma discordância de
   trinta e duas casas decimais.

   O mesmo acontece, mais de leve, em toda cauda exponencialmente pequena. E há
   uma lei por trás: se as duas amplitudes diferem de \[Delta] em valor
   ABSOLUTO, a probabilidade p = |\[Psi]|^2 difere relativamente de ~2\[Delta]/
   \[Sqrt]p. Medido nos casos de referência, com \[Delta] \[TildeTilde] 1e-15:

     p ~ 1e-4  (completo)    previsto 2e-13    observado 7.5e-13
     p ~ 1e-10 (sbm)         previsto 2e-10    observado 7.3e-11
     p ~ 1e-13 (fechado)     previsto 3e-09    observado 6.2e-10
     p ~ 1e-23 (hipercubo)   previsto 2e-04    observado 9.3e-06

   Cinco ordens de grandeza no erro relativo, todas explicadas pela MESMA
   concordância de ~1e-15 na amplitude. O erro relativo da probabilidade não
   está medindo discordância entre os programas; está medindo o quão pequena é
   a probabilidade.

   Então a verificação reporta TRÊS números, e nenhum deles sozinho:

     PiorAbsoluto   — max |a-b| na tabela de escalares. É o critério que pega
                      erro de física: se o núcleo dá 1e-23 e nós damos 1e-3, a
                      diferença absoluta é 1e-3 e nada disso passa.
     PiorRelativo   — max |a-b|/|a| RESTRITO a |a| >= PisoRelativo, que é a
                      faixa onde exigir algarismos significativos faz sentido.
     PiorAmplitude  — max |\[Psi]_núcleo - \[Psi]_nossa| no estado final. É a
                      grandeza de precisão uniforme, e é o número que de fato
                      caracteriza os dois programas.

   O piso NÃO afrouxa o teste: abaixo dele o critério absoluto continua valendo,
   e é ele que reprova qualquer discordância real. O piso só desliga uma
   pergunta que não tem resposta em dupla precisão.

   Os dois critérios se dividem o trabalho, e nenhum sobra: norm, ipr e coh_l1
   estão sempre muito acima do piso e são governados pelo relativo; p_alvo e as
   populações por módulo passam por caudas de 1e-30 e são governadas pelo
   absoluto.

   ---------------------------------------------------------------------------
   PIOR CASO OBSERVADO, MEDIDO (não estimado), em 12 casos de referência

     absoluto   9.4e-12  linha,  coh_l1 — que é O(n) e vale ~30 ali, ou seja,
                         3e-13 em termos relativos
     relativo   9.6e-11  microtubulo-seam0, pmod1 ~ 1.3e-9
     amplitude  1.5e-14  grade 2D (n = 240)

   A grade 2D é mesmo o pior em amplitude, como se espera de um grafo com
   espectro muito degenerado: \[Kappa](V) é maior e a soma espectral cancela
   mais. O piso de ruído da soma espectral escala com n:

     grade  n = 240 \[Rule] 1.5e-14        linha  n = 300 \[Rule] 1.3e-14
     lei aproximada: n \[Epsilon]/4, com \[Epsilon] = 2.2e-16

   As tolerâncias declaradas abaixo são uma ordem de grandeza acima do pior
   observado (a de amplitude, quase duas). Quem rodar N muito maior que 300 deve
   esperar o piso subir junto com n e reajustar a de amplitude — e reajustar
   sabendo por quê, não para fazer passar.

   ---------------------------------------------------------------------------
   ONDE OS DOIS MAIS DISCORDAM, E POR QUÊ ISSO É BOM

   O maior desvio de amplitude do conjunto todo está no caso linha, no sítio
   j = 250, a 250 saltos da origem em t curto. Ali o núcleo devolve 1.1e-68 e
   este pacote devolve ~1.3e-14.

   Nenhum dos dois está errado, e a diferença é estrutural: Chebyshev constrói o
   estado por multiplicações esparsas sucessivas, então um sítio a 250 saltos só
   é alcançado depois de 250 termos e sai com a cauda exponencial correta. A
   soma espectral é uma soma de n termos de módulo O(1) que se cancelam, e o que
   sobra do cancelamento é \[Epsilon] de máquina.

   Isto é a evidência funcionando como projetada: os dois métodos concordam onde
   há física e discordam onde só há aritmética, e a discordância aparece na
   estrutura de erro esperada de cada um. Se concordassem também em 1e-68, o
   motivo mais provável seria estarem rodando o mesmo algoritmo.
   ========================================================================== *)

daeCompararCaso[nome_String, csvRef_String, espec0_Association,
                tolAbs_, tolRel_, piso_, tolAmp_] := Module[
           {espec = espec0, cabRef, tabRef, estRef, r, csvNosso, tabNossa, estNossa,
            digitalRef, falhasRef, linhasRef, colunas, escalaNucleo,
            piorAbs = 0., ondeAbs = "", piorRel = 0., ondeRel = "",
            piorAmp = 0., ondeAmp = "", comparados = 0},
      cabRef = leCabecalhoCSV[csvRef];
      tabRef = leTabelaCSV[csvRef];
      estRef = leEstadoCSV[csvRef];
      espec["core_hash"] = cabRef["core_hash"];

      (* A escala do núcleo é IMPOSTA aqui. Ela é estimativa de Lanczos, a nossa
         é o raio espectral exato, e as duas discordam por construção. Sem impor,
         o H seria outro e TODAS as colunas divergiriam por um motivo só — o que
         esconderia qualquer discordância de propagador atrás de uma de norma.
         A diferença entre as duas escalas sai reportada logo abaixo. *)
      escalaNucleo = paraNumero[cabRef["scale"]];
      r = DaedalusRodar[espec, "EscalaImposta" -> escalaNucleo];
      csvNosso = DaedalusCSV[r["Espec"], r["Rede"], r["Metricas"], r["Propagacao"], r["Observaveis"]];
      tabNossa = leTabelaCSV[csvNosso];
      estNossa = Last[r["Propagacao"]["Estados"]];

      (* ANTI-VACUIDADE: sem linhas, tudo bate. *)
      linhasRef = If[KeyExistsQ[tabRef, "t"], Length[tabRef["t"]], 0];
      If[linhasRef < 2,
        Message[DaedalusVerificarReferencia::vazio, nome];
        Return[<|"Caso" -> nome, "Veredito" -> "VAZIO"|>, Module]];
      If[Length[estRef] =!= Length[estNossa],
        Message[DaedalusVerificarReferencia::semestado, nome];
        Return[<|"Caso" -> nome, "Veredito" -> "SEM ESTADO"|>, Module]];

      (* 2. DIGITAL, antes de qualquer observável e sem tolerância nenhuma. *)
      digitalRef = ToExpression[cabRef["graph_fingerprint"]];
      If[digitalRef =!= r["Rede"]["Digital"],
        Return[<|"Caso" -> nome, "Veredito" -> "DIGITAL DIVERGIU",
          "DigitalRef" -> digitalRef, "DigitalNossa" -> r["Rede"]["Digital"],
          "Nota" -> "os dois construíram grafos diferentes: é o gerador ou o \
PRNG, não o propagador"|>, Module]];

      (* 3. RELIGAÇÕES QUE FALHARAM. A digital não vê isto: uma tentativa que
            colide consome 100 sorteios e não muda aresta nenhuma. Em K_N toda
            tentativa colide, o grafo sai idêntico ao original, e este contador
            é a única testemunha de que os dois gastaram o mesmo fluxo. *)
      falhasRef = ToExpression[Lookup[cabRef, "rewire_failed", "0"]];
      If[falhasRef =!= r["Rede"]["ReligacoesFalhas"],
        Return[<|"Caso" -> nome, "Veredito" -> "RELIGAÇÕES FALHAS DIVERGIRAM",
          "FalhasRef" -> falhasRef, "FalhasNossa" -> r["Rede"]["ReligacoesFalhas"],
          "Nota" -> "mesmo grafo, fluxos de sorteio diferentes: confira o limite \
de 100 tentativas e a regra de b == i consumir tentativa"|>, Module]];

      (* 4. Observáveis, só as colunas presentes nos dois lados. *)
      colunas = Intersection[Keys[tabRef], Keys[tabNossa]];
      Do[
        Module[{a = tabRef[c], b = tabNossa[c]},
          Do[
            Module[{va = a[[k]], vb = b[[k]], d, onde},
              onde = c <> " linha " <> ToString[k] <> ": " <>
                ToString[CForm[va]] <> " contra " <> ToString[CForm[vb]];
              Which[
                MissingQ[va] && MissingQ[vb], Null,
                (* nan de um lado só é discordância de ESTRUTURA, não de número:
                   p_alvo ausente e p_alvo igual a zero são coisas diferentes. *)
                MissingQ[va] || MissingQ[vb],
                  piorAbs = Infinity; ondeAbs = c <> " linha " <> ToString[k] <>
                    ": um lado é nan e o outro não",
                True,
                  ++comparados;
                  d = Abs[va - vb];
                  If[d > piorAbs, piorAbs = d; ondeAbs = onde];
                  If[Abs[va] >= piso,
                    d = Abs[va - vb]/Abs[va];
                    If[d > piorRel, piorRel = d; ondeRel = onde]]]],
            {k, Min[Length[a], Length[b]]}]],
        {c, colunas}];

      (* 5. Amplitudes do estado final: precisão uniforme, sem o piso. *)
      Do[
        Module[{d = Abs[estRef[[j]] - estNossa[[j]]]},
          If[d > piorAmp, piorAmp = d;
            ondeAmp = "j=" <> ToString[j - 1] <> ": " <>
              ToString[CForm[estRef[[j]]]] <> " contra " <> ToString[CForm[estNossa[[j]]]]]],
        {j, Length[estRef]}];

      (* ANTI-VACUIDADE: nenhuma célula comparada também "passa". *)
      If[comparados < linhasRef,
        Return[<|"Caso" -> nome, "Veredito" -> "NADA COMPARADO",
          "Celulas" -> comparados|>, Module]];

      <|"Caso" -> nome,
        "Veredito" -> If[piorAbs <= tolAbs && piorRel <= tolRel && piorAmp <= tolAmp,
          "OK", "FORA DA TOLERÂNCIA"],
        "N" -> r["Rede"]["N"], "Linhas" -> linhasRef, "Celulas" -> comparados,
        "Metodo" -> r["Propagacao"]["Metodo"],
        "PiorAbsoluto" -> piorAbs, "OndeAbsoluto" -> ondeAbs,
        "PiorRelativo" -> piorRel, "OndeRelativo" -> ondeRel,
        "PiorAmplitude" -> piorAmp, "OndeAmplitude" -> ondeAmp,
        "EscalaNucleo" -> escalaNucleo,
        "EscalaExataNossa" -> r["Hamiltoniano"]["EscalaExata"],
        "Lambda2Nucleo" -> paraNumero[cabRef["lambda2"]],
        "Lambda2Nossa" -> r["Metricas"]["Lambda2"],
        "ReligacoesFalhas" -> r["Rede"]["ReligacoesFalhas"]|>];

Options[DaedalusVerificarReferencia] = {
  "TolAbsoluta"  -> 1.*^-10,   (* pior observado 9.4e-12, margem 10x  *)
  "TolRelativa"  -> 1.*^-9,    (* pior observado 9.6e-11, margem 10x  *)
  "PisoRelativo" -> 1.*^-10,
  "TolAmplitude" -> 1.*^-12,   (* pior observado 1.5e-14, margem 68x  *)
  "Verboso" -> True};

DaedalusVerificarReferencia[pasta_String, OptionsPattern[]] := Module[
  {tabPRNG, prngOK, casos, relatorio,
   tolAbs = OptionValue["TolAbsoluta"], tolRel = OptionValue["TolRelativa"],
   piso = OptionValue["PisoRelativo"], tolAmp = OptionValue["TolAmplitude"],
   verboso = OptionValue["Verboso"]},

  (* 1. PRNG. Verificação barata e diagnóstica: se ele divergir, a falha aparece
        AQUI, e não como digital divergente lá na frente, onde não se sabe se a
        culpa é do gerador, da ordem dos sorteios ou da religação. *)
  tabPRNG = Import[FileNameJoin[{pasta, "prng.json"}], "RawJSON"];
  prngOK = AllTrue[tabPRNG["fluxos"],
    DaedalusFluxoAleatorio[#["semente"], Length[#["u64"]]] === (ToExpression /@ #["u64"]) &];
  If[verboso, Print["PRNG contra a tabela do núcleo: ",
    If[prngOK, "IDÊNTICO", "DIVERGIU — pare aqui, nada abaixo faz sentido"]]];

  casos = Select[FileNames["*.json", pasta], ! StringEndsQ[#, "prng.json"] &];

  relatorio = Table[
    daeCompararCaso[FileBaseName[caso],
      Import[StringDrop[caso, -5] <> ".csv", "Text"],
      DaedalusEspecLer[caso], tolAbs, tolRel, piso, tolAmp],
    {caso, casos}];

  If[verboso, Print[Dataset[relatorio]]];
  <|"PRNG" -> prngOK, "Casos" -> relatorio,
    "Tolerancias" -> <|"Absoluta" -> tolAbs, "Relativa" -> tolRel,
      "PisoRelativo" -> piso, "Amplitude" -> tolAmp|>,
    "PiorAbsoluto" -> Max[Cases[relatorio, a_Association :> Lookup[a, "PiorAbsoluto", 0.]]],
    "PiorRelativo" -> Max[Cases[relatorio, a_Association :> Lookup[a, "PiorRelativo", 0.]]],
    "PiorAmplitude" -> Max[Cases[relatorio, a_Association :> Lookup[a, "PiorAmplitude", 0.]]],
    "Veredito" -> If[prngOK && AllTrue[relatorio, Lookup[#, "Veredito", ""] === "OK" &],
      "TUDO OK", "HÁ CASOS FORA"]|>];


(* ==========================================================================
   LINDBLAD EXATO, POR SUPEROPERADOR VETORIZADO

   Defasagem pura, L_j = sqrt(gamma)|j><j|. Vetorizando por LINHAS, que é como
   Flatten empilha em Wolfram, vec(A rho B) = (A ⊗ B^T) vec(rho), e daí

     L = -i (H ⊗ I - I ⊗ H^T) + gamma (Sum_j P_j ⊗ P_j) - gamma I

   O último par de termos é o dissipador: Sum_j P_j rho P_j mantém só a diagonal,
   e -gamma/2 {Sum_j L_j^dag L_j, rho} = -gamma rho porque Sum_j L_j^dag L_j =
   gamma I. Ou seja, as coerências decaem a taxa gamma e as populações não.

   Podia-se pular tudo isso e escrever direto "off-diagonal decai a exp(-gamma t)".
   NÃO se pula de propósito: o atalho analítico usa o mesmo raciocínio que se
   quer conferir, e duas contas com o mesmo raciocínio concordam mesmo quando o
   raciocínio está errado. Montar o superoperador é o caminho independente.
   ========================================================================== *)

Options[DaedalusLindblad] = {"EscalaImposta" -> Automatic, "Amostras" -> Automatic};

DaedalusLindblad[espec_Association, OptionsPattern[]] := Module[
  {rede, ham, n, H, Id, proj, L, rho0, dt, stride, nt, tempos, vec},
  rede = DaedalusRede[espec];
  ham = DaedalusHamiltoniano[rede, espec,
          "EscalaImposta" -> OptionValue["EscalaImposta"]];
  H = N[Normal[ham["H"]]];
  n = Length[H];
  Id = IdentityMatrix[n, SparseArray];

  (* Sum_j P_j ⊗ P_j: 1 em cada posicao ((j,j),(j,j)) da matriz n^2 x n^2 *)
  proj = SparseArray[
    Table[{(j - 1) n + j, (j - 1) n + j} -> 1., {j, n}], {n^2, n^2}];

  L = -I (KroneckerProduct[H, Id] - KroneckerProduct[Id, Transpose[H]]) +
      espec["trajectories"]["gamma_deph"] (proj - IdentityMatrix[n^2, SparseArray]);

  rho0 = ConstantArray[0., {n, n}];
  rho0[[espec["initial"]["site"] + 1, espec["initial"]["site"] + 1]] = 1.;
  vec = Flatten[rho0];

  nt = espec["time"]["nt"];
  dt = espec["time"]["t1"]/nt;
  stride = espec["trajectories"]["rho_stride"];
  tempos = Table[k stride dt, {k, 0, Floor[nt/stride]}];

  <|"Tempos" -> tempos,
    "Rho" -> Table[ArrayReshape[MatrixExp[L t, vec], {n, n}], {t, tempos}],
    "Escala" -> ham["Escala"]|>];

(* Le o rho medio que `daedalus traj` escreveu. *)
leRhoCSV[texto_String] := Module[{linhas, cab, n, na, dados, rho},
  cab = leCabecalhoCSV[texto];
  n = ToExpression[cab["n"]];
  na = ToExpression[cab["n_amostras"]];
  linhas = StringSplit[texto, "\n"];
  dados = Select[Drop[linhas, Position[linhas, "s,t,i,j,re,im", 1][[1, 1]]],
    StringContainsQ[#, ","] &];
  rho = ConstantArray[0. + 0. I, {na, n, n}];
  Do[Module[{c = StringSplit[l, ","]},
      rho[[ToExpression[c[[1]]] + 1, ToExpression[c[[3]]] + 1,
           ToExpression[c[[4]]] + 1]] = paraNumero[c[[5]]] + I paraNumero[c[[6]]]],
    {l, dados}];
  <|"Cabecalho" -> cab, "N" -> n, "NAmostras" -> na, "Rho" -> rho,
    "NTraj" -> ToExpression[cab["n_traj"]],
    "Escala" -> paraNumero[cab["scale"]]|>];

DaedalusVerificarTrajetorias::poucos = "Sao precisos ao menos tres arquivos com n_traj diferentes: a LEI de convergencia so aparece com tres pontos, e um desvio isolado nao distingue \"convergiu\" de \"convergiu para outro lugar\".";

Options[DaedalusVerificarTrajetorias] = {"Verboso" -> True, "GammaSabotado" -> Automatic};

DaedalusVerificarTrajetorias[arquivos_List, OptionsPattern[]] := Module[
  {lidos, espec, exato, desvios, maximos, ns, razoes, esperado, expoente, veredito,
   verboso = OptionValue["Verboso"]},
  If[Length[arquivos] < 3, Message[DaedalusVerificarTrajetorias::poucos]; Return[$Failed]];

  lidos = leRhoCSV[Import[#, "Text"]] & /@ arquivos;
  lidos = SortBy[lidos, #["NTraj"] &];
  espec = DaedalusEspecLer[lidos[[1]]["Cabecalho"]["spec"]];

  (* A escala do NUCLEO e imposta: ela e estimativa de Lanczos, a nossa seria o
     raio exato, e sem impor os dois compararao dinamicas de H diferentes — a
     discordancia apareceria como desdobramento errado. *)
  (* Sabotagem opcional: resolve Lindblad com OUTRA taxa de defasagem. As
     trajetorias entao convergem para um alvo que nao e este, e o desvio tem de
     ESTABILIZAR em vez de cair. E a companheira que separa "o teste detecta
     convergencia" de "o teste detecta qualquer coisa". *)
  If[NumericQ[OptionValue["GammaSabotado"]],
    espec["trajectories"]["gamma_deph"] = OptionValue["GammaSabotado"]];

  exato = DaedalusLindblad[espec, "EscalaImposta" -> lidos[[1]]["Escala"]];

  (* RMS sobre as entradas, nao o maximo. O maximo sobre 2500 entradas e
     estatistica de EXTREMO: ele tambem cai como 1/sqrt(n), mas com dispersao
     grande, e um portao construido sobre ele reprovaria por sorte. O RMS mede a
     mesma lei com muito menos ruido. O maximo continua reportado, porque e ele
     que responde "qual o pior erro em qualquer entrada". *)
  desvios = Table[
    Sqrt@Mean@Flatten@Table[
      Abs[r["Rho"][[k]] - exato["Rho"][[k]]]^2, {k, r["NAmostras"]}],
    {r, lidos}];
  maximos = Table[
    Max@Table[Max@Abs@Flatten[r["Rho"][[k]] - exato["Rho"][[k]]], {k, r["NAmostras"]}],
    {r, lidos}];
  ns = #["NTraj"] & /@ lidos;

  (* A LEI, e ela e ajustada sobre TODOS os pontos: desvio ~ n^(-p), com p = 1/2
     para Monte Carlo nao enviesado. Uma razao entre dois pontos consecutivos tem
     variancia propria, e um limiar sobre ela reprova o caso CORRETO por sorte
     (1.29 contra um limiar de 1.3, por exemplo). O expoente usa a faixa inteira
     e separa os dois casos por duas ordens de grandeza: 0.49 quando converge,
     0.006 quando o alvo esta errado. *)
  razoes = Table[desvios[[i]]/desvios[[i + 1]], {i, Length[desvios] - 1}];
  esperado = Table[Sqrt[ns[[i + 1]]/ns[[i]]], {i, Length[ns] - 1}];
  expoente = -Coefficient[
    Fit[Transpose[{Log[N@ns], Log[desvios]}], {1, \[FormalX]}, \[FormalX]], \[FormalX]];

  If[verboso,
    Print["Lindblad exato contra media de trajetorias do nucleo em C"];
    Do[Print["  n_traj = ", StringPadLeft[ToString[ns[[i]]], 6],
             "   rms = ", ScientificForm[desvios[[i]], 3],
             "   maximo = ", ScientificForm[maximos[[i]], 3]],
       {i, Length[ns]}];
    Do[Print["  ", ns[[i]], " -> ", ns[[i + 1]], ": caiu ",
             NumberForm[razoes[[i]], 3], "x  (a lei 1/sqrt(n) preve ",
             NumberForm[esperado[[i]], 3], "x)"], {i, Length[razoes]}]];

  (* LIMIAR CALIBRADO PELA MEDIDA, nao pela teoria. Um limiar de 0.35 — "0.5
     menos uma folga" — da FALSO ALARME em microtubulo g = 0.09, que mede
     0.331. Medido em 10 casos (2 topologias x 5 valores de gamma), com a
     companheira de gamma dobrado em cada um:

       pior caso CORRETO   p = 0.331      (microtubulo, gamma = 0.09)
       pior caso SABOTADO  p = 0.017      (microtubulo, gamma = 0.09)
       separacao                 20x

     0.15 fica 2.2x abaixo do pior correto e 8.8x acima do pior sabotado. Com
     quatro pontos o expoente tem variancia propria, e um limiar encostado no
     valor teorico reprova por sorte, como as razoes consecutivas.
     Ver CONVENTIONS.md 10.1.7. *)
  veredito = If[expoente > 0.15,
    "CONVERGE para Lindblad", "NAO CONVERGE — o desvio nao cai com n_traj"];
  If[verboso,
    Print["  expoente ajustado: desvio ~ n^(-", NumberForm[expoente, 3],
          "), a lei preve n^(-0.5)"];
    Print["  ", veredito]];

  <|"NTraj" -> ns, "Desvios" -> desvios, "Maximos" -> maximos, "Expoente" -> expoente,
    "Razoes" -> razoes,
    "RazoesEsperadas" -> esperado, "Veredito" -> veredito,
    "Exato" -> exato|>];

(* ==========================================================================
   ANTI-VACUIDADE

   Toda afirmação de invariância precisa da companheira que prova que a
   invariância PODERIA ter quebrado. "Os 12 casos bateram" só quer dizer alguma
   coisa depois que se mostra que um 13º, sabotado, não bate — e que a mensagem
   aponta para a sabotagem certa, porque "o gerador divergiu" e "o propagador
   divergiu" se corrigem em lugares diferentes.
   ========================================================================== *)

(* Perturba a coluna coh_l1 da primeira linha de dados por uma fração relativa.
   coh_l1 porque é O(n) e está muito acima do piso: a perturbação atravessa os
   dois critérios, o absoluto e o relativo. *)
daePerturbar[texto_String, deltaRel_] := Module[{linhas, i, c},
  linhas = StringSplit[texto, "\n"];
  i = Position[linhas, _String?(StringStartsQ[#, "t,norm,"] &), 1][[1, 1]];
  c = StringSplit[linhas[[i + 1]], ","];
  c[[4]] = ToString[CForm[paraNumero[c[[4]]] (1 + deltaRel)]];
  linhas[[i + 1]] = StringRiffle[c, ","];
  StringRiffle[linhas, "\n"]];

daeCortarApos[texto_String, marca_String] := Module[{linhas, i},
  linhas = StringSplit[texto, "\n"];
  i = Position[linhas, _String?(StringStartsQ[#, marca] &), 1];
  If[i === {}, texto, StringRiffle[Take[linhas, i[[1, 1]]], "\n"]]];

DaedalusAntiVacuidade[pasta_String] := Module[
  {arq, csv, espec, sab, resultados, ident},

  sab[nome_, csvS_, especS_] := Quiet[
    Lookup[daeCompararCaso[nome, csvS, especS, 1.*^-10, 1.*^-9, 1.*^-10, 1.*^-12],
      "Veredito", "SEM VEREDITO"],
    {DaedalusVerificarReferencia::vazio, DaedalusVerificarReferencia::semestado}];

  resultados = {};

  (* 1. GRAFO DIFERENTE. Muda o seam_shift e compara contra o CSV do original:
        o passo da costura muda a topologia sem mudar |V| nem |E|, que é
        exatamente o tipo de divergência que passaria despercebida em qualquer
        conferência de tamanho. *)
  arq = FileNameJoin[{pasta, "microtubulo-seam3.json"}];
  csv = Import[FileNameJoin[{pasta, "microtubulo-seam3.csv"}], "Text"];
  espec = DaedalusEspecLer[arq];
  AppendTo[resultados, <|
    "Sabotagem" -> "seam_shift 3 -> 2",
    "Esperado" -> "DIGITAL DIVERGIU",
    "Obtido" -> sab["sabotado", csv,
      Append[espec, "graph" -> Append[espec["graph"],
        "params" -> Append[espec["graph"]["params"], "seam_shift" -> 2]]]]|>];

  (* 2. FLUXO DE SORTEIOS DIFERENTE, MESMO GRAFO. Em K_12 toda tentativa de
        religação colide, o grafo sai idêntico e a digital não vê nada. Só o
        contador de falhas testemunha os 3000 sorteios gastos. *)
  arq = FileNameJoin[{pasta, "completo-religado.json"}];
  csv = Import[FileNameJoin[{pasta, "completo-religado.csv"}], "Text"];
  AppendTo[resultados, <|
    "Sabotagem" -> "rewire_failed 30 -> 29, digital intacta",
    "Esperado" -> "RELIGAÇÕES FALHAS DIVERGIRAM",
    "Obtido" -> sab["sabotado",
      StringReplace[csv, "#! rewire_failed 30" -> "#! rewire_failed 29"],
      DaedalusEspecLer[arq]]|>];

  (* 3 e 4. O LIMIAR TEM DOIS LADOS. Uma tolerância que reprova tudo é tão
     inútil quanto uma que aprova tudo. *)
  arq = FileNameJoin[{pasta, "linha.json"}];
  csv = Import[FileNameJoin[{pasta, "linha.csv"}], "Text"];
  espec = DaedalusEspecLer[arq];
  AppendTo[resultados, <|
    "Sabotagem" -> "coh_l1 deslocada de 1e-6 relativo",
    "Esperado" -> "FORA DA TOLERÂNCIA",
    "Obtido" -> sab["sabotado", daePerturbar[csv, 1.*^-6], espec]|>];
  AppendTo[resultados, <|
    "Sabotagem" -> "coh_l1 deslocada de 1e-15 relativo",
    "Esperado" -> "OK",
    "Obtido" -> sab["sabotado", daePerturbar[csv, 1.*^-15], espec]|>];

  (* 5. ARQUIVO TRUNCADO. Sem linhas, toda comparação passa: é o modo de falha
        mais perigoso, porque o resultado é verde. *)
  AppendTo[resultados, <|
    "Sabotagem" -> "CSV cortado no cabeçalho da tabela",
    "Esperado" -> "VAZIO",
    "Obtido" -> sab["sabotado", daeCortarApos[csv, "t,norm,"], espec]|>];

  (* 6. SEM O BLOCO DE ESTADO. A comparação de amplitude é a mais sensível das
        três; se ela sumir em silêncio, o veredito melhora sem que nada tenha
        melhorado. *)
  AppendTo[resultados, <|
    "Sabotagem" -> "CSV sem o bloco de estado final",
    "Esperado" -> "SEM ESTADO",
    "Obtido" -> sab["sabotado", daeCortarApos[csv, "# estado final"], espec]|>];

  (* 7. IDENTIDADE DA CONCURRENCE. C_MM = s_M^2 - q_M e C_MN = 2 s_M s_N somam,
        sobre M <= N, exatamente (Σ|ψ_j|)^2 - 1 = coerência ℓ1. Isto não é
        comparação com o núcleo: é a álgebra da própria convenção, e se ela não
        fechar, a decomposição por módulo não é decomposição de coisa nenhuma. *)
  ident = Module[{r, obs, nmod, somas, coh},
    r = DaedalusRodar[DaedalusEspecLer[FileNameJoin[{pasta, "sbm.json"}]]];
    obs = r["Observaveis"]; nmod = r["Rede"]["NModulos"];
    somas = Table[
      Total[Flatten[Table[If[a <= b, c[[a, b]], 0.], {a, nmod}, {b, nmod}]]],
      {c, obs["ConcurrencePorModulo"]}];
    coh = obs["CoerenciaL1"];
    Max[Abs[somas - coh]/Map[Max[Abs[#], 1.] &, coh]]];
  AppendTo[resultados, <|
    "Sabotagem" -> "identidade Σ_{M<=N} C_MN = C_ℓ1 (álgebra, não comparação)",
    "Esperado" -> "OK",
    "Obtido" -> If[ident <= 1.*^-12, "OK",
      "QUEBROU em " <> ToString[CForm[ident]]]|>];

  resultados = (Append[#, "Veredito" ->
    If[#["Obtido"] === #["Esperado"], "detectou", "NÃO DETECTOU"]] & /@ resultados);

  <|"Companheiras" -> resultados,
    "DesvioIdentidade" -> ident,
    "Veredito" -> If[AllTrue[resultados, #["Veredito"] === "detectou" &],
      "a verificação pode falhar, e falha pelo motivo certo",
      "HÁ SABOTAGEM NÃO DETECTADA — a verificação está parcialmente vazia"]|>];


End[];
EndPackage[];
