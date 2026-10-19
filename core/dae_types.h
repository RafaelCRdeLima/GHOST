/* dae_types.h — tipos base do núcleo do Daedalus.
 *
 * REGRAS DO NÚCLEO (valem para todo arquivo em core/):
 *   1. C99 que também é C++17 válido. O mesmo fonte compila com gcc -std=c99 e
 *      com g++ -std=c++17, porque o .cpp exportado é a amalgamação deste
 *      diretório. Sem `_Complex`, sem VLA, sem inicializador designado, sem
 *      literal composto, todo malloc explicitamente convertido.
 *   2. Zero estado global mutável. O template exportado roda realizações sob
 *      `#pragma omp parallel for`; qualquer estado compartilhado seria corrida.
 *   3. Sem <stdio.h>. Impressão só em native/ e nos templates.
 *   4. Nunca compilar com -ffast-math: quebra o NaN de p_alvo e a
 *      reprodutibilidade bit a bit.
 *
 * Comentários em português, identificadores em inglês.
 */
#ifndef DAE_TYPES_H
#define DAE_TYPES_H

#include <stddef.h>
#include <stdint.h>

/* -ffast-math é PROIBIDO, e a proibição é verificada aqui em vez de depender de
 * um teste. `t91_conventions` pega o sintoma — p_alvo deixa de ser NaN — mas
 * esse canário depende de um observável opcional continuar existindo com essa
 * forma. Se alguém trocar p_alvo por -1.0 para a interface ter algo plotável, o
 * portão sumiria junto e ninguém perceberia. Isto não some.
 *
 * -ffinite-math-only entra separado porque sozinho ele já colapsa o NaN sem
 * definir __FAST_MATH__. */
#if defined(__FAST_MATH__)
#error "Daedalus: -ffast-math (ou -Ofast) esta ligado. Ele quebra o NaN de p_alvo e a reprodutibilidade bit a bit entre alvos. Ver CONVENTIONS.md, parte 1."
#endif
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__
#error "Daedalus: -ffinite-math-only esta ligado. Ele colapsa o NaN de p_alvo, que e como a ausencia de sitio-alvo se distingue de p_alvo = 0. Ver CONVENTIONS.md, parte 7.2."
#endif

/* Versão do formato do spec.json. Muda quando o significado de um campo muda.
 *
 * 2 — bloco "trajectories", com "output_mode" obrigatório quando o
 *     bloco existe. Subiu porque a regra de dae_spec.h manda subir em campo
 *     novo, e o custo é justamente o ponto: um arquivo da versão 1 não
 *     descreve o formato 2, e recusá-lo é preferível a interpretá-lo por
 *     omissão. */
/* 3 — bloco "neutrino" do GHOST: PMNS, massas, energia, materia. Com ele,
 *     `graph` e `hamiltonian` sao ignorados (a fisica monta H) e `time.t1`
 *     passa a ser lido em quilometros.
 *
 * FAIXA DE VERSOES, e ela e deliberada: o nucleo aceita a FAIXA
 * [MIN, DAE_FORMAT_VERSION], e nao uma versao so por igualdade, e cada recurso
 * novo checa a sua propria versao minima (o bloco neutrino exige 3, em
 * dae_spec.c). O motivo e concreto: subir a versao por causa de um bloco novo
 * nao deve invalidar os specs existentes, e o repositorio tem dezenas deles em
 * specs/. A regra de dae_spec.h continua valendo — campo novo SOBE a versao;
 * arquivo antigo continua legivel, mas campo novo nao sai de graca. */
/* 4 — blocos "magnetic" e "trajectory" do setor magnético: conversão de
 *     helicidade ao longo de um raio já integrado, com o campo projetado na
 *     base transversal transportada. Exclusivos com o bloco "neutrino". */
#define DAE_FORMAT_VERSION     4
#define DAE_FORMAT_VERSION_MIN 2

/* Teto da concurrence par a par completa: N² entradas em memória. */
#define DAE_FULL_CONC_MAX 2000

/* Tolerância da rede de segurança da norma em dae_cheb_step. */
#define DAE_NORM_TOL 1e-10

typedef enum {
  DAE_OK = 0,
  DAE_ERR_ALLOC,          /* falta de memória                                */
  DAE_ERR_PARAM,          /* argumento inválido                              */
  DAE_ERR_JSON,           /* spec.json malformado (ver dae_error)            */
  DAE_ERR_VERSION,        /* format_version incompatível                     */
  DAE_ERR_TOOBIG,         /* pedido acima de um teto explícito               */
  DAE_ERR_CANCELLED,      /* callback de progresso pediu parada              */
  DAE_ERR_NORM,           /* norma deixou de se conservar — resultado é lixo  */
  DAE_ERR_UNIMPLEMENTED   /* caminho previsto na API, ainda não escrito      */
} dae_status;

const char *dae_strerror(dae_status s);

/* Matriz real simétrica em CSR.
 *
 * H é real (adjacência ou laplaciana com pesos reais) e só o estado é
 * complexo — por isso o SpMV aplica a MESMA matriz a re[] e a im[], o que
 * corta pela metade o tráfego de memória do laço quente.
 *
 * Dentro de cada linha, colind é estritamente crescente (ordenado e sem
 * duplicatas). dae_csr_from_edges garante essa invariante.
 */
typedef struct {
  int32_t  n;        /* ordem                                    */
  int32_t  nnz;      /* entradas armazenadas                     */
  int32_t *rowptr;   /* n+1                                      */
  int32_t *colind;   /* nnz                                      */
  double  *val;      /* nnz                                      */
} dae_csr;

/* Estado como dois arrays reais separados (structure of arrays). Nunca
 * intercalado: o SIMD do WASM só rende neste layout. */
typedef struct {
  int32_t  n;
  double  *re;
  double  *im;
} dae_vec;

/* Diagnóstico de erro textual — usado pelo parser do spec e exibido na
 * interface. Uma única fonte de validação, do lado C. */
typedef struct {
  int32_t line;
  int32_t col;
  char    msg[192];
} dae_error;

#endif /* DAE_TYPES_H */
