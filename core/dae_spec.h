/* dae_spec.h — o spec.json, interpretado SÓ aqui.
 *
 * O TypeScript monta o texto; quem interpreta é este arquivo. O .cpp exportado
 * embute o mesmo texto num raw string literal e usa este mesmo parser. Não há
 * validação paralela em JavaScript: a interface chama isto via WASM e exibe
 * dae_error.line/col. Ver CONVENTIONS.md, parte 2.
 *
 * ESTRITO DE PROPÓSITO: chave desconhecida é ERRO, não campo ignorado. Um
 * emissor que escreva "seam" onde o formato diz "seam_shift" produziria, com um
 * parser tolerante, um arquivo internamente consistente que resolve o problema
 * ERRADO — e o oráculo de verificação não pegaria, porque ele valida o
 * propagador, não o emissor. Com parser estrito, o erro aparece na linha e na
 * coluna. Campo novo exige subir format_version.
 */
#ifndef DAE_SPEC_H
#define DAE_SPEC_H

#include "dae_graph.h"
#include "dae_neutrino.h"
#include "dae_numag.h"
#include "dae_version.generated.h"
#include "dae_ham.h"
#include "dae_obs.h"
#include "dae_traj.h"
#include "dae_types.h"

typedef struct {
  int32_t        format_version;
  char           core_hash[24];      /* "" se ausente no arquivo */
  uint64_t       seed;

  dae_gen_params gen;                /* modo procedimental                    */
  int32_t        n_edges;            /* modo explícito: 0 = procedimental     */
  int32_t       *ei, *ej;
  double        *w;

  dae_ham_kind   ham;
  double         gamma;
  dae_norm_kind  norm;
  int32_t        lanczos_steps;

  int32_t        init_site;          /* -1 quando há vetor explícito          */
  double        *init_re, *init_im;
  int32_t        init_n;

  double         t1;
  int32_t        nt;

  int32_t        target;
  int            want_pop, want_conc_mod, want_conc_full;
  int32_t        pop_stride;
  int32_t        realizations;       /* varredura de semente; 1 = execução única */

  /* Trajetórias. `n_traj == 0` significa execução unitária, sem trajetórias,
     e o bloco "trajectories" pode estar ausente. Se o bloco EXISTE,
     `output_mode` é obrigatório: ver dae_traj.h para por que não há padrão
     implícito. */
  int32_t        n_traj;
  double         gamma_deph;
  int32_t        rho_stride;
  dae_saida_traj saida_traj;

  /* GHOST: oscilacao de sabor. `nu_ativo == 0` significa que o bloco
     "neutrino" nao existe e o programa se comporta como o nucleo generico —
     specs de format_version 2 continuam valendo, sem mudanca.
     Com o bloco presente:
       - `graph` e `hamiltonian` sao IGNORADOS: quem monta H e a fisica;
       - `initial.site` e IGNORADO: o sabor inicial vem de `neutrino.flavor`;
       - `time.t1` passa a ser lido em QUILOMETROS, e convertido por
         DAE_NU_T_POR_KM. E o unico campo que muda de unidade, e muda porque
         linha de base em km e o que um experimento reporta.
     Exige format_version >= 3. */
  int            nu_ativo;
  dae_nu_params  nu;
  double         nu_rho;             /* g/cm3; atalho para perfil constante   */
  double         nu_ye;              /* fracao eletronica                     */
  /* Perfil de materia. Sem o bloco `matter`, vale `nu_rho` constante. Com ele,
     a densidade varia ao longo do caminho e a propagacao passa a ser por
     fatias com controle de passo adaptativo (dae_nu_propaga). */
  int            nu_perfil_ativo;
  dae_nu_perfil  nu_perfil;
  /* A tabela do perfil vive AQUI e nao no dae_nu_perfil, porque o spec e dono
     da memoria: dae_nu_perfil so guarda ponteiros const para estes arrays. */
  double        *perf_r, *perf_rho, *perf_ye;
  double        *ruido_amp, *ruido_fase;
  int32_t        perf_n;
  double         nu_tol;             /* tolerancia do guarda de adiabaticidade */

  /* VARREDURA EM ENERGIA. `scan_pontos == 0` significa ausente.
     Existe porque experimento mede ESPECTRO, nao ponto — e porque um ponto so
     nao separa fase CP de efeito de materia: com delta_cp = 0 e materia, a
     assimetria nu/antinu ja vale +0,42, e o que distingue as duas causas e a
     dependencia em E. `realizations` nao serve: ele varre SEMENTE, e a
     dinamica unitaria de sabor nao tem nada de estocastico. */
  int32_t        scan_pontos;
  double         scan_min, scan_max;   /* MeV */
  int            scan_log;             /* espacamento logaritmico            */
  int            scan_ambos;           /* roda nu E antinu, e emite A_CP     */
  int32_t        nu_sabor;           /* 0 = e, 1 = mu, 2 = tau                */

  /* SETOR MAGNÉTICO (format_version >= 4). `numag_ativo == 0` significa que os
     blocos "magnetic" e "trajectory" não existem. Com eles, quem roda é
     dae_numag_run, não dae_run: a física é outra (helicidade, não sabor) e o
     resultado tem o seu próprio CSV. Exclusivo com o bloco "neutrino".
       - `time.t1` é o comprimento a propagar, em km, medido a partir de
         s_km[0] no comprimento local do observador estático;
       - `time.nt` é o número de fatias, e de linhas do CSV.
     As tabelas vivem aqui (o spec é dono da memória); `numag_traj` só guarda
     ponteiros const para elas. */
  int                   numag_ativo;
  int                   numag_tem_traj;
  dae_numag_params      numag;
  double                numag_tol;
  int32_t               numag_helicidade;   /* 0 = left, 1 = right          */
  int                   numag_base_declarada;
  int                   numag_integrador;   /* 0 = brute_force, 1 = windows */
  int                   numag_coerencia;    /* -1 = não declarada, 0 = coherent,
                                               1 = incoherent               */
  dae_numag_janelas_cfg numag_jan;
  /* Colisões: taxa Gamma de interação do nu_L, em eV, constante ou por
     nó da trajetória; com ela, dae_numag_run roda o ENSEMBLE de trajetórias
     (dae_numag_ensemble) com `numag_n_traj` membros e semente `seed`. */
  double                numag_gamma_eV;
  double               *traj_gamma;          /* opcional, eV por nó         */
  int32_t               numag_n_traj;
  double               *traj_s, *traj_rho, *traj_ye, *traj_b1, *traj_b2, *traj_e;
  double               *traj_ynu;            /* opcional: Y_nue por nó       */
  int32_t               traj_n;
  dae_numag_trajetoria  numag_traj;
} dae_spec;

void       dae_spec_default(dae_spec *S);
dae_status dae_spec_parse(dae_spec *S, const char *json, dae_error *err);
void       dae_spec_free(dae_spec *S);

/* 1 se o core_hash do arquivo bate com o núcleo que está executando. Ausente
 * conta como "não confere" e o chamador decide o que fazer: um resultado sem
 * hash não é rastreável, mas também não está errado. */
int        dae_spec_hash_confere(const dae_spec *S);

/* JSON canônico: ordem de chaves fixa, números com 17 dígitos. É o que entra no
 * cabeçalho do CSV e no bloco de reprodutibilidade dos arquivos exportados.
 * Devolve quantos bytes escreveria (como snprintf); trunca se não couber. */
int32_t    dae_spec_canonical(const dae_spec *S, char *buf, int32_t cap);

#endif /* DAE_SPEC_H */
