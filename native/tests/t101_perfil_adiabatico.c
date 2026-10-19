/* t101_perfil_adiabatico.c — ACEITAÇÃO 101: conversão MSW num perfil de
 * supernova, e o guarda que impede o resultado errado em silêncio.
 *
 * O CASO. Um envelope rho ~ r^-3, nu_e saindo do centro, atravessando as
 * ressonâncias H (Dm31^2) e L (Dm21^2). Com a escada de densidade grossa
 * demais na ressonância, P(nu_e->nu_e) sai 0,0027 contra os 0,0224 corretos —
 * errado por um fator 8, com a curva na FORMA CERTA. Um expm de referência
 * erra igual: não é bug de propagador, é erro de discretização, e por isso
 * comparar com outro código não pega.
 *
 * A MÉTRICA NÃO É P(nu_e). Essa oscila e engana — ela vale 0,0094, 0,0027,
 * 0,0225 e 0,0225 para N crescente, sem monotonia. A métrica é
 * |<nu_3|psi>|^2 -> 1: se a conversão foi adiabática, o estado sai como nu_3
 * puro, e aí a probabilidade INCOERENTE que chega à Terra é |U_e3|^2 = sin^2
 * theta13, sem depender de fase nenhuma.
 *
 * ANTI-VACUIDADE (t101b): com tolerância frouxa o teste tem de FALHAR. Sem
 * isso, um guarda que aceitasse qualquer coisa passaria no teste principal, e
 * a única coisa provada seria que 6000 fatias bastam.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

/* Roda o envelope inteiro em `nfatias` fatias logarítmicas, cada uma resolvida
 * com controle de passo adaptativo até `tol`. Devolve |<nu_i|psi>|^2 e o total
 * de sub-passos que o guarda julgou necessário. */
static dae_status envelope(dae_test *T, double tol, int32_t nfatias,
                           double *p_massa, int32_t *subpassos, double *p_inc)
{
  dae_nu_params P;
  dae_nu_perfil M;
  double hre[9], him[9], ure[9], uim[9];
  double re[6], im[6], psi_re[3], psi_im[3];
  /* Comeca em 1000 km (rho ~ 1e7 g/cm3), nao no centro. Nao e conveniencia:
     mais para dentro V domina Dm2/2E por ~1e7 e a fase acumulada por km fica
     grande demais para QUALQUER discretizacao pratica — e fisicamente aquele
     regime nao faz nada, porque com V dominante o nu_e ja E o autoestado de
     materia. As duas ressonancias (H em ~18 500 km, L em ~80 000 km) estao
     inteiras dentro da janela. */
  const double r0 = 1000.0, r1 = 1.0e6;         /* km */
  int merg = 0;
  int32_t i, a, m;
  dae_status st;

  dae_nu_params_default(&P);
  P.delta_cp = -PI / 2.0;                        /* exige o mergulho */
  P.energy_eV = 20.0e6;

  dae_nu_perfil_default(&M);
  M.kind = DAE_NU_PERFIL_POTENCIA;
  M.ye = 0.5;
  M.rho_ref = 1.0e10;                            /* g/cm3 em r_ref */
  M.r_ref = 100.0;
  M.expoente = 3.0;
  M.r_start = r0;

  /* decide o mergulho UMA vez: Im(H) não depende de rho */
  st = dae_nu_hamiltoniano(&P, 1.0, M.ye, hre, him);
  if (st != DAE_OK) return st;
  merg = dae_nu_precisa_mergulho(him, P.n_flavors);
  dae_test_ok(T, merg == 1, "delta_cp=-90 exige mergulho no perfil");

  dae_nu_inicial(0, P.n_flavors, merg, re, im);   /* nu_e */
  *subpassos = 0;
  for (i = 0; i < nfatias; ++i) {
    const double f = (double)i / (double)nfatias;
    const double g = (double)(i + 1) / (double)nfatias;
    const double ra = r0 * pow(r1 / r0, f);
    const double rb = r0 * pow(r1 / r0, g);
    int32_t sub = 0;
    st = dae_nu_propaga(&P, &M, ra, rb, tol, merg, re, im, &sub);
    if (st != DAE_OK) return st;
    *subpassos += sub;
  }

  dae_nu_projeta(re, im, P.n_flavors, merg, psi_re, psi_im);
  st = dae_nu_pmns(&P, ure, uim);
  if (st != DAE_OK) return st;

  /* |<nu_m|psi>|^2 = |sum_a conj(U_am) psi_a|^2 */
  for (m = 0; m < 3; ++m) {
    double sr = 0.0, si = 0.0;
    for (a = 0; a < 3; ++a) {
      const double cr = ure[a * 3 + m], ci = -uim[a * 3 + m];
      sr += cr * psi_re[a] - ci * psi_im[a];
      si += cr * psi_im[a] + ci * psi_re[a];
    }
    p_massa[m] = sr * sr + si * si;
  }
  /* O que a Terra recebe: a coerência morre em ~32 km, muito antes de sair da
     galáxia, então o observável é a soma INCOERENTE sobre autoestados de massa. */
  *p_inc = 0.0;
  for (m = 0; m < 3; ++m) {
    const double ue = ure[0 * 3 + m], ui = uim[0 * 3 + m];
    *p_inc += (ue * ue + ui * ui) * p_massa[m];
  }
  return DAE_OK;
}

int main(void)
{
  dae_test T;
  dae_nu_params P;
  double p_massa[3], p_inc = 0.0;
  int32_t sub = 0;
  double sin2t13;
  dae_status st;

  dae_test_begin(&T, "t101 perfil de supernova: conversao adiabatica e o guarda");

  dae_nu_params_default(&P);
  sin2t13 = sin(P.theta13) * sin(P.theta13);

  /* ------------------------------------------- o teste principal, tol apertada */
  st = envelope(&T, 1e-6, 120, p_massa, &sub, &p_inc);
  dae_test_ok(&T, st == DAE_OK, "o envelope rodou: %s", dae_strerror(st));
  if (st == DAE_OK) {
    dae_test_ok(&T, p_massa[2] > 1.0 - 1e-4,
                "a conversao e adiabatica: |<nu3|psi>|^2 = %.9f", p_massa[2]);
    dae_test_near(&T, p_inc, sin2t13, 1e-3,
                  "P incoerente contra sin^2(theta13)");
    dae_test_note("   . |<nu1|psi>|^2=%.3e  |<nu2|psi>|^2=%.3e  |<nu3|psi>|^2=%.9f",
                  p_massa[0], p_massa[1], p_massa[2]);
    dae_test_note("   . P(nu_e->nu_e) incoerente = %.9f   sin^2(theta13) = %.9f",
                  p_inc, sin2t13);
    dae_test_note("   . o guarda pediu %d sub-passos ao todo", (int)sub);
  }

  /* ---------------- ANTI-VACUIDADE: com tolerancia frouxa o teste tem de FALHAR */
  {
    double pm[3], pi_frouxo = 0.0;
    int32_t sub2 = 0;
    /* tol = 0.5 aceita a primeira tentativa em qualquer fatia: e a escada
       grossa que produz o resultado errado por um fator 8. */
    st = envelope(&T, 0.5, 8, pm, &sub2, &pi_frouxo);
    dae_test_ok(&T, st == DAE_OK, "a corrida frouxa tambem roda (o teste e sobre o VALOR)");
    if (st == DAE_OK) {
      dae_test_ok(&T, pm[2] < 1.0 - 1e-4,
                  "com tol frouxa a adiabaticidade SE PERDE: |<nu3|psi>|^2 = %.6f", pm[2]);
      /* NAO se cobra P incoerente aqui, e a razao e a licao deste arquivo:
         P e a metrica que engana. Na corrida frouxa ela errou por 4,7e-4 —
         pequeno o bastante para passar num limiar de 1e-3, enquanto a
         adiabaticidade ja tinha se perdido em 1,7e-3. Cobrar P faria a
         anti-vacuidade falhar em detectar exatamente o defeito que ela existe
         para detectar. */
      dae_test_note("   . anti-vacuidade: tol=0.5, 8 fatias -> |<nu3|psi>|^2=%.6f, "
                    "P_inc=%.6f (%d sub-passos)", pm[2], pi_frouxo, (int)sub2);
      dae_test_note("   . e o motivo do guarda: a curva teria a forma certa e o valor errado");
    }
  }

  /* ------------------------------- o guarda RECUSA quando nao consegue convergir */
  {
    dae_nu_perfil M;
    double re[6], im[6];
    dae_nu_params Q;
    int32_t s2 = 0;
    dae_nu_params_default(&Q);
    Q.delta_cp = -PI / 2.0;
    Q.energy_eV = 20.0e6;
    dae_nu_perfil_default(&M);
    M.kind = DAE_NU_PERFIL_POTENCIA;
    M.ye = 0.5; M.rho_ref = 1.0e10; M.r_ref = 100.0; M.expoente = 3.0;
    dae_nu_inicial(0, 3, 1, re, im);
    /* tolerancia impossivel numa fatia enorme: tem de RECUSAR, nao devolver */
    st = dae_nu_propaga(&Q, &M, 100.0, 1.0e6, 1e-15, 1, re, im, &s2);
    dae_test_ok(&T, st == DAE_ERR_TOOBIG,
                "tolerancia inatingivel e RECUSADA (%s), nao aproximada",
                dae_strerror(st));
    dae_test_note("   . o guarda recusa em vez de devolver numero: %s", dae_strerror(st));
  }

  /* ------------------------------------------ o perfil devolve o que promete */
  {
    dae_nu_perfil M;
    dae_nu_perfil_default(&M);
    M.kind = DAE_NU_PERFIL_POTENCIA;
    M.rho_ref = 1.0e10; M.r_ref = 100.0; M.expoente = 3.0;
    dae_test_near(&T, dae_nu_densidade(&M, 100.0), 1.0e10, 1.0, "rho(r_ref) = rho_ref");
    dae_test_near(&T, dae_nu_densidade(&M, 200.0), 1.25e9, 1.0, "rho cai com r^-3");
    M.kind = DAE_NU_PERFIL_CONSTANTE; M.rho = 42.0;
    dae_test_near(&T, dae_nu_densidade(&M, 1.0e9), 42.0, 1e-12,
                  "perfil constante ignora r");
  }

  /* ------------------------------ o perfil de TABELA interpola e segura ---- */
  {
    /* Uma tabela que E uma lei de potencia exata: a interpolacao log-log tem
       de reproduzi-la em qualquer ponto, nao so nos nos. E o teste certo —
       comparar tabela com tabela nao pegaria um erro de interpolacao. */
    static double rt[5], pt[5];
    dae_nu_perfil M, Q;
    int32_t i;
    double pior = 0.0;
    for (i = 0; i < 5; ++i) {
      rt[i] = 100.0 * pow(10.0, (double)i);
      pt[i] = 1.0e10 * pow(100.0 / rt[i], 3.0);
    }
    dae_nu_perfil_default(&M);
    M.kind = DAE_NU_PERFIL_TABELA;
    M.r_tab = rt; M.rho_tab = pt; M.n_tab = 5; M.ye = 0.5;
    dae_nu_perfil_default(&Q);
    Q.kind = DAE_NU_PERFIL_POTENCIA;
    Q.rho_ref = 1.0e10; Q.r_ref = 100.0; Q.expoente = 3.0;
    for (i = 0; i < 60; ++i) {
      const double r = 100.0 * pow(1.0e4, (double)i / 59.0);
      const double a = dae_nu_densidade(&M, r), b = dae_nu_densidade(&Q, r);
      const double d = fabs(a - b) / b;
      if (d > pior) pior = d;
    }
    dae_test_ok(&T, pior < 1e-12,
                "tabela log-log reproduz a lei de potencia entre os nos (%.3e)", pior);

    /* Fora da tabela SEGURA na borda, nao extrapola: extrapolar lei de
       potencia num perfil estelar erra por ordens de grandeza, em silencio. */
    dae_test_near(&T, dae_nu_densidade(&M, 1.0), pt[0], 1e-12, "abaixo da tabela: segura");
    dae_test_near(&T, dae_nu_densidade(&M, 1.0e9), pt[4], 1e-12, "acima da tabela: segura");
    dae_test_ok(&T, dae_nu_densidade(&M, 1.0e9) > 0.0, "e o valor segurado e positivo");
    dae_test_note("   . perfil de tabela: interpolacao log-log a %.1e, e segura nas bordas", pior);
  }

  return dae_test_end(&T);
}
