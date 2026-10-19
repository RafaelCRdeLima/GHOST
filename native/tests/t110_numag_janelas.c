/* t110_numag_janelas.c — ACEITAÇÃO 110: o integrador por janelas contra a
 * FORÇA BRUTA nos mesmos extremos, nos casos que o perfil de um só cruzamento
 * não cobre.
 *
 * A comparação é sempre do OPERADOR DE TRANSFERÊNCIA COMPLEXO U (2x2, as duas
 * helicidades iniciais), não só de P: populações certas podem esconder fase
 * errada, e a fase relativa decide a interferência entre cruzamentos.
 *
 *  1. DOIS CRUZAMENTOS (rampa, platô, rampa). U contra a força bruta. E o modo
 *     INCOERENTE contra a MÉDIA do coerente sobre 4 comprimentos de platô que
 *     varrem um período da fase relativa chi. As populações ADIABÁTICAS no fim
 *     são A + B cos(chi + c) — só o 1o harmônico —, e a média de N >= 2 pontos
 *     igualmente espaçados é exata. A comparação é nas populações adiabáticas,
 *     não em P_R: P_R coerente tem um termo cruzado de borda, ~2 theta_fim,
 *     que não depende de chi e não some na média; o incoerente o descarta por
 *     definição. O mesmo vale para a ponta INICIAL: nu_L é superposição dos
 *     dois autoestados (mistura ~1/K), e o termo cruzado dela também não
 *     depende de chi. Por isso o caso parte do AUTOESTADO adiabático e mede a
 *     população adiabática no fim. Esses termos são erro de BORDA.
 *  2. CAMPO GIRANTE: Delta < 0 em todo o domínio, mas Delta_ef = Delta +
 *     phi'/T cruza zero. As janelas têm de achar o cruzamento ONDE Delta_ef se
 *     anula, e U tem de bater com a força bruta. A rotação tabelada por
 *     componentes dá phi' CONTÍNUO nos nós (cada intervalo gira o mesmo
 *     ângulo), então `saltos_tabela` tem de ser ~0: os saltos de theta~ só
 *     aparecem com rotação não uniforme.
 *  3. |Delta| GRANDE NÃO BASTA. Omega modulado em ressonância com o gap, com
 *     |Delta| = 300 x 2 Omega, fora da janela K = 100. Com o critério completo:
 *     U da força bruta. Com eps_max desligado (só a cota de variação total):
 *     ainda U da força bruta. Com os dois desligados: o integrador trata o
 *     trecho como adiabático e erra grosseiramente — é a companheira que prova
 *     que o critério está fazendo o trabalho.
 *  4. RESOLUÇÃO DA TABELA, separada do integrador: um perfil CURVO com N, 2N e
 *     4N nós, partindo do autoestado e medindo a população adiabática no fim.
 *     A diferença entre tabulações cai ~4x por duplicação (interpolação linear,
 *     O(h^2)). Não se usa P_R no fim: a fase dinâmica (~1e4 rad) muda O(1) rad
 *     com a tabulação e move o termo oscilante de borda, que domina P_R.
 */
#include "harness.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define TK DAE_NU_T_POR_KM
#define RHO 1.0e14
#define VR (DAE_NU_V_POR_RHO * RHO)

typedef struct { double *s, *rho, *ye, *b1, *b2; int32_t n; } tab;

static void tab_aloca(tab *t, int32_t n)
{
  t->n = n;
  t->s = (double *)calloc((size_t)n, sizeof(double));
  t->rho = (double *)calloc((size_t)n, sizeof(double));
  t->ye = (double *)calloc((size_t)n, sizeof(double));
  t->b1 = (double *)calloc((size_t)n, sizeof(double));
  t->b2 = (double *)calloc((size_t)n, sizeof(double));
}

static void tab_free(tab *t)
{
  free(t->s); free(t->rho); free(t->ye); free(t->b1); free(t->b2);
  memset(t, 0, sizeof(*t));
}

static dae_numag_trajetoria traj_de(const tab *t)
{
  dae_numag_trajetoria T;
  T.s_km = t->s; T.rho = t->rho; T.ye = t->ye; T.b1 = t->b1; T.b2 = t->b2;
  T.razao_energia = NULL; T.ynu = NULL; T.n = t->n;
  return T;
}

/* Y_e que dá Delta (eV) com rho = RHO */
static double ye_de(double delta) { return (2.0 * delta / VR + 1.0) / 3.0; }

/* U[linha*2+coluna] por força bruta, em `nf` fatias */
static dae_status U_bruta(const dae_numag_params *P, dae_numag_trajetoria *T, double s0,
                          double s1, int32_t nf, double ur[4], double ui[4])
{
  int col, k;
  for (col = 0; col < 2; ++col) {
    double re[4], im[4], pr[2], pi[2];
    dae_nu_inicial(col, 2, 1, re, im);
    for (k = 0; k < nf; ++k) {
      const dae_status st = dae_numag_propaga_tabela(P, T, s0 + (s1 - s0) * k / nf,
                                                     s0 + (s1 - s0) * (k + 1) / nf,
                                                     1e-11, re, im, NULL);
      if (st != DAE_OK) return st;
    }
    dae_nu_projeta(re, im, 2, 1, pr, pi);
    ur[col] = pr[0]; ui[col] = pi[0]; ur[2 + col] = pr[1]; ui[2 + col] = pi[1];
  }
  return DAE_OK;
}

static dae_status U_jan(const dae_numag_params *P, const dae_numag_trajetoria *T, double s0,
                        double s1, const dae_numag_janelas_cfg *c, double ur[4], double ui[4],
                        dae_numag_janelas_info *inf)
{
  int col;
  for (col = 0; col < 2; ++col) {
    double re[4], im[4], pr[2], pi[2];
    dae_status st;
    dae_nu_inicial(col, 2, 1, re, im);
    st = dae_numag_propaga_janelas(P, T, s0, s1, c, re, im, col == 0 ? inf : NULL);
    if (st != DAE_OK) return st;
    dae_nu_projeta(re, im, 2, 1, pr, pi);
    ur[col] = pr[0]; ui[col] = pi[0]; ur[2 + col] = pr[1]; ui[2 + col] = pi[1];
  }
  return DAE_OK;
}

static double dU(const double ar[4], const double ai[4], const double br[4], const double bi[4])
{
  double m = 0.0;
  int k;
  for (k = 0; k < 4; ++k) { const double d = hypot(ar[k] - br[k], ai[k] - bi[k]); if (d > m) m = d; }
  return m;
}

/* População adiabática (+) no fim, modo incoerente, partindo de nu_L */
static dae_status pop_inc(const dae_numag_params *P, const dae_numag_trajetoria *T, double s0,
                          double s1, const dae_numag_janelas_cfg *c, double *p_mais)
{
  double th, fi, pop[2];
  dae_status st;
  (void)th; (void)fi;
  pop[0] = 0.0; pop[1] = 1.0;           /* o autoestado (-), como no coerente */
  st = dae_numag_propaga_janelas_inc(P, T, s0, s1, c, pop, NULL, NULL);
  if (st != DAE_OK) return st;
  *p_mais = pop[0];
  return DAE_OK;
}

/* Estado físico -> mergulhado (psi, -i psi)/sqrt(2) */
static void para_merg(const double *pr, const double *pi, double *re, double *im)
{
  const double c = 1.0 / sqrt(2.0);
  int a;
  for (a = 0; a < 2; ++a) {
    re[a] = pr[a] * c; im[a] = pi[a] * c;
    re[2 + a] = pi[a] * c; im[2 + a] = -pr[a] * c;
  }
}

/* O autoestado (-) de H~ em s0, no laboratório: R(s0) (-sin th, cos th) */
static void autoestado_menos(const dae_numag_params *P, const dae_numag_trajetoria *T,
                             double s0, double *re, double *im)
{
  double th, fi, pr[2], pi[2];
  dae_numag_angulos(P, T, s0, 0, &th, &fi);
  pr[0] = -sin(th) * cos(0.5 * fi); pi[0] = -sin(th) * sin(0.5 * fi);
  pr[1] = cos(th) * cos(0.5 * fi);  pi[1] = -cos(th) * sin(0.5 * fi);
  para_merg(pr, pi, re, im);
}

/* população adiabática (+) no fim, coerente, partindo do autoestado (-) */
static double pop_mais_coer(const dae_numag_params *P, const dae_numag_trajetoria *T,
                            double s0, double s1, const dae_numag_janelas_cfg *c,
                            dae_numag_janelas_info *inf, dae_status *st);

/* |<+~(s1)| R(s1)^dag psi|^2 para o estado coerente psi = coluna 0 de U */
static double pop_mais(const dae_numag_params *P, const dae_numag_trajetoria *T, double s1,
                       const double ur[4], const double ui[4])
{
  double th, fi, c, s, t0r, t0i, t1r, t1i;
  dae_numag_angulos(P, T, s1, 1, &th, &fi);
  c = cos(0.5 * fi); s = sin(0.5 * fi);
  t0r = c * ur[0] + s * ui[0]; t0i = c * ui[0] - s * ur[0];   /* e^{-i fi/2} psi0 */
  t1r = c * ur[2] - s * ui[2]; t1i = c * ui[2] + s * ur[2];   /* e^{+i fi/2} psi1 */
  c = cos(th) * t0r + sin(th) * t1r;
  s = cos(th) * t0i + sin(th) * t1i;
  return c * c + s * s;
}

static double pop_mais_coer(const dae_numag_params *P, const dae_numag_trajetoria *T,
                            double s0, double s1, const dae_numag_janelas_cfg *c,
                            dae_numag_janelas_info *inf, dae_status *st)
{
  double re[4], im[4], ur[4], ui[4], pr[2], pi[2];
  autoestado_menos(P, T, s0, re, im);
  *st = dae_numag_propaga_janelas(P, T, s0, s1, c, re, im, inf);
  dae_nu_projeta(re, im, 2, 1, pr, pi);
  ur[0] = pr[0]; ui[0] = pi[0]; ur[2] = pr[1]; ui[2] = pi[1];
  ur[1] = ui[1] = ur[3] = ui[3] = 0.0;
  return pop_mais(P, T, s1, ur, ui);
}

/* rampa-platô-rampa: Y_e sobe com inclinação g de 1/3 - d a 1/3 + d, fica
   em 1/3 + d por `lp` km e desce. Nós: nr por rampa, 3 no platô. */
static void dois_cruzamentos(tab *t, double g, double d, double lp, int32_t nr)
{
  const double lr = 2.0 * d / g;
  int32_t i, k = 0;
  tab_aloca(t, 2 * nr + 2);
  for (i = 0; i < nr; ++i) { t->s[k] = lr * i / (nr - 1); t->ye[k] = 1.0 / 3.0 - d + g * t->s[k]; ++k; }
  t->s[k] = lr + 0.5 * lp; t->ye[k] = 1.0 / 3.0 + d; ++k;
  t->s[k] = lr + lp; t->ye[k] = 1.0 / 3.0 + d; ++k;
  for (i = 1; i < nr; ++i) {
    t->s[k] = lr + lp + lr * i / (nr - 1);
    t->ye[k] = 1.0 / 3.0 + d - g * (t->s[k] - lr - lp); ++k;
  }
  t->n = k;
  for (i = 0; i < k; ++i) { t->rho[i] = RHO; t->b1[i] = 1e15; t->b2[i] = 0.0; }
}

int main(void)
{
  dae_test T;
  dae_numag_params P;
  dae_numag_janelas_cfg cfg;
  const double om12 = 1e-12 * DAE_NUMAG_MU_B_EV_POR_T * 1e15 * DAE_NUMAG_T_POR_GAUSS;
  double br[4], bi[4], jr[4], ji[4];
  dae_status st;

  dae_test_begin(&T, "t110 integrador por janelas contra forca bruta nos mesmos extremos");
  dae_numag_params_default(&P);
  dae_numag_janelas_default(&cfg);

  /* ------------------------------------------------ 1. dois cruzamentos -- */
  {
    const double g = 0.1, d = 200.0 * 2.0 * om12 / (1.5 * VR);
    const double Ep = 1.5 * VR * d;                       /* Delta no platô, eV */
    const double per = 2.0 * PI / (Ep * TK);              /* km: um período de fase */
    const double gamma = 2.0 * PI * om12 * om12 / (1.5 * VR * g / TK);
    const double q = 1.0 - exp(-gamma);
    double media = 0.0, pinc = 0.0, pior = 0.0, pr0 = 0.0;
    tab t;
    dae_numag_trajetoria Tr;
    dae_numag_janelas_info inf;
    int j;
    dae_numag_janelas_cfg c1 = cfg;
    /* K = 50 < 200 do platô: os cruzamentos ficam em janelas SEPARADAS, com
       trecho adiabático entre eles — sem isso não há fase relativa a
       descartar no modo incoerente. */
    c1.K = 50.0;
    P.mu_re[0] = 1e-12;
    for (j = 0; j < 4; ++j) {
      dois_cruzamentos(&t, g, d, 1e-3 + per * j / 4.0, 41);
      Tr = traj_de(&t);
      memset(&inf, 0, sizeof inf);
      st = U_jan(&P, &Tr, t.s[0], t.s[t.n - 1], &c1, jr, ji, &inf);
      if (st != DAE_OK) { dae_test_ok(&T, 0, "1. janelas: %s", dae_strerror(st)); return dae_test_end(&T); }
      media += pop_mais_coer(&P, &Tr, t.s[0], t.s[t.n - 1], &c1, NULL, &st) / 4.0;
      if (j == 0) {
        pr0 = jr[2] * jr[2] + ji[2] * ji[2];
        st = U_bruta(&P, &Tr, t.s[0], t.s[t.n - 1], 40, br, bi);
        pior = dU(jr, ji, br, bi);
        dae_test_ok(&T, inf.n_cruzamentos == 2 && inf.n_janelas == 2,
                    "1. dois cruzamentos, duas janelas (%d, %d)", inf.n_cruzamentos, inf.n_janelas);
        dae_test_ok(&T, pior < 1e-5, "1. coerente: ||U_jan - U_bf|| = %.1e (cota %.1e)",
                    pior, inf.cota_adiabatica);
        st = pop_inc(&P, &Tr, t.s[0], t.s[t.n - 1], &c1, &pinc);
      }
      tab_free(&t);
    }
    /* comecando no autoestado (-), o (+) no fim e o caminho com EXATAMENTE um
       salto nao adiabatico: q(1-q) + (1-q)q */
    dae_test_ok(&T, fabs(pinc - media) < 1e-5,
                "1. incoerente = media do coerente sobre a fase relativa (populacao adiabatica): "
                "%.8f contra %.8f (dif %.1e)", pinc, media, fabs(pinc - media));
    dae_test_ok(&T, fabs(pinc - 2.0 * q * (1.0 - q)) < 2.0 / 200.0,
                "1. e a formula incoerente 2q(1-q) de LZ: %.5f contra %.5f (resolucao 2/K)",
                pinc, 2.0 * q * (1.0 - q));
    dae_test_note("   . P_R coerente partindo de nu_L, uma fase: %.5f (os termos de borda nao entram acima)", pr0);
  }

  /* ---------------------------------------------- 2. campo girante ---- */
  {
    /* Omega menor deixa a rotação tabelável: mu = 1e-13 */
    const double om = om12 * 0.1;
    const double gamma = 1.0;
    const double dD = 2.0 * PI * om * om * TK / gamma;   /* eV/km */
    const double D0 = 60.0 * 2.0 * om;                   /* Delta(0) = -D0 */
    const double L = 0.9 * D0 / dD;                      /* Delta < 0 em [0, L] */
    const double c = 30.0 * 2.0 * om * TK;               /* phi' (rad/km): Delta_ef = Delta + 30 x 2 Omega */
    const double s_ef = (D0 - c / TK) / dD;              /* zero de Delta_ef */
    const double rot = c * L;                            /* rad no domínio */
    double sal[2] = { 0.0, 0.0 };
    int rep;
    dae_numag_janelas_cfg c2 = cfg;
    c2.K = 10.0;
    P.mu_re[0] = 1e-13;
    for (rep = 0; rep < 2; ++rep) {
      const int32_t n = (int32_t)(rot / (rep ? 0.025 : 0.05)) + 2;
      tab t;
      dae_numag_trajetoria Tr;
      dae_numag_janelas_info inf;
      int32_t i, troca = 0;
      tab_aloca(&t, n);
      for (i = 0; i < n; ++i) {
        t.s[i] = L * i / (n - 1);
        t.rho[i] = RHO;
        t.ye[i] = ye_de(-D0 + dD * t.s[i]);
        t.b1[i] = 1e15 * cos(c * t.s[i]);
        t.b2[i] = 1e15 * sin(c * t.s[i]);
        if (i && ((-D0 + dD * t.s[i]) < 0.0) != ((-D0 + dD * t.s[i - 1]) < 0.0)) ++troca;
      }
      Tr = traj_de(&t);
      memset(&inf, 0, sizeof inf);
      st = U_jan(&P, &Tr, 0.0, L, &c2, jr, ji, &inf);
      sal[rep] = inf.saltos_tabela;
      if (rep == 0) {
        dae_test_ok(&T, st == DAE_OK, "2. janelas no campo girante (%s)", dae_strerror(st));
        dae_test_ok(&T, troca == 0, "2. Delta < 0 em todo o dominio (%d trocas de sinal)", (int)troca);
        /* Com a rotação tabelada por componentes, |b| encolhe no meio de cada
           intervalo e phi' oscila de leve: perto do zero, Delta_ef pode trocar
           de sinal mais de uma vez num aglomerado estreito. É o meio
           interpolado, e o aglomerado é UMA ressonância. */
        dae_test_ok(&T, inf.n_cruzamentos >= 1 && fabs(inf.cruzamento_s_km - s_ef) < 1e-3 * L,
                    "2. o cruzamento achado e o zero de Delta_ef: s = %.6e km (esperado %.6e; "
                    "%d troca(s) de sinal no aglomerado)", inf.cruzamento_s_km, s_ef, inf.n_cruzamentos);
        st = U_bruta(&P, &Tr, 0.0, L, 200, br, bi);
        dae_test_ok(&T, st == DAE_OK && dU(jr, ji, br, bi) < 1e-5,
                    "2. ||U_jan - U_bf|| = %.1e (%d nos, %.0f rad de rotacao)",
                    dU(jr, ji, br, bi), (int)n, rot);
      }
      tab_free(&t);
    }
    dae_test_ok(&T, sal[0] < 1e-10 && sal[1] < 1e-10,
                "2. rotacao uniforme por cordas: phi' continuo nos nos, saltos de theta~ ~0 (%.1e, %.1e)",
                sal[0], sal[1]);
  }

  /* --------------------------------- 3. |Delta| grande nao basta ---- */
  {
    const double a = 0.5;
    const double Dc = 300.0 * 2.0 * om12;                /* fora de K = 100 */
    const double Ec = sqrt(Dc * Dc + 4.0 * om12 * om12);
    const double k = Ec * TK;                            /* modulação no gap */
    const double L = PI / (a * om12) / TK;               /* ~ uma transferência */
    const int32_t n = (int32_t)(L * k / (2.0 * PI) * 24.0) + 2;
    dae_numag_janelas_cfg c1 = cfg, c2 = cfg, c3 = cfg;
    double pbf;
    tab t;
    dae_numag_trajetoria Tr;
    int32_t i;
    P.mu_re[0] = 1e-12;
    tab_aloca(&t, n);
    for (i = 0; i < n; ++i) {
      t.s[i] = L * i / (n - 1);
      t.rho[i] = RHO;
      t.ye[i] = ye_de(Dc);
      t.b1[i] = 1e15 * (1.0 + a * cos(k * t.s[i]));
      t.b2[i] = 0.0;
    }
    Tr = traj_de(&t);
    c1.K = c2.K = c3.K = 100.0;
    c2.eps_max = 1e9;                       /* só a cota de variação total */
    c3.eps_max = 1e9; c3.eps_tol = 1e9;     /* tudo desligado */
    st = U_bruta(&P, &Tr, 0.0, L, 100, br, bi);
    pbf = br[2] * br[2] + bi[2] * bi[2];
    dae_test_ok(&T, st == DAE_OK && pbf > 0.3,
                "3. a forca bruta mostra a conversao ressonante: P_R = %.3f com |Delta| = 300 x 2 Omega", pbf);
    U_jan(&P, &Tr, 0.0, L, &c1, jr, ji, NULL);
    dae_test_ok(&T, dU(jr, ji, br, bi) < 1e-5, "3. criterio completo: ||U_jan - U_bf|| = %.1e", dU(jr, ji, br, bi));
    U_jan(&P, &Tr, 0.0, L, &c2, jr, ji, NULL);
    dae_test_ok(&T, dU(jr, ji, br, bi) < 1e-5,
                "3. so a cota de variacao total (eps_max desligado): ||U_jan - U_bf|| = %.1e", dU(jr, ji, br, bi));
    U_jan(&P, &Tr, 0.0, L, &c3, jr, ji, NULL);
    dae_test_ok(&T, fabs(jr[2] * jr[2] + ji[2] * ji[2] - pbf) > 0.2,
                "3. companheira: sem criterio, o trecho 'adiabatico' erra P_R (%.3f contra %.3f)",
                jr[2] * jr[2] + ji[2] * ji[2], pbf);
    tab_free(&t);
  }

  /* ------------------------------------ 4. resolução da tabela ---- */
  {
    /* Y_e = 1/3 + d tanh((s - s*)/w): um cruzamento, perfil curvo */
    const double om = om12, d = 300.0 * 2.0 * om / (1.5 * VR), w = 2.0e-3;
    const double L = 1.2e-2, s0c = 0.5 * L;
    static const int32_t NS[3] = { 81, 161, 321 };
    double p[3];
    int r;
    dae_numag_janelas_cfg c4 = cfg;
    c4.K = 50.0;
    P.mu_re[0] = 1e-12;
    for (r = 0; r < 3; ++r) {
      tab t;
      dae_numag_trajetoria Tr;
      int32_t i;
      tab_aloca(&t, NS[r]);
      for (i = 0; i < NS[r]; ++i) {
        t.s[i] = L * i / (NS[r] - 1);
        t.rho[i] = RHO;
        t.ye[i] = 1.0 / 3.0 + d * tanh((t.s[i] - s0c) / w);
        t.b1[i] = 1e15; t.b2[i] = 0.0;
      }
      Tr = traj_de(&t);
      p[r] = pop_mais_coer(&P, &Tr, 0.0, L, &c4, NULL, &st);
      tab_free(&t);
    }
    dae_test_ok(&T, st == DAE_OK && fabs(p[1] - p[0]) > 1e-6 && fabs(p[2] - p[1]) < 0.35 * fabs(p[1] - p[0]),
                "4. resolucao da tabela (populacao adiabatica): dif N->2N %.1e -> %.1e (cai ~4x, O(h^2))",
                fabs(p[1] - p[0]), fabs(p[2] - p[1]));
  }

  return dae_test_end(&T);
}
