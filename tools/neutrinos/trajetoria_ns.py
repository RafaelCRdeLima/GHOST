"""trajetoria_ns.py — o raio nulo numa estrela estática esférica, e o que o
setor magnético do GHOST precisa saber ao longo dele.

Este script NÃO propaga neutrino nenhum: integra a geometria e escreve as
tabelas do bloco "trajectory" do spec. Quem propaga é o núcleo
(dae_numag_run). A divisão é a de sempre no GHOST: a física de propagação
vive no C; a preparação das entradas vive aqui, e vai inteira para o spec.

MÉTRICA. ds^2 = -e^nu dt^2 + e^lambda dr^2 + r^2 dOmega^2, com G = c = 1 e
comprimentos em km (M em km: M_sol = 1,4766 km).

O COMPRIMENTO s. É o comprimento espacial LOCAL medido pelo observador
estático ao longo do raio, dl^2 = e^lambda dr^2 + r^2 dOmega^2. O intervalo
próprio da geodésica nula é zero e não serve de parâmetro; o parâmetro afim
não tem unidade física fixa. É em l que H (em eV, no referencial do
observador estático) gera a evolução: i dpsi/dl = H psi.

O RAIO. Com psi o ângulo local entre o momento e a direção radial, medido
pelo observador estático, a constante de movimento b = L/E dá

    sin psi = b e^{nu/2} / r,

e o raio obedece, em l,

    dr/dl   = e^{-lambda/2} cos psi
    dphi/dl = sin psi / r
    dpsi/dl = sin psi (nu'/2 - 1/r) e^{-lambda/2}

A terceira equação é a derivada da constante de movimento; com ela o ponto de
retorno (psi = pi/2) não é singular, ao contrário de escrever dr/dl com a raiz.

A BASE TRANSVERSAL TRANSPORTADA. Numa métrica estática esférica o raio fica
num plano, e por simetria de reflexão nesse plano a NORMAL à órbita L é
transportada paralelamente ao longo dele; o outro vetor transversal é então
fixado, e2 = k x e1. É a ausência de rotação de Faraday gravitacional em
Schwarzschild (ela aparece com rotação, em Kerr). Então

    e1 = L (constante),   e2 = k x e1,   e1 x e2 = k,

a mesma orientação de dae_numag_base_transversa. Uma base reconstruída por
regra local (z x k em cada ponto, por exemplo) gira em torno de k quando o
raio se curva, e essa rotação entraria em Delta_ef como um dphi_B/ds que não
existe — ver `rotacao_base_ingenua`.

O CAMPO. `campo(x, n)` devolve B (gauss) em componentes do referencial
ortonormal do observador estático, escritas nos eixos do "embedding" plano
(r-hat = n, e as direções tangenciais). O dipolo daqui é NEWTONIANO, com
interior uniforme: serve de campo de teste e de primeira aproximação, não de
modelo de protoestrela.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

M_SOL_KM = 1.4766250614046494       # G M_sol / c^2
T_POR_KM = 5067730727.952153        # DAE_NU_T_POR_KM


# ---------------------------------------------------------------- métricas --

class Plana:
    nome = "plana"

    def nu(self, r):
        return np.zeros_like(np.asarray(r, dtype=float))

    def dnu(self, r):
        return np.zeros_like(np.asarray(r, dtype=float))

    def lam(self, r):
        return np.zeros_like(np.asarray(r, dtype=float))

    def dlam(self, r):
        return np.zeros_like(np.asarray(r, dtype=float))


class Schwarzschild:
    """Exterior de massa M (km). Só vale para r > 2M."""
    nome = "schwarzschild"

    def __init__(self, M):
        self.M = float(M)

    def nu(self, r):
        return np.log(1.0 - 2.0 * self.M / np.asarray(r, dtype=float))

    def dnu(self, r):
        r = np.asarray(r, dtype=float)
        return 2.0 * self.M / (r * r * (1.0 - 2.0 * self.M / r))

    def lam(self, r):
        return -self.nu(r)

    def dlam(self, r):
        return -self.dnu(r)


class InteriorConstante:
    """Estrela de densidade constante (solução interior de Schwarzschild),
    casada com o exterior em R. Analítica, e por isso o teste do interior."""
    nome = "interior_constante"

    def __init__(self, M, R):
        self.M, self.R = float(M), float(R)
        self.k = 2.0 * self.M / self.R ** 3
        self.ext = Schwarzschild(M)
        if not 2.0 * M / R < 8.0 / 9.0:
            raise ValueError("Buchdahl: 2M/R precisa ser < 8/9")

    def _a(self, r):
        return (1.5 * np.sqrt(1.0 - 2.0 * self.M / self.R)
                - 0.5 * np.sqrt(1.0 - self.k * r * r))

    def nu(self, r):
        r = np.asarray(r, dtype=float)
        return np.where(r < self.R, 2.0 * np.log(self._a(np.minimum(r, self.R))),
                        self.ext.nu(np.maximum(r, self.R)))

    def dnu(self, r):
        r = np.asarray(r, dtype=float)
        ri = np.minimum(r, self.R)
        dentro = 2.0 * (0.5 * self.k * ri / np.sqrt(1.0 - self.k * ri * ri)) / self._a(ri)
        return np.where(r < self.R, dentro, self.ext.dnu(np.maximum(r, self.R)))

    def lam(self, r):
        r = np.asarray(r, dtype=float)
        return np.where(r < self.R, -np.log(1.0 - self.k * np.minimum(r, self.R) ** 2),
                        self.ext.lam(np.maximum(r, self.R)))

    def dlam(self, r):
        """lambda' SALTA em R (a densidade salta); nu' não."""
        r = np.asarray(r, dtype=float)
        ri = np.minimum(r, self.R)
        return np.where(r < self.R, 2.0 * self.k * ri / (1.0 - self.k * ri * ri),
                        self.ext.dlam(np.maximum(r, self.R)))


class Tabelada:
    """Interior tabelado (r, nu, lambda, por exemplo de uma integração TOV)
    casado com Schwarzschild de massa M a partir do último ponto R. nu e lambda
    por spline cúbica; as derivadas são as da spline. Toda métrica nova deve
    passar pelo teste de transporte covariante (`transporte_covariante`) antes de
    produzir trajetória."""
    nome = "tabelada"

    def __init__(self, r, nu, lam, M):
        from scipy.interpolate import CubicSpline
        r, nu, lam = map(np.asarray, (r, nu, lam))
        self.R, self.M = float(r[-1]), float(M)
        self.ext = Schwarzschild(M)
        self._nu = CubicSpline(r, nu)
        self._lam = CubicSpline(r, lam)

    def _sel(self, r, dentro, fora):
        r = np.asarray(r, dtype=float)
        return np.where(r < self.R, dentro(np.minimum(r, self.R)), fora(np.maximum(r, self.R)))

    def nu(self, r):
        return self._sel(r, self._nu, self.ext.nu)

    def dnu(self, r):
        return self._sel(r, self._nu.derivative(), self.ext.dnu)

    def lam(self, r):
        return self._sel(r, self._lam, self.ext.lam)

    def dlam(self, r):
        return self._sel(r, self._lam.derivative(), self.ext.dlam)


# ------------------------------------------------------------------ campos --

def campo_uniforme(B0, eixo=(0.0, 0.0, 1.0)):
    e = np.asarray(eixo, dtype=float)
    e = e / np.linalg.norm(e)
    return lambda x, n: B0 * e


def campo_dipolo_newtoniano(B0, R, eixo=(0.0, 0.0, 1.0)):
    """Interior uniforme B0 e, exterior dipolar com B_r contínuo em R
    (m = B0 R^3 / 2). Newtoniano: sem a correção relativística do dipolo."""
    e = np.asarray(eixo, dtype=float)
    e = e / np.linalg.norm(e)
    m = 0.5 * B0 * R ** 3

    def B(x, n):
        r = np.linalg.norm(x)
        if r < R:
            return B0 * e
        return (m / r ** 3) * (3.0 * np.dot(e, n) * n - e)
    return B


# -------------------------------------------------------------------- raio --

def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def raio(metrica, x0, k0, l_max, n_amostras=2001, campo=None, rtol=1e-12,
         max_passo=np.inf, r_saida=None):
    """Integra o raio nulo de x0 (km) na direção local k0 por l_max (km).

    Devolve um dict de arrays ao longo de l: r, phi, psi, x, k, e1, e2, n,
    razao_energia = e^{-nu/2}, b_const = r sin psi e^{-nu/2} (tem de ficar
    constante) e, com `campo`, b1, b2 e b_par.

    Com `r_saida`, para quando o raio cruza r = r_saida SAINDO, e devolve o
    ponto exato do cruzamento em `evento` = (l, r, phi, psi).
    """
    x0 = np.asarray(x0, dtype=float)
    r0 = np.linalg.norm(x0)
    n0 = x0 / r0
    k0 = _unit(k0)
    cpsi0 = np.clip(np.dot(k0, n0), -1.0, 1.0)
    psi0 = np.arccos(cpsi0)
    L = np.cross(n0, k0)
    if np.linalg.norm(L) < 1e-14:
        # raio radial: o plano é indefinido; qualquer normal fixa serve
        aux = np.array([0.0, 0.0, 1.0]) if abs(n0[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        L = np.cross(aux, n0)
    L = _unit(L)
    t0 = np.cross(L, n0)                       # k0 = cos psi0 n0 + sin psi0 t0
    b = r0 * np.sin(psi0) * np.exp(-0.5 * float(metrica.nu(r0)))

    def f(l, y):
        r, phi, psi = y
        el = np.exp(-0.5 * float(metrica.lam(r)))
        return [el * np.cos(psi),
                np.sin(psi) / r,
                np.sin(psi) * (0.5 * float(metrica.dnu(r)) - 1.0 / r) * el]

    def centro(l, y):                          # r -> 0 num raio radial para dentro
        return y[0] - 1e-9 * r0
    centro.terminal = True

    eventos = [centro]
    if r_saida is not None:
        def saida(l, y):
            return y[0] - r_saida
        saida.terminal = True
        saida.direction = 1
        eventos.append(saida)
    ls = np.linspace(0.0, l_max, n_amostras)
    sol = solve_ivp(f, (0.0, l_max), [r0, 0.0, psi0], method="DOP853",
                    t_eval=ls, rtol=rtol, atol=rtol * r0, events=eventos,
                    max_step=max_passo)
    if sol.status == 1 and len(sol.t_events[0]):
        raise RuntimeError("raio chegou ao centro")
    if sol.status < 0:
        raise RuntimeError(f"raio interrompido: {sol.message}")
    r, phi, psi = sol.y
    n = np.outer(np.cos(phi), n0) + np.outer(np.sin(phi), t0)
    tt = -np.outer(np.sin(phi), n0) + np.outer(np.cos(phi), t0)
    k = np.cos(psi)[:, None] * n + np.sin(psi)[:, None] * tt
    e1 = np.tile(L, (len(sol.t), 1))      # sol.t, não ls: o evento r_saida encurta o raio
    e2 = np.cross(k, e1)
    out = dict(l=sol.t, r=r, phi=phi, psi=psi, x=r[:, None] * n, n=n, k=k,
               e1=e1, e2=e2, b=b,
               razao_energia=np.exp(-0.5 * metrica.nu(r)),
               b_const=r * np.sin(psi) * np.exp(-0.5 * metrica.nu(r)))
    if r_saida is not None and len(sol.t_events[1]):
        out["evento"] = (sol.t_events[1][0],) + tuple(sol.y_events[1][0])
    if campo is not None:
        Bv = np.array([campo(out["x"][i], n[i]) for i in range(len(sol.t))])
        out["B"] = Bv
        out["b1"] = np.einsum("ij,ij->i", Bv, e1)
        out["b2"] = np.einsum("ij,ij->i", Bv, e2)
        out["b_par"] = np.einsum("ij,ij->i", Bv, k)
    return out


def base_ingenua(k):
    """A regra LOCAL de dae_numag_base_transversa aplicada ponto a ponto:
    e1 = z x k (ou x x k perto de z). Existe só para medir o erro de usá-la
    num raio curvo."""
    z, xh = np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])
    e1 = np.cross(z, k)
    perto = np.linalg.norm(e1, axis=1) < 1e-3
    e1[perto] = np.cross(xh, k[perto])
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    return e1, np.cross(k, e1)


def rotacao_base_ingenua(tr):
    """dphi/dl (rad/km) da base ingênua em torno de k, relativo à transportada.
    É o termo espúrio que ela poria em Delta_ef: dividido por T_POR_KM, em eV."""
    e1n, _ = base_ingenua(tr["k"])
    ang = np.unwrap(np.arctan2(np.einsum("ij,ij->i", e1n, tr["e2"]),
                               np.einsum("ij,ij->i", e1n, tr["e1"])))
    return np.gradient(ang, tr["l"])


# ---------------------------------------------------------------- o spec --

def spec_magnetico(tr, rho, ye, mu_bohr, energia_MeV, nt, tol=1e-10,
                   helicidade="left", t1=None, ynu=None):
    """Monta o spec (format_version 4) a partir de um raio com campo.
    `rho`, `ye` e (opcional) `ynu` = Y_nue são arrays ao longo de tr["l"].
    Sem `ynu`, o spec não leva ynu_table e Delta é o de matéria neutra."""
    l = tr["l"]
    sp = {
        "format_version": 4, "seed": 1,
        "magnetic": {"nature": "dirac", "mu_nu_bohr": float(mu_bohr),
                     "energy_MeV": float(energia_MeV),
                     "initial_helicity": helicidade, "tolerance": float(tol)},
        "trajectory": {"transverse_basis": "parallel_transported",
                       "s_km": [float(v) for v in l],
                       "rho_table": [float(v) for v in rho],
                       "ye_table": [float(v) for v in ye],
                       "b1_gauss": [float(v) for v in tr["b1"]],
                       "b2_gauss": [float(v) for v in tr["b2"]],
                       "energy_ratio": [float(v) for v in tr["razao_energia"]]},
        "time": {"t1": float(t1 if t1 is not None else l[-1] - l[0]), "nt": int(nt)},
        "realizations": 1,
    }
    if ynu is not None:
        sp["trajectory"]["ynu_table"] = [float(v) for v in ynu]
    return sp


def escreve_spec(caminho, spec):
    with open(caminho, "w") as fh:
        json.dump(spec, fh)


# ------------------------------------------------- transporte covariante --

def _coord(x, v, metrica):
    """Vetor espacial v (componentes ortonormais do observador estático, nos
    eixos do embedding) -> componentes coordenadas (t, r, theta, phi)."""
    r = np.linalg.norm(x, axis=1)
    th = np.arccos(np.clip(x[:, 2] / r, -1.0, 1.0))
    ph = np.arctan2(x[:, 1], x[:, 0])
    rh = x / r[:, None]
    thh = np.stack([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)], axis=1)
    phh = np.stack([-np.sin(ph), np.cos(ph), np.zeros_like(ph)], axis=1)
    el = np.exp(-0.5 * metrica.lam(r))
    out = np.zeros((len(r), 4))
    out[:, 1] = el * np.einsum("ij,ij->i", v, rh)
    out[:, 2] = np.einsum("ij,ij->i", v, thh) / r
    out[:, 3] = np.einsum("ij,ij->i", v, phh) / (r * np.sin(th))
    return out, r, th


def _christoffel(metrica, r, th, a, b):
    """c^mu = Gamma^mu_{alpha beta} a^alpha b^beta, métrica estática esférica
    em (t, r, theta, phi)."""
    nu, dnu = metrica.nu(r), metrica.dnu(r)
    lam, dlam = metrica.lam(r), metrica.dlam(r)
    s, c = np.sin(th), np.cos(th)
    out = np.zeros_like(a)
    out[:, 0] = 0.5 * dnu * (a[:, 0] * b[:, 1] + a[:, 1] * b[:, 0])
    out[:, 1] = (0.5 * dnu * np.exp(nu - lam) * a[:, 0] * b[:, 0]
                 + 0.5 * dlam * a[:, 1] * b[:, 1]
                 - r * np.exp(-lam) * a[:, 2] * b[:, 2]
                 - r * s * s * np.exp(-lam) * a[:, 3] * b[:, 3])
    out[:, 2] = (a[:, 1] * b[:, 2] + a[:, 2] * b[:, 1]) / r - s * c * a[:, 3] * b[:, 3]
    out[:, 3] = ((a[:, 1] * b[:, 3] + a[:, 3] * b[:, 1]) / r
                 + (c / s) * (a[:, 2] * b[:, 3] + a[:, 3] * b[:, 2]))
    return out


def transporte_covariante(metrica, tr, e1=None, e2=None):
    """Derivada covariante da base transversal ao longo do raio, em l.

    De/dl = de/dl + Gamma(p, e), com p^mu = dx^mu/dl (p^t = e^{-nu/2}: ao longo
    do raio nulo, o tempo local do observador estático anda junto com l). O
    transporte paralelo não depende da parametrização, então l serve.

    Devolve, ao longo do raio:
      De1_max   maior componente ORTONORMAL de De1/dl (tem de ser 0);
      omega     g(e1, De2/dl), a rotação no plano transversal (rad/km);
    e2 transversal pode ganhar componente ao longo de p sem que isso seja
    rotação; com De1 = 0 e a base ortonormal e ortogonal a p, De2 é
    proporcional a p, e omega é o que sobra para medir.

    A derivada de e é numérica (diferença central de 2a ordem) sobre a
    amostragem do raio: a resolução é ~h^2 e piora onde lambda' salta.
    """
    e1 = tr["e1"] if e1 is None else e1
    e2 = tr["e2"] if e2 is None else e2
    x, l = tr["x"], tr["l"]
    p, r, th = _coord(x, tr["k"], metrica)
    p[:, 0] = np.exp(-0.5 * metrica.nu(r))
    E1, _, _ = _coord(x, e1, metrica)
    E2, _, _ = _coord(x, e2, metrica)
    dE1 = np.gradient(E1, l, axis=0, edge_order=2) + _christoffel(metrica, r, th, p, E1)
    dE2 = np.gradient(E2, l, axis=0, edge_order=2) + _christoffel(metrica, r, th, p, E2)
    g = np.stack([-np.exp(metrica.nu(r)), np.exp(metrica.lam(r)), r * r,
                  (r * np.sin(th)) ** 2], axis=1)
    # componentes ortonormais de De1: e^{nu/2}, e^{lambda/2}, r, r sin(theta)
    escala = np.sqrt(np.abs(g))
    return dict(De1_max=np.max(np.abs(dE1 * escala), axis=1),
                omega=np.sum(g * E1 * dE2, axis=1),
                r=r)
