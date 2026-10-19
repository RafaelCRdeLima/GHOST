"""O GHOST como prescricao de transformacao de sabor no SNEWPY.

POR QUE ISTO EXISTE
-------------------
O pipeline padrao da comunidade (SNEWPY -> SNOwGLoBES -> taxas no detector) tem
uma lacuna estrutural no passo de oscilacao. A interface `FlavorTransformation`
recebe `(t, E)`, mas das dez prescricoes implementadas na 1.5.1 **nenhuma depende
do tempo**, e as duas que dependem da energia sao fisica alem do Modelo Padrao
(decaimento e decoerencia quantica). Para MSW o pipeline oferece duas constantes:

    AdiabaticMSW      P_ee = sin^2(theta13)          (P_H = 0)
    NonAdiabaticMSWH  o outro limite                 (P_H = 1)

A passagem do choque pela casca da ressonancia H, entre 1 e 10 s, faz `P_H` varrer
o intervalo entre esses dois limites em funcao de `t` E de `E`. Isso nao cabe em
constante nenhuma, e e justamente o regime em que o DUNE — o unico detector com
sensibilidade dominante ao canal `nu_e` — teria algo a ver.

Esta classe propaga de verdade: para cada `(t, E)` roda o GHOST atraves do perfil
de densidade daquele instante e devolve a probabilidade **incoerente**, que e o
observavel correto para fonte astrofisica (os pacotes dos autoestados de massa se
separam muito antes de sair da galaxia, entao a fase que define o P coerente ja
nao existe no detector).

USO
---
    from snewpy_ghost import GhostMSW, perfil_lei_potencia

    xform = GhostMSW(perfil_lei_potencia(), energias_MeV=[5,10,20,40],
                     tempos_s=[0.0])
    xform.prob_ee(0.0, 20.0)      # -> P(nu_e -> nu_e) propagado

O perfil e um chamavel `t_s -> (r_km, rho_g_cm3)`. Perfis estaticos podem ser
passados como a propria tupla.
"""
import hashlib
import json
import os
import subprocess
import tempfile

import numpy as np

from snewpy.flavor_transformation import FlavorTransformation, MassHierarchy, MixingParameters

BINARIO_PADRAO = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "native", "build", "bin", "ghost")


def identidade_binario(caminho):
    """(hash do nucleo, sha256 do executavel) do binario que vai propagar.

    Entra na chave do cache. Sem isso uma mudanca no C reaproveita
    probabilidades antigas sem aviso nenhum: o hash que cada CSV carrega nao
    resolve o problema se a execucao foi PULADA por cache, porque nesse caso
    nenhum CSV novo e' escrito.

    `ghost --version` sai com codigo 2 (imprime o uso junto), entao o
    returncode nao e' conferido; o que importa e' a linha do nucleo.
    """
    nucleo = "?"
    try:
        p = subprocess.run([caminho, "--version"],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for l in p.stdout.decode(errors="replace").splitlines():
            if "nucleo" in l:
                nucleo = l.split("nucleo")[1].strip(" ()\t")
                break
    except OSError:
        pass
    h = hashlib.sha256()
    try:
        with open(caminho, "rb") as f:
            for b in iter(lambda: f.read(1 << 20), b""):
                h.update(b)
        exe = h.hexdigest()[:16]
    except OSError:
        exe = "?"
    return nucleo, exe


def versoes_pacotes():
    """Versoes que podem mudar o resultado sem mudar nada no nosso codigo."""
    v = {}
    for nome in ("snewpy", "snowglobes", "numpy", "astropy"):
        try:
            mod = __import__(nome)
            v[nome] = str(getattr(mod, "__version__", "?"))
        except Exception:
            v[nome] = "ausente"
    return v


def _em_graus(x):
    """Angulo -> graus. O SNEWPY usa astropy Quantity; numeros nus sao radianos."""
    if hasattr(x, "to"):
        return float(x.to("deg").value)
    return float(np.degrees(x))


# ----------------------------------------------------------------- perfis ---

def perfil_lei_potencia(rho_ref=1.0e10, r_ref=100.0, expoente=3.0,
                        r0=1.0e3, r1=5.0e5, n=700):
    """Perfil suave `rho ~ r^-3`, o de referencia da literatura de turbulencia.

    Devolve um chamavel que ignora `t`: e o caso adiabatico puro, e serve de
    portao de validacao contra a `AdiabaticMSW` do SNEWPY.

    O alcance importa e nao e arbitrario. O perfil tem de comecar acima da
    ressonancia H e terminar bem abaixo da L PARA TODA a faixa de energia usada
    — senao o estado nao termina de converter e `P_ee` sai errado, com desvio
    que CRESCE com a energia. Com `rho ~ r^-3`:

        rho_H = 1,6e3 g/cm3 (20 MeV)  ->  r = 1,85e4 km
        rho_L =   20 g/cm3  (20 MeV)  ->  r = 7,9e4  km

    e ambas escalam como `1/E`, entao a 80 MeV a L ja esta em `rho ~ 5` g/cm3.
    O padrao termina em `5e5` km, onde `rho = 8e-2` g/cm3 — duas ordens de folga.

    O inicio tem a mesma exigencia do outro lado, e ela morde na energia BAIXA:
    a formula adiabatica supoe producao em densidade infinita, e comecar so um
    fator ~12 acima da ressonancia deixa o estado inicial longe do autoestado de
    materia. A 5 MeV a ressonancia H sobe para `rho = 6,3e3` g/cm3; comecando em
    `1e3` km (`rho = 1e7`) a folga volta a ser 1600.
    """
    r = np.geomspace(r0, r1, n)
    rho = rho_ref * (r_ref / r) ** expoente
    return lambda t: (r, rho)


# ------------------------------------------------------------- propagador ---

class GhostMSW(FlavorTransformation):
    """Transformacao de sabor calculada propagando, nao por formula fechada.

    Parameters
    ----------
    perfil : callable ou tupla
        `t_s -> (r_km, rho_g_cm3)`. Uma tupla e tratada como perfil estatico.
    mh : MassHierarchy
        Ordenamento de massa.
    mix_angles : tupla ou None
        `(theta12, theta13, theta23)` em radianos; None usa os do SNEWPY.
    energias_MeV, tempos_s : sequencias
        A grade em que o GHOST e de fato rodado. Consultas fora dela sao
        interpoladas (log-log em energia, linear em tempo) e **seguradas** nas
        bordas em vez de extrapoladas.
    Ye : float
        Fracao de eletrons.
    tolerancia : float
        Tolerancia da propagacao adaptativa por fatias.
    cache : str ou None
        Arquivo JSON onde a grade calculada e guardada e relida.
    """

    def __init__(self, perfil, mh=MassHierarchy.NORMAL, mix_angles=None,
                 energias_MeV=(5, 10, 20, 30, 50, 80), tempos_s=(0.0,),
                 Ye=0.5, tolerancia=1e-6, nt=300, binario=BINARIO_PADRAO,
                 cache=None, verbose=False, so_neutrino=False):
        if not isinstance(mh, MassHierarchy):
            raise TypeError("mh deve ser MassHierarchy")
        self.mass_order = mh
        self.perfil = perfil if callable(perfil) else (lambda t: perfil)
        self.Ye, self.tol, self.nt = float(Ye), float(tolerancia), int(nt)
        self.binario, self.verbose = binario, verbose

        pars = MixingParameters(mh)
        if mix_angles is not None:
            th12, th13, th23 = mix_angles
        else:
            th12, th13, th23 = pars.get_mixing_angles()
        # ATENCAO: o SNEWPY devolve Quantity em GRAUS. Aplicar np.degrees ao
        # .value multiplica por 57,3 e destroi tudo em silencio.
        self.th = [_em_graus(x) for x in (th12, th13, th23)]
        self.dm21 = float(pars.dm21_2.to("eV2").value)
        self.dm31 = float(pars.dm31_2.to("eV2").value)
        self.dcp = float(_em_graus(pars.deltaCP)) if hasattr(pars, "deltaCP") else -90.0
        # os mesmos De_i da AdiabaticMSW, para o canal que nao propagamos
        c13, s13 = np.cos(np.radians(self.th[1])), np.sin(np.radians(self.th[1]))
        c12, s12 = np.cos(np.radians(self.th[0])), np.sin(np.radians(self.th[0]))
        self.De1, self.De2, self.De3 = (c12 * c13) ** 2, (s12 * c13) ** 2, s13 ** 2

        self.E = np.asarray(sorted(float(e) for e in energias_MeV), dtype=float)
        self.t = np.asarray(sorted(float(x) for x in tempos_s), dtype=float)
        self._cache = cache
        self._grade = None
        self._proveniencia = {}
        # `so_neutrino` preenche o canal de antineutrino com a constante
        # adiabatica em vez de propagar. Corta o custo pela metade e e legitimo
        # QUANDO o observavel nao depende dele — o caso do canal nue_Ar40, que e
        # alimentado so por prob_ee e prob_ex. Nao usar fora disso.
        self.so_neutrino = bool(so_neutrino)

    # -- a grade ------------------------------------------------------------

    def _chave(self):
        """Chave do cache: fisica, perfil E A IDENTIDADE DO QUE VAI RODAR.

        Na CHAVE entram o hash do nucleo e o sha256 do executavel: a grade
        e' P_ee calculado pelo C, e uma grade cacheada por uma versao anterior
        do nucleo e' resultado de outro programa. As versoes dos pacotes NAO
        entram na chave — nao tocam nesta grade, e um bump de patch do numpy
        invalidaria tudo sem motivo — mas ficam registradas na proveniencia,
        porque afetam o fluxo e as contagens adiante.
        """
        m = hashlib.sha256()
        nucleo, exe = identidade_binario(self.binario)
        self._proveniencia = dict(nucleo=nucleo, executavel=exe,
                                  **versoes_pacotes())
        m.update(json.dumps([self.th, self.dm21, self.dm31, self.dcp,
                             int(self.so_neutrino),
                             self.Ye, self.tol, self.nt,
                             list(self.E), list(self.t),
                             int(self.mass_order == MassHierarchy.NORMAL),
                             nucleo, exe],
                            sort_keys=True).encode())
        for tt in self.t:
            r, rho = self.perfil(tt)
            m.update(np.asarray(r, dtype=float).tobytes())
            m.update(np.asarray(rho, dtype=float).tobytes())
        return m.hexdigest()[:16]

    def _rodar(self, r, rho, E_MeV, antinu):
        spec = {
            "format_version": 3, "seed": 1,
            "neutrino": {"n_flavors": 3,
                         "theta12_deg": self.th[0], "theta13_deg": self.th[1],
                         "theta23_deg": self.th[2], "delta_cp_deg": self.dcp,
                         "dm21_sq_eV2": self.dm21,
                         "dm31_sq_eV2": self.dm31,
                         "energy_MeV": float(E_MeV), "flavor": "e",
                         "antineutrino": bool(antinu)},
            "matter": {"profile": "table", "Ye": self.Ye, "tolerance": self.tol,
                       "r_start_km": float(r[0]),
                       "r_km": [float(x) for x in r],
                       "rho_table": [float(x) for x in rho]},
            "time": {"t1": float(r[-1] - r[0]), "nt": self.nt},
            "observables": {"target": 0, "population": False},
            "realizations": 1,
        }
        with tempfile.TemporaryDirectory() as d:
            js, csv = os.path.join(d, "s.json"), os.path.join(d, "s.csv")
            json.dump(spec, open(js, "w"))
            p = subprocess.run([self.binario, "run", js, "--saida", csv],
                               capture_output=True)
            if p.returncode:
                raise RuntimeError("GHOST falhou: " + p.stderr.decode()[-400:])
            linhas = [l for l in open(csv).read().splitlines()
                      if not l.startswith("#!")]
        col = linhas[0].split(",")
        return float(linhas[-1].split(",")[col.index("inc_e")])

    def grade(self):
        """Roda o GHOST em toda a grade `(t, E)`, para nu e antinu. Cacheia."""
        if self._grade is not None:
            return self._grade
        if self._cache and os.path.exists(self._cache):
            d = json.load(open(self._cache))
            if d.get("chave") == self._chave():
                self._grade = (np.array(d["Pee"]), np.array(d["Peebar"]))
                if self.verbose:
                    q = d.get("proveniencia", {})
                    print("  cache: nucleo %s, snewpy %s"
                          % (q.get("nucleo", "?"), q.get("snewpy", "?")))
                return self._grade
        Pee = np.empty((self.t.size, self.E.size))
        Peb = np.empty_like(Pee)
        for i, tt in enumerate(self.t):
            r, rho = self.perfil(tt)
            for j, ee in enumerate(self.E):
                Pee[i, j] = self._rodar(r, rho, ee, False)
                Peb[i, j] = (self.De1 if self.mass_order == MassHierarchy.NORMAL
                             else self.De3) if self.so_neutrino \
                    else self._rodar(r, rho, ee, True)
                if self.verbose:
                    print(f"  t={tt:8.4f} s  E={ee:6.1f} MeV   "
                          f"P_ee={Pee[i,j]:.6f}  P_eebar={Peb[i,j]:.6f}",
                          flush=True)
        self._grade = (Pee, Peb)
        if self._cache:
            # a proveniencia vai LEGIVEL junto da chave: a chave e' um
            # digest e nao diz qual nucleo produziu a grade.
            json.dump({"chave": self._chave(),
                       "proveniencia": self._proveniencia,
                       "Pee": Pee.tolist(),
                       "Peebar": Peb.tolist()}, open(self._cache, "w"))
        return self._grade

    # -- consulta -----------------------------------------------------------

    @staticmethod
    def _nu(x):
        return np.asarray(getattr(x, "value", x), dtype=float)

    def _consulta(self, M, t, E):
        ta, Ea = self._nu(t), self._nu(E)
        te, Ee = np.atleast_1d(ta), np.atleast_1d(Ea)
        # log-log em energia, segurando nas bordas; linear em tempo
        def em_t(k):
            return np.exp(np.interp(np.log(np.clip(Ee, self.E[0], self.E[-1])),
                                    np.log(self.E), np.log(M[k])))
        if self.t.size == 1:
            saida = em_t(0)
        else:
            tc = np.clip(te, self.t[0], self.t[-1])
            porE = np.array([np.interp(tc, self.t, M[:, j])
                             for j in range(self.E.size)]).T   # (Nt, NE)
            saida = np.array([
                np.exp(np.interp(np.log(np.clip(Ee, self.E[0], self.E[-1])),
                                 np.log(self.E), np.log(porE[i])))
                for i in range(tc.size)])
            if te.size == 1:
                saida = saida[0]
        if ta.ndim == 0 and Ea.ndim == 0:
            return float(np.ravel(saida)[0])
        return np.squeeze(saida)

    def prob_ee(self, t, E):
        return self._consulta(self.grade()[0], t, E)

    def prob_eebar(self, t, E):
        return self._consulta(self.grade()[1], t, E)

    def prob_ex(self, t, E):
        return 1. - self.prob_ee(t, E)

    def prob_xe(self, t, E):
        return (1. - self.prob_ee(t, E)) / 2.

    def prob_xx(self, t, E):
        return (1. + self.prob_ee(t, E)) / 2.

    def prob_exbar(self, t, E):
        return 1. - self.prob_eebar(t, E)

    def prob_xebar(self, t, E):
        return (1. - self.prob_eebar(t, E)) / 2.

    def prob_xxbar(self, t, E):
        return (1. + self.prob_eebar(t, E)) / 2.


def perfil_choque(rho_ref=1.0e10, r_ref=100.0, expoente=3.0,
                  r0=6.0e2, r1=5.0e5, v_frente=2.0e4, compressao=10.0,
                  rarefacao=0.15, espessura=0.15, largura=None,
                  n_global=700, n_fino=1500):
    """Progenitor em lei de potencia com uma frente de choque que se propaga.

    Parametrizacao no espirito de Schirato & Fuller (2002) e Fogli et al. (2003):
    a montante do choque o perfil e o do progenitor; atras da frente ha uma casca
    comprimida; mais para dentro, a regiao rarefeita deixada pela expansao.
    `r_f(t) = v_frente * t`.

    A LARGURA DA FRENTE e o parametro delicado: `gamma ~ H`, entao `P_H`
    depende criticamente da altura de escala no ponto de ressonancia. Um degrau
    resolvido pela malha tem largura igual ao passo da malha — numero
    arbitrario, exatamente como as ~4 zonas por que um esquema numerico de
    captura de choque o espalha: em simulacao, a largura do choque e definida
    pela numerica, nao pela fisica.

    Fisicamente, dentro da estrela, o choque e mediado por radiacao (Levinson &
    Nakar 2020) e tem largura de ordem de milimetros, cerca de oito ordens de
    grandeza abaixo de `L_osc = 67` km. O criterio que importa e
    `largura << L_osc`: abaixo disso a travessia e SUBITA e `P_H -> 1`, um
    limite bem definido e maximo, e o resultado deixa de depender da largura. `largura=None` usa o degrau puro; um
    valor em km suaviza as duas bordas com `tanh` e adensa a malha para resolver.
    """
    def perfil(t):
        rf = max(v_frente * float(t), 1.05 * r0)
        r = np.geomspace(r0, r1, n_global)
        if r0 < rf < r1:
            if largura is None:
                fino = np.linspace(max(r0, 0.7 * rf), min(r1, 1.3 * rf), n_fino)
            else:
                # resolve CADA borda com passo ~ largura/20
                w = float(largura)
                fino = np.concatenate([
                    np.linspace(max(r0, b - 12 * w), min(r1, b + 12 * w), 500)
                    for b in ((1.0 - espessura) * rf, rf)]
                    + [np.linspace(max(r0, 0.7 * rf), min(r1, 1.3 * rf), 400)])
            r = np.unique(np.concatenate([r, fino]))
            # junta as duas grades sem deixar pontos quase coincidentes, que
            # estragariam a interpolacao log-log do nucleo
            manter = np.concatenate([[True], np.diff(np.log(r)) > 1e-9])
            r = r[manter]
        rho = rho_ref * (r_ref / r) ** expoente
        r_in = (1.0 - espessura) * rf
        if largura is None:
            x = r / rf
            fator = np.where(x > 1.0, 1.0,
                             np.where(x > 1.0 - espessura, compressao, rarefacao))
        else:
            w = float(largura)
            A = 0.5 * (1.0 + np.tanh((r - r_in) / w))    # 0 dentro, 1 fora de r_in
            B = 0.5 * (1.0 + np.tanh((r - rf) / w))      # 0 dentro, 1 fora de r_f
            fator = rarefacao * (1.0 - A) + compressao * A * (1.0 - B) + B
        return r, rho * fator
    return perfil
