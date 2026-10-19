#__DAE_BANNER_PY__

"""Reprodução em Python, só com numpy e scipy.

Como o .wl, este arquivo RECEBE a lista de arestas do núcleo em vez de gerar o
grafo: ver a nota longa no template Wolfram. Ele valida o propagador, não o
emissor.

    pip install numpy scipy
    python3 __ARQUIVO__ > saida.csv
"""
import sys

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

n = __DAE_N__
arestas = __DAE_ARESTAS__          # [[i, j, w], ...] com i, j em 0..n-1
tipo = "__DAE_HAM__"
gamma = __DAE_GAMMA__
escala = __DAE_ESCALA__            # medida pelo núcleo
sitio = __DAE_SITIO__
alvo = __DAE_ALVO__                # -1 = sem alvo
t1 = __DAE_T1__
nt = __DAE_NT__
modulos = np.array(__DAE_MODULOS__)
nmod = __DAE_NMOD__

li = np.array([e[0] for e in arestas] + [e[1] for e in arestas])
lj = np.array([e[1] for e in arestas] + [e[0] for e in arestas])
lw = np.array([e[2] for e in arestas] * 2, dtype=float)
A = csr_matrix((lw, (li, lj)), shape=(n, n))

if tipo == "adjacency":
    H = -A
else:
    H = diags(np.asarray(A.sum(axis=1)).ravel()) - A
H = (gamma / escala) * H

psi = np.zeros(n, dtype=complex)
psi[sitio] = 1.0
dt = t1 / nt

print("#! oracle python")
print("#! core_hash __DAE_HASH__")
print(f"#! n {n}")
print(f"#! nmod {nmod}")
print(f"#! dt {dt!r}")
print("#! spec __DAE_SPEC_LINHA__")

cab = "t,norm,ipr,coh_l1,p_target" + "".join(f",pmod{m}" for m in range(nmod))
print(cab)

for s in range(1, nt + 1):
    psi = expm_multiply(-1j * H * dt, psi)
    p = np.abs(psi) ** 2
    alvo_txt = repr(float(p[alvo])) if alvo >= 0 else "nan"
    pmods = [repr(float(p[modulos == m].sum())) for m in range(nmod)]
    print(",".join([repr(s * dt), repr(float(p.sum())), repr(float((p ** 2).sum())),
                    repr(float(np.abs(psi).sum() ** 2 - 1)), alvo_txt] + pmods))

print()
print("# estado final")
print("j,re,im")
for j in range(n):
    print(f"{j},{psi[j].real!r},{psi[j].imag!r}")
