#!/bin/sh
# teste_cli.sh — o CLI, que os testes de aceitação NÃO cobrem.
#
# A suíte de `native/tests/` exercita o núcleo. O CLI é casca, e a casca pode
# divergir entre subcomandos sobre o que é um grafo válido: `run` passa por
# dae_run, que trata `edges`; um `traj` que montasse o grafo sozinho e só
# chamasse dae_graph_build rejeitaria DAE_G_EDGELIST, com
# "grafo: argumento invalido" num arquivo que `run` aceita sem reclamar.
#
# Um bug de casca não aparece em nenhum teste de núcleo, por melhor que ele
# seja. Este arquivo existe para essa classe.
set -eu
RAIZ=$(cd "$(dirname "$0")/.." && pwd)
BIN="$RAIZ/native/build/bin/ghost"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT INT TERM

falhou=0
ok()   { printf '   . %s\n' "$1"; }
ruim() { printf '   X %s\n' "$1"; falhou=1; }

# Dois sabores em matéria: auto-laços na diagonal, mistura na aresta.
cat > "$TMP/edgelist.json" <<'JSON'
{ "format_version": 2, "seed": 1,
  "graph": { "generator": "edgelist", "n": 2,
             "edges": [[0,0,0.9397], [1,1,-0.9397], [0,1,-0.3420]] },
  "hamiltonian": { "kind": "adjacency", "gamma": 1.0, "normalization": "none" },
  "initial": { "site": 0 },
  "time": { "t1": 6.2832, "nt": 20 },
  "observables": { "target": 1, "population": false },
  "realizations": 1,
  "trajectories": { "n_traj": 20, "gamma_deph": 0.05, "rho_stride": 5,
                    "output_mode": "accumulate_rho" } }
JSON

cat > "$TMP/gerador.json" <<'JSON'
{ "format_version": 2, "seed": 3,
  "graph": { "generator": "cycle", "params": { "n": 16 } },
  "hamiltonian": { "kind": "adjacency", "gamma": 1.0, "normalization": "none" },
  "initial": { "site": 0 },
  "time": { "t1": 4.0, "nt": 8 },
  "observables": { "population": false }, "realizations": 1,
  "trajectories": { "n_traj": 10, "gamma_deph": 0.05, "rho_stride": 4,
                    "output_mode": "accumulate_rho" } }
JSON

# Aresta apontando para fora do grafo: dae_csr_from_edges recusa com
# DAE_ERR_PARAM. Serve de anti-vacuidade — e tambem de segunda metade da
# pergunta: run e traj precisam concordar no que RECUSAM, nao so no que aceitam.
#
# (Nao serve `edges` sem `n`: o default de dae_gen_params_default e 64, entao
#  run aceita com n=64, e traj aceitar tambem e o comportamento CERTO. Testar
#  isso cobraria uma promessa que o programa nao faz.)
cat > "$TMP/fora.json" <<'JSON'
{ "format_version": 2, "seed": 1,
  "graph": { "generator": "edgelist", "n": 2, "edges": [[0,99,1.0]] },
  "hamiltonian": { "kind": "adjacency", "gamma": 1.0, "normalization": "none" },
  "initial": { "site": 0 }, "time": { "t1": 1.0, "nt": 2 },
  "observables": { "population": false }, "realizations": 1,
  "trajectories": { "n_traj": 4, "gamma_deph": 0.1, "rho_stride": 1,
                    "output_mode": "accumulate_rho" } }
JSON

echo "== CLI: run e traj concordam sobre o que e um grafo"

# 1. `run` com edgelist — serve de controle.
if "$BIN" run "$TMP/edgelist.json" --saida "$TMP/a.csv" >/dev/null 2>&1; then
  ok "run aceita edgelist"
else
  ruim "run REJEITOU edgelist (controle quebrado)"
fi

# 2. `traj` com edgelist — tem de aceitar o mesmo que `run` (regressao).
if "$BIN" traj "$TMP/edgelist.json" --saida "$TMP/b.csv" >/dev/null 2>&1; then
  ok "traj aceita edgelist"
else
  ruim "traj REJEITOU edgelist — a divergencia com run voltou"
fi

# 3. A saída de traj é utilizável, não só um arquivo que existe.
if [ -s "$TMP/b.csv" ] \
   && grep -q '^#! core_hash ' "$TMP/b.csv" \
   && grep -q '^s,t,i,j,re,im$' "$TMP/b.csv"; then
  ok "o rho de traj tem procedencia e colunas"
else
  ruim "a saida de traj nao tem cabecalho de procedencia ou colunas esperadas"
fi

# 4. `traj` com gerador continua funcionando.
if "$BIN" traj "$TMP/gerador.json" --saida "$TMP/c.csv" >/dev/null 2>&1; then
  ok "traj continua aceitando gerador"
else
  ruim "traj QUEBROU para gerador — o fix derrubou o caminho antigo"
fi

# 5. ANTI-VACUIDADE: o teste enxerga um spec ruim?
#    Sem isto, um `traj` que aceitasse qualquer coisa passaria em 2 e 4.
if "$BIN" traj "$TMP/fora.json" --saida "$TMP/d.csv" >/dev/null 2>&1; then
  ruim "traj ACEITOU aresta fora do grafo — o guarda nao existe"
else
  ok "traj recusa aresta fora do grafo"
fi

# 6. E `run` recusa o MESMO arquivo? A concordancia vale nos dois sentidos.
if "$BIN" run "$TMP/fora.json" --saida "$TMP/e.csv" >/dev/null 2>&1; then
  ruim "run ACEITOU o que traj recusou — divergencia, na direcao oposta"
else
  ok "run recusa o mesmo arquivo que traj recusa"
fi

# 7. O espectro: `run` com bloco `scan` muda o formato de saida, e a suite de
#    aceitacao (que exercita o NUCLEO) nao ve isso — o laco de varredura mora
#    no CLI. Mesma classe de defeito de casca dos testes 2 a 6.
cat > "$TMP/scan.json" <<'JSON'
{ "format_version": 3, "seed": 1,
  "neutrino": { "n_flavors": 3, "theta12_deg": 33.4, "theta13_deg": 8.6,
                "theta23_deg": 49.0, "delta_cp_deg": -90.0,
                "dm21_sq_eV2": 7.42e-5, "dm31_sq_eV2": 2.51e-3,
                "energy_MeV": 2500.0, "flavor": "mu",
                "rho_g_cm3": 2.848, "Ye": 0.4957 },
  "scan": { "energy_min_MeV": 500.0, "energy_max_MeV": 5000.0,
            "points": 20, "log": false, "both_signs": true },
  "time": { "t1": 1300.0, "nt": 40 },
  "observables": { "target": 0, "population": false }, "realizations": 1 }
JSON

if "$BIN" run "$TMP/scan.json" --saida "$TMP/s.csv" >/dev/null 2>&1 \
   && grep -q '^E_MeV,p_e,.*,acp_tau$' "$TMP/s.csv" \
   && [ "$(grep -vc '^#!' "$TMP/s.csv")" -eq 21 ]; then
  ok "scan emite espectro com colunas nomeadas e 20 energias"
else
  ruim "o scan nao produziu o espectro esperado"
fi

# 8. ANTI-VACUIDADE do scan: as energias VARIAM de linha para linha.
#    Um laco que esquecesse de atualizar a energia daria 20 linhas identicas,
#    e o teste 7 passaria.
if [ "$(grep -v '^#!' "$TMP/s.csv" | tail -n +2 | cut -d, -f2 | sort -u | wc -l)" -gt 15 ]; then
  ok "as 20 energias dao probabilidades distintas"
else
  ruim "o espectro repete valores — a energia nao esta variando"
fi

if [ "$falhou" -ne 0 ]; then
  echo "   TESTE DE CLI FALHOU"
  exit 1
fi
echo "   PASSOU  8 verificacoes"
