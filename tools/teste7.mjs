#!/usr/bin/env node
/* teste7.mjs — ACEITAÇÃO 7: navegador contra .cpp exportado, mesmo spec.json.
 *
 * Três lados, um núcleo: o WASM (o que a página roda), o binário nativo e o
 * .cpp emitido pelo exportador. O .cpp REGENERA o grafo a partir do spec, o que
 * é o que torna este teste um teste do EMISSOR e não só do propagador — o
 * oráculo Wolfram não pode fazer isso, porque recebe o grafo pronto.
 *
 * VARREDURA EM alpha, NÃO EM dt. Um defeito na escolha da ordem de Chebyshev
 * pode ser não-monotônico em alpha = a·dt (CONVENTIONS.md 6.2), e a grade
 * exportada tem dt livre com `a` dependendo da normalização: grade fina com H
 * normalizada dá alpha pequeno, grade grossa sem normalização dá alpha grande.
 * Varrer dt e reportar dt esconderia isso.
 */
import { execFileSync } from 'node:child_process';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

const RAIZ = new URL('..', import.meta.url).pathname;
const TMP = process.env.DAE_TMP ?? '/tmp/daedalus-teste7';
mkdirSync(TMP, { recursive: true });
process.env.DAE_DATA = '1970-01-01T00:00:00.000Z';

/* Regenera o nucleo amalgamado antes de emitir: sem isto, um core mexido e um
   recursos.gerado.ts velho fazem o teste comparar duas versoes diferentes do
   mesmo programa. */
execFileSync(process.execPath, [join(RAIZ, 'tools/amalgamate.mjs'), '--ts'], { stdio: 'ignore' });

const { Daedalus } = await import('../wasm/build/daedalus.mjs');
const dae = await Daedalus.criar();

/* --- casos: os specs do repositório mais uma varredura de alpha --- */
const base = JSON.parse(readFileSync(join(RAIZ, 'specs/microtubulo.json'), 'utf8'));
const casos = [
  ['microtubulo religado', base],
  ['linha longa', {
    ...base,
    graph: { generator: 'path', params: { n: 401 } },
    hamiltonian: { kind: 'adjacency', gamma: 1, normalization: 'none', lanczos_steps: 0 },
    initial: { site: 200 }, time: { t1: 30, nt: 30 },
    observables: { target: 260, module_concurrence: false },
  }],
  ['SBM laplaciana', {
    ...base, seed: 99,
    graph: { generator: 'sbm', params: { n: 240, n_modules: 4, p_in: 0.3, p_out: 0.01 } },
    hamiltonian: { kind: 'laplacian', gamma: 1, normalization: 'mean_degree', lanczos_steps: 40 },
    initial: { site: 0 }, time: { t1: 25, nt: 40 },
    observables: { target: 239, module_concurrence: true },
  }],
  ['lista de arestas explicita', {
    format_version: 2, seed: 7,
    graph: { generator: 'edgelist', n: 6,
             edges: [[0, 1, 1.0], [1, 2, 1.0], [2, 3, 0.5], [3, 4, 1.0], [4, 5, 1.0], [5, 0, 0.5]] },
    hamiltonian: { kind: 'adjacency', gamma: 1, normalization: 'none', lanczos_steps: 0 },
    initial: { site: 0 }, time: { t1: 6, nt: 20 },
    observables: { target: 3, module_concurrence: true },
  }],
];

/* alpha = a·dt; com normalizacao espectral a ~ 1, entao t1/nt controla alpha. */
for (const alvo of [0.2, 5, 25, 60, 150]) {
  casos.push([`alpha ~ ${alvo}`, {
    ...base,
    graph: { generator: 'microtubule',
             params: { n_par: 40, n_perp: 13, seam_shift: 3, n_modules: 4, ws_p: 0 } },
    hamiltonian: { kind: 'adjacency', gamma: 1, normalization: 'spectral', lanczos_steps: 40 },
    initial: { site: 260 }, time: { t1: alvo * 20, nt: 20 },
    observables: { target: 0, module_concurrence: true },
  }]);
}

const alphas = [];
let falhas = 0;

for (const [nome, spec] of casos) {
  const specPath = join(TMP, 'spec.json');
  writeFileSync(specPath, JSON.stringify(spec, null, 2));

  /* 1. WASM — é o que a página roda */
  const rede = dae.carregar(readFileSync(specPath, 'utf8'));
  dae.avancar(rede.nt);
  const csvWasm = join(TMP, 'wasm.csv');
  writeFileSync(csvWasm, dae.csv(true));
  alphas.push(rede.alpha);

  /* 2. nativo */
  const csvNat = join(TMP, 'nativo.csv');
  execFileSync(join(RAIZ, 'native/build/bin/ghost'),
               ['run', specPath, '--estado', '--saida', csvNat]);

  /* 3. .cpp exportado: emitido, compilado, executado */
  const cpp = join(TMP, 'run.cpp'), bin = join(TMP, 'run'), csvCpp = join(TMP, 'cpp.csv');
  execFileSync(process.execPath, [join(RAIZ, 'tools/exportar.mjs'), specPath, 'cpp', cpp],
               { stdio: 'ignore' });
  execFileSync('g++', ['-O2', '-fopenmp', cpp, '-o', bin]);
  execFileSync(bin, ['--estado', '--saida', csvCpp], { stdio: ['ignore', 'ignore', 'ignore'] });

  console.log(`\n== ${nome}   n=${rede.n}  alpha=${rede.alpha.toFixed(4)}  digital=${rede.fingerprint}`);
  for (const [rot, a, b] of [['navegador x .cpp exportado', csvWasm, csvCpp],
                             ['navegador x nativo', csvWasm, csvNat]]) {
    console.log(`  ${rot}`);
    try {
      execFileSync(process.execPath, [join(RAIZ, 'tools/compare_csv.mjs'), a, b, '1e-14'],
                   { stdio: 'inherit' });
    } catch { ++falhas; }
  }
}

/* Anti-vacuidade da própria varredura: se todos os casos caíssem no mesmo
   alpha, a varredura seria decorativa e um defeito dependente de alpha
   passaria intacto. */
const lo = Math.min(...alphas), hi = Math.max(...alphas);
console.log(`\n  alpha coberto: ${lo.toFixed(3)} a ${hi.toFixed(1)}  (${alphas.length} casos)`);
if (!(hi / Math.max(lo, 1e-12) > 100)) {
  console.error('  FALHOU: a varredura de alpha nao cobre faixa suficiente');
  ++falhas;
}
console.log(falhas === 0 ? '\n  aceitacao 7: PASSOU' : `\n  aceitacao 7: FALHOU (${falhas})`);
process.exit(falhas === 0 ? 0 : 1);
