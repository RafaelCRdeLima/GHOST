#!/usr/bin/env node
/* amalgamate.mjs — concatena core/ num único tradutor.
 *
 * É a peça que faz o teste 7 valer alguma coisa. O texto emitido aqui é o
 * MESMO que vai para dentro do WASM e para dentro do .cpp exportado, então o
 * navegador e o cluster não rodam duas implementações que concordam: rodam o
 * mesmo código sob dois compiladores.
 *
 * Uso: node tools/amalgamate.mjs [saida.c]
 *      node tools/amalgamate.mjs --ts    (gera src/export/recursos.gerado.ts)
 *
 * O modo --ts embute o núcleo E os três templates num módulo TypeScript, para
 * que o emissor funcione igual no navegador (empacotado pelo Vite) e no Node.
 * Ler os arquivos do disco funcionaria só no Node, e aí o exportador do
 * navegador precisaria de uma segunda cópia dos templates.
 */
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const CORE = join(ROOT, 'core');

const files = readFileSync(join(CORE, 'amalgam.list'), 'utf8')
  .split('\n')
  .map((l) => l.trim())
  .filter((l) => l && !l.startsWith('#'));

/* Os #include internos somem: no tradutor único não há o que incluir. Os
   #include <...> do sistema ficam, e as guardas de inclusão tornam a
   repetição inofensiva. */
const INTERNAL = /^\s*#\s*include\s+"dae_[A-Za-z0-9_.]+"\s*$/;

function hashDoNucleo() {
  const h = readFileSync(join(CORE, 'dae_version.generated.h'), 'utf8');
  const m = h.match(/DAE_CORE_HASH\s+"([0-9a-f]+)"/);
  return m ? m[1] : '?';
}

const parts = [
  '/* GERADO por tools/amalgamate.mjs — não editar.',
  ' *',
  ' * Núcleo do Daedalus concatenado na ordem de core/amalgam.list. Compila',
  ' * como C99 e como C++17 a partir do mesmo texto.',
  ' */',
  '',
];
for (const f of files) {
  const body = readFileSync(join(CORE, f), 'utf8')
    .split('\n')
    .filter((l) => !INTERNAL.test(l))
    .join('\n');
  parts.push(`/* ===== ${f} ===== */`, body, '');
}

const out = parts.join('\n');

if (process.argv[2] === '--ts') {
  const lit = (t) => '`' + t.replace(/\\/g, '\\\\').replace(/`/g, '\\`')
                            .replace(/\$\{/g, '\\${') + '`';
  const tpl = (f) => readFileSync(join(ROOT, 'templates', f), 'utf8');
  const destino = join(ROOT, 'src', 'export', 'recursos.gerado.ts');
  const ts = `/* GERADO por tools/amalgamate.mjs --ts — não editar.
 *
 * O núcleo amalgamado e os templates dos exportadores, embutidos como texto.
 * É o MESMO texto que compila no WASM: o .cpp exportado não é uma reescrita.
 */
export const HASH_NUCLEO = ${JSON.stringify(hashDoNucleo())};
export const NUCLEO_C = ${lit(out)};
export const TPL_CPP = ${lit(tpl('cpp_main.tpl.cpp'))};
export const TPL_WL_PKG = ${lit(tpl('wolfram/Daedalus.wl'))};
export const TPL_WL_NB = ${lit(tpl('wolfram/DaedalusDemo.nb'))};
export const TPL_PY = ${lit(tpl('python.tpl.py'))};
`;
  mkdirSync(dirname(destino), { recursive: true });
  writeFileSync(destino, ts);
  console.error(`recursos: ${(ts.length / 1024).toFixed(0)} KB -> ${destino}`);
  process.exit(0);
}

const target = process.argv[2];
if (target) {
  writeFileSync(target, out);
  console.error(`amalgamado: ${files.length} arquivos, ${out.length} bytes -> ${target}`);
} else {
  process.stdout.write(out);
}
