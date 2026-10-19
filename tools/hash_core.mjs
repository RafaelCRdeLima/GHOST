#!/usr/bin/env node
/* hash_core.mjs — identidade do núcleo.
 *
 * O hash cobre todos os arquivos de core/amalgam.list EXCETO o próprio
 * dae_version.generated.h, que não pode conter o hash de si mesmo. Cada figura
 * de artigo precisa ser rastreável até o par (spec.json, core_hash), então
 * este número entra no cabeçalho de todo arquivo exportado e de todo CSV.
 */
import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const CORE = join(ROOT, 'core');
const GENERATED = 'dae_version.generated.h';
const VERSION = '0.1.0';

const files = readFileSync(join(CORE, 'amalgam.list'), 'utf8')
  .split('\n')
  .map((l) => l.trim())
  .filter((l) => l && !l.startsWith('#') && l !== GENERATED);

const h = createHash('sha256');
for (const f of files) {
  h.update(f, 'utf8');
  h.update('\0', 'utf8');
  h.update(readFileSync(join(CORE, f)));
}
const hash = h.digest('hex').slice(0, 16);

const out = `/* GERADO por tools/hash_core.mjs — não editar à mão.
 *
 * DAE_CORE_HASH cobre todos os fontes de core/amalgam.list, na ordem em que
 * a amalgamação os concatena. Um spec.json mais este hash determinam o
 * resultado bit a bit. */
#ifndef DAE_VERSION_GENERATED_H
#define DAE_VERSION_GENERATED_H

#define DAE_VERSION   "${VERSION}"
#define DAE_CORE_HASH "${hash}"

#endif /* DAE_VERSION_GENERATED_H */
`;

const target = join(CORE, GENERATED);
let previous = null;
try { previous = readFileSync(target, 'utf8'); } catch { /* primeira geração */ }
if (previous !== out) writeFileSync(target, out);
console.log(`ghost ${VERSION}  core ${hash}  (${files.length} arquivos)`);
