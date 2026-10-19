/* catalog.test.ts — o catálogo é a única parte do programa que fala, e os dois
 * modos de falha dele são silenciosos:
 *
 *   1. chave que existe numa língua e não em outra — a tela mostra o
 *      identificador cru, ou pior, uma string vazia;
 *   2. chave usada no código e ausente de todos os catálogos.
 *
 * Nos dois casos nada quebra e o programa continua rodando. São dez linhas de
 * teste e a CI pega antes de o usuário ver.
 */
import { readFileSync, readdirSync } from 'node:fs';
import { describe, expect, it as caso } from 'vitest';

import { en } from './en.ts';
import { fr } from './fr.ts';
import { it as catalogoItaliano } from './it.ts';
import { pt } from './pt.ts';
import { LINGUAS_PRONTAS, LOCALES, NOME_DA_LINGUA, t, type Locale } from './index.ts';

const CATALOGOS: Record<Locale, Record<string, unknown>> =
  { pt, en, fr, it: catalogoItaliano };

/** Todo `t('chave')` que aparece na interface, lido da fonte. */
const chavesUsadas = (): string[] => {
  const dirs = ['src/ui', 'src/nucleo'];
  const achadas = new Set<string>();
  for (const d of dirs) {
    for (const f of readdirSync(d)) {
      if (!/\.(ts|tsx)$/.test(f) || f.endsWith('.test.ts') || f.endsWith('.test.tsx')) continue;
      const texto = readFileSync(`${d}/${f}`, 'utf8');
      for (const m of texto.matchAll(/\bt\(\s*'([a-z0-9_]+)'/g)) achadas.add(m[1] as string);
    }
  }
  return [...achadas].sort();
};

describe('catalogo de linguas', () => {
  caso('as quatro linguas tem EXATAMENTE as mesmas chaves', () => {
    const referencia = Object.keys(pt).sort();
    for (const l of LOCALES) {
      expect(Object.keys(CATALOGOS[l]).sort(), l).toEqual(referencia);
    }
  });

  caso('o seletor oferece exatamente as linguas com catalogo', () => {
    expect([...LINGUAS_PRONTAS].sort()).toEqual(
      LOCALES.filter((l) => Object.keys(CATALOGOS[l]).length > 0).sort(),
    );
  });

  caso('toda chave usada na interface existe no catalogo', () => {
    const usadas = chavesUsadas();
    expect(usadas.length).toBeGreaterThan(10);   /* anti-vacuidade: achou mesmo? */
    for (const chave of usadas) {
      expect(Object.keys(pt), `chave usada e nao traduzida: ${chave}`).toContain(chave);
    }
  });

  caso('nenhuma traducao devolve string vazia', () => {
    for (const l of LOCALES) {
      for (const chave of Object.keys(CATALOGOS[l])) {
        const frase = t(chave, { j: 1, n: 2, m: 3, k: 4, i: 5, arquivo: 'a', atual: 'b',
                                 linhas: 6, erro: 'e' }, l);
        expect(frase.length, `${l}.${chave}`).toBeGreaterThan(0);
      }
    }
  });

  caso('cada lingua se nomeia na propria lingua', () => {
    for (const l of LOCALES) expect(NOME_DA_LINGUA[l].length).toBeGreaterThan(2);
    expect(new Set(Object.values(NOME_DA_LINGUA)).size).toBe(LOCALES.length);
  });
});
