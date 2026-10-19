import { en } from './en.ts';
import { fr } from './fr.ts';
import { it } from './it.ts';
import { pt } from './pt.ts';
import { LOCALES, type Catalog, type Locale, type Params } from './types.ts';

export * from './types.ts';

const CATALOGOS: Readonly<Record<Locale, Catalog>> = { pt, en, fr, it };

/** Só oferece o que está traduzido: oferecer uma língua que responde em outra
 *  é pior do que não oferecer. */
export const LINGUAS_PRONTAS: readonly Locale[] =
  LOCALES.filter((l) => Object.keys(CATALOGOS[l]).length > 0);

const PADRAO: Locale = 'pt';
const CHAVE = 'daedalus.lingua';

const ehLocale = (v: string): v is Locale =>
  (LINGUAS_PRONTAS as readonly string[]).includes(v);

/** A língua do navegador, caindo no padrão. `navigator.languages` vem em ordem
 *  de preferência e traz variantes regionais: compara-se só o prefixo. */
export const detectarLingua = (): Locale => {
  try {
    const guardada = globalThis.localStorage?.getItem(CHAVE);
    if (guardada && ehLocale(guardada)) return guardada;
  } catch { /* modo privado, armazenamento bloqueado */ }
  const preferidas = typeof navigator === 'undefined'
    ? [] : [...(navigator.languages ?? [navigator.language])];
  for (const tag of preferidas) {
    const base = (tag ?? '').toLowerCase().split('-')[0] ?? '';
    if (ehLocale(base)) return base;
  }
  return PADRAO;
};

let atual: Locale = detectarLingua();
const ouvintes = new Set<() => void>();

export const lingua = (): Locale => atual;

export const definirLingua = (l: Locale): void => {
  if (l === atual) return;
  atual = l;
  try { globalThis.localStorage?.setItem(CHAVE, l); } catch { /* idem */ }
  for (const f of ouvintes) f();
};

export const observarLingua = (f: () => void): (() => void) => {
  ouvintes.add(f);
  return () => { ouvintes.delete(f); };
};

/**
 * Chave ausente devolve a PRÓPRIA chave, e isso é deliberado: o modo de falha
 * de um catálogo é silencioso — uma frase que some, ou vira string vazia, e
 * ninguém nota. Devolver o identificador cru faz o buraco aparecer na tela.
 * O teste de paridade em catalog.test.ts existe para que isso nunca chegue lá.
 */
export const t = (chave: string, params: Params = {}, l: Locale = atual): string => {
  const entrada = CATALOGOS[l]?.[chave] ?? CATALOGOS[PADRAO]?.[chave];
  return entrada ? entrada(params) : chave;
};
