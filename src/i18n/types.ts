/* types.ts — as quatro línguas e o contrato do catálogo.
 *
 * O núcleo NÃO fala nenhuma delas: ele devolve códigos e números, e quem vira
 * isso em frase é o catálogo, que vive aqui fora.
 */
export const LOCALES = ['pt', 'en', 'fr', 'it'] as const;
export type Locale = (typeof LOCALES)[number];

/** O nome de cada língua NA PRÓPRIA língua: quem não lê a atual precisa se achar. */
export const NOME_DA_LINGUA: Readonly<Record<Locale, string>> = {
  pt: 'Português', en: 'English', fr: 'Français', it: 'Italiano',
};

export type Params = Readonly<Record<string, string | number>>;

/**
 * Cada entrada é uma FUNÇÃO, não um texto com marcadores.
 *
 * Português, inglês, francês e italiano concordam em número e gênero de jeitos
 * diferentes. Com função, cada língua resolve a própria gramática no seu
 * arquivo, sem biblioteca de pluralização no meio.
 */
export type Catalog = Readonly<Record<string, (p: Params) => string>>;
