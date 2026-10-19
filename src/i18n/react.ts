import { useSyncExternalStore } from 'react';
import { lingua, observarLingua, t as traduz, type Locale, type Params } from './index.ts';

export const useLingua = (): Locale =>
  useSyncExternalStore(observarLingua, lingua, () => 'pt' as Locale);

/**
 * A função de tradução amarrada à língua atual.
 *
 * Recebe a língua por parâmetro em vez de ler a global: assim o React sabe que
 * o componente depende dela e redesenha quando ela muda.
 */
export const useT = (): ((chave: string, params?: Params) => string) => {
  const l = useLingua();
  return (chave, params) => traduz(chave, params, l);
};
