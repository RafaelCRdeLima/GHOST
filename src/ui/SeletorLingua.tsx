import { LINGUAS_PRONTAS, NOME_DA_LINGUA, definirLingua } from '../i18n/index.ts';
import { useLingua } from '../i18n/react.ts';

/** Cada língua aparece NA PRÓPRIA língua: quem não lê a atual precisa se achar. */
export function SeletorLingua() {
  const atual = useLingua();
  return (
    <select className="lingua" value={atual} aria-label="idioma"
            onChange={(e) => definirLingua(e.target.value as never)}>
      {LINGUAS_PRONTAS.map((l) => (
        <option key={l} value={l}>{NOME_DA_LINGUA[l]}</option>
      ))}
    </select>
  );
}
