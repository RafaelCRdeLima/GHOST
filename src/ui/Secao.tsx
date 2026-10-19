import { useState, type ReactNode } from 'react';

/* Secao.tsx — seção recolhível com resumo de uma linha.
 *
 * Recolhida, ela mostra o que está configurado ali dentro. Sem o resumo,
 * recolher esconderia o estado e obrigaria a abrir cada seção para lembrar o
 * que está valendo — trocaria rolagem por cliques, sem ganhar nada.
 */
export function Secao({ titulo, resumo, aberta = false, children }: {
  titulo: string;
  resumo: string;
  aberta?: boolean;
  children: ReactNode;
}) {
  const [ab, setAb] = useState(aberta);
  return (
    <section className={ab ? 'secao aberta' : 'secao'}>
      <button className="cabeca" onClick={() => setAb((v) => !v)} aria-expanded={ab}>
        <span className="seta">{ab ? '▾' : '▸'}</span>
        <span className="titulo">{titulo}</span>
        {!ab && <span className="resumo mono">{resumo}</span>}
      </button>
      {ab && <div className="corpo-secao">{children}</div>}
    </section>
  );
}
