/* Metricas.tsx — o painel de leitura, junto ao mapa.
 *
 * São os números que se anota a cada varredura: λ₂, Q, ⟨d⟩, α, N, |E|. Ficam
 * separados do diagnóstico, que é outra coisa. Diagnóstico (norma,
 * procedência) responde "posso confiar nisto?" e fica permanente no rodapé;
 * estes respondem "o que eu medi?" e precisam de corpo legível.
 */
export interface Metrica {
  rotulo: string;
  valor: string;
  nota?: string;
  alerta?: boolean;
}

export function Metricas({ itens }: { itens: Metrica[] }) {
  return (
    <div className="metricas">
      {itens.map((m) => (
        <div key={m.rotulo} className={m.alerta ? 'metrica alerta' : 'metrica'}>
          <span className="rot">{m.rotulo}</span>
          <span className="val">{m.valor}</span>
          {m.nota && <span className="nota">{m.nota}</span>}
        </div>
      ))}
    </div>
  );
}
