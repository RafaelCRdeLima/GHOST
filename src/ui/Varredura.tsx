/* Varredura.tsx — observável médio contra p, com barras de erro.
 *
 * O parâmetro de controle é p, mas o resultado deve poder ser reportado contra
 * λ₂ e Q (CONVENTIONS.md). Por isso o eixo x é rotulado com p E o ponto
 * carrega quantas realizações entraram na média: uma barra de erro sem n é
 * ornamento.
 */
export interface Ponto { p: number; media: number; desvio: number; n: number; }

export function Varredura({ pontos, rotuloX, rotuloY }: {
  pontos: Ponto[]; rotuloX: string; rotuloY: string;
}) {
  const L = 720, A = 200, ml = 58, mr = 12, mt = 12, mb = 30;
  const largura = L - ml - mr, altura = A - mt - mb;
  if (pontos.length === 0) return null;

  const xs = pontos.map((q) => q.p);
  const lo = Math.min(...pontos.map((q) => q.media - q.desvio), 0);
  const hi = Math.max(...pontos.map((q) => q.media + q.desvio));
  const faixa = hi - lo || 1;
  const px = (p: number) => ml + (largura * (p - Math.min(...xs)))
    / Math.max(1e-12, Math.max(...xs) - Math.min(...xs));
  const py = (v: number) => mt + altura - (altura * (v - lo)) / faixa;

  const caminho = pontos.map((q, i) => `${i ? 'L' : 'M'}${px(q.p).toFixed(1)},${py(q.media).toFixed(1)}`).join('');

  return (
    <svg viewBox={`0 0 ${L} ${A}`} className="series" role="img">
      <rect x={ml} y={mt} width={largura} height={altura} className="quadro" />
      <path d={caminho} fill="none" stroke="#8FC2E4" strokeWidth={1.6} />
      {pontos.map((q) => (
        <g key={q.p}>
          <line x1={px(q.p)} x2={px(q.p)} y1={py(q.media - q.desvio)} y2={py(q.media + q.desvio)}
                stroke="#2C77A8" strokeWidth={1.2} />
          <circle cx={px(q.p)} cy={py(q.media)} r={3} fill="#E5A83F">
            <title>p = {q.p.toFixed(3)} · {q.media.toPrecision(4)} ± {q.desvio.toPrecision(2)} (n = {q.n})</title>
          </circle>
        </g>
      ))}
      <text x={ml} y={A - 8} className="eixo">{Math.min(...xs).toFixed(2)}</text>
      <text x={ml + largura} y={A - 8} textAnchor="end" className="eixo">{Math.max(...xs).toFixed(2)}</text>
      <text x={ml + largura / 2} y={A - 8} textAnchor="middle" className="eixo">{rotuloX}</text>
      <text x={4} y={mt + 10} className="legenda">{rotuloY}</text>
      <text x={4} y={mt + 22} className="legenda">{hi.toPrecision(3)}</text>
      <text x={4} y={mt + altura} className="legenda">{lo.toPrecision(3)}</text>
    </svg>
  );
}
