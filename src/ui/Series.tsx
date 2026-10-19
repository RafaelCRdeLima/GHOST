/* Series.tsx — séries temporais com eixos de verdade.
 *
 * O eixo x é TEMPO (γt), não número de passo: passo é detalhe da discretização,
 * e comparar duas corridas com grades diferentes exige a mesma abscissa física.
 *
 * Dois eixos y, porque as grandezas não compartilham faixa: IPR vive em [0,1] e
 * a coerência ℓ₁ vai até N−1. Empilhar as duas num eixo só achataria uma delas
 * contra o zero. Cada eixo é rotulado NA COR da sua curva — sem isso, dois eixos
 * viram adivinhação.
 */
export interface Serie {
  nome: string;
  cor: string;
  valores: ArrayLike<number>;
  eixo: 'esq' | 'dir';
}

const faixaDe = (ss: Serie[], n: number): [number, number] => {
  let lo = Infinity, hi = -Infinity;
  for (const s of ss) {
    for (let i = 0; i < n; ++i) {
      const v = s.valores[i];
      if (!Number.isFinite(v)) continue;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
  }
  if (!Number.isFinite(lo)) return [0, 1];
  if (hi - lo < 1e-12) return [lo - 0.5, lo + 0.5];
  const folga = (hi - lo) * 0.06;
  /* Não descer abaixo de zero numa grandeza que não é negativa: um eixo que
     começa em −0,046 sugere que a curva poderia ir para lá, e IPR, coerência e
     probabilidade não podem. */
  return [lo >= 0 ? Math.max(0, lo - folga) : lo - folga, hi + folga];
};

const fmt = (v: number): string => {
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 1e4 || a < 1e-3) return v.toExponential(1);
  return v.toPrecision(3).replace(/\.?0+$/, '');
};

export function Series({ series, tempo, cursor, aoClicar, rotuloX }: {
  series: Serie[];
  tempo: ArrayLike<number>;
  cursor: number;
  aoClicar?: (passo: number) => void;
  rotuloX: string;
}) {
  const n = tempo.length;
  if (n === 0 || series.length === 0) return null;

  const L = 1000, A = 300, ml = 62, mr = 62, mt = 14, mb = 34;
  const largura = L - ml - mr, altura = A - mt - mb;

  const esq = series.filter((s) => s.eixo === 'esq');
  const dir = series.filter((s) => s.eixo === 'dir');
  const [loE, hiE] = faixaDe(esq, n);
  const [loD, hiD] = faixaDe(dir, n);

  const t0 = tempo[0] ?? 0, t1 = tempo[n - 1] ?? 1;
  const px = (i: number) => ml + (largura * ((tempo[i] ?? 0) - t0)) / Math.max(1e-12, t1 - t0);
  const py = (v: number, e: 'esq' | 'dir') => {
    const [lo, hi] = e === 'esq' ? [loE, hiE] : [loD, hiD];
    return mt + altura - (altura * (v - lo)) / Math.max(1e-12, hi - lo);
  };

  const caminho = (s: Serie) => {
    let d = '', abriu = false;
    for (let i = 0; i < n; ++i) {
      const v = s.valores[i];
      if (!Number.isFinite(v)) { abriu = false; continue; }  /* NaN: buraco, não zero */
      d += `${abriu ? 'L' : 'M'}${px(i).toFixed(1)},${py(v, s.eixo).toFixed(1)}`;
      abriu = true;
    }
    return d;
  };

  const marcasY = [0, 0.25, 0.5, 0.75, 1];
  const marcasX = [0, 0.25, 0.5, 0.75, 1];
  const xCursor = px(Math.max(0, Math.min(n - 1, cursor - 1)));
  const corE = esq[0]?.cor ?? '#9AA6B0';
  const corD = dir[0]?.cor ?? '#9AA6B0';

  return (
    <div className="painel-series">
      {/* legenda FORA da área de plotagem: dentro, ela cobre justamente a parte
          da curva que costuma interessar */}
      <div className="legenda-series">
        {series.map((s) => (
          <span key={s.nome} style={{ color: s.cor }}>
            <i style={{ background: s.cor }} />{s.nome}
            <em>{s.eixo === 'esq' ? '↤' : '↦'}</em>
          </span>
        ))}
      </div>
      <svg viewBox={`0 0 ${L} ${A}`} className="series" preserveAspectRatio="none"
           onClick={(e) => {
             if (!aoClicar) return;
             const r = (e.currentTarget as SVGSVGElement).getBoundingClientRect();
             const f = (((e.clientX - r.left) / r.width) * L - ml) / largura;
             aoClicar(Math.round(Math.min(1, Math.max(0, f)) * (n - 1)) + 1);
           }}>
        {marcasY.map((f) => (
          <line key={f} className="grade" x1={ml} x2={ml + largura}
                y1={mt + altura * (1 - f)} y2={mt + altura * (1 - f)} />
        ))}
        {marcasX.map((f) => (
          <line key={`x${f}`} className="grade" y1={mt} y2={mt + altura}
                x1={ml + largura * f} x2={ml + largura * f} />
        ))}
        <rect x={ml} y={mt} width={largura} height={altura} className="quadro" />

        {series.map((s) => (
          <path key={s.nome} d={caminho(s)} fill="none" stroke={s.cor} strokeWidth={1.8}
                vectorEffect="non-scaling-stroke" />
        ))}
        <line x1={xCursor} x2={xCursor} y1={mt} y2={mt + altura} className="cursor" />

        {/* eixo y esquerdo, na cor da sua curva */}
        {marcasY.map((f) => (
          <text key={`e${f}`} x={ml - 8} y={mt + altura * (1 - f) + 3.5}
                textAnchor="end" className="marca" fill={corE}>
            {fmt(loE + (hiE - loE) * f)}
          </text>
        ))}
        {dir.length > 0 && marcasY.map((f) => (
          <text key={`d${f}`} x={ml + largura + 8} y={mt + altura * (1 - f) + 3.5}
                textAnchor="start" className="marca" fill={corD}>
            {fmt(loD + (hiD - loD) * f)}
          </text>
        ))}
        {marcasX.map((f) => (
          <text key={`tx${f}`} x={ml + largura * f} y={A - 14} textAnchor="middle" className="marca">
            {fmt(t0 + (t1 - t0) * f)}
          </text>
        ))}
        <text x={ml + largura / 2} y={A - 2} textAnchor="middle" className="rotulo">{rotuloX}</text>
      </svg>
    </div>
  );
}
