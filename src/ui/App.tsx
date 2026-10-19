import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { useT } from '../i18n/react.ts';
import { Heatmap } from '../gl/heatmap';
import { MAX_ARESTAS_DESENHADAS, Nuvem } from '../gl/nuvem';
import { ELEVACAO_PADRAO, Tubo } from '../gl/tubo';
import { compressaoParaCaber, coordenadasTubo } from '../nucleo/tubo.ts';
import { empacotarLattice, sitioNoTexel } from '../nucleo/indices';
import { gradienteCss } from '../nucleo/paleta';
import { bytesDosQuadros, passoParaCaber, passosDosQuadros, quadroDoPasso } from '../nucleo/quadros.ts';
import { origemImplementacao, reimportar, type Reimportado } from '../nucleo/reimportar.ts';
import { amostrasDeTempo, bytesAcumularRho, bytesArquivarPsi, formatarBytes,
         type ModoSaida } from '../nucleo/trajetorias.ts';
import { SeletorLingua } from './SeletorLingua';
import { Metricas } from './Metricas';
import { Secao } from './Secao';
import { Series, type Serie } from './Series';
import { Varredura, type Ponto } from './Varredura';

const GERADORES = ['microtubule', 'sbm', 'path', 'cycle', 'grid2d', 'hypercube', 'complete'];
const TETO_QUADROS = 6_000_000;

interface Rede {
  n: number; nnz: number; nmod: number; nPar: number; nPerp: number;
  arestasDescartadas: number; religacoesFalhas: number; fingerprint: bigint;
  geom: number;
  escala: number; dt: number; alpha: number;
  lambda2: number; lambda2Residuo: number; lambda2Convergiu: boolean;
  Q: number; grauMedio: number; caminhoMedio: number; caminhoExato: boolean;
  componentes: number; arestas: number; nt: number;
}

function baixar(nome: string, texto: string) {
  const url = URL.createObjectURL(new Blob([texto], { type: 'text/plain' }));
  const a = document.createElement('a');
  a.href = url; a.download = nome; a.click();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

export default function App() {
  const t = useT();
  /* A classificacao mora em reimportar.ts e tem teste proprio; aqui so o
     rotulo. As chaves ficam como literais para o teste de paridade do catalogo
     enxerga-las na varredura por t('...'). */
  const rotuloOrigem = (v: string | undefined): string => {
    const o = origemImplementacao(v);
    return o.tipo === 'c' ? t('impl_c')
         : o.tipo === 'wolfram' ? t('impl_wolfram')
         : o.tipo === 'literal' ? o.literal : t('impl_desconhecida');
  };
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const heatRef = useRef<Heatmap | null>(null);
  const nuvemRef = useRef<Nuvem | null>(null);
  const tuboRef = useRef<Tubo | null>(null);
  const canvasTubo = useRef<HTMLCanvasElement>(null);
  /* Orientação em REF, não em estado: arrastar dispara dezenas de eventos por
     segundo, e re-renderizar a interface inteira a cada um deles trava o giro
     justamente quando o usuário está tentando entender a forma. */
  const orientRef = useRef<[number, number]>([0, ELEVACAO_PADRAO]);
  const arrasto = useRef<{ x: number; y: number; moveu: boolean } | null>(null);
  const canvasNuvem = useRef<HTMLCanvasElement>(null);
  const arestasRef = useRef<Array<[number, number, number]> | null>(null);
  const workerRef = useRef<Worker | null>(null);
  /* Cada quadro carrega o PASSO em que foi tirado: sem isso o deslizador conta
     quadros e o rótulo mente sobre o tempo. */
  const quadrosRef = useRef<Array<{ pop: Float32Array; passo: number }>>([]);
  const desenharRef = useRef<(p: Float32Array) => void>(() => {});
  const redeRef = useRef<Rede | null>(null);
  const primeira = useRef(true);
  const pendenteCsv = useRef<{ nome: string; texto: string } | null>(null);

  const [gerador, setGerador] = useState('microtubule');
  const [nPar, setNPar] = useState(160);
  const [nPerp, setNPerp] = useState(13);
  const [seam, setSeam] = useState(3);
  const [fechado, setFechado] = useState(false);
  const [modulos, setModulos] = useState(8);
  const [jPerp, setJPerp] = useState(1);
  const [wsP, setWsP] = useState(0);
  const [religar, setReligar] = useState(true);
  const [semente, setSemente] = useState(2026);
  const [pIn, setPIn] = useState(0.3);
  const [pOut, setPOut] = useState(0.01);
  const [nGen, setNGen] = useState(400);

  const [ham, setHam] = useState('adjacency');
  const [gamma, setGamma] = useState(1);
  const [norm, setNorm] = useState('spectral');
  const [t1, setT1] = useState(80);
  const [nt, setNt] = useState(400);
  const [sitio, setSitio] = useState(0);
  const [alvo, setAlvo] = useState(-1);
  const [modoClique, setModoClique] = useState<'inicial' | 'alvo'>('inicial');

  const [vPassos, setVPassos] = useState(9);
  const [vReal, setVReal] = useState(4);
  const [pontos, setPontos] = useState<Ponto[]>([]);
  const [varrendo, setVarrendo] = useState<{ i: number; n: number } | null>(null);

  const [rede, setRede] = useState<Rede | null>(null);
  const [versao, setVersao] = useState('');
  const [series, setSeries] = useState<Float64Array | null>(null);
  const [cursor, setCursor] = useState(0);
  const [total, setTotal] = useState(0);
  const [rodando, setRodando] = useState(false);
  const [tocando, setTocando] = useState(false);
  const [nQuadros, setNQuadros] = useState(1);
  const [passoQuadro, setPassoQuadro] = useState(1);
  const [diagQuadros, setDiagQuadros] = useState<{ tem: number; esperado: number } | null>(null);
  const [erro, setErro] = useState<string | null>(null);
  const [fase, setFase] = useState('f_iniciando');
  const [autoPropagar, setAutoPropagar] = useState(false);
  const [escalaFixa, setEscalaFixa] = useState(false);
  const [importado, setImportado] = useState<Reimportado | null>(null);
  const [maxQuadro, setMaxQuadro] = useState(0);
  /* 'desenrolada' só existe quando o gerador deu a rede desenrolada; 'propria'
     quando ele deu geometria de verdade (linha, ciclo, grade); 'espectral' vale
     para qualquer grafo, e é o padrão de quem não tem geometria nenhuma. */
  const [vista, setVista] = useState<'desenrolada' | 'propria' | 'espectral' | 'tubo'>('desenrolada');
  /* seam_shift vem do spec CANÔNICO que o núcleo devolveu, não do controle na
     tela: é o número que construiu a rede desenhada. Lendo o controle, mexer
     nele sem regerar deixaria o cilindro descrevendo um grafo que não está ali,
     e a discordância não apareceria como erro — apareceria como uma hélice
     ligeiramente errada, que ninguém confere de olho. */
  const [seamShift, setSeamShift] = useState(0);
  /* Falso = eixo encurtado para caber; verdadeiro = escala métrica. O fator
     aparece na tela nos dois casos: uma distorção silenciosa seria pior que
     nenhuma vista. */
  const [tuboMetrico, setTuboMetrico] = useState(false);
  /* Trajetórias. `nTraj = 0` é execução unitária e o bloco nem entra no spec.
     `modoSaida` começa NULO de propósito: os dois modos resolvem problemas
     diferentes, e escolher por omissão escolheria por quem não sabe que estava
     escolhendo. Quem valida isso é o parser em C, não esta tela. */
  const [nTraj, setNTraj] = useState(0);
  const [gammaDeph, setGammaDeph] = useState(0.1);
  const [rhoStride, setRhoStride] = useState(10);
  const [modoSaida, setModoSaida] = useState<ModoSaida | null>(null);
  const [xyEspectral, setXyEspectral] = useState<Float32Array | null>(null);
  const [xyProprio, setXyProprio] = useState<Float32Array | null>(null);

  /* O spec.json É a interface: a mesma moeda que o WASM lê, que o CSV carrega e
     que o .cpp exportado embute. Não existe caminho paralelo de parâmetros. */
  const spec = useMemo(() => {
    const conn = { ws_p: wsP, conn_mode: religar ? 'rewire' : 'add' };
    const params: Record<string, unknown> =
      gerador === 'microtubule'
        ? { n_par: nPar, n_perp: nPerp, seam_shift: seam, longitudinal_closed: fechado,
            j_par: 1, j_perp: jPerp, n_modules: modulos, ...conn }
        : gerador === 'sbm' ? { n: nGen, n_modules: modulos, p_in: pIn, p_out: pOut, ...conn }
        : gerador === 'grid2d' ? { rows: 20, cols: 20, ...conn }
        : gerador === 'hypercube' ? { dim: 9, ...conn }
        : { n: nGen, ...conn };
    return {
      format_version: 2, seed: semente,
      graph: { generator: gerador, params },
      hamiltonian: { kind: ham, gamma, normalization: norm, lanczos_steps: 40 },
      initial: { site: sitio },
      time: { t1, nt },
      observables: { target: alvo, population: false, module_concurrence: true, pop_stride: 1 },
      realizations: 1,
      /* O bloco só existe quando há trajetórias, e `output_mode` só entra
         quando foi escolhido: um spec com o campo preenchido por omissão
         mentiria sobre ter havido escolha. Sem ele, o parser em C recusa com a
         mensagem certa — a regra mora lá, não aqui. */
      ...(nTraj > 0
        ? { trajectories: {
              n_traj: nTraj, gamma_deph: gammaDeph, rho_stride: rhoStride,
              ...(modoSaida ? { output_mode: modoSaida } : {}),
            } }
        : {}),
    };
  }, [gerador, nPar, nPerp, seam, fechado, jPerp, modulos, wsP, religar, semente,
      pIn, pOut, nGen, ham, gamma, norm, t1, nt, sitio, alvo,
      nTraj, gammaDeph, rhoStride, modoSaida]);

  /* O custo dos dois modos, dos parâmetros ATUAIS. Ele muda quando o usuário
     mexe em N, no número de trajetórias ou na amostragem — que é o ponto: um
     número que reage é lido, uma nota de rodapé não. */
  const custoTraj = useMemo(() => {
    const n = rede?.n ?? 0;
    const amostras = amostrasDeTempo(nt, rhoStride);
    return {
      amostras,
      rho: formatarBytes(bytesAcumularRho(n, amostras)),
      psi: formatarBytes(bytesArquivarPsi(n, amostras, nTraj)),
      psiPlano: formatarBytes(bytesArquivarPsi(n, amostras, nTraj) * 100),
    };
  }, [rede, nt, rhoStride, nTraj]);

  /* O passo pedido, elevado ao mínimo que cabe no teto de memória. O custo
     aparece ao lado do controle: guardar um quadro por passo em rede grande
     estoura, e é melhor o usuário ver o número do que descobrir travando. */
  const passoEfetivo = useMemo(() => {
    if (!rede) return Math.max(1, passoQuadro);
    return Math.max(passoQuadro, passoParaCaber(nt, rede.n, TETO_QUADROS));
  }, [passoQuadro, rede, nt]);
  const custoQuadros = useMemo(() => {
    if (!rede) return { k: 0, mb: '0' };
    const k = passosDosQuadros(nt, passoEfetivo).length;
    return { k, mb: (bytesDosQuadros(nt, passoEfetivo, rede.n) / 1e6).toFixed(1) };
  }, [rede, nt, passoEfetivo]);

  const forma = useMemo(() => ({
    largura: rede?.nPar || 1, altura: rede?.nPerp || 1,
  }), [rede]);
  const xyVista = vista === 'espectral' ? xyEspectral : xyProprio;
  /* O cilindro é derivado de (m, q), como o núcleo previu: o grafo tem
     topologia, não geometria tridimensional, e o tubo é uma leitura dela. */
  const compressaoTubo = useMemo(
    () => (rede && !tuboMetrico ? compressaoParaCaber(rede.nPar, rede.nPerp) : 1),
    [rede, tuboMetrico]);
  const xyzTubo = useMemo(
    () => (rede && rede.geom === 2
      ? coordenadasTubo(rede.nPar, rede.nPerp, seamShift, compressaoTubo) : null),
    [rede, seamShift, compressaoTubo]);

  const desenhar = useCallback((pop: Float32Array) => {
    if (!rede) return;
    let max = 0;
    for (let i = 0; i < pop.length; ++i) if (pop[i] > max) max = pop[i];
    const usado = escalaFixa ? Math.min(1, 8 / rede.n) : max;
    setMaxQuadro(usado);
    if (vista === 'desenrolada') {
      const h = heatRef.current;
      if (!h) return;
      h.desenhar(empacotarLattice(pop, forma.largura, forma.altura),
                 forma.largura, forma.altura, usado);
    } else if (vista === 'tubo') {
      tuboRef.current?.desenhar(pop, usado);
    } else {
      const nu = nuvemRef.current;
      if (!nu || !xyVista) return;
      nu.desenhar(pop, usado);
    }
  }, [forma, rede, escalaFixa, vista, xyVista]);
  desenharRef.current = desenhar;


  const worker = useCallback(() => {
    if (!workerRef.current) {
      const w = new Worker(new URL('../nucleo/daedalus.worker.ts', import.meta.url),
                           { type: 'module' });
      /* Erro na CARGA do módulo acontece antes de qualquer onmessage e não vira
         mensagem: sem isto o worker morre e a interface espera para sempre. */
      w.onerror = (ev: ErrorEvent) => {
        setErro(`worker: ${ev.message || '?'}`); setRodando(false); setVarrendo(null);
      };
      w.onmessage = (ev) => {
        const m = ev.data;
        if (m.tipo === 'erro') { setErro(m.mensagem); setRodando(false); setVarrendo(null); setFase('f_erro'); return; }
        if (m.tipo === 'rede') {
          redeRef.current = m.rede;
          setRede(m.rede); setVersao(`${m.versao} · ${m.hashNucleo}`);
          setSeries(null); setCursor(0); setTotal(m.rede.nt);
          quadrosRef.current = [{ pop: m.pop, passo: 0 }];
          setNQuadros(1); setDiagQuadros(null);
          arestasRef.current = m.arestas ?? null;
          setXyProprio(m.xy);
          setXyEspectral(null);
          /* Do spec canônico: é a conta que o NÚCLEO fez, depois dos padrões. */
          try {
            const c = JSON.parse(m.canonico);
            setSeamShift(Number(c?.graph?.params?.seam_shift ?? 0));
          } catch { setSeamShift(0); }
          /* rede desenrolada é o padrão quando existe; layout espectral quando
             o gerador não sabe dar geometria nenhuma. */
          setVista(m.rede.geom === 2 ? 'desenrolada'
                 : m.rede.geom === 1 ? 'propria' : 'espectral');
          setErro(null); setFase('f_pronta');
          if (primeira.current) { primeira.current = false; setAutoPropagar(true); }
        } else if (m.tipo === 'reiniciado') {
          /* Os quadros da corrida anterior nao se somam a esta. */
          quadrosRef.current = [{ pop: m.pop, passo: 0 }];
          setNQuadros(1); setCursor(0); setDiagQuadros(null);
        } else if (m.tipo === 'quadro') {
          /* Cada quadro guarda o PASSO em que foi tirado: e ele que rotula o
             tempo e posiciona o cursor das series. O deslizador conta QUADROS. */
          quadrosRef.current.push({ pop: m.pop, passo: m.passo });
          setNQuadros(quadrosRef.current.length);
          setCursor(quadrosRef.current.length - 1);
          desenharRef.current(m.pop);
        } else if (m.tipo === 'pronto' || m.tipo === 'cancelado') {
          setSeries(m.series); setRodando(false);
          setFase(m.tipo === 'pronto' ? 'f_pronto' : 'f_cancelado');
          const ultimo = quadrosRef.current[quadrosRef.current.length - 1];
          if (ultimo) desenharRef.current(ultimo.pop);
          /* Instrumentacao: quantos quadros o worker disse que emitiria, quantos
             a interface tem. Se divergirem, a animacao para antes do fim e o
             usuario precisa saber POR QUE, nao so ver o deslizador travado. */
          if (m.tipo === 'pronto') {
            const tem = quadrosRef.current.length - 1;
            setDiagQuadros(tem === m.previstos ? null : { tem, esperado: m.previstos });
          }
        } else if (m.tipo === 'exportado') {
          /* Espacados: o navegador trata um lote simultaneo de downloads como
             suspeito e engole os seguintes em silencio. O Wolfram sai em tres
             arquivos, e dois deles sumirem sem aviso seria o pior modo de
             falhar — o usuario ficaria com um pacote que nao roda. */
          (m.arquivos as Array<{ nome: string; texto: string }>).forEach((a, i) =>
            setTimeout(() => baixar(a.nome, a.texto), i * 250));
        } else if (m.tipo === 'csv') {
          baixar('daedalus.csv', m.texto);
        } else if (m.tipo === 'varredura_passo') {
          setPontos(m.pontos); setVarrendo({ i: m.i, n: m.total });
        } else if (m.tipo === 'varredura_fim') {
          setPontos(m.pontos); setVarrendo(null); setFase('f_pronto');
        } else if (m.tipo === 'espectral') {
          setXyEspectral(m.xy);
        } else if (m.tipo === 'validado') {
          /* O spec do CSV passou (ou não) pelo parser estrito em C. */
          const pend = pendenteCsv.current;
          pendenteCsv.current = null;
          if (!pend) return;
          const r = reimportar(pend.nome, pend.texto, m.hash, () => m.erro ?? null);
          setImportado(r);
        }
      };
      workerRef.current = w;
    }
    return workerRef.current;
  }, []);

  /* Um canvas por vez, MONTADO, não escondido com display:none.
     Manter os dois no DOM e alternar a visibilidade parece mais barato, e
     custa caro: um canvas que nunca foi disposto na página não tem tamanho,
     o contexto WebGL nasce com um buffer de outra dimensão, e o que aparece ao
     trocar de vista é conteúdo antigo esticado. Montar só o ativo custa recriar
     o contexto na troca — milissegundos — e elimina a classe inteira. */
  useEffect(() => {
    try {
      if (vista === 'desenrolada') {
        nuvemRef.current = null; tuboRef.current = null;
        if (canvasRef.current) heatRef.current = new Heatmap(canvasRef.current);
      } else if (vista === 'tubo') {
        heatRef.current = null; nuvemRef.current = null;
        if (canvasTubo.current) {
          tuboRef.current = new Tubo(canvasTubo.current);
          tuboRef.current.orientar(orientRef.current[0], orientRef.current[1]);
        }
      } else {
        heatRef.current = null; tuboRef.current = null;
        if (canvasNuvem.current) nuvemRef.current = new Nuvem(canvasNuvem.current);
      }
    } catch (e: any) { setErro(String(e?.message ?? e)); }
  }, [vista]);

  /* REDESENHAR AO REDIMENSIONAR. O conteúdo de um canvas WebGL não se
     re-renderiza sozinho: mudou o tamanho, o que está lá vira imagem antiga
     esticada. Acontece na troca de vista — o canvas nasce com um tamanho, é
     desenhado, e o layout o redimensiona depois — e também ao mexer na
     janela. */
  useEffect(() => {
    const alvos = [canvasRef.current, canvasNuvem.current, canvasTubo.current]
      .filter(Boolean) as HTMLCanvasElement[];
    if (alvos.length === 0 || typeof ResizeObserver === 'undefined') return;
    const ro = new ResizeObserver(() => {
      const q = quadrosRef.current[Math.min(cursor, quadrosRef.current.length - 1)];
      if (q) desenharRef.current(q.pop);
    });
    for (const a of alvos) ro.observe(a);
    return () => ro.disconnect();
  }, [cursor]);

  /* O layout espectral é pedido só quando a vista é escolhida. */
  useEffect(() => {
    if (vista === 'espectral' && !xyEspectral && rede) {
      worker().postMessage({ tipo: 'espectral' });
    }
  }, [vista, xyEspectral, rede, worker]);

  useEffect(() => {
    if (vista === 'desenrolada' || vista === 'tubo' || !nuvemRef.current || !xyVista || !rede) return;
    nuvemRef.current.rede(xyVista, rede.n, arestasRef.current);
    const q = quadrosRef.current[Math.min(cursor, quadrosRef.current.length - 1)];
    if (q) desenharRef.current(q.pop);
    /* eslint-disable-next-line */
  }, [vista, xyVista, rede]);

  useEffect(() => {
    if (vista !== 'tubo' || !tuboRef.current || !xyzTubo || !rede) return;
    tuboRef.current.rede(xyzTubo, rede.n, arestasRef.current);
    const q = quadrosRef.current[Math.min(cursor, quadrosRef.current.length - 1)];
    if (q) desenharRef.current(q.pop);
    /* eslint-disable-next-line */
  }, [vista, xyzTubo, rede]);

  const carregar = () => {
    setErro(null); setFase('f_gerando'); setImportado(null);
    worker().postMessage({ tipo: 'carregar', spec: JSON.stringify(spec) });
  };
  const propagar = () => {
    setErro(null); setRodando(true); setSeries(null); setCursor(0); setFase('f_propagando');
    worker().postMessage({ tipo: 'propagar', passoQuadro: passoEfetivo });
  };
  const varrer = () => {
    setPontos([]); setVarrendo({ i: 0, n: vPassos });
    worker().postMessage({ tipo: 'varrer', spec: JSON.stringify(spec),
                           pMin: 0, pMax: 1, passos: vPassos, realizacoes: vReal });
  };

  useEffect(() => { carregar(); /* eslint-disable-next-line */ }, []);
  useEffect(() => {
    if (autoPropagar && rede) { setAutoPropagar(false); propagar(); }
    /* eslint-disable-next-line */
  }, [autoPropagar, rede]);

  useEffect(() => {
    if (!tocando || quadrosRef.current.length < 2) return;
    const id = setInterval(() => {
      setCursor((c) => {
        const prox = c + 1 > quadrosRef.current.length - 1 ? 0 : c + 1;
        const q = quadrosRef.current[prox];
        if (q) desenharRef.current(q.pop);
        return prox;
      });
    }, 40);
    return () => clearInterval(id);
  }, [tocando]);

  /** Índice de QUADRO, não de passo. */
  const irParaQuadro = (i: number) => {
    setTocando(false);
    const k = Math.max(0, Math.min(quadrosRef.current.length - 1, i));
    setCursor(k);
    const q = quadrosRef.current[k];
    if (q) desenharRef.current(q.pop);
  };
  /** Recebe um PASSO (vindo do clique nas séries) e vai ao quadro mais próximo. */
  const irParaPasso = (passo: number) => {
    irParaQuadro(quadroDoPasso(quadrosRef.current.map((q) => q.passo), passo));
  };

  const marcar = (j: number) => {
    if (j < 0 || !rede || j >= rede.n) return;
    if (modoClique === 'inicial') setSitio(j); else setAlvo(j);
  };
  const clicarLattice = (ev: React.MouseEvent<HTMLCanvasElement>) => {
    const h = heatRef.current;
    if (!h || !rede) return;
    const p = h.texelDoEvento(ev);
    if (p) marcar(sitioNoTexel(p.m, p.q, forma.largura, forma.altura));
  };
  const clicarNuvem = (ev: React.MouseEvent<HTMLCanvasElement>) => {
    if (nuvemRef.current) marcar(nuvemRef.current.sitioNoEvento(ev));
  };

  /* GIRAR E CLICAR NO MESMO GESTO. Arrastar gira o tubo; soltar sem ter
     arrastado marca o vértice. Sem o limiar de deslocamento, todo giro
     terminaria escolhendo um sítio novo — o usuário perderia a origem da
     caminhada só por ter olhado o outro lado do cilindro. */
  const tuboDesce = (ev: React.PointerEvent<HTMLCanvasElement>) => {
    arrasto.current = { x: ev.clientX, y: ev.clientY, moveu: false };
    /* A captura é conveniência — o giro continua se o ponteiro sair do canvas —
       e não é garantida: um ponteiro sintético não existe para o navegador e
       setPointerCapture lança. Sem o try, o giro morreria na primeira linha do
       gesto e a vista pareceria travada. */
    try { ev.currentTarget.setPointerCapture(ev.pointerId); } catch { /* segue */ }
  };
  const tuboMove = (ev: React.PointerEvent<HTMLCanvasElement>) => {
    const a = arrasto.current, tb = tuboRef.current;
    if (!a || !tb) return;
    const dx = ev.clientX - a.x, dy = ev.clientY - a.y;
    if (!a.moveu && Math.hypot(dx, dy) < 4) return;
    a.moveu = true; a.x = ev.clientX; a.y = ev.clientY;
    /* horizontal aponta o EIXO para dentro ou para fora; vertical GIRA o tubo
       em torno do próprio eixo, trocando que protofilamento fica de frente */
    const [az, el] = orientRef.current;
    orientRef.current = [az + dy * 0.007, el + dx * 0.007];
    tb.orientar(orientRef.current[0], orientRef.current[1]);
    orientRef.current = tb.orientacao();          /* a trava do polo é dele */
    const q = quadrosRef.current[Math.min(cursor, quadrosRef.current.length - 1)];
    if (q) desenharRef.current(q.pop);
  };
  const tuboSobe = (ev: React.PointerEvent<HTMLCanvasElement>) => {
    const a = arrasto.current;
    arrasto.current = null;
    if (a && !a.moveu && tuboRef.current) marcar(tuboRef.current.sitioNoEvento(ev));
  };

  /* A reimportação é a única entrada que não nasce aqui: o spec embutido no CSV
     vai ao parser estrito em C antes de qualquer coisa ser plotada. */
  const abrirCsv = async (arquivo: File | undefined) => {
    if (!arquivo) return;
    const texto = await arquivo.text();
    const linhas = texto.split('\n');
    const spec = linhas.find((l) => l.startsWith('#! spec '))?.slice(8) ?? '';
    pendenteCsv.current = { nome: arquivo.name, texto };
    if (!spec) {
      const r = reimportar(arquivo.name, texto, '', () => null);
      pendenteCsv.current = null;
      setImportado(r);
      return;
    }
    worker().postMessage({ tipo: 'validar', spec });
  };

  /* O passo do quadro atual: e ele que posiciona o cursor nas series e que
     rotula o tempo. `cursor` conta QUADROS. */
  const passoAtual = quadrosRef.current[cursor]?.passo ?? 0;
  const passo = Math.max(0, Math.min(total - 1, passoAtual - 1));
  const norma = series && total > 0 ? series[passo * 4 + 0] : NaN;
  const desvio = Number.isFinite(norma) ? Math.abs(norma - 1) : NaN;
  const normaOk = !Number.isFinite(desvio) || desvio < 1e-9;

  /* Eixo x em tempo, nao em passo: passo e detalhe da discretizacao, e comparar
     duas corridas com grades diferentes exige a mesma abscissa fisica. */
  const tempo = useMemo(() => {
    if (importado?.utilizavel) {
      const tt = importado.dados.get('t');
      let g = 1;
      try { g = JSON.parse(importado.meta.get('spec') ?? '{}').hamiltonian?.gamma ?? 1; } catch { g = 1; }
      if (tt) { const v = new Float64Array(tt.length); for (let i = 0; i < tt.length; ++i) v[i] = g * tt[i]; return v; }
      return new Float64Array(0);
    }
    if (!rede || total === 0) return new Float64Array(0);
    const v = new Float64Array(total);
    for (let i = 0; i < total; ++i) v[i] = gamma * (i + 1) * rede.dt;
    return v;
  }, [importado, rede, total, gamma]);

  const listaSeries = useMemo(() => {
    /* Reimportado tem prioridade e vem em BRONZE: o número deixou de ser
       calculado aqui, e a cor diz isso. */
    /* Reimportado vem em BRONZE: o numero deixou de ser calculado aqui. */
    const bronze = !!importado?.utilizavel;
    const cIpr = bronze ? '#B5813C' : '#E5A83F';
    const cCoh = bronze ? '#C99A56' : '#3F7C74';
    const out: Serie[] = [];
    if (bronze) {
      const col = (n: string) => importado!.dados.get(n);
      if (col('ipr')) out.push({ nome: t('se_ipr'), cor: cIpr, valores: col('ipr')!, eixo: 'esq' });
      if (col('coh_l1')) out.push({ nome: t('se_coerencia'), cor: cCoh, valores: col('coh_l1')!, eixo: 'dir' });
      const pa = col('p_target');
      if (pa && Number.isFinite(pa[0])) out.push({ nome: t('se_palvo'), cor: '#A8452C', valores: pa, eixo: 'esq' });
      return out;
    }
    if (!series || total === 0) return out;
    const col = (k: number) => {
      const v = new Float64Array(total);
      for (let i = 0; i < total; ++i) v[i] = series[i * 4 + k];
      return v;
    };
    out.push({ nome: t('se_ipr'), cor: cIpr, valores: col(1), eixo: 'esq' });
    out.push({ nome: t('se_coerencia'), cor: cCoh, valores: col(2), eixo: 'dir' });
    const pa = col(3);
    /* p no alvo divide o eixo esquerdo com o IPR: as duas vivem em [0,1]. */
    if (Number.isFinite(pa[0])) out.push({ nome: t('se_palvo'), cor: '#A8452C', valores: pa, eixo: 'esq' });
    return out;
  }, [series, total, importado, t]);

  const exportar = (alvoExp: string) =>
    worker().postMessage({ tipo: 'exportar', alvo: alvoExp, spec: JSON.stringify(spec) });

  const metricas = useMemo(() => {
    if (!rede) return [];
    return [
      { rotulo: 'N', valor: String(rede.n) },
      { rotulo: '|E|', valor: String(rede.arestas) },
      { rotulo: '⟨d⟩', valor: rede.grauMedio.toFixed(3) },
      { rotulo: 'λ₂', valor: rede.lambda2.toPrecision(5),
        nota: rede.lambda2Convergiu ? undefined : t('r_lambda2_alto'),
        alerta: !rede.lambda2Convergiu },
      { rotulo: 'Q', valor: rede.Q.toFixed(4) },
      { rotulo: 'α', valor: rede.alpha.toFixed(3) },
    ];
  }, [rede, t]);

  const faixa = (rot: string, valor: string, campo: React.ReactNode) => (
    <>
      <div className="campo"><label>{rot}</label><span className="num">{valor}</span></div>
      {campo}
    </>
  );

  return (
    <div className="app" onDragOver={(e) => e.preventDefault()}
         onDrop={(e) => { e.preventDefault(); void abrirCsv(e.dataTransfer.files[0]); }}>
      <header className="topo">
        <svg width="22" height="22" viewBox="0 0 100 100" className="marca" aria-hidden="true">
          <path d="M19 19 V81 M19 19 H81 V81 H19 M19 34.5 H65.5 V65.5 H19 M19 50 H50"
                fill="none" stroke="currentColor" strokeWidth="9"
                strokeLinecap="square" strokeLinejoin="miter" />
        </svg>
        <span className="nome">DAEDALUS</span>
        {/* Azul Egeu: calculado aqui, agora. Bronze: veio de fora.
            E QUEM calculou, que e outra pergunta: o mesmo core_hash pode ter
            produzido o numero pelo nucleo em C ou pelo pacote Wolfram, que usa
            outro metodo. Um resultado vindo do Wolfram nao deve se passar por
            um do nucleo, e um CSV que nao declara origem tambem nao. */}
        <span className={importado?.utilizavel ? 'selo exportado' : 'selo'}>
          {importado?.utilizavel
            ? `${t('modo_reimportado')} · ${rotuloOrigem(importado.meta.get('implementacao'))}`
            : t('modo_local')}
        </span>
        <span className="espaco" />
        <span className="arquivo mono">
          {importado?.utilizavel ? importado.nome
            : `${t(`ger_${gerador}`)}${rede ? ` · ${rede.n}` : ''}`}
        </span>
        <SeletorLingua />
      </header>

      <div className="corpo">
        <main className="palco">
          {erro && <div className="erro">{erro}</div>}
          {importado?.avisos.map((a, i) => (
            <div key={i} className={a.grave ? 'erro' : 'aviso'}>{t(a.chave, a.params)}</div>
          ))}
          <div className="mapa">
            {vista === 'desenrolada'
              ? <canvas ref={canvasRef} onClick={clicarLattice} />
              : vista === 'tubo'
              ? <canvas ref={canvasTubo} className="giravel"
                        onPointerDown={tuboDesce} onPointerMove={tuboMove}
                        onPointerUp={tuboSobe} onPointerCancel={tuboSobe} />
              : <canvas ref={canvasNuvem} onClick={clicarNuvem} />}
          </div>
          {/* O nome do layout em uso fica visível: a figura muda de significado
              conforme ele, e adivinhar qual está ativo é o começo de ler errado. */}
          <div className="vistas">
            {rede?.geom === 2 && (
              <button className={vista === 'desenrolada' ? 'sel' : ''}
                      onClick={() => setVista('desenrolada')}>{t('v_desenrolada')}</button>)}
            {rede?.geom === 2 && (
              <button className={vista === 'tubo' ? 'sel' : ''}
                      onClick={() => setVista('tubo')}>{t('v_tubo')}</button>)}
            {rede && rede.geom === 1 && (
              <button className={vista === 'propria' ? 'sel' : ''}
                      onClick={() => setVista('propria')}>{t('v_propria')}</button>)}
            <button className={vista === 'espectral' ? 'sel' : ''}
                    onClick={() => setVista('espectral')}>{t('v_espectral')}</button>
            <span className="espaco" />
            {vista === 'tubo' && (
              <label className="caixa"><input type="checkbox" checked={tuboMetrico}
                     onChange={(e) => setTuboMetrico(e.target.checked)} /> {t('v_tubo_metrico')}</label>)}
            {vista === 'espectral' && rede && rede.n > MAX_ARESTAS_DESENHADAS &&
              <span className="dica">{t('v_sem_arestas')}</span>}
            {vista === 'espectral' && !xyEspectral && rede && rede.n > 20000 &&
              <span className="aviso">{t('v_sem_espectral')}</span>}
          </div>
          <div className="escala">
            <span>0.00</span>
            <div className="rampa" style={{ background: gradienteCss() }} />
            <span>{maxQuadro ? maxQuadro.toPrecision(3) : '—'}</span>
            <span>|ψⱼ|²</span>
          </div>
          {vista === 'espectral' && <p className="dica">{t('v_espectral_dica')}</p>}
          {vista === 'tubo' && <p className="dica">{t('v_tubo_dica')}</p>}
          {vista === 'tubo' && (
            <p className={compressaoTubo < 1 ? 'aviso' : 'dica'}>
              {t('v_tubo_escala', { fator: compressaoTubo.toFixed(2) })}
            </p>)}
          {vista === 'tubo' && fechado && <p className="dica">{t('v_tubo_fechado')}</p>}
          <Metricas itens={metricas} />
          <div className="transporte">
            <button onClick={() => setTocando((v) => !v)}
                    disabled={quadrosRef.current.length < 2}>{tocando ? '❚❚' : '▶'}</button>
            <input type="range" min={0} max={Math.max(1, nQuadros - 1)} value={cursor}
                   onChange={(e) => irParaQuadro(+e.target.value)} />
            {/* Rotulo em TEMPO, nao em indice de quadro. */}
            <span className="num">γt = {(gamma * passoAtual * (rede?.dt ?? 0)).toFixed(2)}</span>
            <span className="num fraco">{cursor + 1}/{nQuadros}</span>
            {diagQuadros && <span className="aviso">{t('q_incompleto', diagQuadros)}</span>}
            <label className="caixa"><input type="checkbox" checked={escalaFixa}
                   onChange={(e) => setEscalaFixa(e.target.checked)} /> {t('a_escala_fixa')}</label>
          </div>
          {listaSeries.length > 0 && (
            <Series series={listaSeries} tempo={tempo} cursor={cursor} rotuloX="γt"
                    aoClicar={importado?.utilizavel ? undefined : irParaPasso} />
          )}
          {pontos.length > 0 && (
            <Varredura pontos={pontos} rotuloX={t('v_eixo')} rotuloY={t('v_media')} />
          )}
        </main>

        <aside className="painel">
          <Secao titulo={t('sec_gerador')} aberta
                 resumo={`${t(`ger_${gerador}`)}${gerador === 'microtubule' ? ` ${nPar}×${nPerp}` : ` N=${nGen}`} · p=${wsP.toFixed(2)}`}>
            <select value={gerador} onChange={(e) => setGerador(e.target.value)}>
              {GERADORES.map((g) => <option key={g} value={g}>{t(`ger_${g}`)}</option>)}
            </select>

            {gerador === 'microtubule' && (<>
              {faixa(t('c_npar'), String(nPar),
                <input type="range" min={4} max={600} value={nPar} onChange={(e) => setNPar(+e.target.value)} />)}
              {faixa(t('c_nperp'), String(nPerp),
                <input type="range" min={3} max={26} value={nPerp} onChange={(e) => setNPerp(+e.target.value)} />)}
              {faixa(t('c_costura'), String(seam),
                <input type="range" min={0} max={12} value={seam} onChange={(e) => setSeam(+e.target.value)} />)}
              <label className="caixa"><input type="checkbox" checked={fechado}
                onChange={(e) => setFechado(e.target.checked)} /> {t('c_fechar')}</label>
              {faixa(t('c_acoplamento'), `1.00 / ${jPerp.toFixed(2)}`,
                <input type="range" min={0.1} max={2} step={0.05} value={jPerp}
                       onChange={(e) => setJPerp(+e.target.value)} />)}
            </>)}
            {gerador === 'sbm' && (<>
              {faixa(t('c_n'), String(nGen),
                <input type="range" min={20} max={2000} step={10} value={nGen}
                       onChange={(e) => setNGen(+e.target.value)} />)}
              {faixa(t('c_pinpout'), `${pIn.toFixed(2)} / ${pOut.toFixed(3)}`, <>
                <input type="range" min={0} max={1} step={0.01} value={pIn} onChange={(e) => setPIn(+e.target.value)} />
                <input type="range" min={0} max={0.2} step={0.002} value={pOut} onChange={(e) => setPOut(+e.target.value)} />
              </>)}
            </>)}
            {(gerador === 'path' || gerador === 'cycle' || gerador === 'complete') &&
              faixa(t('c_n'), String(nGen),
                <input type="range" min={4} max={2000} step={2} value={nGen}
                       onChange={(e) => setNGen(+e.target.value)} />)}
            {(gerador === 'microtubule' || gerador === 'sbm') &&
              faixa(t('c_modulos'), String(modulos),
                <input type="range" min={1} max={24} value={modulos}
                       onChange={(e) => setModulos(+e.target.value)} />)}

            {faixa(t('c_religacao'), wsP.toFixed(2),
              <input type="range" min={0} max={1} step={0.01} value={wsP}
                     onChange={(e) => setWsP(+e.target.value)} />)}
            <label className="caixa"><input type="checkbox" checked={religar}
              onChange={(e) => setReligar(e.target.checked)} /> {t('c_religar_fixo')}</label>
            {!religar && <p className="aviso">{t('av_acrescentar')}</p>}
            {faixa(t('c_semente'), String(semente),
              <input type="number" value={semente} onChange={(e) => setSemente(+e.target.value)} />)}
            <button className="primario larga" onClick={carregar}>{t('a_gerar')}</button>
          </Secao>

          <Secao titulo={t('sec_hamiltoniano')}
                 resumo={`${ham === 'adjacency' ? '−γA' : 'γL'} · γ=${gamma.toFixed(2)} · ${t(`norm_${norm}`)}`}>
            <div className="segmentado">
              <button className={ham === 'adjacency' ? 'sel' : ''}
                      onClick={() => setHam('adjacency')}>{t('ham_adjacency')}</button>
              <button className={ham === 'laplacian' ? 'sel' : ''}
                      onClick={() => setHam('laplacian')}>{t('ham_laplacian')}</button>
            </div>
            {faixa(t('c_gamma'), gamma.toFixed(2),
              <input type="range" min={0.05} max={4} step={0.05} value={gamma}
                     onChange={(e) => setGamma(+e.target.value)} />)}
            <div className="campo"><label>{t('c_normalizacao')}</label></div>
            <select value={norm} onChange={(e) => setNorm(e.target.value)}>
              <option value="spectral">{t('norm_spectral')}</option>
              <option value="mean_degree">{t('norm_mean_degree')}</option>
              <option value="none">{t('norm_none')}</option>
            </select>
            {norm === 'none' && <p className="aviso">{t('av_sem_normalizar')}</p>}
          </Secao>

          <Secao titulo={t('sec_tempo')} aberta
                 resumo={`t=${t1.toFixed(0)} · ${nt} ${t('c_pontos')}`}>
            {faixa(t('c_tfinal'), t1.toFixed(1),
              <input type="range" min={1} max={400} step={1} value={t1} onChange={(e) => setT1(+e.target.value)} />)}
            {faixa(t('c_pontos'), String(nt),
              <input type="range" min={20} max={2000} step={20} value={nt}
                     onChange={(e) => setNt(+e.target.value)} />)}
            {faixa(t('c_passo_quadro'), t('q_quadros', custoQuadros),
              <input type="range" min={1} max={50} value={passoQuadro}
                     onChange={(e) => setPassoQuadro(+e.target.value)} />)}
            {passoEfetivo > passoQuadro &&
              <p className="dica">{t('q_quadros', custoQuadros)}</p>}
            <div className="transporte">
              <button className="primario" style={{ flex: 1 }} onClick={propagar}
                      disabled={rodando || !rede}>{rodando ? t('a_propagando') : t('a_propagar')}</button>
              <button onClick={() => worker().postMessage({ tipo: 'cancelar' })}
                      disabled={!rodando && !varrendo}>{t('a_parar')}</button>
            </div>
          </Secao>

          <Secao titulo={t('sec_sitios')} resumo={`${sitio} → ${alvo < 0 ? '—' : alvo}`}>
            <div className="segmentado">
              <button className={modoClique === 'inicial' ? 'sel' : ''}
                      onClick={() => setModoClique('inicial')}>{t('s_inicial', { j: sitio })}</button>
              <button className={modoClique === 'alvo' ? 'sel' : ''}
                      onClick={() => setModoClique('alvo')}>{t('s_alvo', { j: alvo < 0 ? '—' : alvo })}</button>
            </div>
            <p className="dica">{t('s_dica')}{alvo >= 0 && <> <a href="#" style={{ color: '#8FC2E4' }}
              onClick={(e) => { e.preventDefault(); setAlvo(-1); }}>{t('a_limpar_alvo')}</a></>}</p>
          </Secao>

          <Secao titulo={t('sec_varredura')} resumo={`${vPassos} × ${vReal}`}>
            {faixa(t('c_passos'), String(vPassos),
              <input type="range" min={2} max={21} value={vPassos} onChange={(e) => setVPassos(+e.target.value)} />)}
            {faixa(t('c_realizacoes'), String(vReal),
              <input type="range" min={1} max={16} value={vReal} onChange={(e) => setVReal(+e.target.value)} />)}
            <button className="primario larga" onClick={varrer} disabled={!!varrendo || alvo < 0}>
              {varrendo ? t('f_varrendo', { i: varrendo.i, n: varrendo.n }) : t('a_varrer')}</button>
            <p className="dica">{t('v_dica')}</p>
          </Secao>

          {/* Bronze: tudo que sai daqui e volta de fora. Distinção funcional. */}
          <Secao titulo={t('sec_trajetorias')}
                 resumo={nTraj > 0
                   ? `${nTraj} · ${modoSaida ? t(`tr_${modoSaida === 'accumulate_rho' ? 'rho' : 'psi'}`) : t('tr_sem_modo_curto')}`
                   : t('tr_desligado')}>
            {faixa(t('tr_n'), String(nTraj),
              <input type="range" min={0} max={500} step={10} value={nTraj}
                     onChange={(e) => setNTraj(+e.target.value)} />)}
            {nTraj > 0 && <>
              {faixa(t('tr_gamma'), gammaDeph.toFixed(2),
                <input type="range" min={0} max={2} step={0.01} value={gammaDeph}
                       onChange={(e) => setGammaDeph(+e.target.value)} />)}
              {faixa(t('tr_stride'), `${custoTraj.amostras} ${t('tr_amostras')}`,
                <input type="range" min={1} max={50} value={rhoStride}
                       onChange={(e) => setRhoStride(+e.target.value)} />)}
              {/* RADIO, nao caixa: os modos sao mutuamente exclusivos, e nenhum
                  e padrao. O custo fica ao lado de cada um, calculado dos
                  parametros de agora. */}
              <div className="modos">
                <label className={modoSaida === 'accumulate_rho' ? 'modo sel' : 'modo'}>
                  <input type="radio" name="modo_saida" checked={modoSaida === 'accumulate_rho'}
                         onChange={() => setModoSaida('accumulate_rho')} />
                  <span className="modo_nome">{t('tr_rho')}</span>
                  <span className="modo_custo mono">{custoTraj.rho} {t('tr_por_thread')}</span>
                  <span className="modo_quando">{t('tr_rho_quando')}</span>
                </label>
                <label className={modoSaida === 'archive_psi' ? 'modo sel' : 'modo'}>
                  <input type="radio" name="modo_saida" checked={modoSaida === 'archive_psi'}
                         onChange={() => setModoSaida('archive_psi')} />
                  <span className="modo_nome">{t('tr_psi')}</span>
                  <span className="modo_custo mono">{custoTraj.psi} {t('tr_por_celula')}
                    {' · '}{custoTraj.psiPlano} {t('tr_plano')}</span>
                  <span className="modo_quando">{t('tr_psi_quando')}</span>
                </label>
              </div>
              {!modoSaida && <p className="aviso">{t('tr_sem_modo')}</p>}
            </>}
          </Secao>

          <Secao titulo={t('sec_exportar')} resumo="C++ · Wolfram · Python · CSV">
            <div className="exportar">
              <button onClick={() => exportar('cpp')} disabled={!rede}>C++</button>
              <button onClick={() => exportar('wl')} disabled={!rede}>Wolfram</button>
              <button onClick={() => exportar('py')} disabled={!rede}>Python</button>
            </div>
            <div className="exportar">
              <button onClick={() => worker().postMessage({ tipo: 'csv' })}
                      disabled={!series}>{t('a_csv')}</button>
            </div>
            <p className="dica">{t('av_exportar')}</p>
            <div className="exportar">
              <label className="botao-arquivo">
                {t('a_abrir_csv')}
                <input type="file" accept=".csv,text/csv,text/plain" style={{ display: 'none' }}
                       onChange={(e) => void abrirCsv(e.target.files?.[0])} />
              </label>
            </div>
          </Secao>
        </aside>
      </div>

      <footer className={normaOk ? 'rodape' : 'rodape ruim'}>
        <span>{t('r_norma')} = <b>{Number.isFinite(norma) ? norma.toFixed(12) : '—'}</b>
          {Number.isFinite(desvio) && ` (Δ ${desvio.toExponential(1)})`}</span>
        {rede && <>
          <span>hash <b>{rede.fingerprint.toString(16).slice(0, 8)}</b></span>
          {rede.componentes > 1 && <span className="alerta">{t('r_componentes', { k: rede.componentes })}</span>}
          {rede.arestasDescartadas > 0 && <span className="alerta">{t('r_duplicadas', { k: rede.arestasDescartadas })}</span>}
          {rede.religacoesFalhas > 0 && <span className="alerta">{t('r_religacoes', { k: rede.religacoesFalhas })}</span>}
        </>}
        <span className="espaco" />
        <span>{t(fase)}</span>
        <span>{versao}</span>
      </footer>
    </div>
  );
}
