/* daedalus.worker.ts — o núcleo roda aqui; a interface nunca bloqueia.
 *
 * A entrada é o TEXTO do spec.json — o mesmo que o .cpp exportado embute. A
 * propagação avança em BLOCOS e o worker volta ao laço de eventos entre eles:
 * é assim que "cancelar" funciona sem callback atravessando a fronteira do
 * WASM. O que sai daqui é sempre CÓPIA (`slice()`): uma view do heap do WASM
 * não sobrevive a um postMessage nem ao crescimento do heap.
 */
import { emitirCpp, emitirPython, emitirWolfram } from '../export/emissor.ts';
import { passosDosQuadros } from './quadros.ts';
import { Daedalus } from '../../wasm/build/daedalus.mjs';

let dae: any = null;
let cancelado = false;

async function garante() {
  if (!dae) dae = await Daedalus.criar();
  return dae;
}

self.onmessage = async (ev: MessageEvent) => {
  const msg = ev.data;
  try {
    switch (msg.tipo) {
      case 'carregar': {
        const d = await garante();
        const rede = d.carregar(msg.spec);
        /* Arestas só quando cabem na tela. O teto é mais alto que os 2000 da
           nuvem espectral, onde acima disso elas viram mancha. No tubo é o
           contrário: a estrutura é regular, e são as arestas que fazem o
           cilindro parecer cilindro — sem elas sobra uma poeira de pontos e o
           3D não se anuncia. O microtúbulo padrão tem 2080 vértices. Cada
           renderizador tem o seu próprio limite de DESENHO. */
        const arestas = rede.n <= 8000 ? d.arestas() : null;
        self.postMessage({
          tipo: 'rede', rede,
          xy: d.posicoes().slice(), modulos: d.modulos().slice(),
          arestas,
          pop: d.populacao().slice(), canonico: d.specCanonico(),
          versao: d.versao, hashNucleo: d.hashNucleo,
        });
        break;
      }
      case 'propagar': {
        const d = await garante();
        cancelado = false;
        /* Propagar de novo RECOMECA. Sem isto a segunda corrida encontraria o
           cursor no fim, nao avancaria nada, e emitiria quadros do mesmo estado
           final — a animacao "andaria" sobre uma figura parada. */
        d.reiniciar();
        self.postMessage({ tipo: 'reiniciado', pop: d.populacao().slice() });
        const nt = d.rede().nt;
        /* O passo de amostragem do MAPA vem de quem chama, e a interface o
           conhece. Se fosse escolhido aqui dentro, a animação pararia antes do
           fim — não por ter travado, mas porque os quadros acabariam ali
           enquanto o deslizador contaria passos. */
        const passo = Math.max(1, Math.floor(msg.passoQuadro ?? 1));
        const plano = passosDosQuadros(nt, passo);
        let emitidos = 0;
        let respiro = performance.now();
        for (let k = 1; k < plano.length && !cancelado; ++k) {
          d.avancar(plano[k] - plano[k - 1]);
          self.postMessage({ tipo: 'quadro', passo: d.cursor(), nt,
                             quadro: k, quadros: plano.length,
                             pop: d.populacao().slice() });
          ++emitidos;
          /* Respira por TEMPO, não por quadro: com passo 1 e nt grande, um
             setTimeout por quadro faria a propagação levar segundos só de
             ida-e-volta ao laço de eventos. */
          if (performance.now() - respiro > 16) {
            await new Promise((r) => setTimeout(r, 0));
            respiro = performance.now();
          }
        }
        self.postMessage({
          tipo: cancelado ? 'cancelado' : 'pronto',
          series: d.series().slice(), seriesModulo: d.seriesModulo().slice(),
          /* diagnóstico: quantos quadros o worker entregou, e até onde chegou */
          emitidos, previstos: plano.length - 1, cursorFinal: d.cursor(), nt,
        });
        break;
      }
      case 'cancelar':
        cancelado = true;
        break;
      case 'exportar': {
        const d = await garante();
        const canonico = d.specCanonico();
        const rede = d.rede();
        const quando = new Date().toISOString();
        const grafo = {
          n: rede.n as number,
          arestas: d.arestas() as Array<[number, number, number]>,
          modulos: Array.from(d.modulos()) as number[],
          nmod: rede.nmod as number,
          escala: rede.escala as number,
          fingerprint: String(rede.fingerprint),
        };
        const spec = JSON.parse(msg.spec);
        /* Wolfram sai como PACOTE (pacote + caderno + spec); os outros como um
           arquivo só. A lista uniformiza os dois casos para a interface. */
        const arquivos = msg.alvo === 'cpp'
          ? [{ nome: 'daedalus_run.cpp', texto: emitirCpp(canonico, quando, grafo, spec) }]
          : msg.alvo === 'wl'
          ? emitirWolfram(spec, grafo, canonico, quando)
          : [{ nome: 'daedalus_oraculo.py', texto: emitirPython(spec, grafo, canonico, quando) }];
        self.postMessage({ tipo: 'exportado', alvo: msg.alvo, arquivos });
        break;
      }
      /* O layout espectral custa um Lanczos deflacionado a mais, então só é
         calculado quando a vista é de fato pedida. */
      case 'espectral': {
        const d = await garante();
        const xy = d.espectral();
        self.postMessage({ tipo: 'espectral', xy: xy ? xy.slice() : null });
        break;
      }
      case 'validar': {
        const d = await garante();
        self.postMessage({ tipo: 'validado', erro: d.validar(msg.spec), hash: d.hashNucleo });
        break;
      }
      /* Varredura pequena: roda aqui, um ponto por vez, devolvendo progresso
         entre eles para a interface não travar. Varredura longa é o que o .cpp
         exportado com realizações existe para fazer. */
      case 'varrer': {
        const d = await garante();
        cancelado = false;
        const { pMin, pMax, passos, realizacoes } = msg;
        const pontos: Array<{ p: number; media: number; desvio: number; n: number }> = [];
        const base = JSON.parse(msg.spec);
        for (let i = 0; i < passos && !cancelado; ++i) {
          const p = passos === 1 ? pMin : pMin + ((pMax - pMin) * i) / (passos - 1);
          const amostras: number[] = [];
          for (let r = 0; r < realizacoes && !cancelado; ++r) {
            const spec = JSON.parse(JSON.stringify(base));
            spec.graph.params.ws_p = p;
            spec.seed = (base.seed ?? 0) + r;
            const rede = d.carregar(JSON.stringify(spec));
            d.avancar(rede.nt);
            const s = d.series();
            /* média TEMPORAL de p_alvo — p̄, não o valor final, que oscila */
            let soma = 0, cont = 0;
            for (let k = 0; k < rede.nt; ++k) {
              const v = s[k * 4 + 3];
              if (Number.isFinite(v)) { soma += v; ++cont; }
            }
            amostras.push(cont > 0 ? soma / cont : Number.NaN);
            await new Promise((res) => setTimeout(res, 0));
          }
          const bons = amostras.filter(Number.isFinite);
          const media = bons.reduce((a, b) => a + b, 0) / Math.max(1, bons.length);
          const desvio = bons.length > 1
            ? Math.sqrt(bons.reduce((a, b) => a + (b - media) ** 2, 0) / (bons.length - 1))
            : 0;
          pontos.push({ p, media, desvio, n: bons.length });
          self.postMessage({ tipo: 'varredura_passo', i: i + 1, total: passos, pontos });
        }
        self.postMessage({ tipo: 'varredura_fim', pontos, cancelado });
        break;
      }
      case 'csv': {
        const d = await garante();
        self.postMessage({ tipo: 'csv', texto: d.csv(true) });
        break;
      }
      default:
        self.postMessage({ tipo: 'erro', mensagem: `tipo desconhecido: ${msg.tipo}` });
    }
  } catch (e: any) {
    self.postMessage({ tipo: 'erro', mensagem: String(e?.message ?? e) });
  }
};
