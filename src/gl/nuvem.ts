/* nuvem.ts — vértices como pontos, para grafos sem geometria própria.
 *
 * O layout vem do núcleo (autovetores 2 e 3 da laplaciana). Aqui só se desenha,
 * e a única regra é NORMALIZAR PARA A CAIXA: as coordenadas cruas do embedding
 * têm escala arbitrária, e sem isso a nuvem ocuparia uma fração do canvas.
 *
 * A escala é UNIFORME nos dois eixos, não uma por eixo. Esticar x e y de forma
 * independente encheria a caixa por completo, mas faria a distância no plano
 * mentir — e o ponto do embedding espectral é justamente que distância ali
 * significa alguma coisa. Com escala uniforme, o eixo limitante encosta na
 * margem e o outro fica centrado.
 *
 * A rede desenrolada (src/gl/heatmap.ts) faz o oposto e estica por eixo, porque
 * lá o plano é uma grade de ÍNDICES e distância no desenho não mede nada. A
 * diferença entre os dois é física, não estética.
 *
 * O que se pode fazer sem mentir é GIRAR. Rotação é isometria: não muda
 * distância nenhuma. Um embedding quase unidimensional — o caso de um tubo
 * longo, onde os dois autovetores mais baixos quase se alinham — sai na
 * diagonal e ocupa uma fração da caixa; alinhar o eixo principal da nuvem com o
 * eixo maior do canvas usa a área toda e preserva a métrica.
 *
 * O raio do vértice também codifica |ψ_j|², redundância deliberada para leitura
 * sem cor.
 */
import { lut } from '../nucleo/paleta';

const VS = `#version 300 es
in vec2 aPos;
in float aVal;
uniform vec2 uEscala;
uniform vec2 uCentro;
uniform float uMax;
uniform float uRaio;
out float vT;
void main() {
  vec2 p = (aPos - uCentro) * uEscala;
  gl_Position = vec4(p.x, -p.y, 0.0, 1.0);
  vT = uMax > 0.0 ? clamp(aVal / uMax, 0.0, 1.0) : 0.0;
  /* raio cresce com a raiz da probabilidade: a ÁREA fica proporcional a ela */
  gl_PointSize = uRaio * (0.55 + 1.35 * sqrt(vT));
}`;

const FS = `#version 300 es
precision highp float;
uniform sampler2D uLut;
in float vT;
out vec4 cor;
void main() {
  vec2 d = gl_PointCoord - vec2(0.5);
  float r = length(d);
  if (r > 0.5) discard;                       /* ponto redondo, não quadrado */
  float borda = smoothstep(0.5, 0.42, r);
  cor = vec4(texture(uLut, vec2(vT, 0.5)).rgb, borda);
}`;

const VS_ARESTA = `#version 300 es
in vec2 aPos;
uniform vec2 uEscala;
uniform vec2 uCentro;
void main() {
  vec2 p = (aPos - uCentro) * uEscala;
  gl_Position = vec4(p.x, -p.y, 0.0, 1.0);
}`;

const FS_ARESTA = `#version 300 es
precision highp float;
out vec4 cor;
/* Alfa muito baixo de propósito: com ~10 mil arestas sobrepostas, 0.10 satura
   em branco e a estrutura vira névoa. O tom é o cinza-cera da identidade, não
   o mármore — as arestas são contexto, não conteúdo. */
void main() { cor = vec4(0.486, 0.533, 0.573, 0.055); }`;

function programa(gl: WebGL2RenderingContext, vs: string, fs: string): WebGLProgram {
  const compilar = (tipo: number, fonte: string) => {
    const s = gl.createShader(tipo)!;
    gl.shaderSource(s, fonte);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(`shader: ${gl.getShaderInfoLog(s)}`);
    return s;
  };
  const p = gl.createProgram()!;
  gl.attachShader(p, compilar(gl.VERTEX_SHADER, vs));
  gl.attachShader(p, compilar(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(`programa: ${gl.getProgramInfoLog(p)}`);
  return p;
}

/** Acima disto as arestas viram uma mancha e custam caro; só os vértices. */
export const MAX_ARESTAS_DESENHADAS = 2000;

export class Nuvem {
  private gl: WebGL2RenderingContext;
  private progPto: WebGLProgram;
  private progAre: WebGLProgram;
  private bufPos: WebGLBuffer;
  private bufVal: WebGLBuffer;
  private bufAre: WebGLBuffer;
  private texLut: WebGLTexture;
  private nAre = 0;
  private n = 0;
  private centro: [number, number] = [0, 0];
  private meia: [number, number] = [1, 1];
  private cosT = 1;
  private senT = 0;
  private posicoes: Float32Array = new Float32Array(0);

  constructor(private canvas: HTMLCanvasElement) {
    const gl = canvas.getContext('webgl2',
      { antialias: true, alpha: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error('WebGL2 indisponivel neste navegador');
    this.gl = gl;
    this.progPto = programa(gl, VS, FS);
    this.progAre = programa(gl, VS_ARESTA, FS_ARESTA);
    this.bufPos = gl.createBuffer()!;
    this.bufVal = gl.createBuffer()!;
    this.bufAre = gl.createBuffer()!;
    this.texLut = gl.createTexture()!;
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.texLut);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB8, 256, 1, 0, gl.RGB, gl.UNSIGNED_BYTE, lut(256));
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  }

  /** Recebe o layout e a topologia. Chamar quando o grafo muda, não por quadro.
   *  `arestas` é a mesma lista que o exportador usa — não uma segunda travessia
   *  da CSR aqui dentro. */
  rede(xy: Float32Array, n: number, arestas?: Array<[number, number, number]> | null) {
    const gl = this.gl;
    this.n = n;

    /* Eixo principal da nuvem, pela covariância. Girar para alinhá-lo é a
       única transformação que ganha área sem tocar em nenhuma distância. */
    let mx = 0, my = 0;
    for (let i = 0; i < n; ++i) { mx += xy[2 * i]; my += xy[2 * i + 1]; }
    mx /= n; my /= n;
    let cxx = 0, cyy = 0, cxy = 0;
    for (let i = 0; i < n; ++i) {
      const dx = xy[2 * i] - mx, dy = xy[2 * i + 1] - my;
      cxx += dx * dx; cyy += dy * dy; cxy += dx * dy;
    }
    let teta = 0.5 * Math.atan2(2 * cxy, cxx - cyy);
    /* Qual das duas orientações põe o eixo longo na horizontal é questão de
       convenção de sinal, e convenção de sinal é o tipo de coisa que se erra em
       silêncio. Em vez de confiar na dedução, mede-se: gira, compara as
       extensões, e acrescenta 90° se ficou de pé. O canvas é largo. */
    for (let tent = 0; tent < 2; ++tent) {
      const c = Math.cos(teta), sn = Math.sin(teta);
      let ex = 0, ey = 0;
      let xa = Infinity, xb = -Infinity, ya = Infinity, yb = -Infinity;
      for (let i = 0; i < n; ++i) {
        const dx = xy[2 * i] - mx, dy = xy[2 * i + 1] - my;
        const x = c * dx + sn * dy, y = -sn * dx + c * dy;
        if (x < xa) xa = x; if (x > xb) xb = x;
        if (y < ya) ya = y; if (y > yb) yb = y;
      }
      ex = xb - xa; ey = yb - ya;
      if (ex >= ey) break;
      teta += Math.PI / 2;
    }
    this.cosT = Math.cos(teta); this.senT = Math.sin(teta);

    const girado = new Float32Array(2 * n);
    let xlo = Infinity, xhi = -Infinity, ylo = Infinity, yhi = -Infinity;
    for (let i = 0; i < n; ++i) {
      const dx = xy[2 * i] - mx, dy = xy[2 * i + 1] - my;
      const x = this.cosT * dx + this.senT * dy;
      const y = -this.senT * dx + this.cosT * dy;
      girado[2 * i] = x; girado[2 * i + 1] = y;
      if (x < xlo) xlo = x; if (x > xhi) xhi = x;
      if (y < ylo) ylo = y; if (y > yhi) yhi = y;
    }
    if (!Number.isFinite(xlo)) { xlo = ylo = 0; xhi = yhi = 1; }
    this.centro = [(xlo + xhi) / 2, (ylo + yhi) / 2];
    this.meia = [Math.max((xhi - xlo) / 2, 1e-9), Math.max((yhi - ylo) / 2, 1e-9)];
    this.posicoes = girado;
    xy = girado;

    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufPos);
    gl.bufferData(gl.ARRAY_BUFFER, girado, gl.STATIC_DRAW);

    this.nAre = 0;
    if (arestas && n <= MAX_ARESTAS_DESENHADAS) {
      const linhas: number[] = [];
      for (const [i, j] of arestas) {
        if (i === j) continue;
        linhas.push(xy[2 * i], xy[2 * i + 1], xy[2 * j], xy[2 * j + 1]);
      }
      this.nAre = linhas.length / 2;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufAre);
      gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(linhas), gl.STATIC_DRAW);
    }
  }

  desenhar(valores: Float32Array, max: number) {
    const gl = this.gl;
    const cw = this.canvas.clientWidth || 1, ch = this.canvas.clientHeight || 1;
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    if (this.canvas.width !== Math.round(cw * dpr) || this.canvas.height !== Math.round(ch * dpr)) {
      this.canvas.width = Math.round(cw * dpr);
      this.canvas.height = Math.round(ch * dpr);
    }
    const L = this.canvas.width, A = this.canvas.height;
    gl.viewport(0, 0, L, A);
    /* --dd-canvas, e NÃO --dd-ink: #101A24 é o extremo zero do mapa de cor, e
       limpar o fundo com ele faria todo vértice de população desprezível ficar
       exatamente da cor do fundo. A rede sumiria e só o pacote apareceria — o
       que parece um layout que não preenche, quando o layout está certo. O fundo
       precisa ser distinguível de "zero". */
    gl.clearColor(0.043, 0.071, 0.098, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
    if (this.n === 0) return;

    /* Margem FIXA em pixels, escala UNIFORME: o eixo limitante encosta na
       margem e o outro fica centrado. Ver o cabeçalho deste arquivo. */
    const margem = 18 * dpr;
    const k = Math.min((L / 2 - margem) / this.meia[0], (A / 2 - margem) / this.meia[1]);
    const escala: [number, number] = [(2 * k) / L, (2 * k) / A];


    if (this.nAre > 0) {
      gl.useProgram(this.progAre);
      gl.uniform2f(gl.getUniformLocation(this.progAre, 'uEscala'), escala[0], escala[1]);
      gl.uniform2f(gl.getUniformLocation(this.progAre, 'uCentro'), this.centro[0], this.centro[1]);
      const a = gl.getAttribLocation(this.progAre, 'aPos');
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufAre);
      gl.enableVertexAttribArray(a);
      gl.vertexAttribPointer(a, 2, gl.FLOAT, false, 0, 0);
      gl.drawArrays(gl.LINES, 0, this.nAre);
    }

    gl.useProgram(this.progPto);
    gl.uniform2f(gl.getUniformLocation(this.progPto, 'uEscala'), escala[0], escala[1]);
    gl.uniform2f(gl.getUniformLocation(this.progPto, 'uCentro'), this.centro[0], this.centro[1]);
    gl.uniform1f(gl.getUniformLocation(this.progPto, 'uMax'), max);
    gl.uniform1f(gl.getUniformLocation(this.progPto, 'uRaio'),
                 Math.max(3, Math.min(16, 340 / Math.sqrt(this.n))) * dpr);
    gl.uniform1i(gl.getUniformLocation(this.progPto, 'uLut'), 0);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.texLut);

    const aPos = gl.getAttribLocation(this.progPto, 'aPos');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufPos);
    gl.enableVertexAttribArray(aPos);
    gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);

    const aVal = gl.getAttribLocation(this.progPto, 'aVal');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.bufVal);
    gl.bufferData(gl.ARRAY_BUFFER, valores, gl.DYNAMIC_DRAW);
    gl.enableVertexAttribArray(aVal);
    gl.vertexAttribPointer(aVal, 1, gl.FLOAT, false, 0, 0);

    gl.drawArrays(gl.POINTS, 0, this.n);
  }

  /** Vértice mais próximo do clique, ou -1. */
  sitioNoEvento(ev: { clientX: number; clientY: number }): number {
    const r = this.canvas.getBoundingClientRect();
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    const L = this.canvas.width, A = this.canvas.height;
    const px = ((ev.clientX - r.left) / r.width) * L;
    const py = ((ev.clientY - r.top) / r.height) * A;
    const margem = 18 * dpr;
    const k = Math.min((L / 2 - margem) / this.meia[0], (A / 2 - margem) / this.meia[1]);
    let melhor = -1, dist = Infinity;
    for (let i = 0; i < this.n; ++i) {
      const x = L / 2 + (this.posicoes[2 * i] - this.centro[0]) * k;
      const y = A / 2 + (this.posicoes[2 * i + 1] - this.centro[1]) * k;
      const d = (x - px) ** 2 + (y - py) ** 2;
      if (d < dist) { dist = d; melhor = i; }
    }
    return dist < (28 * dpr) ** 2 ? melhor : -1;
  }
}
