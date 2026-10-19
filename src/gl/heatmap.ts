/* heatmap.ts — a rede desenrolada em WebGL.
 *
 * ESCALA POR EIXO, de propósito, e a razão é física. A rede desenrolada é uma
 * GRADE DE ÍNDICES: (m, q) são posição ao longo do eixo e número do
 * protofilamento, e a distância no plano do desenho não corresponde a distância
 * nenhuma — nada é medido nela. Então esticar cada eixo até encher a caixa não
 * distorce informação, só usa o espaço.
 *
 * O layout espectral (src/gl/nuvem.ts) faz o oposto e mantém escala uniforme,
 * porque ali distância É a informação: é o embedding que separa os módulos.
 *
 * A escala de cor vem de src/nucleo/paleta.ts, a mesma que os testes verificam:
 * o shader monta a LUT a partir dela em vez de ter uma cópia própria. Um
 * colormap invertido é um dos dois jeitos de o heatmap ficar bonito e mentir; o
 * outro é a transposição, que mora em src/nucleo/indices.ts.
 */
import { lut } from '../nucleo/paleta';

const VS = `#version 300 es
const vec2 P[3] = vec2[3](vec2(-1.0,-1.0), vec2(3.0,-1.0), vec2(-1.0,3.0));
out vec2 vUV;
void main() {
  vec2 p = P[gl_VertexID];
  gl_Position = vec4(p, 0.0, 1.0);
  vUV = vec2((p.x + 1.0) * 0.5, 1.0 - (p.y + 1.0) * 0.5);
}`;

const FS = `#version 300 es
precision highp float;
uniform sampler2D uDados;
uniform sampler2D uLut;
uniform float uMax;
in vec2 vUV;
out vec4 cor;
void main() {
  float v = texture(uDados, vUV).r;
  if (v < 0.0) { cor = vec4(0.10, 0.11, 0.13, 1.0); return; }  /* fora do grafo */
  float t = uMax > 0.0 ? clamp(v / uMax, 0.0, 1.0) : 0.0;
  cor = vec4(texture(uLut, vec2(t, 0.5)).rgb, 1.0);
}`;

function compilar(gl: WebGL2RenderingContext, tipo: number, fonte: string) {
  const s = gl.createShader(tipo)!;
  gl.shaderSource(s, fonte);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
    throw new Error(`shader: ${gl.getShaderInfoLog(s)}`);
  }
  return s;
}

export class Heatmap {
  private gl: WebGL2RenderingContext;
  private prog: WebGLProgram;
  private texDados: WebGLTexture;
  private texLut: WebGLTexture;
  private uMax: WebGLUniformLocation;
  largura = 1;
  altura = 1;

  constructor(private canvas: HTMLCanvasElement) {
    /* preserveDrawingBuffer: o teste de fumaça precisa ler os pixels de volta
       para provar que o mapa não está preto — e sem isto o buffer vem limpo.
       De quebra, o usuário consegue salvar a figura pelo menu do navegador. */
    const gl = canvas.getContext('webgl2',
      { antialias: false, alpha: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error('WebGL2 indisponivel neste navegador');
    this.gl = gl;
    this.prog = gl.createProgram()!;
    gl.attachShader(this.prog, compilar(gl, gl.VERTEX_SHADER, VS));
    gl.attachShader(this.prog, compilar(gl, gl.FRAGMENT_SHADER, FS));
    gl.linkProgram(this.prog);
    if (!gl.getProgramParameter(this.prog, gl.LINK_STATUS)) {
      throw new Error(`programa: ${gl.getProgramInfoLog(this.prog)}`);
    }
    gl.useProgram(this.prog);
    gl.uniform1i(gl.getUniformLocation(this.prog, 'uDados'), 0);
    gl.uniform1i(gl.getUniformLocation(this.prog, 'uLut'), 1);
    this.uMax = gl.getUniformLocation(this.prog, 'uMax')!;

    this.texDados = gl.createTexture()!;
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.texDados);
    /* NEAREST de propósito: interpolar inventaria valores entre sítios, e um
       sítio é uma coisa discreta. */
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    this.texLut = gl.createTexture()!;
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.texLut);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    const tabela = lut(256);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB8, 256, 1, 0, gl.RGB, gl.UNSIGNED_BYTE, tabela);
  }

  desenhar(dados: Float32Array, largura: number, altura: number, max: number) {
    const gl = this.gl;
    const cw = this.canvas.clientWidth || 1;
    const ch = this.canvas.clientHeight || 1;
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    if (this.canvas.width !== Math.round(cw * dpr) || this.canvas.height !== Math.round(ch * dpr)) {
      this.canvas.width = Math.round(cw * dpr);
      this.canvas.height = Math.round(ch * dpr);
    }
    this.largura = largura;
    this.altura = altura;
    gl.viewport(0, 0, this.canvas.width, this.canvas.height);
    gl.useProgram(this.prog);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.texDados);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, largura, altura, 0, gl.RED, gl.FLOAT, dados);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.texLut);
    gl.uniform1f(this.uMax, max);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }

  /** Pixel do canvas → (m, q). Usa a MESMA orientação do shader. */
  texelDoEvento(ev: { clientX: number; clientY: number }): { m: number; q: number } | null {
    const r = this.canvas.getBoundingClientRect();
    const x = (ev.clientX - r.left) / r.width;
    const y = (ev.clientY - r.top) / r.height;
    if (x < 0 || x >= 1 || y < 0 || y >= 1) return null;
    return { m: Math.floor(x * this.largura), q: Math.floor(y * this.altura) };
  }
}
