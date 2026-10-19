#!/bin/sh
# Constrói, sobe o servidor de pré-visualização, roda o teste de fumaça e
# derruba o servidor — inclusive se o teste falhar.
set -eu
RAIZ=$(cd "$(dirname "$0")/.." && pwd)
cd "$RAIZ"
npx vite build > /dev/null
npx vite preview --port 4173 --strictPort > /dev/null 2>&1 &
SERVIDOR=$!
trap 'kill $SERVIDOR 2>/dev/null || true' EXIT INT TERM
for _ in 1 2 3 4 5 6 7 8 9 10; do
  if curl -sf -o /dev/null http://localhost:4173/; then break; fi
  sleep 1
done
node tools/fumaca.mjs http://localhost:4173/
# Mesmo servidor, segunda pergunta: a vista do tubo desenha, e GIRAR muda a
# figura. A segunda metade e que e o teste — "o tubo aparece" passaria tambem
# com uma projecao sem profundidade.
node tools/fumaca_tubo.mjs http://localhost:4173/
