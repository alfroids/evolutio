panel convert index.py \
  --to pyodide-worker \
  --out docs/ \
  --resources src/evolution_models/*.py tabs/*.py \
  --requirements numpy plotly pandas