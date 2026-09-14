$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
Set-Location $repoRoot
$python = Join-Path $repoRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { $python = 'python' }

& $python -m pytest -q `
  tests/test_mu3.py `
  tests/test_mu3_trace.py `
  tests/test_mu3_color_bounds.py `
  tests/test_mu3_q126_k6_minus_edge.py `
  tests/test_mu3_q126_triangle_components.py `
  tests/test_mu3_q126_k4_four_forest.py `
  tests/test_mu3_q126_pure_k5.py `
  tests/test_mu3_q141_energy36_minimal_cross.py `
  tests/test_mu3_q150_boundary.py
if ($LASTEXITCODE -ne 0) { throw "fast μ3 replay failed with exit code $LASTEXITCODE" }
