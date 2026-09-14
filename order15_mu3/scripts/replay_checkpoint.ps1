$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
Set-Location $repoRoot
$python = Join-Path $repoRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { $python = 'python' }

& $python -m order15_mu3.scripts.verify_benchmark
if ($LASTEXITCODE -ne 0) { throw 'benchmark certificate failed' }
& $python -m order15_mu3.scripts.verify_bounds
if ($LASTEXITCODE -ne 0) { throw 'baseline bounds certificate failed' }
& $python -m order15_mu3.scripts.verify_pre_q144_audit
if ($LASTEXITCODE -ne 0) { throw 'pre-Q144 exact audit failed' }
& $python -m order15_mu3.scripts.verify_q144_boundary
if ($LASTEXITCODE -ne 0) { throw 'Q=144 boundary certificate failed' }
& $python -m order15_mu3.scripts.verify_q150_boundary
if ($LASTEXITCODE -ne 0) { throw 'Q=150 boundary certificate failed' }
& $python -m order15_mu3.scripts.verify_aggregate_upper
if ($LASTEXITCODE -ne 0) { throw 'aggregate integer sieve failed' }
