$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
Set-Location $repoRoot
$python = Join-Path $repoRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $python)) { $python = 'python' }

# This is intentionally the complete deterministic test/replay path.  It is
# compute-bound and must not be replaced by a sampled or floating-point run.
& $python -m order15_mu3.scripts.verify_pre_q144_audit
if ($LASTEXITCODE -ne 0) { throw 'pre-Q144 exact audit failed' }
& $python -m order15_mu3.scripts.verify_aggregate_upper
if ($LASTEXITCODE -ne 0) { throw 'aggregate integer sieve failed' }
& $python -m pytest -q
if ($LASTEXITCODE -ne 0) { throw 'full pytest suite failed' }
