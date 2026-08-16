# install-local-msa.ps1 — LocalColabFold local MSA databases installer
# Purpose: end dependence on the colabfold MMseqs2 server (unreachable 2026-08:
# official API blocked, backup site abandoned). Installs mmseqs2 (in WSL) plus
# UniRef30 + colabfold_envdb (~70 GB total) so colabfold_search runs locally.
#
# Usage:
#   pwsh -File install-local-msa.ps1 -InstallMmseqs
#   pwsh -File install-local-msa.ps1 -DownloadDb              # hours; resume-safe
#   pwsh -File install-local-msa.ps1 -DownloadDb -Proxy http://127.0.0.1:7897
#   pwsh -File install-local-msa.ps1 -Verify
#
# Source status tested 2026-08 from this machine:
#   GWDG     https://wwwuser.gwdg.de/~compbiol/colabfold/    REACHABLE (file-level OK)
#   upstream https://colabfold.steineggerlab.workers.dev/    TIMEOUT (canonical, keep as fallback)
# NOTE: keep this file pure ASCII.
param(
    [string]$DataDir = 'D:\bioai\msa-db',
    [switch]$InstallMmseqs,
    [switch]$DownloadDb,
    [switch]$Verify,
    [string]$Proxy = '',
    [int]$Chunks = 8,
    [string]$ParallelScript = 'D:\bioai\dsh-bioinfo\deploy\parallel-download.ps1'
)

$ErrorActionPreference = 'Stop'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$GWDG   = 'https://wwwuser.gwdg.de/~compbiol/colabfold'
$UPSTRM = 'https://colabfold.steineggerlab.workers.dev'
$EnvDb  = 'colabfold_envdb_202108.tar.gz'   # ~9 GB
$UniRef = 'uniref30_2302.tar.gz'            # ~60 GB
$CondaForgeTUNA = 'https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge'

function Invoke-Download([string]$url, [string]$out) {
    $pargs = @{ Url = $url; Out = $out; Chunks = $Chunks; ChunkDir = "$DataDir\.chunks"; TimeoutSec = 7200 }
    if ($Proxy) { $pargs.Proxy = $Proxy }
    & $ParallelScript @pargs
    if ($LASTEXITCODE -ne 0) { throw "download failed: $url" }
}

if ($InstallMmseqs) {
    $bash = @'
set -e
if command -v mmseqs >/dev/null 2>&1; then
  echo "mmseqs already present:"; mmseqs version | head -1; exit 0
fi
CONDA=""
for c in conda ~/miniforge3/bin/conda ~/miniconda3/bin/conda; do
  if command -v "$c" >/dev/null 2>&1 || [ -x "$c" ]; then CONDA="$c"; break; fi
done
if [ -n "$CONDA" ]; then
  echo "installing mmseqs2 via $CONDA (TUNA conda-forge mirror)"
  "$CONDA" install -y -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge mmseqs2
else
  echo "no conda found; trying apt (may prompt for sudo password)"
  sudo apt-get update && sudo apt-get install -y mmseqs2
fi
mmseqs version | head -1
'@
    wsl -d Ubuntu -- bash -lc $bash
}

if ($DownloadDb) {
    New-Item -ItemType Directory -Force -Path $DataDir | Out-Null
    foreach ($file in @($EnvDb, $UniRef)) {
        $out = Join-Path $DataDir $file
        if ((Test-Path $out) -and (Get-Item $out).Length -gt 1MB) {
            Write-Output "resume: $file"
        }
        try {
            Invoke-Download "$GWDG/$file" $out
        } catch {
            Write-Output "GWDG failed ($($_.Exception.Message)); trying upstream (may timeout)"
            Invoke-Download "$UPSTRM/$file" $out
        }
    }
    Write-Output "downloads done. extract in WSL:"
    Write-Output "  wsl -d Ubuntu -- bash -lc 'cd /mnt/d/bioai/msa-db && tar xzf $EnvDb && tar xzf $UniRef'"
}

if ($Verify) {
    Write-Output "--- DB directories ---"
    wsl -d Ubuntu -- bash -lc "test -d /mnt/d/bioai/msa-db/uniref30_2302 && echo uniref30_2302: OK || echo uniref30_2302: MISSING; test -d /mnt/d/bioai/msa-db/colabfold_envdb_202108 && echo colabfold_envdb_202108: OK || echo colabfold_envdb_202108: MISSING"
    Write-Output "--- mmseqs ---"
    wsl -d Ubuntu -- bash -lc "mmseqs version 2>/dev/null | head -1 || echo 'mmseqs: NOT FOUND (run -InstallMmseqs)'"
}

Write-Output @'

Two-step local-MSA recipe (see protein-modeling SKILL, section 1):
  1) search:   colabfold_search --db1 <DataDir>/uniref30_2302 --db2 <DataDir>/colabfold_envdb_202108 query.fasta out_msa
  2) predict:  colabfold_batch out_msa result_dir --model-type alphafold2_multimer_v3
(WSL paths: /mnt/d/bioai/msa-db/... ; activate the colabfold conda env first)
'@
