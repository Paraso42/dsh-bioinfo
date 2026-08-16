# download-af2-params.ps1 — AF2/AF2-Multimer params -> COLABFOLDDIR (D:\bioai\models\colabfold)
# Usage:
#   pwsh -File D:\bioai\bin\download-af2-params.ps1                 # direct GCS (blocked on this network)
#   pwsh -File D:\bioai\bin\download-af2-params.ps1 -Proxy http://127.0.0.1:7897   # via local proxy (enable your proxy first)
#   pwsh -File D:\bioai\bin\download-af2-params.ps1 -MonomerOnly | -MultimerOnly
# NOTE: ASCII-only comments.
param([string]$Proxy = '', [switch]$MonomerOnly, [switch]$MultimerOnly)

$ErrorActionPreference = 'Stop'
$dest = 'D:\bioai\models\colabfold'
New-Item -ItemType Directory -Force -Path $dest | Out-Null
$dl = 'D:\bioai\bin\parallel-download.ps1'

$files = @(
    @{ url = 'https://storage.googleapis.com/alphafold-colabfold/alphafold2_ptm_2.3.1'; name = 'alphafold2_ptm_2.3.1' },
    @{ url = 'https://storage.googleapis.com/alphafold-colabfold/multimer_v3'; name = 'multimer_v3' }
)

foreach ($f in $files) {
    if ($MonomerOnly -and $f.name -ne 'alphafold2_ptm_2.3.1') { continue }
    if ($MultimerOnly -and $f.name -ne 'multimer_v3') { continue }
    $out = Join-Path $dest $f.name
    if ((Test-Path $out) -and ((Get-Item $out).Length -gt 1GB)) {
        Write-Output ("skip existing ({0:N1} GB): {1}" -f ((Get-Item $out).Length/1GB), $out)
        continue
    }
    Write-Output ("downloading {0} -> {1} (proxy: {2})" -f $f.name, $out, $(if ($Proxy) { $Proxy } else { 'direct' }))
    if ($Proxy) { & $dl -Url $f.url -Out $out -Chunks 8 -TimeoutSec 14400 -Proxy $Proxy }
    else { & $dl -Url $f.url -Out $out -Chunks 8 -TimeoutSec 14400 }
}
Write-Output 'AF2 PARAMS DONE'
