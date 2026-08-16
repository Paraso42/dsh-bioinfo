# run_colabfold.ps1 — LocalColabFold 1.5.5 (AF2 / AF2-Multimer) runner on WSL2 GPU
# Prereq: D:\bioai deployment (WSL2 colabfold env + params provisioned in BOTH
#   COLABFOLDDIR/params and XDG_CACHE_HOME/colabfold/params — see deploy-plan.md).
# Usage:
#   pwsh -File run_colabfold.ps1 -Fasta query.fasta -OutDir D:\bioai\jobs\af2_out
#   pwsh -File run_colabfold.ps1 -Fasta complex.fasta -ModelType alphafold2_multimer_v3
# Multi-chain complex: join chains with ':' in one fasta record.
# NOTE: keep this file pure ASCII — PowerShell on this box mis-decodes UTF-8 ps1 comments.
param(
    [Parameter(Mandatory = $true)][string]$Fasta,
    [string]$OutDir = 'D:\bioai\jobs\af2_out',
    [ValidateSet('auto', 'alphafold2', 'alphafold2_ptm', 'alphafold2_multimer_v1', 'alphafold2_multimer_v2', 'alphafold2_multimer_v3', 'deepfold_v1')][string]$ModelType = 'alphafold2_multimer_v3',
    [int]$NumModels = 1,
    [int]$NumRecycle = 3,
    [string]$MsaMode = 'mmseqs2_uniref_env',
    [switch]$SavePAE,
    [switch]$DryRun
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path -LiteralPath $Fasta)) { Write-Error "fasta not found: $Fasta"; exit 1 }

$distros = @((wsl -l -q 2>$null | Out-String).Trim() -split "`r?`n" | ForEach-Object { $_.Trim() } | Where-Object { $_ -and $_ -notmatch '^\x00' })
$distro = $distros | Select-Object -First 1
if (-not $distro) {
    Write-Error @"
WSL distro not found. Deployment guide: D:\bioai\deploy-plan.md
(admin) wsl --install -d Ubuntu --no-launch  ->  reboot  ->  pwsh -File D:\bioai\bin\wsl-setup.ps1
"@
    exit 1
}

function To-WslPath([string]$p) {
    if ($p -match '^([A-Za-z]):(.*)$') { return '/mnt/' + $Matches[1].ToLower() + ($Matches[2] -replace '\\', '/') }
    return $p
}

$wFasta = To-WslPath (Resolve-Path $Fasta)
$wOut = To-WslPath $OutDir
$paeFlag = $(if ($SavePAE) { '--save-pae-plot' } else { '' })

$inner = @(
    'export HOME=/root',
    'cd /root',
    'unset PYTHONPATH',
    'export COLABFOLDDIR=/mnt/d/bioai/models/colabfold',
    'export XDG_CACHE_HOME=/mnt/d/bioai/models/cache',
    'export XLA_PYTHON_CLIENT_PREALLOCATE=false',
    'export XLA_PYTHON_CLIENT_MEM_FRACTION=0.85',
    'source /root/miniforge3/etc/profile.d/conda.sh',
    'conda activate colabfold',
    ("mkdir -p `"{0}`"" -f $wOut),
    ("colabfold_batch --model-type {0} --num-models {1} --num-recycle {2} --disable-unified-memory --msa-mode `"{3}`" {4} `"{5}`" `"{6}`"" -f $ModelType, $NumModels, $NumRecycle, $MsaMode, $paeFlag, $wFasta, $wOut)
) -join '; '

Write-Host "[wsl:$distro] colabfold_batch --model-type $ModelType --msa-mode $MsaMode -> $OutDir"
if ($DryRun) { Write-Host $inner; exit 0 }
wsl -d $distro -- bash -lc $inner
$rc = $LASTEXITCODE
Write-Host "colabfold(wsl) exit: $rc (output: $OutDir)"
exit $rc
