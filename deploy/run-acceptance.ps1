# run-acceptance.ps1 — AF2-Multimer end-to-end acceptance (barnase-barstar)
# Prereqs: colabfold env installed in WSL + params at <BioaiRoot>\models\colabfold\params
# Usage: pwsh -File run-acceptance.ps1            # canonical layout (D:\bioai)
#        pwsh -File run-acceptance.ps1 -BioaiRoot E:\bioai -PresetDir C:\my\dsh\presets\bioinfo
# NOTE: ASCII-only comments.
param(
    [string]$MsaMode = 'mmseqs2_uniref_env',
    [int]$NumModels = 1,
    [int]$NumRecycle = 3,
    [string]$BioaiRoot = 'D:\bioai',
    [string]$PresetDir = 'C:\deepseek-harness\.dsh\.agent-presets\bioinfo',
    [string]$Python313 = 'C:\Program Files\Python313\python.exe',
    [string]$BioPythonDir = 'D:\biopython'
)

$ErrorActionPreference = 'Continue'

# WSL mount path for the BioaiRoot drive (D:\bioai -> /mnt/d/bioai)
$drive = $BioaiRoot.Substring(0, 1).ToLowerInvariant()
$wslRoot = '/mnt/' + $drive + '/' + ($BioaiRoot.Substring(2).Trim('\').Replace('\', '/'))

$fasta = Join-Path $BioaiRoot 'jobs\acceptance\1brs_complex.fasta'
$outdir = Join-Path $BioaiRoot 'jobs\acceptance\af2_out'
$params = Join-Path $BioaiRoot 'models\colabfold\params'

if (-not (Test-Path $fasta)) { Write-Error "missing fasta: $fasta"; exit 1 }
if (-not (Test-Path (Join-Path $params 'alphafold_params_colab_2022-03-02.tar'))) {
    Write-Error "multimer params not downloaded yet ($params) — wait for wsl-download-params.ps1"
    exit 1
}

$cmd = @(
  "export HOME=/root",
  "export COLABFOLDDIR=$wslRoot/models/colabfold",
  "export XDG_CACHE_HOME=$wslRoot/models/cache",
  "cd /root",
  "unset PYTHONPATH",
  "export XLA_PYTHON_CLIENT_PREALLOCATE=false",
  "export XLA_PYTHON_CLIENT_MEM_FRACTION=0.85",
  "source /root/miniforge3/etc/profile.d/conda.sh",
  "conda activate colabfold",
  "cd $wslRoot/models/colabfold/params",
  "if [ ! -f download_complexes_multimer_v3_finished.txt ]; then tar -xf alphafold_params_colab_2022-12-06.tar 2>/dev/null; touch download_complexes_multimer_v3_finished.txt; fi",
  "cd $wslRoot/jobs/acceptance",
  "colabfold_batch --model-type alphafold2_multimer_v3 --num-models $NumModels --num-recycle $NumRecycle --disable-unified-memory --msa-mode `"$MsaMode`" $wslRoot/jobs/acceptance/1brs_complex.fasta $wslRoot/jobs/acceptance/af2_out 2>&1 | tail -30"
) -join "; "

Write-Output "running colabfold_batch (multimer) inside WSL..."
wsl -d Ubuntu -- bash -lc $cmd
$rc = $LASTEXITCODE
Write-Output "colabfold exit: $rc"

$ranked = Get-ChildItem (Join-Path $outdir '*ranked*.pdb') -ErrorAction SilentlyContinue | Select-Object -First 1
if (-not $ranked) { Write-Output "NO ranked PDB produced"; exit 1 }

Write-Output "predicted complex: $($ranked.FullName)"

# Biopython interface analysis on the prediction (chain A = barnase, B = barstar)
$pp = Join-Path $PresetDir 'skills\protein-modeling\resources\pp_interact.py'
$jsonOut = Join-Path $outdir 'predicted_interface.json'
$env:PYTHONPATH = $BioPythonDir
& $Python313 $pp --complex $ranked.FullName --chains A B --cutoff 5.0 --out $jsonOut 2>&1 | Select-Object -Last 15

Write-Output '--- comparison vs crystal baseline (A-D: 55 contacts, BSA 1280.9 A^2) ---'
$j = Get-Content $jsonOut -Raw | ConvertFrom-Json
Write-Output ("predicted contacts={0} BSA={1} A^2" -f $j.n_contacts, $j.bsa_total)
Write-Output 'ACCEPTANCE RUN DONE'
