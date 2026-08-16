# verify-layout.ps1 — self-consistency check for an installed dsh-bioinfo replica.
# Verifies that every file the preset composition and the protein-tools plugin
# reference exists under the canonical deployment layout, then (optionally) runs
# the schema validator from this repository.
# Usage after docs/INSTALL.md:
#   pwsh -File scripts\verify-layout.ps1
#   pwsh -File scripts\verify-layout.ps1 -PresetDir C:\my\dsh\.agent-presets\bioinfo -BioaiRoot E:\bioai
# NOTE: ASCII-only comments.
param(
    [string]$PresetDir = 'C:\deepseek-harness\.dsh\.agent-presets\bioinfo',
    [string]$BioaiRoot = 'D:\bioai',
    [string]$Python313 = 'C:\Program Files\Python313\python.exe',
    [string]$BioPythonDir = 'D:\biopython',
    [switch]$SkipNodeTest
)

$ErrorActionPreference = 'Continue'
$script:fail = 0

function Check([string]$label, [string]$path) {
    if ($path -and (Test-Path -LiteralPath $path)) {
        Write-Output ("[ok]   " + $label)
    } else {
        Write-Output ("[FAIL] " + $label + "  ->  " + $path)
        $script:fail++
    }
}

Write-Output '== preset layer =='
Check 'preset.yml' (Join-Path $PresetDir 'preset.yml')
Check 'agent.cordis.yml' (Join-Path $PresetDir 'agent.cordis.yml')
Check 'plugins\protein-tools.js' (Join-Path $PresetDir 'plugins\protein-tools.js')

Write-Output '== skills layer =='
$skills = 'biopython', 'biopython-analyses', 'protein-modeling', 'protein-quality', 'chem-informatics', 'bio-data-hub', 'bio-visualization'
foreach ($s in $skills) {
    Check ("skills\$s\SKILL.md") (Join-Path $PresetDir ("skills\" + $s + "\SKILL.md"))
}

$res = Join-Path $PresetDir 'skills\protein-modeling\resources'
$resPq = Join-Path $PresetDir 'skills\protein-quality\resources'
$resCi = Join-Path $PresetDir 'skills\chem-informatics\resources'

Write-Output '== plugin backends (default BIO_TOOLS_* paths) =='
Check 'esmfold_api.py' (Join-Path $res 'esmfold_api.py')
Check 'pp_interact.py' (Join-Path $res 'pp_interact.py')
Check 'vina_dock.py' (Join-Path $res 'vina_dock.py')
Check 'run_colabfold.ps1' (Join-Path $res 'run_colabfold.ps1')
Check 'md_mmgbsa.py' (Join-Path $res 'md_mmgbsa.py')
Check 'pdb_to_pdbqt.py' (Join-Path $res 'pdb_to_pdbqt.py')
Check 'struct_eval.py' (Join-Path $resPq 'struct_eval.py')
Check 'prodigy_affinity.py' (Join-Path $resPq 'prodigy_affinity.py')
Check 'virtual_screen.py' (Join-Path $resCi 'virtual_screen.py')
Check 'mol_tools.py' (Join-Path $resCi 'mol_tools.py')

Write-Output '== environment layer =='
Check 'python 3.13 interpreter' $Python313
Check 'venv interpreter' (Join-Path $BioaiRoot 'venv\Scripts\python.exe')
Check 'venv-esm interpreter' (Join-Path $BioaiRoot 'venv-esm\Scripts\python.exe')
Check 'deploy scripts in <BioaiRoot>\bin' (Join-Path $BioaiRoot 'bin\wsl-setup.ps1')
Check 'biopython package dir' (Join-Path $BioPythonDir 'Bio')
Check 'jobs dir' (Join-Path $BioaiRoot 'jobs')
Check 'models dir' (Join-Path $BioaiRoot 'models')

$vina = Get-ChildItem (Join-Path $BioaiRoot 'bin') -Filter 'vina*.exe' -ErrorAction SilentlyContinue | Select-Object -First 1
if ($vina) { Write-Output ("[ok]   vina binary (" + $vina.Name + ")") }
else { Write-Output ("[FAIL] vina binary (vina*.exe under " + (Join-Path $BioaiRoot 'bin') + ")"); $script:fail++ }

if (-not $SkipNodeTest) {
    $validator = Join-Path $PSScriptRoot '..\validate-plugin-schemas.js'
    if (Test-Path $validator) {
        Write-Output '== schema validator =='
        $node = Get-Command node -ErrorAction SilentlyContinue
        if ($node) {
            $testPath = (Join-Path $PSScriptRoot '..\plugins\protein-tools.js')
            node --check $testPath
            if ($LASTEXITCODE -ne 0) { $script:fail++ }
            node $validator
            if ($LASTEXITCODE -ne 0) { $script:fail++ }
        } else {
            Write-Output '[FAIL] node not found — cannot run schema validator'
            $script:fail++
        }
    } else {
        Write-Output ('[warn] validator not found: ' + $validator + ' (run from a dsh-bioinfo checkout)')
    }
}

Write-Output ''
if ($script:fail -eq 0) { Write-Output 'LAYOUT OK - replica is self-consistent.'; exit 0 }
Write-Output ("LAYOUT INCOMPLETE — $script:fail check(s) failed."); exit 1
