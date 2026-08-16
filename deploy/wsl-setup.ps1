# wsl-setup.ps1 — one-shot: import pre-downloaded Ubuntu + LocalColabFold + LightDock
# PREREQ (admin, once, then REBOOT — no download needed, local component store):
#   dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
#   dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
# Usage after reboot (no admin needed):
#   pwsh -File D:\bioai\bin\wsl-setup.ps1
# NOTE: ASCII-only comments (PowerShell on this box mis-decodes UTF-8 ps1 comments).
param([string]$Distro = 'Ubuntu')

$ErrorActionPreference = 'Stop'

# 1) WSL features enabled? (wsl prints an install hint when they are not)
$wslOut = wsl -l -v 2>&1 | Out-String
if ($wslOut -match 'i\s*n\s*s\s*t\s*a\s*l\s*l' -and $wslOut -notmatch 'docker-desktop') {
    Write-Error @"
WSL features are not enabled yet. In an ADMIN PowerShell run:
  dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
  dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
Then REBOOT and re-run this script. (No slow download involved: DISM enables from the local component store.)
"@
    exit 1
}

$distros = @((wsl -l -q 2>$null | Out-String).Trim() -split "`r?`n" | ForEach-Object { $_.Trim() } | Where-Object { $_ -and $_ -notmatch '^\x00' })

# 2) import pre-downloaded rootfs (TUNA mirror) when the distro is missing
if (-not ($distros -contains $Distro)) {
    $rootfs = 'D:\bioai\wsl\ubuntu-jammy.rootfs.tar.gz'
    if (-not (Test-Path $rootfs)) {
        Write-Error "missing $rootfs (re-download: parallel-download.ps1 -Url https://mirrors.tuna.tsinghua.edu.cn/ubuntu-cloud-images/wsl/jammy/current/ubuntu-jammy-wsl-amd64-ubuntu22.04lts.rootfs.tar.gz -Out $rootfs -Chunks 8)"
        exit 1
    }
    $installDir = "D:\bioai\wsl\$Distro"
    Write-Output "importing $Distro from $rootfs -> $installDir (vhdx on D: drive)"
    wsl --import $Distro $installDir $rootfs
    if ($LASTEXITCODE -ne 0) { Write-Error 'wsl --import failed (maybe reboot still pending?)'; exit 1 }
}

# 3) run the Linux-side bootstrap (Miniconda + colabfold env + LightDock)
$boot = 'D:\bioai\bin\wsl-bootstrap.sh'
if (-not (Test-Path $boot)) { Write-Error "missing $boot"; exit 1 }
Write-Output "running wsl-bootstrap.sh in $Distro (log: D:\bioai\pip-logs\wsl-bootstrap.log)"
New-Item -ItemType Directory -Force -Path 'D:\bioai\pip-logs' | Out-Null
wsl -d $Distro -- bash -lc "bash /mnt/d/bioai/bin/wsl-bootstrap.sh 2>&1 | tee /mnt/d/bioai/pip-logs/wsl-bootstrap.log; exit `${PIPESTATUS[0]}"
$rc = $LASTEXITCODE
$logOk = Select-String -Path 'D:\bioai\pip-logs\wsl-bootstrap.log' -Pattern 'BOOTSTRAP DONE' -Quiet -ErrorAction SilentlyContinue
Write-Output "bootstrap exit: $rc ; BOOTSTRAP DONE marker: $logOk"
if (-not $logOk) { exit 1 }
exit 0
