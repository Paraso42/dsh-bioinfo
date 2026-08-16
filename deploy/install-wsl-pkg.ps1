# install-wsl-pkg.ps1 — install the WSL component bundle (ASCII only)
$ErrorActionPreference = 'Stop'
$pkg = 'D:\bioai\wsl\Microsoft.WSL_2.7.11.0_x64_ARM64.msixbundle'
if (-not (Test-Path $pkg)) { Write-Error "missing $pkg"; exit 1 }
Write-Output "installing $pkg ..."
Add-AppxPackage -Path $pkg
Write-Output "Add-AppxPackage completed"
