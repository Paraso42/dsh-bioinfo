# parallel-download.ps1 v2 — ranged parallel downloader with resume support
# Per-connection throttle bypass via parallel Range requests; chunks persist in
# -ChunkDir across restarts and resume with curl -C -.
# Usage: pwsh -File parallel-download.ps1 -Url <url> -Out <file> [-Chunks 8] [-TimeoutSec 3600] [-ChunkDir <dir>] [-Base <name>]
# NOTE: keep this file pure ASCII.
param(
    [Parameter(Mandatory = $true)][string]$Url,
    [Parameter(Mandatory = $true)][string]$Out,
    [int]$Chunks = 8,
    [int]$TimeoutSec = 3600,
    [string]$Proxy = '',
    [string]$ChunkDir = '',
    [string]$Base = '',
    [string]$ResolveIp = '',
    [string]$ResolveIps = ''
)

$ErrorActionPreference = 'Stop'
if (-not (Get-Command curl.exe -ErrorAction SilentlyContinue)) { Write-Error 'curl.exe required'; exit 1 }
$px = @()
if ($Proxy) { $px = @('-x', $Proxy) }
$hostPart = ([System.Uri]$Url).Host
$portPart = ([System.Uri]$Url).Port
$ips = @()
if ($ResolveIps) { $ips = @($ResolveIps -split ',') } elseif ($ResolveIp) { $ips = @($ResolveIp) }

# 1) probe total size
$head = curl.exe -sI --max-time 30 $Url 2>$null
$len = 0
foreach ($line in @($head)) {
    if ($line -match '(?im)^content-length:\s*(\d+)') { $len = [long]$Matches[1]; break }
}
if ($len -le 0) {
    Write-Output 'no content-length; single-stream fallback'
    curl.exe @px -s -L -o $Out --max-time $TimeoutSec $Url
    exit $LASTEXITCODE
}

# already complete?
if ((Test-Path $Out) -and ((Get-Item $Out).Length -eq $len)) {
    Write-Output "already complete: $Out"
    exit 0
}

Write-Output ("total {0:N1} MB | {1} chunks" -f ($len / 1MB), $Chunks)

# 2) chunk storage (persistent when -ChunkDir given)
if ($ChunkDir) {
    New-Item -ItemType Directory -Force -Path $ChunkDir | Out-Null
    $dir = $ChunkDir
} else {
    $dir = Join-Path $env:TEMP ('pdl_' + [guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $dir | Out-Null
}
if (-not $Base) { $Base = [System.IO.Path]::GetFileName($Out) }

$chunkSize = [math]::Ceiling($len / $Chunks)
$jobs = @()
$expect = @{}
for ($i = 0; $i -lt $Chunks; $i++) {
    $start = $i * $chunkSize
    if ($start -ge $len) { break }
    $end = [math]::Min($len - 1, $start + $chunkSize - 1)
    $want = $end - $start + 1
    $part = Join-Path $dir ($Base + '.chunk' + ('{0:D3}' -f $i))
    $expect[$part] = $want
    $chunkRsv = @()
    if ($ips.Count -gt 0) {
        $ip = $ips[$i % $ips.Count]
        $chunkRsv = @('--resolve', "$hostPart`:$portPart`:$ip")
    }
    $jobs += Start-Job -ScriptBlock {
        param($u, $s, $e, $p, $t, $px, $rsv, $want)
        for ($a = 1; $a -le 60; $a++) {
            $got = 0
            if (Test-Path $p) { $got = (Get-Item $p).Length }
            if ($got -eq $want) { return }
            if ($got -gt $want) { Remove-Item -LiteralPath $p -Force -ErrorAction SilentlyContinue; $got = 0 }
            $rs = $s + $got
            curl.exe @px @rsv -s -r "$rs-$e" -o $p --max-time $t $u 2>$null
            if ((Test-Path $p) -and ((Get-Item $p).Length -eq $want)) { return }
            Start-Sleep -Seconds 3
        }
    } -ArgumentList $Url, $start, $end, $part, $TimeoutSec, $px, $chunkRsv, $want
}
$null = Wait-Job $jobs -Timeout $TimeoutSec
$failed = @($jobs | Where-Object { $_.State -ne 'Completed' })
Get-Job | Remove-Job -Force -ErrorAction SilentlyContinue
if ($failed.Count -gt 0) {
    Write-Output "$($failed.Count) chunk(s) still incomplete; run the same command again to resume"
    exit 1
}

# 3) assemble
$outDir = Split-Path $Out -Parent
if ($outDir) { New-Item -ItemType Directory -Force -Path $outDir | Out-Null }
$fs = [System.IO.File]::Create($Out)
try {
    for ($i = 0; $i -lt $Chunks; $i++) {
        $part = Join-Path $dir ($Base + '.chunk' + ('{0:D3}' -f $i))
        if (Test-Path $part) {
            $bytes = [System.IO.File]::ReadAllBytes($part)
            $fs.Write($bytes, 0, $bytes.Length)
        }
    }
} finally { $fs.Dispose() }

$final = (Get-Item $Out).Length
if ($final -ne $len) {
    Write-Output ("size mismatch: got {0}, expected {1}; chunks kept — rerun to resume" -f $final, $len)
    exit 1
}
Write-Output ("done: {0:N1} MB -> {1}" -f ($final / 1MB), $Out)
# clean used chunks only when they are ours (persistent dir: remove by base)
for ($i = 0; $i -lt $Chunks; $i++) {
    $part = Join-Path $dir ($Base + '.chunk' + ('{0:D3}' -f $i))
    Remove-Item -LiteralPath $part -Force -ErrorAction SilentlyContinue
}
