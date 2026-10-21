# CI-only functional smoke test of the same Windows Triangle binary PHASE installs.
$ErrorActionPreference = 'Stop'
$probeDir = Join-Path $env:RUNNER_TEMP 'phase-triangle-smoke'
New-Item -ItemType Directory -Path $probeDir -Force | Out-Null
$zip = Join-Path $probeDir 'stamps-win64-binaries.zip'
Invoke-WebRequest -Uri 'https://github.com/pyccino/StaMPS/releases/download/windows-port-bins-v1/stamps-win64-binaries.zip' -OutFile $zip
$extract = Join-Path $probeDir 'archive'
Expand-Archive -LiteralPath $zip -DestinationPath $extract -Force
$triangleSource = Get-ChildItem -LiteralPath $extract -Recurse -Filter triangle.exe -File | Select-Object -First 1
if (-not $triangleSource) { throw 'Windows StaMPS archive has no triangle.exe.' }
$runtimeDir = Join-Path $probeDir 'runtime'
New-Item -ItemType Directory -Path $runtimeDir -Force | Out-Null
Copy-Item -LiteralPath $triangleSource.FullName -Destination $runtimeDir
$gitCommand = (Get-Command git.exe).Source
$gitRoot = Split-Path (Split-Path $gitCommand -Parent) -Parent
$mingwBin = Join-Path $gitRoot 'mingw64\bin'
foreach ($dll in @('libstdc++-6.dll', 'libgcc_s_seh-1.dll', 'libwinpthread-1.dll')) {
    $source = Join-Path $mingwBin $dll
    if (-not (Test-Path -LiteralPath $source)) { throw "Missing Git for Windows runtime: $source" }
    Copy-Item -LiteralPath $source -Destination $runtimeDir
}
$node = Join-Path $probeDir 'phase-probe.1.node'
[System.IO.File]::WriteAllText($node, "3 2 0 0`n1 0 0`n2 1 0`n3 0 1`n")
$triangle = Join-Path $runtimeDir 'triangle.exe'
$stdout = Join-Path $probeDir 'triangle.out'
$stderr = Join-Path $probeDir 'triangle.err'
$probe = Start-Process -FilePath $triangle -ArgumentList '-e', 'phase-probe.1.node' `
    -WorkingDirectory $probeDir -Wait -PassThru -NoNewWindow `
    -RedirectStandardOutput $stdout -RedirectStandardError $stderr
$edge = Join-Path $probeDir 'phase-probe.2.edge'
if ($probe.ExitCode -ne 0 -or -not (Test-Path -LiteralPath $edge)) {
    Get-Content -LiteralPath $stdout -ErrorAction SilentlyContinue
    Get-Content -LiteralPath $stderr -ErrorAction SilentlyContinue
    throw "Windows Triangle functional probe failed (exit $($probe.ExitCode))."
}
Write-Output 'Windows Triangle created a valid edge file.'
