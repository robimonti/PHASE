# Functional regression of production installer helpers, without launching WPF.
$ErrorActionPreference = 'Stop'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
Add-Type -AssemblyName System.Windows.Forms
$installer = Join-Path $PSScriptRoot '..\installer\install-phase.ps1'
$tokens = $null; $errors = $null
$ast = [System.Management.Automation.Language.Parser]::ParseFile(
    (Resolve-Path $installer).Path, [ref]$tokens, [ref]$errors)
if ($errors.Count -gt 0) { throw ($errors | Out-String) }
$names = @('Find-Git', 'Get-RemoteFile', 'Install-PortableGit',
           'Install-TriangleRuntime', 'Invoke-StampsBinariesDownload')
$definitions = $ast.FindAll({ param($node)
    $node -is [System.Management.Automation.Language.FunctionDefinitionAst] -and
    $names -contains $node.Name
}, $true)
if ($definitions.Count -ne $names.Count) { throw 'Production installer helpers not found.' }
foreach ($definition in $definitions) { Invoke-Expression $definition.Extent.Text }
$urlAssignment = $ast.Find({ param($node)
    $node -is [System.Management.Automation.Language.AssignmentStatementAst] -and
    $node.Left.Extent.Text -eq '$Script:GitPortableUrl'
}, $true)
if (-not $urlAssignment) { throw 'Pinned Portable Git URL not found.' }
Invoke-Expression $urlAssignment.Extent.Text

$probeDir = Join-Path $env:RUNNER_TEMP 'phase-triangle-smoke'
New-Item -ItemType Directory -Path $probeDir -Force | Out-Null
$savedLocalAppData = $env:LOCALAPPDATA
$env:LOCALAPPDATA = Join-Path $probeDir ([guid]::NewGuid().ToString('N'))
$gitRoot = Join-Path $env:LOCALAPPDATA 'PHASE\portable-git'
$stamps = Join-Path $probeDir ([guid]::NewGuid().ToString('N'))
$runtimeDir = Join-Path $stamps 'external\triangle\bin'
$status = { param($message) Write-Host $message }
try {
    $systemGit = Find-Git
    if (-not $systemGit -or $systemGit.StartsWith($env:LOCALAPPDATA)) {
        throw 'Regression fixture requires an existing system Git and no PHASE private cache.'
    }
    if (Test-Path -LiteralPath $gitRoot) { throw 'Expected an empty private runtime cache.' }
    Write-Host "Regression: system Git already exists at $systemGit; PHASE private runtime absent."
    if (-not (Invoke-StampsBinariesDownload -StampsRoot $stamps -StatusCallback $status)) {
        throw 'Fresh runtime preparation failed with existing system Git.'
    }
    $mingwBin = Join-Path $gitRoot 'mingw64\bin'
    $dlls = @('libstdc++-6.dll', 'libgcc_s_seh-1.dll', 'libwinpthread-1.dll')
    foreach ($dll in $dlls) {
        $a = (Get-FileHash -Algorithm SHA256 (Join-Path $mingwBin $dll)).Hash
        $b = (Get-FileHash -Algorithm SHA256 (Join-Path $runtimeDir $dll)).Hash
        if ($a -ne $b) { throw "Unexpected installed Triangle DLL: $dll" }
    }

    # Reproduce an incomplete cache with git.exe still present and an incorrect
    # DLL retained beside triangle.exe by an earlier/manual repair attempt.
    Remove-Item -LiteralPath (Join-Path $mingwBin 'libgcc_s_seh-1.dll') -Force
    [System.IO.File]::WriteAllText((Join-Path $runtimeDir 'libstdc++-6.dll'), 'invalid copied DLL')
    if (-not (Test-Path -LiteralPath (Join-Path $gitRoot 'bin\git.exe'))) {
        throw 'Partial-cache fixture lost git.exe.'
    }
    if (-not (Invoke-StampsBinariesDownload -StampsRoot $stamps -StatusCallback $status)) {
        throw 'Repair of incomplete cache or manually copied DLL failed.'
    }
    foreach ($dll in $dlls) {
        $a = (Get-FileHash -Algorithm SHA256 (Join-Path $mingwBin $dll)).Hash
        $b = (Get-FileHash -Algorithm SHA256 (Join-Path $runtimeDir $dll)).Hash
        if ($a -ne $b) { throw "Runtime repair left a mixed DLL set: $dll" }
    }

    # A complete cache must work without another Portable Git download.
    $Script:GitPortableUrl = 'https://invalid.example/should-not-be-requested.exe'
    if (-not (Invoke-StampsBinariesDownload -StampsRoot $stamps -StatusCallback $status)) {
        throw 'Complete runtime reuse unexpectedly required a download.'
    }
    if ((Find-Git) -ne $systemGit) { throw 'Installer changed system Git selection.' }
    Write-Output 'Production Windows installer passed: system Git, absent/partial/reused cache, DLL repair and real Triangle edge-file probes.'
} finally {
    $env:LOCALAPPDATA = $savedLocalAppData
}
