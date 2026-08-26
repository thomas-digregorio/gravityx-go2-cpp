[CmdletBinding()]
param(
    [switch] $Execute
)

$ErrorActionPreference = 'Stop'

$repositoryRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$runsRoot = (Resolve-Path -LiteralPath (Join-Path $repositoryRoot 'runs')).Path
$diagnosticsRoot = (Resolve-Path -LiteralPath (Join-Path $runsRoot 'diagnostics')).Path
$keptParentName = 'c2fen19402_b5f65c9_20260826'
$keptRunName = 'C2FEN19402_s006_cold_r18_measured_finalization'
$failedSiblingName = 'C2FEN19402_s010_cold_r1'
$keptParent = (Resolve-Path -LiteralPath (Join-Path $runsRoot $keptParentName)).Path
$keptRun = (Resolve-Path -LiteralPath (Join-Path $keptParent $keptRunName)).Path
$failedSibling = (Resolve-Path -LiteralPath (Join-Path $keptParent $failedSiblingName)).Path

if ($repositoryRoot.IndexOf('OneDrive', [StringComparison]::OrdinalIgnoreCase) -ge 0) {
    throw "OneDrive path is forbidden: $repositoryRoot"
}

$topLevelTargets = @(
    Get-ChildItem -LiteralPath $runsRoot -Directory -Filter 'c2fen19402_*' |
        Where-Object { $_.Name -ne $keptParentName }
)
$diagnosticTargets = @(
    Get-ChildItem -LiteralPath $diagnosticsRoot -Directory -Filter 'c2fen19402_*'
)

if ($topLevelTargets.Count -ne 38) {
    throw "Expected exactly 38 failed top-level run trees; found $($topLevelTargets.Count)."
}
if ($diagnosticTargets.Count -ne 220) {
    throw "Expected exactly 220 diagnostic trees; found $($diagnosticTargets.Count)."
}

function Assert-ChildTarget {
    param(
        [Parameter(Mandatory)] [string] $Root,
        [Parameter(Mandatory)] [string] $Target,
        [Parameter(Mandatory)] [string] $RequiredPrefix
    )

    $expectedRootPrefix = $Root + [IO.Path]::DirectorySeparatorChar
    if (-not $Target.StartsWith($expectedRootPrefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Target escaped its approved root: $Target"
    }
    $name = Split-Path -Leaf $Target
    if (-not $name.StartsWith($RequiredPrefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Target does not have the required prefix '$RequiredPrefix': $Target"
    }
}

foreach ($target in $topLevelTargets) {
    Assert-ChildTarget -Root $runsRoot -Target $target.FullName -RequiredPrefix 'c2fen19402_'
}
foreach ($target in $diagnosticTargets) {
    Assert-ChildTarget -Root $diagnosticsRoot -Target $target.FullName -RequiredPrefix 'c2fen19402_'
}
Assert-ChildTarget -Root $keptParent -Target $failedSibling -RequiredPrefix 'C2FEN19402_s010_'

$keptStatus = Join-Path $keptRun 'run_status.json'
$keptStatusObject = Get-Content -LiteralPath $keptStatus -Raw | ConvertFrom-Json
if ($keptStatusObject.success -ne $true -or $keptStatusObject.official_infeasibility -ne 0) {
    throw "The retained scenario-006 result no longer passes its success guard: $keptStatus"
}

$summary = [pscustomobject]@{
    RepositoryRoot = $repositoryRoot
    RetainedVerifiedRun = $keptRun
    FailedTopLevelTrees = $topLevelTargets.Count
    DiagnosticTrees = $diagnosticTargets.Count
    FailedSiblingTrees = 1
    Execute = [bool] $Execute
}
$summary | Format-List

if (-not $Execute) {
    Write-Host 'Dry run only. Re-run with -Execute to permanently remove these generated artifacts.'
    exit 0
}

foreach ($target in $diagnosticTargets) {
    Remove-Item -LiteralPath $target.FullName -Recurse -Force
}
foreach ($target in $topLevelTargets) {
    Remove-Item -LiteralPath $target.FullName -Recurse -Force
}
Remove-Item -LiteralPath $failedSibling -Recurse -Force

if (-not (Test-Path -LiteralPath $keptStatus)) {
    throw "The retained verified result was unexpectedly removed: $keptStatus"
}

$remainingTopLevel = @(
    Get-ChildItem -LiteralPath $runsRoot -Directory -Filter 'c2fen19402_*'
)
$remainingDiagnostics = @(
    Get-ChildItem -LiteralPath $diagnosticsRoot -Directory -Filter 'c2fen19402_*'
)
$drive = Get-PSDrive -Name ([IO.Path]::GetPathRoot($repositoryRoot).TrimEnd(':', '\'))

[pscustomobject]@{
    RetainedVerifiedRun = $keptRun
    RemainingC2FEN19402TopLevelTrees = $remainingTopLevel.Count
    RemainingC2FEN19402DiagnosticTrees = $remainingDiagnostics.Count
    FreeGB = [math]::Round($drive.Free / 1GB, 2)
} | Format-List
