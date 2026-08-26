[CmdletBinding()]
param(
    [switch] $Execute
)

$ErrorActionPreference = 'Stop'

$repositoryRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$runsRoot = (Resolve-Path -LiteralPath (Join-Path $repositoryRoot 'runs')).Path
$diagnosticsRoot = (Resolve-Path -LiteralPath (Join-Path $runsRoot 'diagnostics')).Path

if ($repositoryRoot.IndexOf('OneDrive', [StringComparison]::OrdinalIgnoreCase) -ge 0) {
    throw "OneDrive path is forbidden: $repositoryRoot"
}

$retainedRunSpecs = @(
    [pscustomobject]@{
        ParentName = 'c2fen19402_b5f65c9_20260826'
        RunName = 'C2FEN19402_s006_cold_r18_measured_finalization'
    },
    [pscustomobject]@{
        ParentName = 'c2fen19402_52103f8_20260826'
        RunName = 's010_cold_r24_h11'
    },
    [pscustomobject]@{
        ParentName = 'c2fen19402_1dd8972_20260826'
        RunName = 's069_cold_r3_corrected'
    }
)

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

$retainedRuns = @()
$retainedParentNames = @()
foreach ($spec in $retainedRunSpecs) {
    $parent = (Resolve-Path -LiteralPath (Join-Path $runsRoot $spec.ParentName)).Path
    $run = (Resolve-Path -LiteralPath (Join-Path $parent $spec.RunName)).Path
    $statusPath = Join-Path $run 'run_status.json'
    $status = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
    if ($status.success -ne $true -or $status.official_infeasibility -ne 0) {
        throw "A retained result no longer passes its success guard: $statusPath"
    }
    $retainedRuns += [pscustomobject]@{
        Parent = $parent
        Run = $run
        StatusPath = $statusPath
        ParentName = $spec.ParentName
        RunName = $spec.RunName
    }
    $retainedParentNames += $spec.ParentName
}

$allTopLevel = @(
    Get-ChildItem -LiteralPath $runsRoot -Directory -Filter 'c2fen19402_*'
)

# Refuse to prune any newly discovered successful result. A future verified run
# must be registered explicitly before this maintenance script can remove data.
$unexpectedVerifiedRuns = @()
foreach ($root in $allTopLevel) {
    foreach ($child in @(Get-ChildItem -LiteralPath $root.FullName -Directory -ErrorAction SilentlyContinue)) {
        $statusPath = Join-Path $child.FullName 'run_status.json'
        if (-not (Test-Path -LiteralPath $statusPath)) {
            continue
        }
        try {
            $status = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
        } catch {
            continue
        }
        if ($status.success -eq $true -and $status.official_infeasibility -eq 0) {
            $isExpected = $false
            foreach ($retained in $retainedRuns) {
                if ($child.FullName.Equals($retained.Run, [StringComparison]::OrdinalIgnoreCase)) {
                    $isExpected = $true
                    break
                }
            }
            if (-not $isExpected) {
                $unexpectedVerifiedRuns += $child.FullName
            }
        }
    }
}
if ($unexpectedVerifiedRuns.Count -gt 0) {
    throw "Found unregistered verified result(s); add them to retainedRunSpecs before pruning: $($unexpectedVerifiedRuns -join ', ')"
}

$topLevelTargets = @(
    $allTopLevel | Where-Object { $_.Name -notin $retainedParentNames }
)
$failedSiblingTargets = @()
foreach ($retained in $retainedRuns) {
    $failedSiblingTargets += @(
        Get-ChildItem -LiteralPath $retained.Parent -Directory -ErrorAction SilentlyContinue |
            Where-Object { -not $_.FullName.Equals($retained.Run, [StringComparison]::OrdinalIgnoreCase) }
    )
}
$diagnosticTargets = @(
    Get-ChildItem -LiteralPath $diagnosticsRoot -Directory -Filter 'c2fen19402_*'
)

foreach ($target in $topLevelTargets) {
    Assert-ChildTarget -Root $runsRoot -Target $target.FullName -RequiredPrefix 'c2fen19402_'
}
foreach ($target in $diagnosticTargets) {
    Assert-ChildTarget -Root $diagnosticsRoot -Target $target.FullName -RequiredPrefix 'c2fen19402_'
}
foreach ($target in $failedSiblingTargets) {
    $matchingRetained = $retainedRuns | Where-Object {
        $target.FullName.StartsWith($_.Parent + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)
    } | Select-Object -First 1
    if ($null -eq $matchingRetained) {
        throw "Failed sibling escaped every retained parent: $($target.FullName)"
    }
}

[pscustomobject]@{
    RepositoryRoot = $repositoryRoot
    RetainedVerifiedRuns = ($retainedRuns.Run -join [Environment]::NewLine)
    FailedTopLevelTrees = $topLevelTargets.Count
    DiagnosticTrees = $diagnosticTargets.Count
    FailedSiblingTrees = $failedSiblingTargets.Count
    Execute = [bool] $Execute
} | Format-List

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
foreach ($target in $failedSiblingTargets) {
    Remove-Item -LiteralPath $target.FullName -Recurse -Force
}

foreach ($retained in $retainedRuns) {
    if (-not (Test-Path -LiteralPath $retained.StatusPath)) {
        throw "A retained verified result was unexpectedly removed: $($retained.StatusPath)"
    }
}

$remainingTopLevel = @(
    Get-ChildItem -LiteralPath $runsRoot -Directory -Filter 'c2fen19402_*'
)
$remainingDiagnostics = @(
    Get-ChildItem -LiteralPath $diagnosticsRoot -Directory -Filter 'c2fen19402_*'
)
$drive = Get-PSDrive -Name ([IO.Path]::GetPathRoot($repositoryRoot).TrimEnd(':', '\'))

[pscustomobject]@{
    RetainedVerifiedRuns = ($retainedRuns.Run -join [Environment]::NewLine)
    RemainingC2FEN19402TopLevelTrees = $remainingTopLevel.Count
    RemainingC2FEN19402DiagnosticTrees = $remainingDiagnostics.Count
    FreeGB = [math]::Round($drive.Free / 1GB, 2)
} | Format-List
