param(
    [Parameter(Mandatory = $true)]
    [string]$CampaignRoot,

    [Parameter(Mandatory = $true)]
    [string]$RunDirectory,

    [switch]$AllowFailedTimedOutRun
)

$ErrorActionPreference = 'Stop'

$resolvedCampaign = [System.IO.Path]::GetFullPath($CampaignRoot).TrimEnd(
    [System.IO.Path]::DirectorySeparatorChar
)
$resolvedRun = [System.IO.Path]::GetFullPath($RunDirectory).TrimEnd(
    [System.IO.Path]::DirectorySeparatorChar
)
$requiredPrefix = $resolvedCampaign + [System.IO.Path]::DirectorySeparatorChar

if (-not $resolvedRun.StartsWith(
        $requiredPrefix,
        [System.StringComparison]::OrdinalIgnoreCase
    )) {
    throw "Refusing to prune outside campaign root: $resolvedRun"
}

$statusPath = Join-Path $resolvedRun 'run_status.json'
$certificatePath = Join-Path $resolvedRun `
    'internal\official_evaluation_certificate.json'
if (-not (Test-Path -LiteralPath $statusPath -PathType Leaf)) {
    throw "Missing run status: $statusPath"
}
$status = Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json
if ($status.success) {
    if (-not (Test-Path -LiteralPath $certificatePath -PathType Leaf)) {
        throw "Missing official evaluation certificate: $certificatePath"
    }
    $certificate = Get-Content -LiteralPath $certificatePath -Raw |
        ConvertFrom-Json
    if (-not $certificate.complete_label_set) {
        throw 'Refusing to prune an incomplete official label set'
    }
    if ([double]$certificate.reported_infeasibility -ne 0.0) {
        throw 'Refusing to prune a run with nonzero official infeasibility'
    }
    if (
        [int]$certificate.expected_detail_count -ne
        [int]$certificate.observed_unique_detail_count
    ) {
        throw 'Refusing to prune a run with incomplete official detail coverage'
    }
} elseif ($AllowFailedTimedOutRun) {
    $workerLogDirectory = Join-Path $resolvedRun `
        'internal\fast_screen_worker_logs'
    $workerLogs = @(
        Get-ChildItem -LiteralPath $workerLogDirectory -File `
            -Filter 'worker_*.log'
    )
    if (
        $status.stage -ne 'code2' -or
        -not $status.code2_timed_out -or
        [string]::IsNullOrWhiteSpace([string]$status.error) -or
        $workerLogs.Count -eq 0
    ) {
        throw 'Refusing failed-run pruning without preserved timeout evidence'
    }
} else {
    throw 'Refusing to prune a run without success=true'
}

$payloads = @(
    Get-ChildItem -LiteralPath $resolvedRun -File -Recurse |
        Where-Object { $_.Name -like 'solution_*.txt' }
)
$removedBytes = ($payloads | Measure-Object -Property Length -Sum).Sum
foreach ($payload in $payloads) {
    Remove-Item -LiteralPath $payload.FullName -Force
}

[pscustomobject]@{
    RunDirectory = $resolvedRun
    RemovedFiles = $payloads.Count
    RemovedGiB = [math]::Round($removedBytes / 1GB, 3)
    CertificateRetained = Test-Path -LiteralPath $certificatePath -PathType Leaf
}
