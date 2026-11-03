param(
    [ValidateSet('Plan','Execute','Verify')][string]$Mode='Plan'
)
$ErrorActionPreference='Stop'
$taskRepo='C:\Users\thoma\Documents\gravityx-go2-cpp'
$taskReceipts=Join-Path $taskRepo 'runs\cleanup_receipts'
$taskPlanPath=Join-Path $taskReceipts '20260916_paper_only_plan.json'
$taskReceiptPath=Join-Path $taskReceipts '20260916_paper_only_completed.json'
$taskRevision='3b335872647296c2d0e2c8ae5acbe4359ae72f0a'
$taskReportRevision='b5433e9d8a8fd4a22ab7df33e5c78198ad6a1f7d'
$taskKeepCampaign=Join-Path $taskRepo 'runs\reliability_v34_20260916'
$taskFallback=Join-Path $taskRepo 'runs\reliability_v33_20260916\19402_s077'

function Hash([string]$Path) {
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}
function Read-Json([string]$Path) {
    return Get-Content -LiteralPath $Path -Raw | ConvertFrom-Json
}
function Write-Json([string]$Path,$Value) {
    if (-not $Path.StartsWith($taskReceipts+'\',[StringComparison]::OrdinalIgnoreCase)) {throw 'Unsafe receipt path'}
    $taskText=$Value | ConvertTo-Json -Depth 14
    [IO.File]::WriteAllText($Path,$taskText+[Environment]::NewLine,[Text.UTF8Encoding]::new($false))
}
function Assert-SafePath([string]$Path) {
    $taskResolved=(Resolve-Path -LiteralPath $Path).ProviderPath.TrimEnd('\')
    if ($taskResolved -match '(?i)onedrive' -or -not $taskResolved.StartsWith($taskRepo+'\',[StringComparison]::OrdinalIgnoreCase)) {throw "Unsafe path: $taskResolved"}
    $taskNode=Get-Item -LiteralPath $taskResolved -Force
    while($null -ne $taskNode) {
        if (($taskNode.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {throw "Reparse ancestor: $($taskNode.FullName)"}
        if ($taskNode -is [IO.FileInfo]) {$taskNode=$taskNode.Directory} else {$taskNode=$taskNode.Parent}
    }
    return $taskResolved
}
function Assert-Idle {
    $taskActive=@(Get-CimInstance Win32_Process | Where-Object {
        $_.Name -match '^(gravityx.*|python.*|wsl|bash|julia|highs|ipopt|cmake|ctest)\.exe$' -and
        ($_.Name -match '^gravityx' -or $_.CommandLine -match 'gravityx-go2-cpp|run_experiment|run_reliability')
    })
    if ($taskActive.Count) {throw 'Relevant process active; stopping cleanup'}
    $taskWsl=((wsl.exe --list --running --quiet) -join '').Replace([string][char]0,'').Trim()
    if($taskWsl) {throw "WSL is running; stopping cleanup: $taskWsl"}
    if ((git -C $taskRepo rev-parse HEAD).Trim() -ne $taskReportRevision) {throw 'Repository revision changed'}
    git -C $taskRepo diff --quiet HEAD --
    if($LASTEXITCODE -ne 0) {throw 'Tracked worktree changed'}
}
function Get-TreeInfo([string]$Path) {
    $taskResolved=Assert-SafePath $Path
    $taskFirst=Get-Item -LiteralPath $taskResolved -Force
    [long]$taskFiles=0;[long]$taskBytes=0;[long]$taskDirectories=0;[long]$taskMaxWrite=0
    if ($taskFirst -is [IO.FileInfo]) {
        $taskFiles=1;$taskBytes=$taskFirst.Length;$taskMaxWrite=$taskFirst.LastWriteTimeUtc.Ticks
    } else {
        $taskStack=[Collections.Generic.Stack[IO.DirectoryInfo]]::new()
        $taskStack.Push($taskFirst)
        while($taskStack.Count) {
            $taskDirectory=$taskStack.Pop();$taskDirectories++
            foreach($taskEntry in $taskDirectory.EnumerateFileSystemInfos()) {
                if (($taskEntry.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {throw "Reparse descendant: $($taskEntry.FullName)"}
                if(-not $taskEntry.FullName.StartsWith($taskResolved+'\',[StringComparison]::OrdinalIgnoreCase)) {throw 'Tree escaped target'}
                if($taskEntry -is [IO.DirectoryInfo]) {$taskStack.Push($taskEntry)} else {
                    $taskFiles++;$taskBytes+=$taskEntry.Length
                    if($taskEntry.LastWriteTimeUtc.Ticks -gt $taskMaxWrite) {$taskMaxWrite=$taskEntry.LastWriteTimeUtc.Ticks}
                }
            }
        }
    }
    return [pscustomobject]@{path=$taskResolved;fileCount=$taskFiles;directoryCount=$taskDirectories;apparentBytes=$taskBytes;maximumFileWriteTicks=$taskMaxWrite}
}
function Assert-Target([string]$Path) {
    $taskResolved=Assert-SafePath $Path
    $taskParent=[IO.Path]::GetDirectoryName($taskResolved)
    $taskAllowedParents=@((Join-Path $taskRepo 'runs'),(Join-Path $taskRepo 'diagnostics'),(Join-Path $taskRepo 'runs\reliability_v33_20260916'))
    if($taskParent -notin $taskAllowedParents) {throw "Target is not an approved direct child: $taskResolved"}
    foreach($taskProtected in @($taskKeepCampaign,$taskReceipts,$taskFallback,(Join-Path $taskRepo 'diagnostics\v34_admin'),(Join-Path $taskRepo 'diagnostics\artifact_read'))) {
        if($taskResolved -eq $taskProtected -or $taskResolved.StartsWith($taskProtected+'\',[StringComparison]::OrdinalIgnoreCase) -or $taskProtected.StartsWith($taskResolved+'\',[StringComparison]::OrdinalIgnoreCase)) {throw "Protected overlap: $taskResolved"}
    }
    return $taskResolved
}
function Get-PaperAudit {
    $taskData=Read-Json (Join-Path $taskRepo 'paper\data\evidence.json')
    $taskQa=Read-Json (Join-Path $taskRepo 'paper\qa_report.json')
    if($taskData.algorithm_revision -ne $taskRevision -or $taskData.cases.Count -ne 37) {throw 'Unexpected paper identity'}
    $taskEvidenceHash=Hash (Join-Path $taskRepo 'paper\data\evidence.json')
    if($taskEvidenceHash -ne $taskQa.evidence_sha256) {throw 'Paper evidence hash mismatch'}
    if((Hash (Join-Path $taskRepo 'output\pdf\verified_ac_dispatch.pdf')) -ne $taskQa.pdf_sha256) {throw 'Paper PDF hash mismatch'}
    if((Hash (Join-Path $taskRepo 'output\pdf\verified_ac_dispatch_latex.zip')) -ne $taskQa.latex_zip_sha256) {throw 'Paper source bundle hash mismatch'}
    $taskSummaryPath=Join-Path $taskRepo 'docs\evidence\V34_ALL_SCENARIOS_SUMMARY_20260916.json'
    if((Hash $taskSummaryPath) -ne $taskData.summary_sha256) {throw 'Campaign summary hash mismatch'}
    $taskConfigPath=Join-Path $taskRepo 'config\v34_all_scenarios_20260916.json'
    if((Hash $taskConfigPath) -ne $taskData.config_sha256) {throw 'Frozen configuration mismatch'}
    if((Hash (Join-Path $taskRepo '.build-native\gravityx_go2')) -ne $taskData.executable_sha256) {throw 'Frozen executable mismatch'}
    $taskConfig=Read-Json $taskConfigPath
    $taskRecords=@()
    foreach($taskCase in $taskData.cases) {
        $taskRun=Assert-SafePath (Join-Path $taskRepo $taskCase.run_root)
        if(-not $taskRun.StartsWith($taskKeepCampaign+'\',[StringComparison]::OrdinalIgnoreCase)) {throw 'Paper references another campaign'}
        $taskStatusPath=Join-Path $taskRun 'run_status.json'
        $taskRunSummaryPath=Join-Path $taskRun 'run_summary.json'
        if((Hash $taskStatusPath) -ne $taskCase.status_sha256 -or (Hash $taskRunSummaryPath) -ne $taskCase.summary_sha256) {throw "Paper run hashes mismatch: $taskRun"}
        $taskStatus=Read-Json $taskStatusPath
        $taskRunSummary=Read-Json $taskRunSummaryPath
        $taskCert=Read-Json (Join-Path $taskRun 'internal\official_evaluation_certificate.json')
        if(-not $taskStatus.success -or $taskStatus.stage -ne 'complete' -or $taskStatus.git_revision -ne $taskRevision -or -not $taskCert.complete_label_set -or $taskCert.reported_infeasibility -ne 0 -or $taskCert.expected_detail_count -ne ($taskCase.contingencies+1) -or $taskCert.observed_unique_detail_count -ne ($taskCase.contingencies+1) -or $taskCert.reported_objective -ne $taskCase.objective) {throw "Paper acceptance mismatch: $taskRun"}
        if((Hash (Join-Path $taskRun $taskRunSummary.base.artifact)) -ne $taskRunSummary.base.sha256) {throw 'Paper base artifact mismatch'}
        $taskVectors=@(Get-ChildItem -LiteralPath $taskRun -File -Filter 'solution_*.txt')
        if($taskVectors.Count -ne $taskCase.full_root_vectors -or @($taskVectors | Where-Object {$_.Length -eq 0}).Count) {throw "Paper root vectors changed: $taskRun"}
        $taskFamily=$taskCase.family;$taskScenario=$taskCase.scenario
        $taskPrefix=if($taskFamily -eq '00617' -and $taskScenario -eq '005') {'617'} else {$taskFamily}
        $taskModel=if($taskFamily -in @('17700','19402')) {Join-Path $taskRepo ".data\final_$taskFamily\scenario_$taskScenario\model.json"} else {Join-Path $taskRepo ".data\final_$($taskPrefix)_$taskScenario\model.json"}
        if((Hash $taskModel) -ne $taskCase.input_hashes.normalized_case_sha256) {throw 'Normalized input mismatch'}
        $taskSource=Join-Path $taskConfig.source_root "C2FEN$taskFamily\scenario_$taskScenario"
        foreach($taskSourceFile in @(@('case.raw','raw_sha256'),@('case.json','json_sha256'),@('case.con','con_sha256'))) {
            if((Hash (Join-Path $taskSource $taskSourceFile[0])) -ne $taskCase.input_hashes.($taskSourceFile[1])) {throw 'Raw source input mismatch'}
        }
        $taskRecords+=[pscustomobject]@{run=$taskRun;statusSha256=$taskCase.status_sha256;summarySha256=$taskCase.summary_sha256;rootVectors=$taskVectors.Count;sourceContingencies=$taskCase.contingencies}
    }
    $taskFallbackRun=Join-Path $taskFallback 'C2FEN19402_s077_cold'
    $taskFallbackCert=Read-Json (Join-Path $taskFallbackRun 'internal\official_evaluation_certificate.json')
    $taskFallbackVectors=@(Get-ChildItem -LiteralPath $taskFallbackRun -File -Filter 'solution_*.txt').Count
    if(-not $taskFallbackCert.complete_label_set -or $taskFallbackCert.reported_infeasibility -ne 0 -or $taskFallbackVectors -ne 6585) {throw 'Fallback solution is incomplete'}
    return [pscustomobject]@{algorithmRevision=$taskRevision;paperEvidenceSha256=$taskEvidenceHash;sourceContingencies=($taskData.cases | Measure-Object -Property contingencies -Sum).Sum;runs=$taskRecords;fallbackRun=$taskFallbackRun;fallbackRootVectors=$taskFallbackVectors;missingPaperVectorCase='19402/077: pruned before this cleanup; its V34 records are preserved'}
}

Set-Location -LiteralPath $taskRepo
Assert-Idle
$taskTracked=@(git -C $taskRepo ls-files -- runs diagnostics)
if($taskTracked.Count) {throw 'Tracked files found in generated-data scope'}
if($Mode -eq 'Plan') {
    if(Test-Path -LiteralPath $taskPlanPath) {throw 'Plan already exists; refusing overwrite'}
    Write-Host 'Checking the paper, all 37 saved runs, source hashes, and frozen artifacts.'
    $taskAudit=Get-PaperAudit
    $taskTargets=@()
    foreach($taskItem in Get-ChildItem -LiteralPath (Join-Path $taskRepo 'runs') -Force) {
        if($taskItem.Name -in @('reliability_v34_20260916','cleanup_receipts')) {continue}
        if($taskItem.Name -eq 'reliability_v33_20260916') {
            foreach($taskSub in Get-ChildItem -LiteralPath $taskItem.FullName -Force) {if($taskSub.Name -ne '19402_s077') {$taskTargets+=$taskSub.FullName}}
        } else {$taskTargets+=$taskItem.FullName}
    }
    foreach($taskItem in Get-ChildItem -LiteralPath (Join-Path $taskRepo 'diagnostics') -Force) {if($taskItem.Name -notin @('v34_admin','artifact_read')) {$taskTargets+=$taskItem.FullName}}
    $taskInventory=@();$taskIndex=0
    foreach($taskTarget in $taskTargets) {
        $taskSafe=Assert-Target $taskTarget
        $taskInventory+=Get-TreeInfo $taskSafe
        $taskIndex++
        if($taskIndex%20 -eq 0) {Write-Host "Inventoried $taskIndex / $($taskTargets.Count) old generated targets."}
    }
    $taskPlan=[ordered]@{schemaVersion=1;createdUtc=[datetime]::UtcNow.ToString('o');repository=$taskRepo;operation='Remove obsolete generated run and diagnostic trees, retaining every paper V34 run and the V33/077 fallback pending explicit decision';protectedPaperAudit=$taskAudit;protectedCampaignInventory=(Get-TreeInfo $taskKeepCampaign);protectedFallbackInventory=(Get-TreeInfo $taskFallback);targets=$taskInventory;targetCount=$taskInventory.Count;fileCount=($taskInventory | Measure-Object -Property fileCount -Sum).Sum;apparentBytes=($taskInventory | Measure-Object -Property apparentBytes -Sum).Sum;freeBytesAtPlan=[IO.DriveInfo]::new('C:\').AvailableFreeSpace}
    Write-Json $taskPlanPath $taskPlan
    Write-Host "PLAN COMPLETE: $($taskPlan.targetCount) targets, $($taskPlan.fileCount) files; apparent GiB $([math]::Round($taskPlan.apparentBytes/1GB,2)) (hard links can double-count). No deletion performed."
    exit 0
}
$taskPlan=Read-Json $taskPlanPath
if($taskPlan.repository -ne $taskRepo -or $taskPlan.protectedPaperAudit.algorithmRevision -ne $taskRevision) {throw 'Plan identity mismatch'}
if($Mode -eq 'Execute') {
    if(Test-Path -LiteralPath $taskReceiptPath) {throw 'Execution receipt exists; refusing repeat'}
    $taskAudit=Get-PaperAudit
    if($taskAudit.paperEvidenceSha256 -ne $taskPlan.protectedPaperAudit.paperEvidenceSha256) {throw 'Paper changed since plan'}
    Write-Host 'Checking every exact deletion target against the saved plan.'
    $taskIndex=0
    foreach($taskTarget in $taskPlan.targets) {
        $taskSafe=Assert-Target $taskTarget.path
        $taskNow=Get-TreeInfo $taskSafe
        foreach($taskField in @('fileCount','directoryCount','apparentBytes','maximumFileWriteTicks')) {if($taskNow.$taskField -ne $taskTarget.$taskField) {throw "Target changed: $taskSafe ($taskField)"}}
        $taskIndex++
        if($taskIndex%25 -eq 0) {Write-Host "Rechecked $taskIndex / $($taskPlan.targetCount) targets."}
    }
    Assert-Idle
    $taskReceipt=[ordered]@{schemaVersion=1;startedUtc=[datetime]::UtcNow.ToString('o');planSha256=(Hash $taskPlanPath);freeBytesBefore=[IO.DriveInfo]::new('C:\').AvailableFreeSpace;completed=$false;removed=@();error=$null}
    Write-Json $taskReceiptPath $taskReceipt
    try {
        foreach($taskTarget in $taskPlan.targets) {
            Assert-Idle
            $taskSafe=Assert-Target $taskTarget.path
            if(Test-Path -LiteralPath (Join-Path $taskReceipts 'STOP_PAPER_CLEANUP')) {throw 'User stop marker detected'}
            Remove-Item -LiteralPath $taskSafe -Recurse -Force
            if(Test-Path -LiteralPath $taskSafe) {throw "Deletion incomplete: $taskSafe"}
            $taskReceipt.removed+=@{path=$taskSafe;fileCount=$taskTarget.fileCount;apparentBytes=$taskTarget.apparentBytes;finishedUtc=[datetime]::UtcNow.ToString('o')}
            $taskReceipt.freeBytesLatest=[IO.DriveInfo]::new('C:\').AvailableFreeSpace
            Write-Json $taskReceiptPath $taskReceipt
            if($taskReceipt.removed.Count%10 -eq 0 -or $taskTarget.apparentBytes -gt 5GB) {Write-Host "Removed $($taskReceipt.removed.Count) / $($taskPlan.targetCount) targets; C: free $([math]::Round($taskReceipt.freeBytesLatest/1GB,2)) GiB."}
        }
        $taskReceipt.completed=$true
        $taskReceipt.finishedUtc=[datetime]::UtcNow.ToString('o')
        $taskReceipt.freeBytesAfter=[IO.DriveInfo]::new('C:\').AvailableFreeSpace
        $taskReceipt.recoveredPhysicalBytes=$taskReceipt.freeBytesAfter-$taskReceipt.freeBytesBefore
        Write-Json $taskReceiptPath $taskReceipt
    } catch {
        $taskReceipt.error=$_.Exception.Message
        $taskReceipt.freeBytesLatest=[IO.DriveInfo]::new('C:\').AvailableFreeSpace
        Write-Json $taskReceiptPath $taskReceipt
        throw
    }
}
$taskReceipt=Read-Json $taskReceiptPath
if(-not $taskReceipt.completed -or $taskReceipt.removed.Count -ne $taskPlan.targetCount) {throw 'Cleanup did not fully complete'}
foreach($taskTarget in $taskPlan.targets) {if(Test-Path -LiteralPath $taskTarget.path) {throw "Removed target reappeared: $($taskTarget.path)"}}
Write-Host 'Rechecking the paper and all protected run/input hashes after deletion.'
$taskAudit=Get-PaperAudit
foreach($taskPair in @(@($taskKeepCampaign,$taskPlan.protectedCampaignInventory),@($taskFallback,$taskPlan.protectedFallbackInventory))) {
    $taskAfter=Get-TreeInfo $taskPair[0]
    foreach($taskField in @('fileCount','directoryCount','apparentBytes','maximumFileWriteTicks')) {if($taskAfter.$taskField -ne $taskPair[1].$taskField) {throw "Protected inventory changed: $($taskPair[0])"}}
}
$taskReceipt | Add-Member -NotePropertyName postcheck -NotePropertyValue @{passed=$true;checkedUtc=[datetime]::UtcNow.ToString('o');paperRunCount=37;sourceContingencies=$taskAudit.sourceContingencies;paperHashesUnchanged=$true;sourceHashesUnchanged=$true;protectedInventoriesUnchanged=$true;algorithmChanged=$false;fallbackKept=$true} -Force
Write-Json $taskReceiptPath $taskReceipt
Write-Host "COMPLETE: removed $($taskReceipt.removed.Count) obsolete targets; recovered $([math]::Round($taskReceipt.recoveredPhysicalBytes/1GB,3)) GiB; C: free $([math]::Round([IO.DriveInfo]::new('C:\').AvailableFreeSpace/1GB,3)) GiB. Paper and protected data verified."
