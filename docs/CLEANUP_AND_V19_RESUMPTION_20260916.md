# Superseded output cleanup and V19 resumption

The user authorized deleting unnecessary old files in this local repository
and resuming the cold 19,402-bus improvement study. No source cases, code,
Git history, environment, configuration, detailed contingency JSON, worker
logs, status, official detail records, or verification certificates were
deleted. Only `solution_*.txt` files inside the five exact run directories
below were removed with the existing PowerShell payload-pruning helper.

Campaign: `runs/c2fen19402_6f0afaa_20260915_postcleanup`.

| Run directory | Removed file entries, including hard links | Reason |
|---|---:|---|
| C2FEN19402_s006_cold | 12,682 | Incomplete deadline-limited run; later complete run retained |
| C2FEN19402_s010_cold | 13,410 | Incomplete deadline-limited run; later complete run retained |
| C2FEN19402_s069_cold | 13,288 | Superseded by higher-scoring complete run |
| C2FEN19402_s077_cold | 13,216 | Superseded by higher-scoring complete run |
| C2FEN19402_s095_cold | 13,206 | Superseded by complete V16 |

Total: 65,802 file entries. The deduplicated candidate data size was
55.16864911932498 GiB. C: free space measured immediately before/after the
operation changed from 31,934,111,744 to 90,974,695,424 bytes: a net gain of
54.986 GiB, leaving 84.727 GiB free. This is a measured net volume change,
not a claim that summed hard-link lengths are recoverable storage.

Before deletion, every resolved target and ancestor was checked for reparse
points, every target was constrained to the named local campaign, no related
controller/native worker was running, and each status hash matched the
pre-existing compact archive. The archive SHA256 was
`8508e4f3f72585ef67ddb4fcb53209d76ce3db6f172c17d3d843c2addb8e4bcb`:
`docs/evidence/C2FEN19402_6F0AFAA_20260915_POSTCLEANUP.json`.
The three completed runs also passed the helper's complete-label-set,
official-zero-infeasibility and detail-count gates; the two incomplete runs
retained timeout status and worker-log evidence. After deletion, the status
hashes still matched, completed certificates remained, and no solution-text
entry remained in any of the five targets.

Complete newer solution-text sets remain protected:

| Scenario | Retained directory below runs/ | Root solution files including base |
|---|---|---:|
| 006 | reliability_20260915/target_19402_s006_v6/C2FEN19402_s006_cold | 6,694 |
| 010 | reliability_20260915/frozen_a1bebd5_19402_s010/C2FEN19402_s010_cold | 6,694 |
| 069 | reliability_20260915/frozen_a1bebd5_19402_s069/C2FEN19402_s069_cold | 6,621 |
| 077 | reliability_20260915/frozen_a1bebd5_19402_s077/C2FEN19402_s077_cold | 6,585 |
| 095 | top5_cold_20260915/v16_s095/C2FEN19402_s095_cold | 6,580 |

These counts and certificates were checked before and after cleanup. The
deleted text outputs cannot be independently re-evaluated without recreating
them; retained status/certificate records are verification history, not a
backup of those exact solutions. V18 and every other campaign were untouched.

## Next experiment

The already implemented V19 strict early-rejection change from `55cc0e8`
passed all four native CTest groups and all 89 Python tests again. Binary
SHA256: `e05958fdd35dc077defd083ab39f67f1eb923cda96ac2d25bacd044eaf84a8c1`.
Configuration SHA256:
`292f39175a66005cd9099a76690f43a5de253ac7f914b00cd12b1eafd1c61c4e`.

Run exactly once, cold, on scenario 095 with
`config/c2fen19402_top5_strict_trial_rejection_v19.json`, into the previously
nonexistent `runs/top5_cold_20260915/v19_s095`. Retain the unchanged 300-second
end-to-end deadline, 30 GiB pre-run floor, source constraints, tolerances,
full contingency set and official evaluation. No prior solution initializes
this run. This documentation/archive milestone changes no solver code or
configuration. No result is claimed until the experiment is complete.

## Superseded V19 payload cleanup after V20

V19 completed all native tasks but failed during worker shutdown before a
complete official evaluation. V20 later retained a complete independently
evaluated candidate set (audited after its own finalization timeout). V19's
solution-text payloads were therefore pruned after checking its compact
archive/status hashes, exact target/ancestor/subtree paths, no active related
processes, timeout evidence, and all 6,580 retained V20 root solutions.

Exact target:
`runs/top5_cold_20260915/v19_s095/C2FEN19402_s095_cold`.
Only 13,206 `solution_*.txt` entries were deleted. V19's status, worker logs,
partial official detail files, attribution and committed compact archive
remain. No V19 solution-text file remains; those exact text outputs cannot
be re-evaluated without recreation. V16 and V20 vectors remain intact.

C: free bytes changed from 56,074,645,504 to 67,826,470,912 during this cleanup:
a net gain of 10.945 GiB, leaving 63.168 GiB free. Across the two authorized
cleanup operations, 79,008 file entries were removed and measured net gains
totaled approximately 65.93 GiB. The new experiments and their working data
consumed some space in between; total reclaimed bytes should not be confused
with the change from the beginning-of-turn free-space reading.

## Superseded V16 payload cleanup after accepted V21

V21 now replaces V16 as the accepted on-time scenario-095 study incumbent:
375925.97863455745 versus 356907.86098102527, with full source-label coverage,
zero official infeasibility and unchanged residual/deadline gates. V21's
6,580 root vectors remain intact, as does the V20 candidate set.

After verifying the V16/V21 archive and status hashes, exact local target,
ancestor/subtree reparse checks and absence of active related workers, only
`runs/top5_cold_20260915/v16_s095/C2FEN19402_s095_cold/**/solution_*.txt`
was pruned. All 13,206 text-file entries were removed; no matching file remains.
Their summed apparent size was 21.950 GiB, which double-counts hard links and
is NOT a physical-space recovery measurement. The subsequent C: free-space
reading was 82,313,555,968 bytes (about 76.66 GiB).

V16 certificates, logs, evaluator details and committed compact audit survive;
they cannot recreate the deleted exact solution vectors. Source data, source
code, Git history, unrelated files, V18/V20/V21 vectors and other scenarios'
retained solutions were not removed. This supersedes earlier V16 retention
statements in this report.

## Superseded V21 payload cleanup after accepted V23

V23 passed all gates in 298.0589177999973 seconds with objective
379236.3363775767, replacing V21's accepted objective 375925.97863455745.
Its full raw-label set, source/model/config/binary/revision hashes, independent
residuals and official score were rechecked; all 6,580 root vectors remain.

After verifying the V21/V23 compact archive and status hashes, no active
related controller/native processes, exact contained paths, and no ancestor
or subtree reparse points, only
`runs/top5_cold_20260915/v21_s095/C2FEN19402_s095_cold/**/solution_*.txt`
was pruned. The operation removed 13,206 entries; zero matching text files
remain. C: free bytes changed from 57,207,144,448 to 68,958,728,192, a measured
net recovery of **10.9445152283 GiB**, leaving about **64.22 GiB free**.
The 21.953-GiB sum of removed file lengths includes hard-linked duplicates and
is not the physical recovery figure.

V21 logs, official detail records, certificates and committed audits remain.
Its deleted exact outputs cannot be independently re-evaluated without
recreating them. V20/V22 candidate vectors, V23's new accepted set, every
other scenario's retained solutions, raw sources, code and Git history were
untouched. This supersedes earlier V21 full-vector retention statements.

## Superseded V23 payload cleanup after accepted V25

V25 passed in 295.47726110000076 s with objective 381248.30942716985,
improving on V23's 379236.3363775767. Its independent objective audit,
exact native source-label set, complete official coverage, zero official
infeasibility and residual/deadline gates passed. All 6,580 V25 root vectors
are retained at `runs/top5_cold_20260915/v25_s095/C2FEN19402_s095_cold`.
The result and compact provenance were committed and pushed before pruning.

After verifying old/new archive and status hashes, absence of active relevant
Python/native workers, the exact contained target, and no ancestor/subtree
reparse points, only
`runs/top5_cold_20260915/v23_s095/C2FEN19402_s095_cold/**/solution_*.txt`
was deleted with the existing guarded PowerShell helper. All 13,206 entries
were removed; zero matching old payloads remain. The old status hash is
unchanged and its certificate, worker logs, official details and audit remain.
Code, Git history, source inputs, unrelated files and all V25 vectors remain.

C: free bytes changed from 43,280,080,896 to 55,023,595,520: a measured net
recovery of **10.9370002747 GiB**, leaving **51.2447166443 GiB free**. The
21.953-GiB sum of removed file lengths double-counts hard links and is not a
physical-space recovery figure. Across the cleanup operations recorded here,
118,626 solution-text directory entries have been removed, including links.

The deleted exact V23 solution texts cannot be recovered from retained
certificates or objective audits; those are verification records, not vector
backups. V20/V22/V24 candidate vectors and every other scenario's protected
solutions remain untouched. This supersedes earlier V23 full-vector retention
statements; V25 is now the protected accepted scenario-095 set.

## Failed V26 work-deadline payload cleanup

V26 failed in Code2 after 295.2568164000004 seconds, with 6,428/6,579
contingencies on time. Its compact timeout archive was checked against the
unchanged status hash and pushed in `7095c57` before cleanup. It has no final
official certificate or accepted score. V25's successful certificate and
all 6,580 root solution vectors were checked and protected.

After exact local containment, ancestor/subtree reparse checks, and absence
of relevant Windows/native workers, the existing failed-Code2-timeout helper
removed only
`runs/top5_cold_20260915/v26_s095/C2FEN19402_s095_cold/**/solution_*.txt`.
All 12,901 text entries were removed; zero remain. The status hash stayed
unchanged. C: free bytes rose from 42,515,759,104 to 54,015,488,000, a net
recovery of **10.7099571228 GiB**, leaving about **50.31 GiB free**. The
21.445-GiB sum of removed lengths double-counts hard links.

V26 logs, status, partial official detail records and committed archive
remain. Its deleted exact partial vectors cannot be recovered from these
records. V25's 6,580 vectors, other retained scenario outputs, raw cases,
source code, Git history and unrelated files were not deleted.

## Superseded V25 payload cleanup after accepted V27/095

The independently audited V27/095 result improved the accepted objective to
382016.64361338085 in 295.39416020000135 seconds. Its complete certificate,
source hashes, native checks and all 6,580 root vectors were verified, and
its compact archive/audit were pushed in `c6bcb2a` before pruning V25.

After checking old/new archive and status hashes, exact local containment,
all ancestor/subtree reparse points and absence of active Windows/native
workers, only
`runs/top5_cold_20260915/v25_s095/C2FEN19402_s095_cold/**/solution_*.txt`
was removed. The helper deleted 13,206 entries; zero remain. The V25
certificate, logs and audits remain, as do all 6,580 new V27/095 root vectors.
C: free bytes increased from 41,211,498,496 to 52,962,840,576: a net
**10.9442901611 GiB** recovered. Apparent lengths of 21.953 GiB include hard
links and are not physical recovered bytes.

The deleted exact V25 vectors cannot be recovered from its certificates.
This supersedes earlier V25 full-vector-retention statements. Source data,
code, Git history, unrelated files, other scenarios' complete reference
solutions and the newer V27/095 full vectors remain protected.

## V27 failed 006/010/069 payload cleanup

All three timed out in Code2. Their individual zero-diagnostic provenance
archives, hashes, exact coverage and failure report were pushed in `81c3740`
before deleting any payloads. The previous verified full references were
checked: 006 retained 6,694 root vectors, 010 retained 6,694, and 069 retained
6,621, with complete official certificates and zero reported infeasibility.
V27/077 and V27/095 also remained intact with 6,585 and 6,580 root vectors.

The same exact-path, containment, no-OneDrive, ancestor/subtree no-reparse,
archive/status hash, and no-active-worker checks passed before applying the
existing failed-Code2-timeout helper. Only solution text files were deleted:

| Run below `runs/top5_cold_20260915/` | Removed entries | Apparent GiB, includes hard links |
|---|---:|---:|
| v27_s006/C2FEN19402_s006_cold | 1919 | 3.291 |
| v27_s010/C2FEN19402_s010_cold | 3160 | 5.420 |
| v27_s069/C2FEN19402_s069_cold | 4584 | 7.625 |

Zero matching payloads remain in these three targets. Their original status
hashes are unchanged. C: free bytes rose from 31,981,330,432 to
40,858,066,944: a net **8.2671051025 GiB** recovered. Logs, failed-state JSON,
partial official records and committed timeout archives remain. These runs
never had complete official certificates. Deleted partial vectors cannot be
recovered from the retained records.

## Superseded original 077 reference payload cleanup

V27/077 passed in 177.62342740000167 seconds with objective
373417.7496024141 and complete independent verification. The superseded
campaign reference was archived separately before pruning:
`docs/evidence/C2FEN19402_BASELINE_S077_PREPRUNE_20260916.json`, SHA256
`3d814ca8b4c1770b25799e1b4260f360ccd91811205d87b02b49c99912624f7a`.
Archive diagnostics were zero, and it was pushed in `6aba067`. Its status
hash is `0838946948483b11bdefa1d00daf2d69c0c613c9fb041c15be6265a7de028db0`.
The replacement V27/077 archive/status hashes, full certificate and 6,585
root vectors were checked again.

After the same contained-path, no-reparse and no-active-worker checks, only
`runs/reliability_20260915/frozen_a1bebd5_19402_s077/C2FEN19402_s077_cold/**/solution_*.txt`
was pruned. It removed 13,216 entries; zero matching old payloads remain.
The old complete certificate, logs, audit/provenance and official details
remain, but they cannot recreate its deleted exact vectors. V27/077 retains
all 6,585 full root vectors and is now the protected 077 campaign solution.

C: free bytes rose from 40,849,612,800 to 52,636,852,224, a net recovery of
**10.9777221680 GiB**. Apparent deleted lengths of 21.998 GiB include links.
Across the V26, V25, three failed V27 and old-077 cleanup operations in this
iteration, net recovery totaled **40.8990745544 GiB**, with 48,986 text-file
directory entries removed. About **49.0 GiB** was free after retaining the
new experiments; recovery total is not the net free-space change across
intervening runs. Code, Git history, sources, environments and unrelated
files remain untouched. Earlier original-077 full-retention statements are
superseded; 006/010/069 retain their earlier complete reference vectors.
