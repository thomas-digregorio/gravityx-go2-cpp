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
