# Debug plan: est_haps_var reproducibility at chr3R:19610000

Goal: Prove that production and test runs call the same function with the same inputs and produce the same outputs; if not, identify exactly what differs.

What we check (in order):
- Input equality (hashes): df4 data, founders, parameters (pos, sample, h_cutoff, method, chr)
- Code identity (hash): SHA-1 of body(est_haps_var)
- Core diagnostic: trace(Err) = sum of variances
- Full comparison: Groups, Names, aligned Err matrix diffs

Artifacts added (no production edits):
- scripts/ErrMatrix/est_haps_var_probe.R
  - Loads payload RDS (df4 + metadata), sources scripts/ErrMatrix/BASE_VAR_WIDE.R, runs est_haps_var
  - Prints SHA-1 digests for df4, founders, parameters; and function-body digest
  - Prints trace(Err) and Frobenius norm
  - Saves a __probe.txt next to the payload for record
- scripts/ErrMatrix/compare_single_position_against_prod.R
  - Loads production R.haps.chr3R.out.rds, extracts entry for chr3R:19610000, Rep01_W_F
  - Runs the probe on the extracted payload (position_data_chr3R_19610000_Rep01_W_F.RDS)
  - Prints probe report and compares production vs test (Groups/Names/Err)

Run commands on cluster:
- Rscript scripts/ErrMatrix/compare_single_position_against_prod.R

Interpretation:
- If input hashes and code hash match but Err differs: investigate numerical path (window chosen, A,y used for LSEI). Extend probe to log chosen window and A,y.
- If hashes differ: fix the pipeline so the exact same payload reaches the function.

Notes:
- In adaptive mode, h_cutoff MUST come from CLI/SLURM (third argument). Parameter files must not define h_cutoff.
- We avoid changing production logic; debug uses explicit dplyr::select at call-sites to avoid MASS masking.
- Keep trace(Err) as first-line diagnostic; it’s the signal we care about.
