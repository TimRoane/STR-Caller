# STR Genotype Caller

A compact STR genotyping algorithm that separates true allele signal from PCR stutter in read-count data. NGS STRs usually suffer from PCR stutter causing repeat counts above and below the true peak to be predictably boosted artificially. To reliably call STR genotypes we need to correct for this stutter to see the true peaks with a method of confidence of the call, that is what this caller does.

## What the algorithm does
Models expected STR stutter from a parent allele using an asymmetric polymerase-slippage model, then evaluates residual signal for evidence of a second allele.

For each locus, the caller:

1. collapses reads into repeat-length peaks;
2. selects the dominant peak as the primary allele;
3. predicts asymmetric PCR stutter from that peak;
4. subtracts the expected stutter contribution from neighboring peaks;
5. evaluates residual peaks using distance-aware evidence thresholds;
6. accepts the strongest supported secondary allele, if present;
7. removes stutter attributable to the secondary allele;
8. reports the genotype together with a confidence score and a peak-level audit table.

The stutter model is directional. Contraction stutter decays from 40% per repeat, while expansion stutter decays from 10% per repeat. Secondary-allele thresholds are stricter near the dominant peak, where artifact signal is expected to be strongest.

## Why it is structured this way

The calling core is independent of workflow engines and file-system conventions. `call_locus()` is the main analytical entry point, while the command-line wrapper only handles tabular I/O. This makes the algorithm easy to test, inspect, and reuse inside a larger workflow.

## Input

Read-count TSV:

```text
locus  repeat  forward_reads  reverse_reads
L1     10      620            380
L1      9      250            150
```

Allele-map TSV:

```text
locus  repeat  allele
L1     10      A
L1     11      B
```

## Output

The caller produces a locus-level genotype table and a peak-level audit table containing raw and stutter-corrected read counts. The audit output is intended to make every genotype decision traceable to the underlying peak evidence.

## Tests

The included tests cover the directional stutter model, residual-signal calculation, homozygous calls, heterozygous calls, low-depth rejection, and low-confidence rejection.
