import pandas as pd
import pytest

from str_caller import (
    CallerParameters,
    call_locus,
    remove_stutter,
    stutter_fraction,
)


PARAMS = CallerParameters(
    min_allele_reads=20,
    min_locus_reads=100,
    min_confidence=70,
)
ALLELES = {repeat: str(repeat) for repeat in range(1, 30)}


def reads(repeats, counts):
    return pd.DataFrame(
        {
            "locus": ["L1"] * len(repeats),
            "repeat": repeats,
            "forward_reads": counts,
            "reverse_reads": [0] * len(repeats),
        }
    )


def test_directional_stutter_model():
    assert stutter_fraction(-1, PARAMS) == 0.40
    assert stutter_fraction(-2, PARAMS) == pytest.approx(0.16)
    assert stutter_fraction(1, PARAMS) == 0.10
    assert stutter_fraction(2, PARAMS) == pytest.approx(0.01)


def test_stutter_subtraction_clamps_at_zero():
    assert remove_stutter(400, 1000, -1, PARAMS) == 0
    assert remove_stutter(450, 1000, -1, PARAMS) == 50


def test_clean_homozygous_pattern_calls_homozygous():
    locus = reads(
        [10, 9, 8, 11, 12],
        [1000, 400, 160, 100, 10],
    )
    call, audit = call_locus("L1", locus, ALLELES, PARAMS)

    assert call.genotype == "10/10"
    assert call.status == "PASS"
    assert audit.loc[audit["repeat"] == 9, "corrected_reads"].iat[0] == 0
    assert audit.loc[audit["repeat"] == 11, "corrected_reads"].iat[0] == 0


def test_true_secondary_signal_survives_stutter_correction():
    locus = reads(
        [10, 9, 8, 11, 12],
        [1000, 400, 160, 100, 350],
    )
    call, audit = call_locus("L1", locus, ALLELES, PARAMS)

    assert call.genotype == "10/12"
    assert call.status == "PASS"
    assert call.secondary_repeat == 12
    assert set(audit.loc[audit["selected"], "repeat"]) == {10, 12}


def test_low_depth_returns_no_call():
    call, _ = call_locus("L1", reads([10], [40]), ALLELES, PARAMS)
    assert call.genotype == "NoCall"
    assert call.status == "LOW_DEPTH"


def test_low_confidence_returns_candidate_genotype_as_context():
    strict = CallerParameters(min_confidence=95)
    locus = reads([10, 12, 13], [500, 150, 73])
    call, _ = call_locus("L1", locus, ALLELES, strict)

    assert call.status == "LOW_CONFIDENCE"
    assert call.genotype.startswith("NoCall(")
