import pytest
from pydantic import ValidationError

from anarci_toolz.api.models import NumberRequest, Scheme, Species

VALID_SEQ = "QVQLVQSGAEVKKPGASVKVSCKASGYTFT"  # 30 aa, valid alphabet


def test_valid_request_passes():
    request = NumberRequest(
        sequences=[VALID_SEQ],
        scheme=Scheme.imgt,
        allowed_species=[Species.human],
    )
    assert request.sequences == [VALID_SEQ]
    assert request.scheme == Scheme.imgt


def test_invalid_character_raises():
    bad_seq = VALID_SEQ[:-1] + "X"  # X is not one of the 20 standard aa
    with pytest.raises(ValidationError):
        NumberRequest(sequences=[bad_seq], scheme=Scheme.kabat)


def test_sequence_too_short_raises():
    with pytest.raises(ValidationError):
        NumberRequest(sequences=["ACDEFG"], scheme=Scheme.imgt)  # 6 aa, below 10 min


def test_sequence_too_long_raises():
    with pytest.raises(ValidationError):
        NumberRequest(sequences=["A" * 501], scheme=Scheme.imgt)  # above 500 max


def test_mismatched_seq_dna_length_raises():
    with pytest.raises(ValidationError):
        NumberRequest(
            sequences=[VALID_SEQ, VALID_SEQ],
            seq_dna=["ATG"],
            scheme=Scheme.imgt,
        )
