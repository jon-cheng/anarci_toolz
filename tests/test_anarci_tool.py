import pytest

from anarci_toolz.anarci_tool import ANARCI_RESULTS_COLS, extract_useful_info, get_anarci_alignment
import anarci_toolz.anarci_tool as anarci_tool_module


def _make_alignment_details():
    return [
        [
            {
                "query_start": 0,
                "query_end": 120,
                "chain_type": "H",
                "germlines": {
                    "v_gene": [("human", "IGHV3-66*01"), ("mouse", "other")],
                    "j_gene": [("human", "IGHJ4*01")],
                },
                "evalue": 1e-50,
                "bitscore": 200.5,
                "bias": 0.1,
            }
        ]
    ]


class TestExtractUsefulInfo:
    def test_unpacks_expected_fields(self):
        alignment_details = _make_alignment_details()
        result = extract_useful_info(alignment_details)

        assert result == (0, 120, "H", "IGHV3-66*01", "IGHJ4*01", 1e-50, 200.5, 0.1)


class TestGetAnarciAlignment:
    def test_type_error_fallback_returns_correct_shape(self, monkeypatch):
        def raise_type_error(*args, **kwargs):
            raise TypeError("boom")

        monkeypatch.setattr(anarci_tool_module, "get_alignment_details", raise_type_error)

        result = get_anarci_alignment(
            "EVQLVESGGGLVQPGG", "seq1", "imgt", allowed_species=["human"]
        )

        assert len(result) == len(ANARCI_RESULTS_COLS) + 1
        assert result[0] == "seq1"
        assert result[1] is False
        assert all(v is None for v in result[2:])

    def test_success_path_returns_passed_true(self, monkeypatch):
        monkeypatch.setattr(
            anarci_tool_module,
            "get_alignment_details",
            lambda seq, name, scheme, allowed_species: _make_alignment_details(),
        )

        result = get_anarci_alignment(
            "EVQLVESGGGLVQPGG", "seq1", "imgt", allowed_species=["human"]
        )

        assert result[0] == "seq1"
        assert result[1] is True
        assert result[2:] == (0, 120, "H", "IGHV3-66*01", "IGHJ4*01", 1e-50, 200.5, 0.1)
