import re
from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field, field_validator

_AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
_MIN_SEQ_LEN = 10
_MAX_SEQ_LEN = 500


class Scheme(str, Enum):
    # TODO(jon): confirm which schemes are actually wired up in
    # pipeline.py/anarci_tool.py/abnumber_tool.py — README's Arguments table
    # lists imgt/chothia/kabat/aho but Jon believes only kabat+imgt are
    # implemented. Check before filling in.
    imgt = "imgt"
    kabat = "kabat"


class Species(str, Enum):
    human = "human"
    mouse = "mouse"
    rat = "rat"
    rabbit = "rabbit"
    rhesus = "rhesus"
    pig = "pig"
    alpaca = "alpaca"


class NumberRequest(BaseModel):
    sequences: list[str] = Field(..., min_length=1)
    scheme: Scheme
    allowed_species: Optional[list[Species]] = None
    # Parallel, index-aligned list to `sequences` — the JSON-API-native
    # replacement for the CLI's seq_dna_header column-name concept, since
    # there is no DataFrame here.
    seq_dna: Optional[list[str]] = None
    skip_run_base_anarci: bool = False
    retain_indices: bool = False

    @field_validator("sequences")
    @classmethod
    def validate_sequences(cls, sequences: list[str]) -> list[str]:
        for i, seq in enumerate(sequences):
            if not (_MIN_SEQ_LEN <= len(seq) <= _MAX_SEQ_LEN):
                raise ValueError(
                    f"sequences[{i}] has length {len(seq)}, which is outside "
                    f"the allowed range of {_MIN_SEQ_LEN}-{_MAX_SEQ_LEN} amino acids"
                )
            if not _AA_PATTERN.match(seq):
                raise ValueError(
                    f"sequences[{i}] contains characters outside the standard "
                    "20 amino acids (expected only ACDEFGHIKLMNPQRSTVWY)"
                )
        return sequences

    @field_validator("seq_dna")
    @classmethod
    def validate_seq_dna_length(
        cls, seq_dna: Optional[list[str]], info
    ) -> Optional[list[str]]:
        if seq_dna is None:
            return seq_dna
        sequences = info.data.get("sequences")
        if sequences is not None and len(seq_dna) != len(sequences):
            raise ValueError(
                f"seq_dna has {len(seq_dna)} entries but sequences has "
                f"{len(sequences)} entries — seq_dna must be index-aligned "
                "with sequences"
            )
        return seq_dna


class NumberResponse(BaseModel):
    # TODO(jon): replace dict with a typed model once AIRR column mapping is
    # confirmed against run_anarci_toolz()'s actual output columns.
    results: list[dict]
