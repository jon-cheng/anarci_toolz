import pandas as pd
from Bio.Seq import Seq

from anarci_toolz.pipeline import file_or_dir_input, validate_aa_translation

# TODO: pipeline.validate_aa_translation is nearly identical to
# validation.py's copy of the same function, but only this pipeline.py copy
# is actually called by run_anarci_toolz — the validation.py copy appears to
# be dead code. Not resolved here per CLAUDE.md's known-issues list.


class TestFileOrDirInput:
    def test_directory_returns_only_csv_files(self, tmp_path):
        (tmp_path / "a.csv").write_text("a")
        (tmp_path / "b.csv").write_text("b")
        (tmp_path / "c.txt").write_text("c")

        result = file_or_dir_input([str(tmp_path)])

        assert sorted(result) == sorted(
            [str(tmp_path / "a.csv"), str(tmp_path / "b.csv")]
        )

    def test_direct_file_path_passed_through(self, tmp_path):
        file_path = tmp_path / "data.csv"
        file_path.write_text("data")

        result = file_or_dir_input([str(file_path)])

        assert result == [str(file_path)]

    def test_mix_of_dir_and_file(self, tmp_path):
        sub_dir = tmp_path / "sub"
        sub_dir.mkdir()
        (sub_dir / "a.csv").write_text("a")
        direct_file = tmp_path / "direct.csv"
        direct_file.write_text("direct")

        result = file_or_dir_input([str(sub_dir), str(direct_file)])

        assert sorted(result) == sorted([str(sub_dir / "a.csv"), str(direct_file)])


class TestValidateAaTranslation:
    def test_correct_translation_is_true(self):
        aa_seq = "EVQ"
        dna_seq = "GAAGTTCAG"  # translates to E V Q
        assert str(Seq(dna_seq).translate()) == aa_seq

        df = pd.DataFrame({"sequence_aa": [aa_seq], "sequence_dna": [dna_seq]})
        result = validate_aa_translation(df, "sequence_aa", "sequence_dna")

        assert result["is_correct_translation"].tolist() == [True]
        assert f"sequence_aa_translation" not in result.columns

    def test_incorrect_translation_is_false(self):
        aa_seq = "EVQ"
        dna_seq = "GGGGGGGGG"  # translates to G G G, not EVQ

        df = pd.DataFrame({"sequence_aa": [aa_seq], "sequence_dna": [dna_seq]})
        result = validate_aa_translation(df, "sequence_aa", "sequence_dna")

        assert result["is_correct_translation"].tolist() == [False]
