"""Regression tests for parsing, validation, MAFFT calls, and FASTA output."""

import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from msa import MSAInputError, align_sequences, format_fasta, parse_sequences


class MSATest(unittest.TestCase):
    """Small examples that protect the app's public alignment helpers."""

    def test_parse_fasta(self):
        """FASTA names and sequences are preserved."""
        records = parse_sequences(">a\nACD\n>b\nACDE\n")

        self.assertEqual([record.name for record in records], ["a", "b"])
        self.assertEqual([record.sequence for record in records], ["ACD", "ACDE"])

    def test_parse_plain(self):
        """Plain input lines receive stable generated names."""
        records = parse_sequences("ACD\nACDE\n")

        self.assertEqual(records[0].name, "sequence_1")
        self.assertEqual(records[1].sequence, "ACDE")

    def test_rejects_invalid_characters(self):
        """Non-amino-acid symbols fail before alignment."""
        with self.assertRaises(MSAInputError):
            parse_sequences(">a\nACD1\n>b\nACDE\n")

    def test_sequence_length_limit(self):
        """Sequences may be 1000 aa long, but not longer."""
        records = parse_sequences(f">a\n{'A' * 1000}\n>b\nACD\n")

        self.assertEqual(len(records[0].sequence), 1000)
        with self.assertRaises(MSAInputError):
            parse_sequences(f">a\n{'A' * 1001}\n>b\nACD\n")

    def test_alignment_runs_mafft(self):
        """Alignment delegates to MAFFT and parses aligned FASTA."""
        def fake_mafft(command, **kwargs):
            """Confirm MAFFT receives FASTA input and return aligned output."""
            self.assertIn("--auto", command)
            self.assertIn("--amino", command)
            self.assertIn(">a\nACD", Path(command[-1]).read_text())
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=">a\nAC-D\n>b\nACED\n",
                stderr="",
            )

        records = parse_sequences(">a\nACD\n>b\nACED\n")
        with patch("msa.shutil.which", return_value="/usr/bin/mafft"):
            with patch("msa.subprocess.run", side_effect=fake_mafft):
                alignment = align_sequences(records)

        self.assertEqual(alignment, [("a", "AC-D"), ("b", "ACED")])

    def test_missing_mafft_is_reported(self):
        """A missing local MAFFT binary becomes a user-facing input error."""
        records = parse_sequences(">a\nACD\n>b\nACDE\n")

        with patch("msa.shutil.which", return_value=None):
            with self.assertRaisesRegex(MSAInputError, "MAFFT is not installed"):
                align_sequences(records)

    def test_format_fasta(self):
        """FASTA output uses one header and wrapped sequence per record."""
        output = format_fasta([("a", "ACD-"), ("b", "ACDE")])

        self.assertEqual(output, ">a\nACD-\n>b\nACDE\n")


if __name__ == "__main__":
    unittest.main()
