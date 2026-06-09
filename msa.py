"""Parsing, validation, and MAFFT-backed protein MSA."""

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

ALLOWED_AA = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")
MAX_SEQUENCES = 10
MAX_SEQUENCE_LENGTH = 1000
MAFFT_CMD = "mafft"
MAFFT_TIMEOUT_SECONDS = 120


class MSAInputError(ValueError):
    """Raised when submitted sequences cannot be aligned."""


@dataclass(frozen=True)
class SequenceRecord:
    """Named amino acid sequence after validation."""

    name: str
    sequence: str


def parse_sequences(text: str) -> list[SequenceRecord]:
    """Parse FASTA or one-sequence-per-line text and validate app limits."""
    text = text.strip()
    if not text:
        raise MSAInputError("Paste at least two amino acid sequences.")

    records = _parse_fasta(text) if text.lstrip().startswith(">") else _parse_plain(text)
    if len(records) < 2:
        raise MSAInputError("Provide at least two sequences to align.")
    if len(records) > MAX_SEQUENCES:
        raise MSAInputError(f"Provide no more than {MAX_SEQUENCES} sequences.")

    clean_records = []
    for record in records:
        sequence = "".join(record.sequence.upper().split()).replace("-", "")
        if not sequence:
            raise MSAInputError(f"{record.name} has no sequence.")
        if len(sequence) > MAX_SEQUENCE_LENGTH:
            raise MSAInputError(
                f"{record.name} is {len(sequence)} aa long; maximum is {MAX_SEQUENCE_LENGTH}."
            )
        invalid = sorted(set(sequence) - ALLOWED_AA)
        if invalid:
            chars = ", ".join(invalid)
            raise MSAInputError(f"{record.name} contains invalid amino acid character(s): {chars}")
        clean_records.append(SequenceRecord(record.name, sequence))
    return clean_records


def align_sequences(records: list[SequenceRecord]) -> list[tuple[str, str]]:
    """Run MAFFT and return aligned `(name, sequence)` pairs."""
    if not shutil.which(MAFFT_CMD):
        raise MSAInputError(
            "MAFFT is not installed. Deploy with packages.txt, or install it locally first."
        )

    aligned = _parse_fasta(_run_mafft(format_fasta((r.name, r.sequence) for r in records)))
    if len(aligned) != len(records):
        raise MSAInputError("MAFFT returned an unexpected number of aligned sequences.")
    return [(record.name, record.sequence) for record in aligned]


def format_fasta(alignment: Iterable[tuple[str, str]]) -> str:
    """Format records as wrapped FASTA text."""
    lines = [line for name, seq in alignment for line in [f">{name}", *_wrap(seq, 80)]]
    return "\n".join(lines) + "\n"


def _run_mafft(fasta: str) -> str:
    """Write FASTA to a temporary file, execute MAFFT, and return stdout."""
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as handle:
        handle.write(fasta)
        input_path = Path(handle.name)

    try:
        result = subprocess.run(
            [MAFFT_CMD, "--auto", "--amino", "--thread", "1", "--quiet", str(input_path)],
            capture_output=True,
            check=False,
            text=True,
            timeout=MAFFT_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as exc:
        raise MSAInputError("MAFFT timed out while aligning these sequences.") from exc
    finally:
        input_path.unlink(missing_ok=True)

    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip() or "unknown MAFFT error"
        raise MSAInputError(f"MAFFT failed: {detail}")
    return result.stdout


def _parse_fasta(text: str) -> list[SequenceRecord]:
    """Read FASTA records, allowing multiline sequences."""
    records = []
    name = None
    chunks = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append(SequenceRecord(name, "".join(chunks)))
            name = line[1:].strip() or f"sequence_{len(records) + 1}"
            chunks = []
        else:
            if name is None:
                raise MSAInputError("FASTA sequence data must appear after a header line.")
            chunks.append(line)

    if name is not None:
        records.append(SequenceRecord(name, "".join(chunks)))
    return records


def _parse_plain(text: str) -> list[SequenceRecord]:
    """Treat each nonblank line as an unnamed sequence."""
    sequences = [line.strip() for line in text.splitlines() if line.strip()]
    return [
        SequenceRecord(f"sequence_{index}", seq)
        for index, seq in enumerate(sequences, start=1)
    ]


def _wrap(text: str, width: int) -> list[str]:
    """Split long FASTA sequence lines to the requested width."""
    return [text[index : index + width] for index in range(0, len(text), width)]
