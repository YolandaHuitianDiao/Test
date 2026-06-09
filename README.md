# Protein MSA App

A minimalist Streamlit app for aligning a small number of amino acid sequences.

The app accepts FASTA records or one sequence per line, validates the input, runs
MAFFT, and returns aligned FASTA.

## Why Streamlit Community Cloud?

For this use case, Streamlit Community Cloud is the simplest free host:

- It deploys directly from a GitHub repository.
- It is built for small Python web apps.
- It manages the web server/container layer for you.
- Your workload is tiny: fewer than 10 sequences, each up to 1000 amino acids.

Hugging Face Spaces is also a good free option, especially if you want the app to
sit near other ML/bio demos. Render is more general-purpose, but its free web
services can spin down after inactivity, so it is less convenient for this app.

## Local Development

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Tests

```bash
python -m unittest discover -s tests
```

## Deploy On Streamlit Community Cloud

1. Push this folder to a GitHub repository.
2. Go to Streamlit Community Cloud.
3. Connect your GitHub account.
4. Create a new app from the repository.
5. Set the main file path to `app.py`.
6. Deploy. Streamlit will install Python packages from `requirements.txt` and
   MAFFT from `packages.txt`.

## Input Format

FASTA:

```text
>seq1
MKWVTFISLLFLFSSAYS
>seq2
MKWVTFISLLFLFSSAYT
```

Plain text, one sequence per line:

```text
MKWVTFISLLFLFSSAYS
MKWVTFISLLFLFSSAYT
```

## Notes

This implementation delegates alignment to MAFFT with `--auto --amino`, which is
appropriate for the small protein alignments targeted here. To run locally,
install MAFFT first, for example with Homebrew, Conda, or apt.
