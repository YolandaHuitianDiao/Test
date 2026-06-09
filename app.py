"""Streamlit front end for the small protein MSA tool."""

import streamlit as st

from msa import MSAInputError, align_sequences, format_fasta, parse_sequences


EXAMPLE_INPUT = """>human
MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE
>mouse
MKWVTFISLLFLFSSAYSREVFRRDAHKSEVAHRFKDLGE
>rat
MKWVTFISLLFLFSSAYSRGVFRRDAHKSEIAHRFKDLGE
"""


def main() -> None:
    """Render the UI, run MSA on button click, and display aligned FASTA."""
    st.set_page_config(page_title="Protein MSA", page_icon=":material/science:", layout="centered")
    st.title("Protein MSA")
    st.caption("Align up to 10 amino acid sequences, each up to 1000 residues.")

    text = st.text_area(
        "Sequences",
        value=EXAMPLE_INPUT,
        height=240,
        help="Paste FASTA records or one amino acid sequence per line.",
    )

    with st.sidebar:
        st.header("Alignment")
        st.write("Tool: MAFFT")
        st.write("Mode: auto, amino acid")
        st.header("Limits")
        st.write("Maximum sequences: 10")
        st.write("Maximum sequence length: 1000 aa")

    if not st.button("Align sequences", type="primary", use_container_width=True):
        return

    try:
        alignment = align_sequences(parse_sequences(text))
    except MSAInputError as exc:
        st.error(str(exc))
        return

    fasta = format_fasta(alignment)
    st.subheader("Aligned sequences")
    st.code(fasta, language="text")
    st.download_button("Download FASTA", fasta, "aligned_sequences.fasta", "text/plain")
    st.subheader("Alignment view")
    st.dataframe(
        [{"Name": name, "Aligned sequence": seq} for name, seq in alignment],
        hide_index=True,
        use_container_width=True,
    )


if __name__ == "__main__":
    main()
