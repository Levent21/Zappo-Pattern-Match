# Zappo-Pattern-Match
A chemical logic alignment tool for protein sequences in the 'Twilight Zone'. Detects conserved Zappo patterns (5-mers to 8-mers) across divergent species, focusing on SARS-CoV-2 ORF10.
Zappo-Pattern-Matcher: Decoding the Protein "Twilight Zone"
🧬 The Challenge

Standard alignment tools (like MAFFT or Clustal Omega) often fail when sequence identity drops below 25% in certain regions (the "Twilight Zone"). In these cases, global alignments introduce artificial gaps and fail to align functional modules that have been shifted by insertions or deletions (indels) over millions of years of evolution.
💡 The Solution: Chemical Logic Alignment

This tool bypasses traditional letter-based matching by transforming amino acid sequences into Zappo Color Scheme integers (1-7). By focusing on chemical properties (hydrophobicity, charge, polarity) rather than specific residues, we can identify Hidden Homology.

Key Feature: The algorithm detects conserved n-mer patterns (5-mers to 8-mers) that are chemically identical but spatially shifted across different organisms (e.g., SARS-CoV-2 ORF10 vs. Plant/Fungal proteins).
🚀 Key Results with ORF10

Using this approach on the mysterious ORF10 protein (SARS-CoV-2), we identified conserved structural "anchors" that are invisible in standard Jalview/MAFFT alignments:

    Pattern 12511173 (8-mer): Found in ORF10 and a Platyhelminth protein, shifted by 16 residues.

    Pattern 5111 (4-mer): A universal structural motif found in 4 out of 5 highly divergent species, indicating a conserved alpha-helical turn.

🛠 How to Use

    Place your protein sequences in a FASTA format.

    Run python zappo_analyzer.py.

    The tool will output:

        Identity Boost: Shows how chemical similarity outweighs sequence identity.

        HSP Detector: Finds the longest identical chemical blocks and their "Shift" positions.

        Cluster Scanner: Identifies hydrophobic and charged regions.
