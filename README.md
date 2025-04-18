# BiT – Bioinformatics Tools for biologist

BiT is a tool for researchers who need to quickly analyze the results of their work: from obtaining complementary nucleic acid to filtering sequences

## Table of contents

[Installation](#installation)

[Available options in main script](#available-options-in-main-script)

[Available options in bio_files_processor script](#available-options-in-bio_files_processor-script)

[Usage](#usage)

[System requirements](#system-requirements)

[Contacts](#contacts)

## Installation

To install the BiT, use the following commands:

`git clone https://github.com/opalinn/bit.git`

`cd bit`

## Available options in main script

### Classes

1. BiologicalSequence(ABC)

An abstract base class representing a biological sequence with essential methods for length, indexing, string representation, and checking correctness.

Methods:

__len__(self) -> int: Returns the length of the biological sequence.

__getitem__(self) -> Self: Returns a slice or an item at a given index in the sequence.

__repr__(self) -> str: Returns an "official" representation of the sequence.

__str__(self) -> str: Returns a string version of the sequence.

check_correctness(self) -> bool: Checks if the biological sequence contains only valid characters based on the defined alphabet.

2. BioSequenceFunctions(BiologicalSequence)

This class manages basic biological sequence functions such as length, indexing, and correctness checking.

Methods:

__init__(self, sequence: str): Initializes the BioSequenceFunctions object with a given biological sequence.

__len__(self) -> int: Returns the length of the sequence.

__getitem__(self, index: int | tuple) -> str: Returns the nucleotide or amino acid at the specified index in the sequence.

__repr__(self) -> str: Returns a formal string representation of the sequence.

__str__(self) -> str: Returns the string version of the sequence.

check_correctness(self) -> bool: Validates the sequence against a set of characters defined by the alphabet.

3. NucleicAcidSequence(BioSequenceFunctions)

This class extends BioSequenceFunctions and provides functionality for working with nucleic acid sequences (DNA/RNA) using operations like complementing, reversing, and reverse complementing the sequence.

Methods:

complement(self) -> Self: Returns the complement of the nucleotide sequence.

reverse(self) -> Self: Returns the reversed sequence.

reverse_complement(self) -> Self: Returns the reverse complement sequence.

4. DNASequence(NucleicAcidSequence)

This class is specifically for DNA sequences and extends NucleicAcidSequence. It supports operations like transcribing DNA to RNA.

Methods:

__init__(self, sequence: str): Initializes the DNASequence object with a given DNA sequence.

transcribe(self) -> Self: Transcribes the DNA sequence into an RNA sequence by replacing "T" with "U".

5. RNASequence(NucleicAcidSequence)

This class is for RNA sequences and extends NucleicAcidSequence. It represents an RNA sequence and supports sequence manipulations.

Methods:

__init__(self, sequence: str): Initializes the RNASequence object with a given RNA sequence.

6. AminoAcidSequence(BioSequenceFunctions)

This class extends BioSequenceFunctions and provides functionality for amino acid sequences, including counting occurrences of each amino acid in the sequence.

Methods:
__init__(self, sequence: str): Initializes the AminoAcidSequence object with a given amino acid sequence.

count_aa(self) -> dict: Counts the occurrences of each amino acid in the sequence and returns the results in a dictionary format.

### Functions

`filter_fastq()`

This function filters .fastq files based on various criteria including GC-content, read length, and PHRED quality scores.

Parameters:

- input_fastq (str): Path to the input .fastq file.

- output_fastq (str): Path to the output .fastq file.

- gc_bounds (int | tuple, default (0, 100)): The range of GC content in percentage to filter by.

- length_bounds (int | tuple, default (0, 2^32)): The range of read lengths to filter by.

- quality_threshold (int, default 0): Minimum PHRED quality score for reads to be filtered.

Returns:

None: The filtered sequences are written to the output .fastq file.

## Available options in bio_files_processor script

`convert_multiline_fasta_to_oneline`

Makes one long string out of several short strings.

`parse_blast_output`

Searches in the file with BLAST output in the 'Description' column for the name of the first protein and writes it to a new file in alphabet order


## Usage

```python
import bit
# Create a DNA sequence object
dna_sequence = DNASequence("ATCGATTAG")

# Get the length of the sequence
print(len(dna_sequence))  # Output: 9

# Transcribe DNA to RNA
rna_sequence = dna_sequence.transcribe()
print(rna_sequence)  # Output: AUCGAUUAG

# Check sequence correctness
print(dna_sequence.check_correctness())  # Output: True
```

```python
convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
```

## System requirements

**Operating System**: Linux/Mac/Windows

**Python**: Version 3.12 or higher

**Biopython**: Version 1.85

Also take a healthy nervous system to figure out someone else's code 

## Contacts

If you have any ideas or you encounter a bug, please contact me: pp.malysheva@gmail.com

