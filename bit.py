import os
import click
import logging
from abc import ABC, abstractmethod
from typing import Self
from Bio import SeqIO, SeqRecord, Seq
from Bio.SeqUtils import gc_fraction


class BiologicalSequense(ABC):
    """
    Abstract base class representing a biological sequence with methods for length, indexing,
    string representation and checking correctness
    """

    @abstractmethod
    def __len__(self: Self) -> int:
        """
        Returns the length of the biological sequence
        """
        pass

    @abstractmethod
    def __getitem__(self: Self) -> Self:
        """
        Returns the slice or item at a given index in the sequence
        """
        pass

    @abstractmethod
    def __repr__(self: Self) -> str:
        """
        Returns a “official” representation of the sequence
        """
        pass

    @abstractmethod
    def __str__(self: Self) -> str:
        """
        Returns a string version of the sequence
        """
        pass

    @abstractmethod
    def check_correctness(self: Self) -> bool:
        """
        Checks if the biological sequence is correct
        """
        pass


class BioSequenceFunctions(BiologicalSequense):
    """
    Class for managing basic biological sequence functions, including length, indexing,
    and correctness checking

    """

    def __init__(self: Self, sequence: str) -> None:
        """
        Initializes the BioSequenceFunctions object

        Args:
            sequence (str): the biological sequence
        """
        self.sequence = sequence

    def __len__(self: Self) -> int:
        """
        Returns the length of object
        """
        return len(self.sequence)

    def __getitem__(self: Self, index: int | tuple) -> str:
        """
        Returns slice of sequence
        Args:
            index (int | tuple): the index to retrieve from the sequence

        Returns:
            str: The nucleotide or amino acid at the specified index
        """
        return self.sequence.__getitem__(index)

    def __repr__(self: Self) -> str:
        """
        Returns a “official” representation of the sequence
        """
        return self.sequence

    def __str__(self: Self) -> str:
        """
        Returns the string version of the sequence
        """
        return str(self.sequence)

    def check_correctness(self: Self) -> bool:
        """
        Validates the sequence against a set of characters

        Returns:
            bool: True if the sequence contains only valid characters, False otherwise
        """
        return set(self.sequence) <= set(self.alphabet)


class NucleicAcidSequence(BioSequenceFunctions):
    """
    Class for working with nucleic acid sequences using complementary and inverse operations
    """

    alphabet = None
    complement_dict = None

    def __init__(self: Self, sequence: str) -> None:
        """
        Initializes a NucleicAcidSequence object

        Args:
            sequence (str): a string with nucleic acid sequence
        """
        self.sequence = sequence

    def complement(self: Self) -> Self:
        """
        Returns the complement sequence

        Returns:
            NucleicAcidSequence: a new object containing the complement of the original sequence
        """
        complement_seq = "".join([self.complement_dict[bp] for bp in self.sequence])
        return type(self)(complement_seq)

    def reverse(self: Self) -> Self:
        """
        Returns the reversed sequence

        Returns:
            NucleicAcidSequence: a new sequence object with the reversed sequence
        """
        return type(self)(self.sequence[::-1])

    def reverse_complement(self) -> Self:
        """
        Returns the reverse complement copy of sequence

        Returns:
            NucleicAcidSequence: a new sequence object containing the reversed complement sequence
        """
        self.reverse()
        return self.complement()


class DNASequence(NucleicAcidSequence):
    """
    Class for handling DNA sequences with the transcribe method
    """

    alphabet = "ATGCatgc"
    complement_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
    }

    def __init__(self: Self, sequence: str):
        """
        Initializes a DNASequence object

        Args:
            sequence (str): a string with the DNA sequence
        """
        super().__init__(sequence)

    def transcribe(self: Self) -> Self:
        """
        Transcribes DNA into RNA (replaces "T" with "U")

        Returns:
            RNASequence: a new RNA sequence object
        """
        return RNASequence(self.sequence.replace("T", "U").replace("t", "u"))


class RNASequence(NucleicAcidSequence):
    """
    Class for representing RNA sequences
    """

    alphabet = "AUGCaugc"
    complement_dict = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
        "a": "u",
        "u": "a",
        "g": "c",
        "c": "g",
    }

    def __init__(self: Self, sequence: str):
        """
        Initializes an RNASequence object.

        Args:
            sequence (str): a string with the RNA sequence
        """
        super().__init__(sequence)


class AminoAcidSequence(BioSequenceFunctions):
    """
    Class for representing amino acid sequences with methods for counting amino acids in sequence
    """

    alphabet = "GAVLITSMCPFYWHKRDENQgavlitsmcpfywhkrdenq"

    def __init__(self: Self, sequence: str):
        """
        Initializes an AminoAcidSequence object

        Args:
            sequence (str): a string representing the amino acid sequence
        """
        super().__init__(sequence)

    def count_aa(self: Self) -> dict:
        """
        Counts the occurrences of each amino acid in the sequence

        Returns:
            dict: a dictionary with amino acids as keys and their counts as values
        """
        counter = dict()
        for aa in self.sequence:
            if aa not in counter.keys():
                counter[aa] = 0
            counter[aa] += 1
        return counter


logging.basicConfig(
    filename="fastq_filtrator.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

logger = logging.getLogger()


def make_bounds(bounds: str):
    try:
        return tuple(map(int, bounds.split(","))) if "," in bounds else (0, int(bounds))
    except Exception:
        raise click.BadParameter(
            f"Invalid bounds format: '{bounds}'. Use 50 or two numbers like 30,60"
        )


@click.command()
@click.argument("input_fastq", type=click.Path(exists=True))
@click.argument("output_fastq", type=click.Path())
@click.argument("gc_bounds", type=str)
@click.argument("length_bounds", type=str)
@click.argument("quality_threshold", type=int)
def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: str,
    length_bounds: str,
    quality_threshold: int,
):
    """
    Filters .fastq files by GC-content, lenght of reads and PHRED quality threshold

    Parameters:

    input_fastq : str
        Path to the input .fastq file

    output_fastq : str
        Path to the output .fastq file

    gc_bounds : str
        Range of GC content (in percentage)

    length_bounds : str
        Range of read length

    quality_threshold : int, default 0
        Minimum PHRED quality score for reads to be filtered

    Returns:
        None
    """
    try:
        logger.info(
            f"Start filtering with parameters: "
            f"input_fastq={input_fastq}, output_fastq={output_fastq}, "
            f"gc_bounds={gc_bounds}, length_bounds={length_bounds}, "
            f"quality_threshold={quality_threshold}"
        )

        gc_bounds = make_bounds(gc_bounds)
        length_bounds = make_bounds(length_bounds)

        records = SeqIO.parse(input_fastq, "fastq")

        filter_by_qual = (
            record
            for record in records
            if min(record.letter_annotations["phred_quality"]) >= quality_threshold
        )
        filter_by_length = (
            record
            for record in filter_by_qual
            if length_bounds[0] <= len(record.seq) <= length_bounds[1]
        )
        filter_by_gc = (
            record
            for record in filter_by_length
            if gc_bounds[0] <= gc_fraction(record.seq) * 100 <= gc_bounds[1]
        )
        SeqIO.write(filter_by_gc, output_fastq, format="fastq")
        logger.info("Filtering completed successfully")

    except Exception as error:
        logger.error(f"Error during run fastq_filtrator.py: {error}")
        raise click.ClickException(f"Error: {str(error)}")


if __name__ == "__main__":
    filter_fastq()
