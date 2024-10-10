from .modules.run_dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
)
from .modules.filter_fastq import (
    filter_by_quality,
    filter_by_length,
    filter_by_gc_content,
)


def run_dna_rna_tools(*args: str) -> str | list:
    """
    The main function of the module dna_rna_tools.
    This function accepts an arbitrary number of arguments with DNA or RNA
    sequences (str) as input and the name of the procedure to be performed
    (this is always the last argument, str).

    Parameters:
    *args (str): strings with DNA or RNA sequences and name of procedure, which
    should be last

    Returns:
    str|list: modified strings under the entered procedure
    """
    lst_input = list(args)
    operation = lst_input.pop()

    def check_input(lst_input: list) -> list:
        """
        The helper function for run_dna_rna_tools().
        This function checks the list of entered strings

        Parameters:
        lst_input (list): list with input strings

        Returns:
        list: list with DNA or RNA sequences
        """
        valid_dna = set("ATGCatgc")
        valid_rna = set("AUGCaugc")
        lst_nucl_seq = []
        for seq in lst_input:
            if set(seq) <= valid_dna or set(seq) <= valid_rna:
                lst_nucl_seq.append(seq)
            else:
                lst_nucl_seq.append(None)
        return lst_nucl_seq

    function_map = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    res = function_map[operation](lst_nucl_seq)
    if len(res) == 1:
        return res[0]
    return res


def filter_fastq(
    seqs: dict,
    gc_bounds: int | tuple = (0, 100),
    length_bounds: int | tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
    """
    The main function of the module filter_fastq.
    This function takes arguments to filter fastq sequnces and return
    filtered seqs.

    Parameters:
    - seqs (dict): contains names of sequences (keys, str) and sequence
    with quality string (merged into tuple, values)
    - gc_bounds (int|tuple): argument used to filter sequences by their
    GC content, if int, it is converted to a tuple using check_input_args()
    - lenght_bounds (int|tuple): argument used to filter sequences by
    their length, if int, it is converted to a tuple using check_input_args()
    - quality_threshold (int): the argument for filtering reads by quality
    (the lower bound)

    Returns:
    dict: filtered sequences based on passed parameters

    """

    def check_input_args(gc_bounds: int | tuple, length_bounds: int | tuple) -> tuple:
        """
        The helper function for filter_fastq().
        This function converts arguments to tuples
        if they were passed as int

        Parameters:
        - gc_bounds (int | tuple)
        - length_bounds (int | tuple)

        Returns:
        tuples: gc_bounds and length_bounds
        """
        if isinstance(gc_bounds, int):
            gc_bounds = (0, gc_bounds)
        if isinstance(length_bounds, int):
            length_bounds = (0, length_bounds)
        return gc_bounds, length_bounds

    gc_bounds, length_bounds = check_input_args(gc_bounds, length_bounds)
    seqs_output = dict()
    seqs_output = filter_by_quality(seqs, quality_threshold)
    seqs_output = filter_by_length(seqs_output, length_bounds)
    seqs_output = filter_by_gc_content(seqs_output, gc_bounds)
    return seqs_output

