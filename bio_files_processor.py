from typing import TextIO


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> TextIO:
    """
    The function of bio_files_processor.py
    It makes one long string out of several short strings

    Parameters:
    input_fasta (str): path to .fasta file with sequences
    output_fasta (str = None): path to file with output of function.
    Optional argument

    Returns:
    TextIO: file with long strings with sequences
    """
    if output_fasta is None:
        output_fasta = "converted_" + input_fasta
    else:
        output_fasta = output_fasta + ".fasta"
    with open(input_fasta, "r") as file, open(output_fasta, "w") as output_file:
        output_seqs = {}
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                name = line
                output_seqs[name] = ""
            elif name is not None:
                output_seqs[name] += line
            for key, value in output_seqs.items():
                output_file.writelines(f"\n{key}\n{value}")


def parse_blast_output(input_file: str, output_file: None) -> TextIO:
    """
    The function of bio_files_processor.py
    Searches in the file with BLAST output in the 'Description' column
    for the name of the first protein and writes it to a new file in
    alphabet order

    Parameters:
    input_file (str): path to .txt file with sequences
    output_file (str = None): path to file with output of function

    Returns:
    TextIO: file with protein sequences
    """
    with open(input_file, "r") as file, open(output_file, "w") as out_file:
        if output_file is None:
            output_file = "proteins_blast_output.txt"
        lst_output = []
        flag = 0
        for line in file:
            if line.startswith("Description"):
                while flag == 0:
                    prot_line = file.readline().strip()
                    lst_output.append(prot_line[:66])
                    flag += 1
            flag = 0
        lst_output = sorted(lst_output)
        for line in lst_output:
            out_file.write(line)
            out_file.write("\n")
