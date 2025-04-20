import os
from Bio import SeqIO
from click.testing import CliRunner
from bit import filter_fastq, make_bounds

INPUT_FASTQ = os.path.join(os.path.dirname(__file__), "data_in.fastq")
OUTPUT_FASTQ = os.path.join(os.path.dirname(__file__), "data_out.fastq")


def test_make_bounds():
    assert make_bounds("30,70") == (30, 70)
    assert make_bounds("10") == (0, 10)


def test_input_file_exists():
    assert os.path.exists(INPUT_FASTQ), f"Input file {INPUT_FASTQ} does not exist"


def test_input_is_fastq():
    with open(INPUT_FASTQ, "r") as file:
        first_line = file.readline()
        assert first_line.startswith("@"), "Input file is not .fastq"


def test_gc_filtration_cli():
    if os.path.exists(OUTPUT_FASTQ):
        os.remove(OUTPUT_FASTQ)
    runner = CliRunner()
    runner.invoke(filter_fastq, [INPUT_FASTQ, OUTPUT_FASTQ, "0", "10", "30"])

    records = list(SeqIO.parse(OUTPUT_FASTQ, "fastq"))
    assert len(records) == 1, "Wrong filtration by GC content"


def test_lenght_filtration_cli():
    if os.path.exists(OUTPUT_FASTQ):
        os.remove(OUTPUT_FASTQ)
    runner = CliRunner()
    runner.invoke(filter_fastq, [INPUT_FASTQ, OUTPUT_FASTQ, "50", "20", "30"])

    records = list(SeqIO.parse(OUTPUT_FASTQ, "fastq"))
    assert len(records) == 2, "Wrong filtration by lenght"


def test_quality_filtration_cli():
    if os.path.exists(OUTPUT_FASTQ):
        os.remove(OUTPUT_FASTQ)
    runner = CliRunner()
    runner.invoke(filter_fastq, [INPUT_FASTQ, OUTPUT_FASTQ, "50", "20", "30"])

    records = list(SeqIO.parse(OUTPUT_FASTQ, "fastq"))
    assert len(records) == 2, "Wrong filtration by quality"


def test_multicondition_filtration_cli():
    if os.path.exists(OUTPUT_FASTQ):
        os.remove(OUTPUT_FASTQ)
    runner = CliRunner()
    runner.invoke(filter_fastq, [INPUT_FASTQ, OUTPUT_FASTQ, "0,100", "0,20", "30"])

    records = list(SeqIO.parse(OUTPUT_FASTQ, "fastq"))
    assert len(records) == 3, "Wrong filtration by all conditions"


def test_filtration_empty_output_cli():
    if os.path.exists(OUTPUT_FASTQ):
        os.remove(OUTPUT_FASTQ)
    runner = CliRunner()
    runner.invoke(filter_fastq, [INPUT_FASTQ, OUTPUT_FASTQ, "75", "0,30", "100"])

    records = list(SeqIO.parse(OUTPUT_FASTQ, "fastq"))
    assert len(records) == 0, "Wrong filtration by all conditions"
