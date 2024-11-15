import os


def filter_by_quality(seqs: dict, quality_threshold: int = 0) -> dict:
    seqs_output = {}
    for key, value in seqs.items():
        sum_of_scores = 0
        quality_string = value[1]
        for el in quality_string:
            sum_of_scores += ord(el) - 33
        mean_quality = sum_of_scores / len(quality_string)
        if mean_quality > quality_threshold:
            seqs_output[key] = value
    return seqs_output


def filter_by_length(seqs_output: dict, length_bounds: tuple) -> dict:
    filtered_by_length = {}
    for key, value in seqs_output.items():
        if length_bounds[0] <= len(value[0]) <= length_bounds[1]:
            filtered_by_length[key] = value
    return filtered_by_length


def filter_by_gc_content(filtered_by_length: dict, gc_bounds: tuple) -> dict:
    filtered_by_gc = {}
    for key, value in filtered_by_length.items():
        gc_sum = value[0].count("G") + value[0].count("C")
        gc_content = (gc_sum / len(value[0])) * 100
        if gc_bounds[0] <= gc_content <= gc_bounds[1]:
            filtered_by_gc[key] = value
    return filtered_by_gc


def convert_fastq(input_fastq: str) -> dict:
    with open(input_fastq, "r") as file:
        seqs = {}
        for line in file:
            line = line.strip()
            if line.startswith("@S"):
                name = line
                seqs[name] = list()
            elif name is not None:
                seqs[name].append(line)
            elif line.startswith("+"):
                quality_info = line
            elif quality_info is not None:
                seqs[name].append(line)
        for key in seqs.keys():
            if len(seqs[key]) != 0:
                del seqs[key][1]
    return seqs


def write_output_file(output_fastq: str, output: dict) -> None:
    output_folder = "filtered"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_fastq = os.path.join("filtered", output_fastq)
    with open(output_fastq, "w") as file_output:
        for key, value in output.items():
            file_output.writelines(f"\n{key}\n{value[0]}\n{'+'}{key[1:]}\n{value[1]}")
