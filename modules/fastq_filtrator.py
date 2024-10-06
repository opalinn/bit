def filter_by_quality(
    seqs: dict, seqs_output: dict, quality_threshold: int = 0
) -> dict:
    for key, value in seqs.items():
        sum_of_scores = 0
        quality_string = value[1]
        for el in quality_string:
            sum_of_scores += ord(el) - 33
        mean_quality = sum_of_scores / len(quality_string)
        if mean_quality > quality_threshold:
            seqs_output[key] = value
    return seqs_output


def filter_by_length(seqs: dict, length_bounds: tuple, seqs_output: dict) -> dict:
    for key, value in seqs.items():
        if length_bounds[0] <= len(value[0]) <= length_bounds[1]:
            seqs_output[key] = value
    return seqs_output


def filter_by_gc_content(seqs: dict, gc_bounds: tuple, seqs_output: dict) -> dict:
    for key, value in seqs.items():
        gc_sum = value[0].count("G") + value[0].count("C")
        gc_content = (gc_sum / len(value[0])) * 100
        if gc_bounds[0] <= gc_content <= gc_bounds[1]:
            seqs_output[key] = value
    return seqs_output
