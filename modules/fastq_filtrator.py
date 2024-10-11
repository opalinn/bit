def filter_by_quality(
    seqs: dict, quality_threshold: int = 0
) -> dict:
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
