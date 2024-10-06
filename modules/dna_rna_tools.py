complement_dict_dna = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
}

complement_dict_rna = {
    "A": "U",
    "U": "A",
    "G": "C",
    "C": "G",
    "a": "u",
    "u": "a",
    "g": "c",
    "c": "g",
}


def transcribe(lst_nucl_seq: str|list) -> str|list:
    lst_trscr_out = []
    for seq in lst_nucl_seq:
        if "U" in seq or "u" in seq:
            return "It is RNA"
        else:
            seq = seq.replace("T", "U").replace("t", "u")
            lst_trscr_out.append(seq)
    return lst_trscr_out if len(lst_trscr_out) > 1 else lst_trscr_out[0]


def reverse(lst_nucl_seq: str|list) -> str|list:
    lst_rev_out = []
    for seq in range(len(lst_nucl_seq)):
        lst_rev_out.append(lst_nucl_seq[seq][::-1])
    return lst_rev_out if len(lst_rev_out) > 1 else lst_rev_out[0]


def complement(lst_nucl_seq: str|list) -> str|list:
    lst_compl_out = []
    for seq in lst_nucl_seq:
        complement_output = ""
        if set(seq) <= set("ATGCatgc"):
            for bp in seq:
                complement_output += complement_dict_dna[bp]
            lst_compl_out.append(complement_output)
        elif set(seq) <= set("AUGCaugc"):
            for bp in seq:
                complement_output += complement_dict_rna[bp]
            lst_compl_out.append(complement_output)
    return lst_compl_out if len(lst_compl_out) > 1 else lst_compl_out[0]


def reverse_complement(lst_nucl_seq: str|list) -> str|list:
    lst_revcom_out = []
    for seq in lst_nucl_seq:
        if set(seq) <= set("ATGCatgc"):
            seq = seq[::-1]
            reverse_compl_out = "".join(complement_dict_dna[bp] for bp in seq)
            lst_revcom_out.append(reverse_compl_out)
        elif set(seq) <= set("AUGCaugc"):
            seq = seq[::-1]
            reverse_compl_out = "".join(complement_dict_rna[bp] for bp in seq)
            lst_revcom_out.append(reverse_compl_out)
    return lst_revcom_out if len(lst_revcom_out) > 1 else lst_revcom_out[0]
