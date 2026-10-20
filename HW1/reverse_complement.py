def reverse_complement(pattern: str) -> str:
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    return "".join(
        complement[nucleotide]
        for nucleotide in reversed(pattern)
    )