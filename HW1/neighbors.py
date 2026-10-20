def hamming_distance(first: str, second: str) -> int:
    if len(first) != len(second):
        raise ValueError("Strings must have equal lengths")

    return sum(
        a != b
        for a, b in zip(first, second)
    )


def neighbors(pattern: str, d: int) -> set[str]:
    if d < 0:
        return set()

    if d == 0:
        return {pattern}

    if len(pattern) == 0:
        return {""}

    if d >= len(pattern):
        return {
            "".join(chars)
            for chars in __import__("itertools").product(
                "ACGT", repeat=len(pattern)
            )
        }

    neighborhood = set()

    suffix = pattern[1:]
    suffix_neighbors = neighbors(suffix, d)

    for text in suffix_neighbors:
        if hamming_distance(suffix, text) < d:
            for nucleotide in "ACGT":
                neighborhood.add(nucleotide + text)
        else:
            neighborhood.add(pattern[0] + text)

    return neighborhood
