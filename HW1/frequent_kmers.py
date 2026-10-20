from collections import defaultdict
from itertools import product

from reverse_complement import reverse_complement
from neighbors import neighbors


def frequent_words_with_mismatches_and_reverse_complements(
    text: str,
    k: int,
    d: int
) -> set[str]:
    if k <= 0 or k > len(text) or d < 0:
        return set()

    counts = defaultdict(int)

    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]

        for neighbor in neighbors(pattern, d):
            counts[neighbor] += 1
            counts[reverse_complement(neighbor)] += 1

    if not counts:
        return set()

    max_count = max(counts.values())

    return {
        pattern
        for pattern, count in counts.items()
        if count == max_count
    }