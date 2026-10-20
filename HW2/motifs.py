import random


def random_motifs(
    dna: list[str],
    k: int,
    rng: random.Random
) -> list[str]:
    """Choose one random k-mer from each DNA sequence."""
    if k <= 0:
        raise ValueError("k must be positive")

    motifs = []

    for sequence in dna:
        if len(sequence) < k:
            raise ValueError("Every DNA sequence must have length >= k")

        start = rng.randrange(len(sequence) - k + 1)
        motifs.append(sequence[start:start + k])

    return motifs

def profile_randomly_generated_kmer(
    text: str,
    k: int,
    profile: dict[str, list[float]],
    rng: random.Random
) -> str:
    """Sample a k-mer with probability proportional to its profile probability."""
    if k <= 0 or len(text) < k:
        raise ValueError("Invalid k or text shorter than k")

    windows = [
        text[i:i + k]
        for i in range(len(text) - k + 1)
    ]

    weights = []

    for window in windows:
        probability = 1.0

        for position, nucleotide in enumerate(window):
            probability *= profile[nucleotide][position]

        weights.append(probability)

    if sum(weights) == 0:
        return rng.choice(windows)

    return rng.choices(windows, weights=weights, k=1)[0]

def build_profile(
    motifs: list[str],
    k: int
) -> dict[str, list[float]]:
    """Build a profile with pseudocounts."""
    nucleotides = "ACGT"
    counts = {nucleotide: [1] * k for nucleotide in nucleotides}

    for motif in motifs:
        for position, nucleotide in enumerate(motif):
            counts[nucleotide][position] += 1

    denominator = len(motifs) + 4

    return {
        nucleotide: [
            count / denominator
            for count in counts[nucleotide]
        ]
        for nucleotide in nucleotides
    }


def score_motifs(motifs: list[str]) -> int:
    """Calculate the number of non-consensus symbols across columns."""
    if not motifs:
        return 0

    k = len(motifs[0])
    score = 0

    for position in range(k):
        counts = {
            nucleotide: sum(
                motif[position] == nucleotide
                for motif in motifs
            )
            for nucleotide in "ACGT"
        }

        score += len(motifs) - max(counts.values())

    return score

def gibbs_sampler(
    dna: list[str],
    k: int,
    n_iterations: int,
    rng: random.Random
) -> list[str]:
    """Run Gibbs sampling and return the best motifs found."""
    if not dna:
        raise ValueError("dna must not be empty")
    if k <= 0 or any(len(sequence) < k for sequence in dna):
        raise ValueError("Invalid k or DNA sequence shorter than k")
    if n_iterations < 0:
        raise ValueError("n_iterations must be non-negative")

    motifs = random_motifs(dna, k, rng)

    best_motifs = motifs.copy()
    best_score = score_motifs(best_motifs)

    for _ in range(n_iterations):
        i = rng.randrange(len(dna))

        # Exclude the motif belonging to the selected DNA sequence.
        remaining_motifs = motifs[:i] + motifs[i + 1:]

        profile = build_profile(remaining_motifs, k)

        new_motif = profile_randomly_generated_kmer(
            dna[i], k, profile, rng
        )
        motifs[i] = new_motif

        current_score = score_motifs(motifs)

        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs.copy()

    return best_motifs

def repeated_gibbs_sampler(
    dna: list[str],
    k: int,
    n_iterations: int,
    n_restarts: int,
    seed: int
) -> tuple[list[str], list[int]]:
    """Run Gibbs sampler repeatedly and return the best motifs and all scores."""
    if n_restarts <= 0:
        raise ValueError("n_restarts must be positive")

    rng = random.Random(seed)

    best_motifs = []
    best_score = float("inf")
    scores = []

    for _ in range(n_restarts):
        motifs = gibbs_sampler(dna, k, n_iterations, rng)
        current_score = score_motifs(motifs)
        scores.append(current_score)

        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs.copy()

    return best_motifs, scores