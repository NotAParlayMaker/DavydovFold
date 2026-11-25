"""Generate model sequences with exact amino-acid composition."""

from collections.abc import Mapping
from random import Random

__all__ = ["design_sequences"]

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def design_sequences(
    composition: Mapping[str, int],
    *,
    fixed: Mapping[int, str] | None = None,
    count: int = 1,
    seed: int | None = None,
) -> list[str]:
    """Sample sequences with exact residue counts and one-based fixed positions.

    Candidates are sampled independently and may repeat. A seed reproduces
    results for the same inputs and Python version.
    """
    if type(count) is not int or count < 1:
        raise ValueError("count must be a positive integer")

    remaining = dict(composition)
    for residue, amount in remaining.items():
        if not isinstance(residue, str) or len(residue) != 1 or residue not in AMINO_ACIDS:
            raise ValueError(f"invalid amino-acid code: {residue!r}")
        if type(amount) is not int or amount < 0:
            raise ValueError(f"count for {residue} must be a nonnegative integer")

    length = sum(remaining.values())
    if length == 0:
        raise ValueError("composition must define at least one residue")

    template = [""] * length
    for position, residue in (fixed or {}).items():
        if type(position) is not int or not 1 <= position <= length:
            raise ValueError(f"fixed position must be between 1 and {length}: {position!r}")
        if not isinstance(residue, str) or remaining.get(residue, 0) == 0:
            raise ValueError(f"composition cannot supply fixed residue {residue!r} at {position}")
        template[position - 1] = residue
        remaining[residue] -= 1

    pool = [residue for residue in AMINO_ACIDS for _ in range(remaining.get(residue, 0))]
    random = Random(seed)
    sequences = []
    for _ in range(count):
        shuffled = iter(random.sample(pool, len(pool)))
        sequences.append("".join(residue or next(shuffled) for residue in template))
    return sequences
