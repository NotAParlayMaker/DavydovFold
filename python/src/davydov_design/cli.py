"""Command-line FASTA output for the sequence design API."""

import argparse
import textwrap

from . import design_sequences


def _assignments(items, parse_key, parse_value):
    result = {}
    for item in items:
        key, value = item.split("=", 1)
        key = parse_key(key)
        if key in result:
            raise ValueError(f"duplicate assignment for {key!r}")
        result[key] = parse_value(value)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("composition", nargs="+", metavar="AA=COUNT", help="exact residue counts")
    parser.add_argument("--fixed", action="append", default=[], metavar="POSITION=AA",
                        help="one-based fixed position, repeatable")
    parser.add_argument("--count", type=int, default=1, help="number of candidates")
    parser.add_argument("--seed", type=int, help="seed for reproducible sampling")
    args = parser.parse_args(argv)
    try:
        composition = _assignments(args.composition, str, int)
        fixed = _assignments(args.fixed, int, str)
        sequences = design_sequences(composition, fixed=fixed, count=args.count, seed=args.seed)
    except ValueError as exc:
        parser.error(str(exc))

    for index, sequence in enumerate(sequences, start=1):
        print(f">candidate_{index}")
        print(textwrap.fill(sequence, width=70))
    return 0
