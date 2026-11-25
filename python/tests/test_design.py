from collections import Counter
import subprocess
import sys
import unittest

from davydov_design import design_sequences


class SequenceDesignTests(unittest.TestCase):
    def test_seeded_design_preserves_composition_and_fixed_positions(self):
        composition = {"G": 3, "A": 4, "S": 3}
        args = {"fixed": {1: "A", 10: "S"}, "count": 8, "seed": 42}
        candidates = design_sequences(composition, **args)
        self.assertEqual(len(candidates), 8)
        self.assertEqual(candidates, design_sequences({"A": 4, "S": 3, "G": 3}, **args))
        for candidate in candidates:
            self.assertEqual(Counter(candidate), composition)
            self.assertEqual((candidate[0], candidate[-1]), ("A", "S"))

    def test_fully_fixed_composition(self):
        self.assertEqual(design_sequences({"A": 2}, fixed={1: "A", 2: "A"}), ["AA"])

    def test_impossible_constraints(self):
        for fixed in ({1: "A", 2: "A"}, {0: "A"}, {1: "G"}):
            with self.subTest(fixed=fixed), self.assertRaises(ValueError):
                design_sequences({"A": 1, "S": 1}, fixed=fixed)

    def test_cli_emits_fasta_and_reports_duplicate_assignments(self):
        command = [sys.executable, "-c", "from davydov_design.cli import main; main()"]
        result = subprocess.run(
            command + ["A=2", "S=1", "--fixed", "1=S"],
            capture_output=True, text=True, check=True,
        )
        self.assertEqual(result.stdout, ">candidate_1\nSAA\n")
        rejected = subprocess.run(
            command + ["A=2", "A=3"], capture_output=True, text=True,
        )
        self.assertEqual(rejected.returncode, 2)
        self.assertIn("duplicate assignment", rejected.stderr)
