from unittest import TestCase, TestLoader, TextTestRunner

import sys
import cobra
import probanno_standalone.probanno as probanno

from contextlib import contextmanager

from six import StringIO

# Test Data Files
SIMPLE_1_MODEL = 'data/simple1_model.json'
SIMPLE_1_UNIVERSE = 'data/simple1_universe.json'
SIMPLE_1_REACTION_PROBABILITIES = 'data/simple1_reaction_probs.json'
SIMPLE_1B_REACTION_PROBABILITIES = 'data/simple1b_reaction_probs.json'

SIMPLE_2_MODEL = 'data/simple2_model.json'
SIMPLE_2_UNIVERSE = 'data/simple2_universe.json'
SIMPLE_2_REACTION_PROBABILITIES = 'data/simple2_reaction_probs.json'

if __name__ == "__main__":
    sys.path.insert(0, "../..")
    from cobra.core import Model, Reaction, Metabolite
    from cobra.flux_analysis import *
    sys.path.pop(0)
else:
    from cobra.core import Model, Reaction, Metabolite
    from cobra.solvers import get_solver_name
    from cobra.flux_analysis import *


@contextmanager
def captured_output():
    """ A context manager to test the IO summary methods """
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestProbabilisticAnnotation(TestCase):
    """Test the simulation functions of probabilistic annotation"""

    def setUp(self):
        pass

    def test_probabilistic_gapfill(self):
        """
        Test the probabilistic gap-filling approach
        :return:
        """
        # Adapted from test_gapfilling in test.cobra.flux_analysis
        try:
            solver = get_solver_name(mip=True)
        except:
            self.skipTest("no MILP solver found")
        # Simple Test Case 1
        model1 = cobra.io.json.load_json_model(SIMPLE_1_MODEL)
        universe1 = cobra.io.json.load_json_model(SIMPLE_1_UNIVERSE)
        rxn_probs1 = probanno.ReactionProbabilities.from_json_file(SIMPLE_1_REACTION_PROBABILITIES)
        reactions = probanno.probabilistic_gapfill(model1, universe1, rxn_probs1)
        reaction_ids = [r.id for r in reactions[0]]
        self.assertTrue('a2b' not in reaction_ids)
        self.assertTrue('a2d' in reaction_ids)
        self.assertTrue('d2b' in reaction_ids)

        # Simple Test Case 1b
        rxn_probs1b = probanno.ReactionProbabilities.from_json_file(SIMPLE_1B_REACTION_PROBABILITIES)
        reactions = probanno.probabilistic_gapfill(model1, universe1, rxn_probs1b)
        reaction_ids = [r.id for r in reactions[0]]
        self.assertTrue('a2b' in reaction_ids)
        self.assertTrue('a2d' not in reaction_ids)
        self.assertTrue('d2b' not in reaction_ids)

        # Simple Test Case 2
        model1 = cobra.io.json.load_json_model(SIMPLE_2_MODEL)
        universe1 = cobra.io.json.load_json_model(SIMPLE_2_UNIVERSE)
        rxn_probs2 = probanno.ReactionProbabilities.from_json_file(SIMPLE_2_REACTION_PROBABILITIES)
        reactions = probanno.probabilistic_gapfill(model1, universe1, rxn_probs2)
        reaction_ids = [r.id for r in reactions[0]]
        self.assertTrue('a2e' in reaction_ids)
        self.assertTrue('e2f' in reaction_ids)
        self.assertTrue('f2d' in reaction_ids)
        self.assertTrue('e2b' not in reaction_ids)
        self.assertTrue('a2b' not in reaction_ids)
        self.assertTrue('c2f' not in reaction_ids)
        self.assertTrue('c2d' not in reaction_ids)
        self.assertTrue('e2c' not in reaction_ids)








# make a test suite to run all of the tests
loader = TestLoader()
suite = loader.loadTestsFromModule(sys.modules[__name__])


def test_all():
    TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    test_all()