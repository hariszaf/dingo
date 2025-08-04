# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2022 Apostolos Chalkis
# Copyright (c) 2022 Vissarion Fisikopoulos
# Copyright (c) 2022 Haris Zafeiropoulos
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file


import sys
import unittest
import numpy as np
from pathlib import Path

from dingo import MetabolicNetwork, PolytopeSampler
from dingo.pyoptinterface_based_impl import set_default_solver

root_dir = Path(__file__).parent.parent
ext_data = root_dir / "ext_data"

def sampling(model, testing_class):

    sampler = PolytopeSampler(model)

    test_cases = [
        ("billiard_walk", {"n": 500, "burn_in": 20, "thinning":100}),
        ("gaussian_hmc_walk", {"n": 500}),
        ("exponential_hmc_walk", {"n": 500, "variance": 50}),
        ("hmc_leapfrog_gaussian", {"n": 500}),
        ("hmc_leapfrog_exponential", {"n": 500, "variance": 50}),
    ]

    for method, kwargs in test_cases:

        steady_states = sampler.generate_steady_states_no_multiphase(method=method, **kwargs)
        try:
            testing_class.assertEqual(steady_states.shape[0], 95)
            testing_class.assertFalse(np.all(steady_states == 0))
        except AssertionError as e:
            print(f"❌ Test failed for method: {method}")
            raise

class TestSampling(unittest.TestCase):

    def test_sample_json(self):

        input_file_json = ext_data / "e_coli_core.json"
        model = MetabolicNetwork.from_json(input_file_json.as_posix())
        sampling(model, self)

    def test_sample_mat(self):

        input_file_mat = ext_data / "e_coli_core.mat"
        model = MetabolicNetwork.from_mat(input_file_mat.as_posix())
        sampling(model, self)

    def test_sample_sbml(self):

        input_file_sbml = ext_data / "e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml.as_posix())
        sampling(model, self)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main()
