# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2022 Apostolos Chalkis
# Copyright (c) 2022 Vissarion Fisikopoulos
# Copyright (c) 2022 Haris Zafeiropoulos
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.pyoptinterface_based_impl import set_default_solver

def sampling(model, testing_class):

    sampler = PolytopeSampler(model)

    test_cases = [
        ("mmcs", {"ess": 1000}),
        ("gaussian_hmc_walk", {"n": 500}),
        ("exponential_hmc_walk", {"n": 500, "variance": 50}),
        ("hmc_leapfrog_gaussian", {"n": 500}),
        ("hmc_leapfrog_exponential", {"n": 500, "variance": 50}),
    ]

    for method, kwargs in test_cases:
        steady_states = sampler.generate_steady_states_no_multiphase(method=method, **kwargs)

        testing_class.assertEqual(steady_states.shape[0], 95)
        testing_class.assertFalse(np.all(steady_states == 0))

class TestSampling(unittest.TestCase):

    def test_sample_json(self):

        input_file_json = os.getcwd() + "/ext_data/e_coli_core.json"
        model = MetabolicNetwork.from_json(input_file_json)
        sampling(model, self)

    def test_sample_mat(self):

        input_file_mat = os.getcwd() + "/ext_data/e_coli_core.mat"
        model = MetabolicNetwork.from_mat(input_file_mat)
        sampling(model, self)

    def test_sample_sbml(self):

        input_file_sbml = os.getcwd() + "/ext_data/e_coli_core.xml"
        model = MetabolicNetwork.from_sbml(input_file_sbml)
        sampling(model, self)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main()
