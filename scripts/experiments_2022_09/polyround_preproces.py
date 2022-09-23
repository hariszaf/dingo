#!/usr/bin/python3

import os, sys, time, getopt
import numpy as np
import pickle
from PolyRound.api import PolyRoundApi
from PolyRound.static_classes.lp_utils import ChebyshevFinder
from PolyRound.settings import PolyRoundSettings
import hopsy


def polyround_preprocess(model_path):

    name = model_path.split("/")[-1]

    # Import model and create Polytope object
    polytope = PolyRoundApi.sbml_to_polytope(model_path)
    print("Polyope for network " + name + " was built.")

    # Make a settings object for the polyround library - optional
    settings = PolyRoundSettings()

    # Simplify the polytope
    start = time.time()
    simplified_polytope = PolyRoundApi.simplify_polytope(polytope)
    end   = time.time()
    time_for_simplification = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_simplification) + " sec to get simplified.")

    # Polytope transformation
    start = time.time()
    transformed_polytope = PolyRoundApi.transform_polytope(simplified_polytope)
    end   = time.time()
    time_for_transformation = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_transformation) + " sec to get transformed.")

    # Export simplified and transformed polytope as pickle file
    polytope_info = (
        transformed_polytope,
        name,
    )

    with open(
        "simpl_transf_polytopes/polytope_" + name + ".pckl", "wb"
    ) as polyround_polytope_file:
        pickle.dump(polytope_info, polyround_polytope_file)

    # Rounding
    start = time.time()
    rounded_polytope = PolyRoundApi.round_polytope(transformed_polytope)
    end   = time.time()
    time_for_rounding = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_rounding) + " sec to get rounded.")

    # Export rounded polytope as pickle file
    polytope_info = (
        rounded_polytope,
        name,
    )

    with open(
        "polyround_polytopes/polytope_" + name + ".pckl", "wb"
    ) as polyround_polytope_file:
        pickle.dump(polytope_info, polyround_polytope_file)

    return rounded_polytope, name


if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], "n:p:s:w:", ["network=", "psrf=", "polytope=", "walk="])
    current_directory = os.getcwd()
    file_name = None ; polytope = None
    print(opts)
    for i in opts:

        if i[0] == "--network":
            file_name = i[1]
            path_to_net = current_directory + "/" + file_name

        if i[0]  == "--psrf":
            psrf = i[1]

        if i[0] == "--walk":
            walk = i[1]

    if file_name:
        rpolytope, name = polyround_preprocess(path_to_net)

