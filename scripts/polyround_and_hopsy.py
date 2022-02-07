#!/usr/bin/python

import os, sys, getopt
import numpy as np
import pickle
from PolyRound.api import PolyRoundApi
from PolyRound.static_classes.lp_utils import ChebyshevFinder
from PolyRound.settings import PolyRoundSettings
import hopsy
from time import process_time

def polyround_preprocess(model_path):

    name = model_path.split("/")[-1]

    # Import model and create Polytope object
    polytope = PolyRoundApi.sbml_to_polytope(model_path)
    print("Polyope for network " + name + " was built.")

    # Make a settings object for the polyround library - optional
    settings = PolyRoundSettings()

    # Simplify the polytope
    start = process_time()
    simplified_polytope = PolyRoundApi.simplify_polytope(polytope)
    end   = process_time()
    time_for_simplification = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_simplification) + " sec to get simplified.")

    # Polytope transformation
    start = process_time()
    transformed_polytope = PolyRoundApi.transform_polytope(simplified_polytope)
    end   = process_time()
    time_for_transformation = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_transformation) + " sec to get transformed.")

    # Export simplified and transformed polytope as pickle file
    polytope_info = (
        transformed_polytope,
        name,
    )

    with open(
        "simpified_polytope_" + name + ".pckl", "wb"
    ) as polyround_polytope_file:
        pickle.dump(polytope_info, polyround_polytope_file)


    # Rounding
    start = process_time()
    rounded_polytope = PolyRoundApi.round_polytope(transformed_polytope)
    end   = process_time()
    time_for_rounding = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_rounding) + " sec to get rounded.")

    # Export rounded polytope as pickle file
    polytope_info = (
        rounded_polytope,
        name,
    )

    with open(
        "polyrounded_polytope_" + name + ".pckl", "wb"
    ) as polyround_polytope_file:
        pickle.dump(polytope_info, polyround_polytope_file)

    return rounded_polytope, name


def sample_on_rounded_polytope(name, rounded_polytope, psrf, walk_length):

    # Set the target distribution
    model = hopsy.UniformModel()

    # Build the problem using the rounded polytope and the model 
    start   = process_time()
    problem = hopsy.Problem(rounded_polytope.A, rounded_polytope.b, model)
    end     = process_time()

    # The run object contains and constructs the markov chains. 
    start                        = process_time()
    run                          = hopsy.Run(problem)

    # fixed for the case of Recon3D as mentioned in HOPS paper
    run.number_of_chains         = 5

    run.sample_until_convergence = True
    run.diagnostics_threshold    = float(psrf)

    if walk_length == "100d":
	    run.thinning                 = 100 * rounded_polytope.A.shape[1]
    elif walk_length == "8d2":
            run.thinning                 = 8 * rounded_polytope.A.shape[1] ^ 2 

    elif walk_length == "200d":
            run.thinning                 = 200 * rounded_polytope.A.shape[1]

    end                          = process_time()
    time_to_build_the_run        = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_to_build_the_run) + " sec for Run hopsy building.")


    # We finally sample
    start    = process_time()
    sampling = run.sample()
    end                          = process_time()
    time_for_run                 = end - start
    print("Polytope derived from the " + name + " network, took " + str(time_for_run) + " sec for hopsy sampling.")

    # From the run, we can now extract the produced data
    #data = run.data

    # Samples is a list of lists of numpy.ndarrays, which can be casted to a numpy.ndarray
    # which then has the shape (m,n,d), where m is the number of chains, n the number of samples and d the dimenion
    #samples  = data.states
    run.data.states = np.concatenate(run.data.states)
    #print("Number of samples retrieved: " + str(len(samples[0])))
    run.data.states = rounded_polytope.back_transform(run.data.states.T)
    #non_mapped_samples = np.concatenate(samples)
    #mapped_states = rounded_polytope.back_transform(non_mapped_samples.T)

    samples_info = (
        samples,
        name,
    )

    with open(
        "hopsy_samples/hopsy_samples_of_" + name + "_" + str(psrf) + "_" + str(walk_length) + ".pckl", "wb"
    ) as hopsy_samples_file:
        pickle.dump(samples_info, hopsy_samples_file)

    steady_states_info = (
            mapped_states,
            name,
    )

    with open(
        "steady_states/states_of_" + name + "_" + str(psrf) + "_" + str(walk_length) +  ".pckl", "wb"
    ) as hopsy_states_file:
        pickle.dump(steady_states_info, hopsy_states_file)

    return samples, mapped_states


if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], "n:p:s:w:", ["network=", "psrf=", "polytope=", "walk="])
    current_directory = os.getcwd()
    file_name = None ; polytope = None
    print(opts)
    for i in opts:

        if i[0] == "--network":
            file_name = i[1]
            #path_to_net = current_directory + "/" + file_name
            path_to_net = file_name

        if i[0]  == "--psrf":
            psrf = i[1]

        if i[0] == "--walk":
            walk = i[1]

        if i[0] == "--polytope":
            polytope_file = i[1]
            #polytope_file = current_directory + "/" + i[1]
            name = polytope_file.split("/")[-1]
            name = name.split(".")[0][9:]
            #print("here is a name: " + name)
            with open(polytope_file, "rb") as f:
                 obj = pickle.load(f)
            polytope = obj[0]

    if file_name:
        rpolytope, name = polyround_preprocess(path_to_net)
        sys.exit()
        sample_on_rounded_polytope(name, rpolytope, psrf, walk)

    if polytope:
        sample_on_rounded_polytope(name, polytope, psrf, walk)


