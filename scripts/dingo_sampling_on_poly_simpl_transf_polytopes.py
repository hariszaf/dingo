import os, sys
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo import plot_histogram
from dingo import gmscale, apply_scaling
from time import process_time
import pickle
import multiprocessing

def sample_on_polyround_processed_polytope(input_polytope):

   current_directory = os.getcwd()
   name              = input_polytope.split("/")[-1]
   polytope_file     = input_polytope

   with open(polytope_file, "rb") as f:
      obj = pickle.load(f)
   simpl_polytope = obj[0]

   A = simpl_polytope.A.to_numpy()
   b = simpl_polytope.b.to_numpy()

   res = gmscale(A, 0.99)
   res = apply_scaling(A, b, res[0], res[1])

   scaled_A = res[0]
   scaled_b = res[1]


   # Sample from the polytope built using the parrallel MMCS 
   start         = process_time() 
   dingo_samples = PolytopeSampler.sample_from_polytope(scaled_A, 
							scaled_b,
                                                   	psrf = True, 
                                                   	parallel_mmcs = False
                                                  	)
   end = process_time() 
   sampling_runtime = end - start
   
   print("Dingo for the " + name + " model took " + str(sampling_runtime) + " sec, using time.procesess_time(), to sample using polyround processed polytope")

   with open("../steady_states/dingo_samples_" + name + ".pckl", "wb") as dingo_samples_file: 
         pickle.dump(dingo_samples, dingo_samples_file)


   steady_states = simpl_polytope.back_transform(dingo_samples)
   with open("../steady_states/dingo_steady_states_" + name + ".pckl", "wb") as dingo_steadystates_file:
      pickle.dump(steady_states, dingo_steadystates_file)



if __name__ == '__main__':

   file_name = sys.argv[1]
   current_directory = os.getcwd()
   path_to_net = current_directory + "/" + file_name
   sample_on_polyround_processed_polytope(path_to_net)



