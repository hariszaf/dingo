<!--  Bio-projects for GSoC2024  -->


# Leveraging the increased statistical value of flux sampling.

## Overview

The project will provide routines for furthere statistical analysis of the samples and visualizations 
statistical evidence and types of statistics highlight added value of sampling 


graph with associates fluxes scores

random env ball 


All methods will be implemented in Python and merged in the `dingo` library.

The contributor will also run experiments with several metabolic networks to investigate the scaling of their findings.


## Related work

Flux sampling is 





## Details of your coding project

The contributor will have to initiate the a post-process analysis of the samples return and extend the Python [`illustrations.py`]() to visualize reaction pairs found correlated.
Then they will have to run a few experiments on benchmark metabolic networks to assess how their methods scale in real-world metabolic networks and write a brief report with the results.

Difficulty: Medium





## Expected impact




## Mentors

- [Apostolos Chalkis](https://tolischal.github.io) <tolis.chal at gmail.com> is an expert in statistical software, computational geometry, and optimization, and has previous GSoC student experience (2018 & 2019) and mentoring experience with GeomScale (2020 & 2021).
- [Haris Zafeiropoulos](https://hariszaf.github.io) <haris.zafeiropoulos at kuleuven.be> is working on metabolic modeling software development and applications as a post-doc in the [Lab of Systems Biology](http://msysbiology.com) at KU Leuven and has previous GSoC student experience (2021) and mentoring experience with GeomScale (2022) and NRNB (2023).


## Tests

Students, please do one or more of the following tests before contacting the mentors above.

- Easy: compile and run [`dingo`](https://github.com/GeomScale/dingo). Use the documentation to sample from the flux space of the [e_coli](https://github.com/GeomScale/dingo/tree/develop/ext_data) model.
- Medium: Compare FBA solution with your samples. Identify reaction pairs with statistical significance.
- Hard: Generate a random polytope in Python and use `dingo` to sample from it using MMCS.









# Predicting metabolic interactions using sampling 





