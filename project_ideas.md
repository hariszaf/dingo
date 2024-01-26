<!--  Bio-projects for GSoC2024  -->

# Leveraging the increased statistical value of flux sampling

## Overview

The project will introduce new features for further statistical analysis of the samples and visualizations generated by `dingo`.

These will include at least the following features: 
- inference of pairwise correlated reactions
- construction of a weighted graph of the model's reactions with the correlation coefficients as weights
- annotation of these weights to the metabolic model and extraction to an annotated SBML file


All methods will be implemented in Python and merged in the `dingo` library.

The contributor will also run experiments with several metabolic networks to investigate the scaling of their findings.


## Related work

In constraint-based metabolic modelling, physical and biochemical constraints define a polyhedral convex set of feasible flux vectors. 

Contrary to Flux Balance Analysis (FBA), flux sampling is not dependent on an objective function [1]. 
Thus, sampling of this set provides an unbiased characterization of the metabolic capabilities of a biochemical network. 
Moreover, by sampling a sufficient number of samples one can study the properties of certain components of the whole network and deduce significant biological insights such as correlated reactions and/or pathways and more.

To sample uniformly from a convex polytope, `dingo` uses Multiphase Monte Carlo Sampling method based on Billiard walk [2].

Metabolic models are usually in SBML format and can be visualized through [Cytoscape](https://cytoscape.org).

For more information contact the mentors.


**References:**

[1] Apostolos Chalkis, Vissarion Fisikopoulos, Elias Tsigaridas, Haris Zafeiropoulos, [Geometric algorithms for sampling the flux space of metabolic networks](https://drops.dagstuhl.de/storage/00lipics/lipics-vol189-socg2021/LIPIcs.SoCG.2021.21/LIPIcs.SoCG.2021.21.pdf), 2021. 

[2] Wedmark, Y. K., Vik, J. O., & Øyås, O. (2023). [A hierarchy of metabolite exchanges in metabolic models of microbial species and communities](https://doi.org/10.1101/2023.09.05.556413). bioRxiv, 2023-09.


## Details of your coding project

The contributor will have to initiate the a post-process statistical analysis of the returned samples.
Then, they will extend dingo's [`illustrations.py`](https://github.com/GeomScale/dingo/blob/develop/dingo/illustrations.py) to visualize reaction pairs found correlated.
Both [Plotly](https://plotly.com) and [plotnine](https://plotnine.readthedocs.io/en/v0.12.4/index.html) libraries will be considered. 

Then they will have to run a few experiments on benchmark metabolic networks to assess how their methods scale in real-world metabolic networks and write a brief report with the results.

Difficulty: Medium


## Expected impact

The project will provide great help in the interpretation of the sampling findings. 
This benefits both the biologists community as they would gain novel insight and the geometry community highlighting the added value of the random sampling methods.
Also, it brings together GeomScale Org. with the NRNB community supporting Cytoscape.


## Mentors

- [Haris Zafeiropoulos](https://hariszaf.github.io) <haris.zafeiropoulos at kuleuven.be> is working on metabolic modeling software development and applications as a post-doc in the [Lab of Systems Biology](http://msysbiology.com) at KU Leuven and has previous GSoC student experience (2021) and mentoring experience with GeomScale (2022) and NRNB (2023).


- [Apostolos Chalkis](https://tolischal.github.io) <tolis.chal at gmail.com> is an expert in statistical software, computational geometry, and optimization, and has previous GSoC student experience (2018 & 2019) and mentoring experience with GeomScale (2020 & 2021).


## Tests

Students, please do one or more of the following tests before contacting the mentors above.

- **Easy:** compile and run [`dingo`](https://github.com/GeomScale/dingo). Use the documentation to sample from the flux space of the [e_coli](https://github.com/GeomScale/dingo/tree/develop/ext_data) model.

- **Medium:** Compare FBA solution with your samples. Choose radom pairs of reactions and check if they are correlated or not.

- **Hard:** For a pair of correlated reactions, show whether all the reactions of the metabolic pathway they belong to are also correlated. 



<!-- ======================================================= -->




# Metabolic interactions inference using community flux sampling

## Overview
 
The project will enhance `dingo` by incorporating various approaches for modeling microbial communities and conducting sampling in their flux space.

These will include at least: 
- a bag-of-reactions model: the union of all reactions and metabolites within at least one model of the community members 
- a bag of genomes model: each member of the community is depicted as a separate compartment, and their interactions are facilitated by shuttle reactions.


Flux sampling will be then conducted using the supported samplers within the `dingo` framework with or without a 
ginen set of elementary conversion modes.
The latter, would focus on the exchange reactions of the community model.

Cross-feeding interactions will be inferred from the samples based on the 
dependencies between fluxes of the pair.

The results obtained will undergo validation against experimental data. 
A comparative analysis of the employed modeling approaches will be performed.


## Related work

Community Genome-Scale models have been used for some time now [1].
Yet, there are several ways to model a community; each of whom having their own pros and cons [2]. 
Random flux sampling has become feasible for GEMs in recent years [3], but to scale pathway analysis it is necessary to focus on a subset of pathways or a subnetwork of the full metabolic network.

All features build is going to be in Puthon.
Methods such as `ecmtool` will be investigated to get elementary conversion modes.
All methods and analyses will be performn in Python.



**References:**


[1] Apostolos Chalkis, Vissarion Fisikopoulos, Elias Tsigaridas, Haris Zafeiropoulos, [Geometric algorithms for sampling the flux space of metabolic networks](https://drops.dagstuhl.de/storage/00lipics/lipics-vol189-socg2021/LIPIcs.SoCG.2021.21/LIPIcs.SoCG.2021.21.pdf), 2021. 


[2] Khandelwal, Ruchir A., et al. "Community flux balance analysis for microbial consortia at balanced growth." PloS one 8.5 (2013): e64567.


[2] Khandelwal, R. A., Olivier, B. G., Röling, W. F., Teusink, B., & Bruggeman, F. J. (2013). Community flux balance analysis for microbial consortia at balanced growth. PloS one, 8(5), e64567.






[3] Wedmark, Y. K., Vik, J. O., & Øyås, O. (2023). [A hierarchy of metabolite exchanges in metabolic models of microbial species and communities](https://doi.org/10.1101/2023.09.05.556413). bioRxiv, 2023-09.


<!-- https://doi.org/10.1016/j.csbj.2020.12.003
http://dx.doi.org/10.1098/rsif.2016.0627
cFBA: https://doi.org/10.1371/journal.pone.0064567  -->

<!-- Implementations -->
<!-- https://github.com/manuelgloeckler/ncmw/blob/main/ncmw/community/community_models.py -->
<!-- https://github.com/micom-dev/micom/blob/main/micom/community.py 
Sampling on EFM: https://gitlab.com/YlvaKaW/exchange-enumeration/-/blob/main/multicellular_analysis_ove.ipynb -->


<!-- [2] https://github.com/micom-dev/micom/blob/main/micom/community.py 
Sampling on EFM: https://gitlab.com/YlvaKaW/exchange-enumeration/-/blob/main/multicellular_analysis_ove.ipynb -->




## Details of your coding project

ecmtool can be used for get the elemntary flux (conversion) modes. 
the minimal set of conformal generators of a conversion cone (Urbanczik and Wagner, 2005; Clement et al., 2021). Like EFPs, ECMs can be

Methods such as `ecmtool` will be investigated to get elementary conversion modes.


Already implemented algorithms, regarding both the creation of the community model and the enumeration of the conversion modes will be considered for that, 
e.g. (/exchange-enumeration/multicellular_analysis_ove.ipynb)(https://gitlab.com/YlvaKaW/exchange-enumeration/-/blob/main/multicellular_analysis_ove.ipynb)



Difficulty: Hard



## Expected impact

The scientific impact of this project could be fundamental.
Only recently was explained how flux sampling when focusing on exchange reactions, could be implemented in community models.
This will be the first implementation of this feature in Python and in large scale.


## Mentors

- [Haris Zafeiropoulos](https://hariszaf.github.io) <haris.zafeiropoulos at kuleuven.be> is working on metabolic modeling software development and applications as a post-doc in the [Lab of Systems Biology](http://msysbiology.com) at KU Leuven and has previous GSoC student experience (2021) and mentoring experience with GeomScale (2022) and NRNB (2023).


- [Apostolos Chalkis](https://tolischal.github.io) <tolis.chal at gmail.com> is an expert in statistical software, computational geometry, and optimization, and has previous GSoC student experience (2018 & 2019) and mentoring experience with GeomScale (2020 & 2021).


## Tests

Students, please do one or more of the following tests before contacting the mentors above.

- **Easy:** compile and run [`dingo`](https://github.com/GeomScale/dingo). Use the documentation to sample from the flux space of the [e_coli](https://github.com/GeomScale/dingo/tree/develop/ext_data) model.

- **Medium:** Sample the flux space of each species of a pair of models known from the literature. Describe potential cross-feeding reactions based on these individual samples. 

- **Hard:** Enumerate Elementary Conversion Modes on e_colie_Core. 

