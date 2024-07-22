# comm_detection_via_redundancy_optimization

Here you will find MATLAB code to perform community detection in human brain networks via redundancy maximization. Results are reported in the paper "Leveraging multivariate information for community detection in functional brain networks"

## Paper abstract

Embedded in neuroscience is the concept that brain functioning is underpinned by specialized systems whose integration enables cognition and behavior. Modeling the brain as a network of interconnected brain regions, allowed us to capitalize on network science tools and identify these segregated systems (modules, or communities) by optimizing the weights of pairwise connections within them. However, just knowing how strongly two brain areas are connected does not paint the whole picture. Brain dynamics is also engendered by interactions involving more areas at the same time, namely, higher-order interactions. In this paper, we propose a community detection algorithm that accounts for higher-order interactions and finds modules of brain regions whose brain activity is maximally redundant. Compared to modules identified with methods based on bivariate interactions, our redundancy-dominated modules are more symmetrical between the hemispheres, they overlap with canonical systems at the level of the sensory cortex, but describe a new organization of the transmodal cortex. By detecting redundant modules across spatial scales, we identified a sweet spot of maximum balance between segregation and integration of information as that scale where redundancy within modules and synergy between modules peaked. Moreover, we defined a local index that distinguishes brain regions in segregators and integrators based on how much they participate in the redundancy of their modules versus the redundancy of the whole system. Finally, we applied the algorithm to a lifespan dataset and tracked how redundant subsystems change across time. The results of this paper serve as educated guesses on how the brain organizes itself into modules accounting for higher-order interactions of its fundamental units, and pave the way for further investigation that could link them to cognition, behavior, and disease.

## Libraries

This project relies on external libraries:
- BCT toolbox: to extract measures on brain networks (https://sites.google.com/site/bctnet/).
- To plot communities on the human neocortex: https://github.com/faskowit/parc_plotter.
- To pick nice colors: https://it.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps.
- other nice colors: https://it.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps.
