# NeuroFdpToolbox
MATLAB toolbox to automatically detect and analyse fractional LÃ©vy motions of propagating neural activity patterns, developed by Dr. Pulin Gong's group at University of Sydney. This toolbox includes the generation of a spiking neural network model, and the detailed processes for analysing our simulated and experimental data.
If you use our code in your research, please cite us as follows:

Liu Y., Long X., Martin PR, Solomon SG and Gong P., LÃ©vy walk dynamics explain gamma burst patterns in primate cerebral cortex. Under Review, 2020. 

## Network generation
A spiking neural circuit simulation model edited by Yifan Gu, Yuxi Liu, James Henderson, Guozhang Chen. The instruction is detailed in the readme file in the sub-folder.

SpikeNet is a software that has three stand-alone components.
1. User interface for configuring spiking neuronal networks
2. A c++ simulator 
3. User interface for parsing and post-analyzing the simulation results.

The design of SpikeNet provides the following four main features.

* **Configurability** SpikeNet supports any user-defined structure of synaptic connectivity topologies, coupling strengths and conduction delays. It can be easily extended by developers to support any variations of integrate-and-fire neuron and synapse models.

* **Performance**  Simulation of spiking neuronal network quickly becomes computationally intensive if the number of neurons in the network exceeds a few thousand. To achieve superior performance, various measures have been taken at both algorithmic and implementation level. 

* **User-friendly interface** In SpikeNet, although c++ is used for heavy-duty computation, its user-interface is written in  high-level programming language (Matlab) for user-friendliness and fast prototyping. This means SpikeNet does not require non-developer users to be familiar with c++.

* **Scalability** The design of the SpikeNet c++ simulator readily supports parallel computing using Message Passing Interface (MPI). Additionally, the HDF5-based I/O file format provides big data handling capability. Also Portable Batch System (PBS) scripts for array jobs are provided if the you have access to a cluster.

## Simulation data analysis
The matlab function for generating the neuron microcircuit is main_Gamma.m 
The detection, analysis and visualization of the fractional propagating patterns in the simulation data can be found in the folder model_data_analysis, such as GetBurst2.m; LFPAmpPattern5.m; MSDAnalysis.m; SpikeMatchLFPPhase.m; SpikesVSLFPPattern.m; Visualization3DLFPBurst.m; get_MSD_PBC.m, etc.
The code to generate the figure 4-6 in the paper "LÃ©vy walk dynamics explain gamma burst patterns in primate cerebral cortex" can be found in the sub-folder experimental_data_analysis/New/GammaPaperFig1-4.

## Experimental data analysis
The detection, analysis and visualization of the fractional propagating patterns in the experimental data can be found in the main function Project1.m in the sub-folder experimental_data_analysis/Toolbox_CSC. An example movie of the experimental data can be found:

![Example superdiffusive gamma burst pattern movie](https://github.com/longxian319/PhD_XL/blob/master/GammaDynaPatt/example%20movies/GammaBurstPatterns.avi)

The full data can be shared upon requested. The code to generate the figure 2 and figure 3 in the paper "LÃ©vy walk dynamics explain gamma burst patterns in primate cerebral cortex" can be found in the sub-folder experimental_data_analysis/Toolbox_CSC/p1.


## Authors

* **Yuxi Liu** - *Model generation and analysis* - [yliu2521](https://github.com/yliu2521)
* **Xian Long** - *Experimental data analysis* - [Xian Long](https://github.com/longxian319)
* **Pulin Gong** - *Coordinator* - pulin.gong@sydney.edu.au
