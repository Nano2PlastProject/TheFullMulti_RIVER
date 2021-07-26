# The Full Multi Plastic Particles River Model
=======================================

River version of the Full Multi framework

This repository contains one of the models developed in the ECO48-Nano2plast Project: [EXTENDING NANOPARTICLE MODELS TO OPEN SOURCE MODELS OF THE FATE AND TRANSPORT OF MICROPLASTIC IN AQUATIC SYSTEMS](https://bit.ly/35WxEPF).

The Full Multi Team: Proff. Matthew MacLeod (@MacLeodMatt), Dr. Antonia Praetorius (@apraetorius) and Dr. Maria del Prado Domercq (@PradoDomercq)

The code contained in this repository is intended as a framework to study fate and transport of plastic particles along a generic river system.

The model can be parameterised for different plastic particle types (composition, shape, size, etc.) and river characteristics (flow velocity, suspended particulates, dimensions, etc.). Furthermore, this framework is intended as a flexible tool where new or updated fate processes can be included or reparameterised easily. 

This model follows a modular multimedia mass-balance modelling approach in which a system of coupled mass balance equations is built describing the transformation and transport processes of the plastic particles within the modelled system with specific rate constants. All processes are described using first-order kinetics based on our own review of the literature. A database of the modelled processess is included in the documentation folder of this respository.

The modelled system consists of a generic one directional river structure where the river is subdivided in a set of river sections (RS) connected horizontally. Each river section is, at the same time, subdivided into four compartments of different types and dimensions representing a surface water layer, a flowing water body, a stagnant water body compartment and the top sediment layer of the river bed. The system of differential equations is parameterized according to each compartment properties and dimensions and solved in a dynamic mode so that it yields the particles’s concentration (in number or mass of particles per volume) in each river section and compartment as a function of time.


![Alt text](https://github.com/PradoDomercq/Nano2Plast_RiverModel/blob/main/FigureGenericRiver.png "Generic River")


# Guide for Users

This model can be run in the cloud with Binder: you can run the FullMulti_RiverModel_Notebook Jupyter notebook directly from your web browser without installing anything. Just click on the launch binder button above. A temporary Jupyter Notebook server with all dependencies will be automatically launched in the cloud. It is not persistent: all your changes will be lost after some time. (#To be implemented#)

For a more interactive and modifiable version of this model follow the instructions below:

## Getting started for Windows

- First download and install Anaconda python distribution for your system from [here](https://www.anaconda.com/products/individual#windows) (scroll down to "Anaconda installers"), Python 3 is the recomended version for this model. If installation is sucessful, you will be able to work with the Jupyter Notebooks or use the Python Coding Environment Spyder to run the model. 
- Secondly, download the repository to your computer by clicking on the green CODE button on the rigth of the repository screen.
- Jupyter Notebooks: 
  - You can launch Jupyter notebook from visual Anaconda interface, or by opening `Anaconda prompt`, typing "jupyter notebook" and press Enter
  - Open the Jupyter notebook FullMulti_RiverModel_Notebook.ipynb and follow the instructions provided in the notebook to parameterise you scenario analysis and run the model.
  - It is recomended to make a copy of the notebook and rename it for each new model run/analysis
- Spyder:
  - You can launch Spyder from visual Anaconda interface, or by opening `Anaconda prompt`, typing "spyder" and press Enter
  - Open the MAIN_MPsRiver.py file in Spyder, follow the instructions on the code and run the model by pressing the green Play button

## Making changes in the model
- Reparameterization of imput files can be done by creating new imput files following the same structure as the original files in the imput folder or by adding new entries to the alaready existing imput files:

  - microplasticsSizeClass : contains all characteristics of MP particles to be modelled (i.e. composition, density, size and shape)
  - process_paramRiver: contains relevant MPs fate process descriptors organised by river compartment and MP aggregation state (e.g. biofilm growth rate, fragmentation rate constant, attachment efficiency)
  - compartmentsGenericRiverSec_prop: contains the parameters describing the Generic River dimensions and properties of each of its river sections and compartments (i.e. flow velocity, depth, width, length, etc.). A file to generate new river parameterizations is provided in the same folder called GenerateGenericRiverImputFile.
  
- Reformulation of process descriptions: the fate process descriptions for plastic particles included in this model are formulated in the file called RC_GeneratorRiver. New processes or alternative process descriptions can be added here and should be also added to the list of processes "processList" defined in the MAIN_MPsRiver.py file. Simultaneusly, in order for the model to take on account the new processes, new entries should be added in the estimation of rate constants section (starting line #140 of the MAIN_MpsRiver.py file).

### Author
===========

Maria del Prado Domercq

prado.domercq@gmail.com
