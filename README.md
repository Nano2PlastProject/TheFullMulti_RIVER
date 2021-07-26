This Jupyter Notebook is associated to the Repository Full_Multi_Microplastics_River_Model_Test

Test before making the Full-Multi River Model available of the Full Multi Microplastics River Model open source version.

Please bare in mind this is a repository under developement and some sections migth be still incomplete and its update is an ongoing process.

Comments and feedback are wellcomed and will be adressed as time permits

# The Full Multi Microplastics River Model
=======================================

Microplastics River Model from Nano2Plast project

This repository contains one of the models developed in the ECO48-Nano2plast Project: [EXTENDING NANOPARTICLE MODELS TO OPEN SOURCE MODELS OF THE FATE AND TRANSPORT OF MICROPLASTIC IN AQUATIC SYSTEMS](https://bit.ly/35WxEPF).

The Full Multi Team: Proff. Matthew MacLeod (@MacLeodMatt), Dr. Antonia Praetorius (@apraetorius) and Dr. Maria del Prado Domercq (@PradoDomercq)

The code contained in this repository is intended as a framework to study microplastics fate and transport along a generic river system.

The model can be parameterised for different microplastic types (composition, shape, size, etc.) and river characteristics (flow velocity, suspended particulates, dimensions, etc.). Furthermore, this framework is intended as a flexible tool where new or updated fate processes can be included or reparameterised easily. 

This model follows a modular multimedia mass-balance modelling approach in which a system of coupled mass balance equations is built describing the transformation and transport processes of the plastic particles within the modelled system with specific rate constants. All processes are described using first-order kinetics based on our own review of the literature.

The modelled system consists of a generic one directional river structure where the river is subdivided in a set of river sections (RS) connected horizontally. Each river section is, at the same time, subdivided into four compartments of different types and dimensions representing a surface water layer (comp1), a flowing water body (comp2), a stagnant water body compartment (comp3) and the top sediment layer of the river bed (comp4). The system of differential equations is parameterized according to each compartment properties and dimensions and solved in a dynamic mode so that it yields the particles’s concentration (in number or mass of particles per volume) in each river section and compartment as a function of time.


![Alt text](https://github.com/PradoDomercq/Nano2Plast_RiverModel/blob/main/FigureGenericRiver.png "Generic River")


# Guide for Users

Upcoming

### Author
===========

Maria del Prado Domercq

prado.domercq@gmail.com
