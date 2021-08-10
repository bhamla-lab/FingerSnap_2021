# Finger Snap 2021 Github Guide

## Overview

This is the github created to house all the data files and model files used and created for the paper "The ultrafast snap of a finger is mediated by skin friction." 

## Repository Structure


- Data Files
    - High Speed Kinematics Experiment Data
    - Force Dynamics Experiment Data
    - Force Dynamics Experiment Data
- Fingersnap LaMSA Models

## Data Files

_High Speed Kinematics Experiment Data_

This branch houses the data files created from analyzing the high speed videos we took. The files are further organized based on the snapper (three different people provided high speed footage of the snap) and are named through the following naming convention:
```sh
s(snapper number) (material) (trial number) (video fps) fps.xlsx
```

Each file contains x and y data for each of the five circular reflectors used in each frame. Therefore, every five rows represents data for each frame. From first to fifth row, the x and y position data in each frame correspond to the wrist joint, the base of the fingers (metacarpophalangeal joint), the first knuckle (proximal interphalangeal joint) of the middle finger, the second knuckle (distal interphalangeal joint) of the middle finger, and the tip of the middle finger of the snapper. See Fig 1 in the main paper for additional detail. Additional details on the experiment performed can be found in the Materials and Methods seciton of the main manuscript.

_Force Dynamics Experiment Data_

This branch contains the data files created during the experiments measuring changes in force dynamics when different materials cover the fingers. These files are named as follows:
```sh
(material) (trial number).xlsx
```
Each file contains the force measured by a Force Sensitive Resistor (See SI Figure 1) every 0.1 ms. Additional details on the experiment performed can be found in the Materials and Methods seciton of the main manuscript.

_Friction Experiment Data_

This branch contains the data files created during experiments designed to measure the friction coefficient of various material combinations. The files are named as follows:
```sh
(material) - (trial number).xlsx
```
Each file contains the force measured by a load cell (See SI Figure 2) at every point in time. Additional details on the experiment performed can be found in the Materials and Methods section of the main manuscript. 

#Finger snap LaMSA Model

Here you can find the components of the compressible, friction-based LaMSA model we developed to reflect and capture the trends observed experimentally. Additional details on the development of this model can be found in the Results and Discussion section of the main manuscript. This model includes two example files and eleven files which make up the overall model. Additional information on these files and the workings of the model can be found in the readme.md file within this branch. 
