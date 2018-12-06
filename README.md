# APTTipCarving
A tool to generate synthetic Atom Probe Tomography (APT) tips for evaporating them with the TAPSim simulation package.
For details see the following research papers:

"Building a Library of Simulated Atom Probe Data for Different Crystal Structures and Tip Orientations Using TAPSim"
M. Kühbach, A. Breen, M. Herbig, and B. Gault: in Microscopy & Microanalysis, 2019

"3D nanostructural characterisation of grain boundaries in atom probe data utilising machine learning techniques"
Y. Wei, Z. Peng, M. Kühbach, A. Breen, M. Legros, D. Raabe, and B. Gault, Frontiers, 2019

The input is currently tailored to the TAPSim simulation package:
C. Oberdorfer, S. M. Eich, and G. Schmitz: "A full-scale simulation approach for atom probe tomography"
Ultramicroscopy, 2013, 128, 55-67, dxdoi: 10.1016/j.ultramic.2013.01.005

The examples folder contains all scripts and key input xml files which exemplify how we used this tool for 
generating single-crystalline pillar datasets for the Microscopy and Microanalysis paper.
Further, it contains the xml settings file for generating a bicrystal setup with tailored 
boundary position and inclination angle that we used for the work in the Frontiers paper.
For generating bicrystal setups the reading of triangle patches is in beta stage and has not been tested.
We have NOT used it for our publication in Frontiers from January 2019. In this paper we used only the 
explicit planar boundary.

# Compiling
Modify the CMakeLists.txt file to allow the CMake-powered compilation of the tool with a C/C++11 compatible compiler.
The run_queue_bicarving shell scripts specify how to execute the tool. 
So far only fcc is supported, key class function is SiXHdl.cpp were either an axis-angle "ax_dummy" 
or an orientation matrix "om_dummy" matrix can be passed. 

# Questions
In case of questions of how to use the tool do not hesitate contacting me.
