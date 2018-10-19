# APTTipCarving
A tool to generate synthetic Atom Probe Tomography (APT) tips in preparation for evaporating them with the TAPSim.
For details see the following research papers:

"Building a Library of Simulated Atom Probe Data for Different Crystal Structures and Tip Orientations Using TAPSim"
M. Kühbach, A. Breen, M. Herbig, and B. Gault: Microscopy & Microanalysis, minor revision about to be submitted 2018/10/20

The input is currently tailored to the TAPSim simulation package:
C. Oberdorfer, S. M. Eich, and G. Schmitz: "A full-scale simulation approach for atom probe tomography"
Ultramicroscopy, 2013, 128, 55-67, dxdoi: 10.1016/j.ultramic.2013.01.005

The examples folder contains all scripts and key input xml files which exemplify how we used this tool for generating datasets.
Modify the CMakeLists.txt file to allow the CMake-powered compilation of the tool with a C/C++11 compatible compiler.
The run_queue_bicarving shell scripts specify how to execute the tool. 
So far only fcc is supported, key class function is SiXHdl.cpp were either an axis-angle "ax_dummy" or an orientation matrix "om_dummy" matrix can be passed. In case of questions of how to use the tool do not hesitate contacting me.
