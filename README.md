# EXSAN: An (E)NDF (C)ross (S)ection (An)alysis and visualization application for nuclear sciences and engineering

## Welcome

Thank you for your interest in EXSAN. If you like the code, please send feedback to hasani@lanl.gov.

If EXSAN helps you in your research, please cite this [**open source article**](https://goo.gl/MPt6eY) from the Annals of Nuclear Energy:

Wooten, H., *"An application for streamlined and automated ENDF Cross Section Analysis and visualization,"* Ann Nuc Energy,     **129**, 482-486, 2019


## Getting Started

### Requirements
1. EXSAN is written in Python 2.7, and requires the [Python (2.7) Anaconda distribution](https://www.anaconda.com/distribution/).
2. Processing raw ENDF files requires an [NJOY 2016](https://www.njoy21.io/NJOY2016/) executable. Instructions for downloading and compiling the source code are [here](http://www.njoy21.io/Build/index.html)
3. After cloning the repository, create a subdirectory called `working`.  In `working`, create a softlink to your NJOY 2016 executable.

### Caveats
EXSAN has been tested on Mac OSX 10.11 and higher. It has **not** yet been tested on any other operating systems!
