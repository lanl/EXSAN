# EXSAN: An (E)NDF (C)ross (S)ection (An)alysis and visualization application for nuclear sciences and engineering

## Welcome

Thank you for your interest in EXSAN. If you like the code, please send feedback to hasani@lanl.gov.

If EXSAN helps you in your research, please cite this [**open source article**](https://goo.gl/MPt6eY) from the Annals of Nuclear Energy:

Wooten, H., *"An application for streamlined and automated ENDF Cross Section Analysis and visualization,"* Ann Nuc Energy,     **129**, 482-486, 2019.


## Getting Started

**1. Install Anaconda**

EXSAN is written in Python 2.7, and requires the [Python (2.7) Anaconda distribution](https://www.anaconda.com/distribution/).  

**2. Add Anaconda to your path**

After Anaconda installation, make sure that Anaconda is pre-pended to your path. On Mac OSX, this can be done with by typing the following command in the terminal, and/or adding it to your `$HOME/.bash_profile` file:

`export PATH="$HOME/anaconda2/bin:$PATH"`

**3. Install NJOY2016**

Processing raw ENDF files requires an [NJOY 2016](https://www.njoy21.io/NJOY2016/) executable. Instructions for downloading and compiling NJOY are [here](http://www.njoy21.io/Build/index.html).

**4. Create a softlink for the NJOY2016 executable in the `EXSAN/working` subdirectory and name it `njoy`.**

E.g. `ln -s Path_to_NJOY/NJOY2016/bin/njoy working/njoy`

### Caveats
EXSAN has been tested on Mac OSX 10.11 and higher. It has **not** yet been tested on any other operating systems!

### Release
This software has been approved for open source release by Los Alamos National Laboratory as **LA-CC-19-002**.
