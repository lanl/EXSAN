# EXSAN: An (E)NDF (C)ross (S)ection (An)alysis and visualization application for nuclear sciences and engineering

## Welcome

Thank you for your interest in EXSAN. If you like the code, please send feedback to hasani@lanl.gov.

If EXSAN helps you in your research, please cite this [**open source article**](https://goo.gl/MPt6eY) from the Annals of Nuclear Energy:

Wooten, H., *"An application for streamlined and automated ENDF Cross Section Analysis and visualization,"* Ann Nuc Energy,     **129**, 482-486, 2019.


## Getting Started

**1. Install Anaconda**

EXSAN is written in Python 2.7, and requires the [Python (2.7) Anaconda distribution](https://www.anaconda.com/distribution/).  

**2. Add Anaconda to your path**

After Anaconda installation, make sure that Anaconda is pre-pended to your path. On Mac OSX, this can be done by typing the following command in the terminal, and/or adding it to your `$HOME/.bash_profile` file:

`export PATH="$HOME/anaconda2/bin:$PATH"`

**3. Install NJOY2016**

Processing raw ENDF files requires an [NJOY 2016](https://www.njoy21.io/NJOY2016/) executable. Instructions for downloading and compiling NJOY are [here](http://www.njoy21.io/Build/index.html). The minor modifications required to obtain NJOY2016 are as follows:
```
# Download the source code
git clone https://github.com/njoy/NJOY2016.git

# Get the desired version of NJOY21 (1.0.0 in this example)
cd NJOY2016
wget https://raw.githubusercontent.com/njoy/signatures/master/NJOY21/1.0.0-NJOY21.json
./metaconfigure/fetch_subprojects.py 1.0.0-NJOY21.json

# Configure the build process
mkdir bin
cd bin
cmake -D fetched_subprojects=true -D CMAKE_BUILD_TYPE=release ../

# Build NJOY1
make

# Test NJOY1
make test
```
Make note of the absolute path for the NJOY2016 executable, as this will be needed in the next steps.


**4. View EXSAN's help menu**

```
python exsan.py -h
```
One of the options, `-n`, lets you tell EXSAN where the NJOY2016 executable file is located (i.e., provide the absolute path to the NJOY2016 executable)

**5. Launch EXSAN and specify the NJOY2016 executable absolute path**

```
python exsan.py -n /Absolute/Path/To/NJOY2016/bin/njoy
```


### Tutorial
After you have NJOY installed and linked, this [tutorial](https://docs.google.com/presentation/d/1P2XkV6NQazMvSrXv0gHJSfl8rgEkq9fNbID7oCnMgA8/edit?usp=sharing) will walk you through how to use the main features of EXSAN.

### Caveats
EXSAN has been tested on Mac OSX 10.11 and higher. It has **not** yet been tested on any other operating systems!

### Release
This software has been approved for open source release by Los Alamos National Laboratory as **LA-CC-19-002**.
