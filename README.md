[![CircleCI](https://dl.circleci.com/status-badge/img/gh/MDU-PHL/tbtamr/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/MDU-PHL/tbtamr/tree/master)
[![Python 3.9](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/downloads/release/python-390/)


# tbtAMR  <img src="https://github.com/MDU-PHL/tbtamr/blob/master/tbtamr_logo_transparent.png" width="100" height="70">

`tbtAMR` implements TB-profiler and custom logic developed at MDU to identify mutations linked to AMR mechanisms in _M. tuberculosis_ and generate reports suitable for public health in Victoria. It may also be suitable for use in research settings.

**`tbtAMR` is now accredited to ISO15189 standard by NATA for use in Victoria Australia.** 


## tbtAMR installation

In order to install `tbtAMR` conda is strongly recommended - installation instructions can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). 

```
conda create -n tbtamr python=3.9
```

Once you have installed `tbtAMR` you will need to setup the environment. This will install the validated version of `TB-profiler`, if you wish to use another version of `TB-profiler` you may get unexpected behaviour or errors.

```
conda activate tbtamr
tbtamr setup
tbtamr check
```

## Usage

See our [wiki](https://github.com/MDU-PHL/tbtamr/wiki) page for further information on `tbtamr` usage.
