[![Documentation Status](https://readthedocs.org/projects/streaminginstability-yj14/badge/?version=latest)](https://streaminginstability-yj14.readthedocs.io/en/latest/?badge=latest)
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/LGPL-3.0)

# Effects of Optically Thick Regions in Protoplanetary Disks 

We explore the observational implications of optically thick regions in protoplanetary disks. Our work includes a radiative transfer analysis, for details see the [documentation](https://streaminginstability-yj14.readthedocs.io/en/latest/).

# Analysis

We analyzed the numerical simulations of planet formation model via streaming instability presented by [Yang & Johansen (2014)](https://arxiv.org/pdf/1407.5995.pdf). This GitHub repo contains the density cube at one particular instant, which we used to construct our analysis pipeline. This data cube is from their simulations of the streaming instability, which employed a shearing-box of length 1.6 scale heights. 

This snapshot is available in the data folder and can be loaded as such:

```python

from StreamingInstability_YJ14 import shearing_box

cube = shearing_box.density_cube()
```



# How to Contribute?

Want to contribute? Bug detections? Comments? Suggestions? Please email us : godines@nmsu.edu, wlyra@nmsu.edu
