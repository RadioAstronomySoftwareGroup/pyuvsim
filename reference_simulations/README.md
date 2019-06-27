# Reference Simulations

``pyuvsim`` reference simulations are intended as a reference point for comparison
with existing simulators. These reference simulations span axes chosen to
highlight important physical directions (time, frequency, baselines, sources).
The total number of data points in each reference was constrained by the current
level of optimization and availability of computing resources. As performance
improves, we plan to add more reference simulations to cover multiple axes of
interest at once. Some of the reference simulations were run multiple times using
different beams of interest.

``pyuvsim`` is primarily run on the Oscar cluster at Brown University, but it has
been tested on other machines and clusters.
Our future plans are to restructure and optimize our code to allow for more
sources/times/frequencies/baselines to be simulated at a given time.
More documentation on ``pyuvsim`` can be found at [ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/).

Organization of this folder is as follows:
 - Information about the reference simulations can be found in the **Description** folder.
 - Catalog files can be found in the **catalog_files** folder.
 - Telescope configuration files can be found in the **catalog_config** folder.
 - Memos detailing the comparisons of reference simulations with other simulators
 and with HERA data can be found in the **Memos** folder.

Information about how to run the simulations may be found at
[ReadTheDocs](https://pyuvsim.readthedocs.io/en/latest/usage.html).
