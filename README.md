# BottleRockets

Experiments in bottle rocket simulations.


# Getting  started

1. Install [mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install)
1. `git clone https://github.com/avisagie/BottleRockets.git`
1. `cd BottleRockets`
1. `mamba env create -f environment.yml`
1. `mamba activate bottlerockets`
1. `panel serve panel_gui/gui.py` to run simulations using the gui
1. You could also run `python -u fly.py`

Then check out fly.py, edit it with new parameters and have fun!
also check out rocket_architectures.py. There are two (at the time of writing):
one with a center stage and three boosters and one with a single stage.

fly*_ga.py try to optimise rocket parameters with a Genetic
Algorithm.
