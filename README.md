# BottleRockets

Experiments in bottle rocket simulations.


# Getting  started

1. Install [mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install)
1. `git clone https://github.com/avisagie/BottleRockets.git`
1. `cd BottleRockets`
1. `mamba env create -f environment.yml`
1. `mamba activate bottlerockets`
1. `panel serve panel_gui/gui.py` to run simulations using the gui. (zoom out on small screens with `ctrl -`)
1. You could also run `python -u fly.py`

Then check out `fly.py`, edit it with new parameters and have fun!
also check out rocket_architectures.py. There are two (at the time of writing):
one with a center stage and three boosters and one with a single stage.

# Some interesting ideas

- `fly*_ga.py` try to optimise rocket parameters with a Genetic
Algorithm.
- See `grid_search.py` to see how changing different parameters affect distance.

- If you have a slow motion video of a rocket launch you can play around with the programs in `calibrate/`.
You can run `python calibrate/<name_of_prog> --help` to see a description and some options.
It uses movement and color to detect the rocket.
It will most likely not properly track the rocket out of the box so you will have to jump into the code.

# Code formatting
The `black` code formatter with line length set to 100 was used.
Run `black --line-length 100 *.py`
