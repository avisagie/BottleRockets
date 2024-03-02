import inspect
from collections import namedtuple

import rocket_architectures as ra

SimulatorTemplate = namedtuple("simulator", ["name", "arguments", "function"])


def get_simulators() -> list[SimulatorTemplate]:
    simulators_functions = [
        ra.sim_single_bottle,
        ra.sim_rocket_with_boosters,
        ra.sim_3_boosters_bullet,
        ra.sim_3_stage,
    ]
    simulators = [
        SimulatorTemplate(sim.__name__, get_signature(sim), sim) for sim in simulators_functions
    ]
    return simulators


def get_signature(function):
    func_parameters = inspect.signature(function).parameters
    return {key: value.default for key, value in func_parameters.items()}
