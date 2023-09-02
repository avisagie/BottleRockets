import panel as pn
from plotly import express as px
from simulation_templates import (
    SimulatorTemplate,
    get_simulators
)

simulation_results = []
simulation_display = pn.pane.Plotly()

def update_display():
    print(simulation_results)

def run_simulation(arguments,function):
    simulation_results.append(function(**arguments))
    update_display()

def create_simulator_pane(simulator : SimulatorTemplate):

    arguments = pn.widgets.JSONEditor(
        mode = "form",
        value = simulator.arguments
    )
    simulate_button = pn.widgets.Button(name = "run")
    simulate_button.on_click(lambda e:run_simulation(arguments.value,simulator.function))
    pane = pn.Column(
        simulator.name,
        arguments,
        simulate_button
    )
    return simulator.name,pane
tabs = [tuple(create_simulator_pane(sim)) for sim in get_simulators()]
tabs = pn.Tabs(*tabs)
tabs.servable()


