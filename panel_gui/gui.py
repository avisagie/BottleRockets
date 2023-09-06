import panel as pn
from plotly import express as px
from simulation_templates import (
    SimulatorTemplate,
    get_simulators
)
from panel_gui.plotting import plot_single_trace,plot_multiple_traces

simulation_results = []
simulation_display = pn.pane.Plotly()
stacked_display = pn.pane.Plotly()

def update_display():
    if len(simulation_results) >= 1:
        simulation_display.object = plot_single_trace(simulation_results[-1][0],simulation_results[-1][1])
    
    if len(simulation_results) >= 2:
        unzipped = zip(*simulation_results)
        names = next(unzipped)
        traces = next(unzipped)
        stacked_display.object = plot_multiple_traces(names,traces)

def run_simulation(name,arguments,function):
    simulation_results.append((name,function(**arguments)))
    update_display()

def create_simulator_pane(simulator : SimulatorTemplate):

    arguments = pn.widgets.JSONEditor(
        mode = "form",
        value = simulator.arguments
    )
    simulate_button = pn.widgets.Button(name = "run")
    simulate_button.on_click(lambda e:run_simulation(simulator.name,arguments.value,simulator.function))
    pane = pn.Column(
        simulator.name,
        arguments,
        simulate_button
    )
    return simulator.name,pane
tabs = [tuple(create_simulator_pane(sim)) for sim in get_simulators()]
tabs = pn.Tabs(*tabs)
plots = pn.Column(simulation_display,stacked_display)
main_row = pn.Row(
    tabs,plots
)
main_row.servable()


