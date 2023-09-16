if not __name__.startswith("bokeh_app"):
    raise RuntimeError("Start the gui with the command `panel serve panel_gui/gui.py`")

# add the root directory to the python path
import sys
from pathlib import Path
from icecream import ic
ic.disable()

parent_dir = Path(".").parent.absolute()
sys.path.append(str(parent_dir))

import panel as pn
from panel_gui.simulation_templates import SimulatorTemplate, get_simulators
from panel_gui.plotting import plot_single_trace, plot_multiple_traces

pn.extension("plotly")
pn.extension("jsoneditor")

simulation_results = []
simulation_display = pn.pane.Plotly(sizing_mode="stretch_both")
error_display = pn.pane.Alert("Ready",alert_type = "success")


def update_display():
    if len(simulation_results) == 1:
        simulation_display.object = plot_single_trace(
            simulation_results[-1][0], simulation_results[-1][1]
        )
    elif len(simulation_results) > 1:
        unzipped = zip(*simulation_results)
        names = next(unzipped)
        traces = next(unzipped)
        simulation_display.object = plot_multiple_traces(names, traces)

    simulation_display.object.layout.autosize = True


def run_simulation(name, arguments, function):
    try :
        simulation_results.append((name, function(**arguments)))
    except Exception as e:
        error_display.alert_type = "danger"
        error_display.object = str(e)
    else:
        error_display.alert_type = "success"
        error_display.object = "done"
        update_display()


def create_simulator_pane(simulator: SimulatorTemplate):
    arguments = pn.widgets.JSONEditor(mode="form", value=simulator.arguments)
    simulate_button = pn.widgets.Button(name="run")
    simulate_button.on_click(
        lambda e: run_simulation(simulator.name, arguments.value, simulator.function)
    )
    pane = pn.Column(simulator.name, arguments, simulate_button)
    return simulator.name, pane


tabs = [tuple(create_simulator_pane(sim)) for sim in get_simulators()]
tabs = pn.Tabs(*tabs)
main_row = pn.Row(tabs, pn.Column(error_display,simulation_display))
main_row.servable()
