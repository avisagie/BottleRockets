from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
from rocket_architectures import Traces
import pandas as pd
import numpy as np
from rocket_architectures import sim_single_bottle


def create_scatter_plots_from_trace(name, trace: Traces):
    pos = go.Scatter(
        x=trace.position[:, 0],
        y=trace.position[:, 1],
        hovertext=[f"{time:.2f}s" for time in trace.time],
        name=name + " position",
    )
    spd = np.linalg.norm(trace.velocity, axis=1, keepdims=True)
    speed = go.Scatter(
        x=trace.time,
        y=spd[:, 0],
        name=name + " speed",
    )
    # direction_of_movement = trace.velocity/spd
    magnitude_of_acceleration = np.linalg.norm(trace.acceleration, axis=1, keepdims=True)
    # direction_of_acceleration = trace.acceleration/magnitude_of_acceleration
    # forward_acceleration = np.array(magnitude_of_acceleration)
    # for i in range(len(direction_of_acceleration)):
    #     dot_v_a = direction_of_movement[i].T@direction_of_acceleration[i]
    #     forward_acceleration[i]*=dot_v_a

    acc = go.Scatter(
        x=trace.time,
        y=magnitude_of_acceleration[:, 0] / 9.8,
        name=name + " acceleration",
    )
    return pos, speed, acc


def plot_single_trace(test_name: str, trace: Traces) -> go.Figure:
    fig = make_subplots(
        rows=2,
        cols=1,
        subplot_titles=("Position", "Speed/Acceleration"),
        specs=[[{}], [{"secondary_y": True}]],
    )
    pos, speed, acc = create_scatter_plots_from_trace(test_name, trace)

    fig.add_trace(pos, row=1, col=1)

    fig.add_trace(speed, row=2, col=1, secondary_y=False)

    fig.add_trace(
        acc,
        row=2,
        col=1,
        secondary_y=True,
    )

    fig.update_xaxes(title_text="Horizontal Distance (m)", row=1, col=1)
    fig.update_yaxes(title_text="Altitude (m)", row=1, col=1)

    fig.update_xaxes(title_text="Time", row=2, col=1)
    fig.update_yaxes(title="speed (m/s)", secondary_y=False, row=2, col=1)
    fig.update_yaxes(title="acceleration (g)", secondary_y=True, row=2, col=1)

    return fig


def plot_multiple_traces(test_names: list[str], traces: list[Traces]) -> go.Figure:
    fig = make_subplots(
        rows=3,
        cols=1,
        subplot_titles=("Position", "Speed", "Acceleration"),
    )

    for test_number, (name, trace) in enumerate(zip(test_names, traces)):
        scatters = create_scatter_plots_from_trace(name, trace)
        for i, scatter in enumerate(scatters):
            fig.add_trace(scatter, row=i + 1, col=1)
            fig.data[-1].update(marker=dict(color=px.colors.qualitative.Plotly[test_number]))

    fig.update_xaxes(title_text="Horizontal Distance (m)", row=1, col=1)
    fig.update_yaxes(title_text="Altitude (m)", row=1, col=1)

    fig.update_xaxes(title_text="Time", row=2, col=1)
    fig.update_yaxes(title="speed (m/s)", row=2, col=1)
    fig.update_yaxes(title="forward acceleration (g)", row=3, col=1)

    return fig


if __name__ == "__main__":
    n = 100
    test_trace = sim_single_bottle()

    fig = plot_single_trace("sim 1", test_trace)
    fig.show()

    fig = plot_multiple_traces(["sim 1"] * 2, [test_trace] * 2)
    fig.layout.autosize = True
    fig.show()
