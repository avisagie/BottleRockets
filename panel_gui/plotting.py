from plotly.subplots import make_subplots
import plotly.graph_objects as go
from rocket_architectures import Traces
import pandas as pd
import numpy as np
from rocket_architectures import sim1


def plot_single_trace(test_name : str,traces : Traces) -> go.Figure:
    fig = make_subplots(
        rows = 2, 
        cols = 1,
        subplot_titles=("Position","Speed/Acceleration"),
        specs=[[{}],[{"secondary_y" : True}]],
    )

    fig.add_trace(
        go.Scatter(
            x = traces.position[:,0],
            y = traces.position[:,1],
            name = "position",
        ),
        row = 1,
        col = 1
    )

    fig.add_trace(
        go.Scatter(
            x = traces.time,
            y = np.linalg.norm(traces.velocity,axis=1),
            name = "speed",
        ),
        row = 2,
        col = 1,
        secondary_y = False
    )

    fig.add_trace(
        go.Scatter(
            x = traces.time,
            y = np.linalg.norm(traces.acceleration,axis=1)/9.8,
            name = "acceleration",
        ),
        row = 2,
        col = 1,
        secondary_y=True,
    )

    fig.update_xaxes(title_text = "Horizontal Distance (m)",row=1,col=1)
    fig.update_yaxes(title_text = "Altitude (m)",row=1,col=1)

    fig.update_xaxes(title_text="Time",row=2,col=1)
    fig.update_yaxes(title="speed (m/s)",secondary_y = False,row=2,col=1)
    fig.update_yaxes(title="acceleration (g)",secondary_y = True,row=2,col=1)

    return fig



if __name__ == "__main__":
    n = 100
    test_trace = sim1()

    fig = plot_single_trace("sim 1",test_trace)
    fig.show()
