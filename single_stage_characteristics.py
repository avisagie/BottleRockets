import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from rocket_architectures import plot_basic, sim_single_bottle


def grid_search():
    shared_config = dict(
        radius=0.043,  # meters
        nozzle_radius=0.0105,  # meters
        launch_tube_length=0.25,  # meters
        theta=30,  # degrees
        rail_length=1.5,  # meters
        dry_mass=0.269,  # kg
        volume=2.36,  # l
        windspeed=-2.0,  # m/s
        extra_frontal_surface=0.0,
        timestep=0.0001,
        bottle_shape="naive",
        water_l=2.36 / 3,
        C_drag=0.245,
    )
    pressures = [p for p in np.linspace(2, 8, 4)]

    fig = go.Figure()
    sims_to_run = [dict(**shared_config, pressure=p) for p in pressures]
    traces = [sim_single_bottle(**sim_cf) for sim_cf in sims_to_run]
    distances = [trace.position[-1][0] for trace in traces]
    fig.add_trace(go.Scatter(mode="lines", x=pressures, y=distances, name=f"C"))

    fig.add_trace(go.Scatter(mode="markers", x=[8, 6, 4], y=[182, 139, 94], name="measurements"))
    fig.show()


def grid_search_weight():
    shared_config = dict(
        radius=0.052875,
        nozzle_radius=0.0105,
        launch_tube_length=0.25,
        theta=42,
        rail_length=3.0,
        extra_frontal_surface=0.001,
        timestep=0.0001,
        bottle_shape="naive",
        pressure=6,
        volume=2.25,
        water_l=2.25 / 3.0,
    )

    mass = np.linspace(0.01, 0.6, 100)
    fig = go.Figure()
    for w in np.array([0, -4.0, -8.0, -16.0, -32.0]) * -1:
        for c in [0.3]:
            sims_to_run = [
                dict(
                    **shared_config,
                    C_drag=c,
                    dry_mass=m,
                    windspeed=w,
                )
                for m in mass
            ]
            traces = [sim_single_bottle(**sim_cf) for sim_cf in sims_to_run]
            distances = [trace.position[-1][0] for trace in traces]
            fig.add_trace(go.Scatter(mode="lines", x=mass, y=distances, name=f"Wind :{w} m/s"))
            fig.update_layout(
                xaxis_title="mass (kg)",
                yaxis_title="distance (m)",
                title="Distance vs Mass (tail wind)",
            )

    # fig.add_trace(go.Scatter(mode="markers", x=[8, 6, 4], y=[182, 139, 94], name="measurements"))
    fig.show()


def grid_search_angle():  # units wil be in either meters,square meters,seconds,bars,litres or kilograms
    shared_config = dict(
        radius=0.052875,
        nozzle_radius=0.0105,
        launch_tube_length=0.25,
        rail_length=3.0,
        extra_frontal_surface=0.001,
        timestep=0.0001,
        bottle_shape="naive",
        pressure=6,
        volume=2.25,
        water_l=2.25 / 3.0,
        dry_mass=0.3,
    )

    angle = np.linspace(15, 60, 100)
    fig = go.Figure()
    for w in np.array([0, -4.0, -8.0, -16.0, -32.0]) * -1:
        for c in [0.3]:
            sims_to_run = [
                dict(
                    **shared_config,
                    C_drag=c,
                    theta=a,
                    windspeed=w,
                )
                for a in angle
            ]
            traces = [sim_single_bottle(**sim_cf) for sim_cf in sims_to_run]
            distances = [trace.position[-1][0] for trace in traces]
            fig.add_trace(go.Scatter(mode="lines", x=angle, y=distances, name=f"Wind :{w} m/s"))
            fig.update_layout(
                xaxis_title="angle (degrees)",
                yaxis_title="distance (m)",
                title="Distance vs Angle (tail wind)",
            )

    # fig.add_trace(go.Scatter(mode="markers", x=[8, 6, 4], y=[182, 139, 94], name="measurements"))
    fig.show()


if __name__ == "__main__":
    # grid_search()
    grid_search_weight()
    grid_search_angle()
