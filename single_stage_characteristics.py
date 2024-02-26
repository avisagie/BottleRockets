from rocket_architectures import sim_single_bottle,plot_basic
import numpy as np
import plotly.express as px


def rename_me():
    shared_config = dict(
        radius=0.04,
        dry_mass=0.269,
        volume=2.36,
        water_l=2.36 / 3.0,
        nozzle_radius=0.0105,
        launch_tube_length=0.025,
        theta=30,
        rail_length=1.5,
        extra_frontal_surface=0.0,
        timestep=0.0001,
        bottle_shape="naive",
        windspeed=-2.5,
   )
    #pressures = [p for p in np.linspace(2,8,100)]
    pressures = []
    dragC = []
    sims_to_run = []
    for p in [8,6,4]:
        for c in [0.1]:
            sims_to_run.append(dict(**shared_config,pressure=p,C_drag=c))
            pressures.append(p)
            dragC.append(c)
    traces = [sim_single_bottle(**sim_cf) for sim_cf in sims_to_run]
    distances = [trace.position[-1][0] for trace in traces]
    fig = px.scatter(dict(dist = distances,bar = pressures,c = dragC),x="bar",y="dist",hover_data="c")
    fig.show()
    # for trace in traces:
    #    plot_basic(trace)


 



if __name__ == "__main__":
    rename_me()
