from rocket_architectures import sim_single_bottle


def rename_me():
    sim_config = dict(
        radius=0.04,
        C_drag=0.3,
        dry_mass=0.269,
        volume=2.36,
        water_l=2.36 / 3.0,
        nozzle_radius=0.0105,
        launch_tube_length=0.025,
        theta=30,
        rail_length=1.5,
        extra_frontal_surface=0.0,
        timestep=0.001,
        bottle_shape="naive",
        windspeed=-2.5,
    )


if __name__ == "__main__":
    rename_me()
