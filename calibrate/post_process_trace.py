import argparse

import numpy as np
from plotly import express as px
from scipy import signal

from rocket import Traces, load_traces, save_traces

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Post Process Trace",
        description="""
        Shifts position trace to start at 0.
        Converts velocity to speed and applies a low pass filter.
        Its important to filter out noisy sections by using the start and end index arguments.
        Stores the smoothed speed in the x direction of the velocity trace.
        Recalculates acceleration (scalar) from the smoothed speed.
        """,
    )

    parser.add_argument(
        "input_trace", help="Name of the json file in which the input trace is stored."
    )
    parser.add_argument(
        "--start_index",
        dest="start_index",
        type=int,
        default=0,
        help="Trim the trace by only reading data from the start index onward.",
    )
    parser.add_argument(
        "--end_index",
        dest="end_index",
        type=int,
        default=None,
        help="Trim the trace by not reading data at or after the end index.",
    )
    parser.add_argument(
        "--output_trace",
        dest="output_trace",
        type=str,
        default="smoothed.json",
        help="Name of the smoothed trace output file.help",
    )
    args = parser.parse_args()
    traces = load_traces(args.input_trace)
    start = args.start_index
    end = args.end_index if args.end_index is not None else len(traces.velocity)

    subset_traces = {k: v[start:end] for k, v in traces._asdict().items()}

    subset_speed = np.linalg.norm(subset_traces["velocity"], axis=1)
    padding = (20, 20)
    subset_speed = np.pad(subset_speed, padding, "edge")
    a, b = signal.butter(N=2, Wn=0.2)
    smoothed_speed = signal.filtfilt(a, b, subset_speed)
    px.line(y=[smoothed_speed, subset_speed], title="smoothed velocity").show()

    smoothed_speed = smoothed_speed[padding[0] : -padding[1]]
    smoothed_acceleration = np.diff(smoothed_speed) / np.diff(subset_traces["time"])

    px.line(smoothed_acceleration, title="smoothed acceleration").show()

    velocity = np.zeros(subset_traces["velocity"].shape)
    acceleration = np.zeros(subset_traces["acceleration"].shape)

    velocity[:, 0] = smoothed_speed
    acceleration[:-1, 0] = smoothed_acceleration
    position = subset_traces["position"] - np.min(subset_traces["position"], axis=0)
    smoothed_traces = Traces(
        time=traces.time[: len(smoothed_speed)],
        position=position,
        velocity=velocity,
        acceleration=acceleration,
    )
    save_traces(smoothed_traces, args.output_trace)
