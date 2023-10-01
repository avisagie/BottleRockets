#%%
from rocket import (
    Traces,
    save_traces,
    load_traces,
)
import numpy as np
from scipy import signal
from plotly import express as px
#%%
if __name__ == "__main__":
    #%%
    traces = load_traces("traces_from_final.mp4.json")
    #%%
    start = 16
    end = 106
    subset_traces = {k:v[start:end] for k,v in traces._asdict().items()}

    subset_speed = np.linalg.norm(subset_traces["velocity"],axis= 1)
    padding = (20,20)
    subset_speed = np.pad(subset_speed,padding,'edge')
    a,b = signal.butter(N=2,Wn=0.2)
    smoothed_speed = signal.filtfilt(a,b,subset_speed)
    px.line(y=[smoothed_speed,subset_speed],title = "smoothed velocity").show()

    smoothed_speed = smoothed_speed[padding[0]:-padding[1]]
    smoothed_acceleration = np.diff(smoothed_speed)

    #%%
    px.line(smoothed_acceleration,title = "smoothed acceleration").show()

    #%%
    smoothed_traces = Traces(time=subset_traces["time"],position=subset_traces["position"],velocity=np.expand_dims(smoothed_speed,-1),acceleration=smoothed_acceleration)
    save_traces(smoothed_traces,"smoothed.json")
