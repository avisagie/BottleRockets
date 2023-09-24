#%%
import cv2 as cv
import numpy as np
import plotly.express as px
from icecream import ic
from pathlib import Path

def read_video(video_file : str) -> np.ndarray:
    cap = cv.VideoCapture(video_file)
    frames = []
    hsv_frames = []
    while cap.isOpened():
        ret, frame = cap.read()

        if not ret:
            break
        frames.append(frame) 
        hsv_frames.append(cv.cvtColor(frame,cv.COLOR_BGR2HSV))

    cap.release()
    return np.array(frames)

def bgr_to_hsv_colorspace(video_frames: np.ndarray) -> np.ndarray:
    hsv_frames = np.zeros(video_frames.shape)
    for i in enumerate(video_frames):
        hsv_frames[i,:] = cv.cvtColor(frames[i],cv.COLOR_BGR2HSV)

    return hsv_frames
    
def motion_detect_hsv(video_frames : np.ndarray):
    diff = video_frames.astype(int)
    # value_mask = diff[:,:,:,2] > 10
    # diff = diff*np.expand_dims(value_mask,axis=3)
    diff = np.diff(diff,axis=0)
    # the hue in hsv is circular and goes from 0 to 179
    diff[:,:,:,0] = ((diff[:,:,:,0] + 90) % 180 - 90)*2
    diff = np.abs(diff)  
    diff = np.sum(diff,axis=3)/3.0
    return diff.astype(np.uint8)

def motion_detect_bgr(video_frames : np.ndarray):
    diff = video_frames.astype(int)
    diff = np.diff(diff,axis=0)
    diff = np.abs(diff)
    diff = np.sum(diff,axis=3)/3
    return diff.astype(np.uint8)

def hue_filter(frames_to_mask : np.ndarray,hsv_video_frames : np.ndarray, hsv_start : np.ndarray, hsv_end : np.ndarray,invert_hue_range = False) -> np.ndarray:
    frames = np.array(hsv_video_frames)
    for i,frame in enumerate(hsv_video_frames):
        if not invert_hue_range  :
            mask = cv.inRange(frame,hsv_start,hsv_end)
        else:
            start_A = np.array(hsv_start)
            start_B = np.array(hsv_start)
            end_A = np.array(hsv_end)
            end_B = np.array(hsv_end)
            
            start_A[0] = 0
            end_A[0] = hsv_start[0]

            start_B[0] = hsv_end[0]
            end_B[0] = 179

            mask_A = cv.inRange(frame,start_A,end_A)
            mask_B = cv.inRange(frame,start_B,end_B)
            mask = cv.bitwise_or(mask_A,mask_B)

        frames[i,:] = cv.bitwise_and(frames_to_mask[i,:],frames_to_mask[i,:],mask = mask)

    return frames

def clamp(video_frames : np.ndarray,sensitivity = 0.5):
    video_frames = video_frames.astype(float)
    max_value = video_frames.max()
    video_frames = np.clip(video_frames - max_value*sensitivity,a_min = 0, a_max=255)
    video_frames = video_frames/(max_value - sensitivity*max_value)*255
    return video_frames.astype(np.uint8)

def save_mono_channel_video(video_file : str,video_frames) -> None:
    path = Path(video_file)
    if path.exists() and path.is_file():
        path.unlink()
    height,width = video_frames[0].shape
    print(width,height)
    writer = cv.VideoWriter(
        filename=video_file,
        fourcc=cv.VideoWriter.fourcc(*"mp4v"),
        frameSize=(width,height),
        fps=10,
        isColor= False,
    )
    for frame in video_frames:
        writer.write(frame)

    writer.release()

def save_color_video(video_file : str,video_frames) -> None:
    path = Path(video_file)
    if path.exists() and path.is_file():
        path.unlink()
    height,width,_ = video_frames[0].shape
    print(width,height)
    writer = cv.VideoWriter(
        filename=video_file,
        fourcc=cv.VideoWriter.fourcc(*"mp4v"),
        frameSize=(width,height),
        fps=10,
        isColor= True,
    )
    for frame in video_frames:
        writer.write(frame)

    writer.release()

def show(img,is_bgr = False):
    if is_bgr:
        img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

    px.imshow(img,binary_string=not is_bgr).show()

def track_largest_blob(video_frames : np.ndarray) -> tuple[list,np.ndarray] :
    frame_contours = []
    largest_contour_centers = np.zeros([len(video_frames),2])
    for i,frame in enumerate(video_frames):
        contours,_ = cv.findContours(frame,cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
        contours = sorted(contours,key=lambda x:cv.contourArea(x))
        if len(contours) > 0:
            largest_contour = contours[-1]
            M = cv.moments(largest_contour)
            M['m00'] += 1E-6
            cx = int(M['m10']/M['m00'])
            cy = int(M['m01']/M['m00'])
            frame_contours.append(contours)
            largest_contour_centers[i,:] = [cx,cy]
        else:
            frame_contours.append([])
            largest_contour_centers[i,:] = np.inf

    return frame_contours,largest_contour_centers

def add_track_and_contours(video_frames : np.ndarray, track : np.ndarray, contours_s : list):
    for i,(frame,contours) in enumerate(zip(video_frames[1:],contours_s)):
        n_contours = len(contours)
        if n_contours > 0:
            frame = cv.drawContours(frame, contours[-1], -1, (0,0,255), 3)
        if n_contours > 1:
            frame = cv.drawContours(frame, contours[:-1], -1, (0,255,255), 3)

        for point in track[:i]:
            frame = cv.circle(frame,point.astype(np.uint),3,(0,0,255))

def estimate_realworld_motion(track : np.ndarray, delta_time_seconds : float, pixel_width_meters : float):
    time = np.arange(len(track))*delta_time_seconds
    track = track*pixel_width_meters

    is_inf = np.isinf(track[:,0])
    time = time[~is_inf]
    track = track[~is_inf]

    dt = np.diff(time)
    dxdy = np.diff(track,axis=0)
    velocity = np.linalg.norm(dxdy,axis=1)/dt
    acceleration = np.diff(velocity)

    return track,time,velocity,acceleration

#%%
        
if __name__ == "__main__":
    # %%
    media_folder = "/mnt/c/Users/bottlerocket/"
    video_file = "single_bottle.mp4"
    frames= read_video(media_folder + video_file)

    #%%
    #red diff
    red_frames = frames.astype(int)
    red_frames = np.clip(red_frames[:,:,:,2] - red_frames[:,:,:,0] - red_frames[:,:,:,1],a_min=0,a_max=255)
    red_diff = np.diff(red_frames,axis=0)
    red_clamped = clamp(red_diff,0.1)
    #%%
    test_frame = red_clamped[50]
    contours,hierarchy = cv.findContours(test_frame,cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    show(cv.drawContours(np.array(test_frame), contours, -1, (255,255,255), 1))
    #%%
    contours,track = track_largest_blob(red_clamped)

    #%%
    add_track_and_contours(frames,track,contours)

    output_file = "final.mp4"
    save_color_video(media_folder + output_file,frames)
    #%%
    track,time,velocity,acceleration = estimate_realworld_motion(track,0.01,0.01)

