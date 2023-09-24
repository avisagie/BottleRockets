import argparse
from pathlib import Path
from from_video import read_video,show

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Interactive Video Frame',
        description='Display a single video frame in an interactive window using plotly',
    )
    parser.add_argument('filename')
    parser.add_argument('--frame_index',dest='frame_index',type=int,help = 'index of the frame to display')
    args = parser.parse_args()

    path_to_video = Path(args.filename)
    if not (path_to_video.exists() and path_to_video.is_file()):
        raise ValueError(f"Video file {path_to_video} not found.")
    print("Loading video")
    frames = read_video(str(path_to_video.absolute()))
    frame_index = args.frame_index if args.frame_index is not None else len(frames//2)
    if frame_index < 0 or frame_index > len(frames):
        raise ValueError(f"Invalid frame index of {frame_index} was selected for video with {len(frames)} frames.")
    show(frames[len(frames)//2],is_bgr=True)


