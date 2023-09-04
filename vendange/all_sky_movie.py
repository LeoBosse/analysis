#!/usr/bin/python3
# -*-coding:utf-8 -*

################################################################################
# This script creates a movie from a set of images.
# It loads the images from a folder (with a filter on the names), transform them, and stack them in a movie.
################################################################################

import cv2
from astropy.io import fits
import argparse
import os
from vendange_configuration import Config
import numpy as np
import datetime as dt

config = Config()

image_folder = config.data_path + 'allsky_camera/KIL20200224/'  #plip/Temp"

# Construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-ext", "--extension", required=False, default='png', help="extension name. default is 'png'.")
ap.add_argument("-o", "--output", required=False, default='output.mp4', help="output video file")
ap.add_argument("-fps", "--framepersec", required=False, default=24, help="Number of frames per seconds for the output video.")
args = vars(ap.parse_args())

# Arguments
dir_path = image_folder
ext = args['extension']
output = image_folder + args['output']
fps = float(args['framepersec'])

start = dt.time(hour=0, minute=45, second=0)
end   = dt.time(hour=1, minute=20, second=0)

def TransformFrame(frame, time = None):
    '''Transform a frame to fit the desired appearence of the video.
    '''
    ###Transform the black-and white to blue
    frame = np.interp(frame, (1, 3), (0, 255))
    frame[:, :, 1] = frame[:, :, 0] / 2.
    frame[:, :, 2] = frame[:, :, 0] / 2.

    ### For fits images
    # frame = np.interp(frame, np.percentile(frame, (10, 90)), (0, 255))

    ### For Kilpisjarvi all sky images. Find the pixel corresponding to the obs direction and draw a red cross
    x, y = GetPixelFromAzEl(-45.412617, 52.50683357, frame.shape)
    frame = DrawCross(x, y, frame)

    ### Add the text label to the image
    if time is not None:
        cv2.putText(frame, 'UT ' + time.strftime("%H:%M"), (int(frame.shape[0]*0.5), int(frame.shape[1]*0.9)), cv2.FONT_HERSHEY_SIMPLEX, 1., (255, 255, 255), 2, cv2.LINE_AA)
        cv2.putText(frame, '427.8 nm', (int(frame.shape[0]*0.05), int(frame.shape[1]*0.9)), cv2.FONT_HERSHEY_SIMPLEX, 1., (255, 255, 255), 2, cv2.LINE_AA)

    frame = frame.astype('uint8')
    return frame

def GetPixelFromAzEl(az, el, frame_shape):
    '''For Kilpisjarvi images, returns the pixel position corresponding to a given az/el observation coordinates'''
    height, width, channels = frame_shape
    center_height = int(height / 2)
    center_width = int(width / 2)

    radius = (90 - el) * (width / 2) / 90

    x = center_width + int(radius * np.sin(az*np.pi/180.))
    y = center_height - int(radius * np.cos(az*np.pi/180.))

    return x, y

def DrawCross(center_x, center_y, frame, in_rad = 10, out_rad = 20):
    
    cv2.line(frame, (center_y, center_x + in_rad), (center_y, center_x + out_rad), (0, 0, 255), 3) 
    cv2.line(frame, (center_y, center_x - in_rad), (center_y, center_x - out_rad), (0, 0, 255), 3) 
    cv2.line(frame, (center_y + in_rad, center_x), (center_y + out_rad, center_x), (0, 0, 255), 3) 
    cv2.line(frame, (center_y - in_rad, center_x), (center_y - out_rad, center_x), (0, 0, 255), 3) 
    
    return frame


def ReadFrame(name):
    if name.endswith('fits'):
        with fits.open(name) as hdu_list:
            bw_image = hdu_list[0].data
            print(bw_image.shape)
        frame = np.zeros((bw_image.shape[0], bw_image.shape[1], 3))
        frame[:, :, 0] += bw_image
        frame[:, :, 1] += bw_image
        frame[:, :, 2] += bw_image
    else:
        frame = cv2.imread(image_path)

    return frame

#List all image names to be used in the movie. You can hardcode conditions on the names here.
images = []
times  = []
for f in os.listdir(dir_path):
    if f.endswith(ext) : #and ('_01' in f):
        time = dt.datetime.strptime(f.split('_')[1][:6], "%H%M%S").time()
        if start <= time <= end:
            times.append(time)
            images.append(f)

images.sort()
times.sort()

#Sort images. Hard code the sort function depending on the name formats.
# images.sort(key = lambda s: s.split('_')[2][:15]) #For PLIP fits


# Determine the width and height from the first image
image_path = os.path.join(dir_path, images[0])
frame = ReadFrame(image_path)
frame = TransformFrame(frame)
cv2.imshow('video',frame)
height, width, channels = frame.shape
print(frame.shape)


# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v') # For .mp4 format. Be sure to use lower case
# fourcc = cv2.VideoWriter_fourcc('F', 'M', 'P', '4') # For .avi format
out = cv2.VideoWriter(output, fourcc, fps, (width, height))


for image, time in zip(images, times):
    image_path = os.path.join(dir_path, image)
    frame = ReadFrame(image_path)

    # Transform each frame as defined in the function
    frame = TransformFrame(frame, time = time)

    out.write(frame) # Write out frame to video

    cv2.imshow('video', frame)
    if (cv2.waitKey(1) & 0xFF) == ord('q'): # Hit `q` to exit
        break

# Release everything if job is finished
out.release()
cv2.destroyAllWindows()

print("The output video is {}".format(output))
