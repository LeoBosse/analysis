#!/usr/bin/python3
# -*-coding:utf-8 -*

################################################################################
# This script creates a movie from a set of images.
# It loads the images from a folder (with a filter on the names), transform them, and stack them in a movie.
################################################################################

import cv2
import argparse
import os
from vendange_configuration import Config
import numpy as np

config = Config()

image_folder = config.data_path + "allsky_camera/KIL20200224/"

# Construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-ext", "--extension", required=False, default='png', help="extension name. default is 'png'.")
ap.add_argument("-o", "--output", required=False, default='output.mp4', help="output video file")
args = vars(ap.parse_args())

# Arguments
dir_path = f'{config.data_path}allsky_camera/KIL20200224/'
ext = args['extension']
output = f'{config.data_path}allsky_camera/KIL20200224/' + args['output']



def TransformFrame(frame):
    '''Transform a frame to fit the desired appearence of the video.
    Here it transform the black-and white to blue, and draw a red cross where ptcu was looking.
    '''
    frame = np.interp(frame, (1, 3), (0, 255))
    frame[:, :, 1] = 0
    frame[:, :, 2] = 0

    x, y = GetPixelFromAzEl(-45.412617, 52.50683357, frame.shape)

    frame = DrawCross(x, y, frame)

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

def DrawCross(center_x, center_y, frame):
    for i in range(10):
        frame[center_y, center_x + i + 10] = 0, 0, 255
        frame[center_y, center_x - i - 10] = 0, 0, 255
        frame[center_y + i + 10, center_x] = 0, 0, 255
        frame[center_y - i - 10, center_x] = 0, 0, 255
    return frame


#Load and filter the images
images = []
for f in os.listdir(dir_path):
    if f.endswith(ext) and ('25_01' in f):
        images.append(f)
images.sort()


# Determine the width and height from the first image
image_path = os.path.join(dir_path, images[0])
frame = cv2.imread(image_path)
frame = TransformFrame(frame)
cv2.imshow('video',frame)
height, width, channels = frame.shape
print(frame.shape)


# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v') # Be sure to use lower case
out = cv2.VideoWriter(output, fourcc, 20.0, (width, height))


for image in images:
    image_path = os.path.join(dir_path, image)
    frame = cv2.imread(image_path)

    # Transform each frame as defined in the function
    frame = TransformFrame(frame)

    out.write(frame) # Write out frame to video

    cv2.imshow('video', frame)
    if (cv2.waitKey(1) & 0xFF) == ord('q'): # Hit `q` to exit
        break

# Release everything if job is finished
out.release()
cv2.destroyAllWindows()

print("The output video is {}".format(output))
