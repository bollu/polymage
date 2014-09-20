#!/usr/bin/env python

'''
Multithreaded video processing sample.
Usage:
   video_threaded.py {<video device number>|<video file name>}

   Shows how python threading capabilities can be used
   to organize parallel captured frame processing pipeline
   for smoother playback.

Keyboard shortcuts:

   ESC - exit
   space - switch between multi and single threaded processing
'''


import ctypes
import numpy as np
import cv2

from multiprocessing.pool import ThreadPool
from collections import deque

from common import clock, draw_str, StatValue
import video

# load polymage shared libraries
libmask = ctypes.cdll.LoadLibrary("./unsharp.so")
libharris = ctypes.cdll.LoadLibrary("./harris.so")

unsharp_par = libmask.unsharp_mask_par
unsharp = libmask.unsharp_mask
harris = libharris.pipeline_harris


class DummyTask:
    def __init__(self, data):
        self.data = data
    def ready(self):
        return True
    def get(self):
        return self.data

if __name__ == '__main__':
    import sys

    print __doc__

    try: fn = sys.argv[1]
    except: fn = 0
    cap = video.create_capture(fn)


    def process_frame(frame, t0):
        # some intensive computation...
        frame = cv2.medianBlur(frame, 19)
        frame = cv2.medianBlur(frame, 19)
        frame = cv2.medianBlur(frame, 19)
        return frame, t0

    def unsharp_mask(frame, t0):
        #GaussianBlur(img_region, blury, Size(5, 5), sigma, sigma);
        #mask = (abs(img_region - blury) < threshold);
        #ocv_sharpened = img_region*(1+weight) + blury*(-weight);
        w = 1
        blur = cv2.GaussianBlur(frame ,(5,5), -0.75)
        sharp = frame*(1 + w) - blur * w
        return frame, t0

    def harris_cv(frame, t0):
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        gray = np.float32(gray)
        coarsity = cv2.cornerHarris(gray, 3, 3, 0.04)
        return coarsity, t0

    def harris_polymage(frame, t0):
        rows = frame.shape[0]
        cols = frame.shape[1]
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        gray = np.float32(gray)
        res = np.zeros((rows, cols), np.float32) 
        harris(ctypes.c_int(cols-2), ctypes.c_int(rows-2), \
               ctypes.c_void_p(gray.ctypes.data), ctypes.c_void_p(res.ctypes.data))
        return res, t

    threadn = cv2.getNumberOfCPUs()
    pool = ThreadPool(processes = threadn)
    pending = deque()

    threaded_mode = True

    latency = StatValue()
    frame_interval = StatValue()
    last_frame_time = clock()
    while True:
        while len(pending) > 0 and pending[0].ready():
            res, t0 = pending.popleft().get()
            latency.update(clock() - t0)
            draw_str(res, (40, 40), "threaded      :  " + str(threaded_mode))
            draw_str(res, (40, 80), "latency        :  %.1f ms" % (latency.value*1000))
            draw_str(res, (40, 120), "frame interval :  %.1f ms" % (frame_interval.value*1000))
            cv2.imshow('threaded video', res)
        if len(pending) < threadn:
            ret, frame = cap.read()
            t = clock()
            frame_interval.update(t - last_frame_time)
            last_frame_time = t
            if threaded_mode:
                #task = pool.apply_async(process_frame, (frame.copy(), t))
                #inframe = frame.transpose(2, 0, 1).copy()
                inframe = frame.copy()
                #task = pool.apply_async(unsharp_mask, (inframe, t))
                task = pool.apply_async(harris_cv, (inframe, t))
                #task = pool.apply_async(harris_polymage, (inframe, t))
            else:
                task = DummyTask(process_frame(frame, t))
            pending.append(task)
        ch = 0xFF & cv2.waitKey(1)
        if ch == ord(' '):
            threaded_mode = not threaded_mode
        if ch == 27:
            break
cv2.destroyAllWindows()
