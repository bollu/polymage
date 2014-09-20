import ctypes
import numpy as np
import time
import cv2

'''
libgray = ctypes.cdll.LoadLibrary("./gray.so")
gray = libgray.gray
'''

# load our shared library
libmask = ctypes.cdll.LoadLibrary("./unsharp.so")
# nickname our library gray function
unsharp = libmask.unsharp_mask


#Ahem I have this file do you! 
#cap = cv2.VideoCapture("/big-data/Interstellar-TLR-3-51ch-1080p-HDTN.mp4")
cap = cv2.VideoCapture("/home/ravi/Videos/Interstellar-TLR-3-51ch-4K-HDTN.mp4")


frames = 0
#startTime = time.clock()
startTime = time.time()

while(cap.isOpened()):
    ret, frame = cap.read()

    # OpenCV's Gray
    # =============
    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) #convert to grayscale
    #gray = np.float32(gray)
    #cv2.imshow('frame', gray)
    #print gray.dtype

    '''
    # Our Gray
    # ========
    rows = frame.shape[0]
    cols = frame.shape[1]
    blank = np.zeros((rows, cols), frame.dtype)
    # Shared library function call here
    gray(ctypes.c_int(cols), ctypes.c_int(rows), ctypes.c_void_p(frame.ctypes.data), ctypes.c_void_p(blank.ctypes.data))
    cv2.imshow('frame', blank)
    '''
 

    # OpenCV's bluuuurrrr
    # ===================
    #blur = frame
    #for i in xrange(0, 3):
    #    blur = cv2.blur(blur, (5, 5))
    #cv2.imshow('frame', blur)

    # Our Unsharp Mask
    # ========
    rows = frame.shape[0]
    cols = frame.shape[1]
    # Shared library function call here
    frame = np.float32(frame)
    inframe = frame.transpose(2, 0, 1).copy()
    blank = np.zeros((3, rows-4, cols-4), np.float32)
    unsharp(ctypes.c_int(cols-4), ctypes.c_int(rows-4), ctypes.c_float(0.01), ctypes.c_float(1), ctypes.c_void_p(inframe.ctypes.data), ctypes.c_void_p(blank.ctypes.data))
    #frame = frame.transpose(1, 2, 0)
    out = blank.transpose(1, 2, 0).copy()
    #cv2.imshow('frame', frame/255)
    cv2.imshow('frame1', out/255)
    cv2.imshow('frame2', frame/255)
 

    if cv2.waitKey(1) & 0xFF == ord('q'):
        break
    frames += 1
      
    if frames % 10 == 0:
        #currTime = time.clock()
        currTime = time.time()
        numsecs = currTime - startTime
        fps = 10 / numsecs
        startTime = currTime
        print "average FPS:", fps

cap.release()
cv2.destroyAllWindows()
