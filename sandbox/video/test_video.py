import ctypes
import numpy as np
import time
import cv2

# load our shared library
libgray = ctypes.cdll.LoadLibrary("./gray.so")
# nickname our library gray function
gray = libgray.gray


#Ahem I have this file do you! 
cap = cv2.VideoCapture("/big-data/Interstellar-TLR-3-51ch-1080p-HDTN.mp4")


frames = 0
#startTime = time.clock()
startTime = time.time()

while(cap.isOpened()):
    ret, frame = cap.read()

    # OpenCV's Gray
    # =============
    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) #convert to grayscale
    #cv2.imshow('frame', gray)

    # Our Gray
    # ========
    rows = frame.shape[0]
    cols = frame.shape[1]
    blank = np.zeros((rows, cols), frame.dtype)
    # Shared library function call here
    gray(ctypes.c_int(cols), ctypes.c_int(rows), ctypes.c_void_p(frame.ctypes.data), ctypes.c_void_p(blank.ctypes.data))
    cv2.imshow('frame', blank)
 

    # OpenCV's bluuuurrrr
    # ===================
    #blur = frame
    #for i in xrange(0, 3):
    #    blur = cv2.blur(blur, (5, 5))
    #cv2.imshow('frame', blur)


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
