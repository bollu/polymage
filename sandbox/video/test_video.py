import numpy as np
import time
import cv2

#Ahem I have this file do you! 
cap = cv2.VideoCapture("/home/ravi/Videos/Oblivion.2013.1080p.BluRay.x264.POISON.mp4")

frames = 0
startTime = time.clock()

while(cap.isOpened()):
    ret, frame = cap.read()

    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) #convert to grayscale
    #cv2.imshow('frame', gray)
 
    blur = frame 
    for i in xrange(0, 3):
        blur = cv2.blur(blur, (5, 5))
    cv2.imshow('frame', blur)

    #cv2.imshow('frame', frame)
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break
    frames += 1
      
    if frames % 10 == 0:
        currTime = time.clock()
        numsecs = currTime - startTime
        fps = 10 / numsecs
        startTime = currTime
        print "average FPS:", fps

cap.release()
cv2.destroyAllWindows()
