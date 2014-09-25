import numpy as np
import time
import cv2

frame = cv2.imread("/home/vinay/papers/image-proc-opt/images/baboon.jpg")

# OpenCV's Gray
# =============
gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) #convert to grayscale
#gray = np.float32(gray)
cv2.imshow('frame', gray)
