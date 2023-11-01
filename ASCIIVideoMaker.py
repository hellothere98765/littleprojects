import math
import numpy as np
from moviepy.editor import *
import pyautogui
import os
from numpy import asarray
import time
import subprocess
ascii="$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\|()1{}[]?-_+~<>i!lI;:,\"`'.  "
#ascii="@%#*+=-:.  "
#Uses Pillow version 9.5.0
Width=450
path = input("What's the file? ")
writeto=input("What's the name of the file you want this to be written to? ")
clip = VideoFileClip(path)
w=2*clip.w*130/clip.h
h=130
clip = clip.resize((w,h))
clip=clip.subclip(0,20)
n=1/clip.fps
print(clip.size)
a=open("a.txt","w").close()
def make_frame(t):
    a=open("a.txt","r+")
    i=clip.get_frame(t)
    for x in range(0,len(i)):
        for s in range(0,len(i[x])):
            a.write(ascii[math.floor((len(ascii)-1)*(i[x][s][0]*.21+.72*i[x][s][1]+.07*i[x][s][2])/255)])
        a.write("\n")
    a.close()
    os.popen("a.txt","r")
    time.sleep(.7)
    imageder = pyautogui.screenshot(region=(0,50,3*w,6*h))
    time.sleep(.7)
    subprocess.call("TASKKILL /F /IM notepad.exe", shell=True)
    return asarray(imageder)
Concater=VideoClip(make_frame,duration=clip.duration)
Concater = Concater.set_audio(clip.audio)
Concater.write_videofile(writeto,clip.fps)
