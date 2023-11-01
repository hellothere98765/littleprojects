import numpy
import PIL
from moviepy.editor import *

path = input("What's the file? ")
writeto=input("What's the name of the file you want this to be written to? ")
comp=int(input("How compressed do you want it as a percent?"))
a=float(input("What time do you want it to start at?"))
b=float(input("What time do you want it to end at?"))
clip = VideoFileClip(path)
clip=clip.subclip(a,b)
def fftcompress(r):
    x=numpy.fft.fft2(r)# Originally: x=numpy.fft.fftshift(numpy.fft.fft2(r))
    a=[[numpy.abs(x[i][j]) for j in range(0,len(x[1]))] for i in range(0,len(x))]
    s=[a[int(i/len(a[1]))][int(len(a[1])*(i/len(a[1])-int(i/len(a[1]))))] for i in range(0,len(a)*len(a[1]))]
    s.sort()
    magn=s[int(len(s)*comp/100)]
    g=[[x[i][j] if a[i][j]>magn else 0 for j in range(0,len(x[1]))] for i in range(0,len(x))]
    return(numpy.clip((numpy.fft.ifft2(g)),0,255))#Originally: return(numpy.clip((numpy.fft.ifft2(numpy.fft.ifftshift(g))),0,255))
    #[int(widthgap):int((clip.w-widthgap))][int(heightgap):int((clip.h-heightgap))]
def make_frame(t):
    i=clip.get_frame(t)
    red_array=i[:,:,0]
    green_array=i[:,:,1]
    blue_array=i[:,:,2]
    return numpy.stack((fftcompress(red_array),fftcompress(green_array),fftcompress(blue_array)),axis=-1)

    #x=fftcompress(red_array)
    #y=fftcompress(green_array)
    #z=fftcompress(blue_array)

video=VideoClip(make_frame,duration=clip.duration)
video = video.set_audio(clip.audio)
video.write_videofile(writeto,clip.fps)
