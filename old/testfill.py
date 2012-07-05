#!/usr/bin/env python

from pylab import *
x1 = arange (0 , 2 , 0.01)
y1 = sin(2*pi*x1)
y2 = sin(4*pi*x1) + 2

x = concatenate((x1,x1[::-1]))
y = concatenate((y1,y2[::-1]))

p = fill(x,y,facecolor='g')
hold(True)
fill(array([0,0,1,1,0]),array([1,0,0,-1,1]),'r',alpha=0.9)
hold(False)

show()
