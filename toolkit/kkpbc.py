#!/usr/bin/python2.7

import numpy
import copy

class PBC(object):
    def __init__(self, x, y, z, a=90.0, b=90.0, c=90.0):
        self.L = numpy.array([x,y,z])
        self.angle = numpy.array([a,b,c])
        self.L_half = self.L * 0.5
        self.L_inv = 1.0 / self.L 
        self.origin = numpy.array([0.0, 0.0, 0.0])
        self.set_center(numpy.array([0.0, 0.0, 0.0]))
    def diff_crd_minim_image(self, crd1, crd2):
        """
        calculate deviation between the two coordinates
        crd1 and crd2 must be numpy.array 
        """
        d = (crd1 - crd2) 
        d = d - self.L * numpy.array([int(x) for x in d * self.L_inv])
        return d
    def info_txt(self):
        info = ""
        info += str(self.x)  + ", "
        info += str(self.y)  + ", "
        info += str(self.z)  + ", "
        info += str(self.a)  + ", "
        info += str(self.b)  + ", "
        info += str(self.c)  + "\N"
        return info
    def set_center(self, center):
        self.origin = center - (self.L_half)
