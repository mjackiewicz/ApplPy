######################################################################
# ApplPy Software 2012 Matthew Robinson, Matthew Jackiewicz          #
# Version 0.5, last updated 13 May 2012                              #
######################################################################

"""
Plotting Module

Defines procedures for plotting random variables

"""
from sympy import *
from applpy import *

try:
    from pylab import *
except:
    print 'WARNING: Plotting not currently enabled'
    print 'Download matplotlib to enable plotting.'
    print ''

def mat_plot(funclist,suplist,lab1=None,lab2=None,ftype='continuous'):
    """
    Procedure Name: mat_plot
    Purpose: Create a matplotlib plot of a random variable
    Arguments:  1. RVar: A random variable
                2. suplist: The support of the plot
    Output:     1. A plot of the random variable
    """
    # if the random variable is continuous, plot the function
    if ftype=='continuous':
        for i in range(len(funclist)):
            x=arange(suplist[i],suplist[i+1],0.01)
            s=eval(funclist[i])
            plot(x,s,linewidth=1.0,color='green')
        if lab1=='idf':
            xlabel('s')
        else:
            xlabel('x')
        if lab1!=None:
            ylabel(lab1)
        if lab2!=None:
            title(lab2)
        grid(True)
    # If the random variable is discrete, plot the function
    if ftype=='discrete':
        plot(suplist,funclist,'ro')
        if lab1=='F-1(s)':
            xlabel('s')
        else:
            xlabel('x')
        if lab1!=None:
            ylabel(lab1)
        if lab2!=None:
            title(lab2)
        grid(True)


