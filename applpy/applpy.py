#############################################################
# ApplPy Software 2012 Matthew Robinson, Matthew Jackiewicz #
# Version 0.5, last updated 30 March 2012                   #
#############################################################

"""
Main ApplPy Module

Imports supporting ApplPy Modules

"""
from __future__ import division
from rv import *
from dist_type import *
from stats import *

def Menu():
    print 'ApplPy Procedures'
    print ""
    print 'Procedure Notation'
    print ""
    print 'Capital letters are random variables'
    print 'Lower case letters are number'
    print 'Greek letters are parameters'
    print 'gX indicates a function'
    print 'n and r are positive integers where n>=r'
    print 'Square brackets [] denote a list'
    print 'Curly bracks {} denote an optional variable'
    print ""
    print ""

    print 'RV Class Procedures'
    print 'X.variate(n,x),X.verifyPDF()'
    print ""

    print 'Functional Form Conversion'
    print 'CDF(X,{x}),CHF(X,{x}),HF(X,{x}),IDF(X,{x})'
    print 'PDF(X,{x}),SF(X,{x}),BootstrapRV([data])'
    print 'Convert(X,{x})'
    print ""    

    print 'Procedures on One Random Variable'
    print 'ConvolutionIID(X,n),CoefOfVar(X),ExpectedValue(X,gx)'
    print 'Kurtosis(X),MaximumIID(X,n),Mean(X),MGF(X)'
    print 'MinimumIID(X,n),OrderStat(X,n,r),ProductIID(X,n)'
    print 'Skewness(X),Transform(X,gX),Truncate(X,[x1,x2])'
    print 'Variance(X)'
    print ""

    print 'Procedures on Two Random Variables'
    print 'Convolution(X,Y),Maximum(X,Y),Minimum(X,Y)'
    print 'Mixture([p1,p2],[X,Y]),Product(X,Y)'
    print ""

    print 'Utilities'
    print 'PlotDist(X,{[x1,x2]}),PlotDisplay([plotlist],{[x1,x2]})'
    print ""

    print 'Continuous Distributions'
    print 'BetaRV(alpha,beta),CauchyRV(a,alpha),ChiRV(N),ChiSquareRV(N)'
    print 'ErlangRV(theta,N),ExponentialRV(theta)'
    print 'ExponentialPowerRV(theta,kappa),ExtremeValueRV(alpha,beta)'
    print 'GammaRV(theta,kappa),GompertzRV(theta,kappa)'
    print 'InverseGaussianRV(theta,mu),InverseGammaRV(alpha,beta)'
    print 'LogGammaRV(alpha,beta), LogisticRV(kappa,theta)'
    print 'LogLogisticRV(theta,kappa),LogNormalRV(mu,sigma)'
    print 'LomaxRV(kappa,theta),MuthRV(kappa),NormalRV(mu,sigma)'
    print 'ParetoRV(theta,kappa),RayleighRV(theta),TriangularRV(a,b,c)'
    print 'TRV(N),UniformRV(a,b),WeibullRV(theta,kappa)'
    print ""

    print 'Discrete Distributions'
    print 'BenfordRV(),BinomialRV(n,p),GeometricRV(p),PoissonRV(theta)'

print '-----------------'
print 'Welcome to ApplPy'
print '-----------------'
Menu()

