#############################################################
# ApplPy Software 2012 Matthew Robinson, Matthew Jackiewicz #
# Version 0.5, last updated 01 April 2012                   #
#############################################################

"""
Distribution Subclass Module

Defines commonly used distributions as subclasses of the
    RV class

"""

from rv import *

"""
Distributions included:

ArcSinRV(), BetaRV(alpha,beta), CauchyRV(a,alpha), ChiSquareSRV(N),
ErlangRV(theta,N), ExponentialRV(theta), GammaRV(theta,kappa), NormalRV(mu,sigma),
TriangularRV(a,b,c), UniformRV(a,b)

"""

def list_dist():
    print('ArcSinRV(), BetaRV(alpha,beta), CauchyRV(a,alpha), ChiSquareSRV(N),')
    print('ErlangRV(theta,N), ExponentialRV(theta), GammaRV(theta,kappa),')
    print('NormalRV(mu,sigma),TriangularRV(a,b,c), UniformRV(a,b)') 

def param_check(param):
    for i in range(len(param)):
        if type(param[i])!=int and type(param[i])!=float:
            return False
        else:
            return True

class ArcSinRV(RV):
    def __init__(self):
        X_dummy=RV([1/(pi*sqrt(x*(1-x)))],[0,1])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype

    def variate(self,n=1,s='sim'):
        # Generate arc sin variates
        idf_func=(1/2)-((1/2)*cos(pi*t))
        varlist=[]
        for i in range(n):
            val=random()
            var=idf_func.subs(t,val).evalf()
            varlist.append(var)
        varlist.sort()
        return varlist


class BetaRV(RV):
    def __init__(self,alpha=Symbol('alpha'),beta=Symbol('beta')):
        X_dummy=RV((gamma(alpha+beta)*(x**(alpha-1))*(1-x)**(beta-1))/
               (gamma(alpha)*gamma(beta)),[0,1])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[alpha,beta]

class CauchyRV(RV):
    def __init__(self,a=Symbol('a'),alpha=Symbol('alpha')):
        X_dummy=RV((1)/(alpha*pi*(1+((x-alpha)**2/alpha**2))),[-oo,oo])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[a,alpha]

    def variate(self,n=1,s='sim'):
        # If no parameter is specified, return an error
        if param_check(self.parameter)==False:
            raise RVError('Not all parameters specified')

        # Generate exponential variates
        idf_func=self.parameter[0]-cot(pi*t)*self.parameter[1]
        varlist=[]
        for i in range(n):
            if s=='sim':
                val=random()
            else:
                val=s
            var=idf_func.subs(t,val).evalf()
            varlist.append(var)
        varlist.sort()
        return varlist

class ChiSquareRV(RV):
    def __init__(self,N=Symbol('N')):
        X_dummy=RV((x**(N-1)*exp(-(x**2/2)))/
                   (2**((N/2)-1)*gamma(N/2)),[0,oo])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[N]

class ErlangRV(RV):
    def __init__(self,N=Symbol('N'),theta=Symbol('theta')):
        X_dummy=RV((N*(N*x)**(theta-1)*exp(-N*x))/(factorial(theta-1)),[0,oo])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[N,theta]

class ExponentialRV(RV):
    def __init__(self,theta=Symbol('theta')):
        X_dummy=RV([theta*exp(-theta*x)],[0,oo])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[theta]

    def variate(self,n=1,s='sim'):
        # If no parameter is specified, return an error
        if param_check(self.parameter)==False:
            raise RVError('Not all parameters specified')

        # Generate exponential variates
        idf_func=-(ln(1-t))/(self.parameter[0])
        varlist=[]
        for i in range(n):
            if s=='sim':
                val=random()
            else:
                val=s
            var=idf_func.subs(t,val)
            varlist.append(var)
        varlist.sort()
        return varlist

class GammaRV(RV):
    def __init__(self,theta=Symbol('theta'),kappa=Symbol('kappa')):
        X_dummy=RV((theta*(theta*x)**(kappa-1)*exp(-theta*x))/(gamma(kappa)),
                   [0,oo])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[theta,kappa]

class NormalRV(RV):
    def __init__(self,mu=Symbol('mu'),sigma=Symbol('sigma')):
        X_dummy=RV((exp((-(x-mu)**2)/(2*sigma**2))*sqrt(2))/(2*sigma*sqrt(pi))
                   ,[-oo,oo])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[mu,sigma]


class TriangularRV(RV):
    def __init__(self,a=Symbol('a'),b=Symbol('b'),c=Symbol('c')):
        X_dummy=RV([(2*(x-a))/((c-a)*(b-a)),(2*(c-x))/((c-a)*(c-b))],[a,b,c])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[a,b,c]

class UniformRV(RV):
    def __init__(self,a=Symbol('a'),b=Symbol('b')):
        X_dummy=RV((b-a)**(-1),[a,b])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[a,b]

    def variate(self,n=1,sim='sim'):
        # If no parameter is specified, return an error
        if param_check(self.parameter)==False:
            raise RVError('Not all parameters specified')

        # Generate uniform variates
        idf_func=-t*self.parameter[0]+t*self.parameter[1]+self.parameter[0]
        varlist=[]
        for i in range(n):
            if s=='sim':
                val=random()
            else:
                val=s
            var=idf_func.subs(t,val)
            varlist.append(var)
        varlist.sort()
        return varlist

class WeibullRV(RV):   
    def __init__(self,theta=Symbol('theta'),kappa=Symbol('kappa')):
        X_dummy=RV(kappa*theta**(kappa)*x**(kappa-1)*exp(-(theta*x)**kappa),
                   [0,oo])
        self.func=X_dummy.func
        self.support=X_dummy.support
        self.ftype=X_dummy.ftype
        self.parameter=[theta,kappa]


    def variate(self,n=1,s='sim'):
        # If no parameter is specified, return an error
        if param_check(self.parameter)==False:
            raise RVError('Not all parameters specified')

        # Generate weibull variates
        idf_func=exp(-(-ln(-ln(1-t))+self.parameter[1]*ln(self.parameter[0]))/
                     self.parameter[1])
        varlist=[]
        for i in range(n):
            if s=='sim':
                val=random()
            else:
                val=s
            var=idf_func.subs(t,val)
            varlist.append(var)
        varlist.sort()
        return varlist













        

