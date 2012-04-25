######################################################################
# ApplPy Software 2012 Matthew Robinson, Matthew Jackiewicz          #
# Version 0.5, last updated 18 April 2012                            #
######################################################################

"""
Main Random Variable Module

Defines the random variable class
Defines procedures for chaning functional form
Defines procedures on one random variable
Defines procudures on two random variables

"""

from __future__ import division
from sympy import *
import plot as plt
from random import random
x,y,z,t=symbols('x y z t')

class RVError(Exception):
    """
    RVError Class
    Defines a custom error message for exceptions relating
    to the random variable class
    """
    def __init__(self,value):
        self.value=value
    def __str__(self):
        return repr(self.value)

class RV:
    """
    RV Class
    Defines the data structure of ApplPy random variables
    Defines procedures relating to ApplPy random variables
    """

    def __init__(self,func,support,ftype=['continuous','pdf']):
        """
        Creates an instance of the random variable class
            The random variable default is to produce a continuous pdf
        Checks the random variable for errors
        """

        # Check for errors in the data structure of the random
        #   variable

        # Check to make sure that the given function is in the
        #   form of a list
        # If it is not in the form of a list, place it in a list
        if isinstance(func,list)!=True:
            func1=func
            func=[func1]
        # Check to make sure that the given support is in the form of
        #   a list
        if isinstance(support,list)!=True:
            raise RVError('Support must be a list')
        # Check to make sure that the support list has the correct
        #   length
        # The support list should be one element larger than the
        #   function list for continuous distributions, and the same
        #   size for discrete
        if ftype[0]=='continuous':
            if len(support)-len(func)!=1:
                raise RVError('Support has incorrect number of elements')
        if ftype[0]=='discrete':
            if len(support)-len(func)!=0:
                raise RVError('Support has incorrect number of elements')
        # Check to make sure that the elements of the support list are
        #   in ascending order
        for i in range(len(support)-1):
            if support[i]>support[i+1]:
                raise RVError('Support is not in ascending order')
        # Check to make sure that the random variable is either
        #   discrete or continuous
        if ftype[0] not in ['continuous','discrete']:
            raise RVError('Random variable must either be discrete or continuous')
        # Initialize the random variable
        self.func=func
        self.support=support
        self.ftype=ftype

    """
    Special Class Methods

    Procedures:
        1. display(self)
        2. __repr__(self)
        3. __len__(self)
        4. __add__(self,other)
    """

    def display(self):
        """
        Creates a default print setting for the random variable class
        """
        if self.ftype[0]=='continuous':
            print '%s %s with support %s:'%(self.ftype[0],self.ftype[1],self.support)
            return self.func
        if self.ftype[0]=='discrete':
            print '%s %s where {x->f(x)}:'%(self.ftype[0],self.ftype[1])
            for i in range(len(self.support)):
                if i!=(len(self.support)-1):
                    print '{%s -> %s}, '%(self.support[i],self.func[i]),
                else:
                    print '{%s -> %s}'%(self.support[i],self.func[i])

    def __repr__(self):
        """
        Sets the default string display setting for the random variable class
        """
        return repr(self.display())

    def __len__(self):
        """
        Sets the behavior for the len() procedure when an instance of the
            random variable class is given as input
        """
        return len(self.func)

    # Set the behavior for the operators '+,-,*,/'

    def __add__(self,other):
        """
        Sets the behavior of the '+' operator
        """
        return Convolution(self,other)

    def __sub__(self,other):
        """
        Sets the behavior of the '-' operator
        """
        gX=[[-x],[-oo,oo]]
        RVar=Transform(other,gX)
        return Convolution(self,RVar)

    def __mul__(self,other):
        """
        Sets the behavior of the '*' operator
        """
        return Product(self,other)

    def __truediv__(self,other):
        """
        Sets the behavior of the '/' operator
        """
        gX=[[1/x],[-oo,oo]]
        RVar=Transform(other,gX)
        return Product(self,RVar)

        # To add later:
        #   1. self.__getitem__() / self.__delitem__ for use in dictionaries
        #   2. self.__cmp__() to allow for comparison of two random variables

    """
    Utility Methods

    Procedures:
        1. verifyPDF(self)
        2. variate(self,n)
    """
    def verifyPDF(self):
        """
        Checks whether of not the pdf of a random variable is valid
        """
        if self.ftype[0]=='continuous':
            # Convert the random variable to PDF form
            X_dummy=PDF(self)
            # Check to ensure that the area under the PDF is 1
            print 'Now checking for area...'
            area=0
            for i in range(len(X_dummy.func)):
                val=integrate(X_dummy.func[i],(x,X_dummy.support[i],
                                               X_dummy.support[i+1]))
                area+=val
            print 'The area under f(x) is: %s'%(area)
            
            # Check to ensure that the PDF is positive over the entire support of the PDF
            #print 'Now checking for absolute value...'

            # Need an algorithm to check the function for non-negative
            # integrate command does not seem to work for functions like Abs(x**2)
            # finding max of function may be best bet, no max function
            # Could use numerical maximization procedure coded in APPLBayes for Maple            

            print 'The pdf of the random variable:'
            print '%s'%(X_dummy.func)
            print 'continuous pdf with support %s'%(X_dummy.support)
            if area>.9999 and area<1.00001:
                print 'is valid'
            else:
                print 'is not valid'
        if self.ftype[0]=='discrete':
            # Convert the random variable to PDF form
            X_dummy=PDF(self)
            # Check to ensure that the area under the PDF is 1
            print 'Now checking for area...'
            area=0
            for i in range(len(self.support)):
                area+=self.func[i]
            print 'The area under f(x) is: %s'%(area)
            print 'The pdf of the random variable'
            if area>.9999 and area<1.0001:
                print 'is valid'
            else:
                print 'is not valid'

    def variate(self,n=1,s='sim',sensitivity=.00001):
        """
        Generates a list of n random variates from the random variable
            using the Newton-Raphson Method
        """   
        # Find the cdf and pdf functions (to avoid integrating for each
        #   variate
        cdf=CDF(self)
        pdf=PDF(self)
        mean=Mean(self)
        # Create a list of variates
        varlist=[]
        for i in range(n):
            guess=mean
            if s=='sim':
                val=random()
            else:
                val=s
            for i in range(10):
                try:
                    if len(self.func)==1:
                        guess=(guess-((cdf.func[0].subs(x,guess)-val)/
                                     pdf.func[0].subs(x,guess))).evalf()
                    else:
                        guess=(guess-((CDF(cdf,guess)-val)/PDF(pdf,guess))).evalf()
                except:
                    if guess>self.support[len(self.support)-1]:
                        cfunc=cdf.func[len(self.func)-1].subs(x,guess)
                        pfunc=pdf.func[len(self.func)-1].subs(x,guess)
                        guess=(guess-((cfunc-val)/pfunc)).evalf()
                    if guess<self.support[0]:
                        cfunc=cdf.func[0].subs(x,guess)
                        pfunc=pdf.func[0].subs(x,guess)
                        guess=(guess-((cfunc-val)/pfunc)).evalf()
            varlist.append(guess)            
        varlist.sort()
        return varlist


"""
Procedures for converting functional form

Procedures:
    1. CDF(RVar,value)
    2. CHF(RVar,value)
    3. HF(RVar,value)
    4. IDF(RVar,value)
    5. PDF(RVar,value)
    6. SF(RVar,value)
    7. BootstrapRV(varlist)
"""

def check_value(value,sup):
    # Not intended for use by end user
    """
    Procedure Name: check_value
    Purpose: Check to see if a value passed to CDF,CHF,HF,PDF or SF is in the
                support of the random variable
    Arguments:  1. value: The value passed to RV procedure
                2. sup: The support of the RV in the procedure
    Output:     1. True if the value given is within the support
                2. False otherwise
    """
    if value==x:
        return True
    else:
        max_idx=len(sup)-1
        if value<sup[0] or value>sup[max_idx]:
            return False
        else:
            return True

def CDF(RVar,value=x):
    """
    Procedure Name: CDF
    Purpose: Compute the cdf of a random variable
    Arguments:  1. RVar: A random variable
                2. value: An integer or floating point number (optional)
    Output:     1. CDF of a random variable (if value not specified)
                2. Value of the CDF at a given point (if value is specified)
    """

    # Check to make sure the value given is within the random variable's support
    if check_value(value,RVar.support)!=True:
        raise RVError('Value is not within the support of the random variable')

    # If the distribution is continous, find and return the distribution
    #   of the random variable
    if RVar.ftype[0]=='continuous':
        # If the random variable is already a cdf, nothing needs to be done
        if RVar.ftype[1]=='cdf':
            if value==x:
                return RVar
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        cdfvalue=RVar.func[i].subs(x,value)
                        return simplify(cdfvalue)
        # If the random variable is a sf, find and return the cdf of the random
        #   variable
        if RVar.ftype[0]=='sf':
            X_dummy=SF(RVar)
            # Compute the sf for each segment
            cdflist=[]
            for i in range(len(X_dummy.func)):
                cdflist.append(1-X_dummy.funcCD[i])
            # If no value is specified, return the sf function
            if value==x:
                return RV(cdflist,X_dummy.support,['continuous','cdf'])
            # If not, return the value of the cdf at the specified value
            else:
                for i in range(len(X_dummy.support)):
                    if value>=X_dummy.support[i] and value<=X_dummy.support[i+1]:
                        cdfvalue=cdflist[i].subs(x,value)
                        return simplify(cdfvalue)
        # If the random variable is not a cdf or sf, compute the pdf of the random
        #   variable, and then compute the cdf by integrating over each
        #   segment of the random variable
        else:
            X_dummy=PDF(RVar)
            # Substitue the dummy variable 't' into the dummy rv
            funclist=[]
            for i in range(len(X_dummy.func)):
                newfunc=X_dummy.func[i].subs(x,t)
                funclist.append(newfunc)
            # Integrate to find the cdf
            cdflist=[]
            for i in range(len(funclist)):
                cdffunc=integrate(funclist[i],(t,X_dummy.support[i],x))
                # Adjust the constant of integration
                if i!=0:
                    const=cdflist[i-1].subs(x,X_dummy.support[i])-cdffunc.subs(x,X_dummy.support[i])
                    cdffunc=cdffunc+const
                if i==0:
                    const=0-cdffunc.subs(x,X_dummy.support[i])
                    cdffunc=cdffunc+const
                cdflist.append(cdffunc)
            # If no value is specified, return the cdf
            if value==x:
                return RV(cdflist,X_dummy.support,['continuous','cdf'])
            # If a value is specified, return the value of the cdf
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        cdfvalue=cdflist[i].subs(x,value)
                        return simplify(cdfvalue)
                    
    # If the distribution is discrete, find and return the cdf of the random variable
    if RVar.ftype[0]=='discrete':
        # If the distribution is already a cdf, nothing needs to be done
        if RVar.ftype[1]=='cdf':
            if value==x:
                return RVar
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return RVar.func[RVar.support.index(value)]
        # If the distribution is a sf, find the cdf by reversing the function list
        if RVar.ftype[1] in ['sf','chf','hf']:
            X_dummy=SF(RVar)
            newfunc=[]
            for i in reversed(range(len(X_dummy.func))):
                newfunc.append(X_dummy.func[i])
            Xsf=RV(newfunc,X_dummy.support,['discrete','cdf'])
            if value==x:
                return Xsf
            if value!=x:
                if value not in Xsf.support:
                    return 0
                else:
                    return Xsf.func[Xsf.support.index(value)]
        # If the distribution is not a cdf or sf, find the pdf and then compute the cdf
        #   by summation
        else:
            X_dummy=PDF(RVar)
            cdffunc=[]
            area=0
            for i in range(len(X_dummy.support)):
                area+=X_dummy.func[i]
                cdffunc.append(area)
            if value==x:
                return RV(cdffunc,X_dummy.support,['discrete','cdf'])
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return cdffunc.func[RVar.support.index(value)]
                

def CHF(RVar,value=x):
    """
    Procedure Name: CHF
    Purpose: Compute the chf of a random variable
    Arguments:  1. RVar: A random variable
                2. value: An integer or floating point number (optional)
    Output:     1. CHF of a random variable (if value not specified)
                2. Value of the CHF at a given point (if value is specified)
    """
    
    # Check to make sure the value given is within the random variable's support
    if check_value(value,RVar.support)!=True:
        raise RVError('Value is not within the support of the random variable')

    # If the distribution is continuous, find and return the chf of the random variable
    if RVar.ftype[0]=='continuous':
        # If the distribution is already a chf, nothing needs to be done
        if RVar.ftype[1]=='chf':
            if value==x:
                return RVar
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        chfvalue=RVar.func[i].subs(x,value)
                        return simplify(chfvalue)
        # Otherwise, find and return the chf
        else:
            X_dummy=SF(RVar)
            # Generate a list of sf functions
            sflist=[]
            for i in range(len(X_dummy.func)):
                sflist.append(X_dummy.func[i])
            # Generate chf functions
            chffunc=[]
            for i in range(len(sflist)):
                newfunc=-ln(sflist[i])
                chffunc.append(simplify(newfunc))
            # If a value is not specified, return the chf of the random variable
            if value==x:
                return RV(chffunc,X_dummy.support,['continuous','chf'])
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.func[i] and value<=RVar.support[i+1]:
                        chfvalue=chffunc[i].subs(x,value)
                        return simplify(chfvalue)
                    
    # If the random variable is discrete, find and return the chf
    if RVar.ftype[0]=='discrete':
        # If the distribution is already a chf, nothing needs to be done
        if RVar.ftype[1]=='chf':
            if value==x:
                return RVar
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return RVar.func[RVar.support.index(value)]
        # Otherwise, use the survivor function to find the chf
        else:
            X_sf=SF(RVar)
            chffunc=[]
            for i in range(len(X_sf.func)):
                chffunc.append(-log(X_sf.func[i]))
            if value==x:
                return RV(chffunc,X_sf.support,['discrete','chf'])
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return chffunc[RVar.support.index(value)]
                    

def HF(RVar,value=x):
    """
    Procedure Name: HF
    Purpose: Compute the hf of a random variable
    Arguments:  1. RVar: A random variable
                2. value: An integer or floating point number (optional)
    Output:     1. HF of a random variable (if value not specified)
                2. Value of the HF at a given point (if value is specified)
    """
    
    # Check to make sure the value given is within the random variable's support
    if check_value(value,RVar.support)!=True:
        raise RVError('Value is not within the support of the random variable')

    # If the distribution is continuous, find and return the hf of the random variable
    if RVar.ftype[0]=='continuous':
        # If the distribution is already a hf, nothing needs to be done
        if RVar.ftype[1]=='hf':
            if value==x:
                return RVar
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        hfvalue=RVar.func[i].subs(x,value)
                        return simplify(hfvalue)
        # If the distribution is in chf form, use differentiation to find the hf
        if RVar.ftype[1]=='chf':
            X_dummy=CHF(RVar)
            # Generate a list of hf functions
            hflist=[]
            for i in range(len(X_dummy.func)):
                newfunc=diff(X_dummy.func[i],x)
                hflist.append(newfunc)
            if value==x:
                return RV(hflist,X_dummy.support,['continuous','hf'])
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        hfvalue=hflist[i].subs(x,value)
                        return simplify(hfvalue)
        # In all other cases, use the pdf and the sf to find the hf
        else:
            X_pdf=PDF(RVar).func
            X_sf=SF(RVar).func
            # Create a list of hf functions
            hflist=[]
            for i in range(len(RVar.func)):
                hfunc=(X_pdf[i])/(X_sf[i])
                hflist.append(simplify(hfunc))
            if value==x:
                return RV(hflist,RVar.support,['continuous','hf'])
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        hfvalue=hflist[i].subs(x,value)
                        return simplify(hfvalue)

    # If the random variable is discrete, find and return the hf
    if RVar.ftype[0]=='discrete':
        # If the distribution is already a hf, nothing needs to be done
        if RVar.ftype[1]=='hf':
            if value==x:
                return RVar
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return RVar.func[RVar.support.index(value)]
        # Otherwise, use the pdf and sf to find the hf
        else:
            X_pdf=PDF(RVar)
            X_sf=SF(RVar)
            hffunc=[]
            for i in range(len(X_pdf.func)):
                hffunc.append(X_pdf.func[i]/X_sf.func[i])
            if value==x:
                return RV(hffunc,X_pdf.support,['discrete','hf'])
            if value!=x:
                if value not in X_pdf.support:
                    return 0
                else:
                    return hffunc[X_pdf.support.index(value)]


def IDF(RVar,value=x):
    """
    Procedure Name: IDF
    Purpose: Compute the idf of a random variable
    Arguments:  1. RVar: A random variable
                2. value: An integer or floating point number (optional)
    Output:     1. IDF of a random variable (if value not specified)
                2. Value of the IDF at a given point (if value is specified)
    """
    
    # Check to make sure the percentile given is between 0 and 1
    if check_value(value,[0,1])!=True:
        raise RVError('Value is not within the support of the random variable')

    # If the distribution is continuous, find and return the idf of the random variable
    if RVar.ftype[0]=='continuous':
        if value==x:
            if RVar.ftype[1]=='idf':
                return self
            # Convert the random variable to its CDF form
            X_dummy=CDF(RVar)
            # Create values used to check for correct inverse
            check=[]
            for i in range(len(X_dummy.support)-1):
                if X_dummy.support[i]==-oo and X_dummy.support[i+1]==oo:
                    check.append(0)
                elif X_dummy.support[i]==-oo and X_dummy.support[i+1]!=oo:
                    check.append(X_dummy.support[i+1]-1)
                elif X_dummy.support[i]!=-oo and X_dummy.support[i+1]==oo:
                    check.append(X_dummy.support[i]+1)
                else:
                    check.append((X_dummy.support[i]+X_dummy.support[i+1])/2)
            # Use solve to create a list of candidate inverse functions
            # Check to see which of the candidate inverse functions is correct
            idffunc=[]
            for i in range(len(X_dummy.func)):
                invlist=solve(X_dummy.func[i]-t,x)
                if len(invlist)==1:
                    idffunc.append(invlist[0])
                else:
                    # The flag is used to determine if two separate inverses
                    #   could represent the inverse of the CDF. If this is the
                    #   case, an exception is raised
                    flag=False
                    for j in range(len(invlist)):
                        val=invlist[j].subs(t,X_dummy.func[i].subs(x,check[i])).evalf()
                        if abs(val-check[i])<.00001:
                            if flag==True:
                                raise RVError('Could not find the correct inverse')
                            idffunc.append(invlist[j])
                            flag=True
            # Create a list of supports for the IDF
            idfsup=[]
            for i in range(len(X_dummy.support)):
                idfsup.append(CDF(X_dummy,X_dummy.support[i]))
            # Return the IDF
            return RV(idffunc,idfsup,['continuous','idf'])
                    
            
        # If a value is specified, use the newton-raphson method to generate a random variate
        if value!=x:
            X_dummy=IDF(RVar)
            for i in range(len(X_dummy.support)):
                if value>=X_dummy.support[i] and value<=X_dummy.support[i+1]:
                    idfvalue=X_dummy.func[i].subs(t,value)
                    return simplify(idfvalue)
            #varlist=RVar.variate(s=value)
            #return varlist[0]


    # If the distribution is discrete, find and return the idf of the random variable
    if RVar.ftype[0]=='discrete':
        # If the distribution is already an idf, nothing needs to be done
        if RVar.ftype[1]=='idf':
            if value==x:
                return RVar
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return RVar.func[RVar.support.index(value)]
        # Otherwise, find the cdf, and then invert it
        else:
            # If the distribution is a chf or hf, convert to an sf first
            if RVar.ftype[1]=='chf' or RVar.ftype[1]=='hf':
                X_dummy0=SF(RVar)
                X_dummy=CDF(X_dummy0)
            else:
               X_dummy=CDF(RVar)
            if value==x:
                return RV(X_dummy.support,X_dummy.func,['discrete','idf'])
            if value!=x:
                if value not in RVar.func:
                    return 0
                else:
                    return RVar.support[RVar.func.index(value)]
            


def PDF(RVar,value=x):
    """
    Procedure Name: PDF
    Purpose: Compute the pdf of a random variable
    Arguments:  1. RVar: A random variable
                2. value: An integer or floating point number (optional)
    Output:     1. PDF of a random variable (if value not specified)
                2. Value of the PDF at a given point (if value is specified)
    """
    
    # Check to make sure the value given is within the random variable's support
    if check_value(value,RVar.support)!=True:
        raise RVError('Value is not within the support of the random variable')
    
    # If the distribution is continuous, find and return the pdf of the random variable
    if RVar.ftype[0]=='continuous':
        # If the distribution is already a pdf, nothing needs to be done
        if RVar.ftype[1]=='pdf':
            if value==x:
                return RVar
            if value!=x:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        pdfvalue=RVar.func[i].subs(x,value)
                        return simplify(pdfvalue)
        # If the distribution is a hf or chf, use integration to find the pdf
        if RVar.ftype[1]=='hf' or RVar.ftype[1]=='chf':
            X_dummy=HF(RVar)
            # Substitute the dummy variable 't' into the hazard function
            hfsubslist=[]
            for i in range(len(X_dummy.func)):
                newfunc=X_dummy.func[i].subs(x,t)
                hfsubslist.append(newfunc)
            # Integrate the hazard function
            intlist=[]
            for i in range(len(hfsubslist)):
                newfunc=integrate(hfsubslist[i],(t,X_dummy.support[i],x))
                # Correct the constant of integration
                if i!=0:
                    const=intlist[i-1].subs(x,X_dummy.support[i])-newfunc.subs(x,X_dummy.support[i])
                    newfunc=newfunc+const
                if i==0:
                    const=0-newfunc.subs(x,X_dummy.support[i])
                    newfunc=newfunc+const
                intlist.append(simplify(newfunc))
            # Multiply to find the pdf
            pdffunc=[]
            for i in range(len(intlist)):
                newfunc=X_dummy.func[i]*exp(-intlist[i])
                pdffunc.append(simplify(newfunc))
            if value==x:
                return RV(pdffunc,RVar.support,['continuous','pdf'])
            if value!=x:
                for i in range(len(X_dummy.support)):
                    if value>=X_dummy.support[i] and value<=X_dummy.support[i+1]:
                        pdfvalue=pdffunc[i].subs(x,value)
                        return simplify(pdfvalue)
        # In all other cases, find the pdf by differentiating the cdf
        else:
            X_dummy=CDF(RVar)
            if value==x:
                pdflist=[]
                for i in range(len(X_dummy.func)):
                    pdflist.append(diff(X_dummy.func[i],x))
                return RV(pdflist,RVar.support,['continuous','pdf'])
            if value!=x:
                for i in range(len(X_dummy.support)):
                    for i in range(len(X_dummy.support)):
                        if value>=X_dummy.support[i] and value<=X_dummy.support[i+1]:
                            pdffunc=diff(X_dummy.func[i],x)
                            pdfvalue=pdffunc.subs(x,value)
                            return simplify(pdfvalue)
                        
    # If the distribution is discrete, find and return the pdf of the random variable
    if RVar.ftype[0]=='discrete':
        # If the distribution is already a pdf, nothing needs to be done
        if RVar.ftype[1]=='pdf':
            if value==x:
                return RVar
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return RVar.func[RVar.support.index(value)]  
        # Otherwise, find the cdf of the random variable, and compute the pdf
        #   by finding differences
        else:
            X_dummy=CDF(RVar)
            pdffunc=[]
            for i in range(len(X_dummy.func)):
                if i==0:
                    pdffunc.append(X_dummy.func[i])
                else:
                    pdffunc.append(X_dummy.func[i]-X_dummy.func[i-1])
            if value==x:
                return RV(pdffunc,X_dummy.support,['discrete','pdf'])
            if value!=x:
                if value not in X_dummy.support:
                    return 0
                else:
                    return pdffunc.func[X_dummy.support.index(value)]

def SF(RVar,value=x):
    """
    Procedure Name: SF
    Purpose: Compute the SF of a random variable
    Arguments:  1. RVar: A random variable
                2. value: An integer or floating point number (optional)
    Output:     1. SF of a random variable (if value not specified)
                2. Value of the SF at a given point (if value is specified)
    """
    
    # Check to make sure the value given is within the random variable's support
    if check_value(value,RVar.support)!=True:
        raise RVError('Value is not within the support of the random variable')

    # If the distribution is continuous, find and return the sf of the random variable
    if RVar.ftype[0]=='continuous':
        # If the distribution is already a sf, nothing needs to be done
        if RVar.ftype[1]=='sf':
            if value==x:
                return RVar
            else:
                for i in range(len(RVar.support)):
                    if value>=RVar.support[i] and value<=RVar.support[i+1]:
                        sfvalue=RVar.func[i].subs(x,value)
                        return simplify(sfvalue)
        # If not, then use subtraction to find the sf
        else:
            X_dummy=CDF(RVar)
            # Compute the sf for each segment
            sflist=[]
            for i in range(len(X_dummy.func)):
                sflist.append(1-X_dummy.func[i])
            if value==x:
                return RV(sflist,RVar.support,['continuous','sf'])
            if value!=x:
                for i in range(len(X_dummy.support)):
                    if value>=X_dummy.support[i] and value<=X_dummy.support[i+1]:
                        sfvalue=sflist[i].subs(x,value)
                        return simplify(sfvalue)

    # If the distribution is discrete, find and return the sf of the random variable
    if RVar.ftype[0]=='discrete':
        # If the distribution is already an sf, nothing needs to be done
        if RVar.ftype[1]=='sf':
            if value==x:
                return RVar
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return RVar.func[RVar.support.index(value)]
        # If the distribution is a chf use exp(-chf) to find sf
        if RVar.ftype[1]=='chf':
            X_dummy=CHF(RVar)
            sffunc=[]
            for i in range(len(X_dummy.func)):
                sffunc.append(exp(-(X_dummy.func[i])))
            if value==x:
                return RV(sffunc,X_dummy.support,['discrete','sf'])
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return sffunc[RVar.support.index(value)]
        # If the distribution is a hf, use bootstrap rv to find sf:
        if RVar.ftype[1]=='hf':
            X_pdf=BootstrapRV(RVar.support)
            X_hf=RVar
            sffunc=[]
            for i in range(len(RVar.func)):
                sffunc.append(X_pdf.func[i]/X_hf.func[i])
            if value==x:
                return RV(sffunc,RVar.support,['discrete','sf'])
            if value!=x:
                if value not in RVar.support:
                    return 0
                else:
                    return sffunc[RVar.support.index(value)]
        # Otherwise, find the cdf of the random variable, and reverse the function
        #   argument
        else:
            X_dummy=CDF(RVar)
            newfunc=[]
            for i in reversed(range(len(X_dummy.func))):
                newfunc.append(X_dummy.func[i])
            Xsf=RV(newfunc,X_dummy.support,['discrete','sf'])
            if value==x:
                return Xsf
            if value!=x:
                if value not in Xsf.support:
                    return 0
                else:
                    return Xsf.func[Xsf.support.index(value)]
        

def BootstrapRV(varlist):
    """
    Procedure Name: Bootstrap RV
    Purpose: Generate a discrete random variable from a list of variates
    Arguments: 1. varlist: A list of variates
    Output:    1. A discrete random variable, where each element in the given variate
                    list is equally probable
    """
    # Sort the list of variables
    varlist.sort()
    # Find the number of elements in the list of variates
    numel=len(varlist)
    # Use varlist to generate the function and support for the random variable
    #   Count number of times element appears in varlist, divide by number
    #   of elements
    funclist=[]
    supplist=[]
    for i in range(len(varlist)):
        if varlist[i] not in funclist:
            supplist.append(varlist[i])
            funclist.append(varlist.count(varlist[i])/numel)
    # Return the result as a discrete random variable
    return RV(funclist,supplist,['discrete','pdf'])


"""
Procedures on One Random Variable

Procedures:
    1. ConvolutionIID(RVar,n)
    2. MaximumIID(RVar,n)
    3. Mean(RVar)
    4. MinimumIID(RVar,n)
    5. OrderStat(X,n,r)
    6. Transform(RVar,gX)
    7. Truncate(RVar,[lw,up])
    8. Variance(RVar)
"""

def ConvolutionIID(RVar,n):
    """
    Procedure Name: ConvolutionIID
    Purpose: Compute the convolution of n iid random variables
    Arguments:  1. RVar: A random variable
                2. n: an integer
    Output: The convolution of n iid random variables
    """
    # Check to make sure n is an integer
    if type(n)!=int:
        raise RVError('The second argument must be an integer')

    # Compute the iid convolution
    X_dummy=RVar
    for i in range(n-2):
        X_dummy+=X_dummy
    return X_dummy

def MaximumIID(RVar,n):
    """
    Procedure Name: MaximumIID
    Purpose: Comput the maximum of n iid random variables
    Arguments:  1. RVar: A random variable
                2. n: an integer
    Output:     1. The maximum of n iid random variables
    """
    # Check to make sure n is an integer
    if type(n)!=int:
        raise RVError('The second argument must be an integer')

    # Compute the iid maximum
    X_dummy=RVar
    for i in range(n-2):
        X_dummy=Maximum(X_dummy,X_dummy)
    return X_dummy

def Mean(RVar):
    """
    Procedure Name: Mean
    Purpose: Compute the mean of a random variable
    Arguments: 1. RVar: A random variable
    Output:    1. The mean of a random variable
    """
    # Find the PDF of the random variable
    X_dummy=PDF(RVar)
    # If the random variable is continuous, find and return the mean
    if RVar.ftype[0]=='continuous':
        # Create list of x*f(x)
        meanfunc=[]
        for i in range(len(X_dummy.func)):
            meanfunc.append(x*X_dummy.func[i])
        # Integrate to find the mean
        meanval=0
        for i in range(len(X_dummy.func)):
            val=integrate(meanfunc[i],(x,X_dummy.support[i],X_dummy.support[i+1]))
            meanval+=val
        return meanval

    # If the random variable is discrete, find and return the variance
    if RVar.ftype[0]=='discrete':
        # Create a list of x*f(x)
        meanlist=[]
        for i in range(len(X_dummy.func)):
            meanlist.append(X_dummy.func[i]*X_dummy.support[i])
        # Sum to find the mean
        meanval=0
        for i in range(len(meanlist)):
            meanval+=meanlist[i]
        return meanval

def MinimumIID(RVar,n):
    """
    Procedure Name: MinimumIID
    Purpose: Comput the minimum of n iid random variables
    Arguments:  1. RVar: A random variable
                2. n: an integer
    Output:     1. The minimum of n iid random variables
    """
    # Check to make sure n is an integer
    if type(n)!=int:
        raise RVError('The second argument must be an integer')

    # Compute the iid minimum
    X_dummy=RVar
    for i in range(n-2):
        X_dummy=Minimum(X_dummy,X_dummy)
    return X_dummy

def OrderStat(RVar,n,r):
    """
    Procedure Name: OrderStat
    Purpose: Compute the distribution of the rth order statistic
                from a sample puplation of n
    Arguments:  1. RVar: A random variable
                2. n: The number of items randomly drawn from the rv
                3. r: The index of the order statistic
    Output:     1. The desired r out of n OrderStatistic
    """
    if r>n:
        raise RVError('The index cannot be greater than the sample size')

    # If the distribution is continuous, find and return the value of the
    #   order statistic
    if RVar.ftype[0]=='continuous':
        # Compute the PDF, CDF and SF of the random variable
        pdf_dummy=PDF(RVar)
        cdf_dummy=CDF(RVar)
        sf_dummy=SF(RVar)
        # Compute the factorial constant
        const=(factorial(n))/(factorial(r-1)*factorial(n-r))
        # Compute the distribution of the order statistic for each
        #   segment
        ordstat_func=[]
        for i in range(len(RVar.func)):
            fx=pdf_dummy.func[i]
            Fx=cdf_dummy.func[i]
            Sx=sf_dummy.func[i]
            ordfunc=const*(Fx**(r-1))*(Sx**(n-r))*fx
            ordstat_func.append(ordfunc.simplify())
        # Return the distribution of the order statistic
        return RV(ordstat_func,RVar.support,['continuous','pdf'])
            
    

                
def Transform(RVar,gXt):
    """
    Procedure Name: Transform
    Purpose: Compute the transformation of a random variable
                by a a function g(x)
    Arguments:  1. RVar: A random variable
                2. gX: A transformation in list of two lists format
    Output:     1. The transformation of RVar       
    """
    
    # Check to make sure support of transform is in ascending order
    for i in range(len(gXt[1])-1):
        if gXt[1][i]>gXt[1][i+1]:
            raise RVError('Transform support is not in ascending order')

    # Convert the RV to its PDF form
    X_dummy=PDF(RVar)
            
    # If the distribution is continuous, find and return the transformation
    if RVar.ftype[0]=='continuous':
        # Adjust the transformation to include the support of the random
        #   variable
        gXold=[]
        for i in range(len(gXt)):
            gXold.append(gXt[i])
        gXsupp=[]
        for i in range(len(gXold[1])):
            gXsupp.append(gXold[1][i])
        # Add the support of the random variable into the support
        #   of the transformation
        for i in range(len(X_dummy.support)):
            if X_dummy.support[i] not in gXsupp:
                gXsupp.append(X_dummy.support[i])
        gXsupp.sort()
        # Find which segment of the transformation applies, and add it
        #   to the transformation list
        gXfunc=[]
        for i in range(len(gXsupp)-1):
            for j in range(len(gXold[0])):
                if gXsupp[i]>=gXold[1][j]:
                    if gXsupp[i]<=gXold[1][j+1]:
                        gXfunc.append(gXold[0][j])
                        break
        # Set the adjusted transformation as gX
        gX=[]
        gX.append(gXfunc)
        gX.append(gXsupp)
        # If the support of the transformation does not match up with the
        #   support of the RV, adjust the support of the transformation
        
        # Traverse list to find elements that are not within the support
        #   of the rv
        for i in range(len(gX[1])):
            if gX[1][i]<X_dummy.support[0]:
                gX[1][i]=X_dummy.support[0]
            if gX[1][i]>X_dummy.support[len(X_dummy.support)-1]:
                gX[1][i]=X_dummy.support[len(X_dummy.support)-1]
        # Delete segments of the transformation that will not be used
        for i in range(len(gX[0])-1):
            if gX[1][i]==gX[1][i+1]:
                gX[0].remove(gX[0][i])
                gX[1].remove(gX[1][i+1])
        # Create a list of mappings x->g(x)
        mapping=[]
        for i in range(len(gX[0])):
            mapping.append([gX[0][i].subs(x,gX[1][i]),
                            gX[0][i].subs(x,gX[1][i+1])])
        # Create the support for the transformed random variable
        trans_supp=[]
        for i in range(len(mapping)):
            for j in range(2):
                if mapping[i][j] not in trans_supp:
                    trans_supp.append(mapping[i][j])
        trans_supp.sort()
        # Find which segment of the transformation each transformation
        #   function applies to
        applist=[]
        for i in range(len(mapping)):
            temp=[]
            for j in range(len(trans_supp)-1):
                if min(mapping[i])<=trans_supp[j]:
                    if max(mapping[i])>=trans_supp[j+1]:
                        temp.append(j)
            applist.append(temp)
        # Find the appropriate inverse for each g(x)
        ginv=[]
        for i in range(len(gX[0])):
            # Find the 'test point' for the inverse
            if [gX[1][i],gX[1][i+1]]==[-oo,oo]:
                c=0
            elif gX[1][i]==-oo and gX[1][i+1]!=oo:
                c=gX[1][i+1]-1
            elif gX[1][i]!=-oo and gX[1][i+1]==oo:
                c=gX[1][i]+1
            else:
                c=(gX[1][i]+gX[1][i+1])/2
            # Create a list of possible inverses
            invlist=solve(gX[0][i]-t,x)
            # Use the test point to determine the correct inverse
            for j in range(len(invlist)):
                # If g-1(g(c))=c, then the inverse is correct
                if invlist[j].subs(t,gX[0][i].subs(x,c))==c:
                    ginv.append(invlist[j])
        # Find the transformation function for each segment
        seg_func=[]
        for i in range(len(X_dummy.func)):
            # Only find transformation for applicable segments
            for j in range(len(gX[0])):
                if gX[1][j]>=X_dummy.support[i]:
                    if gX[1][j+1]<=X_dummy.support[i+1]:
                        if type(X_dummy.func[i]) not in [float,int]:
                            tran=X_dummy.func[i].subs(x,ginv[j])*diff(ginv[j],t)
                        else:
                            tran=X_dummy.func[i]*diff(ginv[j],t)
                        seg_func.append(tran)
        # Sum the transformations for each piece of the transformed
        #   random variable
        trans_func=[]
        for i in range(len(trans_supp)-1):
            h=0
            for j in range(len(seg_func)):
                if i in applist[j]:
                    if mapping[j][0]<mapping[j][1]:
                        h=h+seg_func[j]
                    else:
                        h=h-seg_func[j]
            trans_func.append(h)
        # Substitute x into the transformed random variable
        trans_func2=[]
        for i in range(len(trans_func)):
            trans_func2.append(trans_func[i].subs(t,x).simplify())
        # Create and return the random variable
        return RV(trans_func2,trans_supp,['continuous','pdf'])

    # If the distribution is discrete, find and return the transformation
    if RVar.ftype[0]=='discrete':
        trans_sup=[]
        # Find the portion of the transformation each element
        #   in the random variable applies to, and then transform it
        for i in range(len(X_dummy.support)):
            for j in range(len(gX[1])-1):
                if X_dummy.support[i]>=gX[1][j]:
                    if X_dummy.support[i]<=gX[1][j+1]:
                        trans_sup.append(gX[0][j].subs(x,X_dummy.support[i]))
        # Sort the function and support lists for the convolution
        sortlist=zip(trans_sup,X_dummy.func)
        sortlist.sort()
        translist=[]
        funclist=[]
        for i in range(len(sortlist)):
            translist.append(sortlist[i][0])
            funclist.append(sortlist[i][1])
        # Combine redundant elements in the list
        translist2=[]
        funclist2=[]
        for i in range(len(translist)):
            if translist[i] not in translist2:
                translist2.append(translist[i])
                funclist2.append(funclist[i])
            elif translist[i] in translist2:
                idx=translist2.index(translist[i])
                funclist2[idx]+=funclist[i]
        # Return the transformed random variable
        return RV(funclist2,translist2,['discrete','pdf'])

def Truncate(RVar,supp):
    """
    Procedure Name: Truncate
    Purpose: Truncate a random variable
    Arguments: 1. RVar: A random variable
               2. supp: The support of the truncated random variable
    Output:    1. A truncated random variable
    """
    # Check to make sure the support of the truncated random
    #   variable is given in ascending order
    if supp[0]>supp[1]:
        raise RVError('The support must be given in ascending order')
    
    # Conver the random variable to its pdf form
    X_dummy=PDF(RVar)
    cdf_dummy=CDF(RVar)

    # If the random variable is continuous, find and return
    #   the truncated random variable
    if RVar.ftype[0]=='continuous':
        # Find the area of the truncated random variable
        area=CDF(cdf_dummy,supp[1])-CDF(cdf_dummy,supp[0])
        # Cut out parts of the distribution that don't fall
        #   within the new limits
        for i in range(len(X_dummy.func)):
            if supp[0]>=X_dummy.support[i]:
                if supp[0]<=X_dummy.support[i+1]:
                    lwindx=i
            if supp[1]>=X_dummy.support[i]:
                if supp[1]<=X_dummy.support[i+1]:
                    upindx=i
        truncfunc=[]
        for i in range(len(X_dummy.func)):
            if i>=lwindx and i<=upindx:
                truncfunc.append(X_dummy.func[i]/area)
        truncsupp=[supp[0]]
        upindx+=1
        for i in range(len(X_dummy.support)):
            if i>lwindx and i<upindx:
                truncsupp.append(X_dummy.support[i])
        truncsupp.append(supp[1])
        # Return the truncated random variable
        return RV(truncfunc,truncsupp,['continuous','pdf'])

    # If the distribution is discrete, find and return the
    #   truncated random variable
    if RVar.ftype[0]=='discrete':
        # Find the area of the truncated random variable
        area=0
        for i in range(len(X_dummy.support)):
            if X_dummy.support[i]>=supp[0]:
                if X_dummy.support[i]<=supp[1]:
                    area+=X_dummy.func[i]
        # Truncate the random variable and find the probability
        #   at each point
        truncfunc=[]
        truncsupp=[]
        for i in range(len(X_dummy.support)):
            if X_dummy.support[i]>=supp[0]:
                if X_dummy.support[i]<=supp[1]:
                    truncfunc.append(X_dummy.func[i]/area)
                    truncsupp.append(X_dummy.support[i])
        # Return the truncated random variable
        return RV(truncfunc,truncsupp,['discrete','pdf'])     


def Variance(RVar):
    """
    Procedure Name: Variance
    Purpose: Compute the variance of a random variable
    Arguments: 1. RVar: A random variable
    Output:    1. The variance of a random variable
    """
    # Find the PDF of the random variable
    X_dummy=PDF(RVar)
    # If the random variable is continuous, find and return the variance
    if RVar.ftype[0]=='continuous':
        # Find the mean of the random variable
        EX=Mean(X_dummy)
        # Find E(X^2)
        # Create list of (x**2)*f(x)
        varfunc=[]
        for i in range(len(X_dummy.func)):
            varfunc.append((x**2)*X_dummy.func[i])
        # Integrate to find E(X^2)
        exxval=0
        for i in range(len(X_dummy.func)):
            val=integrate(varfunc[i],(x,X_dummy.support[i],X_dummy.support[i+1]))
            exxval+=val
        # Find Var(X)=E(X^2)-E(X)^2
        var=exxval-(EX**2)
        return var

    # If the random variable is discrete, find and return the variance
    if RVar.ftype[0]=='discrete':
        # Find the mean of the random variable
        EX=Mean(X_dummy)
        # Find E(X^2)
        # Create a list of (x**2)*f(x)
        exxlist=[]
        for i in range(len(X_dummy.func)):
            exxlist.append(X_dummy.func[i]*(X_dummy.support[i])**2)
        # Sum to find E(X^2)
        exxval=0
        for i in range(len(exxlist)):
            exxval+=exxlist[i]
        # Find Var(X)=E(X^2)-E(X)^2
        var=exxval-(EX**2)
        return var

"""
Procedures on Two Random Variables

Procedures:
    1. Convolution(RVar1,RVar2)
    2. Maximum(RVar1,RVar2)
    3. Minimum(RVar1,RVar2)
    4. Product(RVar1,RVar2)
"""
                    
def Convolution(RVar1,RVar2):
    """
    Procedure Name: Convolution
    Purpose: Compute the convolution of two independent
                random variables
    Arguments:  1. RVar1: A random variable
                2. RVar2: A random variable
    Output:     1. The convolution of RVar1 and RVar2        
    """
    # If the two random variables are not both continuous or
    #   both discrete, return an error
    if RVar1.ftype[0]!=RVar2.ftype[0]:
        raise RVError('Both random variables must have the same type')

    # Convert both random variables to their PDF form
    X1_dummy=PDF(RVar1)
    X2_dummy=PDF(RVar2)

    # If the distributions are continuous, find and return the convolution
    #   of the two random variables
    if RVar1.ftype[0]=='continuous':
        # If the two distributions are both lifetime distributions, treat
        #   as a special case
        if RVar1.support==[0,oo] and RVar2.support==[0,oo] and len(RVar1.func)==1 and len(RVar2.func)==1:
            func1=X1_dummy.func[0]
            func2=X2_dummy.func[0].subs(x,z-x)
            conv=integrate(func1*func2,(x,0,z))
            return RV([conv.subs(z,x)],[0,oo],['continuous','pdf'])
        # Otherwise, compute the convolution using the product method
        else:
            gln=[[ln(x)],[0,oo]]
            ge=[[exp(x),exp(x)],[-oo,0,oo]]
            temp1=Transform(X1_dummy,ge)
            temp2=Transform(X2_dummy,ge)
            temp3=Product(temp1,temp2)
            fz=Transform(temp3,gln)
            convfunc=[]
            for i in range(len(fz.func)):
                convfunc.append(fz.func[i].simplify())
            return RV(convfunc,fz.support,['continuous','pdf'])
            

    # If the distributions are discrete, find and return the convolution
    #   of the two random variables.
    if RVar1.ftype[0]=='discrete':
        # Convert each random variable to its pdf form
        X1_dummy=PDF(RVar1)
        X2_dummy=PDF(RVar2)
        # Create function and support lists for the convolution of the
        #   two random variables
        convlist=[]
        funclist=[]
        for i in range(len(X1_dummy.support)):
            for j in range(len(X2_dummy.support)):
                convlist.append(X1_dummy.support[i]+X2_dummy.support[j])
                funclist.append(X1_dummy.func[i]*X2_dummy.func[j])
        # Sort the function and support lists for the convolution
        sortlist=zip(convlist,funclist)
        sortlist.sort()
        convlist2=[]
        funclist2=[]
        for i in range(len(sortlist)):
            convlist2.append(sortlist[i][0])
            funclist2.append(sortlist[i][1])
        # Remove redundant elements in the support list
        convlist3=[]
        funclist3=[]
        for i in range(len(convlist2)):
            if convlist2[i] not in convlist3:
                convlist3.append(convlist2[i])
                funclist3.append(funclist2[i])
            else:
                funclist3[convlist3.index(convlist2[i])]+=funclist2[i]
        # Create and return the new random variable
        return RV(funclist3,convlist3,['discrete','pdf'])

def Maximum(RVar1,RVar2):
    """
    Procedure Name: Maximum
    Purpose: Compute cdf of the maximum of RVar1 and RVar2
    Arguments:  1. RVar1: A random variable
                2. RVar2: A random variable
    Output:     1. The cdf of the maximum distribution
    """

    # If the two random variables are not of the same type
    #   raise an error
    if RVar1.ftype[0]!=RVar2.ftype[0]:
        raise RVError('The RVs must both be discrete or continuous')

    # If the distributions are continuous, find and return the max
    if RVar1.ftype[0]=='continuous':
        # Special case for lifetime distributions
        if RVar1.support==[0,oo] and RVar2.support==[0,oo]:
            cdf_dummy1=CDF(RVar1)
            cdf_dummy2=CDF(RVar2)
            cdf1=cdf_dummy1.func[0]
            cdf2=cdf_dummy2.func[0]
            maxfunc=cdf1*cdf2
            return RV(maxfunc.simplify(),[0,oo],['continuous','cdf'])
        # Otherwise, compute the min using the full algorithm
        Fx=CDF(RVar1)
        Fy=CDF(RVar2)
        # Create a support list for the 
        max_supp=[]
        for i in range(len(Fx.support)):
            if Fx.support[i] not in max_supp:
                max_supp.append(Fx.support[i])
        for i in range(len(Fy.support)):
            if Fy.support[i] not in max_supp:
                max_supp.append(Fy.support[i])
        max_supp.sort()
        # Remove any elements that are above the lower support max
        lowval=max(min(Fx.support),min(Fy.support))
        max_supp2=[]
        for i in range(len(max_supp)):
            if max_supp[i]>=lowval:
                max_supp2.append(max_supp[i])
        # Compute the minimum function for each segment
        xindx=0
        yindx=0
        max_func=[]
        for i in range(len(max_supp2)-1):
            if max_supp2[i]>Fx.support[0]:
                currFx=0
            elif max_supp2[i]==Fx.support[xindx]:
                currFx=Fx.func[xindx]
                xindx+=1
            if max_supp2[i]>Fy.support[yindx]:
                currFy=0
            elif max_supp2[i]==Fy.support[yindx]:
                currFy=Fy.func[yindx]
                yindx+=1
            Fmax=-(1-currFx)*(1-currFy)
            Fax=Fmax.simplify()
            max_func.append(Fmax)
        # Return the random variable
        return RV(max_func,max_supp2,['continuous','cdf'])

def Minimum(RVar1,RVar2):
    """
    Procedure Name: Minimum
    Purpose: Compute the distribution of the minimum of RVar1 and RVar2
    Arguments:  1. RVar1: A random variable
                2. RVar2: A random variable
    Output:     1. The minimum of the two random variables
    """

    # If the two random variables are not of the same type
    #   raise an error
    if RVar1.ftype[0]!=RVar2.ftype[0]:
        raise RVError('The RVs must both be discrete or continuous')

    # If the distributions are continuous, find and return the min
    if RVar1.ftype[0]=='continuous':
        # Special case for lifetime distributions
        if RVar1.support==[0,oo] and RVar2.support==[0,oo]:
            sf_dummy1=SF(RVar1)
            sf_dummy2=SF(RVar2)
            sf1=sf_dummy1.func[0]
            sf2=sf_dummy2.func[0]
            minfunc=1-(sf1*sf2)
            return RV(minfunc.simplify(),[0,oo],['continuous','cdf'])
        # Otherwise, compute the min using the full algorithm
        Fx=CDF(RVar1)
        Fy=CDF(RVar2)
        # Create a support list for the 
        min_supp=[]
        for i in range(len(Fx.support)):
            if Fx.support[i] not in min_supp:
                min_supp.append(Fx.support[i])
        for i in range(len(Fy.support)):
            if Fy.support[i] not in min_supp:
                min_supp.append(Fy.support[i])
        min_supp.sort()
        # Remove any elements that are above the lower support max
        highval=min(max(Fx.support),max(Fy.support))
        min_supp2=[]
        for i in range(len(min_supp)):
            if min_supp[i]<=highval:
                min_supp2.append(min_supp[i])
        # Compute the minimum function for each segment
        xindx=0
        yindx=0
        min_func=[]
        for i in range(len(min_supp2)-1):
            if min_supp2[i]<Fx.support[0]:
                currFx=0
            elif min_supp2[i]==Fx.support[xindx]:
                currFx=Fx.func[xindx]
                xindx+=1
            if min_supp2[i]<Fy.support[yindx]:
                currFy=0
            elif min_supp2[i]==Fy.support[yindx]:
                currFy=Fy.func[yindx]
                yindx+=1
            Fmin=1-(1-currFx)*(1-currFy)
            Fmin=Fmin.simplify()
            min_func.append(Fmin)
        # Return the random variable
        return RV(min_func,min_supp2,['continuous','cdf'])

def Product(RVar1,RVar2):
    """
    Procedure Name: Product
    Purpose: Compute the product of two independent
                random variables
    Arguments:  1. RVar1: A random variable
                2. RVar2: A random variable
    Output:     1. The product of RVar1 and RVar2        
    """
    #
    # Needs further debugging, as well as quadrants II,III,IV 
    #
    
    # If the random variable is continuous, find and return the
    #   product of the two random variables
    if RVar1.ftype[0]=='continuous':
        v=Symbol('v')
        # Place zero in the support of X if it is not there already
        X1=PDF(RVar1)
        xfunc=[]
        xsupp=[]
        for i in range(len(X1.func)):
            xfunc.append(X1.func[i])
            xsupp.append(X1.support[i])
            if X1.support[i]<0:
                if X1.support[i+1]>0:
                    xfunc.append(X1.func[i])
                    xsupp.append(0)
        xsupp.append(X1.support[len(X1.support)-1])
        X_dummy=RV(xfunc,xsupp,['continuous','pdf'])
        # Place zero in the support of Y if it is not already there
        Y1=PDF(RVar2)
        yfunc=[]
        ysupp=[]
        for i in range(len(Y1.func)):
            yfunc.append(Y1.func[i])
            ysupp.append(Y1.support[i])
            if Y1.support[i]<0:
                if Y1.support[i+1]>0:
                    yfunc.append(Y1.func[i])
                    ysupp.append(0)
        ysupp.append(Y1.support[len(Y1.support)-1])
        Y_dummy=RV(yfunc,ysupp,['continuous','pdf'])
        # Initialize the support list for the product V=X*Y
        vsupp=[]
        for i in range(len(X_dummy.support)):
            for j in range(len(Y_dummy.support)):
                val=X_dummy.support[i]*Y_dummy.support[j]
                if val not in vsupp:
                    vsupp.append(val)
        vsupp.sort()
        # Initialize the pdf segments of v
        vfunc=[]
        for i in range(len(vsupp)-1):
            vfunc.append(0)
        # Loop through each piecewise segment of X
        for i in range(len(X_dummy.func)):
            # Loop through each piecewise segment of Y
            for j in range(len(Y_dummy.func)):
                # Define the corner of the rectangular region
                a=X_dummy.support[i]
                b=X_dummy.support[i+1]
                c=Y_dummy.support[j]
                d=Y_dummy.support[j+1]
                # If the region is in the first quadrant, compute the
                #   required integrals sequentially
                if a>=0 and c>=0:
                    if type(Y_dummy.func[j]) not in [float,int]:
                        gj=Y_dummy.func[j].subs(x,v/x)
                    else:
                        gj=Y_dummy.func[j]
                    fi=X_dummy.func[i]
                    pv=integrate(fi*gj*(1/x),(x,a,b))
                    if d<oo:
                        qv=integrate(fi*gj*(1/x),(x,v/d,b))
                    if c>0:
                        rv=integrate(fi*gj*(1/x),(x,a,v/c))
                    if c>0 and d<oo and a*d<b*c:
                        sv=integrate(fi*gj*(1/x),(x,v/d,v/c))
                    # 1st Qd, Scenario 1
                    if c==0 and d==oo:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=0:
                                vfunc[k]+=pv
                    # 1st Qd, Scenario 2
                    if c==0 and d<oo:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=0 and v[k+1]<=a*d:
                                vfunc[k]+=pv
                            if vsupp[k]>=a*d and v[k+1]<=b*d:
                                vfunc[k]+=qv
                    # 1st Qd, Scenario 3
                    if c>0 and d==oo:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=b*c:
                                vfunc[k]+=pv
                            if vsupp[k]>=a*c and v[k+1]<=b*c:
                                vfunc[k]+=rv
                    # 1st Qd, Scenario 4
                    if c>0 and d<oo:
                        # Case 1
                        if a*d<b*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=a*c and vsupp[k+1]<=a*d:
                                    vfunc[k]+=rv
                                if vsupp[k]>=a*d and vsupp[k+1]<=b*c:
                                    vfunc[k]+=sv
                                if vsupp[k]>=b*c and vsupp[k+1]<=b*d:
                                    vfunc[k]+=qv
                        # Case 2
                        if a*d==b*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=a*c and vsupp[k+1]<=a*d:
                                    vfunc[k]+=rv
                                if vsupp[k]>=b*c and vsupp[k+1]<=b*d:
                                    vfunc[k]+=qv
                        # Case 3
                        if a*d>b*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=a*c and vsupp[k+1]<=b*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=b*c and vsupp[k+1]<=a*d:
                                    vfunc[k]+=pv
                                if vsupp[k]>=a*d and vsupp[k+1]<=b*d:
                                    vfunc[k]+=qv
                # If the region is in the second quadrant, compute
                #   the required integrals sequentially
                if a<0 and c<0:
                    if type(Y_dummy.func[j]) not in [float,int]:
                        gj=Y_dummy.func[j].subs(x,v/x)
                    else:
                        gj=Y_dummy.func[j]
                    fi=X_dummy.func[i]
                    pv=-integrate(fi*gj*(1/x),(x,a,b))
                    if d<0:
                        qv=-integrate(fi*gj*(1/x),(x,(v/d),b))
                    if c>-oo:
                        rv=-integrate(fi*gj*(1/x),(x,a,(v/c)))
                    if c>-oo and d<0:
                        sv=-integrate(fi*gj*(1/x),(x,(v/d),(v/c)))
                    # 2nd Qd, Scenario 1
                    if c==-oo and d==0:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=0:
                                vfunc[k]+=pv
                    # 2nd Qd, Scenario 2
                    if c==-oo and d<0:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=a*d and vsupp[k+1]<=oo:
                                vfunc[k]+=pv
                            if vsupp[k]>=b*d and vsupp[k+1]<=a*d:
                                vfunc[k]+=qv
                    # 2nd Qd, Scenario 3
                    if c>-oo and d==0:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=0 and vsupp[k+1]<=b*c:
                                vfunc[k]+=pv
                            if vsupp[k]>=b*c and vsupp[k+1]<=a*c:
                                vfunc[k]+=rv
                    # 2nd Qd, Scenario 4
                    if c>-oo and d<0:
                        # Case 1
                        if a*d>b*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=a*d and vsupp[k+1]<=a*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=b*c and vsupp[k+1]<=a*d:
                                    vfunc[k]+=sv
                                if vsupp[k]>=b*d and vsupp[k+1]<=b*c:
                                    vfunc[k]+=qv
                        # Case 2
                        if a*d==b*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=a*d and vsupp[k+1]<=a*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=b*d and vsupp[k+1]<=b*c:
                                    vfunc[k]+=qv
                        # Case 3
                        if a*d<b*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=b*c and vsupp[k+1]<=a*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=a*d and vsupp[k+1]<=b*c:
                                    vfunc[k]+=pv
                                if vsupp[k]>=b*d and vsupp[k+1]<=a*d:
                                    vfunc[k]+=qv
                # If the region is in the third quadrant, compute
                #   the required integrals sequentially
                if a<0 and c>=0:
                    if type(Y_dummy.func[j]) not in [float,int]:
                        gj=Y_dummy.func[j].subs(x,v/x)
                    else:
                        gj=Y_dummy.func[j]
                    fi=X_dummy.func[i]
                    pv=-integral(fi*gj*(1/x),(x,a,b))
                    if d<oo:
                        qv=-integral(fi*gj*(1/x),(x,a,(v/d)))
                    if c>0:
                        rv=-integral(fi*gj*(1/x),(x,(v/c),b))
                    if c>0 and d<oo:
                        sv=-integral(fi*gj*(1/x),(x,(v/c),(v/d)))
                    # 3rd Qd, Scenario 1
                    if c==0 and d==oo:
                        for k in range(len(vfunc)):
                            if vsupp[k+1]<=0:
                                vfunc[k]+=pv
                    # 3rd Qd, Scenario 2
                    if c==0 and d<oo:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=b*d and vsupp[k+1]<=0:
                                vfunc[k]+=pv
                            if vsupp[k]>=a*d and vsupp[k+1]<=b*d:
                                vfunc[k]+=qv
                    # 3rd Qd, Scenario 3
                    if c>0 and d<oo:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=-oo and vsupp[k+1]<=a*c:
                                vfunc[k]+=pv
                            if vsupp[k]>=a*c and vsupp[k+1]<=b*c:
                                vfunc[k]+=rv
                    # 3rd Qd, Scenario 4
                    if c>0 and d<oo:
                        # Case 1
                        if b*d>a*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=b*d and vsupp[k+1]<=b*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=a*c and vsupp[k+1]<=b*d:
                                    vfunc[k]+=sv
                                if vsupp[k]>=a*d and vsupp[k+1]<=a*c:
                                    vfunc[k]+=qv
                        # Case 2
                        if a*c==b*d:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=a*d and vsupp[k+1]<=a*c:
                                    vfunc[k]+=qv
                                if vsupp[k]>=b*d and vsupp[k+1]<=b*c:
                                    vfunc[k]+=rv
                        # Case 3
                        if a*c>b*d:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=a*c and vsupp[k+1]<=b*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=b*d and vsupp[k+1]<=a*c:
                                    vfunc[k]+=pv
                                if vsupp[k]>=a*d and vsupp[k+1]<=b*d:
                                    vfunc[k]+=qv
                # If the region is in the fourth quadrant, compute
                #   the required integrals sequentially
                if a>=0 and c<0:
                    if type(Y_dummy.func[j]) not in [float,int]:
                        gj=Y_dummy.func[j].subs(x,v/x)
                    else:
                        gj=Y_dummy.func[j]
                    fi=X_dummy.func[i]
                    pv=integrate(fi*gj*(1/x),(x,a,b))
                    if d<0:
                        qv=integrate(fi*gj*(1/x),(x,a,(v/d)))
                    if c>-oo:
                        rv=integrate(fi*gj*(1/x),(x,(v/c),b))
                    if c>-oo and d<0:
                        sv=integrate(fi*gj*(1/x),(x,(v/c),(v/d)))
                    # 4th Qd, Scenario 1
                    if c==oo and d==0:
                        for k in range(len(vfunc)):
                            if vsupp[k+1]<=0:
                                vfunc[k]+=pv
                    # 4th Qd, Scenario 2
                    if c==oo and d<0:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=-oo and vsupp[k+1]<=b*d:
                                vfunc[k]+=pv
                            if vsupp[k]>=b*d and vsupp[k+1]<=a*d:
                                vfunc[k]+=qv
                    # 4th Qd, Scenario 3
                    if c>-oo and d==0:
                        for k in range(len(vfunc)):
                            if vsupp[k]>=a*c and vsupp[k+1]<=0:
                                vfunc[k]+=pv
                            if vsupp[k]>=b*c and vsupp[k+1]<=a*c:
                                vfunc[k]+=rv
                    # 4th Qd, Scenario 4
                    if c>-oo and d<0:
                        # Case 1
                        if a*c>b*d:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=b*c and vsupp[k+1]<=b*d:
                                    vfunc[k]+=rv
                                if vsupp[k]>=b*d and vsupp[k+1]<=a*c:
                                    vfunc[k]+=sv
                                if vsupp[k]>=a*c and vsupp[k+1]<=a*d:
                                    vfunc[k]+=qv
                        # Case 2
                        if a*d==b*c:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=b*c and vsupp[k+1]<=a*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=a*c and vsupp[k+1]<=a*d:
                                    vfunc[k]+=qv
                        # Case 3
                        if a*c<b*d:
                            for k in range(len(vfunc)):
                                if vsupp[k]>=b*c and vsupp[k+1]<=a*c:
                                    vfunc[k]+=rv
                                if vsupp[k]>=a*c and vsupp[k+1]<=b*d:
                                    vfunc[k]+=pv
                                if vsupp[k]>=b*d and vsupp[k+1]<=a*d:
                                    vfunc[k]+=qv                   
        vfunc_final=[]
        for i in range(len(vfunc)):
            if type(vfunc[i]) not in [int,float]:
                vfunc_final.append(vfunc[i].subs(v,x))
            else:
                vfunc_final.append(vfunc[i])
        return RV(vfunc_final,vsupp,['continuous','pdf'])
    # If the distributions are discrete, find and return the product
    #   of the two random variables.
    if RVar1.ftype[0]=='discrete':
        # Convert each random variable to its pdf form
        X1_dummy=PDF(RVar1)
        X2_dummy=PDF(RVar2)
        # Create function and support lists for the product of the
        #   two random variables
        prodlist=[]
        funclist=[]
        for i in range(len(X1_dummy.support)):
            for j in range(len(X2_dummy.support)):
                prodlist.append(X1_dummy.support[i]*X2_dummy.support[j])
                funclist.append(X1_dummy.func[i]*X2_dummy.func[j])
        # Sort the function and support lists for the convolution
        sortlist=zip(prodlist,funclist)
        sortlist.sort()
        prodlist2=[]
        funclist2=[]
        for i in range(len(sortlist)):
            prodlist2.append(sortlist[i][0])
            funclist2.append(sortlist[i][1])
        # Remove redundant elements in the support list
        prodlist3=[]
        funclist3=[]
        for i in range(len(prodlist2)):
            if prodlist2[i] not in prodlist3:
                prodlist3.append(prodlist2[i])
                funclist3.append(funclist2[i])
            else:
                funclist3[prodlist3.index(prodlist2[i])]+=funclist2[i]
        # Create and return the new random variable
        return RV(funclist3,prodlist3,['discrete','pdf'])

"""
Utilities

Procedures:
    1. PlotDist(RVar,suplist)
    2. PlotDist(plot_list,suplist)
"""

def PlotDist(RVar,suplist=None,opt=None):
    """
    Procedure: Plot Dist
    Purpose: Plot a random variable
    Arguments:  1. RVar: A random variable
                2. suplist: A list of supports for the plot
    Output:     1. A plot of the random variable
    """
    # Create the labels for the plot
    if RVar.ftype[1]=='cdf':
        lab1='F(x)'
        lab2='Cumulative Distribution Function'
    elif RVar.ftype[1]=='chf':
        lab1='H(x)'
        lab2='Cumulative Hazard Function'
    elif RVar.ftype[1]=='hf':
        lab1='h(x)'
        lab2='Hazard Function'
    elif RVar.ftype[1]=='idf':
        lab1='F-1(s)'
        lab2='Inverse Density Function'
    elif RVar.ftype[1]=='pdf':
        lab1='f(x)'
        lab2='Probability Density Function'
    elif RVar.ftype[1]=='sf':
        lab1='S(X)'
        lab2='Survivor Function'

    # If the distribution is continuous, plot the function
    if RVar.ftype[0]=='continuous':
        # Create a list of supports for the plot
        if suplist==None:
            suplist=[RVar.support[0],
                     RVar.support[len(RVar.support)-1]]
            if -oo in suplist:
                suplist[0]=float(RVar.variate(s=.01)[0])
            if oo in suplist:
                suplist[1]=float(RVar.variate(s=.99)[0])
        plot_sup=[]
        if suplist[0]>RVar.support[0]:
            plot_sup.append(float(suplist[0]))
        for i in range(len(RVar.support)):
            if RVar.support[i]>=suplist[0]:
                if RVar.support[i]<=suplist[1]:
                    if RVar.support[i] not in plot_sup:
                        plot_sup.append(float(RVar.support[i]))
        if suplist[1]<RVar.support[len(RVar.support)-1]:
            plot_sup.append(float(suplist[1]))
        print plot_sup
        # Create a list of functions for the plot
        for i in range(len(RVar.func)):
            if suplist[0]>=RVar.support[i]:
                if suplist[0]<=RVar.support[i+1]:
                    lwindx=i
            if suplist[1]>=RVar.support[i]:
                if suplist[1]<=RVar.support[i+1]:
                    upindx=i
        plot_func=[]
        for i in range(len(RVar.func)):
            if i>=lwindx and i<=upindx:
                strfunc=str(RVar.func[i])
                plot_func.append(strfunc)
        plot_sup=[suplist[0]]
        upindx+=1
        for i in range(len(RVar.support)):
            if i>lwindx and i<upindx:
                plot_sup.append(RVar.support[i])
        plot_sup.append(suplist[1])
        print plot_func
        plt.mat_plot(plot_func,plot_sup,lab1,lab2,'continuous')
        if opt!='display':
            plt.show()
    if RVar.ftype[0]=='discrete':
        plt.mat_plot(RVar.func,RVar.support,lab1,lab2,'discrete')
        if opt!='display':
            plt.show()    

def PlotDisplay(plot_list,suplist=None):
    # Create a plot of each random variable in the plot list
    for i in range(len(plot_list)):
        PlotDist(plot_list[i],suplist,'display')
    plt.show()
        
                                         
                                         
        
            
                           
            
            
                    
            
            
                
        
        






























            
