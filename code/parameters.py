#!/usr/bin/env python
# $URL$
# $Rev$
#
# Daniel Rothenberg, 2011-06-07
"""Stub project parameters file. 

Eventually, this will be generated as part of the workflow for using the
homogenizaiton library by the end-user. For now, it can be used by importing
into another workflow.

"""
__docformat__ = "restructuredtext"

#: The following values are constants 

#: The radius of the earth, in kilometers
RADIUS_EARTH = 6371.0 


class Parameters(object):
    
    def __init__(self, **params):
        self.__dict__.update(params)
        
    def __repr__(self):
        return "Parameters (%r)" % self.__dict__

####################################################################

def default_parameters(**user_params):
        
    project = "TEST"
    """This is the reference name for the project. It will be used for
    creating output files such as logs or figures, and should be a string 
    with all capital letters."""
    
    begyr = 1900
    """First year in the interval over which the homogenization analysis 
    will occur. Few if any USHCN data series start earlier than 1890, so this
    makes a good minimum value to consider; if you go to far earlier, then there 
    will be lots of missing data filled into the series."""
    
    endyr = 2000
    """Last year in the interval over which the homogenization analysis will
    occur. The dataset currently available on the CDC FTP site only goes to 2009,
    so this value shouldn't be higher than that!"""
    
    nsnet = 21
    """The total number of stations comprising the network being run through the
    homogenization routine."""
    
    nstns = 21
    """The total number of sub-network or "neighborhood" stations in the analysis.
    Generally, there is only one sub-network in a given set of USHCN stations, so 
    this should be the same as nsnet."""
    
    maxstns = 21
    """The maximum number of stations in the entire network. This should also be
    the same as nsnet."""
    
    maxfound = 2000
    """The maximum number of inhomogenity hits possible to find a single iteration
    of a homogenization analysis. This was used by Menne/Williams to help avoid
    the possibility of entering infinite loops."""
    
    amiss = -9999
    """The value to use in place of missing data."""
    
    numyr = endyr - begyr
    """The total number of years over which the homogenization analysis is
    being performed."""
    
    nmo = numyr*12
    """The total number of months over which the homogenization analysis is 
    being performed."""
    
    pdict = dict(nmo=nmo, numyr=numyr, amiss=amiss, maxfound=maxfound, 
                 maxstns=maxstns, nstns=nstns, nsnet=nsnet, endyr=endyr, 
                 begyr=begyr, project=project)
    pdict.update(user_params)
    
    p = Parameters(**pdict)
    return p
    


