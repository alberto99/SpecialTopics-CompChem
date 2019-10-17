#! /usr/bin/python
#
# Replica.py    

from .Config import *
from . import Chain
from . import Monty
from .Trajectory import *

import random
import string
import math
import sys
import os

class Replica:
    """A container object, to hold the Chain() and Monty() objects"""

    def __init__(self,config,repnum):
        """Initialize the Replica() object."""
	
        print(repnum)
        temp = config.REPLICATEMPS[repnum]
        self.repnum = repnum
        self.repname = 'rep' + string.zfill(str(repnum),2)
        self.repdir = config.DATADIR + self.repname
        self.chain = Chain.Chain(config)
        self.mc = Monty.Monty(config,temp,self.chain)    
        #Implement linear factor scheme for force constants
        if (len(config.REPLICATEMPS) > 1):
            self.mc.factor = 1-float(repnum)/(len(config.REPLICATEMPS)-1)
        else:
            self.mc.factor = 1

