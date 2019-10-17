#! /usr/bin/python

usage = """
mcrex.py <configfile>

Try:  mcrex.py mcrex.conf

This program will read in an HP chain and run parameters specified in the configure file,
and perform a replica exchange Monte Carlo simulation.

For the example "mcrex.conf", an 11-mer sequence is simulated, and the program ends when
the native conformation (contact state) is found.
A directory of results is output to directory ./mcrex_data
"""

import sys

from Config import *
from Chain import *
from Monty import *
from Replica import *
from Trajectory import *

import numpy
import random
import string
import math
import os

g = random.Random(randseed)


if len(sys.argv) < 2:
    print usage
    sys.exit(1)
 
 
VERBOSE = True
    

    
##########################
# Functions 

def setup_restraints_for_replica(replicas,config):
    '''Here we are going to assign a vector with which springs should 
    be active, after that the springs for each replica should be updated'''
    springs = config.POOL_SPRINGS
    for i in range(config.NREPLICAS):
        ene = []
        print "replica %i %s %f %f" % (i, replicas[i].chain.coords, replicas[i].mc.energy(replicas[i].chain),replicas[i].mc.restraint.energy(replicas[i].chain)),replicas[i].mc.active,replicas[i].mc.temp,replicas[i].mc.factor
        for j in range(len(springs)):
            alpha = replicas[i].mc.factor
            s1 = Monty.DistRestraint([springs[j]],config.KSPRING*alpha)
            ene.append(s1.energy(replicas[i].chain))
        for step in range(10*len(springs)**2):
            j = g.randrange(len(springs))
            k = g.randrange(len(springs))
            if replicas[i].mc.active[j] == replicas[i].mc.active[k]:
                continue

            #Flip weights:
            new_j_weight = replicas[i].mc.active[k]
            new_k_weight = replicas[i].mc.active[j]
            #Compute delta energy
            delta_e = (ene[j] * (new_j_weight - replicas[i].mc.active[j]) + 
                      ene[k] * (new_k_weight - replicas[i].mc.active[k]) )

            #Accept or reject
            delta_e =  1.0 / (1.987e-3 * replicas[i].mc.temp) * delta_e
            if delta_e < 0:
                accept = True
            else:
                metrop = math.exp(-delta_e)
                trial = g.random()
                if trial < metrop:
                    accept = True
                else:
                    accept = False

            if accept:
                replicas[i].mc.active[j] = new_j_weight
                replicas[i].mc.active[k] = new_k_weight
        #Now we are done with replica i: setup for run
        active = numpy.array(replicas[i].mc.active)
        assert len(springs) == len(active)
        ind = active == 1
        restraints = [a for a,b in zip(springs,ind) if b]
        print "restraints to put: ",restraints
        print "lastenergy before: ",replicas[i].mc.energy(replicas[i].chain),replicas[i].mc.restraint.energy(replicas[i].chain),replicas[i].mc.lastenergy
        print "restraints had: ",replicas[i].mc.restraint.contacts
        replicas[i].mc.restraint = Monty.DistRestraint(restraints,config.KSPRING*alpha)
        replicas[i].mc.lastenergy = replicas[i].mc.energy(replicas[i].chain) + replicas[i].mc.restraint.energy(replicas[i].chain)
        print "lastenergy after: ",replicas[i].mc.energy(replicas[i].chain),replicas[i].mc.restraint.energy(replicas[i].chain),replicas[i].mc.lastenergy
        print "restraints put: ",replicas[i].mc.restraint.contacts
    return replicas




def attemptswap(replicas,swaps,viableswaps): 

    #Construct matrix with energies:
    U = numpy.zeros( (config.NREPLICAS, config.NREPLICAS) )
    for i in range(config.NREPLICAS):
        pot = replicas[i].mc.energy(replicas[i].chain)

        for j in range(config.NREPLICAS):
            t = replicas[j].mc.temp
            res = replicas[i].mc.restraint.energy(replicas[i].chain) * replicas[j].mc.factor
            U[i,j] = 1.0 / (1.987e-3 * t) * (pot + res)

    for step in xrange(10 * config.NREPLICAS**2):
        # Attempt a swap between replicas    
        if config.SWAPMETHOD == 'random pair':
            # pick pair at random
            r = g.random()
            i = min(int(r*config.NREPLICAS),(config.NREPLICAS-1))
            j = i
            while (j==i):
                s = g.random()
                j = min(int(s*config.NREPLICAS),(config.NREPLICAS-1))
            # make sure j > i (needed for the swap criterion below)
            if j < i:
                tmp = i
                i = j
                j = tmp
        
        elif config.SWAPMETHOD == 'neighbors':    
            # pick neighboring pair at random 
            r = g.random()
            i = min(int(r*(config.NREPLICAS-1)),(config.NREPLICAS-2))
            j = i+1
                    
        else:
            print 'Swap method', config.SWAPMETHOD, 'unknown.  Exiting...'
            sys.exit(1)


        #if (VERBOSE):
        #  print 'REX: Attempting swap between replicas',i,'and',j
         
        randnum = g.random()

        ### if proposing i-->j, 
        ### 
        ###   del = -[U(xj,Ti,Sj(Tifactor)) + U(xi,Tj,Si(Tjfactor))] +
        ###          [U(xi,Ti,Si(Tifactor)) + U(xj,Tj,Sj(Tjfactor))]

        delfactor = -( U[j,i] + U[i,j] ) + (U[i,i] + U[j,j])

        boltzfactor = math.exp(delfactor)
        if randnum < boltzfactor:

            # swap the ***temperatures***
            temp_i = replicas[i].mc.temp
            temp_j = replicas[j].mc.temp
            replicas[i].mc.temp = temp_j
            replicas[j].mc.temp = temp_i
            # and the factors:
            factor_i = replicas[i].mc.factor
            replicas[i].mc.factor =  replicas[j].mc.factor
            replicas[j].mc.factor = factor_i
            retval = 1
        else:
            retval = 0

    for i in range(config.NREPLICAS):
        old_indexfromrep = replicas[i].mc.tempfromrep
        if old_indexfromrep == config.REPLICATEMPS.index(replicas[i].mc.temp):
            #No swap
            retval = 0
        else:
            replicas[i].mc.tempfromrep = config.REPLICATEMPS.index(replicas[i].mc.temp)
            retval = 1
        viableswaps[i] = viableswaps[i] + retval
        swaps[i] = swaps[i] + 1

    return 1


#####################################
# Main Program
    
if __name__ == '__main__':
    
    # load in config file
    configfile = sys.argv[1]
    config = Config()
    config.read_configfile(configfile)
    if VERBOSE:
        config.print_config()

    # look up the native contact state from a pre-compiled flat file     
    if config.STOPATNATIVE == 1:
        nativeclistfile = config.NATIVEDIR + '/' + config.HPSTRING + '.clist'
        fnative = open(nativeclistfile,'r')
        nativeclist = eval(fnative.readline())
        fnative.close()
        if VERBOSE:
            print 'NATIVE CLIST:',nativeclist
    
    
    # Make a list of config.NREPLICAS Replica objects
    replicas = []           # a list of Replica objects
    for i in range(0,config.NREPLICAS):
        replicas.append(Replica(config,i))

    #Configure which restraints are initially active
    n = config.ACTIVE_SPRINGS
    for i in range(config.NREPLICAS):
        x = numpy.zeros(len(config.POOL_SPRINGS))
        a = list(g.sample(numpy.arange(len(config.POOL_SPRINGS)),n))
        for aa in a:
            x[aa] = 1
        replicas[i].mc.active = x
    replicas = setup_restraints_for_replica(replicas,config)

    
    # Each replica is just a container for the objects:
    #     replica[i].chain   <-- the HP chain
    #     replica[i].mc      <-- A Monte Carlo move set object that performs operations on the chain
    
    # a trajectory object to handle the work of writing trajectories
    traj = Trajectory(replicas,config)  
    if VERBOSE:
        print 'Trajectory REPFILES:',traj.repfiles_trj


    # Keep track of statistics for the replica exchange simulation, for each replica:

    ### Monte Carlo stats
    steps = []                 # number of total move steps attempted (viable or not)
    viablesteps = []           # number of viable steps
    acceptedsteps = []         # number of viable steps that were accepted under Metropolis criterion 
    moveacceptance = []        # fraction of 
    acceptance = []            # fraction of MC moves accepted

    ### Replica exchange stats 
    swaps = []                 # the number of attemped replica swaps
    viableswaps = []           # the number of viable swaps
    swap_acceptance = []       # the fraction of swaps that were accepted
    
    for i in range(0,config.NREPLICAS):
         steps.append(0)
         viablesteps.append(0)    
         acceptedsteps.append(0)    
         moveacceptance.append(0)    
         acceptance.append(0)    
         swaps.append(0)
         viableswaps.append(0)
         swap_acceptance.append(0)
    
    
    prodstep = 1
    foundnative = 0
    while (prodstep < (config.MCSTEPS+1))&(foundnative==0):
    
        # Run the replicas for a production cycle...
        for rep in range(0,config.NREPLICAS):

            ### Propose a new MC move
            if string.strip(config.MOVESET) == 'MS1':
                replicas[rep].mc.move1(replicas[rep])
            elif string.strip(config.MOVESET) == 'MS2':
                replicas[rep].mc.move2(replicas[rep])
            elif string.strip(config.MOVESET) == 'MS3':
                replicas[rep].mc.move3(replicas[rep])
            elif string.strip(config.MOVESET) == 'MS4':
                replicas[rep].mc.move4(replicas[rep])
            else: 
                print 'MC MOVESET=',config.MOVESET,'not supported!'
                print 'Exiting....'
                sys.exit(1)

            if replicas[rep].chain.nextviable == 1:     # count this move only if the chain is viable
                ### accept with metroplis criterion
                accept = replicas[rep].mc.metropolis(replicas[rep])
                # if accept is yes, the chain is updated (see Monty.py)
                if (accept):
                    acceptedsteps[rep] = acceptedsteps[rep] + 1
                # keep track of MC steps
                viablesteps[rep] = viablesteps[rep] + 1

            # keep track of total steps
            steps[rep] = steps[rep] + 1

            # write viable steps to the trajectory with specified frequency
            if (viablesteps[rep] % config.TRJEVERY) == 0:
                traj.queue_trj(replicas[rep])
            if (viablesteps[rep] % config.ENEEVERY) == 0:
                traj.queue_ene(replicas[rep])

            # HAVE WE FOUND THE NATIVE YET? If so, stop
            if config.STOPATNATIVE == 1:    
                thisclist  = replicas[rep].chain.contactstate()
                if (nativeclist == thisclist):
                    foundnative = 1


        # After the production cycle,      
        if (prodstep % config.SWAPEVERY) == 0:

            #Before attempting swaps change the restraints 
            replicas = setup_restraints_for_replica(replicas,config)
            ### ...after every production run, attempt a SWAP
            success = attemptswap(replicas,swaps,viableswaps)

            if (VERBOSE):
              if success:   
                outstring = str(prodstep)+'\tSwap successful!\treplica temps: ['
                for rep in range(0,config.NREPLICAS):
                  outstring = outstring + str(replicas[rep].mc.temp) + ' '
                print outstring +']'

              else:  
                print str(prodstep)+'\tUnsuccessful swap attempt.'


        # Print status
        if (prodstep % config.PRINTEVERY) == 0:

            # calc MC move acceptance
            for rep in range(0,len(steps)):
                if steps[rep] > 0:
                    moveacceptance[rep] = float(viablesteps[rep])/float(steps[rep])
                    if viablesteps[rep] > 0:
                        acceptance[rep] = float(acceptedsteps[rep])/float(viablesteps[rep])
                    else:
                        acceptance[rep] = 0.0
                else:
                    moveacceptance[rep] = 0.0
                    acceptance[rep] = 0.0

            # calc replica swap acceptance
            for rep in range(0,len(steps)):
                if swaps[rep] > 0:
                    swap_acceptance[rep] = float(viableswaps[rep])/float(swaps[rep])
                else:
                    swap_acceptance[rep] = 0.0

            # Output the status of the simulation
            print prodstep,'production steps'
            print '%-12s %-12s %-12s %-12s %-12s %-12s %-12s '%('replica','viablesteps','steps','MCaccept','viableswaps','swaps','SWAPaccept')
            for rep in range(0,len(steps)):
                print '%-12d %-12d %-12d %-12s %-12d %-12d %-12s '%(rep,viablesteps[rep],steps[rep],'%1.3f'%moveacceptance[rep],swaps[rep],viableswaps[rep],'%1.3f'%swap_acceptance[rep])
            if config.STOPATNATIVE == 1:
                print 'NATIVE CLIST:', nativeclist
            print '%-8s %-12s %-12s'%('replica','foundnative','contact state') 
            for replica in replicas:
                print '%-8d %-12d %s'%(replica.repnum,foundnative,repr(replica.chain.contactstate()))

        # Continue to the next production cycle!
        prodstep = prodstep + 1

            
    # write the acceptance ratios to the <rundir>/data directory
    faccept = open(config.DATADIR+'/acceptance','w')
    for i in range(0,len(config.REPLICATEMPS)):
        tmp = str(config.REPLICATEMPS[i]) + '\t'
        tmp = tmp +  str(acceptedsteps[i]) + '\t'
        tmp = tmp +  str(viablesteps[i]) + '\t'
        tmp = tmp +  str(steps[i]) + '\n'


    #Output the status of the simulation one last time...
    print prodstep,'production steps'
    print '%-12s %-12s %-12s %-12s %-12s %-12s %-12s '%('replica','viablesteps','steps','MOVEaccept','viableswaps','swaps','SWAPaccept')
    for rep in range(0,len(steps)):
        print '%-12d %-12d %-12d %-12s %-12d %-12d %-12s '%(rep,viablesteps[rep],steps[rep],'%1.3f'%moveacceptance[rep],swaps[rep],viableswaps[rep],'%1.3f'%swap_acceptance[rep])
    if config.STOPATNATIVE == 1:
        print 'NATIVE CLIST:', nativeclist
    print '%-8s %-12s %-12s'%('replica','foundnative','contact state')    
    for replica in replicas:
        print '%-8d %-12d %s'%(replica.repnum,foundnative,repr(replica.chain.contactstate()))

    faccept.write(tmp)  
    faccept.close()


    
    # write the last of the trj and ene buffers
    # and close all the open trajectory file handles
    traj.cleanup(replicas)      

