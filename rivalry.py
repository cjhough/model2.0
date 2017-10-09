##usuage as module with useful functions for rivalry dynamics,
##dominance durations, first epoch, etc.

import pandas as pd
import numpy as np

def gain_pl_sqrt(x):
    """gain function"""
    y = np.sqrt((x>0)*x)
    return y

def stimulus(k,c):
    """divide input stimulus"""
    
    ks = k
    kc = c
    k = (ks + (1+kc))*0.5
    return k

def dynamics(d,k,c,ton,toff,dx,full=None):
    """returns solutions for system variables
    and time as arrays in key,value pairs"""
    
    #unpack fixed parameters
    param = d
    beta = param['beta']
    gamma = param['gamma']
    tau_a = param['tau_a']
    tau_u = param['tau_u']
    mu = param['mu']
    sigma = param['sigma']
    #b = param['b']
    
    #input stimulus
    ks = k
    kc = c
    k1 = stimulus(ks,kc)
    k2 = stimulus(ks,kc)#weaker stimulus
    
    dt = dx
    tf = ton
    equil = toff
    totaltime= tf + equil
    
    #system variables
    U = np.random.uniform(0,1,(2,2))
    A = np.random.uniform(0,1,(2,1))
    Z = np.random.uniform(0,1,(2,1))
    
    u1 = []
    u2 = []
    a1 = []
    a2 = []
    z1 = []
    z2 = []
    t = []
    i=-1
    
    while i < totaltime:
        i += 1
        t.append(i*dt)
        stim_on = int(equil<i)
        S1 = k1*stim_on
        S2 = k2*stim_on
        if full==True:
            Z[0] = sigma*np.random.randn()*(1+mu*U[0][0])
            Z[1] = sigma*np.random.randn()*(1+mu*U[1][0])

            U[0][1] = U[0][0] +dt*(gain_pl_sqrt(S1 - beta*A[1]*U[1][0] + Z[0]))/tau_u
            U[1][1] = U[1][0] +dt*(gain_pl_sqrt(S2 - beta*A[0]*U[0][0] + Z[1]))/tau_u
            U[0][0]= U[0][1]
            U[1][0]= U[1][1]

            A[0] += dt*(1 - A[0] - gamma*A[0]*U[0][1])/tau_a
            A[1] += dt*(1 - A[1] - gamma*A[1]*U[1][1])/tau_a
            
            u1.append(U[0][0])
            u2.append(U[1][0])
            a1.append(A[0])
            a2.append(A[1])
            z1.append(Z[0])
            z2.append(Z[1])
        
        else:
            Z[0] = 1+mu*U[0][0]
            Z[1] = 1+mu*U[1][0]
            
            U[0][1] = gain_pl_sqrt(S1 - beta*A[1]*U[1][0] + Z[0])/tau_u
            U[1][1] = gain_pl_sqrt(S2 - beta*A[0]*U[0][0] + Z[1])/tau_u
            U[0][0]= U[0][1]
            U[1][0]= U[1][1]

            A[0] += dt*(1 - A[0] - gamma*A[0]*U[0][1])/tau_a
            A[1] += dt*(1 - A[1] - gamma*A[1]*U[1][1])/tau_a
            
            u1.append(U[0][0])
            u2.append(U[1][0])
            a1.append(A[0])
            a2.append(A[1])
            z1.append(Z[0])
            z2.append(Z[1])

    solution = {}
    solution['t'] = t
    solution['u1'] = u1
    solution['u2'] = u2
    solution['a1'] = a1
    solution['a2'] = a2
    solution['z1'] = z1
    solution['z2'] = z2

    return solution

def dominance(N,xi,xj,pre=None):
    """return list of dominance time intervals"""
    if pre==None:
        n = N
        i = xi
        j = xj
    else:
        n = N[pre:]
        i = xi[pre:]
        j = xj[pre:]
    
    length_of_time = range(len(n))
    t_n = [increment for increment in length_of_time]
    temp = []
    TD = []
    index=np.greater(i,j)
    
    if np.all(index):
        #all i > j == True
        print "WARNING: NO rivalry dynamics detected"
        TD.append(len(t_n))
        return TD
    else:
        sd = pd.Series(t_n,index=index)
        dom = np.asarray(sd.loc[False])
    
        if dom.ndim == 0:
            print "WARNING: NO rivalry dynamics detected"
            TD.append(len(t_n))
            return TD
        else:
            #dom.ndim=1
            temp = list(dom)
            #convert to list in order to delete elements
            #np.any(temp[:]!=templist[:])
            #should always be the same
            while len(temp) >= 2:
                t2 = temp[-1]
                t1 = temp[-2]
                interval = t2-t1

                if interval <= 1:
                    del temp[-1]
                else:
                    TD.append(interval)
                    del temp[-1]
            
            if len(TD) >= 3:
                temp = []
                del TD[0]
                del TD[-1]
                return TD
            else:
                if len(TD) != 0:
                    temp = []
                    print "rivalry:warning: dom found < 3"
                    return TD
                else:
                    TD.append(1) #small but not zero
                    print "WARNING: dom detected < 1"
                    return TD
                
def first_epoch(x,y,t,pre_t,dt):
    t1 = pre_t
    start = int(1/dt)
    a=[]
    b=[]
    c=[]
    
    a = x[t1:]
    b = y[t1:]
    c = t[t1:]
    
    #find mean of background noise
    bkgnd = np.append(x[start:t1],y[start:t1],axis=0)
    meanbkgnd = np.mean(bkgnd)
    
    #filter for signal above background
    u1corr = np.where(a>meanbkgnd,a,0)
    u2corr = np.where(b>meanbkgnd,b,0)
    
    #set up time for finding intervals of percept dominance
    tz = range(len(c))
    maxu = np.amax(np.append(u2corr,u1corr,axis=0))
    diff = np.absolute(np.subtract(u1corr,u2corr))
    threshold = maxu*0.6
    tfilter = np.where(diff>=threshold,tz,0)
    tsort = [i for i in tfilter if i!=0]
    
    #constrains first epoch to be at least 50ms in duration
    duration = 49
    start=0
    jump=duration
   
    #default if first percept is not found
    default = 0
    u1first=1
    u2first=2
    firstpercept = default
    n=-1
    stop = len(tsort)
    winner=True
    
    while n<stop:
        n=+1
    
        if stop<jump+1:
            n = stop
        else:
            t1 = tsort[start]
            t2 = tsort[jump]
            dom = u1corr[t1]>u2corr[t1]
            if t2-t1 == duration:
                if winner == dom:
                    firstpercept = u1first
                    n = stop
                else:
                    firstpercept = u2first
                    n = stop
            else:
                if t2-t1 > duration:
                    start = start + 1
                    jump = jump + 1
                else:
                    n=stop

    record_winner = firstpercept
    return record_winner

def prob_seq(s,e):
    """degree of bias in probability of being first percept"""
    win_results = s
    winner = float(win_results.count(e))
    loser = float(win_results.count(2))
    count = onseterror(win_results)
    if count == True:
        if float(len(win_results))-float(win_results.count(0)) != 0:
            outcomes = winner+loser
            prob = (winner/outcomes)
        else:
            outcomes = float(len(win_results))
            prob = (winner/outcomes)
        return prob
    else:
        prob = 0.0
        return prob


def onseterror(lis):
    """onset rivarly check"""
    listrials = lis
    total = len(listrials)
    accept = (total/2)
    fails = listrials.count(0)
    passed = True
    failed = False
    if fails < accept:
        return passed
    else:
        return failed


def no_input(x,d,t,dx):
    """solve system without input stimulus
    and return last value for each variable"""
    varlist = x
    U1 = varlist[0]
    U2 = varlist[1]
    A1 = varlist[2]
    A2 = varlist[3]
    Z1 = varlist[4]
    Z2 = varlist[5]

    param = d
    beta = param['beta']
    gamma = param['gamma']
    sigma = param['sigma']
    tau_a = param['tau_a']
    dt = dx

    runtime = t
    S1 = 0
    S2 = 0

    i = -1

    while i < runtime:
        i += 1
        U1 = gain_pl_sqrt(S1-beta*U2-gamma*A1+Z1)
        U2 = gain_pl_sqrt(S2-beta*U1-gamma*A2+Z2)
        A1 = A1 +(dt/tau_a)*(-A1+U1)
        A2 = A2 +(dt/tau_a)*(-A2+U2)
        Z1 = sigma*np.random.randn()
        Z2 = sigma*np.random.randn()

    final = [U1,U2,A1,A2,Z1,Z2]
    return final