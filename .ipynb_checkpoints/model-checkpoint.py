###usage as module to run two pool neural activity rate model
import pandas as pd
import numpy as np
import rivalry as rv

def run(parameters,ntrials,onsettrial):
    """use argument parameters to solve model"""
    theta = {}
    theta = parameters

    condition = [1,2,3,4,5]
    SUSTAINED = np.arange(ntrials)
    ONSET = np.arange(onsettrial)

    dt = 1
    sustainsec = 1000000
    onsetsec = 3000

    prestim = 6000

    PERCENT_IMBALANCE = {}
    PERCENT_IMBALANCE={n:[] for n in condition}
    
    for imbalance in PERCENT_IMBALANCE.keys():
    
        r = imbalance
    
        DOM1 = []
        DOM2 = []
        LEADER = []
    
        run = {}
        run_onset = {}
        U1d = []
        U2d = []
        Td = []
        
        run = rv.dynamics(theta,r,sustainsec,prestim,dt,full=True)
        U1d = np.asarray(run.pop('u1'))
        U2d = np.asarray(run.pop('u2'))
        Td = np.asarray(run.pop('t'))
        
        run.clear()
        
        td1 = rv.dominance(Td,U1d,U2d,prestim)
        td2 = rv.dominance(Td,U2d,U1d,prestim)
        
        DOM1.extend(td1)
        DOM2.extend(td2)
        
        TD1 = np.array(DOM1)
        TD2 = np.array(DOM2)
        
        for trial in ONSET:
            U1p = []
            U2p = []
            Tp = []

            run_onset = rv.dynamics(theta,r,onsetsec,prestim,dt,full=True)
            
            U1p = np.asarray(run_onset.pop('u1'))
            U2p = np.asarray(run_onset.pop('u2'))
            Tp = np.asarray(run_onset.pop('t'))

            run_onset.clear()
        
            first = rv.first_epoch(U1p,U2p,Tp,prestim,dt)
            LEADER.append(first)
    
        FINAL_DATA = {}
        FINAL_DATA['FIRST_EPOCH_PROB'] = rv.prob_seq(LEADER,1)
        FINAL_DATA['MEAN_U1'] = np.mean(TD1)
        FINAL_DATA['MEAN_U2'] = np.mean(TD2)
        FINAL_DATA['SD_U1'] = np.std(TD1)
        FINAL_DATA['SD_U2'] = np.std(TD2)
        
        PERCENT_IMBALANCE[imbalance] = FINAL_DATA

    model_data = pd.DataFrame(PERCENT_IMBALANCE).T
    model_data.index = [0,20,40,60,80]

    return model_data
