###usage as module to run two pool neural activity rate model
import pandas as pd
import numpy as np
import rivalry as rv

def run(parameters,ntrials,onsettrial):
    """use argument parameters to solve model"""
    theta = {}
    theta = parameters

    condition = [0,1,2,3,4]
    #SUSTAINED = np.arange(ntrials)
    ONSET = np.arange(onsettrial)

    dt = 1
    dom_ms = 500000
    onset_ms = 6000

    prestim = 6000

    PERCENT_IMBALANCE = {}
    PERCENT_IMBALANCE={n:[] for n in condition}

    for imbalance in PERCENT_IMBALANCE.keys():

        r = imbalance

        DOM1 = []
        DOM2 = []
        LEADER = []

        run_dom = {}
        U1d = []
        U2d = []
        Td = []

        run_dom = rv.dynamics(theta,r,dom_ms,prestim,dt,full=True)
        U1d = np.asarray(run_dom.pop('u1'))
        U2d = np.asarray(run_dom.pop('u2'))
        Td = np.asarray(run_dom.pop('t'))

        run_dom.clear()

        td1 = rv.dominance(Td,U1d,U2d,prestim)
        td2 = rv.dominance(Td,U2d,U1d,prestim)

        DOM1.extend(td1)
        DOM2.extend(td2)

        TD1 = np.array(DOM1)
        TD2 = np.array(DOM2)

        for trial in ONSET:

            run_onset = {}
            U1p = []
            U2p = []
            Tp = []

            run_onset = rv.dynamics(theta,r,onset_ms,prestim,dt,full=True)

            U1p = np.asarray(run_onset.pop('u1'))
            U2p = np.asarray(run_onset.pop('u2'))
            Tp = np.asarray(run_onset.pop('t'))

            run_onset.clear()

            first = rv.first_epoch(U1p,U2p,Tp,prestim,dt)
            LEADER.append(first)
        
        count1 = LEADER.count(1)
        count2 = LEADER.count(2)

        FINAL_DATA = {}
        FINAL_DATA['FIRST_EPOCH_PROB_U1'] = rv.prob_seq(LEADER,1,2)
        FINAL_DATA['FIRST_EPOCH_PROB_U2'] = rv.prob_seq(LEADER,2,1)
        FINAL_DATA['COUNTS_U1'] = count1
        FINAL_DATA['COUNTS_U2'] = count2
        FINAL_DATA['MEAN_U1'] = np.mean(TD1)
        FINAL_DATA['MEAN_U2'] = np.mean(TD2)
        FINAL_DATA['SD_U1'] = np.std(TD1)
        FINAL_DATA['SD_U2'] = np.std(TD2)

        PERCENT_IMBALANCE[imbalance] = FINAL_DATA

    model_data = pd.DataFrame(PERCENT_IMBALANCE).T
    model_data.index = [0,20,40,60,80]

    return model_data
