#usage mcmc.py 'initialtheta' 'total' 'save intervals' 'outfilename' 'outfilepath'
#Fit only mean dominance times and CV
import pandas as pd
import numpy as np
from sys import argv
import model as model

def chisq(data,model,variance):
    d = data
    m = model
    v = variance

    #for first epoch probability
    #fmax = lambda x:(max(x,0.001)**2)
    #N1 = m['COUNTS_U1']
    #x1 = (m['FIRST_EPOCH_PROB_U1']-d['FIRST_EPOCH_PROB_U1'])**2
    #s1 = m['FIRST_EPOCH_PROB_U1'].map(fmax)
    #N2 = m['COUNTS_U2']
    #x2 = (m['FIRST_EPOCH_PROB_U2']-d['FIRST_EPOCH_PROB_U2'])**2
    #s2 = m['FIRST_EPOCH_PROB_U2'].map(fmax)

    #other chi follow form: chi = ((d-m)**2)/(v**2)
    chi = pd.DataFrame()
    #chi['FIRST_EPOCH_PROB_U1'] = (N1*(x1/s1) + N2*(x2/s2))
    chi['MEAN_U1'] = ((m['MEAN_U1']-d['MEAN_U1'])**2)/(v['MEAN_U1']**2)
    chi['MEAN_U2'] = ((m['MEAN_U2']-d['MEAN_U2'])**2)/(v['MEAN_U2']**2)
    chi['SD_U1'] = ((m['SD_U1']-d['SD_U1'])**2)/(v['SD_U1']**2)
    chi['SD_U2'] = ((m['SD_U2']-d['SD_U2'])**2)/(v['SD_U2']**2)

    return chi

def chisum(df):
    #return scalar
    one = df.sum(axis=0)
    two = one.sum(0)
    return two

def chisave(df):
    s1 = df[df.columns[0]]
    s2 = df[df.columns[1]]
    s3 = df[df.columns[2]]
    s4 = df[df.columns[3]]
    #s5 = df[df.columns[4]]
    #save = [s1,s2,s3,s4,s5]
    save = [s1,s2,s3,s4]
    return save

def taupen(tau_a):
    penalty = (tau_a-2)**2
    return penalty

def kpen(ks,kc):
    delta = 0
    pen = delta*(ks + kc)**2
    return pen

def mapkc(dic,new):
    a = new
    s1 = dic
    for key, value in s1.items():
        # do something with value
        s1[key] = value + a*value*np.random.randn()
    return s1

def save(lis1,lis2,lis3,x,n,name):
    interval = str(n)
    basename = name
    fn = basename+interval
    indx = x

    thetalist = lis1
    chilist = lis2
    sumlist = lis3

    f1 = str(fn+"THETA.pkl")
    d1 = pd.DataFrame(thetalist,index=indx)
    d1.to_pickle(f1)

    f2 = str(fn+"CHISQ.pkl")
    s1 = pd.Series(chilist,index=indx)
    s1.to_pickle(f2)

    f3 = str(fn+"SUM.pkl")
    d2 = pd.DataFrame(sumlist,index=indx)
    d2.to_pickle(f3)
    print "steps written to file:",n

def newtheta(theta):
    a = 0.05
    fa = lambda x:x + a*x*np.random.randn()
    current = theta
    #newkc = mapkc(current.KC.copy())
    #newsig = mapkc(current.sig.copy())
    proposed = current.map(fa)
    #proposed = current.drop(['KC','sig']).map(fa)
    #proposed['KC']=newkc
    #proposed['sigma']=newsig
    return proposed

init_guess = str(argv[1])
burntime = int(argv[2])
savestep = int(argv[3])
outfile = str(argv[4])
savepoints = range(savestep,burntime,savestep)
pathname = str(argv[5])
outfilename = pathname+outfile

thetafirst = pd.read_pickle(init_guess)
theta_cur = pd.Series(thetafirst)
theta_prop = newtheta(theta_cur)
print theta_cur
print theta_prop

index = [0,20,40,60,80]
cv = 0.6
cv_var = 0.1
p1_var = 0.1
fcv = lambda x: x*cv
fcv_var = lambda x:(x*cv_var)
fp1_var = lambda x:(x*p1_var)

#load data
syprob = pd.read_csv('sy_prob.csv',index_col=0)
syU1 = pd.read_csv('syU1.csv',index_col=0)
syU2 = pd.read_csv('syU2.csv',index_col=0)

data = pd.DataFrame()
data['COUNTS_U1'] = syprob['counts_U1']
data['COUNTS_U2'] = syprob['counts_U2']
data['FIRST_EPOCH_PROB_U1'] = syprob['U1']
data['FIRST_EPOCH_PROB_U2'] = syprob['U2']
data['MEAN_U1'] = syU1['MEAN']
data['MEAN_U2'] = syU2['MEAN']
data['SD_U1'] = data['MEAN_U1'].map(fcv)
data['SD_U2'] = data['MEAN_U2'].map(fcv)
data.index=index
col = data.columns

variance = pd.DataFrame()
variance['COUNTS_U1'] = data['COUNTS_U1'].map(fp1_var)
variance['COUNTS_U2'] = data['COUNTS_U2'].map(fp1_var)
variance['FIRST_EPOCH_PROB_U1'] = data['FIRST_EPOCH_PROB_U1'].map(fp1_var)
variance['FIRST_EPOCH_PROB_U2'] = data['FIRST_EPOCH_PROB_U2'].map(fp1_var)
variance['MEAN_U1'] = data['MEAN_U1'].map(fcv_var)
variance['MEAN_U2'] = data['MEAN_U2'].map(fcv_var)
variance['SD_U1'] = data['SD_U1'].map(fcv_var)
variance['SD_U2'] = data['SD_U2'].map(fcv_var)
variance.index = index
variance.columns = col

print"this is the data:"
print data
print "this is variance:"
print variance

print"Fitting only mean dominance times and CV"

print"running current theta"
sustained = 1
onset = 2
m1 = model.run(theta_cur,sustained,onset)
m1.index = index
m1.columns = col

#calculate chi
chicur = pd.DataFrame(chisq(data,m1,variance),index=index)
sumchicur = chisum(chicur)

print"here is Y1:"
print m1
print"here is sum of chi squared:"
print sumchicur
print"all chi values:"
print chicur

THETA = []
CHI = []
SUM = []
x = []
n = -1

while n < burntime:
    n +=1
    print n

    THETA.append(theta_cur)
    CHI.append(chisave(chicur))
    SUM.append(chicur.sum(axis=0))
    x.append(n)

    m2 = model.run(theta_prop,sustained,onset)
    m2.index = index
    m2.columns = col

    chiprop = pd.DataFrame(chisq(data,m2,variance),index=index)
    #chipenalty = taupen(theta_prop['tau_a'])#+kpen(theta_prop['KS'],theta_prop['KC'])
    sumchiprop = chisum(chiprop)

    ratio = np.exp((-sumchiprop+sumchicur)/2)
    random_min = np.random.sample()
    keep = random_min < ratio

    if keep:
        theta_cur = theta_prop
        chicur = chiprop
        sumchicur = sumchiprop
        print "new theta:",theta_cur
        print "model results:"
        print m2
        print"new sum chi squared:",sumchicur
        print"all chi values:"
        print chicur
        theta_prop = newtheta(theta_cur)

    else:
        theta_prop = newtheta(theta_cur)

    if n in savepoints:
        save(THETA,CHI,SUM,x,n,outfilename)
    else:
        pass

save(THETA,CHI,SUM,x,n,outfilename)

print"Fitting only mean dominance times and CV"
print "done"

#file1 = str("THETA"+outfile+".pkl")
#df1 = pd.DataFrame(THETA,index=x)
#df1.to_pickle(file1)

#file2 = str("CHISQ"+outfile+".pkl")
#S1 = pd.Series(CHI,index=x)
#S1.to_pickle(file2)

#to unpickle CHI
#for item in x:
    #dfi = pd.DataFrame(CHI[i])
    #print dfi.T

#file3 = str("SUM"+outfile+".pkl")
#df2 = pd.DataFrame(SUM,index=x)
#df2.to_pickle(file3)
