from numpy import exp, sin,cos, linspace
import matplotlib.pyplot as plt
import os, time, glob
import math
import numpy as np
import pickle

def get_labels(idx,Mol):
    LBL=["R1","R2","R3","R4","R5","R6"]
    if Mol=="HFCO":
        LBL=["CH","CF","CO","HCO","FCO","Phi"]
    if Mol=="HONO":
        LBL=["N=O","O-N","OH","ONO","HON","Phi"]
    if Mol=="CS2":
        LBL=["CS1","CS2","Phi"]

    return LBL[idx]


def pes_6D(hono_r,hono_th,t,t1,ith,ith1,ms):

    
    #r = pickle.load(open('r_params.pkl','rb'))
    #wu = pickle.load(open('w_params.pkl','rb'))
    if ms=="HFCO":
        r = pickle.load(open('r_hfco.pkl','rb'))
        wu = pickle.load(open('w_hfco.pkl','rb'))
        c = 0.219898906 # for HFCO Full 6D PES
    if ms=="HONO":
        r = pickle.load(open('r_hono.pkl','rb'))
        wu = pickle.load(open('w_hono.pkl','rb'))
        c = 0.065600224 # for HONO 6D PES        
    if ms=="CS2":
        r = pickle.load(open('r_cs2.pkl','rb'))
        wu = pickle.load(open('w_cs2.pkl','rb'))
        c = 0.272084269200190 # for CS2 3D PES
    #r=np.loadtxt("r_params.txt")
    #wu=np.loadtxt("w_params.txt")
    #c = 0.219898906 # for HFCO Cis
    
    nofd = len(wu[0]) #6 #Number of degrees of freedom
    nofn = len(r) #100 #Number of Neurons
    H_cm=219474.63 #Hartree to cm-1

    v = []
    #coord = [3.446,2.488,2.4714,0.8619,-0.2519,3.14285] # TRANS HFCO minimum
    #coord = [2.0632,2.534,2.22874,-0.6138,-0.54038,3.14285] #coordinates of EQ HFCO
    coord=[]
    RR=hono_r.split(",")
    RT=hono_th.split(",")
    for ia in RR:
        coord.append(float(ia))
    for ib in RT:
        coord.append(float(ib))
        
    #coord[ith]=t
    coord= coord[:ith]+[t]+coord[ith+1:]
    coord= coord[:ith1]+[t1]+coord[ith1+1:]
    
    for i in range(nofn):
        v.append("0")

    for i in range(nofn):
        v[i] = 1.0
        for j in range(nofd):
            v[i] = v[i]*exp(wu[i][j]*coord[j])

    V = 0.0
    for i in range(nofn):
        V = V + r[i]*v[i] 
    return H_cm*(V+c)



def compute(R,Q, bstr0, bstr1,Mo, resolution=20):
    
    b0=float(bstr0.split(",")[0])
    w0=float(bstr0.split(",")[1])
    T0=int(bstr0.split(",")[2])
    Tf=int(bstr0.split(",")[3])
    b1=float(bstr1.split(",")[0])
    w1=float(bstr1.split(",")[1])
    T1=int(bstr1.split(",")[2])
    Tf1=int(bstr1.split(",")[3])
    
    """Return filename of plot of the NN Fitted Molecular PES."""
    t = linspace(b0, w0, T0)
    t1 = linspace(b1, w1, T1)
    if Tf==Tf1:
        rr=get_labels(Tf,Mo)
        plt.figure()
        plt.title('%s 1D-cut PES'%Mo)
        plt.xlabel("%s"%rr)
        plt.ylabel("V(R$_{%s}$), cm$^{-1}$"%rr)
        #plt.plot(t, hfco_pes(R,Q,t,Tf),'-or')
        if Mo=="CS2" or Mo=="HONO" or Mo=="HFCO":
            plt.plot(t, pes_6D(R,Q,t,t,Tf,Tf,Mo),'-or')
    else:
        
        Z=[]
        for i in t:
            Z0=[]
            for j in t1:
                if Mo=="CS2" or Mo=="HONO" or Mo=="HFCO":
                    Z0.append(pes_6D(R,Q,i,j,Tf,Tf1,Mo))
            Z.append(Z0)
        Z=np.array(Z)

        [X, Y] = np.meshgrid(t, t1)

        rr=get_labels(Tf,Mo)
        rr0=get_labels(Tf1,Mo)
        fig,ax=plt.subplots(1,1)
        ax = plt.axes(projection='3d')
        ax.set_title('%s Contour PES'%Mo)
        ax.set_xlabel('%s'%rr)
        ax.set_ylabel('%s'%rr0)
        plt.contour(X, Y, Z, 100, cmap='RdGy');
        plt.colorbar();

    if not os.path.isdir('static'):
        os.mkdir('static')
    else:
        # Remove old plot files
        for filename in glob.glob(os.path.join('static', '*.png')):
            os.remove(filename)
    # Use time since Jan 1, 1970 in filename in order make
    # a unique filename that the browser has not chached
    plotfile = os.path.join('static', str(time.time()) + '.png')
    plt.savefig(plotfile)
    return plotfile
