# -*- coding: utf-8 -*-

"""The utilities for ``cp_app``."""

import numpy as np
import yaml

Kb=1.3806504e-23  # Boltzmann constant in [J/K]
Ph=1.98644586e-23 # Planck constant in [J.cm]
Avo=6.02214076e23 # 1/mol
J2cal=0.2390
th2cm=33.35641

def read_vibspectrum(filename):
    frequencies=[]
    with open(filename) as fi:
        for line in fi.readlines()[3:-1]:
            frequencies.append(float(line.strip().split()[2]))
    return np.array(frequencies)

def read_frequencies_from_mesh(filename):
    mesh=yaml.load(open(filename))
    w=[fr["frequency"] for fr in mesh["phonon"][0]["band"]]
    w=np.array(w)*th2cm
    return w

def cv_from_pdos(temp, pdos):
    pdos=pdos[np.where(pdos[:,0]>0)]
    x = Ph * pdos[:,0] / Kb / temp
    expVal = np.exp(x)
    cv_contributions= np.sum(pdos[:,1:],axis=1)* Avo*Kb * x ** 2 * expVal / (expVal - 1.0) ** 2
    return np.sum(cv_contributions)

def cv_from_dos(temp, totaldos):
    dos=totaldos[np.where(totaldos[:,0]>0)]
    x = Ph * dos[:,0] / Kb / temp
    expVal = np.exp(x)
    cv_contributions= dos[:,1]* Avo*Kb * x ** 2 * expVal / (expVal - 1.0) ** 2
    return np.sum(cv_contributions)

def cv_from_frequencies(temp, freqs):
    freqs=freqs[freqs>0]
    x = Ph * freqs / Kb / temp
    expVal = np.exp(x)
    cv_contributions=Avo*Kb * x ** 2 * expVal / (expVal - 1.0) ** 2
    return np.sum(cv_contributions)

def read_totaldos(filename):
    data=np.loadtxt(filename,skiprows=1)
    data[:,0]*=th2cm
    return data

def read_pdos(filename):
    data=np.loadtxt(filename,skiprows=1)
    data[:,0]*=th2cm
    return data


def add_type_label(mydict,atomtype,name,label):
    if atomtype in mydict:
        mydict[atomtype][name]=label
    else:
        mydict[atomtype]={name:label}
    return mydict

def read_atoms_from_mesh(filename):
    mesh=yaml.load(open(filename))
    w=[fr["frequency"] for fr in mesh["phonon"][0]["band"]]
    w=np.array(w)*th2cm
    return w

def cv_from_pdos_site(temp, pdos,site):
    pdos=pdos[np.where(pdos[:,0]>0)]
    x = Ph * pdos[:,0] / Kb / temp
    expVal = np.exp(x)
    cv_contributions= pdos[:,site+1]* Avo*Kb * x ** 2 * expVal / (expVal - 1.0) ** 2
    return np.sum(cv_contributions)

def select_structures(nsamples,df):
    selected=set()
    for structure_type in df["structure_type"].unique():
        selected.add(df.loc[df["structure_type"]==structure_type].index.values[0])
        selected.add(df.loc[df["structure_type"]==structure_type].index.values[1])
        if len(selected)>nsamples-1:
            break
    for atom_type in df["atom_types"].unique():
        selected.add(df.loc[df["atom_types"]==atom_type].index.values[0])
        if len(selected)>nsamples-1:
            break
        
    while len(selected)< nsamples:
        selected.add(df.sample(1).index.values[0])
        
    return selected

