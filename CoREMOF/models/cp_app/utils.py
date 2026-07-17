# -*- coding: utf-8 -*-

"""The utilities for ``cp_app``."""

import numpy as np
import yaml

Kb=1.3806504e-23  # Boltzmann constant in [J/K]
Ph=1.98644586e-23 # Planck constant in [J.cm]
Avo=6.02214076e23 # 1/mol
J2cal=0.2390
th2cm=33.35641


def _heat_capacity_factor(x):
    """Return ``x**2 exp(x) / (exp(x)-1)**2`` without overflow."""

    exp_negative = np.exp(-x)
    denominator = -np.expm1(-x)
    return x**2 * exp_negative / denominator**2

def read_vibspectrum(filename):
    frequencies=[]
    with open(filename) as fi:
        for line in fi.readlines()[3:-1]:
            frequencies.append(float(line.strip().split()[2]))
    return np.array(frequencies)

def read_frequencies_from_mesh(filename):
    with open(filename, encoding="utf-8") as stream:
        mesh = yaml.safe_load(stream)
    w=[fr["frequency"] for fr in mesh["phonon"][0]["band"]]
    w=np.array(w)*th2cm
    return w

def cv_from_pdos(temp, pdos):
    pdos=pdos[np.where(pdos[:,0]>0)]
    x = Ph * pdos[:,0] / Kb / temp
    cv_contributions = np.sum(pdos[:, 1:], axis=1) * Avo * Kb * _heat_capacity_factor(x)
    return np.sum(cv_contributions)

def cv_from_dos(temp, totaldos):
    dos=totaldos[np.where(totaldos[:,0]>0)]
    x = Ph * dos[:,0] / Kb / temp
    cv_contributions = dos[:, 1] * Avo * Kb * _heat_capacity_factor(x)
    return np.sum(cv_contributions)

def cv_from_frequencies(temp, freqs):
    freqs=freqs[freqs>0]
    x = Ph * freqs / Kb / temp
    cv_contributions = Avo * Kb * _heat_capacity_factor(x)
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
    with open(filename, encoding="utf-8") as stream:
        mesh = yaml.safe_load(stream)
    w=[fr["frequency"] for fr in mesh["phonon"][0]["band"]]
    w=np.array(w)*th2cm
    return w

def cv_from_pdos_site(temp, pdos,site):
    pdos=pdos[np.where(pdos[:,0]>0)]
    x = Ph * pdos[:,0] / Kb / temp
    cv_contributions = pdos[:, site + 1] * Avo * Kb * _heat_capacity_factor(x)
    return np.sum(cv_contributions)

def select_structures(nsamples, df, random_state=None):
    if nsamples < 0:
        raise ValueError("nsamples must be non-negative")
    if nsamples > len(df):
        raise ValueError(
            f"Cannot select {nsamples} unique structures from a dataframe with {len(df)} rows"
        )
    if nsamples == 0:
        return set()

    selected = []
    selected_set = set()

    def add(index):
        if len(selected) < nsamples and index not in selected_set:
            selected.append(index)
            selected_set.add(index)

    for structure_type in df["structure_type"].unique():
        candidates = df.loc[df["structure_type"] == structure_type].index.values
        for candidate in candidates[:2]:
            add(candidate)
        if len(selected) >= nsamples:
            return set(selected)
    for atom_type in df["atom_types"].unique():
        add(df.loc[df["atom_types"] == atom_type].index.values[0])
        if len(selected) >= nsamples:
            return set(selected)
        
    remaining = df.index.difference(selected)
    needed = nsamples - len(selected)
    if needed:
        for index in df.loc[remaining].sample(n=needed, random_state=random_state).index:
            add(index)
        
    return set(selected)
