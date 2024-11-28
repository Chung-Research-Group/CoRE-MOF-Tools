# -*- coding: utf-8 -*-

"""The methods to compute labels for machine learning the heat capacity."""

from ase.io import read
from .utils import  cv_from_pdos_site 
import pandas as pd
import numpy as np
import phonopy
from phonopy.units import CP2KToTHz

def compute_total_dos_structure(phonopy_params: str, unitfactor=CP2KToTHz, dx:float=0.5, fmax:float=100.0, freq_pitch:float=0.5, saveto: str=None):
    """Compute projected dos from phonopy parameter file

    :param phonopy_params: list of phonopy parameter files (output of DFT)
    :param cif: list of crystal structure in cif format
    :param temperatures: the target temperature 
    :param factor: the unit conversion factor
    :param dx: spacing to compute dos
    :param fmax: max frequency in dos calculations
    :param freq_pitch: pitch frequency in dos calculations
    """
    phonon = None
    phonon = phonopy.load(phonopy_params, factor=unitfactor)
    phonon.run_mesh([1,1,1], with_eigenvectors=True)
    # total dos
    phonon.run_total_dos(sigma=dx,freq_min=0,freq_max=fmax, freq_pitch=freq_pitch)
    dos_dict = phonon.get_total_dos_dict()
    totaldos = np.vstack([dos_dict["frequency_points"],dos_dict["total_dos"]]).T
    if saveto:
        np.savetxt(saveto,totaldos,header="sigma = {}".format(dx))

    return totaldos


def compute_projected_dos_structure(phonopy_params: str, unitfactor=CP2KToTHz, dx=0.5, fmax=100, freq_pitch=0.5, saveto: str=None):
    """Compute projected dos from phonopy parameter file

    :param phonopy_params: list of phonopy parameter files (output of DFT)
    :cif: list of crystal structure in cif format
    :temperatures: the target temperature 
    :factor: the unit conversion factor
    :dx: spacing to compute dos
    :fmax: max frequency in dos calculations
    :freq_pitch: pitch frequency in dos calculations
    :saveto: save the projected dos to a file
    """
    phonon = None
    phonon = phonopy.load(phonopy_params, factor=unitfactor)
    phonon.run_mesh([1,1,1], with_eigenvectors=True)
    mesh = phonon.get_mesh_dict()
    # pdos
    phonon.run_projected_dos(sigma=dx,freq_min=0,freq_max=fmax, freq_pitch=freq_pitch)
    pdos_dict = phonon.get_projected_dos_dict()
    pdos = np.vstack([pdos_dict["frequency_points"],pdos_dict["projected_dos"]]).T
    if saveto:
        np.savetxt(saveto, pdos, header="sigma = {}".format(dx))

    return pdos 

def compute_atomic_cv_dataset(phonopy_params: list[str], cifs: list[str], temperatures: list, verbos=False, saveto: str="labels.csv") -> pd.DataFrame:
    """Compute atomic contribution to total cv from phonopy parameter file

    :param phonopy_params: list of phonopy parameter files (output of DFT)
    :cifs: list of crystal structure in cif format
    :temperatures: the target temperature 
    """
    labels={}
    for phonopy_param,cif in zip(phonopy_params,cifs):
        if verbos:
            print(cif)
        pdos=compute_projected_dos_structure(phonopy_params=phonopy_param, saveto="%s_projected_dos.dat".format(cif.replace(".cif","")))
        pdos[:,0]*=th2cm
        atoms=read(cif)
        if not atoms.get_global_number_of_atoms()==pdos.shape[1]-1:
            print(atoms.get_global_number_of_atoms(),pdos.shape[0])
            print("Warning! number of atoms do not match in pdos and structure for %s"%cif)
            continue
        for atomidx in range(atoms.get_global_number_of_atoms()):
            site_name="%s_%i"%(cif.replace(".cif",""),atomidx)
            for i,temperature in enumerate(temperatures):
                cv_site=cv_from_pdos_site(temperature,pdos,atomidx)
                if i==0:
                    labels[site_name]={"pCv_%05.2f"%temperature:cv_site}
                else:
                    labels[site_name]["pCv_%05.2f"%temperature]=cv_site

    df=pd.DataFrame.from_dict(labels).T
    df.to_csv(saveto)

    return df
