"""Geometery properties calculation based on Zeo++ software
You need install Zeo++ package independently from source (https://www.zeoplusplus.org/download.html) or conda (https://anaconda.org/conda-forge/zeopp-lsmo)
Befor run this class please test "network" commond is works or not.
"""

import os
import subprocess

def ChanDim(structure, probe_radius = 0, high_accuracy = True, prefix="tmp_chan"):

    """Analysis dimension of channel.

    Args:
        structure (str): path to your CIF.
        probe_radius (float): probe of radiu.
        high_accuracy (bool): use high accuracy or not.
        prefix (str): temporary file.

    Returns:
        Dictionary:
            -   unit by ["unit"], always nan.
            -   dimention by ["Dimention"].
    """
    
    results_chan = {}
    results_chan["unit"]="nan"
    
    tmp_file = f"{prefix}.txt"
    
    if high_accuracy:
        cmd = f'network -ha -chan {probe_radius} {tmp_file} {structure}'
        # cmd = f'network -ha S50 -chan {probe_radius} {tmp_file} {structure}'
    else:
        cmd = f'network -chan {probe_radius} {tmp_file} {structure}'
    _ = subprocess.run(
                        cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=True,
                    )

    with open(tmp_file) as f:
        for i, row in enumerate(f):
            if i == 0:
                dim = int(row.split('dimensionality')[1].split()[0])

    results_chan["Dimention"] = dim

    os.remove(tmp_file)

    return results_chan

def FrameworkDim(structure, high_accuracy = True, prefix="tmp_strinfo"):

    """Analysis dimension of framework.

    Args:
        structure (str): path to your CIF.
        high_accuracy (bool): use high accuracy or not.
        prefix (str): temporary file.

    Returns:
        Dictionary:
            -   unit by ["unit"], always nan
            -   dimention by ["Dimention"]
            -   number of 2D framewor by ["N_1D"]
            -   number of 1D framework by ["N_2D"]
            -   number of 3D framewor by ["N_3D"]
    """

    results_strinfo = {}
    results_strinfo["unit"]="nan"

    tmp_file = f"{prefix}.txt"
    
    if high_accuracy:
        cmd = f'network -ha -strinfo {tmp_file} {structure}'
        # cmd = f'network -ha S50 -strinfo {tmp_file} {structure}'
    else:
        cmd = f'network -strinfo {tmp_file} {structure}'
    _ = subprocess.run(
                        cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=True,
                    )

    with open(tmp_file) as f:
        line = f.readline().split()
        try:
            dim = int(line[-1])
            one_dim = int(line[7])
            two_dim = int(line[8])
            three_dim = int(line[9])
        except:
            one_dim = 0
            two_dim = 0
            three_dim = 0
            dim = 0
    results_strinfo["Dimention"] = dim
    results_strinfo["N_1D"] = one_dim
    results_strinfo["N_2D"] = two_dim
    results_strinfo["N_3D"] = three_dim

    os.remove(tmp_file)

    return results_strinfo

def PoreDiameter(structure, high_accuracy = True, prefix="tmp_pd"):

    """Analysis pore diameter of structure.

    Args:
        structure (str): path to your CIF.
        high_accuracy (bool): use high accuracy or not.
        prefix (str): temporary file.

    Returns:
        Dictionary:
            -   unit by ["unit"], always angstrom, Å
            -   largest cavity diameter by ["LCD"]
            -   pore-limiting diameter by ["PLD"]
            -   largest free pore diameter by ["LFPD"]
    """

    results_pd = {}
    results_pd["unit"]="angstrom, Å"

    tmp_file = f"{prefix}.txt"
    
    if high_accuracy:
        cmd = f'network -ha -res {tmp_file} {structure}'
        # cmd = f'network -ha S50 -res {tmp_file} {structure}'
    else:
        cmd = f'network -res {tmp_file} {structure}'
    _ = subprocess.run(
                        cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=True,
                    )
    with open(tmp_file) as f:
        line = f.readline().split()
        results_pd["LCD"], results_pd["PLD"], results_pd["LFPD"] = map(float, line[1:4])
    os.remove(tmp_file)

    return results_pd

def SurfaceArea(structure, chan_radius = 1.655, probe_radius = 1.655, num_samples = 5000, high_accuracy = True, prefix="tmp_sa"):

    """Analysis surface area of structure.

    Args:
        structure (str): path to your CIF.
        chan_radius (float): probe of channel, it is advised to keep chan_radius=probe_radius.
        probe_radius (float): probe of radiu.
        num_samples (int): number of MC samples per atom.
        high_accuracy (bool): use high accuracy or not.
        prefix (str): temporary file.

    Returns:
        Dictionary:
            -   unit by ["unit"], always Å^2, m^2/cm^3, m^2/g
            -   accessible surface area by ["ASA"]
            -   non-accessible surface area by ["NASA"]
    """

    results_sa = {}
    results_sa["unit"]="Å^2, m^2/cm^3, m^2/g"

    tmp_file = f"{prefix}.txt"
    
    if high_accuracy:
        cmd = f'network -ha -sa {chan_radius} {probe_radius} {num_samples} {tmp_file} {structure}'
        # cmd = f'network -ha S50 -sa {chan_radius} {probe_radius} {num_samples} {tmp_file} {structure}'
    else:
        cmd = f'network -sa {chan_radius} {probe_radius} {num_samples} {tmp_file} {structure}'
    _ = subprocess.run(
                        cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=True,
                    )
    with open(tmp_file) as f:
        for i, row in enumerate(f):
            if i == 0:
                ASA = float(row.split('ASA_A^2:')[1].split()[0])
                VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0])
                GSA = float(row.split('ASA_m^2/g:')[1].split()[0])
                NASA = float(row.split('NASA_A^2:')[1].split()[0])
                NVSA = float(row.split('NASA_m^2/cm^3:')[1].split()[0])
                NGSA = float(row.split('NASA_m^2/g:')[1].split()[0])

    results_sa["ASA"] = [ASA, VSA, GSA]
    results_sa["NASA"] = [NASA, NVSA, NGSA]

    os.remove(tmp_file)

    return results_sa

def PoreVolume(structure, chan_radius = 0, probe_radius = 0, num_samples = 5000, high_accuracy = True, prefix="tmp_pv"):

    """Analysis pore volume of structure.

    Args:
        structure (str): path to your CIF.
        chan_radius (float): probe of channel, it is advised to keep chan_radius=probe_radius.
        probe_radius (float): probe of radiu.
        num_samples (int): number of MC samples per atom.
        high_accuracy (bool): use high accuracy or not.
        prefix (str): temporary file.

    Returns:
        Dictionary:
            -   unit by ["unit"], always PV: Å^3, cm^3/g; VF: nan
            -   accessible pore volume by ["PV"]
            -   non-accessible pore volume by ["NPV"]
            -   accessible void fraction by ["VF"]
            -   non-accessible void fraction by ["NVF"]
    """

    results_pv = {}
    results_pv["unit"]="PV: Å^3, cm^3/g; VF: nan"

    tmp_file = f"{prefix}.txt"
    
    if high_accuracy:
        cmd = f'network -ha -volpo {chan_radius} {probe_radius} {num_samples} {tmp_file} {structure}'
        # cmd = f'network -ha S50 -volpo {chan_radius} {probe_radius} {num_samples} {tmp_file} {structure}'
    else:
        cmd = f'network -volpo {chan_radius} {probe_radius} {num_samples} {tmp_file} {structure}'
    _ = subprocess.run(
                        cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        check=True,
                    )
    with open(tmp_file) as f:
        for i, row in enumerate(f):
            if i == 0:
                POAV = float(row.split('POAV_A^3:')[1].split()[0])
                PONAV = float(row.split('PONAV_A^3:')[1].split()[0])
                GPOAV = float(row.split('POAV_cm^3/g:')[1].split()[0])
                GPONAV = float(row.split('PONAV_cm^3/g:')[1].split()[0])
                POAV_volume_fraction = float(row.split('POAV_Volume_fraction:')[1].split()[0])
                PONAV_volume_fraction = float(row.split('PONAV_Volume_fraction:')[1].split()[0])
    results_pv["PV"] = [POAV, GPOAV]
    results_pv["NPV"] = [PONAV, GPONAV]
    results_pv["VF"] = POAV_volume_fraction
    results_pv["NVF"] = PONAV_volume_fraction

    os.remove(tmp_file)

    return results_pv
