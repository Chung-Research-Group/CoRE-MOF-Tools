# import sys
# import os
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# from CoREMOF.curate import clean_pacman
# import glob


# cifs = glob.glob("./ion_pacman/*cif")
# for cif in cifs:
#     print(cif)
#     clean_pacman(structure=cif, initial_skin=0.25, output_folder="ion_pacman/result", saveto="clean_result.csv")

   
# from CoREMOF.calculation.Zeopp import FrameworkDim
# example_cif = "./cal/IRMOF-1.cif"
# results_strinfo = FrameworkDim(structure=example_cif,
#                                         high_accuracy=True,
#                                         prefix="test_strinfo") 
# print(results_strinfo["Dimention"])