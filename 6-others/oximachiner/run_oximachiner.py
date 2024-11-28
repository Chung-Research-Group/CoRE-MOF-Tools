from oximachinerunner import OximachineRunner
import json
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import os
import glob
from tqdm import tqdm
import concurrent

def process_mof(mof):
    runner = OximachineRunner() 
    name = mof.split("/")[-1].replace(".cif", "")
    try:
        result = runner.run_oximachine(mof)
        return (name, {"metal_symbols": result["metal_symbols"], "prediction": result["prediction"]})
    except Exception as e:
        print(f"{name} failed due to {str(e)}")
        return (name, "fail")

def main():
    folder = "./dataset/"
    mofs = glob.glob(os.path.join(folder, '*.cif'))
    results = {}

    with ProcessPoolExecutor() as executor:
        future_to_mof = {executor.submit(process_mof, mof): mof for mof in mofs}
        for future in tqdm(concurrent.futures.as_completed(future_to_mof), total=len(mofs)):
            name, result = future.result()
            results[name] = result

    with open("oximachiner.json", "w") as f:
        json.dump(results, f, indent=4, default=lambda o: int(o) if isinstance(o, np.int64) else o)

if __name__ == "__main__":
    main()
