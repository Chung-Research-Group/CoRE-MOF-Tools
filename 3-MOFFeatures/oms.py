from mof_collection import MofCollection
import warnings

def run(input_folder,output_folder,num_batches=1):
    a_mof_collection = MofCollection.from_folder(collection_folder=input_folder)
    a_mof_collection.analyse_mofs(num_batches=num_batches,overwrite=False)
    print(a_mof_collection.mof_oms_df)
    a_mof_collection.mof_oms_df.to_csv(output_folder + "/oms.csv")
