import threading
import time
import argparse
import numpy as np
import pandas as pd

# Additional imports
from ops.imports import *
from ops.process import Align
import ops.firesnake
from ops.firesnake import Snake
import ops.in_situ
import tifffile
import nd2


from cellpose import models
import torch

class DataProcessor:
    def __init__(self, total_cycle, num_threads):
        self.total_cycle = total_cycle
        self.num_threads = num_threads
        #self.QThreshold = 0.1
        self.THRESHOLD_DAPI = 250
        self.THRESHOLD_CELL = 50
        self.THRESHOLD_STD = 50
        self.QualityThreshold = 0.000
        self.directory_path = '/gpfs/gsfs10/users/sagarm2/OPSdata/230915Data6wellplate/' # seventh data 10x six well
        self.subfolders = [f.path for f in os.scandir(self.directory_path) if f.is_dir()]
        self.subfolders.sort(key=None, reverse=False)



    def inputGenerationfromSingleNd2(self,totalCycle, indexInCycle):

		############# input generation
        dataRPE = np.zeros((totalCycle, 5, 2304, 2304))


		
		## this is number of cycles
        for i in range(totalCycle):
            searchRPE=self.subfolders[i]+"/*.nd2"
            input_files = natsorted(glob(searchRPE))
            print(input_files[indexInCycle])
            nd2file=nd2.imread(input_files[indexInCycle])
            #print(nd2file.shape)
            dataRPE[i,] = nd2file
		
		 
        return dataRPE
    
    def cellposeSegNucG(self,img, min_size=15):
        # Check if GPU is available and set use_gpu accordingly
        use_gpu = torch.cuda.is_available()

        model = models.Cellpose(gpu=use_gpu, model_type="nuclei")
        
        # Include the use_gpu parameter in the eval function
        res, _, _, _ = model.eval(
            img,
            channels=[0, 0],
            diameter=None,
            min_size=min_size,
            #gpu=use_gpu  # Add this line
        )
        return res

    def getBarcodeFromDataCPSeg(self,dataT, THRESHOLD_DAPI = 800, THRESHOLD_CELL = 250, THRESHOLD_STD = 80, qualityIndex=0):
    
		# gets barcode from one FOV- series
        aligned1 = Snake._align_SBS(dataT, method='DAPI')
        loged1 = Snake._transform_log(aligned1, skip_index=0)
        maxed1 = Snake._max_filter(loged1, 3, remove_index=0)
        std1 = Snake._compute_std(loged1, remove_index=0)
        peaks1 = Snake._find_peaks(std1)

		#THRESHOLD_DAPI = 800 # got from ImageJ DAPI image
        NUCLEUS_AREA = 0.25*400, 0.25*35000 #analyze particles imagej
        WILDCARDS = dict(well='A1', tile='7')
		#NUCLEUS_AREA = 0.25*100, 0.25*35000 #analyze particles imagej
		
		#nuclei1 = Snake._segment_nuclei(dataT[0], THRESHOLD_DAPI, area_min=NUCLEUS_AREA[0], area_max=NUCLEUS_AREA[1])

        nuclei_cp1 = self.cellposeSegNucG(dataT[0])
		#THRESHOLD_CELL = 250 #300
		# got excellent result with 250
        cells1 = Snake._segment_cells(dataT[0], nuclei_cp1, THRESHOLD_CELL)

		
        df_bases1 = Snake._extract_bases(maxed1, peaks1, cells1, 
			            THRESHOLD_STD, wildcards=WILDCARDS, bases = 'GTAC')
	   # print(df_bases)
        df_reads1 = Snake._call_reads(df_bases1) # does the quality Q index calculation
		
	   
        df_cells1 = Snake._call_cells2(df_reads1, qualityIndex)
		
		#rt=df_cells['cell_barcode_0'].value_counts() # distribution of barcodes
    
        return df_bases1, df_reads1, df_cells1




    def process_data(self, locationInFolder, positionInSeries):
        # Simulate data processing; replace these lines with actual function calls
        dataRPE_np = self.inputGenerationfromSingleNd2(12, positionInSeries) # new function style reading for comparison

        
        data_read = dataRPE_np.astype(np.uint16)
        df_bases_Train, df_reads_Train, df_cell_Train=self.getBarcodeFromDataCPSeg(data_read, 100, self.THRESHOLD_CELL, self.THRESHOLD_STD, self.QualityThreshold)
        
        # Placeholder for actual processing functions
        # Assuming getBarcodeFromDataCPSeg returns dummy pandas dataframes

        csvname_reads = f"Output20230915/df_reads_Well{locationInFolder}_Position{positionInSeries:02}.csv"
        csvname_cells = f"Output20230915/df_cells_Well{locationInFolder}_Position{positionInSeries:02}.csv"
        
        df_reads_Train.to_csv(csvname_reads)
        df_cell_Train.to_csv(csvname_cells)



    def worker(self, locationInFolder, positionInSeries):
        self.process_data(locationInFolder, positionInSeries)

    def run(self, start_location, end_position):
        threads = []
        for locationInFolder in range(start_location):
            for positionInSeries in range(end_position):
                if len(threads) >= self.num_threads:
                    for thread in threads:
                        thread.join()
                    threads = []
                t = threading.Thread(target=self.worker, args=(locationInFolder, positionInSeries))
                threads.append(t)
                t.start()

        for thread in threads:
            thread.join()

        print("All threads have finished.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process data using multi-threading.")
    parser.add_argument("--start_location", type=int, default=0, help="Start index for location in folder")
    parser.add_argument("--end_position", type=int, default=750, help="End index for position in series")
    parser.add_argument("--total_cycle", type=int, default=12, help="Total cycle count")
    parser.add_argument("--num_threads", type=int, default=6, help="Number of concurrent threads")

    args = parser.parse_args()

    start_time = time.time()
    print(args.total_cycle, args.num_threads)
    processor = DataProcessor(total_cycle=args.total_cycle, num_threads=args.num_threads)

    processor.run(start_location=args.start_location, end_position=args.end_position)
    elapsed_time = time.time() - start_time
    print('Time taken:', elapsed_time)

