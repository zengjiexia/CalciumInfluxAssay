import local_tools 
import os
import pandas as pd
import math as ms
import numpy as np

mainPath = r"C:\Users\zx252\Documents\20191205\20191205_Dimitri_Soaked_brain"


print('Starting...')

if not os.path.isdir(mainPath):
	print('Data folder does not exist. Exit.')
	quit()

resultPath = mainPath + '/results'

### get all samples in the folder ###
sampleNames = [name for name in os.listdir(mainPath) if not name.startswith('.') or name == 'results']
sampleSummary = []

### Loop over all samples ###
for sample in sampleNames:

	print('Running sample: ' + sample)
	sampleErr = 0

	### Set the paths ###
	ionomycinPath = mainPath + '/' + sample + '/Ionomycin/'
	samplePath = mainPath + '/' + sample + '/Sample/'
	blankPath = mainPath + '/' + sample + '/Blank/'

	sampleOutput = pd.DataFrame()

	### Obtain filenames for fields of view ###
	fieldNames = local_tools.extract_filename(ionomycinPath)
	fieldSummary = []

	### Loop over all fields of views ###
	for c, field in enumerate(fieldNames, 1):

		### Average tiff files ###
		ionomycinMean = local_tools.average_frame(ionomycinPath + field)
		sampleMean = local_tools.average_frame(samplePath + field)
		blankMean = local_tools.average_frame(blankPath + field)

		### Align blank and sample images to the ionomycin image ###
		sampleAligned, blankAligned = local_tools.img_alignment(ionomycinMean, sampleMean, blankMean)

		### Locate the peaks on the ionomycin image ###
		peaks = local_tools.peak_locating(ionomycinMean, 200)
		
		### Calculate the intensities of peaks with certain radius (in pixel) ###
		ionInten = local_tools.intensities(ionomycinMean, peaks, 3)
		samInten = local_tools.intensities(sampleAligned, peaks, 3)
		blaInten = local_tools.intensities(blankAligned, peaks, 3)
		
		### Calculate influx of each single liposome and count errors ###
		influx = pd.DataFrame((samInten - blaInten)/(ionInten - blaInten)*100, columns=['Influx'])

		""" 
		if 100% < influx < 200% take as 100%
		if -100% < influx < 0% take as 0
		if influx calculated to be nan or <-100 or >200 count as error
 		""" 
		influx['Influx'] = [100 if i >= 100 and i <= 200 else i for i in influx['Influx']]
		influx['Influx'] = [0 if i <= 0 and i >= -100 else i for i in influx['Influx']]
		influx['Influx'] = ['error' if ms.isnan(np.float(i)) or i < -100 or i > 200 else i for i in influx['Influx']]

		try:
			fieldErr = (influx.Influx.values == 'error').sum()
		except (AttributeError, FutureWarning) as e:
			fieldErr = 0

		### Propagate field Err into sample Err ###
		sampleErr += fieldErr

		### Filter out error data ###
		influx = influx[influx.Influx != 'error']

		### Generate a dataframe which contains the result of current field of view ###
		fieldOutput = pd.concat([
			pd.DataFrame(np.tile(c, (len(peaks), 1)), columns=['Field']),
			pd.DataFrame(peaks, columns=['X', 'Y']),
			influx				
			],axis = 1)

		### Record the mean influx of this field of view ###
		fieldSummary.append([c, fieldOutput.loc[:, 'Influx'].mean(), str(round(fieldErr/len(peaks)*100, 2)) + '%']) 

		### Merge the result of current field into the sample dataframe ###
		sampleOutput = pd.concat([sampleOutput, fieldOutput])

	### Reset the index for sample dataframe ###
	sampleOutput = sampleOutput.reset_index(drop = True)

	### Record the mean influx of this sample ###
	sampleSummary.append([sample, sampleOutput.loc[:, 'Influx'].mean(), str(round(sampleErr/len(sampleOutput)*100, 2)) + '%'])

	print(sample + ' with mean influx: ' + str(sampleOutput.loc[:, 'Influx'].mean()))
	print('Percentage of error in this sample: ' + str(err/len(sampleOutput)*100) + '%')

	### Save the result for current sample ###
	sampleOutput.to_csv(resultPath + '/raw/' + sample + '.csv', index=False)
	fieldSummary_df = pd.DataFrame(fieldSummary, columns=['Field', 'Influx', r'% Error'])
	fieldSummary_df.to_csv(resultPath + '/raw/' + sample + '_field.csv')

### Save sample summaries ###
sampleSummary_df = pd.DataFrame(sampleSummary, columns=['Sample', 'Influx', r"% Error"])
sampleSummary_df.to_csv(resultPath + '/summary.csv')
