import local_tools 
import os


mainPath = os.path.dirname(os.path.abspath(__file__))+'/sample_data',


print('Starting...')

if not os.path.isdir(mainPath)
	print('Data folder does not exist. Exit.')
	quit()

resultPath = mainPath + '/results'

### get all samples in the folder ###
sampleNames = [name for name in os.listdir(mainPath) if not name.startswith('.')]
sampleSummary = {}

### Loop over all samples ###
for sample in sampleNames:

	### Set the paths ###
	ionomycinPath = mainPath + '/' + sample + '/Ionomycin/'
	samplePath = mainPath + '/' + sample + '/Sample/'
	blankPath = mainPath + '/' + sample + '/Blank/'

	sampleOutput = pd.DataFrame()

	### Obtain filenames for fields of view ###
	fieldNames = local_tools.extract_filename(ionomycinPath)
	fieldSummary = {}

	### Loop over all fields of views ###
	for c, field in enumerate(fieldNames, 1):

		### Average tiff files ###
		ionomycinMean = local_tools.average_frame(ionomycinPath + field)
		sampleMean = local_tools.average_frame(samplePath + field)
		blankMean = local_tools.average_frame(blankPath + field)

		### Align blank and sample images to the ionomycin image ###
		sampleAligned, blankAligned = local_tools.img_alignment(ionomycinMean, sampleMean, blankMean)

		### Locate the peaks on the ionomycin image ###
		peaks = local_tools.peak_locating(ionomycinMean, 80)
		
		### Calculate the intensities of peaks with certain radius (in pixel) ###
		ionInten = local_tools.intensities(ionomycinMean, peaks, 3)
		samInten = local_tools.intensities(sampleAligned, peaks, 3)
		blaInten = local_tools.intensities(blankAligned, peaks, 3)
		
		### Generate a dataframe which contains the result of current field of view ###
		fieldOutput = pd.concat([
			pd.DataFrame(np.tile(c, (len(peaks), 1))),
			pd.DataFrame(peaks),
			pd.DataFrame((samInten - blaInten)/(ionInten - blaInten)*100)				
			],axis = 1)

		fieldOutput.columns = ['Field', 'X', 'Y', 'Influx']

		### Record the mean influx of this field of view ###
		fieldSummary[c] = fieldOutput.loc[:, 'Influx'].mean()

		### Merge the result of current field into the sample dataframe ###
		sampleOutput = pd.concat([sampleOutput, fieldOutput])

	### Reset the index for sample dataframe ###
	sampleOutput = sampleOutput.reset_index(drop = True)

	### Record the mean influx of this sample ###
	sampleSummary[sample] = sampleOutput.loc[:, 'Influx'].mean()

	### Save the result for current sample ###
	sampleOutput.to_csv(resultPath + '/raw/' + sample + '.csv', index=False)
	fieldSummary_df = pd.DataFrame.from_dict(fieldSummary, orient='index', columns=['Influx'])
	fieldSummary_df.to_csv(resultPath + '/raw/' + sample + '_field.csv')

### Save sample summaries ###
sampleSummary_df = pd.DataFrame.from_dict(sampleSummary, orient='index', columns=['Influx'])
sampleSummary_df.to_csv(resultPath + '/summary.csv')