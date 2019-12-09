from local_tools import *
import os

Path = r'/Users/zengjiexia/Desktop/2019-12-06 liposome_frac_2_Derya'


result_path = r'/Users/zengjiexia/Desktop/results'

holder = {
'Threshold' : '50',
'Radius' : '3',
'High' : 200,
'Low' : -100
}

def one_sample(root, item):
	sample_path = os.path.join(root, item)
	print('starting '+ sample_path)
	path = {
		'main': sample_path,
		'ionomycin': sample_path + '/Ionomycin/',
		'sample': sample_path + '/Sample/',
		'blank': sample_path + '/Blank/'
	}
	# Error = 0, Safe = 1
	# error_report{'path': [main, ionomycin, sample, blank]}
	error_report = {'path': [os.path.isdir(p) for p in path.values()]} 


	ionomycin = []
	sample = []
	blank = []

	filenames = extract_filename(path['ionomycin'])
	for name in filenames:

		ionomycin_mean = average_frame(path['ionomycin'] + name)
		sample_mean = average_frame(path['sample'] + name)
		blank_mean = average_frame(path['blank'] + name)

		sample_aligned, blank_aligned = img_alignment(ionomycin_mean, sample_mean, blank_mean)

		ionomycin.append(ionomycin_mean)
		sample.append(sample_aligned)
		blank.append(blank_aligned)

	Results = pd.DataFrame()

	for n in range(0, len(ionomycin)):
		print('working on file number '+str(n))

		coordinations = peak_locating(ionomycin[n], 50)
		peak = pd.DataFrame({
			'field': np.repeat(n, len(coordinations)),
			'x': coordinations[:, 0],
			'y': coordinations[:, 1]
			})

		result, error = influx_calculation(
			ionomycin[n],
			sample[n],
			blank[n],
			peak,
			holder['High'],
			holder['Low'],
			int(holder['Radius']))
		error_report['calculation field ' + str(n)] = error

		Results = pd.concat([Results, result])

	Results.to_excel(result_path+'/'+ item+'.xlsx')
	print('Done with '+ sample_path)

	return Results.influx.mean()



samplenames = os.listdir(Path)
for item in samplenames:
	influx = one_sample(Path, item)
	with open(r"/Users/zengjiexia/Desktop/results/summary.txt",'a') as f:
		f.write(item + '\t' + str(influx)+'\n')


