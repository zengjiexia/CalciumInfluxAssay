

from local_tools import *
import os

holder = {
'PATH' : r"D:\Temporary_20191205\20191205_Derya\1",
'Threshold' : '80',
'Radius' : '3',
'High' : 200,
'Low' : -100
}

sample_path = holder['PATH']
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



for n in range(0, len(ionomycin)):
	coordinations = peak_locating(ionomycin[n], 80)
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

	print(result)



