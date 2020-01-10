# -*- coding: utf-8 -*-
"""
LAST UPDATE 22/11/2019

Supporting code for analyze fluorescence images taken by TIRFM
Primarily made for Protein Aggregate Detection Assay project

Must work with UI.py in the same folder.

Created by Zengjie Xia in November 2018.
Copyright © 2018-2019 Zengjie Xia. All rights reserved.
Version: 2.0
"""

import cv2
import os
import tifffile as tiff
from scipy.ndimage import filters
from scipy import ndimage
import pandas as pd
import numpy as np
import math as ms
import time


def extract_filename(path):
	"""
	walk through a directory and put names of all tiff files into an ordered list
	para: path - string
	return: filenames - list of string 
	"""

	filenames = []
	for root, dirs, files in os.walk(path):
		for name in files:
			if name.endswith('.tif'):
				filenames.append(name)
	filenames = sorted(filenames)
	return filenames


def average_frame(path):
	"""
	input 'path' for stacked tiff file and the 'number of images' contained
	separate individual images from a tiff stack.
	para: path - string
	return: ave_img - 2D array
	"""

	ori_img = tiff.imread(path)
	ave_img = np.mean(ori_img, axis=0)
	ave_img = ave_img.astype('uint16')

	return ave_img


def img_alignment(Ionomycin, Sample, Blank):
	"""
	image alignment based on cross-correlation
	Ionomycin image is the reference image
	para: Ionomycin, Sample, Blank - 2D array
	return: Corrected_Sample, Corrected_Blank - 2D array
	"""

	centre_ = (Ionomycin.shape[0]/2, Ionomycin.shape[1]/2)
	# 2d fourier transform of averaged images
	FIonomycin = np.fft.fft2(Ionomycin)
	FSample = np.fft.fft2(Sample)
	FBlank = np.fft.fft2(Blank)

	# Correlation based on Ionomycin image
	FRIS = FIonomycin*np.conj(FSample)
	FRIB = FIonomycin*np.conj(FBlank)

	RIS = np.fft.ifft2(FRIS)
	RIS = np.fft.fftshift(RIS)
	RIB = np.fft.ifft2(FRIB)
	RIB = np.fft.fftshift(RIB)

	[i, j] = np.where(RIS == RIS.max())
	[g, k] = np.where(RIB == RIB.max())

	# offset values
	IS_x_offset = i-centre_[0]
	IS_y_offset = j-centre_[1]
	IB_x_offset = g-centre_[0]
	IB_y_offset = k-centre_[1]

	# Correction
	MIS = np.float64([[1, 0, IS_x_offset], [0, 1, IS_y_offset]])
	Corrected_Sample = cv2.warpAffine(Sample, MIS, Ionomycin.shape)
	MIB = np.float64([[1, 0, IB_x_offset], [0, 1, IB_y_offset]])
	Corrected_Blank = cv2.warpAffine(Blank, MIB, Ionomycin.shape)

	return Corrected_Sample, Corrected_Blank


def peak_locating(data, threshold):
	"""
	Credit to Daniel
	para: data - 2D array
	para: threshold - integer
	return: xy_thresh - 2D array [[x1, y1], [x2, y2]...]
	"""

	data_max = filters.maximum_filter(data, 3)
	maxima = (data == data_max)
	data_min = filters.minimum_filter(data, 3)
	diff = ((data_max - data_min) > threshold)
	maxima[diff == 0] = 0

	labeled, num_objects = ndimage.label(maxima)
	xy = np.array(ndimage.center_of_mass(data, labeled, range(1, num_objects+1)))
	xy_thresh = np.zeros((0, 2))
	for row in xy:
		a = row[0]
		b = row[1]
		if (a > 30) and (a < 480) and (b > 30) and (b < 480):
			ab = np.array([np.uint16(a), np.uint16(b)])
			xy_thresh = np.vstack((xy_thresh, ab))
	xy_thresh = xy_thresh[1:] 

	return xy_thresh


def intensities(image_array, peak_coor, radius):
	"""
	When the local peak is found, extract all the coordinates of pixels in a 'radius'
	para: image_array - 2D array
	para: peak_coor - 2D array [[x1, y1], [x2, y2]]
	para: radius - integer
	return: intensities - 2D array [[I1], [I2]]
	"""

	x_ind, y_ind = np.indices(image_array.shape)
	intensities = np.zeros((0,1))
	#intensities = []
	for (x, y) in peak_coor:
		intensity = 0
		circle_points = ((x_ind - x)**2 + (y_ind - y)**2) <= radius**2
		coor = np.where(circle_points == True)
		coor = np.array(list(zip(coor[0], coor[1])))
		for j in coor:
			intensity += image_array[j[0], j[1]]
		intensities = np.vstack((intensities, intensity))
		#intensities.append(intensity)

	return intensities
	#return np.array(intensities)

def influx_calculation(Ionomycin, Sample, Blank, peak_coor, high, low, radius):
	"""
	Equation is:
		(F_Sample-F_Blank)*100/(F_Ionomycin-F_Blank)
	In the event of photoblench took place, some results are normalized/ignored
				> high limited, < low limited     error
				100%-high limited				  100%
				0-100%							  itself
				low limited-0					  0
	para: Ionomycin, Sample, Blank - 2D array
	para: peak_coor - pd.DataFrame {x:, y:}
	para: radius, high, low - integer
	return: results - pd.DataFrame {x:, y:, ionomycin:, sample:, blank:, influx:}
	return: error - integer
	"""

	error = 0
	results = pd.DataFrame(peak_coor)

	results.columns=['field', 'x', 'y', 'ionomycin', 'sample', 'blank']
	results['influx'] = (results['sample'] - results['blank']) / (results['ionomycin'] - results['blank']) * 100

	results['influx'] = [100 if i >= 100 and i <= high else i for i in results['influx']]
	results['influx'] = [0 if i <= 0 and i >= low else i for i in results['influx']]
	results['influx'] = ['error' if ms.isnan(np.float(i)) or i < low or i > high else i for i in results['influx']]

	try:
		error = (results.influx.values == 'error').sum()
	except (AttributeError, FutureWarning) as e:
		error = 0

	results = results[results.influx != 'error']

	return results, error

if __name__ == '__main__':
	
	Holder = {
	'PATH' : os.path.dirname(os.path.abspath(__file__))+'/sample_data',
	'Threshold' : '80',
	'Radius' : '3',
	'High' : 200,
	'Low' : -100
	}

	print('Running code in test mode')
	if not os.path.isdir(Holder['PATH']):
		print('Sample data not found. Exit test mode.')
		quit()

	resultPath = os.path.dirname(os.path.abspath(__file__))+'/sample_results'

	test_round = 0

	time_list = []
	while test_round <= 0:

		tic = time.clock()

		### get all samples in the folder ###
		sampleNames = [name for name in os.listdir(Holder['PATH']) if not name.startswith('.')]
		sampleSummary = {}

		### Loop over all samples ###
		for sample in sampleNames:
			err = 0

			### Set the paths ###
			ionomycinPath = Holder['PATH'] + '/' + sample + '/Ionomycin/'
			samplePath = Holder['PATH'] + '/' + sample + '/Sample/'
			blankPath = Holder['PATH'] + '/' + sample + '/Blank/'

			sampleOutput = pd.DataFrame()

			### Obtain filenames for fields of view ###
			fieldNames = extract_filename(ionomycinPath)
			fieldSummary = {}

			### Loop over all fields of views ###
			for c, field in enumerate(fieldNames, 1):

				### Average tiff files ###
				ionomycinMean = average_frame(ionomycinPath + field)
				sampleMean = average_frame(samplePath + field)
				blankMean = average_frame(blankPath + field)

				### Align blank and sample images to the ionomycin image ###
				sampleAligned, blankAligned = img_alignment(ionomycinMean, sampleMean, blankMean)

				### Locate the peaks on the ionomycin image ###
				peaks = peak_locating(ionomycinMean, 80)
				
				### Calculate the intensities of peaks with certain radius (in pixel) ###
				ionInten = intensities(ionomycinMean, peaks, 3)
				samInten = intensities(sampleMean, peaks, 3)
				blaInten = intensities(blankMean, peaks, 3)
				
				influx = pd.DataFrame((samInten - blaInten)/(ionInten - blaInten)*100, columns=['Influx'])
				
				### for test ###
				#influx = influx.append({'Influx': float('nan')}, ignore_index=True)

				influx['Influx'] = [100 if i >= 100 and i <= 200 else i for i in influx['Influx']]
				influx['Influx'] = [0 if i <= 0 and i >= -100 else i for i in influx['Influx']]
				influx['Influx'] = ['error' if ms.isnan(np.float(i)) or i < -100 or i > 200 else i for i in influx['Influx']]

				try:
					err += (influx.Influx.values == 'error').sum()
				except (AttributeError, FutureWarning) as e:
					err += 0

				influx = influx[influx.Influx != 'error']


				### Generate a dataframe which contains the result of current field of view ###
				fieldOutput = pd.concat([
					pd.DataFrame(np.tile(c, (len(peaks), 1)), columns=['Field']),
					pd.DataFrame(peaks, columns=['X', 'Y']),
					influx				
					],axis = 1)

				print(fieldOutput)
				quit()
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

		toc = time.clock()
		time_list.append(toc-tic)
		print(str(test_round))
		test_round += 1
	print(np.array(time_list).mean())

	print('Exit test mode')