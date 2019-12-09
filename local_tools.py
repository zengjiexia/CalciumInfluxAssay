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
	para: peak_coor - pd.DataFrame {field:, x:, y:}
	para: radius - integer
	return: intensities - pd.DataFrame  {0:}
	"""

	x_ind, y_ind = np.indices(image_array.shape)
	intensities = []

	for x, y in zip(peak_coor['x'], peak_coor['y']):
		intensity = 0
		circle_points = ((x_ind - x)**2 + (y_ind - y)**2) <= radius**2
		coor = np.where(circle_points == True)
		coor = np.array(list(zip(coor[0], coor[1])))
		for j in coor:
			intensity += image_array[j[0], j[1]]
		intensities.append(intensity)
	intensities = pd.DataFrame(intensities)

	return intensities


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
	results = pd.concat(
		[peak_coor, 
		intensities(Ionomycin, peak_coor, radius), 
		intensities(Sample, peak_coor, radius), 
		intensities(Blank, peak_coor, radius)], 
		axis=1
		)
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


class CalciumSample:

	def __init__(self, sample_path, Holder):
		
		self.Holder = Holder
		self.path = {
		'main': sample_path,
		'ionomycin': sample_path + '/Ionomycin/',
		'sample': sample_path + '/Sample/',
		'blank': sample_path + '/Blank/'
		}
		# Error = 0, Safe = 1
		# error_report{'path': [main, ionomycin, sample, blank]}
		self.error_report = {'path': [os.path.isdir(p) for p in self.path.values()]} 
		
	def img_correction(self):

		self.no_of_field = 0

		self.ionomycin = {}
		self.sample = {}
		self.blank = {}

		filenames = extract_filename(self.path['ionomycin'])
		for name in filenames:

			ionomycin_mean = average_frame(self.path['ionomycin'] + name)
			sample_mean = average_frame(self.path['sample'] + name)
			blank_mean = average_frame(self.path['blank'] + name)

			sample_aligned, blank_aligned = img_alignment(ionomycin_mean, sample_mean, blank_mean)

			self.ionomycin[self.no_of_field] = ionomycin_mean
			self.sample[self.no_of_field] = sample_aligned
			self.blank[self.no_of_field] = blank_aligned
			self.no_of_field += 1
		
		return 1

	def peak_location(self, threshold):
		# return a dataset containing peak located on each field based on inputed threshold
		'''
		if self.peaks in globals():	# clear cache
			del self.peaks
		else:
			pass
		'''
		self.peaks = {} # self.peaks={'field': DataFrame{x:,y:}}
		for n in range(0, self.no_of_field):
			coordinations = peak_locating(self.ionomycin[n], threshold)
			peak = pd.DataFrame({
				'field': np.repeat(n, len(coordinations)),
				'x': coordinations[:, 0],
				'y': coordinations[:, 1]
				})
			self.peaks[n] = peak

		return self.peaks

	def influx(self):
		'''
		if self.results in globals():
			del self.results
		else:
			pass
		'''
		self.results = {}
		for n in range(0, self.no_of_field):
			result, error = influx_calculation(
				self.ionomycin[n],
				self.sample[n],
				self.blank[n],
				self.peaks[n],
				self.Holder['High'],
				self.Holder['Low'],
				int(self.Holder['Radius']))
			self.error_report['calculation field ' + str(n)] = error
			self.results[n] = result

		return self.results
		
'''
if __name__ == '__main__':
	
	Holder = {
	'PATH' : os.path.dirname(os.path.abspath(__file__)),
	'Threshold' : '80',
	'Radius' : '3',
	'High' : 200,
	'Low' : -100
	}

	print('Input -h or -help for mannul')

	stay = 1
	command = None
	while stay:
		if command = '-h' or '-help':
			print("""-info Current setting information
				-path PATH change to a new path""")
'''