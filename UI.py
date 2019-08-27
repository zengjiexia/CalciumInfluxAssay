# -*- coding: utf-8 -*-
"""
LAST UPDATE 26/08/2019

User Interface for analyze fluorescence images taken by TIRFM
Primarily made for Calcium Influx Assay (Liposome Assay) project

This is a User Interface code, which must work with local_tools.py in the same folder.

Created by Zengjie Xia in July 2019.
Copyright © 2019 Zengjie Xia. All rights reserved.
Version: 2.0
"""

import local_tools
from tkinter import *
from tkinter import ttk
import tkinter.messagebox as tkmbox
from tkinter import filedialog
import tkinter.scrolledtext as scrolledtext
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import compress
import datetime
import os
import locale
import time
import shutil
locale.setlocale(locale.LC_ALL, '')


Holder = {
	'PATH' : os.path.dirname(os.path.abspath(__file__)),
	'Threshold' : '80',
	'Radius' : '3',
	'High' : 200,
	'Low' : -100
}

# Holder['PATH'] = r'C:\Users\zx252\Documents\Projects\Calcium_Influx_Image_Analysis\20190315_pd_samples'

class UserInterface(Frame):
	global Holder

	def __init__(self, master):
		tik = time.time()
		
		Frame.__init__(self, master)
		# Master interface
		self.master = master
		master.geometry('660x240')
		master.resizable(True, True) 
		master.title('Calcium Influx Analysis')

		# Control Frame
		self.control_frame = Frame(master)
		self.control_frame.grid(row=0, column=0, sticky=E+W)

			# Directory 'PATH'
		Label(self.control_frame, text='Directory:', anchor=E).grid(row=0, column=0, sticky=E)
		self.dir_var = StringVar()
		self.dir_var.set(Holder['PATH'])
		self.dir_entry = Entry(self.control_frame, textvariable=self.dir_var, width='80')
		self.dir_entry.bind('<Return>', self.updateDirectory_return)
		self.dir_entry.bind('<FocusOut>', self.updateDirectory_return)
		self.dir_entry.grid(row=0, column=1, sticky=E+W)
		self.dir_button = Button(self.control_frame, text='Browse', command=self.updateDirectory_browse, width='15', anchor=CENTER)
		self.dir_button.grid(row=0, column=2, sticky=E+W)

			# Threshold
		Label(self.control_frame, text='Threshold:', anchor=E).grid(row=1, column=0, sticky=E)
		self.thre_var = StringVar()
		self.thre_var.set(Holder['Threshold'])
		self.thre_entry = Entry(self.control_frame, textvariable=self.thre_var)
		self.thre_entry.bind('<Return>', self.updateThreshold)
		self.thre_entry.bind('<FocusOut>', self.updateThreshold)
		self.thre_entry.grid(row=1, column=1, sticky=E+W)

			# Radius
		Label(self.control_frame, text='Radius:', anchor=E).grid(row=2, column=0, sticky=E)
		self.radius_var = StringVar()
		self.radius_var.set(Holder['Radius'])
		self.radius_entry = Entry(self.control_frame, textvariable=self.radius_var)
		self.radius_entry.bind('<Return>', self.updateRadius)
		self.radius_entry.bind('<FocusOut>', self.updateRadius)
		self.radius_entry.grid(row=2, column=1, sticky=E+W)

			# Start
		Button(self.control_frame, text='Start', anchor=CENTER, command=self.start).grid(row=1, rowspan=2, column=2, sticky=E+W+S+N)	

		ttk.Separator(self.control_frame, orient='horizontal').grid(row=3, column=0, columnspan=3, sticky=E+W)

		# Status Frame
		self.status_frame = scrolledtext.ScrolledText(master, height=10, width=80)
		self.status_frame.grid(row=1, column=0, sticky=E+W)
		self.status_frame.insert(END, str(datetime.datetime.now().strftime("%x,%X"))+' Welcome.\n')
		self.status_frame.config(state="disabled")
		self.updateStatus('If you encountered any problem, please try to restart. If can not be solved, contact me by email: zx252@cam.ac.uk')

		toc = time.time()
		print('UI loading time: ' + str(toc-tik) + ' s')

	def updateStatus(self, text):
		"""
		Update message in the status frame
		"""
		self.status_frame.config(state="normal")
		text = str(datetime.datetime.now().strftime("%x,%X")) + '\t' + text + '\n'
		self.status_frame.insert(END, text)
		self.status_frame.config(state="disabled")
		self.status_frame.see(END)
		self.status_frame.update()


	def updateDirectory_return(self, event):
		
		new_directory = self.dir_var.get()
		
		if os.path.isdir(new_directory):
			Holder['PATH'] = new_directory
			self.dir_var.set(Holder['PATH'])
			self.updateStatus('Directory updated.')
			self.master.update()

		else:
			self.updateStatus('Directory does not exist.')
			self.dir_var.set(Holder['PATH'])
			self.master.update()
	
	def updateDirectory_browse(self):
		
		new_directory = filedialog.askdirectory(
			initialdir = Holder['PATH'],
			title = 'Main Directory for Analysis'
			)
		
		if os.path.isdir(new_directory):
			Holder['PATH'] = new_directory
			self.dir_var.set(Holder['PATH'])
			self.updateStatus('Directory updated.')
			self.master.update()

		else:
			self.updateStatus('No directory selected.')

	def updateThreshold(self, event):
		new_thre = self.thre_var.get()
		new_thre = new_thre.split('/')

		for i in new_thre: # remove any non pure numberic term
			try:
				int_i = int(i)
			except ValueError:
				new_thre.remove(i)

		new_thre = list(dict.fromkeys(new_thre)) # remove duplicates

		if len(new_thre) > 5:
			self.updateStatus('WARNING: More than 5 thresholds input. Only first 5 will be used.')
			new_thre = new_thre[:5] # remove extra terms
		

		[self.updateStatus('WARNING: '+i+' might be too small or too large for thresholding') for i in new_thre if int(i) > 200 or int(i) < 20] # warning for small/large threshold

		new_thre = '/'.join(new_thre)
		if new_thre == '':
			Holder['Threshold'] = '80'
		else:
			Holder['Threshold'] = new_thre

		self.thre_var.set(Holder['Threshold'])
		self.updateStatus('Threshold updated.')
		self.master.update()

	def updateRadius(self, event):
		new_radius = self.radius_var.get()
		new_radius = ''.join(c for c in new_radius if c.isdigit())
		Holder['Radius'] = new_radius
		self.radius_var.set(Holder['Radius'])
		self.updateStatus('Radius updated.')
		if int(Holder['Radius']) >= 10:
			self.updateStatus('WARNING: The radius might be too large.')
		self.master.update()

	def start(self):
		result_path = Holder['PATH'] + '/Results'
		if os.path.isdir(result_path):
			answer = tkmbox.askquestion(title='Pre-exist Result Folder', message='Result folder existed. Do you want to replace the old result?')
			if answer == 'yes':
				shutil.rmtree(result_path)
			else:
				self.updateStatus('Programme terminated.')
				return 0

		os.makedirs(result_path)

		thresholds = [int(thre) for thre in Holder['Threshold'].split('/')]
		sample_names = next(os.walk(Holder['PATH']))[1]
		sample_names = [name for name in sample_names if name != 'Results']
		self.summary = pd.DataFrame()

		items = ['Main', 'Ionomycin', 'Sample', 'Blank']
		est = 0
		for name in sample_names:
			if est == 0:
					tic = time.time()

			sample = local_tools.CalciumSample(Holder['PATH'] + '/' + name, Holder) 
			try:
				error = sample.error_report['path'].index(0)
			except ValueError:
				self.updateStatus('Starting '+name)

				sample.img_correction()
				for thre in thresholds:
					self.updateStatus('At threshold:' + str(thre))
					data_file = pd.DataFrame()
					sample.peak_location(thre)
					result = sample.influx()
					for i in result.values():
						data_file = pd.concat([data_file, i]) 
					data_file.to_csv(result_path+'/'+name+'_at_'+str(thre)+'.csv')
					'''
					plt.figure()
					hist_plot = sns.distplot(data_file['influx'])
					hist_plot.figure.savefig(result_path+'/'+name+'_at_'+str(thre)+'.png')
					plt.close()
					'''
					if est == 0:
						toc = time.time()
						est_time = str(datetime.timedelta(seconds=(np.ceil((toc - tic + 15)*len(sample_names)*len(thresholds)))))
						self.updateStatus('Estimated processing time: '+est_time)
						est = 1

					summa = pd.DataFrame({
						'file': [name],
						'threshold': [thre],
						'influx': [round(data_file['influx'].mean(),2)],
						'n': [len(data_file.index)]
						})
					self.summary = pd.concat([self.summary, summa])

					self.updateStatus('Influx for '+name +' at ' + str(thre)+'(threshold):' + str(round(data_file['influx'].mean(),2))+'%')

			else:
				error_items = ', '.join(list(compress(items, np.subtract(1, sample.error_report['path']))))
				self.updateStatus('Path error with ' + error_items + 'folder(s) in ' + name +'.\n'+name+' was skipped.')

		self.summary.to_csv(result_path+'/summary.csv')
		self.updateStatus('Done. Total time cost: '+ str(datetime.timedelta(seconds = (np.ceil(date.time()-tic)))))

if __name__ == "__main__":

	root = Tk()
	UserInterface(root).grid(row=0, column=0, sticky=E+W+S+N)
	root.mainloop()
