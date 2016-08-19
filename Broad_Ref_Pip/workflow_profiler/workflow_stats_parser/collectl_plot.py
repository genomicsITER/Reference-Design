# TODO: polish imports!

import matplotlib as mpl                                                                                                                                                                                           
mpl.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
# import scipy.signal
import numpy

import logging
import argparse
import datetime
import sys
import os
import fnmatch
import subprocess
import re
import collections

import workflow_dictionaries
#from workflow_dictionaries import *


class CollectlColumn:
	def __init__(self, title, axis_text):
		self.title = title
		self.axis_text = axis_text
		self.data = None

	def append(self, values):
		pass


class CpuTotal(CollectlColumn):
	def __init__(self, hdr_index):
		super(CpuTotal, self).__init__('CPU Total Utilization', '%')
		self.data = list()
		self.hdr_index = hdr_index

	def append(self, values):
		value = float(values[self.hdr_index['[CPU]Totl%']])
		self.data.append(value)


class CpuUser(CollectlColumn):
	def __init__(self, hdr_index):
		super(CpuUser, self).__init__('CPU User Utilization', '%')
		self.data = list()
		self.hdr_index = hdr_index

	def append(self, values):
		value = float(values[self.hdr_index['[CPU]User%']])
		self.data.append(value)


class CpuSystem(CollectlColumn):
	def __init__(self, hdr_index):
		super(CpuSystem, self).__init__('CPU System Utilization', '%')
		self.data = list()
		self.hdr_index = hdr_index

	def append(self, values):
		value = float(values[self.hdr_index['[CPU]Sys%']])
		self.data.append(value)


class CpuIOWait(CollectlColumn):
	def __init__(self, hdr_index):
		super(CpuIOWait, self).__init__('CPU IOWait', '%')
		self.data = list()
		self.hdr_index = hdr_index

	def append(self, values):
		value = float(values[self.hdr_index['[CPU]Wait%']])
		self.data.append(value)

class CommitedMem(CollectlColumn):
	def __init__(self, hdr_index):
		super(CommitedMem, self).__init__('Total Commited Memory', 'GBs/s')
		self.data = list()
		self.hdr_index = hdr_index

	def append(self, values):
		value = float(values[self.hdr_index['[MEM]Commit']]) / (1024 * 1024)
		self.data.append(value)


class HddIORead(CollectlColumn):
	def __init__(self, hdr_index):
		super(HddIORead, self).__init__('Total HDD IORead Bandwidth', 'MBs/sec')
		self.data = list()
		self.hdr_index = hdr_index

	def append(self, values):
		value = float(values[self.hdr_index['[DSK]ReadKBTot']])/1024
		self.data.append(value)


class HddIOWrite(CollectlColumn):
	def __init__(self, hdr_index):
		super(HddIOWrite, self).__init__('Total HDD IOWrite Bandwidth', 'MBs/sec')
		self.data = list()
		self.hdr_index = hdr_index

	def append(self, values):
		value = float(values[self.hdr_index['[DSK]WriteKBTot']])/1024
		self.data.append(value)


class HddIOWait(CollectlColumn):
	def __init__(self, hdr_index):
		super(HddIOWait, self).__init__('HDD IOWAIT', 'msec')
		self.hdr_index = hdr_index
		self.data = collections.OrderedDict()
		
		disk_re = re.compile("\[DSK:(?P<disk>\w+)\]Name$")
		for key in self.hdr_index.keys():
			disk = disk_re.match(key)
			if disk:
				disk = disk.group('disk')
				self.data[disk] = list()

	def append(self, values):
		for disk in self.data.keys():
			value = float(values[self.hdr_index['[DSK:' + disk + ']Wait']]) + float(values[self.hdr_index['[DSK:' + disk + ']SvcTim']])
			#TODO: probably not the fastest solution
			self.data[disk].append(value)


class LustreIORead(CollectlColumn):
	def __init__(self, hdr_index):
		super(LustreIORead, self).__init__('LustreFS IORead Bandwidth', 'MBs/sec')
		self.hdr_index = hdr_index
		self.data = collections.OrderedDict()
		
		sub_col_re = re.compile("\[CLT:(?P<sub_col>\w+)\]ReadKB$")
		for col in self.hdr_index.keys():
			sub_col = sub_col_re.match(col)
			if sub_col:
				sub_col = sub_col.group('sub_col')
				self.data[sub_col] = list()

	def append(self, values):
		for key in self.data.keys():
			value = float(values[self.hdr_index['[CLT:' + key + ']ReadKB']]) / 1024
			#TODO: probably not the fastest solution
			self.data[key].append(value)


class LustreIOWrite(CollectlColumn):
	def __init__(self, hdr_index):
		super(LustreIOWrite, self).__init__('LustreFS IOWrite Bandwidth', 'MBs/sec')
		self.hdr_index = hdr_index
		self.data = collections.OrderedDict()
		
		sub_col_re = re.compile("\[CLT:(?P<sub_col>\w+)\]WriteKB$")
		for col in self.hdr_index.keys():
			sub_col = sub_col_re.match(col)
			if sub_col:
				sub_col = sub_col.group('sub_col')
				self.data[sub_col] = list()

	def append(self, values):
		for key in self.data.keys():
			value = float(values[self.hdr_index['[CLT:' + key + ']WriteKB']]) / 1024
			#TODO: probably not the fastest solution
			self.data[key].append(value)

class LustreTotalIORead(CollectlColumn):
	def __init__(self, hdr_index):
		super(LustreTotalIORead, self).__init__('LustreFS Total IORead Bandwidth', 'MBs/sec')
		self.hdr_index = hdr_index
		self.data = list()

	def append(self, values):
		value = float(values[self.hdr_index['[CLT]ReadKB']])/1024
		self.data.append(value)


class LustreTotalIOWrite(CollectlColumn):
	def __init__(self, hdr_index):
		super(LustreTotalIOWrite, self).__init__('LustreFS Total IOWrite Bandwidth', 'MBs/sec')
		self.hdr_index = hdr_index
		self.data = list()

	def append(self, values):
		value = float(values[self.hdr_index['[CLT]WriteKB']])/1024
		self.data.append(value)


class Step:
	def __init__(self, name, path, columns, time_format='%Y%m%d %H:%M:%S'):
		self.name = name
		self.path = path
		self.time = list()
		self.time_text = 'Time(hours)'
		self.columns = collections.OrderedDict()

		logging.info('step load: ' + self.name + ' ' + str([column.__name__ for column in columns]))

		hdr = collections.OrderedDict()
		with open(self.path) as file:
			hdr_line = file.readline().strip()
			hdr_columns = hdr_line.split(' ')
			for i in range(len(hdr_columns)):
				hdr[hdr_columns[i]]=i
			logging.debug('header: ' + self.name + ' ' + str(hdr))

			for column in columns:
				column_object = column(hdr)
				self.columns[column.__name__] = column_object

			for line in file:
				if line[0] == '#': continue
				values = line.split()
				time_value = datetime.datetime.strptime(values[0] + ' ' + values[1], time_format)
				self.time.append(time_value)
				for column in self.columns.values():
					column.append(values)
		loaded_columns = list()
		for key in self.columns.keys():
			if type(self.columns[key].data) is list:
				loaded_columns.append(str(key))
			else:
				loaded_columns.append(str(key) + str(list(self.columns[key].data.keys())))
		
		logging.info('step loaded:' + self.name +': ' + str(loaded_columns))

	def __repr__(self):
		return "Step{{name={name}, path={path}}}".format(name=self.name, path=self.path)

	def __str__(self):
		return "Step{{name={name}, path={path}}}".format(name=self.name, path=self.path)


def plot_column(run_name, title, x_text, y_text, steps_names, x_data, y_data, colors, linewidth):
	class Time2Hours:
		def __init__(self):
			self.offset = None

		def __call__(self, time):
			out = []
			for t in time:
				if not self.offset:
					self.offset = t
				out.append((t - self.offset).total_seconds() / 3600)
			return out

	logging.info('plot:' + run_name + ': ' + title)

	fig = plt.figure()
	fig.suptitle(title, fontsize=20)
	plt.title(run_name, fontsize=18)

	plt.subplot(111)

	plt.ylabel(y_text, fontsize=15)
	plt.xlabel(x_text, fontsize=15)

	t2h = Time2Hours()
	for step_name, x, y, color in zip(steps_names, x_data, y_data, colors):
		time = t2h(x)
		values = y
		#TODO: do we need smoothing?
		# values = scipy.signal.medfilt(values, 13)
		# values = movingaverage(values,12)
		plt.plot(time, values, color, label=step_name, linewidth=linewidth)

	plt.subplots_adjust(bottom=0.35)
	leg = plt.legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.15), fontsize=10, framealpha=0.5, title='Worflow phase(s)')
	for obj in leg.legendHandles:
		obj.set_linewidth(3.0)

def movingaverage(values, window_size):
	return numpy.convolve(values, numpy.ones(window_size) / window_size, 'same')

def plot_steps(sample, steps, linewidth, dpi, output):
	if not len(steps): return
	if not len(steps[0].columns): return
	#TODO: check consistency if all steps have the same columns

	pdf = None
	if not os.path.isdir(output):
		logging.debug('pdf: ' + output)
		pdf = PdfPages(output)

	colors=('red','orange','brown','green','cyan','blue','violet','yellow','green','cyan','blue','violet')
	plt.rc('figure', figsize=(11.69, 8.27))
	
	for column_name in steps[0].columns.keys():
		x_text = steps[0].time_text
		y_text = steps[0].columns[column_name].axis_text

		if type(steps[0].columns[column_name].data) is list:
			title = steps[0].columns[column_name].title
			steps_names = list()
			x_data = list()
			y_data = list()
			for step in steps:
				steps_names.append(step.name)
				x_data.append(step.time)
				y_data.append(step.columns[column_name].data)
			plot_column(sample, title, x_text, y_text, steps_names, x_data, y_data, colors, linewidth)
			if pdf: plt.savefig(pdf, format='pdf', dpi=dpi)
			else:
				file_name = sample + '_' + title
				#TODO: be sure filename is conformant
				file_name = file_name.replace(' ', '_')
				plt.savefig(os.path.join(output, file_name + '.png'), format='png', dpi=dpi, transparent=True)
		else:
			for data_name in steps[0].columns[column_name].data.keys():
				title = steps[0].columns[column_name].title + ' ' + data_name
				steps_names = list()
				x_data = list()
				y_data = list()
				for step in steps:
					steps_names.append(step.name)
					x_data.append(step.time)
					y_data.append(step.columns[column_name].data[data_name])
				plot_column(sample, title, x_text, y_text, steps_names, x_data, y_data, colors, linewidth)
				if pdf: plt.savefig(pdf, format='pdf', dpi=dpi)
				else:
					file_name = sample + '_' + title
					#TODO: be sure filename is conformant
					file_name = file_name.replace(' ', '_')
					plt.savefig(os.path.join(output, file_name + '.png'), format='png', dpi=dpi, transparent=True)
	if pdf:
		pdf.close()


def load_steps(root_dir, steps_files, columns):
	steps = list()
	for step_name, step_file in steps_files.items():
		steps.append(Step(step_name, step_file, columns))
	return steps


def prepare_files(root_dir, steps, sub_dirs=False, lustre=False, force=False):
	root_dir = os.path.normpath(root_dir)
	steps_files = collections.OrderedDict()
	
	def os_walk_onerror(error):
		raise error
	
	all_files = list()
	for dirpath, dirnames, files in os.walk(root_dir, topdown=False, onerror=os_walk_onerror):
		for name in files:
			all_files.append(os.path.join(dirpath,name))
			
	logging.debug('files_in_root:' + str(all_files))

	for step_name, step_srch in steps.items():
		step_files = list()
		if sub_dirs:
			match_str = os.path.join(root_dir,'*' + step_srch + '*','*' + step_srch + '*.collectl-*.raw.gz')
		else:
			match_str = os.path.join(root_dir,'*' + step_srch + '*.collectl-*.raw.gz')
		logging.debug('match_str:' + step_name +': ' + match_str)
		for file in all_files:
			if fnmatch.fnmatch(file, match_str):
				step_files.append(file)
		if len(step_files) > 1:
			raise RuntimeError('duplicates found:', step_name, step_files)
		if len(step_files) < 1:
			raise RuntimeError('no files for:', step_name)
		steps_files[step_name] = step_files[0]
	logging.debug('steps_files: ' + str(steps_files))

	for step_name, path in steps_files.items():
		dir, file_name = os.path.split(path)
		out_path = os.path.join(dir, file_name[:-len('raw.gz')] + 'csv')
		logging.info('create_csv:' +step_name +': ' + path + ' --> ' + out_path)
		if not os.path.exists(out_path) or force:
			with open(out_path, 'w') as out:
				#--import lustreClient.ph,s
				if lustre:
					cmd = ['collectl', '--import', 'lustreClient.ph,s', '-P', '-p', path]
				else:
					cmd = ['collectl', '-P', '-p', path]
				subprocess.check_call(cmd, stdout=out)
		steps_files[step_name]=out_path

	return steps_files


def main():
	LOG_LEVEL_MAP = {'debug': logging.DEBUG, 'info': logging.INFO}

	parser = argparse.ArgumentParser(description='perf stats generator')
	parser.add_argument('-s', '--sample', metavar='SAMPLE', type=str, help='sample name', required=True)
	parser.add_argument('-N', '--workflow_name', metavar='WORKFLOW_NAME', type=str, help='workflow name (workflow dictionaries)', required=True)
	parser.add_argument('-r', '--root_dir', metavar='ROOT DIR', type=str, help='root dir for data', required=True)
	parser.add_argument('-o', '--output', metavar='OUTPUT', type=str, help='pdf output name', required=True)
	parser.add_argument('--lustre', help='collect lustre fs', action='store_true', default=False)
	parser.add_argument('-f', '--force', help='force to generate new csvs', action='store_true', default=False)
	parser.add_argument('-l', '--log', help='Specify the logging level', choices=LOG_LEVEL_MAP.keys(), default='info')
	args = parser.parse_args()

	logging.basicConfig(format='%(levelname)s:[%(asctime)s]:%(message)s', level=LOG_LEVEL_MAP[args.log])

	linewidth = 1
	dpi = 100

	logging.debug('workflow_parse_dict:' + str(workflow_dictionaries.workflow_parse_dict))
	steps_name = workflow_dictionaries.workflow_parse_dict[args.workflow_name]
	steps = getattr(workflow_dictionaries, steps_name)

	if args.lustre:
		columns = (CpuTotal, CpuUser, CpuSystem, CpuIOWait, CommitedMem, HddIORead, HddIOWrite, HddIOWait, LustreTotalIORead, LustreTotalIOWrite)
	else:
		columns = (CpuTotal, CpuUser, CpuSystem, CpuIOWait, CommitedMem, HddIORead, HddIOWrite, HddIOWait)

	steps_files=prepare_files(args.root_dir, steps, sub_dirs=True, lustre=args.lustre, force=args.force)
	steps = load_steps(args.root_dir, steps_files, columns)
	plot_steps(args.sample, steps, linewidth, dpi, args.output)

	#plt.show()
	exit()

if __name__ == "__main__":
	sys.exit(main())
