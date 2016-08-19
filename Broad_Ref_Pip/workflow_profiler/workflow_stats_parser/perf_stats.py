# TODO: polish imports!

import matplotlib.pyplot as plt
# import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
# import scipy.signal
import numpy

import argparse
import datetime
import os
import fnmatch
import subprocess
import re


def read_pcm(root_dir, steps):
	print('read pcm data...')
	out = dict()
	desc = dict()
	desc['time'] = ('Time',)
	desc['Read'] = ('GBs',)
	desc['Write'] = ('GBs',)

	cols_printed = False
	for step in steps:
		file_path = root_dir + '/' + step + '.pcm.txt'
		with open(file_path, 'r') as file:
			values = {}
			header_0 = file.readline()
			header_1 = file.readline()

			if not cols_printed:
				cols = header_1.split(';')
				cols_n = {}
				for i in range(0, len(cols)):
					cols_n[i] = cols[i]
				print(cols_n)
				cols_printed = True

			values['time'] = []
			values['Read'] = []
			values['Write'] = []
			for line in file:
				vals = line.split(';')
				step_time = datetime.datetime.strptime(vals[0] + ' ' + vals[1], "%Y-%m-%d %H:%M:%S.%f")
				values['time'].append(step_time)
				values['Read'].append(float(vals[13]))
				values['Write'].append(float(vals[14]))
			out[step] = values
	return desc, out

def read_collectl(root_dir, steps):
	print('read collectl data...')
	out = dict()
	desc = dict()
	multi = dict()
	desc['time'] = ('Time', '')
	desc['CpuTotl'] = ('%', '[CPU]Totl%')
	desc['CpuSys'] = ('%', '[CPU]Sys%')
	desc['CpuWait'] = ('%', '[CPU]Wait%')  # is it really msec?
	desc['IoRead'] = ('MBs/sec', '[DSK]ReadKBTot')
	desc['IoWrite'] = ('MBs/sec', '[DSK]WriteKBTot')
	desc['Mem'] = ('GBs', '[MEM]Commit')

	# TODO: rewrite; special case for disks
	multi['DSK'] = None

	hdr_processed = False
	hdr_index = dict()
	for step in steps:
		file_path = root_dir + '/' + step + '.collectl'
		with open(file_path, 'r') as file:
			values = {}
			hdr_line = file.readline()

			if not hdr_processed:
				hdr_names = hdr_line.split(' ')
				print('collectl header:')
				for i in range(0, len(hdr_names)):
					hdr_index[hdr_names[i]] = i
					print("{name}: {i} ".format(name=hdr_names[i], i=i), end=', ')
				print('\n')

				dsk_re = re.compile("\[DSK:(?P<dsk>\w+)\]Name$")
				disks = []
				for i in range(len(hdr_names)):
					dsk = dsk_re.match(hdr_names[i])
					if dsk:
						disk = dsk.group('dsk')
						disks.append(disk)
						desc['IOWAIT:' + disk] = ('Time(msec)', '[DSK:' + disk + ']Wait', '[DSK:' + disk + ']SvcTim')
				multi['DSK'] = disks
				print('disks:', disks)

				hdr_processed = True

			values['time'] = []
			values['CpuTotl'] = []
			values['CpuSys'] = []
			values['CpuWait'] = []
			values['IoRead'] = []
			values['IoWrite'] = []
			values['Mem'] = []
			for disk in multi['DSK']:
				values['IOWAIT:' + disk] = []

			for line in file:
				vals = line.split(' ')
				if vals[0][0] == '#':
					continue
				step_time = datetime.datetime.strptime(vals[0] + ' ' + vals[1], "%Y%m%d %H:%M:%S")
				values['time'].append(step_time)
				values['CpuTotl'].append(float(vals[10]))
				values['CpuSys'].append(float(vals[hdr_index[desc['CpuSys'][1]]]))
				values['CpuWait'].append(float(vals[5]))
				values['IoRead'].append(float(vals[63]) / 1024)
				values['IoWrite'].append(float(vals[64]) / 1024)
				values['Mem'].append(float(vals[32]) / (1024 * 1024))
				for disk in multi['DSK']:
					disk_name = 'IOWAIT:' + disk
					values[disk_name].append(float(vals[hdr_index[desc[disk_name][1]]]) +
										float(vals[hdr_index[desc[disk_name][2]]]))
			out[step] = values
	return desc, multi, out


def movingaverage(values, window_size):
	return numpy.convolve(values, numpy.ones(window_size) / window_size, 'same')


def plot_value(run_name, title, value_name, steps, desc, data, colors, linewidth=2):
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

	print('printing:', run_name, title)

	fig = plt.figure()
	fig.suptitle(run_name, fontsize=20)
	plt.title(title, fontsize=18)

	plt.subplot(111)

	plt.ylabel(desc[value_name][0], fontsize=15)
	plt.xlabel('Time(hours)', fontsize=15)

	t2h = Time2Hours()
	for step, color in zip(steps, colors):
		time = t2h(data[step]['time'])
		values = data[step][value_name]
		# values = scipy.signal.medfilt(values, 13)
		# values = movingaverage(values,12)
		plt.plot(time, values, color, label=step, linewidth=linewidth)

	plt.subplots_adjust(bottom=0.35)
	plt.legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.15), fontsize=10)


def plot_time(run_name, steps, colors, root_dir):
	print('printing: ', run_name, 'Steps Time')

	times = []
	times_hours = []
	for step in steps:
		with open(root_dir + '/' + step + '.cmd.txt', 'r') as cmd_file:
			lines = cmd_file.readlines()
			vals = lines[1].rstrip('\n').split(' ')
			start_time = datetime.datetime.strptime(vals[2] + ' ' + vals[3], "%Y-%m-%d %H:%M:%S")
			vals = lines[3].rstrip('\n').split(' ')
			end_time = datetime.datetime.strptime(vals[2] + ' ' + vals[3], "%Y-%m-%d %H:%M:%S")
			times.append(end_time - start_time)
			times_hours.append((end_time - start_time).total_seconds() / (3600))
		# print(times)

	fig = plt.figure()
	fig.suptitle(run_name, fontsize=20)
	plt.title('Steps Time', fontsize=18)

	plt.subplot(111)

	plt.ylabel('Time(hours)', fontsize=15)
	plt.xlabel('Step', fontsize=15)

	plt.xticks([])

	for t in list((range(len(steps)))):
		plt.bar(t, times_hours[t], color=colors[t], align='center', label=steps[t])
		plt.annotate(times[t], (t, times_hours[t]), va='bottom', ha='center', fontsize=10)

	plt.legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.15), fontsize=10)
	plt.subplots_adjust(bottom=0.35)


def prepare_files(root_dir):
	files = os.listdir(root_dir)

	# check if we have all *pcm.txt files
	for step in steps:
		if step + '.pcm.txt' not in files:
			raise RuntimeError('file is absent:', step + '.pcm.txt')

	# check collectl file and no duplicates there
	collectl_files = {}
	for step in steps:
		match_str = step + '.collectl-*.raw.gz'
		# print(step)
		step_files = []
		for file in files:
			# print("\t",file)
			if fnmatch.fnmatch(file, match_str):
				# print('fn---->',file)
				step_files.append(file)
		if len(step_files) > 1:
			raise RuntimeError('duplicates found:', step, step_files)
		if len(step_files) < 1:
			raise RuntimeError('no files for:', step)
		collectl_files[step_files[0]] = step

	# print(collectl_files)
	for file_in, file_out in collectl_files.items():
		in_path = root_dir + '/' + file_in
		out_path = root_dir + '/' + file_out + '.collectl'
		print(in_path, '-->', out_path)
		if not os.path.exists(out_path):
			with open(out_path, 'w') as out:
				subprocess.check_call(['collectl', '-P', '-p', in_path], stdout=out)

# colors=('red','orange','brown','green','cyan','blue','violet','yellow','green','cyan','blue','violet')
colors = ('red', 'gold', 'brown', 'greenyellow', 'cyan', 'blue', 'violet', 'purple', 'lime', 'magenta', 'steelblue',
			'blueviolet')


# bwa, sts2b, sts, md, sti, rtc, ir, br, pr, hc
cls_old = ('Bwa.0',
			'SamtoolsSamToBam.0',
			'SamtoolsSort.0',
			'PicardMarkDuplicates.0',
			'SamtoolsIndex.0',
			'GatkRealignerTargetCreator.0',
			'GatkIndelRealigner.0',
			'GatkBaseRecalibrator.0',
			'GatkPrintReads.0',
			'GatkHaplotypeCaller.0')

sge = ('Bwa',
		'SamtoolsSamToBam',
		'SamtoolsSort',
		'PicardMarkDuplicates',
		'GatkRealignerTargetCreator',
		'GatkIndelRealigner',
		'GatkBaseRecalibrator',
		'GatkPrintReads',
		'GatkHaplotypeCaller')

cls = ('Bwa',
		'SamtoolsSam2Bam',
		'SamtoolsSort',
		'PicardMarkDuplicates',
		'SamtoolsIndex',
		'GatkRealignerTargetCreator',
		'GatkIndelRealigner',
		'GatkBaseRecalibrator',
		'GatkPrintReads',
		'GatkHaplotypeCaller')

sge_old = ('GatkQueue.Bwa',
			'GatkQueue.SamtoolsSamToBam',
			'GatkQueue.SamtoolsSort',
			'GatkQueue.PicardMarkDuplicates',
			'GatkQueue.GatkRealignerTargetCreator',
			'GatkQueue.GatkIndelRealigner',
			'GatkQueue.GatkBaseRecalibrator',
			'GatkQueue.GatkPrintReads',
			'GatkQueue.GatkHaplotypeCaller')

sge_old_old = ('GatkQueue.Bwa',
			'GatkQueue.SamtoolsSamToBam',
			'GatkQueue.SamtoolsSamSort',
			'GatkQueue.PicardMarkDuplicates',
			'GatkQueue.GatkRealignerTargetCreator',
			'GatkQueue.GatkIndelRealigner',
			'GatkQueue.GatkBaseRecalibrator',
			'GatkQueue.GatkPrintReads',
			'GatkQueue.GatkHaplotypeCaller')

parser = argparse.ArgumentParser(description='perf stats generator')
parser.add_argument('--title', metavar='TITLE', type=str, help='sample name with notes', required=True)
parser.add_argument('--root_dir', metavar='ROOT DIR', type=str, help='root dir for data', required=True)
parser.add_argument('--output', metavar='OUTPUT', type=str, help='pdf output name', required=True)
args = parser.parse_args()

root_dir = args.root_dir  # '/home/sergeo/proj/bench/28_08_2015/sge_NA12878_2_t1_i5'
pdf_out = args.output  # '/media/sf_usr/tmp/28_08_2015-sge_NA12878_2_t1_i5.pdf'
title = args.title  # 'sge_NA12878_2_t1_i5'
linewidth = 1
dpi = 100
steps = sge_old_old

print('root_dir:', root_dir)
print('pdf_out:', pdf_out)
print('title:', title)

print('data files preparation...')
prepare_files(root_dir)

(pcm_desc, pcm) = read_pcm(root_dir, steps)
(collectl_desc, multi, collectl) = read_collectl(root_dir, steps)

plt.rc('figure', figsize=(11.69, 8.27))
pdf = PdfPages(pdf_out)

# '''
plot_time(title, steps, colors, root_dir)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'CPU Total Utilization', 'CpuTotl', steps, collectl_desc, collectl, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'CPU System Utilization', 'CpuSys', steps, collectl_desc, collectl, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'CPU I/O Wait', 'CpuWait', steps, collectl_desc, collectl, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'DRAM Read Usage', 'Read', steps, pcm_desc, pcm, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'DRAM Write Usage', 'Write', steps, pcm_desc, pcm, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'Total Commited Memory', 'Mem', steps, collectl_desc, collectl, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'Total HDD I/O Read Bandwidth', 'IoRead', steps, collectl_desc, collectl, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)

plot_value(title, 'Total HDD I/O Write Bandwidth', 'IoWrite', steps, collectl_desc, collectl, colors, linewidth=linewidth)
plt.savefig(pdf, format='pdf', dpi=dpi)
# '''

for i in range(len(multi['DSK'])):
	disk_name = multi['DSK'][i]
	value_name = 'IOWAIT:' + disk_name
	plot_value(title, 'HDD IOWAIT ' + disk_name, value_name, steps, collectl_desc, collectl, colors, linewidth=linewidth)
	plt.savefig(pdf, format='pdf', dpi=dpi)
# '''
# plt.show()
pdf.close()
exit()
