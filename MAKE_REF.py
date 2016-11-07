from AMP_COV import COV_AMP
from sys import argv
import glob
from multiprocessing import Process, Queue
import multiprocessing as mp
import numpy


files = glob.glob(argv[1]+"*.bam")##input folder


queue_arr = [Queue()]*len(files)

def bed_pro(file_name):
	f = open(file_name,"r")
	result_arr = []
	for x in f.readlines():
		#if x.find("Pool#1-Amp#18")==-1 and x.find("Pool#2-Amp#1")==-1:
		temp = x.split("\t")
		if len(temp)<=4:
			result_arr.append([temp[0],temp[1],temp[2],temp[3].strip()])
		elif len(temp) > 4:
			result_arr.append([temp[0],temp[1],temp[2],temp[3].strip(),temp[4].strip()])
	f.close()

	return result_arr

def cov_run(file_name, bed_name, method, output):
	cov_c = COV_AMP(file_name,bed_name) ##input, BED, output
	dat_arr = cov_c.write_proc(method)##method
	output.put(dat_arr)

def barcode_sorting(arr):
	temp_res = []
	for x in range(len(arr)):
		for y in range(len(arr[x])):
			temp_res.append(int(arr[x][y][0]))

	return sorted(list(set(temp_res)))


def normalized_cov(arr):
	res = []
	arr_val = []
	for x in arr:
		val = int(x[4])
		arr_val.append(val)

	arr_min = min(arr_val)
	arr_max = max(arr_val)
	m_diff = arr_max - arr_min

	for x in arr:
		val = int(x[4])
		ch1 = float(val - arr_min)

		norm_v = float(ch1)/float(m_diff)
		res.append([x[0],x[1],x[2],x[3],norm_v])

	return res




def bet_samp_cal(arr):
	res = []
	for x in arr:
		avr = []
		for y in arr:
			temp_a = float(x[4])
			temp_b = float(y[4])

			if temp_a==0.0:
				temp_a = 1.0
			if temp_b==0.0:
				temp_b = 1.0

			af = temp_a/temp_b
			avr.append(af)

		avr_val = numpy.mean(avr)
		res.append([x[0],x[1],x[2],x[3],avr_val])
	return res


def bet_samp_cal_cov(arr):
	res = []
	t_cov = []
	print arr
	for x in arr:
		t_cov.append(int(x[4]))

	t_cov = sum(t_cov)
	print t_cov
	avr = []

	for x in arr:
		#print x
		temp_a = float(x[4])
		#print temp_a

		if temp_a==0.0:
			temp_a = 1.0

		af = float(temp_a)/float(t_cov)
		res.append([x[0],x[1],x[2],x[3],af])

	return res



def external_cal(arr1,ref):

	res = []
	for x in arr1:
		for y in ref:
			if x[1]==y[1] and x[2]==y[2] and x[3]==y[3]:
				temp_a = float(x[4])
				temp_b = float(y[4])

				if temp_a==0.0:
					temp_a = 1.0
				if temp_b==0.0:
					temp_b = 1.0

				af = temp_a/temp_b
				res.append([x[0],x[1],x[2],x[3],af])
	return res


def z_score_cal(arr1, tarr):
	res = []
	for x in arr1:
		t_arr = []
		for sa in tarr:
			for y in sa:
				if x[1]==y[1] and x[2]==y[2] and x[3]==y[3] and x[0]!=y[0]:

					t_arr.append(float(y[4]))

		v_mean = numpy.mean(t_arr)
		v_std = numpy.std(t_arr)

		vf = float(x[4])-v_mean
		z_score =  float(vf)/float(v_std)

		res.append([x[0],x[1],x[2],x[3],x[4],z_score])

	return res


processes = [mp.Process(target=cov_run, args=(files[n],argv[2],argv[3],queue_arr[n])) for n in range(len(files))] #####input, BED, output, method

for p in processes:
	p.start()

res_arr = [None]*len(files)

for a in range(len(files)):
	res_arr[a] = queue_arr[a].get()

for a in range(len(files)):
	queue_arr[a].close()

for p in processes:
	p.join()


rw = open(argv[1]+"ref_analysis.txt","w")#####output
barcode_temp = barcode_sorting(res_arr)

ref_dig = int(argv[4])
ref_arr_avr = []
res_arr_avr = []


for sr in res_arr:
	print "#################"
	#print  sr
		#res_arr_avr.append(bet_samp_cal(sr))
	res_arr_avr.append(bet_samp_cal_cov(sr))


"""
#####p-value calculation
p_val_inc = []
for sr in fin_arr_avr:
	p_val_inc.append(z_score_cal(sr,fin_arr_avr))

fin_arr_avr = p_val_inc
#####p-value calculation
"""

"""
bar_str = "\t"*4
for b in barcode_temp:
	bar_str = bar_str+"Barcode"+str(b)+"\t"

bar_str = bar_str[:len(bar_str)-1]+"\n"
rw.write(bar_str)
"""


for loc in bed_pro(argv[2]):
	d_str = str(loc[0])+"\t"+str(loc[1])+"\t"+str(loc[2])+"\t"
	tavr = []
	for b in barcode_temp:
		for sr in range(len(res_arr_avr)):
			for dat in res_arr_avr[sr]:
				if loc[0]==dat[1] and loc[1]==dat[2] and loc[2]==dat[3] and b==dat[0]:
					d_str = d_str+str(dat[4])+":"
					#d_str = d_str+str(dat[4])+"\t"
					tavr.append(float(dat[4]))


	d_str = d_str.strip()+"\t"+str(numpy.mean(tavr))
	#d_str = d_str.strip()+"\t"+str(numpy.mean(tavr))+"\t"+str(numpy.std(tavr))
	d_str = d_str.strip()+"\n"
	rw.write(d_str)
