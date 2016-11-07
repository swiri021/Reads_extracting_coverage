from AMP_COV import COV_AMP
from sys import argv
import glob
from multiprocessing import Process, Queue
import multiprocessing as mp
import numpy


files = glob.glob(argv[1]+"edit/"+"*.bam")##input folder


queue_arr = [Queue()]*len(files)

def bed_pro(file_name):
	f = open(file_name,"r")
	result_arr = []
	for x in f.readlines():
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
		#print arr
		for y in range(len(arr[x])):
			temp_res.append(int(arr[x][y][0]))

	return sorted(list(set(temp_res)))


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


rw = open(argv[1]+"/edit/coverage/amplicon_speci_cov.txt","w")#####output

barcode_temp = barcode_sorting(res_arr)

if argv[3] != "A":

	bar_str = "\t"*4
	for b in barcode_temp:
		bar_str = bar_str+"Barcode"+str(b)+"\t"

	bar_str = bar_str[:len(bar_str)-1]+"\n"
	rw.write(bar_str)


	for loc in bed_pro(argv[2]):
		d_str = str(loc[0])+"\t"+str(loc[1])+"\t"+str(loc[2])+"\t"+str(loc[3])+"\t"
		for b in barcode_temp:
			for sr in range(len(res_arr)):
				for dat in res_arr[sr]:
					if loc[0]==dat[1] and loc[1]==dat[2] and loc[2]==dat[3] and b==dat[0]:
						d_str = d_str+str(dat[4])+"\t"

		d_str = d_str.strip()+"\n"
		rw.write(d_str)


else:
	bar_str = "\t"*4
	for b in barcode_temp:
		bar_str = bar_str+"Barcode"+str(b)+"\t"+""+"\t"

	bar_str = bar_str[:len(bar_str)-1]+"\n"
	rw.write(bar_str)


	for loc in bed_pro(argv[2]):
		d_str = str(loc[0])+"\t"+str(loc[1])+"\t"+str(loc[2])+"\t"+str(loc[3])+"\t"
		for b in barcode_temp:
			for sr in range(len(res_arr)):
				for dat in res_arr[sr]:
					if loc[0]==dat[1] and loc[1]==dat[2] and loc[2]==dat[3] and b==dat[0]:
						d_str = d_str+str(dat[4])+"\t"+str(dat[5])+"\t"

		d_str = d_str.strip()+"\n"
		rw.write(d_str)

rw.close()



