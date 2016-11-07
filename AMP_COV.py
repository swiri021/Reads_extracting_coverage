import pysam
import numpy

class COV_AMP:

	def bed_pro(self,file_name):
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

	def string_count(self,dat):
		count = 0
		dele = ''
		del_count = 0
		for a in dat:
			if del_count ==0:
				if a.isalpha()==True:
					count += 1
				elif a=='^':
					del_count = 1
					count += 1
			elif del_count ==1:
				if a.isalpha()==True:
					count = count
				elif a.isdigit()==True:
					del_count = 0

		return count

	def file_name_pro(self,f_name):
		ln = f_name.split("/")
		ln = ln[len(ln)-1]

		if ln.find("IonX") > -1:
			temp = ln.split("_")
			bn = temp[1]

			return int(bn)

	def amp_covering_rate(self,amp_st, amp_en, read_st, read_en, strand):
		amp_cov_rate = 0.0
		if strand=="F":
			total_amp = amp_en - amp_st
			if read_st <= amp_st:
				if read_en <= amp_en:
					cov = read_en - amp_st
					amp_cov_rate = float(cov)/float(total_amp)
				else:
					amp_cov_rate = 1.0
			else:
				if read_en <= amp_en:
					cov = read_en - read_st
					amp_cov_rate = float(cov)/float(total_amp)
				else:
					cov = amp_en - read_st
					amp_cov_rate = float(cov)/float(total_amp)
		else:
			total_amp = amp_en - amp_st
			if read_en >= amp_en:
				if read_st >= amp_st:
					cov = amp_en - read_st
					amp_cov_rate = float(cov)/float(total_amp)
				else:
					amp_cov_rate = 1.0

			else:
				if read_st <= amp_st:
					cov = read_en - amp_st
					amp_cov_rate = float(cov)/float(total_amp)
				else:
					cov = read_en - read_st
					amp_cov_rate = float(cov)/float(total_amp)

		return amp_cov_rate


	def distance_measure(self,read, bed_arr):

		result_arr = []
		#if read.is_unmapped==False:
		strand = read.is_reverse
		if strand==False:
			strand = 'F'
		else:
			strand = 'R'
		ref_st = read.reference_start
		ref_end = read.reference_end
		tag = read.get_tag('MD')

		for y in bed_arr:
			if strand == 'F':
				result_arr.append(ref_st-int(y[1]))

			else:
				result_arr.append(int(y[2])-ref_end)

		result_edit_arr = [None]*len(result_arr)

		for y in range(len(result_arr)):
			result_edit_arr[y] = abs(int(result_arr[y]))

		min_value = min(result_edit_arr)
		min_value_index = 0

		if result_edit_arr.count(min_value) == 2:####### count2
			min_value_index = result_edit_arr.index(int(min_value))
		else:
			min_value_index = result_edit_arr.index(int(min_value))


		if result_arr[min_value_index] > 0:
			if self.string_count(tag) > 10:
				result_edit_arr[min_value_index] = 100000
				min_value = min(result_edit_arr)
				min_value_index = result_edit_arr.index(min_value)

		result = bed_arr[min_value_index]

		return result

	def covered_rate_measure(self,read, bed_arr):

		result_arr = []
		#if read.is_unmapped==False:
		strand = read.is_reverse
		if strand==False:
			strand = 'F'
		else:
			strand = 'R'
		ref_st = read.reference_start
		ref_end = read.reference_end
		total_bp = ref_end - ref_st
		tag = read.get_tag('MD')

		for y in bed_arr:
			if strand == 'F':
				if ref_st <= int(y[1]):
					cov_bp = ref_end-int(y[1])
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)

				else:
					cov_bp = int(y[2])-ref_st
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)

			else:
				if ref_end >= int(y[2]):
					cov_bp = int(y[2]) - ref_st
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)

				else:
					cov_bp = ref_end-int(y[1])
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)


		max_value = max(result_arr)
		max_value_index = 0

		if result_arr.count(max_value) == 2:####### count2
			max_value_index = result_arr.index(max_value)
		else:
			max_value_index = result_arr.index(max_value)

		result = bed_arr[max_value_index]
		covering_per = self.amp_covering_rate(int(result[1]),int(result[2]),ref_st, ref_end,strand)

		return result, covering_per

	def amp_covered_rate_measure(self,read, bed_arr):

		result_arr = []
		#if read.is_unmapped==False:
		strand = read.is_reverse
		if strand==False:
			strand = 'F'
		else:
			strand = 'R'
		ref_st = read.reference_start
		ref_end = read.reference_end
		total_bp = ref_end - ref_st
		tag = read.get_tag('MD')

		for y in bed_arr:
			if strand == 'F':
				if ref_st <= int(y[1]):
					cov_bp = ref_end-int(y[1])
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)

				else:
					cov_bp = int(y[2])-ref_st
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)

			else:
				if ref_end >= int(y[2]):
					cov_bp = int(y[2]) - ref_st
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)

				else:
					cov_bp = ref_end-int(y[1])
					cov_rate = float(cov_bp)/float(total_bp)
					result_arr.append(cov_rate)


		max_value = max(result_arr)
		max_value_index = 0

		if result_arr.count(max_value) == 2:####### count2
			max_value_index = result_arr.index(max_value)
		else:
			max_value_index = result_arr.index(max_value)

		result = bed_arr[max_value_index]
		covering_per = self.amp_covering_rate(int(result[1]),int(result[2]),ref_st, ref_end,strand)
		return result, covering_per


	def write_proc(self, method):
		data_result = []

		for x in self.bed:

			read_name_arr = []
			temp = []

			for bc in self.bed:
				if x[0] == bc[0]:
						temp.append(bc)

			covr = []
			covv = 0.0
			for read in self.bamfile.fetch(x[0],int(x[1]),int(x[2])):
				if read.is_unmapped==False:
					result = "";
					if method == "D":
						result = self.distance_measure(read, temp)

					elif method == "C":
						result, covv = self.covered_rate_measure(read, temp)

					elif method == "A":
						result, covv = self.amp_covered_rate_measure(read, temp)

					else:
						print "Please, assign the method for calculation"
						exit()

					if result == x:
						#read_name_arr.append(read.query_name)

						temp_cig = read.cigar
						count = 0
						for y in temp_cig:
							if int(y[0])==0 or int(y[0])==1 or int(y[0])==2:
								count += int(y[1])

						if count >= 40:
						#if read.reference_length >= 40 and covv > 0.8 and read.mapping_quality>60.0:
							read_name_arr.append(read.query_name)
							covr.append(covv)


						#if read.reference_length >= 40:


			read_name_edit_arr = list(set(read_name_arr))
			if covr:
				covr_res = numpy.mean(covr)
			else:
				covr_res = "NA"

			data_result.append([self.barcode,x[0],x[1],x[2],len(read_name_edit_arr), covr_res])

		return data_result

	def __init__(self,fn, bed_fn):
		self.bamfile = pysam.AlignmentFile(fn,"rb")
		self.barcode = self.file_name_pro(fn)
		self.bed = self.bed_pro(bed_fn)
		#self.rw = open(rw,"w")



	def __del__(self):
		self.bamfile.close()

