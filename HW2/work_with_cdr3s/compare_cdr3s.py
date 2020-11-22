import re
import sys
import datetime
#sys.setrecursionlimit(100000000)
begin_time = datetime.datetime.now()

global v
v = 0

index = ''
cdr3s_arr = []

f = open('cdr3s.fasta', 'r')
# read all cdr3s from file
for line in f:
	if line[0] == '>':
		index = line
		cdr3 = ''
		continue
	elif cdr3 == '' and index != '':
		cdr3 = line.rstrip()
		cdr3s_arr.append(cdr3)
		continue
	cdr3s_arr[-1] = cdr3s_arr[-1] + line.rstrip()

# groop cdr3s by length
cdr3s_length_dict = {}
	
for cdr3 in cdr3s_arr:
	cdr3_length = len(cdr3)
	if cdr3_length in cdr3s_length_dict.keys():
		cdr3s_length_dict[cdr3_length].append(cdr3)
	else:
		cdr3s_length_dict[cdr3_length] = [cdr3]

# function, calculating ACGT-content
def ACGT_content(string):
	return (len(re.findall(r'A',string)), len(re.findall(r'C',string)), len(re.findall(r'G',string)), len(re.findall(r'T',string)))
	
# function, calculting difference in ACGT-contents
def ACGT_diff(content1, content2):
	diff = tuple(abs(x-y) for x,y in zip(content1, content2))
	return sum(list(diff)) / sum(list(content1))

# function, returns all elements of a given array,
# whitch are close enough to a given element

def close_cdr3s(cdr3, cdr3s_arr, treshold):
	global v
	res = []
	lg = len(cdr3)
	cdr3s_arr_new = cdr3s_arr.copy()
	for r in cdr3s_arr:
		v=v+1
		if ACGT_diff(ACGT_content(cdr3),ACGT_content(r)) < 2*treshold*lg:
			mismatch_number = 0
			for k in range(lg):
				if cdr3[k] != r[k]:
					mismatch_number = mismatch_number + 1
					if mismatch_number > treshold*lg:
						break
			if mismatch_number <= treshold*lg:
				res.append(r)
				cdr3s_arr_new.remove(r)
	return (cdr3s_arr_new,res)

# function, returns HG with a given element,
# i.e. returns a lineage with a first cdr3 in a list

def f(cdr3, cdr3s_arr, lineage, treshold):
	cdr3s_arr, closeCDR3s = close_cdr3s(cdr3, cdr3s_arr, treshold)
	for r in closeCDR3s:
		lineage.append(r)
		lineage, cdr3s_arr = f(r, cdr3s_arr, lineage, treshold)
	return (lineage, cdr3s_arr)

# compare pairs of similar (by length and ACGT-content) cdr3s
treshold = 0.2
lineages = []
for cdr3s_length in cdr3s_length_dict:
	cdr3s_arr = cdr3s_length_dict[cdr3s_length]
	while len(cdr3s_arr) > 0:
		cdr3 = cdr3s_arr.pop(0)
		lineage, cdr3s_arr = f(cdr3, cdr3s_arr, [cdr3], treshold)
		lineages.append(lineage)

print(v)
print(datetime.datetime.now() - begin_time)
