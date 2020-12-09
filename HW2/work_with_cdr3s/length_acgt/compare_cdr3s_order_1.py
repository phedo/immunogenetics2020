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

# function, calculating ACGT-content
def ACGT_content(string):
	return (len(re.findall(r'A',string)), len(re.findall(r'C',string)), len(re.findall(r'G',string)), len(re.findall(r'T',string)))

# groop cdr3s by length and acgt-content
cdr3s_length_dict = {}
	
for cdr3 in cdr3s_arr:
	cdr3_length = len(cdr3)
	a, c, g, t = ACGT_content(cdr3)
	if cdr3_length in cdr3s_length_dict.keys():
		acgt_cdr3s_dict = cdr3s_length_dict[cdr3_length]
		if (a, c, g, t) in acgt_cdr3s_dict.keys():
			cdr3s_length_dict[cdr3_length][(a, c, g, t)].append(cdr3)
		else:
			cdr3s_length_dict[cdr3_length][(a, c, g, t)] = [cdr3]
	else:
		acgt_cdr3s_dict = {}
		acgt_cdr3s_dict[(a, c, g, t)] = [cdr3]
		cdr3s_length_dict[cdr3_length] = acgt_cdr3s_dict


# function, returns all elements of a given array,
# whitch are close enough to a given element

def close_cdr3s(cdr3, cdr3s_arr, treshold):
	global v
	res = []
	lg = len(cdr3)
	cdr3s_arr_new = cdr3s_arr.copy()
	for r_full in cdr3s_arr:
		v=v+1
		r = r_full[4]
		if ACGT_diff(ACGT_content(cdr3),ACGT_content(r)) < 2*treshold*lg:
			mismatch_number = 0
			for k in range(lg):
				if cdr3[k] != r[k]:
					mismatch_number = mismatch_number + 1
					if mismatch_number > treshold*lg:
						break
			if mismatch_number <= treshold*lg:
				res.append(r)
				cdr3s_arr_new.remove(r_full)
	return (cdr3s_arr_new,res)

# function, returns HG with a given element,
# i.e. returns a lineage with a first cdr3 in a list

def f(cdr3, cdr3s_arr, lineage = [], treshold = 0.2):
	cdr3s_arr, closeCDR3s = close_cdr3s(cdr3, cdr3s_arr, treshold)
	for r in closeCDR3s:
		lineage.append(r)
		lineage, cdr3s_arr = f(r, cdr3s_arr, lineage, treshold)
	return (lineage, cdr3s_arr)

# function, returns ACGT-contents, which are close enough for a given acgt, 
# considering treshold value
# for the realisation simplicity this ACGT-contents are excess.

def deltas(acgt, max_mismatch):
	a,c,g,t = acgt
	variants = []
	for D in range(max_mismatch+2):
		for i in range(D):
			for j in range(max(0,D-i)):
				for k in range(max(0,D-i-j)):
					for l in range(max(0,D-i-j-k)):
						variants.append((i,j,k,l))
		variants = list(dict.fromkeys(variants))
		res = []
		for el in variants:
			i = el[0]
			j = el[1]
			k = el[2]
			l = el[3]
			res.append((a+i,c+j,g+k,t+l))
			res.append((a-i,c+j,g+k,t+l))
			res.append((a+i,c-j,g+k,t+l))
			res.append((a+i,c+j,g-k,t+l))
			res.append((a+i,c+j,g+k,t-l))
			res.append((a-i,c-j,g+k,t+l))
			res.append((a-i,c+j,g-k,t+l))
			res.append((a-i,c+j,g+k,t-l))
			res.append((a+i,c-j,g-k,t+l))
			res.append((a+i,c-j,g+k,t-l))
			res.append((a+i,c+j,g-k,t-l))
			res.append((a-i,c-j,g-k,t+l))
			res.append((a-i,c-j,g+k,t-l))
			res.append((a-i,c+j,g-k,t-l))
			res.append((a+i,c-j,g-k,t-l))
			res.append((a-i,c-j,g-k,t-l))		
	return list(dict.fromkeys(res))



# function, returns all elements of a given array,
# whitch are close enough to a given element

def close_data(cdr3, cdr3s_to_compare, treshold):
	res = []
	lg = len(cdr3)
	cdr3s_to_compare_new = cdr3s_to_compare.copy()
	for r in cdr3s_to_compare:
		mismatch_number = 0
		for k in range(lg):
			if cdr3[k] != r[k]:
				mismatch_number = mismatch_number + 1
				if mismatch_number > treshold*lg:
					break
		if mismatch_number <= treshold*lg:
			res.append(r)
			cdr3s_to_compare_new.remove(r)
	return (cdr3s_to_compare_new,res)

# function, returns HG with a given element,
# i.e. returns a lineage with a first cdr3 in a list

def f(first_cdr3, cdr3s_to_compare, lineage, treshold):
	cdr3s_to_compare, closeCDR3s = close_data(first_cdr3, cdr3s_to_compare, treshold)
	for r in closeCDR3s:
		lineage.append(r)
		lineage, cdr3s_to_compare = f(r, cdr3s_to_compare, lineage, treshold)
	return (lineage, cdr3s_to_compare)


def compare_cdr3s(cdr3s_to_compare, treshold):
	lineages = []
	while len(cdr3s_to_compare)>0:
		first_cdr3 = cdr3s_to_compare.pop(0)
		lineage = [first_cdr3]
		if len(cdr3s_to_compare)>0:
			lineage, cdr3s_to_compare = f(first_cdr3, cdr3s_to_compare, lineage, treshold)
		lineages.append(lineage)
	return lineages

# compare pairs of similar (by length and ACGT-content) cdr3s
treshold = 0.1
lineages = []
for cdr3s_length in cdr3s_length_dict:
	max_mismatch = int(2*cdr3s_length*treshold)
	acgt_cdr3s_dict = cdr3s_length_dict[cdr3s_length]
	while len(acgt_cdr3s_dict)>0:
		acgt = list(acgt_cdr3s_dict.keys())[0]
		dlts = deltas(acgt, max_mismatch)
		cdr3s_to_compare = []
		for arg in dlts:
			if arg in acgt_cdr3s_dict.keys():
				cdr3s_to_compare = cdr3s_to_compare + acgt_cdr3s_dict[arg]
				acgt_cdr3s_dict.pop(arg)
		lineages = lineages + compare_cdr3s(cdr3s_to_compare, treshold)
print('The number of clonal lineages = ' + str(len(lineages)))

print(v)
print(datetime.datetime.now() - begin_time)
