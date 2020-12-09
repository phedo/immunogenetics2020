import re
import sys
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

# groop cdr3s by length
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

def compare_two_strings(string1, string2, max_mismatch):
	l = len(string1)
	mismatch_number = 0
	for i in range(l):
		if string1[i] != string2[i]:
			mismatch_number = mismatch_number + 1
			if mismatch_number > max_mismatch:
				break
	if mismatch_number <= max_mismatch:
		return True
	return False


def compare_two_dicts(dict1, dict2, max_mismatch):
	similar_elements = []
	for el1 in dict1.values():
		for el2 in dict2.values():
			if compare_two_strings(el1, el2, max_mismatch):
				similar_elements.append(el1)
				similar_elements.append(el2)
				dict1[list(dict1.keys())[0]].remove(el1)
				dict2[list(dict2.keys())[0]].remove(el2)
	return (dict1[list(dict1.keys())[0]], dict2[list(dict2.keys())[0]], similar_elements)

# compare pairs of similar (by length and ACGT-content) cdr3s
treshold = 0.2
for cdr3s_length in cdr3s_length_dict:
	max_mismatch = int(2*cdr3s_length*treshold)
	acgt_cdr3s_dict = cdr3s_length_dict[cdr3s_length]
	while len(acgt_cdr3s_dict)>0:
		acgt = list(acgt_cdr3s_dict.keys())[0]
		deltas = deltas(acgt, max_mismatch)
		lineage = []
		for neighbour in deltas:
			if neighbour in acgt_cdr3s_dict:
				acgt_cdr3s_dict[acgt], acgt_cdr3s_dict[neighbour], lineage_add = compare_two_dicts(acgt_cdr3s_dict[acgt], acgt_cdr3s_dict[neighbour], max_mismatch)
				lineage = lineage + lineage_add
				





				for who in acgt_cdr3s_dict[acgt]:
					for whom in acgt_cdr3s_dict[neighbour]:
						mismatch_number = 0
						for i in range(cdr3s_length):
							if who[k] != whom[k]:
								mismatch_number = mismatch_number + 1
								if mismatch_number > max_mismatch:
									break
						if mismatch_number <= max_mismatch:
							res.append(r)
							cdr3s_arr_new.remove(r_full)













	#cdr3s_arr = cdr3s_length_dict[cdr3s_length]
	acgt_cdr3s_dict = cdr3s_length_dict[cdr3s_length]
	first_key = list(acgt_cdr3s_dict.keys())[0]
	first_cdr3 = acgt_cdr3s_dict[first_key].pop(0)
	if len(acgt_cdr3s_dict[first_key]) == 0:
		acgt_cdr3s_dict.pop(first_key)
	lineages = []
	while len(acgt_cdr3s_dict) > 0:
		cdr3 = cdr3s_arr.pop(0)[4]
		lineage, cdr3s_arr = f(cdr3, cdr3s_arr, [cdr3], treshold)
		lineages.append(lineage)

print(v)
print(datetime.datetime.now() - begin_time)