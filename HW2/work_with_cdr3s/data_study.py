import re
import sys

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

print(len(cdr3s_arr))

# function, calculating ACGT-content
def ACGT_content(string):
	return (len(re.findall(r'A',string)), len(re.findall(r'C',string)), len(re.findall(r'G',string)), len(re.findall(r'T',string)))

# function, calculting difference in ACGT-contents
def ACGT_diff(content1, content2):
	diff = tuple(abs(x-y) for x,y in zip(content1, content2))
	return sum(list(diff))

# groop cdr3s by length and CAGT-content
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

treshold = 0.2
number_of_sequences_total = 0
number_of_sequences_i_should_not_compare = 0
for length in cdr3s_length_dict:
	print('_______sequences_length='+str(length)+'____number_of_CAGT_groups='+str(len(cdr3s_length_dict[length]))+'________')
	for acgt1 in cdr3s_length_dict[length]:
		for acgt2 in cdr3s_length_dict[length]:
			number_of_sequences_total = number_of_sequences_total + len(cdr3s_length_dict[length][acgt1]) + len(cdr3s_length_dict[length][acgt2])  
			if ACGT_diff(acgt1, acgt2) > int(2*treshold*length):
				number_of_sequences_i_should_not_compare = number_of_sequences_i_should_not_compare + len(cdr3s_length_dict[length][acgt1]) + len(cdr3s_length_dict[length][acgt2])
	
'''for d in cdr3s_length_dict[el]:
		for l in cdr3s_length_dict[el][d]:
			print (l)'''
			
print(100*number_of_sequences_i_should_not_compare/number_of_sequences_total)		
#print(do_not_compare)
#print(counter)
#		print(str(key)+'    '+str(len(cdr3s_length_dict[el][key])))
