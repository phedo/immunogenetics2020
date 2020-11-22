import re
import sys
import datetime
from operator import itemgetter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
begin_time = datetime.datetime.now()

# read all cdr3s from file

index = ''
cdr3 = ''
# d structure will be: {index: [cdr3, read, v-gene]}
d = {}
f = open('cdr3s.fasta', 'r')
for line in f:
	if line[0] == '>':
		index = re.findall(r'INDEX\:([A-Z0-9\-]*)\|', line)[0]
		cdr3 = ''
		continue
	elif cdr3 == '' and index != '':
		cdr3 = line.rstrip()
		d[index] = [cdr3]
		continue
	d[index] = [d[index][0] + line.rstrip()]
f.close()

# read all reads from file

f = open('PRJNA324093_Dnr4_10k.fasta', 'r')
index = ''
read = ''
for line in f:
	if index == '':
		index = re.findall(r'INDEX\:([A-Z0-9\-]*)\|', line)[0]
		continue
	elif read == '':
		read = line.rstrip()
		d[index].append(read)
	index = ''
	read = ''
f.close()

# read all V genes from file

f = open('v_alignment.fasta', 'r')
index = ''
read_value = ''
gene = ''
gene_value = ''
for line in f:
	if index == '':
		index = re.findall(r'INDEX\:([A-Z0-9\-]*)\|', line)[0]
		continue
	elif read_value == '':
		read_value = line.rstrip()
		continue
	elif gene == '':
		gene = re.findall(r'GENE\:([A-Z0-9\-]*)\*', line)[0]
		d[index].append(gene)
		continue
	elif gene_value == '':
		gene_value = line.rstrip()
	index = ''
	read_value = ''
	gene = ''
	gene_value = ''
f.close()

# groop cdr3s by length

cdr3s_length_dict = {}
for key in d:
	cdr3 = d[key][0]
	cdr3_length = len(cdr3)
	if cdr3_length in cdr3s_length_dict.keys():
		cdr3s_length_dict[cdr3_length].append(d[key])
	else:
		cdr3s_length_dict[cdr3_length] = [d[key]]



# function, returns all elements of a given array,
# whitch are close enough to a given element

def close_data(data, data_with_one_cdr3_length, treshold):
	res = []
	lg = len(data[0])
	data_with_one_cdr3_length_new = data_with_one_cdr3_length.copy()
	for r in data_with_one_cdr3_length:
		mismatch_number = 0
		for k in range(lg):
			if data[0][k] != r[0][k]:
				mismatch_number = mismatch_number + 1
				if mismatch_number > treshold*lg:
					break
		if mismatch_number <= treshold*lg:
			res.append(r)
			data_with_one_cdr3_length_new.remove(r)
	return (data_with_one_cdr3_length_new,res)

# function, returns HG with a given element,
# i.e. returns a lineage with a first cdr3 in a list

def f(data, data_with_one_cdr3_length, lineage, treshold, lineage_number):
	data_with_one_cdr3_length, closeData = close_data(data, data_with_one_cdr3_length, treshold)
	for r in closeData:
		lineage.append(r)
		lineage, data_with_one_cdr3_length = f(r, data_with_one_cdr3_length, lineage, treshold, lineage_number)
	return (lineage, data_with_one_cdr3_length)

# compare pairs of similar (by length) cdr3s
treshold = 0.2
lineages = []
max_lineage = []
max_lineage_length = 0
at_least_ten_counter = 0
lineage_number = 0
lineage_v_gene = {}
for cdr3s_length in cdr3s_length_dict:
	data_with_one_cdr3_length = cdr3s_length_dict[cdr3s_length]
	while len(data_with_one_cdr3_length) > 0:
		data = data_with_one_cdr3_length.pop(0)
		lineage_number = lineage_number+1
		print('lineage_number = ' + str(lineage_number))
		lineage, data_with_one_cdr3_length = f(data, data_with_one_cdr3_length, [data], treshold, lineage_number)
		lineages.append(lineage)
		if lineage[0][2] in lineage_v_gene.keys():
			lineage_v_gene[lineage[0][2]] = lineage_v_gene[lineage[0][2]] + 1
		else:
			lineage_v_gene[lineage[0][2]] = 1
		if len(lineage) > max_lineage_length:
			max_lineage_length = len(lineage)
			max_lineage = lineage
		if len(lineage) > 9:
			at_least_ten_counter = at_least_ten_counter + 1

fasta_output = open('max_lineage_data.fasta', 'a')
for data in max_lineage:
	fasta_output.write('>\n')
	fasta_output.write(data[1] + '\n')
fasta_output.close()

print('The number of clonal lineages = ' + str(len(lineages)))
print('The number of sequences in the largest lineage = ' + str(max_lineage_length))
print('The number of clonal lineages presented by at least 10 sequences = ' + str(at_least_ten_counter))

v_gene_lineages = []
for key, value in lineage_v_gene.items():
    temp = [key,value]
    v_gene_lineages.append(temp)

df = pd.DataFrame(sorted(v_gene_lineages, key=itemgetter(1)))
#sns.boxplot(x=df[0], y=df[1])
ax = sns.boxplot(x=df[0], y=df[1], data=df)
plt.setp(ax.get_xticklabels(), rotation=90)
plt.show()

print(datetime.datetime.now() - begin_time)
