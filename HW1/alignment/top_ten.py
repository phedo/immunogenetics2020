import re
import operator
import needleman_wunsch
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
f = open('plasma.fasta', 'r')

read = ''
read_value = ''
gene = ''
gene_value = ''

gene_reads_dict = {}
gene_gene_dict = {}
gene_vstreca_dict = {}

for line in f:
	if read == '':
		read = line
		continue
	elif read_value == '':
		read_value = line.rstrip()
		continue
	elif gene == '':
		gene = re.findall(r'GENE\:([A-Z0-9\-]*)\*', line)[0]
		continue
	elif gene_value == '':
		gene_value = line.rstrip()

	if gene in gene_vstreca_dict.keys():
		gene_vstreca_dict[gene] = gene_vstreca_dict[gene] + 1
		gene_reads_dict[gene].append(read_value)
	else:
		gene_vstreca_dict[gene] = 1
		gene_reads_dict[gene] = [read_value]
		gene_gene_dict[gene] = gene_value
		
	
	read = ''
	read_value = ''
	gene = ''
	gene_value = ''


sorted_gene_vstreca_dict = sorted(gene_vstreca_dict.items(), key=operator.itemgetter(1), reverse = True)


gene_diff_number = []
for i in range(10):
	gene = sorted_gene_vstreca_dict[i][0]
	gene_string = gene_gene_dict[gene]
	gene_reads = gene_reads_dict[gene]
	print(gene)
	for read in gene_reads:
		gene_diff_number.append([gene, needleman_wunsch.alignments(gene_string,read)[2]])
		print(needleman_wunsch.alignments(gene_string,read)[2])


#print(gene_string)
#print('\n\n')
#print(gene_reads)


arr = np.array([[1,110],[1,120],[1,130],[2,200],[2,220]])
df = pd.DataFrame(gene_diff_number)

sns.boxplot(x=df[0], y=df[1])
plt.show()
