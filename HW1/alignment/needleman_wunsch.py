import numpy as np

d = -1 # gap penalty

def S(a, b): # Similarity_matrix s
	if a == '-' or b == '-':
		return 0
	elif a == b:
		return 1
	else:
		return -1 

#print(S('d','d'))


def F_matrix(str1='GATTACA', str2='GCATGCU'): #F matrix with scores
	len1 = len(str1)
	len2 = len(str2)
	F = np.zeros((len1+1, len2+1))
	for i in range(len1+1):
		F[i,0] = d*i
	for j in range(len2+1):
		F[0,j] = d*j
	for i in range(len1):
		for j in range(len2):
			match = F[i,j] + S(str1[i], str2[j])
			delete = F[i,j+1] + d
			insert = F[i+1, j] + d
			F[i+1,j+1] = max(match,delete,insert)
	return F

#print(F_matrix('GATTACA','GCATGCU'))

	

def alignments(str1, str2): # return alignment(s) with maximum score
	differences_number = 0
	F = F_matrix(str1,str2)
	alignments_arr = []
	alignment1 = ''
	alignment2 = ''
	i = F.shape[0] - 1
	j = F.shape[1] - 1
	while (i>0 or j>0):
		differences_number = differences_number + 1
		if i>0 and j>0 and F[i,j] == F[i-1,j-1] + S(str1[i-1],str2[j-1]):
			alignment1 = str1[i-1] + alignment1
			alignment2 = str2[j-1] + alignment2
			i = i - 1
			j = j - 1
			if str1[i] == str2[j]:
				differences_number = differences_number - 1
		elif i>0 and F[i,j] == F[i-1,j] + d:
			alignment1 = str1[i-1] + alignment1
			alignment2 = '-' + alignment2	
			i = i - 1
		else:
			alignment1 = '-' + alignment1
			alignment2 = str2[j-1] + alignment2
			j = j - 1
	#print(alignment1)				
	#print(alignment2)
	#print(differences_number)
	return [alignment1, alignment2, differences_number]
		

#alignments('GATTACA','GCATGCU')

#alignments('krgagralsrk','dfg-lfa-fdga')
	
	
