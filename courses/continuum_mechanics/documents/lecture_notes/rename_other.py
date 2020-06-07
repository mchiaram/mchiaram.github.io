import os 
import sys
from functools import reduce
if __name__ == '__main__':

	files = os.listdir('./')

	for afile in files:

		p = afile.split('.')
		
		if p[-1] == 'pdf':
			bfile = afile.replace('-','')
			bfile = bfile.replace(' ','_')
			bfile = bfile.replace('__','_')
			bfile.split('_')
			bfile = reduce(lambda y,x: y+'_'+x, bfile.split('_')[:-2] )+'.pdf'
			print('- [documents/lecture_notes/'+bfile, reduce(lambda y,x: y+' '+x, bfile.strip('.pdf').split('_')[2::] ),']' )
			



	
