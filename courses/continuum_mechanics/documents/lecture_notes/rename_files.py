

import os 

files = os.listdir('./')
print(files)
for afile in files[1:]: 
	afile = afile.replace('\& ','')
	print(afile)
	command = 'mv \"%s\" \"%s\"'%( afile, afile.replace(' ','') )
	print( command )
	os.system( command )


