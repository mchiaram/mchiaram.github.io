
from biblatex_dictionary import *
import sys
import re 

if __name__ == "__main__":

    bib = bib_dictionary(sys.argv[1:])
    
    entries = {}

    for key in bib.keys():
        
        for inner_key in bib[key].keys():
            bib[key][inner_key] = bib[key][inner_key].strip(" {},")

        if bib[key]['type'] == 'article':
            entry = '%s %s. /%s/, %s, %s. ' %( bib[key]['author'],\
                bib[key]['title'],
                bib[key]['journal'], 
                ('' + (bib[key]['volume'] if bib[key].has_key('volume')  else '') +\
                     ('('+bib[key]['number']+')'  if bib[key].has_key('number')  else '')+\
                     (   ':'+bib[key]['pages']  if bib[key].has_key('pages') else '') ),
                '%s'%bib[key]['year'])
            if bib[key].has_key('doi'):

                entry += '\[ [%s bib ]| [%s DOI ] | [%s http ] \]'%(
                ('journals_bib.html\#%s'%key),
                ('http://dx.doi.org/%s'%bib[key]['doi']),
                bib[key]['url'] )
            
        if bib[key]['type'] == 'incollection':
            entry = '%s %s. /%s/, %s, %s, %s. ' %( bib[key]['author'],\
                bib[key]['title'],
                bib[key]['booktitle'], 
                bib[key]['publisher'], 
                ('' + (bib[key]['volume'] if bib[key].has_key('volume')  else '') +\
                     ('('+bib[key]['number']+')'  if bib[key].has_key('number')  else '')+\
                     (   ':'+bib[key]['pages']  if bib[key].has_key('pages') else '') ),
                '%s'%bib[key]['year'])
            if bib[key].has_key('doi'):

                entry += '\[ [%s bib ]| [%s DOI ] | [%s http ] \]'%(
                ('journals_bib.html\#%s'%key),
                ('http://dx.doi.org/%s'%bib[key]['doi']),
                bib[key]['url'] )

        if 'Thesis' in bib[key]['type'] :
            entry = '%s /%s/. %s, %s, %s. ' %(
                bib[key]['author'],
                bib[key]['title'],
                bib[key]['type'], 
                bib[key]['school'], 
                '%s'%bib[key]['year'])
            entry += '\[ [%s http ] \]'%( bib[key]['url'] )
        
        if not entries.has_key( int(bib[key]['year']) ) :
            entries[int(bib[key]['year'])] = [] 

        entries[int(bib[key]['year'])].append(entry)

    for key in sorted(entries.iterkeys(),reverse=True):
        for entry in entries[key]:
            print ':{'+entry+'}\n'


    

