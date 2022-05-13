# input: a list of protein names
import pandas
import numpy
# replace 3-based aa to 1
aa ={"Ala":"A","Gly":"G","Ile":"I","Leu":"L","Pro":"P","Val":"V","Phe":"F","Trp":"W","Tyr":"Y","Asp":"D","Glu":"E","Arg":"R","His":"H","Lys":"K","Ser":"S","Thr":"T","Cys":"C","Met":"M","Asn":"N","Gln":"Q"}
with open("/Users/tshe0033/Desktop/change_aa_code/protein_names.txt") as file:
  for line in file:
    line = line.rstrip()
    for word, initial in aa.items():
      line = line.replace(word, initial)
    print(line)
    # returns only variants, get rid of  accession id
    #  only_var = line.split("p.",1)[1]
    #print(only_var)
    


