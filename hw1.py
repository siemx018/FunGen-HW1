import sys
import os

file = sys.argv[1]
text =  open(file).read().split('\n')

print("Number of probes: ", len(text)-2)
vals = [0,0]
line1 = text[1].split('\t')
for i in line1:
    if(i == "0"):
        vals[0] +=1
    else:
        vals[1] += 1
print("Number of Relapse Free: ", vals[0], " Number of Relapses: ", vals[1])


genes = set()

##I'm Suspicous about this part
for line in text[2::]:
    words = line.split("\t")
    if(len(words)> 1):
        genes.add(words[1])
print("Number of unique Genes:", len(genes))
