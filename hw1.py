import sys
import matplotlib.pyplot as plot
import os
import numpy as np
import scipy.stats

##To run from terminal 'python3 hw1.py wang_data.txt'
file = sys.argv[1]
text =  open(file).read().split('\n')
matrix = [i.split("\t") for i in text]

###########################################
##Part 2
###########################################

def part2():
    print("Number of probes: ", len(text)-3)
    vals = [0,0]
    line1 = text[1].split('\t')
    for i in line1:
        if(i == "0"):
            vals[0] +=1
        else:
            vals[1] += 1
    print("Number of Relapse Free: ", vals[0], " Number of Relapses: ", vals[1],
    "Number of Total Patients: ", vals[1]+vals[0])


    genes = set()


    for line in text[2::]:
        words = line.split("\t")
        if(len(words)> 1):
            genes.add(words[1])
    print("Number of unique Genes:", len(genes))

part2()
###########################################
##Part 3
###########################################

def part3():
    ##I'm confused on what is required to do here.

    buckets = []
    for line in matrix[2:]:
        for val in line[2:]:
            buckets.append(float(val))

    # plot.hist(buckets, 10)
    #
    # luckets = [np.log10(i) for i in buckets]
    # plot.hist(luckets, 10)
    # plot.show()

    buckets = [[],[],[],[]]
    for line in range(2, len(matrix)-1):
        for j in range(2, 6):
            val  = np.log10(float(matrix[line][j]))
            buckets[j-2].append(val )
    # for i in range(4):
    #     plot.hist(buckets[i], 20)
    #     plot.show()

    means = []
    for line in matrix[2::]:
        vals = [float(i) for i in line[2:]]
        sum = 0
        for i in vals: sum += i
        if(len(line[2:])> 0):
            mean = sum/len(line[2:])
            means.append(mean)
    means.sort()
    
    cols = [ [] for i in range(288)]
    for i in range(2,len(matrix)):
        if(i % 1000 == 0):
            print(i)
        for j in range(2, len(matrix[i])):
            cols[j].append(float(matrix[i][j]))
        for c in cols:
            c.sort()
            for k in range(len(c)):
                c[k] = means[k]

    for i in range(4):
        plot.hist(buckets[i], 20)
        plot.show()
part3()


###########################################
##Part 4
###########################################

def part4():

    ##This is for the ttest and maybe works
    dflst = []
    ## For each probe we find the t-value
    for probe in range(2,len(matrix)-1):
        laps = []
        unlps = []
        ## Get the expression level for non-lapsed and relapsed samples.
        for g in range(2,len(matrix[probe])):
            if(int(matrix[1][g]) == 1):
                laps.append(float(matrix[probe][g]))
            else:
                unlps.append(float(matrix[probe][g]))
        results = scipy.stats.ttest_ind(laps, unlps, equal_var = False)[1]
        dflst += [(results, probe)]

    scipy.stats.wilcoxon(laps, unlps)
    dflst.sort(reverse=True, key = lambda x : x[0])
    for probe in dflst[:10]:
        print(matrix[probe[1]][1])

    dflst = [-np.log10(i[0]) for i in dflst]
    plot.hist(dflst)
    plot.show()
        ## This is so the comments collapse with the function
        ## It's meaningless
    a = 5


###########################################
##Part 5
###########################################
