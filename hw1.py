import matplotlib.pyplot as plot
import os
import numpy as np
import scipy.stats

##To run from terminal 'python3 hw1.py wang_data.txt'
text =  open('wang_data.txt').read().split('\n')
matrix = [i.split("\t") for i in text]

SHOW_GRAPHS = False
###########################################
##Part 2
###########################################

def part2():
    print("Number of probes: ", len(text)-3)
    vals = [0,0]
    line1 = matrix[1]
    for i in line1:
        if(i == "0"):
            vals[0] +=1
        else:
            vals[1] += 1
    print("Number of Relapse Free: ", vals[0], " Number of Relapses: ", vals[1],
    "Number of Total Patients: ", vals[1]+vals[0])


    genes = set()
    for row in matrix[2::]:
        if(len(row) > 1):
            genes.add(row[1])
    print("Number of unique Genes:", len(genes))

part2()
###########################################
##Part 3
###########################################

def part3():

    ##Part A
    buckets = []
    for line in matrix[2:]:
        for val in line[2:]:
            buckets.append(float(val))


    plot.hist(buckets, 100)
    if(SHOW_GRAPHS):
        plot.show()
    luckets = [np.log10(i) for i in buckets]
    plot.hist(luckets, 100)
    if(SHOW_GRAPHS):
        plot.show()


    ## Part B
    buckets = [[],[],[],[]]
    for line in range(2, len(matrix)-1):
        for j in range(2, 6):
            val  = np.log10(float(matrix[line][j]))
            buckets[j-2].append(val )
    for i in range(4):
         plot.hist(buckets[i], 100)
         if(SHOW_GRAPHS):
             plot.show()

    ## Part C
    means = []
    for line in matrix[2::]:
        vals = [float(i) for i in line[2:]]
        sum = 0
        for i in vals: sum += i
        if(len(line[2:])> 0):
            mean = sum/len(line[2:])
            means.append(mean)
    means.sort()

    trans = scipy.mat(matrix[2::][2::]).T




###########################################
##Part 4
###########################################

def part4():

    ##Part a
    ##This is for the ttest
    dflst = []
    willist = []
    setd = set()
    setw = set()
    ## For each probe we find the t-value

    ##One of these is probably wrong
    for probe in range(2,len(matrix)-1):
        s1 = []
        s2 = []

        dcount = 0
        wcount = 0
        ## Get the expression level for non-lapsed and relapsed samples.
        for g in range(2,len(matrix[probe])):
            if(int(matrix[1][g]) == 1):
                s1.append(float(matrix[probe][g]))
            else:
                s2.append(float(matrix[probe][g]))
        results = scipy.stats.ttest_ind(s1, s2)[1]
        dflst += [(results, probe)]
        if(results < 0.05):
            setd.add(matrix[probe][1])
        results = scipy.stats.mannwhitneyu(s1, s2)[1]
        willist += [(results, probe)]
        if(results < 0.05):
            setw.add(matrix[probe][1])
    ##End Loop

    willist.sort( key = lambda x : x[0])
    dflst.sort( key = lambda x : x[0])
    for probe in dflst[:10]:
        print(probe[0], end=" ")
        print(matrix[probe[1]][1])
    print("\n \n")
    for probe in willist[1:10]:
        print(probe[0], end=" ")
        print(matrix[probe[1]][1])

    print("# of Signifigant from ttest:", len(setd))
    print("# of Signifigant from rank-sum:", len(setw))
    print("Size of Set difference", len(setd.difference(setw)))

    wilog = [-np.log10(i[0]) for i in dflst if i[0] > 0]
    dflst = [-np.log10(i[0]) for i in dflst if i[0] > 0]
    if(SHOW_GRAPHS):
        plot.hist(dflst)
        plot.show()
        ## This is so the comments collapse with the function
        ## It's meaningless
    return willist

ranked = part4()

###########################################
##Part 5
###########################################


def part5(willist):

    ##Part A
    filterd = []
    for ele in willist:
        p, probe = ele[0], ele[1]
        corrected = p*(len(willist))
        if(corrected < 0.05):
            filterd.append((p,probe))
    print(len(filterd))

    #Part B
    FDR = 0.05
    M = len(willist)
    ffilter = 0
    for i in range(M):
        rank = i
        score = (rank/M)*FDR
        if( willist[i][0] < score):
            ffilter += 1
    print(ffilter)
part5(ranked)
