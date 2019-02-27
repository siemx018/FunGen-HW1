import sys
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
    print("Number of probes: ", len(text)-2)
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

###########################################
##Part 3
###########################################

def part3():
    ##I'm confused on what is required to do here.
    hist = dict()
    # for line in text[2::]:
    #     cols = line.split("\t")
    #     if(len(cols)>=1):
    #         for col in cols[2::]:
    #             value = float(col)
    #             if(value in hist):
    #                 hist[value] += 1
    #             else:
    #                 hist[value] = 1


###########################################
##Part 4
###########################################

def part4():
    #### This might be a helpful function to get the first n biggest numbers
    #### For this question we need to look at the 10 most differential probes
    def getNbig(arr, n):
        final_list = []
        for i in range(0, n):
            mx = 0
            for j in arr:
                if j > mx:
                    mx = j

            list1.remove(mx)
            final_list.append(mx)

        return final_list

    dflst = []
    ## For each probe
    for probe in range(len(matrix[2::])):
        laps = []
        unlps = []
        ## Get the expression level for non-lapsed and relapsed samples.
        for g in range(len(matrix[probe][2::])):
            if(matrix[1][g] == 1):
                laps += [g]
            else:
                unlps += [g]
        ## scipy has a built in wilcoxon test if we don't want to implement it
        ##But I am confused here because this method wants to equal sized data sets
        ## Which we don't have. The wiki page of the test also assumes to equal sized
        ##Data sets.

        #score, p  = scipy.stats.wilcoxon(laps, unlps)
        #dflst += [score]
        ## This is so the comments collapse with the function 
    a = 5


###########################################
##Part 5
###########################################
