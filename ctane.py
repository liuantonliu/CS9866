"""------------------------------------------------------------------------------------------
TANE Algorithm for discovery of exact conditional functional dependencies
Author: Nabiha Asghar, nasghar@uwaterloo.ca
March 2015
Use for research purposes only.
Please do not re-distribute without written permission from the author
Any commerical uses strictly forbidden.
Code is provided without any guarantees.
----------------------------------------------------------------------------------------------"""
from pandas import *
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
import sys
import time

def replace_element_in_tuple(tup, elementindex, elementval):
    if type(elementval)==tuple:
        elementval = elementval[0]
    newtup = list(tup)
    newtup[elementindex] = elementval
    newtup = tuple(newtup)
    return newtup

def add_element_in_tuple(spxminusa, ca):
    thelist = list(spxminusa)
    thelist.append(ca[0])
    return tuple(thelist)

#! here
def validApproxcfd(xminusa, x, a, spxminusa, sp, ca):
    global dictpartitions
    global eps
    if xminusa == '' or a == '': 
        return False
    indexofa = x.index(a)
    newsp0 = add_element_in_tuple(spxminusa, ca)
    newsp1 = replace_element_in_tuple(sp, indexofa, ca)   #this is sp, except that in place of value of a we put ca
    if (x, newsp1) in dictpartitions.keys():
        bigLen = max(len(dictpartitions[(xminusa, spxminusa)]), len(dictpartitions[(x, newsp1)]))
        if  abs(len(dictpartitions[(xminusa, spxminusa)]) - len(dictpartitions[(x, newsp1)])) < eps*bigLen: # changed this
            return True    
    return False

def validcfd(xminusa, x, a, spxminusa, sp, ca):
    global dictpartitions
    if xminusa == '' or a == '': 
        return False
    indexofa = x.index(a)
    newsp0 = add_element_in_tuple(spxminusa, ca)
    newsp1 = replace_element_in_tuple(sp, indexofa, ca)   #this is sp, except that in place of value of a we put ca
    if (x, newsp1) in dictpartitions.keys():
        if  len(dictpartitions[(xminusa, spxminusa)]) == len(dictpartitions[(x, newsp1)]):# and twodlen(dictpartitions[(xminusa, spxminusa)]) == twodlen(dictpartitions[(x, newsp1)]):
            return True    
    return False

def twodlen(listoflists):
	summ = 0
	for item in listoflists:
		summ = summ + len(item)
	return summ

def greaterthanorequalto(upxminusa, spxminusa): # this is actually greaterthan or equal to
    if upxminusa == spxminusa: 
        return True
    flag = True
    for index in range(0, len(upxminusa)):
        if not (spxminusa[index]=='--'):
            if (not (upxminusa[index] == spxminusa[index])):
                flag = False
    return flag

def doublegreaterthan(upxminusa, spxminusa): 
    if upxminusa == spxminusa: 
        return False
    flag = True
    for index in range(0, len(upxminusa)):
        if (not spxminusa[index]=='--'):
            if (not (upxminusa[index] == spxminusa[index])):
                flag = False
    return flag

#! here
def compute_dependencies(level, listofcols):
    global dictCplus
    global finallistofCFDs
    global listofcolumns
    global eps
    for (x,sp) in level:
        for a in x:
            for (att, ca) in dictCplus[(x, sp)]:
                if att == a:
                    newtup =  spXminusA(sp, x, a)      ### tuple(y for y in sp if not sp.index(y)==x.index(a)) # this is sp[X\A]                             
                    if eps > 0: #! added this line and changed line below
                        if validApproxcfd( x.replace(a,''), x, a, newtup, sp, ca) and not ([x.replace(a,''), a, [newtup, ca]] in finallistofCFDs):
                            finallistofCFDs.append([x.replace(a,''), a, [newtup, ca]])
                            for (xx, up) in level:
                                if xx==x:
                                    newtup0 =  spXminusA(up, x, a)          ### tuple(y for y in up if not up.index(y)==x.index(a)) # this is up[X\A]
                                    if up[x.index(a)]==ca[0] and greaterthanorequalto(newtup0, newtup) :
                                        if (a, ca) in dictCplus[(x,up)]: dictCplus[(x,up)].remove((a,ca))
                                        #! added below line
                                        if validcfd( x.replace(a,''), x, a, newtup, sp, ca) and not ([x.replace(a,''), a, [newtup, ca]] in finallistofCFDs):
                                            listofcolscopy = listofcols[:]
                                            for j in x: # this loop computes R\X
                                                if j in listofcolscopy: listofcolscopy.remove(j)
                                            for b_att in listofcolscopy: # this loop removes each b in R\X from C+(X,up)
                                                stufftobedeleted = []
                                                for (bbval, sometup) in dictCplus[(x,up)]:
                                                    if b_att == bbval:
                                                        stufftobedeleted.append((bbval,sometup))                        
                                                for item in stufftobedeleted:
                                                    dictCplus[(x,up)].remove(item)
                    else: #! added this line, original code below
                        if validcfd( x.replace(a,''), x, a, newtup, sp, ca) and not ([x.replace(a,''), a, [newtup, ca]] in finallistofCFDs):
                            finallistofCFDs.append([x.replace(a,''), a, [newtup, ca]])
                            for (xx, up) in level:
                                if xx==x:
                                    newtup0 =  spXminusA(up, x, a)          ### tuple(y for y in up if not up.index(y)==x.index(a)) # this is up[X\A]
                                    if up[x.index(a)]==ca[0] and greaterthanorequalto(newtup0, newtup) :
                                        if (a, ca) in dictCplus[(x,up)]: dictCplus[(x,up)].remove((a,ca))
                                        listofcolscopy = listofcols[:]
                                        for j in x: # this loop computes R\X
                                            if j in listofcolscopy: listofcolscopy.remove(j)
                                        for b_att in listofcolscopy: # this loop removes each b in R\X from C+(X,up)
                                            stufftobedeleted = []
                                            for (bbval, sometup) in dictCplus[(x,up)]:
                                                if b_att == bbval:
                                                    stufftobedeleted.append((bbval,sometup))                        
                                            for item in stufftobedeleted:
                                                dictCplus[(x,up)].remove(item)

def prune(level):
    global dictCplus
    stufftobedeleted=[]
    for (x,sp) in level:
        if len(dictCplus[(x,sp)])==0:
            stufftobedeleted.append((x,sp))
    for item in stufftobedeleted:
        level.remove(item)

def computeCplus(level): # for each tuple (x,sp) in the list level, it computes C+(x,sp), which is a list of (attribute, value) tuples) 
    global listofcolumns
    global dictCplus
    listofcols = listofcolumns[:]
    for (x,sp) in level: #sp is a tuple of strings like this: ('aa', 'bb', 'cc') or ('aa', )     
        thesets=[]
        for b in x:
            indx = x.index(b) # the index where b is located in x
            spcopy =  spXminusA(sp, x, b)     ### tuple(y for y in sp if not sp.index(y)==indx)
            spcopy2 = sp[:]            
            if (x.replace(b,''), spcopy ) in dictCplus.keys():
                temp = dictCplus[(x.replace(b,''), spcopy)]
            else: temp = []   # is this correct???? should I put [] here?
            thesets.insert(0, set(temp))
        if list(set.intersection(*thesets)) == []:
            dictCplus[(x,sp)] = []
        else:
            dictCplus[(x,sp)] = list(set.intersection(*thesets))

def initial_Cplus(level):
    global listofcolumns
    global dictCplus
    computeCplus(level)
    for (a,ca) in level:
        stufftobedeleted = []
        for (att, val) in dictCplus[(a,ca)]:
            if att==a and not val==ca:
                stufftobedeleted.append((att,val))
        for item in stufftobedeleted:
            dictCplus[(a,ca)].remove(item)

def populateL1(listofcols):    
    global k_suppthreshold
    l1 = []
    attributepartitions = computeAttributePartitions(listofcols)
    for a in listofcols:
        l1.append((a, ('--',)))
        for eqclass in attributepartitions[a]:
            if len(eqclass)>= k_suppthreshold:
                l1.append( (a, (str(data2D.iloc[eqclass[0]][a]) , ) ) )
    computeInitialPartitions(l1, attributepartitions) # populates the dictpartitions with the initial partitions (X,sp) where X is a single attribute
    return l1

def computeInitialPartitions(level1, attributepartitions):
	global data2D
	global dictpartitions # dictpartitions[(x,sp)] is of the form [[0,1,2]]. So simply a list of lists of indices  
	for (a,sp) in level1:
		dictpartitions[(a,sp)]=[]
		dictpartitions[(a,sp)] = attributepartitions[a]

def old_computeInitialPartitions(level1, attributepartitions):
    global data2D
    global dictpartitions # dictpartitions[(x,sp)] is of the form [[0,1,2]]. So simply a list of lists of indices  
    for (a,sp) in level1:
        dictpartitions[(a,sp)]=[]
        if sp[0]=='--':
            dictpartitions[(a,sp)] = attributepartitions[a]
        else:
            for eqclass in attributepartitions[a]:
                if str(data2D.iloc[eqclass[0]][a])==sp[0]:
                    dictpartitions[(a,sp)].append(eqclass)

def computeAttributePartitions(listofcols): # compute partitions for every attribute 
    global data2D    
    attributepartitions = {}
    for a in listofcols:
        attributepartitions[a]=[]
        for element in list_duplicates(data2D[a].tolist()): # list_duplicates returns 2-tuples, where 1st is a value, and 2nd is a list of indices where that value occurs
            if len(element[1])>0: # if >1, then ignore singleton equivalence classes
                attributepartitions[a].append(element[1])
    return attributepartitions

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)>0)

def sometuplematchesZUP(z,up):
    global dictpartitions
    global k_suppthreshold
    sumofmatches = 0
    for eqclass in dictpartitions[(z, up)]:
        sumofmatches = sumofmatches +  len(eqclass)
    if sumofmatches >= k_suppthreshold:
        return True
    else:
        return False

def generate_next_level(level):
    nextlevel=[]
    for i in range(0,len(level)): # pick an element
        for j in range(i+1, len(level)): # compare it to every element that comes after it. 
            if ((not level[i][0]==level[j][0]) and level[i][0][0:-1]==level[j][0][0:-1] and level[i][1][0:-1]==level[j][1][0:-1]):
                z = level[i][0] + level[j][0][-1]
                up = tuple(list(level[i][1]) + [level[j][1][-1]])
                (z, up) = sortspbasedonx(z, up)
                partition_product((z,up), level[i], level[j])
                if sometuplematchesZUP(z,up):
                    flag = True
                    for att in z:
                        indexofatt = z.index(att) # where is att located in z                        
                        up_zminusa = spXminusA(up, z, att)
                        zminusa = z.replace(att,'')
                        if not ((zminusa, up_zminusa) in level):
                            flag = False
                    if flag:
                        nextlevel.append((z, up))
    return nextlevel

def spXminusA(sp, x, a):
    indexofa = x.index(a)
    mylist=[]
    for i in range(0, len(sp)):
        if not i==indexofa:
            mylist.append(sp[i])
    return tuple(mylist)

def partition_product(zup, xsp, ytp):
    global dictpartitions
    global tableT
    tableS = ['']*len(tableT)
    partitionXSP = dictpartitions[xsp]
    partitionYTP = dictpartitions[ytp]
    partitionZUP = []
    for i in range(len(partitionXSP)):
        for t in partitionXSP[i]:
            tableT[t] = i
        tableS[i]=''
    for i in range(len(partitionYTP)):
        for t in partitionYTP[i]: 
            if ( not (tableT[t] == 'NULL')): 
                tableS[tableT[t]] = sorted(list(set(tableS[tableT[t]]) | set([t]))) 
        for t in partitionYTP[i]: 
            if (not (tableT[t] == 'NULL')) and len(tableS[tableT[t]])>= 1 : 
                partitionZUP.append(tableS[tableT[t]]) 
            if not (tableT[t] == 'NULL'): tableS[tableT[t]]='' 
    for i in range(len(partitionXSP)): 
        for t in partitionXSP[i]: 
            tableT[t]='NULL'
    dictpartitions[zup] = partitionZUP
    dictpartitions[zup] = partitionZUP

def sortspbasedonx(x,sp):
    x = list(x)
    points = zip(x,sp)
    sorted_points = sorted(points)
    new_x = [point[0] for point in sorted_points]
    new_sp = [point[1] for point in sorted_points]
    return (''.join(new_x), tuple(new_sp))

#------------------------------------------------------- START ---------------------------------------------------

def getCFD(data2D, eps, k):
    start_time = time.time()
    global totaltuples
    global listofcolumns
    global tableT
    global k_suppthreshold
    global L0
    
    global dictCplus
    global dictpartitions
    global finallistofCFDs
    global L1
    
    totaltuples = len(data2D.index)
    listofcolumns = list(data2D.columns.values) # returns ['A', 'B', 'C', 'D', .....]
    tableT = ['NULL']*totaltuples # this is for the table T used in the function partition_product
    k_suppthreshold = k
    L0 = []

    dictpartitions = {} # maps 'stringslikethis' to a list of lists, each of which contains indices
    finallistofCFDs=[]
    L1=populateL1(listofcolumns[:])  # L1 is a list of tuples of the form [ ('A', ('val1') ), ('A', ('val2') ), ..., ('B', ('val3') ), ......]
    dictCplus = {('',()): L1[:]}
    l=1
    L = [L0,L1]

    while (not (L[l] == [])):
        if l==1:
            initial_Cplus(L[l])
        else:
            computeCplus(L[l])
        compute_dependencies(L[l],listofcolumns[:])
        prune(L[l])
        temp = generate_next_level(L[l])
        L.append(temp)
        l=l+1
        #print "List of all CFDs: " , finallistofCFDs
        #print "CFDs found: ", len(finallistofCFDs), ", level = ", l-1    

    # print("List of all CFDs: " , finallistofCFDs)
    print("Total number of CFDs found: ", len(finallistofCFDs))
    print("--- %s seconds ---" % round(time.time() - start_time, 5))
    return len(finallistofCFDs), round(time.time() - start_time, 5)


if __name__=="__main__":
    #! rows
    # rows = [5,10,20,50,100,200,500,1000]
    # CFDTime = []
    # ACFDTime = []
    # CFDCount = []
    # ACFDCount = []
    # # baseline: k=5,  eps=0.1, rows=100, columns = 5
    # for row in rows:
    #     data2D = pd.read_csv(r'D:\Western\CS9866\FD_CFD_extraction\adult.csv') #15 columns, 32561 rows
    #     data2D = data2D.head(row) # take first x rows
    #     data2D = data2D.iloc[:, : 4] # take first x+1 columns
        
    #     print("#########################################################") # normal CFD
    #     eps = 0
    #     k = 5
    #     numCFD, timeTaken = getCFD(data2D, eps, k)
    #     CFDTime += [timeTaken]
    #     CFDCount += [numCFD]
        
    #     print("*********************************************************") # Approx CFD
    #     eps = 0.1
    #     k = 5
    #     numCFD, timeTaken = getCFD(data2D, eps, k)
    #     ACFDTime += [timeTaken]
    #     ACFDCount += [numCFD]
        
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    # fig.suptitle('CFD vs ACFD Discovery - Rows')
    # ax1.scatter(rows, CFDTime, c="blue", label="CFD")
    # ax1.scatter(rows, ACFDTime, c="red", label="ACFD")
    # ax1.plot(rows, CFDTime, 'b-', alpha=0.5)  # Blue line for CFD
    # ax1.plot(rows, ACFDTime, 'r-', alpha=0.5)  # Red line for ACFD
    # ax1.set_xscale("log")
    # ax1.set_yscale("log")
    # ax1.legend()
    # ax1.set_xlabel("Number of rows in dataset")
    # ax1.set_ylabel("Time taken")
    
    # ax2.scatter(rows, CFDCount, c="blue", label="CFD")
    # ax2.scatter(rows, ACFDCount, c="red", label="ACFD")
    # ax2.plot(rows, CFDCount, 'b-', alpha=0.5)  # Blue line for CFD
    # ax2.plot(rows, ACFDCount, 'r-', alpha=0.5)  # Red line for ACF
    # ax2.set_xscale("log")
    # # ax2.set_yscale("log")
    # ax2.legend()
    # ax2.set_xlabel("Number of rows in dataset")
    # ax2.set_ylabel("Number of CFDs found")
    # plt.show()
    
    #! columns
    # columns = [2,3,4,5,6,7]
    # CFDTime = []
    # ACFDTime = []
    # CFDCount = []
    # ACFDCount = []
    # # baseline: k=5,  eps=0.1, rows=100, columns = 5
    # for column in columns:
    #     print(column)
    #     data2D = pd.read_csv(r'D:\Western\CS9866\FD_CFD_extraction\adult.csv') #15 columns, 32561 rows
    #     data2D = data2D.head(100) # take first x rows
    #     data2D = data2D.iloc[:, : column-1] # take first column columns
        
    #     print("#########################################################") # normal CFD
    #     eps = 0
    #     k = 5
    #     numCFD, timeTaken = getCFD(data2D, eps, k)
    #     CFDTime += [timeTaken]
    #     CFDCount += [numCFD]
        
    #     print("*********************************************************") # Approx CFD
    #     eps = 0.1
    #     k = 5
    #     numCFD, timeTaken = getCFD(data2D, eps, k)
    #     ACFDTime += [timeTaken]
    #     ACFDCount += [numCFD]
        
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    # fig.suptitle('CFD vs ACFD Discovery - Columns')
    # ax1.scatter(columns, CFDTime, c="blue", label="CFD")
    # ax1.scatter(columns, ACFDTime, c="red", label="ACFD")
    # ax1.plot(columns, CFDTime, 'b-', alpha=0.5)  # Blue line for CFD
    # ax1.plot(columns, ACFDTime, 'r-', alpha=0.5)  # Red line for ACFD
    # # ax1.set_xscale("log")
    # ax1.set_yscale("log")
    # ax1.legend()
    # ax1.set_xlabel("Number of columns in dataset")
    # ax1.set_ylabel("Time taken")
    
    # ax2.scatter(columns, CFDCount, c="blue", label="CFD")
    # ax2.scatter(columns, ACFDCount, c="red", label="ACFD")
    # ax2.plot(columns, CFDCount, 'b-', alpha=0.5)  # Blue line for CFD
    # ax2.plot(columns, ACFDCount, 'r-', alpha=0.5)  # Red line for ACF
    # # ax2.set_xscale("log")
    # # ax2.set_yscale("log")
    # ax2.legend()
    # ax2.set_xlabel("Number of columns in dataset")
    # ax2.set_ylabel("Number of CFDs found")
    # plt.show()
    
    #! k
    # ks = [2,3,4,5,6,7,8,9,10]
    # CFDTime = []
    # ACFDTime = []
    # CFDCount = []
    # ACFDCount = []
    # # baseline: k=5,  eps=0.1, rows=100, columns = 5
    # for k in ks:
    #     print(k)
    #     data2D = pd.read_csv(r'D:\Western\CS9866\FD_CFD_extraction\adult.csv') #15 columns, 32561 rows
    #     data2D = data2D.head(100) # take first x rows
    #     data2D = data2D.iloc[:, : 4] # take first column columns
        
    #     print("#########################################################") # normal CFD
    #     eps = 0
    #     numCFD, timeTaken = getCFD(data2D, eps, k)
    #     CFDTime += [timeTaken]
    #     CFDCount += [numCFD]
        
    #     print("*********************************************************") # Approx CFD
    #     eps = 0.1
    #     numCFD, timeTaken = getCFD(data2D, eps, k)
    #     ACFDTime += [timeTaken]
    #     ACFDCount += [numCFD]
        
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    # fig.suptitle('CFD vs ACFD Discovery - Support k')
    # ax1.scatter(ks, CFDTime, c="blue", label="CFD")
    # ax1.scatter(ks, ACFDTime, c="red", label="ACFD")
    # ax1.plot(ks, CFDTime, 'b-', alpha=0.5)  # Blue line for CFD
    # ax1.plot(ks, ACFDTime, 'r-', alpha=0.5)  # Red line for ACFD
    # # ax1.set_xscale("log")
    # ax1.set_yscale("log")
    # ax1.legend()
    # ax1.set_xlabel("Support threshold k")
    # ax1.set_ylabel("Time taken")
    
    # ax2.scatter(ks, CFDCount, c="blue", label="CFD")
    # ax2.scatter(ks, ACFDCount, c="red", label="ACFD")
    # ax2.plot(ks, CFDCount, 'b-', alpha=0.5)  # Blue line for CFD
    # ax2.plot(ks, ACFDCount, 'r-', alpha=0.5)  # Red line for ACF
    # # ax2.set_xscale("log")
    # # ax2.set_yscale("log")
    # ax2.legend()
    # ax2.set_xlabel("Support threshold k")
    # ax2.set_ylabel("Number of CFDs found")
    # plt.show()
    
    #! eps
    
    # baseline: k=5,  eps=0.1, rows=100, columns = 5
    data2D = pd.read_csv(r'D:\Western\CS9866\FD_CFD_extraction\adult.csv') #15 columns, 32561 rows
    data2D = data2D.head(100) # take first x rows
    data2D = data2D.iloc[:, : 4] # take first column columns
    
    print("#########################################################") # normal CFD
    k=5
    eps = 0
    numCFD, timeTaken = getCFD(data2D, eps, k)
    
    print("*********************************************************") # Approx CFD
    k = 5
    eps = 0.9
    print(eps)
    numCFD, timeTaken = getCFD(data2D, eps, k)
    
    epss = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    CFDTime = [0.146,0.146,0.146,0.146,0.146,0.146,0.146,0.146,0.146]
    ACFDTime = [0.15204,0.143,0.147,0.14503,0.146,0.14197,0.14601,0.142,0.143]
    CFDCount = [19,19,19,19,19,19,19,19,19]
    ACFDCount = [24,25,33,34,34,34,45,45,59]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle('CFD vs ACFD Discovery - Epsilon')
    ax1.scatter(epss, CFDTime, c="blue", label="CFD")
    ax1.scatter(epss, ACFDTime, c="red", label="ACFD")
    ax1.plot(epss, CFDTime, 'b-', alpha=0.5)  # Blue line for CFD
    ax1.plot(epss, ACFDTime, 'r-', alpha=0.5)  # Red line for ACFD
    # ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend()
    ax1.set_xlabel("Error tolerance psilon")
    ax1.set_ylabel("Time taken")
    
    ax2.scatter(epss, CFDCount, c="blue", label="CFD")
    ax2.scatter(epss, ACFDCount, c="red", label="ACFD")
    ax2.plot(epss, CFDCount, 'b-', alpha=0.5)  # Blue line for CFD
    ax2.plot(epss, ACFDCount, 'r-', alpha=0.5)  # Red line for ACF
    # ax2.set_xscale("log")
    # ax2.set_yscale("log")
    ax2.legend()
    ax2.set_xlabel("Error tolerance psilon")
    ax2.set_ylabel("Number of CFDs found")
    plt.show()

#! keep track of time and CFDs found
# increase k
# increase eps
# increase rows
# increase columns