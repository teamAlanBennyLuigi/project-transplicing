"""check sequences"""

#the bug : FBtr0070099[23-285]  S91N6:968:518[80-80]
#29,85

#open files with raw transcripts and ref transcripts
seqfile=open("test1.fasta")
reffile=open("sequence6992.fasta")
#import usegul modules
#import math
import time
from kingsfordlocalaligner import *

#make dictionary for raw transcripts
rawdic=dict()
#initilize variable for sequence data
seq=""
#go through file and assemble sequence data and sequence labels
for line in seqfile:
    #if line is a new name
    if (line[0]==">"):
        # and sequence is not empty
        if(seq!=""):
            #if sequence in dictionary already update the list of labels
            # current sequence has
            if (seq in rawdic):
                # if only one label present make a list of labels
                if type(rawdic[seq][0])!=type([]):newname=[rawdic[seq][0]]+[name[:-1]]
                #otherwise just add new entry to list of labels
                else: newname=rawdic[seq][0]+[name[:-1]]
                #store list of labels
                name=newname
            # store entry of previous sequence data in dictionary
            rawdic[seq]=(name[:-1],[[]])
            # start new sequence to record
            seq=""
        #get label of new sequence
        name=line
    elif (line[0]=="A" or line[0]=="G" or line[0]=="C" or line[0]=="T"):
        #if sequence data add to previous sequence data
        seq=seq+line[:-1]
    elif (line[0]=="\n"):pass
    else: print "something is not right:", line

#make dictionary for raw transcripts
refdic=dict()
#initilize variable for sequence data
seq=""
#counter for number of genes
counter1=0
#go through file and assemble sequence data and sequence labels
for line in reffile:
    #if line is a new name
    if (line[0]==">"):
        # and sequence is not empty
        counter1+=1
        if(seq!=""):
            #print "here1"
            # extract data like name and chromosome
            name=title.split(" ")[0][1:]
            data=title.split(" ")[1].split(":")
            chromo=data[0]
            position=data[1]
            position=position.split("-")
            start=int(position[0])
            length=int(position[1])-start+1
            #print "here1a"
            #use dot upper so that both raw and ref seq are
            #useing the same type of chars
            refdic[seq.upper()]=[name,chromo,start,length]
            #print "here1b"
            #print refdic[seq]
            seq=""
        #get name of new sequence
        title=line
        #print "here2"
    elif (line[0]=="a" or line[0]=="g" or line[0]=="c" or line[0]=="t"):
        #if sequence data add to previous sequence data
        seq=seq+line[:-1]
        #print "here3"
    elif (line[0]=="\n"):pass
    else: print "something is not right:", line

#old testing code
""">>> count =0
>>> for i in refdic:
	print i
	print refdic[i]
	count+=1
	if count ==10:break"""

seqfile.close()
reffile.close()

def alignscore(query,reference):
    #old function to aligning not local not that good
    assert(len(query)<=len(reference))
    score=0
    match=1
    mismatch=0
    for i in xrange(len(query)):
        if (query[i]==reference[i]):score+=match
        else:score+=mismatch
    return score

def simplesearch(x,refdic):
    # brute force search for best alignment , can take hours for 1 transcript
    bestscore=0
    best="non"
    counter=0
    for ref in refdic:
        print counter
        counter+=1
        test=midfastproperalign(x,ref)
        if(test>bestscore):
            bestscore=test
            best=ref
            print bestscore, best
    print bestscore, best

def mode(valueList):
    #got from http://stackoverflow.com/questions/9567687/calculating-the-mode-in-a-multimodal-list-in-python
    # picks the lowest value mode
    frequencies = {}
    mx = None
    for value in valueList:
        if value in frequencies:
            frequencies[value] += 1
        else:
            frequencies[value] = 1
        if not mx or frequencies[value] > mx[1]:
            mx = (value, frequencies[value])

    mode = mx[0]
    return mode

def fastersearchclosealignment(query,refdic):
    #will hold best alignment
    best=()
    #holds "fingerprints", somewhat similar to "seeds" in BLAST
    #which speeds up search, for small pieces that are perfect matches
    fingerprints=[]
    #- len of fingerprint to prevent error
    for i in xrange(0,len(query)-10,10):fingerprints.append(query[i:i+20])
    #print fingerprints[50]
    for ref in refdic:
        #if (ref=="AAAAAG"):
        #    print ref
        # holds ref seq
        reftran=ref
        #holds bestscore over each reference sequence
        bestScore=0
        #holds bestscore over a match between a fingerprint and reference sequence
        bestscore=0
        #position where ref sequence starts matching to raw transcript 
        refposition=0
        #list of best positions that matched to fingerprints
        bestpositions=[]
        #holds current best location
        bestposition=0
        # holds score of best alignment
        score=0
        seq="non"
        # if exact match is found to whole string then that is the best score
        if query in reftran:
            #print "            match"
            bestposition= reftran.index(query)
            #print bestposition, "test"
            bestscore=len(query)
            bestpositions.append(bestposition)
            refposition=bestposition
            seq=reftran
        countprints=0
        #otherwise search for if multiple fingerprint exactly match to
        #the reference sequence, if so then  a match 
        for prin in xrange(len(fingerprints)):
            if fingerprints[prin] in reftran:
                #print "            match with fingerprint"
                positioninquery=prin*10
                bestposition= reftran.index(fingerprints[prin])

                #print bestposition, positioninquery,bestposition-positioninquery
                bestposition=bestposition-positioninquery
                bestpositions.append(bestposition)
                bestscore=len(fingerprints[prin])
                refposition=bestposition
                #print bestposition, "test"
                seq=reftran
                countprints+=1
                # if this reference sequence scroes highly in a local alignment
                # then it might be the match
                if (midfastproperalign(query[:30],reftran[bestposition:bestposition+40])>20):
                    #print "              somewhat good"
                    #print query[:50]
                    #print reftran[bestposition:bestposition+50], bestposition
                    #print midfastproperalign(query[:50],reftran[bestposition:bestposition+50])
                    #print reftran
                    score=midfastproperalign(query[:50],reftran[bestposition:bestposition+50])
                    #print bestscore, score
                    if (bestscore<score):
                        bestScore=score
                        #print bestscore, score, "good"
                        #print bestposition, "test"
                        refposition=bestposition
                        seq=reftran
        if(len(bestpositions)>0):
            # the best position is found by majority because there are false positives
            bestposition=mode(bestpositions)
            #print "HEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHERE"

        if(best==()):
            #print bestscore            
            best=(bestscore,refposition,seq,50)
            #print best
        elif (best[0]<bestScore):
            #if the current mapping to a reference sequence is better than th older best
            #one update the annotation
            #print bestpositions
            #print """




            #bestposition:""", bestposition
            #print best[0] ,bestScore
            best=(bestScore,refposition,seq,50)
            #print best

        """for ref in refdic:
        for i in xrange(len(reftran)-len(query)+1):
            if (i%500==0):print i
            if (fastproperalign(query[:10],reftran[i:i+10])>6):
                print "              somewhat good"
                score=fastproperalign(query,reftran[i:i+len(query)])
            if (bestscore<score):
                bestscore=score
                refposition=i
                seq=reftran
        """
        """if(countprints>9):
            #print "            atleast 9 matches with fingerprint"
            #print best"""
    return best

def searchclosealignment(query,refdic):
    #older version not useable
    best=()
    for ref in refdic:
        #if (ref=="AAAAAG"):
        #    print ref
        reftran=ref
        bestscore=0
        refposition=0
        score=0
        seq="non"
        if query in reftran:
            bestposition= reftran.index(query)
            bestscore=len(query)
            refposition=bestposition
            seq=reftran
        else:
            for i in xrange(len(reftran)-len(query)+1):
                if (i%500==0):print i
                if (fastproperalign(query[:10],reftran[i:i+10])>6):
                    print "              somewhat good"
                    score=fastproperalign(query,reftran[i:i+len(query)])
                if (bestscore<score):
                    bestscore=score
                    refposition=i
                    seq=reftran
        """
        old for reference
        else:
            for i in xrange(len(reftran)-len(query)+1):
                if (alignscore(query[:10],reftran[i:i+10])>8):
                    score=alignscore(query,reftran[i:i+len(query)])
                if (bestscore<score):
                    bestscore=score
                    refposition=i
                    seq=reftran"""
        if(best==()):best=(bestscore,refposition,seq,len(query))
        elif (best[0]<bestscore):best=(bestscore,refposition,seq,len(query))
    return best
            

"""refdic={}
base=["A","G","C","T"]
len(refdic)
for i in xrange(4**6):
	string=""
	if ((i%2==0) or (i%3==0)):continue
	for j in xrange(6):
		string=string+base[i%4]
		i=i/4
	refdic[string]="oh"


"""
def extend10old(transcript,rawposition,ref,refposition):
    #old function do not use
    tran=transcript[rawposition+1:]
    re=ref[refposition+1:]

    return tran[:10],re[:10]


# there are glitches that need to be fixed
trancounter=0
a=time.time()
for transcript in rawdic:
    print "1 Transcript done"
    query=transcript
    #print query, len(query)
    #hit=searchclosealignment(query,refdic)
    b=time.time()
    startinraw=0
    endinraw=0
    stoping=False
    ##########

    while(endinraw<len(query)):
        c=time.time()
        hit2=fastersearchclosealignment(query[startinraw:startinraw+50],refdic)
        #print hit2
        if (hit2[0]>40):
            print hit2[0]
            #print startinraw, len(query[startinraw:]),len(query[startinraw:startinraw+50])
            #print query[startinraw:]
            #print query[startinraw:startinraw+50],  "I'm here!!!!"
            #print
            #print hit2[2][hit2[1]:hit2[1]+50]
            ref=hit2[2]
            refposition=hit2[1]
            #print "position is : ", hit2[1]
            #print ref[refposition:refposition+2*len(query[startinraw:])],"ref"
            print  query[startinraw:],len(query[startinraw:])
            while(True):
                if(not fastproperalign(query[startinraw:startinraw+10],
                                       ref[refposition:refposition+10])>6):
                    startinraw+=10
                    print startinraw
                    refdata=refdic[ref]
                    stuf=fastersearchclosealignment(query[startinraw:],{ref:data})
                    newhit=stuf[1]
                    refposition=newhit
                    print  query[startinraw:],len(query[startinraw:])
                    stoping=True
                else:break
            # careful change 150 to variable
            #cutoff is necessary to prevent overshadowing of later but much larger and longer hits
            data=properalign(query[startinraw:startinraw+150],
                             ref[refposition:refposition+len(query[startinraw:])+5])
            scroingdata=makescorelist(*data)
            end=findreasonablepositions(scroingdata[0],scroingdata[1],6,10)
            #print end
            findendinref=properalign(query[startinraw:end],
                                     ref[refposition:refposition+2*len(query[startinraw:])])
            name= refdic[ref][0]
            start=refposition
            stop=findendinref[1][0]+refposition
            #startinraw=0
            endinraw=end+startinraw
            chromo= refdic[ref][1]
            #if(start<0 or stop<0):break
            rawdic[query][1][-1]=[name,chromo,
                                  str(startinraw)+"-"+str(endinraw),
                                  str(start)+"-"+str(stop)]
            rawdic[query][1].append([])
            print rawdic[query]
            d=time.time()
            #print b-a, "since start"
            #print d-c, "this loop"
            startinraw=endinraw
            if stoping:break
            if ((endinraw)>(len(query)-30)):
                #print rawdic[query]
                break
        else:
            #print startinraw, len(query[startinraw:]),len(query[startinraw:startinraw+50])
            #print query[startinraw:]
            #print query[startinraw:startinraw+50],  "I'm here!!!!"
            #print
            #print hit2[2][hit2[1]:hit2[1]+50]
            #print "position is : ", hit2[1]
            #print ref[refposition:refposition+2*len(query[startinraw:])]
            #print end
            name= "non"
            endinraw=startinraw+50
            chromo= "non"
            rawdic[query][1][-1]=[name,chromo,
                                  str(startinraw)+"-"+str(endinraw),
                                  "non"]
            rawdic[query][1].append([])
            print rawdic[query]
            d=time.time()
            #print d-a, "since start"
            #print d-c, "this loop"""
            startinraw=endinraw
            if stoping:break
            if ((endinraw)>(len(query)-50)):
                #print rawdic[query]
                break
        #contin=extend10(transcript,rawposition,ref,refposition)
        #print contin[0]
        #print contin[1]

    trancounter+=1
    if stoping:break
    if(trancounter>=100):break
"""
>>> test1="CATTACAGTGCCATCTTCAGAAGACAAAAAAGGATTTGCACCAGAGGACCAGGGATCAGGAGCAAAGAAGCAACAGCAAACCATGGCAGTGGAAGTAGTTCAGGAGACGCTGCAACAAGCGGCGTCCAGTA"
>>> test2="ATAGCAACAGTGCCATCTTCAGAAGACAAAAAGGATTTGCACCAGAGGACCAGGGATCAGGAGCAAAGAAGCAACAGCAACCATGGCAGTGGAAGTAGTTCAGGAGACGCTGCAACAAGCGGCGTCCAGTTCGTCGACGACGGTCCTGGGATTCAGTCCTATGTTAACCACCTTAGTGGGCACCCTGGTGGCCATGGCATTGTACGAGTATTGGCGCAGGAATAGCCGGGAATACCGCATGGTTGCCAATATACCATCCCCA"
>>> y=properalign(test1,test2)
>>> y[0]
122

error at
'>S91N6:968:518', [['FBtr0089175', '4', '0-80', '805-885'], ['FBtr0070099', 'X', '90-211', '33-64']
""" 
