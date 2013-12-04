#Local Alignment Python Code
def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in xrange(sizex)]


class ScoreParam:
    """The parameters for an alignment scoring function"""
    def __init__(self, gap, match, mismatch):
        self.gap = gap
        self.match = match
        self.mismatch = mismatch
        
def getalignscore(x,i,y,j,score):
    # a function that safely returns a match or mismatch score if
    # all indexes accessed are safe otherwises returns a number that can not be max
    #in local alignment
    cantalign=-1000000
    if ((i>-1) and (j>-1)):
        if (i<len(x) and j<len(y)):
            return (score.match if x[i] == y[j] else score.mismatch)
        else:
            return cantalign
    else:
        return cantalign


def properalign(x,y):
    #function that properly formats 2 strings to be locally aligned
    #using rpvided algorithm
    # adjusts inputs for the fact first pair is not scored and
    # strings should be of the same length
    if len(x)<len(y):
        x=x+ "X" *(len(y)-len(x))
        #print "print1", (len(y)-len(x)), len(y),len(x)
    if len(x)>len(y):
        y=y+ "Y" *(len(x)-len(y))
        #print "print2", (len(x)-len(y)), len(x),len(y)
    #print "1"+x
    #print "2"+y
    return local_align("1"+x, "2"+y)
        
def midfastproperalign(x,y):
    # makes sures strings are the same length
    # and uses fast aligner that only returns score
    if len(x)<len(y):
        x=x+ "X" *(len(y)-len(x))
        #print "print1", (len(y)-len(x)), len(y),len(x)
    if len(x)>len(y):
        y=y+ "Y" *(len(x)-len(y))
        #print "print2", (len(x)-len(y)), len(x),len(y)
    #print "1"+x
    #print "2"+y
    return fast_local_align("1"+x, "2"+y)

def fastproperalign(x,y):
    # assumes strings are the same length
    # and uses fast aligner that only returns score
    #assert(len(x)==len(y))
    #print "1"+x
    #print "2"+y
    return fast_local_align("1"+x, "2"+y)

def fast_local_align(x, y, score=ScoreParam(-1, 1, -1)):
    """Do a local alignment between x and y and quickly and efficently as possible"""
    # create a zero-filled matrix
    A = make_matrix(len(x) + 1,len(y) + 1)
    
    best = 0
    optloc = (0,0)
    # fill in A in the right order
    for i in xrange(1, len(y)):
        for j in xrange (1, len(x)):
            # the local alignment recurrance rule:
            A[i][j] = max(A[i][j-1] + score.gap,
                        A[i-1][j] + score.gap,
                        A[i-1][j-1] + (score.match if x[i] == y[j] else score.mismatch),
                        0)
            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j)
    #for i in A:print i
    
    # returns only the opt score
    return best

def local_align(x, y, score=ScoreParam(-1, 1, -1)):
    """Do a local alignment between x and y"""
    # create a zero-filled matrix for both scores and pointers for alignment
    A = make_matrix(len(x) + 1,len(y) + 1)
    B = make_matrix(len(x) + 1,len(y) + 1)
    
    best = 0
    optloc = (0,0)
    # fill in A in the right order
    for i in xrange(1, len(y)):
        for j in xrange (1, len(x)):
            # the local alignment recurrance rule:
            A[i][j] = max(A[i][j-1] + score.gap,
                        A[i-1][j] + score.gap,
                        A[i-1][j-1] + getalignscore(x,i,y,j,score),
                        0)
            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j)
            # fills in pointers in matrix to point to optimal path leading up to current
            # alignment

            if(A[i][j]==(A[i][j-1] + score.gap)):B[i][j]=(i,j-1)
            elif(A[i][j]==(A[i-1][j] + score.gap)):B[i][j]=(i-1,j)
            elif(A[i][j]==(A[i-1][j-1] + getalignscore(x,i,y,j,score))):B[i][j]=(i-1,j-1)
            #if(A[i][j-1]==0 and A[i-1][j]==0 and A[i-1][j-1]==0):B[i][j]=="stop"

    #for i in A:print i
    #for i in B:print i

    # return the opt score and the best location snd score matrix and pointer matrix
    return best, optloc, B, A

def makescorelist(best, optloc, B, A):
    #converts score and pointer matrices into
    # 2 lists that hold best alignment scores by position in alignemnt
    # and best position alignments in local alignment
    #lists that hold optimal alignement scores
    scores=[]
    positions=[]
    # start back trace at optimal score position in matrix
    position=optloc
    positions.append(optloc)
    assert(A[position[0]][position[1]]==best)
    while(position!=(0,0)):
        #print position[0], position[1]
        #print A[position[0]]
        #print A[position[0]][position[1]]
        #print position
        # record score at that position
        scores.append(A[position[0]][position[1]])
        #update the current position to point to the next earlier position
        #if (B[position[0]][position[1]]=="stop"):break
        position=B[position[0]][position[1]]
        #this prevents errors for entries that did not get tuple position
        if type(position)!=type((0,0)):
            position=(0,0)
        #appends a zero if the current position is (0,0)
        if (position == (0,0)):scores.append(0)
        #append earlier position to list f scores
        positions.append(position)
    #print scores
    #print positions
    # reverse lists so each list starts at beginning of opt alignment
    positions.reverse()
    scores.reverse()
    #print scores
    #print positions
    return scores,positions

def getshift(positions,position,seqlen):
    # record current postion in query sequences in current alignment position
    oldchar=positions[position][1]
    #initilize variable that will say how many positions in alignemnt
    #must be included so that "seqlen" positions in query are used
    # for next extend alignemnt step
    stop="no"
    for i in xrange(len(positions[position:])):
        #print i, positions[position+i]
        # for each position in rest of positions in positions in alignemnt
        # checks if the current position in the query is "seqlen" more than old position
        # in string, if it is record position and exit loop 
        if (positions[position+i][1]-oldchar==seqlen):
            stop=i
            #print "stopped"
            break
    # if type of stop is not int, means that there are not enough chars left in query go so
    #that "seqlen" more chars can be selected from query in current alignment
    if (type(stop)!=type(i)):stop= len(positions[position:])-1
    return stop

def findreasonablepositions(scores,positions,cutoff,seqlen):
    #the index of the point up until the match meets the cutoff

    
    #initilize position in list
    position=0
    #finds how many positions later in alignemnt will query position will change by
    #seqlen
    shift=getshift(positions,position,seqlen)
    #print "shift", shift
    # while we dont past end of alingment
    while(position<len(scores)-shift):
        #find score of current position in local alignment
        oldscore=scores[position]
        #finds how many positions later in alignemnt will query position will change by
        #seqlen        
        shift=getshift(positions,position,seqlen)
        #print "shift", shift
        # change position by shift
        position+=shift
        #print "position", position
        #get score of position
        newscore=scores[position]
        # find how much score has changed over these "seqlen" positions in query
        delscore=newscore-oldscore
        #print delscore
        # if this change in score is not good enough likely no longer properly aligned
        #to correct transcript
        if (delscore<cutoff):
            #print "not good now"
            # undo last shift, since only up to previous position change was alignment
            #good and save that position
            dposition=position-shift
            # end loop
            break
        # if we reach end of alignment store that position as end of current aligning
        # to be end of current matching between transcripts
        elif(not(position<len(scores)-shift)):
            #print "end"
            dposition=position
    return dposition
x=properalign("ABCDEFXXGHI","ABCDKFGHI")
y=makescorelist(*x)
print y
print x
a=findreasonablepositions(y[0],y[1],0,2)
#print a
# examples for testing
#x=properalign("ABCDEFXXGHI","ABCDKFGHI")
#y=makescorelist(*x)
#a=findresonablepositions(y[0],y[1],0,2)
