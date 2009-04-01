### General Alignment Tools, Erez Lieberman.

### Includes global (Needleman-Wunsch) and local (Smith-Waterman) alignment, following the approach of Durbin, Biological Sequence Analysis.
### This has both linear and affine gap costs. The four critical parameters are listed at the top.

from numarray import array,where,resize,zeros,Float32

matchscore=.5
mismatchscore=-.7
gappenalty=-2.
#gapextend=-.4
gapextend=-1.5

#################################################
def seqs2scorematrix(seq1,seq2,reset=0):

    matchscore=.5
    mismatchscore=-.75
    gappenalty=-2.
#    gapextend=-.4
    gapextend=-1.5

    scorematrix=zeros((len(seq1)+1,len(seq2)+1),type=Float32)
    backtrack=zeros((len(seq1)+1,len(seq2)+1))    

    default=[]
    if reset==1:
        default=[0]

    for row in xrange(len(seq1)+1):
        scorematrix[row][0]=max([gapextend*row]+default)
        backtrack[row][0]=1
    
    for col in xrange(len(seq2)+1):
        scorematrix[0][col]=max([gapextend*col]+default)
        backtrack[0][col]=0        
    
    for row in xrange(1,len(seq1)+1):              
        for col in xrange(1,len(seq2)+1):
            if seq1[row-1]==seq2[col-1]:
                currmatchscore=matchscore
            else:
                currmatchscore=mismatchscore

            if scorematrix[row-1][col]<=scorematrix[row][col-1]:
                bonus1=0
                guy1=scorematrix[row][col-1]+gapextend
                
            else:
                bonus1=1
                guy1=scorematrix[row-1][col]+gapextend
                
            if reset==1 and scorematrix[row-1][col-1]+currmatchscore<0:
                guy2=0
                bonus2=1
                
            else:
                guy2=scorematrix[row-1][col-1]+currmatchscore
                bonus2=0

            if guy1>=guy2:
                scorematrix[row][col]=guy1
                backtrack[row][col]=bonus1
            else:
                scorematrix[row][col]=guy2
                backtrack[row][col]=2+bonus2


    return scorematrix,backtrack

#################################################
def seqs2affinematrix(seq1,seq2,reset=0):

    matchscore=.5
#    mismatchscore=-.7
    mismatchscore=-.75
    gappenalty=-2.
    gapextend=-.4
##    import time
##
##
##    t0=time.time()
    scorematrix=zeros((len(seq1)+1,len(seq2)+1),type=Float32)
    gapmatrix=zeros((len(seq1)+1,len(seq2)+1),type=Float32)
    backtrack=zeros((len(seq1)+1,len(seq2)+1))    
    gapbacktrack=zeros((len(seq1)+1,len(seq2)+1))    

    default=[]
    if reset==1:
        default=[0]

##    t1=time.time()


    
    for row in xrange(len(seq1)+1):
        #scorematrix gets +1 so that opening gap is always tracked back to gapmatrix        
        scorematrix[row][0]=max([gappenalty+gapextend*row]+default)
        gapmatrix[row][0]=max([gappenalty+gapextend*(row-1)]+default)        
        backtrack[row][0]=1
        gapbacktrack[row][0]=3

##    t2=time.time()
    
    for col in xrange(len(seq2)+1):
        #scorematrix gets +1 so that opening gap is always tracked back to gapmatrix
        scorematrix[0][col]=max([gappenalty+gapextend*col]+default)
        gapmatrix[0][col]=max([gappenalty+gapextend*(col-1)]+default)        
        backtrack[0][col]=0
        gapbacktrack[0][col]=2

    #prevents added penalty for starting new alignment at beginning....I think!
    scorematrix[0][0]=0
    gapmatrix[0][0]=0    
    
##
##    t3=time.time()        
##

    for row in xrange(1,len(seq1)+1):
        for col in xrange(1,len(seq2)+1):
            if seq1[row-1]==seq2[col-1]:
                currmatchscore=matchscore
            else:
                currmatchscore=mismatchscore
            #currmatchscore=(seq1[row-1]==seq2[col-1])*(matchscore-mismatchscore)+mismatchscore
            if scorematrix[row-1][col-1]>=gapmatrix[row-1][col-1]:
                backtrack[row][col]=0
                newscore=max(default+[scorematrix[row-1][col-1]+currmatchscore])
            else:
                backtrack[row][col]=1
                newscore=max(default+[gapmatrix[row-1][col-1]+currmatchscore])

            scorematrix[row][col]=newscore
        
            if scorematrix[row][col-1]>=scorematrix[row-1][col]:
                guy1=scorematrix[row][col-1]+gappenalty
                bonus1=0
            else:
                guy1=scorematrix[row-1][col]+gappenalty
                bonus1=1
            
            if gapmatrix[row][col-1]>=gapmatrix[row-1][col]:
                guy2=gapmatrix[row][col-1]+gapextend
                bonus2=0                
            else:
                guy2=gapmatrix[row-1][col]+gapextend
                bonus2=1
                
            if guy1>=guy2:
                gapmatrix[row][col]=max(default+[guy1])
                gapbacktrack[row][col]=bonus1
            else:
                gapmatrix[row][col]=max(default+[guy2])
                gapbacktrack[row][col]=2+bonus2

            
                

##    t4=time.time()

##    print '--'
##    print t1-t0
##    print t2-t1
##    print t3-t2
##    print t4-t3
##    print t4-t0
##    print '--'
##

##    print scorematrix

    return scorematrix,backtrack,gapmatrix,gapbacktrack


    
#################################################    
def dmatrices2global(seq1,seq2,scorematrix,backtrack):

    aln1=''
    aln2=''

    row=len(seq1)
    col=len(seq2)
    
    while row>0 or col>0:

        if backtrack[row][col]==0:
            aln1+='-'
            aln2+=seq2[col-1]
            newcol=col-1
            newrow=row

        elif backtrack[row][col]==1:
            aln1+=seq1[row-1]
            aln2+='-'
            newcol=col
            newrow=row-1

        elif backtrack[row][col]==2:
            aln1+=seq1[row-1]
            aln2+=seq2[col-1]
            newrow=row-1
            newcol=col-1

        row=newrow
        col=newcol

    return [[aln1[::-1],aln2[::-1]],scorematrix[len(seq1)][len(seq2)]]

        
#################################################
def dmatrices2local(seq1,seq2,scorematrix,backtrack,ignorezeros=0):


    topscore=scorematrix.max()

    if ignorezeros==1 and topscore<0:
        alignments=[[['',''],0,[[0,0],[0,0]]]]
        return alignments
    
    maxlocation=where(scorematrix==topscore)

    alignments=[]
    for i in range(len(maxlocation[0])):
        aln1=''
        aln2=''

        row=maxlocation[0][i]
        col=maxlocation[1][i]

        while (row>0 or col>0) and (ignorezeros==1 or scorematrix[row][col]!=0):
            if backtrack[row][col]==0:
                aln1+='-'
                aln2+=seq2[col-1]
                newcol=col-1
                newrow=row

            elif backtrack[row][col]==1:
                aln1+=seq1[row-1]
                aln2+='-'
                newcol=col
                newrow=row-1
                
            elif backtrack[row][col]==2:
                aln1+=seq1[row-1]
                aln2+=seq2[col-1]
                newrow=row-1
                newcol=col-1

            row=newrow
            col=newcol

        alignments.append([[aln1[::-1],aln2[::-1]],topscore,[[row,maxlocation[0][i]],[col,maxlocation[1][i]]]])

    return alignments


#################################################    
def affmatrices2global(seq1,seq2,scorematrix,backtrack,gapmatrix,gapbacktrack):

    aln1=''
    aln2=''

    row=len(seq1)
    col=len(seq2)

    if scorematrix[row][col]>=gapmatrix[row][col]:
        gapped=0
    else:
        gapped=1
    
    while row>0 or col>0:

        if backtrack[row][col]==0 and gapped==0:
            aln1+=seq1[row-1]
            aln2+=seq2[col-1]
            newrow=row-1
            newcol=col-1
            newgapped=gapped            

        if backtrack[row][col]==1 and gapped==0:
            aln1+=seq1[row-1]
            aln2+=seq2[col-1]
            newrow=row-1
            newcol=col-1
            newgapped=1            

        if gapbacktrack[row][col]==0 and gapped==1:
            aln1+='-'
            aln2+=seq2[col-1]
            newcol=col-1
            newrow=row
            newgapped=0

        if gapbacktrack[row][col]==2 and gapped==1:
            aln1+='-'
            aln2+=seq2[col-1]
            newcol=col-1
            newrow=row
            newgapped=gapped

        if gapbacktrack[row][col]==1 and gapped==1:
            aln1+=seq1[row-1]
            aln2+='-'
            newcol=col
            newrow=row-1
            newgapped=0

        if gapbacktrack[row][col]==3 and gapped==1:
            aln1+=seq1[row-1]
            aln2+='-'
            newcol=col
            newrow=row-1
            newgapped=gapped
        
        row=newrow
        col=newcol
        gapped=newgapped

    return [[aln1[::-1],aln2[::-1]],scorematrix[len(seq1)][len(seq2)]]

        
#################################################
def affmatrices2local(seq1,seq2,scorematrix,backtrack,gapmatrix,gapbacktrack,ignorezeros=0):
##
##    print scorematrix
##    print backtrack
##    print gapmatrix
##    print gapbacktrack

    
    topscore=scorematrix.max()
    
    if ignorezeros==1 and topscore<0:
        alignments=[[['',''],0,[[0,0],[0,0]]]]
        return alignments
    
    maxlocation=where(scorematrix==topscore)

##    print topscore
##    print maxlocation


    alignments=[]
    for i in range(len(maxlocation[0])):
        aln1=''
        aln2=''

        row=maxlocation[0][i]
        col=maxlocation[1][i]
        gapped=0


#not (scorematrix[row][col]!=0 and gapped==0) and not (gapmatrix[row][col]!=0 and gapped==1)
        #while (row>0 or col>0) and (ignorezeros==1 or scorematrix[row][col]!=0):            
        while (row>0 or col>0) and (ignorezeros==1 or (not (scorematrix[row][col]==0 and gapped==0) and not (gapmatrix[row][col]==0 and gapped==1))):
            
            if backtrack[row][col]==0 and gapped==0:
                aln1+=seq1[row-1]
                aln2+=seq2[col-1]
                newrow=row-1
                newcol=col-1
                newgapped=gapped            

            if backtrack[row][col]==1 and gapped==0:
                aln1+=seq1[row-1]
                aln2+=seq2[col-1]
                newrow=row-1
                newcol=col-1
                newgapped=1

            if gapbacktrack[row][col]==0 and gapped==1:
                aln1+='-'
                aln2+=seq2[col-1]
                newcol=col-1
                newrow=row
                newgapped=0

            if gapbacktrack[row][col]==2 and gapped==1:
                aln1+='-'
                aln2+=seq2[col-1]
                newcol=col-1
                newrow=row
                newgapped=gapped

            if gapbacktrack[row][col]==1 and gapped==1:
                aln1+=seq1[row-1]
                aln2+='-'
                newcol=col
                newrow=row-1
                newgapped=0

            if gapbacktrack[row][col]==3 and gapped==1:
                aln1+=seq1[row-1]
                aln2+='-'
                newcol=col
                newrow=row-1
                newgapped=gapped


##
##            print row,col,gapped
##            print scorematrix[row][col],gapmatrix[row][col]
##            print backtrack[row][col],gapbacktrack[row][col]
##            print aln1[::-1],aln2[::-1]            
##            print newrow,newcol,newgapped
##            print gapmatrix[row-1][col],gapmatrix[row][col-1],gapmatrix[row-1][col-1]            
##            print scorematrix[row-1][col],scorematrix[row][col-1],scorematrix[row-1][col-1]            
##            print '--'

            
            row=newrow
            col=newcol
            gapped=newgapped


            

        alignments.append([[aln1[::-1],aln2[::-1]],topscore,[[row,maxlocation[0][i]],[col,maxlocation[1][i]]]])


    


    return alignments





#################################################
def align(seq1,seq2,local0orglobal1,linear0oraffine1,reseton=1):

#Still some bugs in reseton==0 mode, which forces the alignment to start at one edge. In particular, opening gaps seem erroneously treated....
#Also a bug with following alignment:
##seq1='CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGTTCGCATCTGTAGTGTGTGGTACTGTGTCAGCGGCACTTAAAATGC'
##seq2='GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCTTCGCAACTGTAGTGTGTGATGGTTACTGTGTCAGCGGCTTAAAATCG'
##I think this bug is now fixed. Note that how the first row/col is treated in backtrack and score matrices is tenuous, as is the resetting.

##    import time
    
    if linear0oraffine1==0:
        if local0orglobal1==0:
            scorematrix,backtrack=seqs2scorematrix(seq1,seq2,reset=reseton)
            alignment=dmatrices2local(seq1,seq2,scorematrix,backtrack,ignorezeros=1-reseton)
        else:
            scorematrix,backtrack=seqs2scorematrix(seq1,seq2)
            alignment=dmatrices2global(seq1,seq2,scorematrix,backtrack)
    else:
        if local0orglobal1==0:
##            a=time.time()
            scorematrix,backtrack,gapmatrix,gapbacktrack=seqs2affinematrix(seq1,seq2,reset=reseton)
##            b=time.time()
            alignment=affmatrices2local(seq1,seq2,scorematrix,backtrack,gapmatrix,gapbacktrack,ignorezeros=1-reseton)
##            c=time.time()            
        else:
##            a=time.time()            
            scorematrix,backtrack,gapmatrix,gapbacktrack=seqs2affinematrix(seq1,seq2)
##            b=time.time()            
            alignment=affmatrices2global(seq1,seq2,scorematrix,backtrack,gapmatrix,gapbacktrack)
##            c=time.time()

##    print 'START'
##    print seq1,seq2
##    print b-a
##    print c-b
##    print 'STOP'            

    return alignment


#################################################
def match(a,b):

    if a==b:
        value=matchscore
    else:
        value=mismatchscore

    return value

#################################################
def complement(a,b):

    if a+b in ['AT','CG','GC','TA']:
        value=palscore
    else:
        value=palmismatch

    return value


#################################################
def gap(g):
    
    return g*gappenalty



#################################################
def affgap(g):
    
    return (g-1)*gapextend+gappenalty


#################################################
def scorealign(aln1,aln2):

    matchscore=.5
    mismatchscore=-.75    
#    mismatchscore=-.7
    gappenalty=-2.
#    gappenalty=-.4
    #gapextend=-.4
    gapextend=-1.5    

    score=0
    gapped=0
    for i in xrange(len(aln1)):        
        if aln1[i]!=aln2[i]:
            if '-'!=aln1[i] and '-'!=aln2[i]:
                score+=mismatchscore
                gapped=0
            else:
                if gapped==0:
                    score+=gappenalty
                    gapped=1
                elif gapped==1:
                    score+=gapextend
                
        else:
            score+=matchscore
            gapped=0               
            


    return score
                
            
