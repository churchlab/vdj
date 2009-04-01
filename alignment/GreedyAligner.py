#Quick Comparison Algorithm From Wu et al. Implementation by Erez Lieberman.

def snakealign(seq1,seq2):

    M=len(seq1)
    N=len(seq2)
    delta=N-M
    
    fp={}
    for number in xrange(-(M+1),N+2):
        fp[number]=-1

    p=-1
    
    while fp[delta]<N:
        p=p+1
        for k in xrange(-p,delta):
            fp[k]=snake(k,max(fp[k-1]+1,fp[k+1]),seq1,seq2)

        for k in xrange(delta+p,delta,-1):                            
            fp[k]=snake(k,max(fp[k-1]+1,fp[k+1]),seq1,seq2)

        fp[delta]=snake(delta,max(fp[delta-1]+1,fp[delta+1]),seq1,seq2)

    return delta+2*p

def snake(k,y,seq1,seq2):

    M=len(seq1)
    N=len(seq2)
    
    x=y-k
    while x<M and y<N and seq1[x]==seq2[y]:
        x+=1
        y+=1

    return y


##from numarray import array,zeros

##def compare2(seq1,seq2)
##
##    M=len(seq1)
##    N=len(seq2)
##
##    #Need so big?
##    R=zeros(M,N)
##
##    i=0
##
##    while i<min(M,N) and seq[i]=seq[i+1]:
##        i=i+1
##
##    R[0][0]=i
##    
