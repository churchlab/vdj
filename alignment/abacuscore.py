# ABACUS VDJ sequence aligner (beta)
# Written by Erez Lieberman
# Modified/clarified/extended by Uri Laserson

from AlignmentCore3 import *
import time

# DEFAULTS.  CAN BE SET IN set_params()
import os
if os.path.exists('/Users/laserson/'):
	indirectory = '/Users/laserson/research/church/vdj-ome/ref-data/align/'
elif os.path.exists('/home/ul2/'):
	indirectory = '/home/ul2/vdj-ome/ref-data/align/'
outdirectory1='./'
vfile='V3.txt'
dfile='D3.txt'
jfile='J3.txt'

def set_params(refdirectory,outputdir,vref,dref,jref):
        global indirectory
        global outdirectory1
        global vfile
        global dfile
        global jfile

        indirectory = refdirectory
        outdirectory1 = outputdir
        vfile = vref
        dfile = dref
        jfile = jref

def fasta2classT2(directory,filename,tag,verbose=0,fastmode=0):
    '''ABACUS classifier
       New Tornado classifier
       Figures out VDJ using structured words. Finishes using SmithWaterman.
    '''

#    import time
    from math import ceil

    packetsize=1000
    
    t0=time.time()      # Start Timing

    # Define Seed Patterns
    patternA='111011001011010111'
    patternB='1111000100010011010111'
    patternC='111111111111'
    patternD='110100001100010101111'
    patternE='1110111010001111'    
#    patterns=[patternA,patternB,patternC]
    patterns=[patternA,patternB,patternC,patternD,patternE]    

    minipats=['111011','110111']
#    minipats=['111011','110111','10111011','11011011','1111','1110011','1100111','11111','111','1110000111']

    # Load germline VDJ database
    vlist=fa2list(indirectory,vfile,1)
    dlist=fa2list(indirectory,dfile,1)
    jlist=fa2list(indirectory,jfile,1)

    vlist=caplist(vlist,1)
    dlist=caplist(dlist,1)
    jlist=caplist(jlist,1)

    vdict=seqlist2seqdict(vlist)
    ddict=seqlist2seqdict(dlist)
    jdict=seqlist2seqdict(jlist)    

    # Generate hashes
    vld,vldsort=seqlist2wordannotlist3(vlist,patterns)
    dld,dldsort=seqlist2wordannotlist3(dlist,patterns)
    jld,jldsort=seqlist2wordannotlist3(jlist,patterns)

    dldmini,dldsortmini=seqlist2wordannotlist3(dlist,minipats)

    
    t1=time.time()      # Record time for database initialization
    
    # NOTE: After this point, runtime is seqlist dependent

    # Load query sequences
    fullseqlist=fa2list(directory,filename,1)
    totseqs=len(fullseqlist[0])
    tfullseqlist=transpose(fullseqlist)

    output1 = open(outdirectory1+'newclass'+tag+'.class', 'w')
    output2 = open(outdirectory1+'trialclass'+tag+'.class', 'w')        
    
    dcount=0
    fail=0
    for packetnum in range(1,1+int(ceil(float(len(tfullseqlist))/packetsize))):

        print "Packet", packetnum
        seqlist=transpose(tfullseqlist[packetsize*(packetnum-1):packetsize*packetnum])

    
        seqdict=seqlist2seqdict(seqlist)

        sld,sldsort=seqlist2wordannotlist3(seqlist,patterns)

        t2=time.time()



        
        vscores={}
        jscores={}
        dscores={}   
        for key in seqdict.keys():   # for each query sequence
            vscores[key]={}
            maxscore=-1        
            for vseg in vdict.keys():   # for each germline V sequence
                score=0
                for pattern in patterns:
                    # for each pattern, query, and germline sequence, how many
                    # shared keys/words are there?
                    score+=len(vldsort[vseg][pattern] & sldsort[key][pattern])
                vscores[key][vseg]=score
                
            
            jscores[key]={}
            maxscore=-1
            for jseg in jdict.keys():
                score=0
                for pattern in patterns:
                    score+=len(jldsort[jseg][pattern] & sldsort[key][pattern])
                jscores[key][jseg]=score
                    


        vthresh=.75
        if fastmode==1:
            vthresh=1
        
        goodvscores={}
        goodjscores={}
        bestv={}
        bestj={}
        bestd={}
        for key in seqdict.keys(): 
            goodvscores[key]={}
            maxscore=-1
            bestv[key]=''
            goodvseglist=transpose(dict2kvsorted(vscores[key],1)[-10:])[0][::-1]

##            if fastmode==1:
##                    goodvseglist=[goodvseglist[0]]
                            
            for goodvseg in goodvseglist:
                if vscores[key][goodvseglist[0]]>0 and float(vscores[key][goodvseg])/vscores[key][goodvseglist[0]]>=vthresh:
                    seqprocdict={}
                    for pattern in patterns:
                        seqprocdict[pattern]=[[sld[key][pattern],sldsort[key][pattern]],[vld[goodvseg][pattern],vldsort[goodvseg][pattern]]]
                    aligned=semidpalign(seqdict[key],vdict[goodvseg],patterns,seqproc=seqprocdict,prewindow=10,postwindow=10)
                    if aligned!=[]:
                        score2=scorealign(aligned[0],aligned[1])
                    else:
                        break
                    
                  
                    if score2>maxscore:
                        best=goodvseg
                        maxscore=score2
                        bestv[key]=best

                else: break                        



            jthresh=.75
            if fastmode==1:
                jthresh=1
            
            goodjscores[key]={}
            maxscore=-1
            bestj[key]=''

            goodjseglist=transpose(dict2kvsorted(jscores[key],1)[-3:])[0][::-1]
            
##            if fastmode==1:
##                    goodjseglist=[goodjseglist[0]]
                
            for goodjseg in goodjseglist:
                if jscores[key][goodjseglist[0]]>0 and float(jscores[key][goodjseg])/jscores[key][goodjseglist[0]]>=jthresh:
                    
            
                
                    seqprocdict={}
                    for pattern in patterns:
                        seqprocdict[pattern]=[[sld[key][pattern],sldsort[key][pattern]],[jld[goodjseg][pattern],jldsort[goodjseg][pattern]]]


                    aligned=semidpalign(seqdict[key],jdict[goodjseg],patterns,seqproc=seqprocdict,prewindow=10,postwindow=10)
                    if aligned!=[]:
                        score2=scorealign(aligned[0],aligned[1])
                    else:
                        break

                    if score2>maxscore:
                        best=goodjseg
                        maxscore=score2
                        bestj[key]=best

                else: break                        



        
#        for key in seqdict.keys():
        for key in seqlist[0]:
            if bestv[key]!='':
                seqprocdict={}
                for pattern in patterns:
                    seqprocdict[pattern]=[[sld[key][pattern],sldsort[key][pattern]],[vld[bestv[key]][pattern],vldsort[bestv[key]][pattern]]]
                va=semidpalign(seqdict[key],vdict[bestv[key]],patterns,seqproc=seqprocdict,prewindow=10,postwindow=10)
            else:
                fail+=len(key.split(','))

            if bestj[key]!='':
                seqprocdict={}
                for pattern in patterns:
                    seqprocdict[pattern]=[[sld[key][pattern],sldsort[key][pattern]],[jld[bestj[key]][pattern],jldsort[bestj[key]][pattern]]]
                ja=semidpalign(seqdict[key],jdict[bestj[key]],patterns,seqproc=seqprocdict,prewindow=10,postwindow=10)
            else:
                fail+=len(key.split(','))

            try:
                # D segment aligner code.
                if bestv[key]!='' and bestj[key]!='':            

                    dpiece=seqdict[key].replace(va[0].replace('-',''),'|').replace(ja[0].replace('-',''),'|').split('|')[1]
                    dpieceannot,dpiecemini=seq2wordannot3(dpiece,minipats)

                    dscores[key]={}
                    maxscore=-1
                    for dseg in ddict.keys():
                        score=0
                        for pattern in minipats:
                            score+=len(dldsortmini[dseg][pattern] & dpiecemini[pattern])
                        dscores[key][dseg]=score




                    gooddseglist=transpose(dict2kvsorted(dscores[key],1))[0][::-1]
                    crudemax=dscores[key][gooddseglist[0]]
                    
                    maxscore=4
                    bestda=[[['']]]
                    dthresh=.75
                    if fastmode==1:
                        dthresh=1

##                    if fastmode==1:
##                        gooddseglist=[gooddseglist[0]]
                

                    for gooddseg in gooddseglist:
                        if dscores[key][gooddseg]>=crudemax*dthresh:         
                            da=align(dpiece,ddict[gooddseg],0,0)
                            if da[0][1]>maxscore:
                                maxscore=da[0][1]
                                bestda=da
                                bestdseq=ddict[gooddseg]
                                bestseg=gooddseg
                        else:
                            break
                        
                    if bestda[0][0][0]!='':
                        bestd[key]=bestseg
                        da=bestda[0][0]
                        dcount+=len(key.split(','))
                    else:
                        bestd[key]=''
            except:
                weirdobug=1
                print key
                print seqdict[key]
                print '--'

            if verbose==1:
                if bestv[key]!='' and bestj[key]!='' and bestd[key]!='':

                    seqv=vdict[bestv[key]].replace(va[1].replace('-',''),va[1])
                    seqd=ddict[bestd[key]].replace(da[1].replace('-',''),da[1])	    
                    seqj=jdict[bestj[key]].replace(ja[1].replace('-',''),ja[1])
                    seqq=seqdict[key].replace(va[0].replace('-',''),va[0]).replace(ja[0].replace('-',''),ja[0]).replace(da[0].replace('-',''),da[0])

                    seqvbar=''
                    seqdbar=''	    
                    seqjbar=''
                    

                    j=0
                    for i in xrange(len(va[0])):            
                        if '-'!=va[0][i] and '-'!=va[1][i]:
                            seqvbar+='|'
                        else:
                            seqvbar+=' '

                    for i in xrange(len(da[0])):
                        if '-'!=da[0][i] and '-'!=da[1][i]:
                            seqdbar+='|'
                        else:
                            seqdbar+=' '
                            
                    for i in xrange(len(ja[0])):
                        if '-'!=ja[0][i] and '-'!=ja[1][i]:
                            seqjbar+='|'
                        else:
                            seqjbar+=' '

                    try:
                        output1.write('>'+key+'\n')
                        output1.write(bestv[key]+', '+bestj[key]+' Alignment:'+'\n')
                        output1.write(' '*(70-seqv[-70:].index(va[1])+seqq.index(va[0]))+seqv[-70:]+'\n')
                        output1.write(' '*(70+seqq.index(va[0]))+seqvbar+'\n')
                        output1.write(' '*(70)+seqq+'\n')
                        output1.write(' '*(70+seqq.index(ja[0]))+seqjbar+'\n')
                        output1.write(' '*(70+seqq.index(ja[0])-seqj.index(ja[1]))+seqj+'\n')
                        output1.write(bestd[key]+' Alignment:'+'\n')
                        output1.write(' '*(70+seqq.index(va[0]))+'V'*len(seqvbar)+' '*(seqq.index(da[0])-seqq.index(va[0])-len(seqvbar))+'D'*len(seqdbar)+' '*(seqq.index(ja[0])-seqq.index(da[0])-len(seqdbar))+'J'*len(seqjbar)+'\n')            
                        output1.write(' '*(70)+seqq+'\n')
                        output1.write(' '*(70+seqq.index(da[0]))+seqdbar+'\n')
                        output1.write(' '*(70+seqq.index(da[0])-seqd.index(da[1]))+seqd+'\n')
                        output1.write('----'+'\n')
                        
                    except:
                        print vdict[bestv[key]]
                        print va
                        print seqdict[key]
        
            output2.write(key+'\t')
            try: output2.write(bestv[key]+'\t')
            except: output2.write('\t')
            try: output2.write(bestd[key]+'\t')
            except: output2.write('\t')
            try: output2.write(bestj[key])
            except: output2.write('')
            output2.write('\n')                   
            
                    

	    
    t3=time.time()



##    t5=time.time()
##
    output1.close()
    output2.close()
## 
##    
####    classdict2file(dictannot,'trialclass'+tag+'.class')
##
##
##
##    t6=time.time()

    print 'Ds Called: '+`dcount`
    print 'Failures: '+`fail`
##    print 'Uniques: '+`len(seqlist[0])`
    print 'Total Seqs: '+`totseqs`    


    

    print t1-t0
    print t2-t1
    print t3-t2
##    print t4-t3
##    print t5-t4
##    print t6-t5
    print '----'
    print t3-t0
    




################################################
def fa2list(directory,filename,fullname):
### Converts a fasta file to a list; either keeping the Fullname or just the first 7 positions
    
    input = open(directory+filename,'r')

    seqlist=[]
    names=[]
    strand=''
    string=''

    while 1:
        string = input.readline()
        if not string:
            strand=strand+lastchar
            seqlist.append(strand.replace('\n',''))
            break
        if string[0]=='>':
            if fullname==1:
                names.append(string[1:len(string)-1])
            else:
                names.append(string[1:7])
                
            seqlist.append(strand)
            string = input.readline()
            strand=string[0:len(string)-1]
        else:
            strand = strand+string[0:len(string)-1]

        lastchar=string[len(string)-1]


    seqlist=seqlist[1:len(seqlist)]
    seqlist=[names,seqlist]
    input.close()
    return seqlist


#################################################
def caplist(seqlist,upper):
### Caps or uncaps a seqlist.

    for i in range(len(seqlist[1])):
        if upper==1:
            seqlist[1][i]=seqlist[1][i].upper()
        else:
            seqlist[1][i]=seqlist[1][i].lower()

    return seqlist


#################################################
def seqlist2seqdict(seqlist):
### Converts a seqlist into a dictionary sorted by sequence name.

    tseqlist=transpose(seqlist)
    seqdict={}
    for seq in tseqlist:
        seqdict[seq[0]]=seq[1]

    return seqdict


#################################################
def seqlist2wordannotlist3(seqlist,patternlist):

    seqlistwannot={}
    seqlistsort={}

    i=0    
    for seq in transpose(seqlist):
        i=i+1
        if i%1000==0:
            print i
##            print seq

        seqannot,sort=seq2wordannot3(seq[1],patternlist)
        seqlistwannot[seq[0]]={}
        seqlistsort[seq[0]]={}
        for pattern in patternlist:
            seqlistwannot[seq[0]][pattern]=seqannot[pattern]
            seqlistsort[seq[0]][pattern]=sort[pattern]
        

    return seqlistwannot,seqlistsort


#################################################
def transpose(lists):
### transposes matrices, lists 

    newlist=[]
    i=0;
    while i<len(lists[0]):
        newentry=[]
        j=0;
        while j<len(lists):
            newentry.append(lists[j][i])
            j=j+1
            
        newlist.append(newentry)
        i=i+1

    return newlist


#####################################
def dict2kvsorted(nmers,keyorvalue):
### Sorts a listed dictionary by key or by value
    
    from operator import itemgetter
    yi=nmers.items()
    yi.sort(key = itemgetter(keyorvalue))

    return yi


#################################################
def semidpalign(seq1,seq2,patternlist,seqproc=[],prewindow=0,postwindow=0,verbose=0):


    a=seqs2anchors(seq1,seq2,patternlist,seqproc=seqproc)
    b=anchors2exanchorsfast(seq1,seq2,a,verbose=verbose)
    try: aligned=exanchors2swanchors(seq1,seq2,b,prewindow,postwindow,verbose=verbose)
    except:
        aligned=[]
    

    t3=time.time()
    

    return aligned


#################################################
def seq2wordannot3(seq,patternlist):

    seqannot={}
    patlens=[]
    for pattern in patternlist:
        patlens.append(len(pattern))
        seqannot[pattern]={}

    maxpat=max(patlens)
        
    for i in xrange(len(seq)):
        word=seq[i:i+maxpat]
        for pattern in patternlist:
            patlen=len(pattern)
            if len(word)>=patlen:
                key=''
                for j in xrange(patlen):
                    if pattern[j]=='1':
                        key+=word[j]

                oldnmers=seqannot[pattern].get(key,[])
                seqannot[pattern][key]=oldnmers+[i]



    sortlist={}
    for pattern in patternlist:
        sortlist[pattern]=set(seqannot[pattern].keys())

    return seqannot,sortlist


#################################################
def seqs2anchors(seq1,seq2,patternlist,seqproc=[]):

    if seqproc==[]:
        seqannots1,sorts1=seq2wordannot3(seq1,patternlist)
        seqannots2,sorts2=seq2wordannot3(seq2,patternlist)
        seqproc={}
        for pattern in patternlist:
            seqproc[pattern]=[[seqannots1[pattern],sorts1[pattern]],[seqannots2[pattern],sorts2[pattern]]]

    anchors=[]
    for pattern in patternlist:
        seqannot1,sort1=seqproc[pattern][0]
        seqannot2,sort2=seqproc[pattern][1] 

        window=len(pattern)
        anchors+=[[[entry,entry+window],[entry2,entry2+window]] for key in sort1 & sort2 for entry in seqannot1[key] for entry2 in seqannot2[key]]


    return anchors


#################################################
def anchors2exanchorsfast(seq1,seq2,anchors,verbose=0):

   # print anchors
    newanchors=[]
    for anchor in anchors:
        needed=1
        for anchornew in newanchors:
            if intervalcontained(anchor[0],anchornew[0])==1 and intervalcontained(anchor[1],anchornew[1])==1:
                needed=0
                break
            
        if needed==1:
            newanchor=anchor2extend(seq1,seq2,anchor)
            newanchors.append(newanchor)

    
    if len(newanchors)>1:
        newanchors=anchors2wellordered(seq1,seq2,newanchors)

    
    return newanchors


#################################################
def exanchors2swanchors(seq1,seq2,exanchors,prewindow,postwindow,verbose=0):

    seqaln1=''
    seqaln2=''



    
    for i in range(len(exanchors)-1):
        alnpiece=align(seq1[exanchors[i][0][1]:exanchors[i+1][0][0]],seq2[exanchors[i][1][1]:exanchors[i+1][1][0]],1,0)
        
        if i>0:

            if exanchors[i][0][0]<exanchors[i-1][0][1] or exanchors[i][1][0]<exanchors[i-1][1][1]:
                diff1=max(exanchors[i-1][0][1]-exanchors[i][0][0],0)
                diff2=max(exanchors[i-1][1][1]-exanchors[i][1][0],0)

                if diff2>diff1:
                    seqaln2+='-'*(diff2-diff1)
                if diff1>diff2:
                    seqaln1+='-'*(diff1-diff2)

            if exanchors[i][0][0]<exanchors[i-1][0][1]:
                seqaln1+='&'+seq1[exanchors[i-1][0][1]:exanchors[i][0][1]]+'|'+alnpiece[0][0]
            else:
                seqaln1+='&'+seq1[exanchors[i][0][0]:exanchors[i][0][1]]+'|'+alnpiece[0][0]

           

            if exanchors[i][1][0]<exanchors[i-1][1][1]:    
                seqaln2+='&'+seq2[exanchors[i-1][1][1]:exanchors[i][1][1]]+'|'+alnpiece[0][1]
            else:
                seqaln2+='&'+seq2[exanchors[i][1][0]:exanchors[i][1][1]]+'|'+alnpiece[0][1]
        else:
            seqaln1+='&'+seq1[exanchors[i][0][0]:exanchors[i][0][1]]+'|'+alnpiece[0][0]
            seqaln2+='&'+seq2[exanchors[i][1][0]:exanchors[i][1][1]]+'|'+alnpiece[0][1]          




    if len(exanchors)>1:
        if exanchors[-1][0][0]<exanchors[-2][0][1] or exanchors[-1][1][0]<exanchors[-2][1][1]:
            diff1=max(exanchors[-2][0][1]-exanchors[-1][0][0],0)
            diff2=max(exanchors[-2][1][1]-exanchors[-1][1][0],0)
            if diff2>diff1:
                seqaln2+='-'*(diff2-diff1)
            if diff1>diff2:
                seqaln1+='-'*(diff1-diff2)  
            

    if len(exanchors)>1 and exanchors[-1][0][0]<exanchors[-2][0][1]:
        seqaln1+='&'+seq1[exanchors[-2][0][1]:exanchors[-1][0][1]]
    else:
        seqaln1+='&'+seq1[exanchors[-1][0][0]:exanchors[-1][0][1]]
    if len(exanchors)>1 and exanchors[-1][1][0]<exanchors[-2][1][1]:    
        seqaln2+='&'+seq2[exanchors[-2][1][1]:exanchors[-1][1][1]]
    else:
        seqaln2+='&'+seq2[exanchors[-1][1][0]:exanchors[-1][1][1]]



    if len(exanchors)>1:
        seqaln1,seqaln2=draftalignment2paredalignment(seqaln1+'|',seqaln2+'|')


    begin1=seq1.index(seqaln1.replace('-','').replace('&','').replace('|',''))
    end1=begin1+len(seqaln1.replace('-','').replace('&','').replace('|',''))
    begin2=seq2.index(seqaln2.replace('-','').replace('&','').replace('|',''))
    end2=begin2+len(seqaln2.replace('-','').replace('&','').replace('|',''))    

    alnstart=align(seq1[max(begin1-prewindow,0):begin1][::-1],seq2[max(begin2-prewindow,0):begin2][::-1],0,0,reseton=0)
    alnend=align(seq1[end1:end1+postwindow],seq2[end2:end2+postwindow],0,0,reseton=0)

    
    seqaln1=alnstart[0][0][0][::-1]+seqaln1+alnend[0][0][0]
    seqaln2=alnstart[0][0][1][::-1]+seqaln2+alnend[0][0][1]

    seqaln1=seqaln1.replace('&','').replace('|','')
    seqaln2=seqaln2.replace('&','').replace('|','')    

    return [seqaln1,seqaln2]


#################################################
def intervalcontained(interval1,interval2):

    contain=0
    if interval1[0]>=interval2[0] and interval1[1]<=interval2[1]:
        contain=1

    return contain


#################################################
def anchor2extend(seq1,seq2,anchor):


    start=min(anchor[0][0],anchor[1][0])
    if start>0:
        for i in xrange(1,start+1):
            if seq1[anchor[0][0]-i]!=seq2[anchor[1][0]-i]:
                i=i-1
                break
    else:
        i=0

    end=min([len(seq1)-anchor[0][1],len(seq2)-anchor[1][1]])
    if end>0:
        for j in xrange(end):
            if seq1[anchor[0][1]+j]!=seq2[anchor[1][1]+j]:
                j=j-1
                break
    else:
        j=-1

    anchornew=[[anchor[0][0]-i,anchor[0][1]+j+1],[anchor[1][0]-i,anchor[1][1]+j+1]]
    
    return anchornew


#################################################
def anchors2wellordered(seq1,seq2,anchors):
# This piece of code makes sure the anchors have a reasonable order. The check is pretty crude and it is imaginable that it occasionally fails...    
    
    anchorlengths=[]
    i=0
    for anchor in anchors:
        anchorlengths.append([i,scorealign(seq1[anchor[0][0]:anchor[0][1]],seq2[anchor[1][0]:anchor[1][1]])])
        i+=1

    priority=transpose(sortbycol(anchorlengths,1,-1))

    newanchors=[]
    for index in priority[0]:
        good=1
        for newanchor in newanchors:
            if intervaloverlap(newanchor[0],anchors[index][0])==1 or intervaloverlap(newanchor[1],anchors[index][1])==1:
                good=0
                break
            elif (newanchor[0][1]<anchors[index][0][1])==(newanchor[1][0]>anchors[index][1][0]):
                good=0
                break
         
        if good==1:
            newanchors.append(anchors[index])
    
    newanchors.sort()
    
    return newanchors


#################################################
def sortbycol(metalist,col,direction):
### Sorts a list by a particular column.
    
    metalist=transpose(metalist)
    newlist=[metalist[col]]
    for alist in metalist:
        newlist.append(alist)

    newlist=transpose(newlist)
    newlist.sort()
    if direction==-1:
        newlist.reverse()
        
    newlist=transpose(newlist)    
    newlist.pop(0)
    newlist=transpose(newlist)
    
    return newlist


#################################################
def intervaloverlap(interval1,interval2):

    overlap=1
    if interval1[1]<interval2[0] or interval1[0]>interval2[1]:
        overlap=0

    return overlap


#################################################
def draftalignment2paredalignment(seqaln1,seqaln2):

    anchorpiecesbef1=seqaln1.split('&')[1:]
    anchorpiecesbef2=seqaln2.split('&')[1:]

    anchorscores={}
    for i in xrange(len(anchorpiecesbef1)):
        anchorscores[i]=scorealign(anchorpiecesbef1[i].split('|')[0],anchorpiecesbef2[i].split('|')[0])

    

    goodorbadbefore={}    
    for i in xrange(len(anchorpiecesbef1)):
        goodorbadbefore[i]=scorealign(anchorpiecesbef1[i].replace('|',''),anchorpiecesbef2[i].replace('|',''))


    scoresbefore={}
    for i in xrange(len(anchorpiecesbef1)):
        scoresbefore[i]={}
        scoresbefore[i][0]=0        
        score=0
        for step in range(i):
            score+=goodorbadbefore[i-step-1]
            scoresbefore[i][step+1]=score



    anchorpiecesaft1=seqaln1.split('|')[:-1]
    anchorpiecesaft2=seqaln2.split('|')[:-1]

    goodorbadafter={}    
    for i in xrange(len(anchorpiecesbef1)):
        goodorbadafter[i]=scorealign(anchorpiecesaft1[i].replace('&',''),anchorpiecesaft2[i].replace('&',''))

    scoresafter={}
    for i in xrange(len(anchorpiecesbef1)):
        scoresafter[i]={}
        scoresafter[i][0]=0
        score=0
        for step in range(len(goodorbadafter)-i-1):
            score+=goodorbadafter[i+step+1]
            scoresafter[i][step+1]=score


    topintervals={}
    topscore={}
    for i in xrange(len(anchorpiecesbef1)):
        bestbefore=dict2kvsorted(scoresbefore[i],1)[-1]
        down=bestbefore[0]
        bestafter=dict2kvsorted(scoresafter[i],1)[-1]
        up=bestafter[0]
        bestscore=anchorscores[i]+bestbefore[1]+bestafter[1]
        topintervals[i]=[down,up,bestscore]
        topscore[i]=bestscore

    bestinterval=dict2kvsorted(topscore,1)[-1]


    topintervals[bestinterval[0]][0]

    bestanchor1=anchorpiecesbef1[bestinterval[0]].split('|')[0]
    bestanchor2=anchorpiecesbef2[bestinterval[0]].split('|')[0]    

    seqbefore1=''
    seqbefore2=''    
    for step in xrange(0,topintervals[bestinterval[0]][0]):
        seqbefore1=anchorpiecesbef1[i-step-1]+seqbefore1
        seqbefore2=anchorpiecesbef2[i-step-1]+seqbefore2       

    seqafter1=''
    seqafter2=''    
    for step in xrange(1,topintervals[bestinterval[0]][1]+1):
        seqafter1+=anchorpiecesaft1[i+step-1]
        seqafter2+=anchorpiecesaft2[i+step-1]        

    bestalign1=seqbefore1+bestanchor1+seqafter1
    bestalign2=seqbefore2+bestanchor2+seqafter2


    return bestalign1,bestalign2

#################################################
def caps(sequence):
### Puts a sequence in caps.
    seqnew=sequence.replace('a','A').replace('c','C').replace('g','G').replace('t','T')

    return seqnew


#################################################
def invcomp(string):
### Computes the inverse complement of a string

    string=caps(string)
    s1=string.replace('A','t')
    s2=s1.replace('T','a')
    s3=s2.replace('C','g')
    s4=s3.replace('G','c')
    s5=s4.replace('a','A')
    s6=s5.replace('c','C')
    s7=s6.replace('g','G')
    s8=s7.replace('t','T')
    i=1;
    s9=''
    while i<=len(s8):
        s9=s9+s8[len(s8)-i]
        i=i+1
        
    string=s9
    return string


#################################################
def fliplist(seqlist):
### Takes a Seqlist and flips the strand around.

    for i in range(len(seqlist[1])):
        seqlist[1][i]=invcomp(seqlist[1][i].upper())

    return seqlist



##############################################
## Uri addendum

from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def uriseqlist2erezseqlist(uriseqlist):
        erezseqlist = []
        erezseqlist.append([])
        erezseqlist.append([])
        for seq in uriseqlist:
                erezseqlist[0].append(seq.description)
                erezseqlist[1].append(seq.seq.data.upper())
        return erezseqlist

def erezseqlist2uriseqlist(erezseqlist):
        tseqlist = transpose(erezseqlist)
        uriseqlist = []
        for seq in tseqlist:
                uriseqlist.append( SeqRecord( Seq(seq[1],generic_dna), description=seq[0] ) )
        return uriseqlist

#################################
def seqlist2positiveseqids(seqlist):

        pattern='111111111111'
        patterns=[pattern]
        
        # NOTE the dependencies on the global variables (change in the future)
        vlist=fa2list(indirectory,vfile,1)
        jlist=fa2list(indirectory,jfile,1)

        vlist=caplist(vlist,1)
        jlist=caplist(jlist,1)

        # Generate hashes
        posvld,posvldsort=seqlist2wordannotlist3(vlist,patterns)
        posjld,posjldsort=seqlist2wordannotlist3(jlist,patterns)
        negvld,negvldsort=seqlist2wordannotlist3(fliplist(vlist),patterns)
        negjld,negjldsort=seqlist2wordannotlist3(fliplist(jlist),patterns)

        # collect possible keys
        posset=set([])
        for key in posvldsort.keys():
                posset=(posset|posvldsort[key][pattern])
        for key in posjldsort.keys():
                posset=(posset|posjldsort[key][pattern])

        negset=set([])
        for key in negvldsort.keys():
                negset=(negset|negvldsort[key][pattern])
        for key in negjldsort.keys():
                negset=(negset|negjldsort[key][pattern])

        # get keys unique to positive or negative versions of reference set
        possetnew=posset-negset
        negsetnew=negset-posset

        posset=possetnew
        negset=negsetnew

        strandid = [1]*len(seqlist[1])
        for i in xrange(len(seqlist[1])):
                seqwords=seq2wordannot3(seqlist[1][i],['111111111111'])[1]['111111111111']
                if len(negset & seqwords) > len(posset & seqwords):
                        strandid[i] = -1

        return strandid
