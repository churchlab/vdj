import rpy2.robjects as robjects
r = robjects.r

qbinom = r['qbinom']
pbinom = r['pbinom']

def binomial(p, lst, prob):
    invBinom = qbinom(p,robjects.IntVector(lst), prob)
    cumulative = pbinom(p,invBinom, prob)
    res = [1]
    for ii in lst:
        if cumulative[ii-1] < p or invBinom[ii-1] == ii:
            res.append(int(invBinom[ii-1] + 1))
        else:
            res.append(int(invBinom[ii-1]))
    return res

