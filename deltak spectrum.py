def eval(elt, tup):
    """Evaluates the permutation elt pointwise in the tuple tup"""
    return [elt.tuple()[i-1] for i in tup]

def todict(lst):
    res = {}
    for i in lst:
        if i in res.keys():
            res[i] += 1
        else:
            res[i] = 1
    return res

def kfixed(elt, k):
    """Returns the number of k-subsets fixed by the permutation elt"""
    cyc_type = todict([len(cyc) for cyc in elt.cycle_tuples(singletons=True)])
    parts = Partitions(k)
    total = 0
    for part in parts:
        d = todict(part)
        if reduce((lambda x,y: x and y), [i in cyc_type.keys() \
                                            for i in d.keys()]):
            total += reduce((lambda x,y: x*y), \
                    [binomial(cyc_type[i], d[i]) for i in d.keys()])
    return total
 
def spectrum(k, n):
    """Returns the spectrum of Delta_k for n elements as a dictionary"""
    g = SymmetricGroup(n)
    cl = g.conjugacy_classes()
    dic = {}
    for cls in cl:
        fixes = kfixed(cls.random_element(), k)
        if fixes in dic.keys():
            dic[fixes] += len(cls)
        else:
            dic[fixes] = len(cls)
    return dic 

def powersums(dic, l):
    """Returns the l-th power sum of the elements in dic.
    The i-th key appears with multiplicity dic[i]"""
    return sum([dic[i]*i**l for i in dic.keys()])

###sample run, where spec3 is a list of dictionaries,
###where the n-th dictionary contains the spectrum of Delta3
###for n elements:

#for n in range (3, 11):
#    for l in range(1, 11):
#        print powersums(spec3[n],l)/fact(n),
#    print ''

#this outputs

#1 1 1 1 1 1 1 1 1 1 
#1 2 5 15 51 187 715 2795 11051 43947 
#1 3 15 107 923 8683 84715 838827 8355243 83420843 
#1 4 28 328 5200 94624 1822528 35909248 713923840 14244657664 
#1 4 35 577 14006 414944 13496125 457605607 15799074116 549739560814 
#1 4 38 749 23846 1000904 48658458 2546625904 138133674626 7620885721484 
#1 4 39 830 30791 1628742 106448149 7836417880 615040519755 49911022879610 
#1 4 39 855 34168 2071063 165527008 15736107435 1658451693436 185368145380911 