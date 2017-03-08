import itertools
from sage.combinat.subset import SubsetsSorted

####basic methods

def fact(n):
    if n == 0:
        return 1
    else:
        return n*fact(n-1)

def basis(n, deg):
    """Returns the standard basis for L(deg-1)(n)"""
    b = []
    for comp in Compositions(n, length=1+deg):    
        for perm in Permutations(n):
            cur = []
            ind = 0
            for i in range(1+deg):
                 cur.append(perm[ind:ind+comp[i]])
                 ind += comp[i]
            b.append(cur)
    return b

def comult(elt):
    """Comultiply a basic element of L0(n)"""
    res = []
    n = len(elt)
    r = range(n)
    for sset in sage.combinat.subset.SubsetsSorted(r):
        if len(sset) != 0 and len(sset) != n:
            comp = [i for i in r if i not in sset]
            res.append([[elt[i] for i in sset], [elt[i] for i in comp]])
    return res
    
def complex(n):
    return [basis(n, i) for i in range(n)]
    
def differentiate(elt, compl):
    """Returns the coordinates of the differential of basic element elt,
    in terms of the basis for the complex given by compl"""
    n = len(elt)
    codomainBasis = compl[n]
    res = [0 for _ in range(len(codomainBasis))]
    sign = 1
    for i in range(n):
        if len(elt[i]) != 1:
            for term in comult(elt[i]):
                cur = elt[:i] + term + elt[i+1:]
                res[codomainBasis.index(cur)] += sign
        sign *= -1
    return res
    
def differential(n, compl):
    """Returns the n-th differential matrix of the complex"""
    rows = []
    for basic_elt in compl[n]:
        rows.append(differentiate(basic_elt, compl))
    return matrix(rows).transpose()

#example: calculating betti numbers for n=3
#n = 3
#ch = ChainComplex([differential(i,complex(n)) for i in range(n-1)])
#ch.betti()


####spectral sequence stuff

def perms_by_fixpoints(n):
    """Lists all permutations of n elements sorted by 
    number of fixpoints, increasingly"""
    perms = [[] for _ in range(n+1)]
    for i in Permutations(n):
        perms[i.number_of_fixed_points()].append(i)
    return perms

#cache all permutations classified by number of fixpoints up to n=10
perms = [perms_by_fixpoints(i) for i in range(1,10)]


def apply(perm, set):   
    return [set[perm[i]-1] for i in range(len(set))]

def filt_basis(n, tensors, fixpoints):
    """Returns a list of basic elements of Ltensors(n) with a 
    prescribed number of fixpoints"""
    b = []
    comps = Compositions(n, length=tensors)
    fixes = Compositions(fixpoints, length=tensors, min_part=0)
    for comp in comps:
        for fix in fixes:
            if reduce(lambda x,y: x and y, [comp[i]>=fix[i] \
            for i in range(tensors)]):
                parts = [[list(i) for i in x] \
                for x in OrderedSetPartitions(range(1,n+1), comp)]
                for part in parts:
                    cur_perms = []
                    for i in range(len(part)):
                        l = len(part[i])
                        cur_perms.append([apply(p, part[i]) \
                                          for p in perms[l-1][fix[i]]])
                    b += [list(elt) for elt in itertools.product(*cur_perms)]
    return b

#example: the dimensions of the associated graded objects for n=5
#n=5
#for parts in range(1,n+1):
#    print [len(basis(n, parts, fixes)) for fixes in range(n+1)]

#after a while, this outputs:
#[44, 45, 20, 10, 0, 1]
#[40, 150, 120, 140, 0, 30]
#[0, 90, 120, 360, 0, 150]
#[0, 0, 0, 240, 0, 240]
#[0, 0, 0, 0, 0, 120]

###d, the first component of the comultiplication map

def d(n):
    pms = Permutations(n).list()
    rows = []
    for p in pms:
        r = [[p[i]]+[x for x in p if x!=p[i]] for i in range(n)]
        rows.append([int(i in r) for i in pms])
    return Matrix(rows)