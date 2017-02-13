import itertools
from sage.combinat.subset import SubsetsSorted

####basic methods

def basis(n, k):
    """Returns a basis for Lk(n)"""
    parts = Compositions(n, length=k).list()
    perms = Permutations(n).list()
    basis = []
    for part in parts:
        for perm in perms:
            basic_elt = []
            i = 0
            for term in part:
                basic_elt.append(perm[i:i+term])
                i += term
            basis.append(basic_elt)
    return basis 
             
def complex(n):
    """Returns an array containing bases of Lk(n) for 0<k<=n"""
    return [basis(n, i) for i in range(1,n+1)]
    
def comult(part):
    return [[list(set), [x for x in part if x not in set]] \
    for set in SubsetsSorted(part)][1:-1]

def differentiate(elt, compl):
    length = len(elt)
    sign = 1
    basis = compl[length]
    diff = [0]*len(basis)
    for i in range(length):
        if len(elt[i]) != 1:
            com = comult(elt[i])
            for j in com:
                new_elt = elt[:i] + j + elt[i+1:]
                diff[basis.index(new_elt)] += sign
        sign *= -1
    return diff
    
def differential(n, compl):
    mat = []
    for elt in compl[n]:
        mat.append(differentiate(elt, compl))
    return matrix(mat).transpose()

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