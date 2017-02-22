def fact(m):
    if m==0: return 1
    else: return m*fact(m-1)

class GroupAlgebraElt:
    def __init__(self, coeffs, order):
        self.coeffs = coeffs
        self.order = order
        self.ncoeffs = fact(order)
    
    def __add__(self, other):
        return GroupAlgebraElt([self.coeffs[i] + other.coeffs[i] \
            for i in range(self.ncoeffs)], self.order)
    
    def __sub__(self, other):
        return GroupAlgebraElt([self.coeffs[i] - other.coeffs[i] \
            for i in range(self.ncoeffs)], self.order)
        
    def __mul__(self, other):
        glist = SymmetricGroup(self.order).list()
        new_coeffs = [0 for _ in range(self.ncoeffs)]
        for indi, coeffi in enumerate(self.coeffs):
            if coeffi != 0:
                for indj, coeffj in enumerate(other.coeffs):
                    if coeffj != 0:
                        new_coeffs[glist.index(glist[indj]*glist[indi])] \
                        += coeffi*coeffj
        return GroupAlgebraElt(new_coeffs, self.order)
    
    def __repr__(self):
        glist = SymmetricGroup(self.order).list()
        s = ''
        flag = False
        for indi, coeffi in enumerate(self.coeffs):
            if coeffi != 0:
                if flag:
                    s += ' + ' + str(coeffi) + ' ' + str(glist[indi])          
                else:
                    s += str(coeffi) + ' ' + str(glist[indi])
                    flag = True
        if flag:
            return s
        else:
            return '0'
            
    def leftMultMatrix(self):
        glist = SymmetricGroup(self.order).list()
        cols = []
        for i in range(self.ncoeffs):
            col = [0 for _ in range(self.ncoeffs)]
            for indi, coeffi in enumerate(self.coeffs):
                 col[glist.index(glist[i] * glist[indi])] += coeffi
            cols.append(col)
        return Matrix(cols).transpose()
        
    def spectrum(self):
        return self.leftMultMatrix().charpoly().roots() 
        
    def eigenvectors(self):
        return self.leftMultMatrix().eigenvectors_right()
        
    def termwiseConjugation(self, elt):
        coeffs = [0 for _ in range(self.ncoeffs)]
        for indi, coeffi in enumerate(self.coeffs):
            if coeffi != 0:
                coeffs[glist.index(elt*glist[indi]*elt.inverse())] = coeffi
        return GroupAlgebraElt(coeffs, self.order)
        
def beta(k,n):
    """Returns the element 1 + (12) + ... + (12...k)
    in the group algebra kSn"""
    g = SymmetricGroup(n)
    glist = g.list()
    coeffs = [0 for _ in range(fact(n))]
    coeffs[0] = 1 #identity
    elt = g('()')
    for i in range(k-1):
        elt =  g('(' + str(i+1) + ',' + str(i+2) + ')') * elt
        coeffs[glist.index(elt)] = 1
    return GroupAlgebraElt(coeffs, n)

#####some ideas from babai's paper, to compute the spectrum in each
#####right irreducible submodule

import itertools
def prodcount(n, k):
    g = SymmetricGroup(n)
    glist = g.list()
    elts = [glist[indi] for indi, coeffi in enumerate(beta(n,n).coeffs) \
            if coeffi != 0]
    classes = dict.fromkeys(g.conjugacy_classes(), 0)
    for factors in itertools.product(*[elts for _ in range(k)]):
        classes[reduce((lambda x,y: x*y), factors).conjugacy_class()] += 1
    return classes
    
def chartally(chr, vals):
    g = SymmetricGroup(n)
    c = g.conjugacy_classes()
    return sum([chr[i]*vals[c[i]] for i in range(len(chr))])


n = 5
g = SymmetricGroup(n)
glist = g.list()
for part in Partitions(n):
    #left/right issues?
    chr = SymmetricGroupRepresentation(part).to_character().values()
    dim = int(chr[0])
    exp = 1
    total = chartally(chr, prodcount(n, exp))
    candidates = [p + [0 for _ in range(dim-len(p))] for p \
                  in Partitions(int(total), max_length=dim)]
    while (len(candidates) > 1) and (exp <= dim):
        exp += 1
        total = chartally(chr, prodcount(n, exp))
        candidates = filter(lambda x: \
                    sum(map(lambda y: y**exp, x)) == total, candidates)
    print part, dict((j, candidates[0].count(j)) for j in candidates[0])

#some computations:

###n=3
#[3] {3: 1}
#[2, 1] {0: 1, 1: 1}
#[1, 1, 1] {1: 1}

###n=4
#[4] {4: 1}
#[3, 1] {0: 1, 1: 1, 2: 1}
#[2, 2] {0: 1, 1: 1}
#[2, 1, 1] {0: 1, 1: 1, 2: 1}
#[1, 1, 1, 1] {0: 1}

###n=5  
#[5] {5: 1}
#[4, 1] {0: 1, 1: 1, 2: 1, 3: 1}
#[3, 2] {0: 2, 1: 2, 2: 1}
#[3, 1, 1] {0: 2, 1: 2, 2: 1, 3: 1}
#[2, 2, 1] {0: 2, 1: 2, 2: 1}
#[2, 1, 1, 1] {0: 2, 1: 2}
#[1, 1, 1, 1, 1] {1: 1}

###n=6
#[6] {6: 1}
#[5, 1] {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
#[4, 2] {0: 3, 1: 3, 2: 2, 3: 1}
#[4, 1, 1] {0: 3, 1: 3, 2: 2, 3: 1, 4: 1}
#[3, 3] {0: 2, 1: 2, 2: 1}
#[3, 2, 1] {0: 6, 1: 6, 2: 3, 3: 1}
#[3, 1, 1, 1] {0: 4, 1: 4, 2: 2}
#[2, 2, 2] {0: 2, 1: 2, 2: 1}
#[2, 2, 1, 1] {0: 4, 1: 4, 2: 1}
#[2, 1, 1, 1, 1] {0: 2, 1: 2, 2: 1}
#[1, 1, 1, 1, 1, 1] {0: 1}

###n=7
#[7] {7: 1}
#[6, 1] {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1}
#[5, 2] {0: 4, 1: 4, 2: 3, 3: 2, 4: 1}
#[5, 1, 1] {0: 4, 1: 4, 2: 3, 3: 2, 4: 1, 5: 1}
#[4, 3] {0: 5, 1: 5, 2: 3, 3: 1}
#[4, 2, 1] {0: 12, 1: 12, 2: 7, 3: 3, 4: 1}
#[4, 1, 1, 1] {0: 7, 1: 7, 2: 4, 3: 2}
#[3, 3, 1] {0: 8, 1: 8, 2: 4, 3: 1}
#[3, 2, 2] {0: 8, 1: 8, 2: 4, 3: 1}
#[3, 2, 1, 1] {0: 14, 1: 14, 2: 6, 3: 1}
#[3, 1, 1, 1, 1] {0: 6, 1: 6, 2: 2, 3: 1}
#[2, 2, 2, 1] {0: 6, 1: 6, 2: 2}
#[2, 2, 1, 1, 1] {0: 6, 1: 6, 2: 2}
#[2, 1, 1, 1, 1, 1] {0: 3, 1: 3}
#[1, 1, 1, 1, 1, 1, 1] {1: 1}

#spectral tally for the robinson schented correspondence
#does not quite work
d = {}
for i in glist:
    p = Permutation(i)
    #print p, p.number_of_fixed_points(), p.robinson_schensted()
    if p.left_tableau().shape() in d.keys():
        d[p.left_tableau().shape()].append(p.number_of_fixed_points())
    else:
        d[p.left_tableau().shape()] = [p.number_of_fixed_points()]
for i in d.keys():
    print i, dict((j, d[i].count(j)) for j in d[i])

#filter out which (left?) young symmetrizers are eigenvectors

from sage.combinat.symmetric_group_algebra import a as aaa
from sage.combinat.symmetric_group_algebra import b as bbb

def is_mult(x,y):
    flag = True
    for i in range(len(x)):
        if x[i] != 0:
            q = y[i]/x[i]
            flag = False
            break
    if flag:
        return True
    else:
        return reduce(lambda u,v: u and v, [y[i]==q*x[i] for i \
                      in range(len(x))])

b = beta(n,n).leftMultMatrix()
for p in Partitions(n):
    for t in p.standard_tableaux():
        x = aaa(t)*bbb(t)
        rcoeffs = []
        for i in glist:
            rcoeffs.append(x.coefficient(list(i.tuple())))
        
        res = (b*(Matrix(rcoeffs).transpose())).transpose()
        if is_mult(rcoeffs, res[0]):
            print p, t
            print rcoeffs
            print res.str()
            print ''

#cayley graph for beta_n
def cayleyGraph(n):
    b = beta(n,n).leftMultMatrix()
    dic = {}
    for i in range(fact(n)):
        dic[i] = [indi for indi, coefi in enumerate(b.rows()[i]) \
                  if coefi == 1]
    graph = DiGraph(dic)
    return graph