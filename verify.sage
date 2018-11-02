import time,sys,os
print version()
load('generators.sage')
print "importing functions ..."


def generate_representation(partition, gen = GEN, extra = False, lookup = True, save = True):
    """
    the Sk representation corresponding to the given partition 
    using the generators (12),(12...k) in gen[partition] 
    """
    
    # first look it up
    filename = "rep-%s" % "".join(map(str,partition))
    if lookup and os.path.isfile(filename):
        return loads(open(filename, 'rb').read())

    # set generators
    k = sum(partition)
    swap, rot = GEN[partition]
    assert swap.is_square() and swap.dimensions() == rot.dimensions()
    rep = {Permutation([2,1] + range(3,k+1)) : matrix(SR,swap).canonicalize_radical(), 
           Permutation(range(2,k+1) + [1]) : matrix(SR,rot).canonicalize_radical()}

    # generate representation
    while len(rep) < factorial(k):
        rk = rep.keys()
        for x in rk:
            for y in rk:
                if x*y not in rep:
                    rep[x*y] = (rep[x] * rep[y]).canonicalize_radical()

    # check Sk relations
    s = [rep[Permutation([(i,i+1),(k,)])] for i in range(1,k-1)] + [rep[Permutation([(k-1,k)])]]
    for i in range(len(s)):
        assert (s[i]^2).is_one()
    for i in range(1,len(s)):
        for j in range(i-1):
            assert s[i]*s[j] == s[j]*s[i]
    for i in range(len(s)-1):
        assert s[i]*s[i+1]*s[i] == s[i+1]*s[i]*s[i+1]
    
    # a more extensive test
    if extra:   
        assert rep[Permutation(range(1,k+1))].is_one()
        for x in rep.keys():
            for y in rep.keys():
                assert rep[x]*rep[y] == rep[x*y]
                
    # done
    if save and not os.path.isfile(filename):
        open(filename, 'wb').write(dumps(rep))
    return rep


def profile_covariance(k, n, lookup = True, save = True):
    """
    the COV matrix of the random k-profile in Sn.
    note: we don't normalize by binomial(n,k)
    """
    
    # first look it up
    filename = "cov-%d-%d" % (k,n)
    if lookup and os.path.isfile(filename):
        return loads(open(filename, 'rb').read())

    # straightforward computation
    cov = matrix(ZZ, factorial(k), factorial(k))
    lex = {p : i for i,p in enumerate(Permutations(k))}
    for P in Permutations(n):
        W = Word(P)
        profile = vector(ZZ, factorial(k))
        for w in Subwords(W,k):
            profile[lex[w.standard_permutation()]] += 1
        cov += profile.outer_product(profile)
    cov /= factorial(n)
    
    # done
    if save and not os.path.isfile(filename):
        open(filename, 'wb').write(dumps(cov))
    return cov


def verbose(st):
    print time.strftime('[%a,%H:%M:%S]'),
    print st
    sys.stdout.flush()

    
def verify_diagonalization(k, steps):
    """
    For the k-profile, test whether the n^(k+s) term is diagonalized
    by the suitable matrix elements of Sk-representations 
    """
    
    # loop computes covariance's coefs of binomial(n,k+s) 
    coefs = [matrix(factorial(k), factorial(k), 1/factorial(k))]
    for s in range(1,steps+1):
        verbose('============ Component %d of %d-profile ============' % (k-s, k))
        
        # find representation
        R = []
        for partition in GEN:
            if sum(partition) == k and partition[0] == s:
                verbose('Finding representation %s ...' % str(partition))
                R.append(generate_representation(partition))
        U = matrix([[rep[p][i,j] for p in Permutations(k)] for rep in R 
                    for i in range(rep.values()[0].dimensions()[0]) 
                    for j in range(rep.values()[0].dimensions()[0])])                    
        
        # find covariance, and subtract lower orders
        verbose("Finding covariance of %d-profile in S%d ..." % (k,k+s))        
        cov = profile_covariance(k, k+s)
        cov -= sum([binomial(k+s, k+j) * coef for j,coef in enumerate(coefs)])
        coefs.append(cov)
        
        # check if projection on representations is diagonal
        verbose("Left projection on %d matrix elements ..." % U.dimensions()[0])
        test = (U.canonicalize_radical() * cov).canonicalize_radical()
        verbose("Right projection on %d matrix elements ..." % U.dimensions()[0])
        test = (test * U.transpose().canonicalize_radical()).canonicalize_radical()
        verbose("Testing restricted covariance coefficients ...")
        if test.nonzero_positions() == [(i,i) for i in range(test.dimensions()[0])]:
            verbose("Order n^%d term is diagonal of full rank!" % (k+s))
        else:
            verbose("There's a problem! \n%s" % str(test.nonzero_positions()))
            return False

    # all orders were diagonal
    return True

    
#======================================================================================
#   optional, distributed computation to avoid long wait

def profile_covariance_part(k, n, prefix):
    cov = matrix(ZZ, factorial(k), factorial(k))
    suffix = [x for x in range(1,n+1) if x not in prefix]
    lex = {p : i for i,p in enumerate(Permutations(k))}
    for P in Permutations(suffix):
        W = Word(prefix) + Word(P)
        profile = vector(ZZ, factorial(k))
        for w in Subwords(W,k):
            profile[lex[w.standard_permutation()]] += 1
        cov += profile.outer_product(profile)
    cov /= factorial(n)
    filename = "cov-%d-%d-%s" % (k,n,'-'.join(map(str,prefix)))
    open(filename, 'wb').write(dumps(cov))
    return cov

def profile_covariance_merge(k, n, pre = 1):
    cov = matrix(ZZ, factorial(k), factorial(k))
    for prefix in Permutations(n,pre):
        filename = "cov-%d-%d-%s" % (k,n,'-'.join(map(str,prefix)))
        cov += loads(open(filename, 'rb').read())
    filename = "cov-%d-%d-merged" % (k,n)
    open(filename, 'wb').write(dumps(cov))
    return cov
