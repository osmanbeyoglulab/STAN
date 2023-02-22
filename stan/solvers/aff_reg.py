import numpy as np

def aff_reg(a,b,c,l1=0, l2=1, svda=None, svdb=None):
    (d1, d2)=a.shape
    (d3, d4)=b.shape


    if svda is None:
        [ua,sa,va]=np.linalg.svd(a, full_matrices=False)
    else:
        [ua,sa,va]=svda
    if svdb is None:
        [ub,sb,vb]=np.linalg.svd(b, full_matrices=False)
    else:
        [ub,sb,vb]=svdb

    scale=np.multiply(sa.reshape((len(sa), 1)), sb.reshape((1,len(sb)))) **2
    scale=np.divide(1, scale+(d1*d4*l2))
    x=va.T.dot(np.multiply(scale, va.dot(a.T.dot(c).dot(b.T)).dot(ub))).dot(ub.T)
    if l1>0:
        x=fista(a,b,c, (d1*d4*l1), (d1*d4*l2), x0=x, L=(max(sa)**2) * (max(sb)**2), m= (min(sa)**2) * (min(sb)**2))

    return x

def fista(a,b,c,l1, l2, x0=None, max_iter=50,L=1,m=0.1):

    ata=a.T.dot(a)
    bbt=b.dot(b.T)
    atcbt=a.T.dot(c).dot(b.T)

    theta=(1-np.sqrt(m/L))/(1+np.sqrt(m/L))
    t=1/L

    t=[]
    xold=x0
    xoldold=x0
    for k in range(max_iter):
        y=xold+theta*(xold-xoldold)
        grad=ata.dot(y).dot(bbt)-atcbt-l2*y
        yminusgrad=y-t*grad
        x=np.sign(x)*np.max(0, abs(x)-(l1/t))
        xoldold=xold
        xold=x
        if np.linalg.norm(xoldold-xold)/np.linalg.norm(xold):
            break
    return x





class AffReg():
    def __init__(self, l1=0, l2=1):
        self.l1=l1
        self.l2=l2

    def fit(self, a, b, c):
        pass
    def perm_test(self, n_trials=100):
        pass