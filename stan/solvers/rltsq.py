import pandas as pd
import numpy as np
from numpy.linalg import lstsq
def rlstsq(a,b,l, intercept=False, return_df=True):
    '''
    Solve the optimization problem:
        min_x  (1/n)*||a*x-b||_F^2 +l ||x||_F^2

    This is a wrapper function for scipy.linalg.lstsq

    Parameters
    ----------
    a : (M, N) array_like
        Left-hand side array
    b : (M, K) array_like
        Right-hand side array
    l : float
        regularization parameter

    Returns
    -------
    x : (N,K) ndarray
    least-squares solution
    '''
    if len(b.shape)==1:
        if type(b) is pd.core.series.Series:
            b=b.to_numpy()
        b=b.reshape(-1,1)

    A=np.concatenate((np.sqrt(1/a.shape[0])*a, np.sqrt(l)*np.eye(a.shape[1])))
    B=np.concatenate((np.sqrt(1/a.shape[0])*b, np.zeros((a.shape[1], b.shape[1]))))
    if intercept:
        A=np.concatenate((A, np.ones((A.shape[0], 1))))

    p, res, rnk, s = lstsq(A,B)

    if intercept:
        intercept_term=p[p.shape[0]-1,:]
        x=p[0:p.shape[0]-2,:]
    else:
        x=p

    if type(a) is pd.core.frame.DataFrame and type(b) is pd.core.frame.DataFrame and return_df:
        x=pd.DataFrame(data=x, index=a.columns, columns=b.columns)
        if intercept:
            intercept_term=pd.Series(data=intercept_term, index=b.columns)

    if intercept:
        return (x, intercept_term)
    else:
        return x