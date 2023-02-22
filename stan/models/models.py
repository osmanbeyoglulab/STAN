import pandas as pd

from .base_model import *
class SpaTraFact_stl(ModelBase):
    def __init__(self, args, **kwargs):
        super().__init__(args, **kwargs)

    def train(self, fold=0, lambda1=1):
        #self.W=rlstsq(self.D_train, self.Y_train, lambda1)
        if self.svdD[fold] is None:
            self.svdD[fold]=np.linalg.svd(self.D_train[fold], full_matrices=False)
        (ua, sa, va)=self.svdD[fold]
        self.W[fold]=va.T.dot(np.diag(1/(lambda1*self.n_genes+sa**2)).dot(va)).dot(self.D_train[fold].T.dot(self.Y_train[fold]))
        #self.kernel=make_kernel(self.X, n_kernel_samples, std=length_scale, im_feats_weight=im_feature_weight), length_scale=3, n_kernel_samples=100, im_feature_weight=0.1, l=1)

class SpaTraFact(ModelBase):
    def __init__(self, args, **kwargs):
        """Initialize the model.
        Parameters
        ----------
        args: arguments for initializing the model.
        """
        super().__init__(args, **kwargs)

    def train(self, fold=0, n_kernel=500, ls=1, lw=5, lL=1, length_scale=2, im_features_weight=0.1):

        if self.svdD[fold] is None:
            self.svdD[fold]=np.linalg.svd(self.D_train[fold], full_matrices=False)
        if self.svdK is None:
            uk,sk,vk=np.linalg.svd(self.kernel, full_matrices=True)
            sk=np.concatenate((sk, [0]*(uk.shape[0]-len(sk))))
            self.svdK=(uk,sk,vk)

        [ud,sd,vd]=self.svdD[fold]
        [uk,sk,vk]=self.svdK

        D=self.D_train[fold]
        Y=self.Y_train[fold]

        scale=np.divide(1, sd.reshape((-1,1))**2+lw*(1-(lw*sk**2/(lw*sk**2+ls))).reshape(1,-1))
        W=vd.T.dot(np.multiply(scale, vd.dot(D.T.dot(Y)).dot(uk))).dot(uk.T)



        colums=self.tf_gene.columns.to_list()
        colums.append('intercept')
        #
        self.W[fold]=W
        #self.W_spat=pd.DataFrame(W.dot(uk.T.dot(uk)), index=colums, columns=self.Y.columns)
        #self.tf_activity=pd.DataFrame(self.W.T, columns=colums, index=self.Y.columns)

