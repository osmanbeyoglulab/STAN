import numpy as np
from sklearn.base import RegressorMixin, BaseEstimator, clone


class ClusteredRegressor(RegressorMixin, BaseEstimator):
    def __init__(self, n_components=2, random_state=1, max_iter=100, tol=1e-10, verbose=False, K=None):
        self.n_components = n_components
        self.random_state = random_state
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose
        self.K=K

    def fit(self, X, y):

        np.random.seed(self.random_state)

        resid=y-X.dot(np.linalg.lstsq(X,y, rcond=None)[0])
        if self.K is not None:
            resid=self.K.dot(resid)/self.K.sum(axis=1, keepdims=True)


        self.class_weights=np.random.rand(len(y),self.n_components)


        for it in range(self.max_iter):
            self.beta = np.stack([np.linalg.lstsq(X*self.class_weights[:,i].reshape((-1,1)),y*self.class_weights[:,i].reshape((-1,1)),rcond=None)[0] for i in range(self.n_components)])
            resid=np.stack([y]*self.n_components)-X.dot(self.beta).transpose((1,0,2))
            resid=resid.squeeze().T
            if self.K is not None:
                resid=self.K.dot(resid.squeeze())/K.sum(axis=1, keepdims=True)

            #w=np.exp(-((resid**2).sum(axis=2)/(resid**2).sum(axis=2).mean())).T

            #self.class_weights=np.exp(-self.n_components*(resid**2)/((resid)**2).mean())
            #self.class_weights=self.class_weights/self.class_weights.sum(axis=1).reshape((-1,1))
            self.class_weights=(resid**2)==(resid**2).min(axis=1, keepdims=True)
        return self

    def score(self, X, y):
        y_fit=self.transform(X,y)
        return pearsonr(y.squeeze(),y_fit.squeeze())[0]

    def transform(self, X, y):
        class_weights=self.predict_proba(X, y, K=self.K)
        return self.predict(X, class_weights)

    def predict(self, X, class_weights):

        """ Calculate a matrix of conditional predictions """

        classes=np.argmax(class_weights, axis=1)

        y_pred=X.dot(self.beta)
        return y_pred[np.arange(y_pred.shape[0]), classes]

    def predict_proba(self, X, y):
        """ Estimate cluster probabilities of labeled data """
        resid=np.stack([y]*self.n_components)-X.dot(self.beta).transpose((1,0,2))
        resid=resid.squeeze().T
        if K is not None:
            resid=K.dot(resid.squeeze())/K.sum(axis=1, keepdims=True)

        #w=np.exp(-((resid**2).sum(axis=2)/(resid**2).sum(axis=2).mean())).T

        class_weights=np.exp(-(resid**2)/((resid)**2).mean())
        class_weights=class_weights/class_weights.sum(axis=1).reshape((-1,1))
        return class_weights

    def score_custers(self, X, y):
        classes=np.argmax(self.predict_proba(X,y), axis=1)
        pred=self.predict(X,self.predict_proba(X,y))


        return [pearsonr(np.asarray(y)[classes==i].squeeze(), pred[classes==i].squeeze()) for i in range(max(classes)+1)]
