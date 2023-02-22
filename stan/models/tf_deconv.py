import numpy as np
from sklearn.neighbors import KNeighborsClassifier


class ClusteredRegressor():
    def __init__(self, n_components=2, random_state=1, max_iter=100, tol=1e-10, verbose=False, K=None):
        self.n_components = n_components
        self.random_state = random_state
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose
        if K is None:
            K= #TODO: Make kernel here...
        self.K=K


        self.clusters=None
        self.knn=KNeighborsClassifier(n_neighbors=6)

    def fit(self, C, X, y):
        self.X=X
        self.y=y

        np.random.seed(self.random_state)

        resid=y-X.dot(np.linalg.lstsq(X,y, rcond=None)[0])
        # if self.K is not None:
        #     resid=self.K.dot(resid)/self.K.sum(axis=1, keepdims=True)


        self.class_weights=np.random.rand(len(y),self.n_components)


        for it in range(self.max_iter):
            self.beta = np.stack([np.linalg.lstsq(X*self.class_weights[:,i].reshape((-1,1)),y*self.class_weights[:,i].reshape((-1,1)),rcond=None)[0] for i in range(self.n_components)])
            resid=np.stack([y]*self.n_components)-X.dot(self.beta).transpose((1,0,2))
            resid=resid.squeeze().T
            # if self.K is not None:
            #     resid=self.K.dot(resid.squeeze())/K.sum(axis=1, keepdims=True)
            self.class_weights=(resid**2)==(resid**2).min(axis=1, keepdims=True)

        self.classes=np.argmax(class_weights, axis=1)
        self.knn.fit(self.C, self.classes)

        return self

    def predict(self,C, X, classes=None):

        """ Calculate a matrix of conditional predictions """
        if class_weights is not None:
            classes=np.argmax(class_weights, axis=1)

        else:
            classes=predict_class(C)

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

    def predict_class(self,C):
        return self.knn.predict(C)

    def score(self, C, X, y, classes=None):
        y_fit=self.predict(C, X, classes=classes)
        return pearsonr(y.squeeze(),y_fit.squeeze())[0]

    def score_custers(self, X, y):
        classes=np.argmax(self.predict_proba(X,y), axis=1)
        pred=self.predict(X,self.predict_proba(X,y))

        return [pearsonr(np.asarray(y)[classes==i].squeeze(), pred[classes==i].squeeze()) for i in range(max(classes)+1)]

    