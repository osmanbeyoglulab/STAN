#
# import os
# import numpy as np
# import pandas as pd
# import scanpy
# from PIL import Image
#
# from tensorflow.keras.applications.resnet50 import ResNet50, preprocess_input as preprocess_resnet, decode_predictions
# from tensorflow.keras import backend as K
# from tensorflow.keras.layers import Dense, GlobalAveragePooling2D, Input, concatenate, Dropout, Lambda
# from tensorflow.keras.applications.inception_v3 import InceptionV3, preprocess_input as preprocess_inception
# from tensorflow.keras.applications.xception import Xception, preprocess_input as preprocess_xception
# from tensorflow.keras.losses import mse, binary_crossentropy
# from tensorflow.keras.models import Model, Sequential
# from tensorflow.keras.optimizers import Adam
#
#
# from tqdm import tqdm
#
# class ResNet:
#     __name__ = "ResNet"
#
#     def __init__(self, batch_size=1):
#         self.model = ResNet50(include_top=False, weights='imagenet', pooling="avg")
#         self.batch_size = batch_size
#         self.data_format = K.image_data_format()
#
#     def predict(self, x):
#         if self.data_format == "channels_first":
#             x = x.transpose(0, 3, 1, 2)
#         x = preprocess_resnet(x.astype(K.floatx()))
#         return self.model.predict(x, batch_size=self.batch_size)
#
# class Inception_V3:
#     """
#     pre-trained Inception_V3 model
#     """
#     __name__ = "Inception_V3"
#
#     def __init__(self, batch_size=1):
#         self.model = InceptionV3(include_top=False, weights='imagenet', pooling="avg")
#         self.batch_size = batch_size
#         self.data_format = K.image_data_format()
#
#     def predict(self, x):
#         if self.data_format == "channels_first":
#             x = x.transpose(0, 3, 1, 2)
#         x = preprocess_inception(x.astype(K.floatx()))
#         return self.model.predict(x, batch_size=self.batch_size)
#
#
# class Xception_imagenet:
#     """
#     pre-trained xception model
#     """
#     __name__ = "xception"
#
#     def __init__(self, batch_size=1):
#         self.model = Xception(include_top=False, weights='imagenet', pooling="avg")
#         self.batch_size = batch_size
#         self.data_format = K.image_data_format()
#
#     def predict(self, x):
#         if self.data_format == "channels_first":
#             x = x.transpose(0, 3, 1, 2)
#         x = preprocess_xception(x.astype(K.floatx()))
#         return self.model.predict(x, batch_size=self.batch_size)
#
#
# def cnn_encode(tiles, model):
#     features = model.predict(tiles)
#     #features = features.ravel()
#     return features
#
# def image_encoder(adata, image, scale=1, window_size=100, method="hist", n_feats= None, n_bins=32):
#     #TODO: this gives a bug when method="hist" and n_feats=8 or 4.  I think it has to do with all of the pixels being the same color.
#     # update it might have just fixed itself?? i cant reproduce the error!
#
#
#     if method=="hist":
#         image_feature_df=pd.DataFrame(None, columns=[str(i) for i in range(0,n_bins)], index=adata.obs_names)
#
#         im = image.convert("P", palette = Image.ADAPTIVE, colors = n_bins)
#         for i in tqdm(range(len(adata.obs_names))):
#             x=round(adata.obsm['spatial'][i,1]*scale)
#             y=round(adata.obsm['spatial'][i,0]*scale)
#             counts=np.array([c[0] for c in im.crop((y-window_size,x-window_size,y+window_size ,x+window_size)).getcolors()])
#             colors=[c[1] for c in im.crop((y-window_size,x-window_size,y+window_size ,x+window_size)).getcolors()]
#             feats=np.zeros((n_bins,))
#             feats[colors]=counts
#             image_feature_df.iloc[i,:]=feats
#
#     elif method=="cnn":
#         image_feature_df=pd.DataFrame(None, columns=[str(i) for i in range(0,2048)], index=adata.obs_names)
#
#         #import tensorflow only if we use the CNN.
#
#         model=Xception_imagenet()
#         im=np.asarray(image)
#         for i in tqdm(range(len(adata.obs_names))):
#             x=round(adata.obsm['spatial'][i,1]*scale)
#             y=round(adata.obsm['spatial'][i,0]*scale)
#
#             tile=im[x-window_size:x+window_size, y-window_size:y+window_size,:]
#             tiles=np.zeros((1,2*window_size,2*window_size,3))
#             tiles[0,:,:,:]=tile
#             feats=cnn_encode(tiles, model)
#             image_feature_df.iloc[i,:]=feats
#
#     if n_feats is not None:
#
#         #todo: dimensionality reduction here
#         pass
#     adata.obsm['hist_feats']=image_feature_df.to_numpy()
#
#     return image_feature_df
