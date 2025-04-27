import cv2, os, timm, torch, torchvision
import numpy as np
import pandas as pd
# import ResNet as ResNet
import torch.nn as nn
from huggingface_hub import login, hf_hub_download
from torch.utils.data import DataLoader, TensorDataset
from torchvision import transforms
from transformers import AutoImageProcessor, AutoModelForImageClassification

from .bigan import BiGAN_train, BiGAN_encode
from .ctran import ctranspath

def get_img_features_pretrained(adata, method, is_model_exist = False, path = ''):
    if method == 'UNI': 
        transform = transforms.Compose([
            transforms.ToPILImage(),
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize(mean=[0.485, 0.456, 0.406],
                                 std=[0.229, 0.224, 0.225])
        ])
        img_data = preprocessing(adata, (112,112))

    
        if img_data.dtype != np.uint8:
            if img_data.max() <= 1.0:
                img_data = (img_data * 255).astype(np.uint8)
            else:
                img_data = img_data.astype(np.uint8)
        processed_patches = torch.stack([transform(patch) for patch in img_data])
        dataset = TensorDataset(processed_patches)
        data_loader = DataLoader(dataset, batch_size=1, shuffle=False)
        return extract_UNI_feature(data_loader, is_model_exist = is_model_exist, path = path)
    elif method == 'Transpath':
        from ctran import ctranspath
        transform = transforms.Compose([
            transforms.ToPILImage(),
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize(mean=[0.485, 0.456, 0.406],
                                 std=[0.229, 0.224, 0.225])
        ])
        img_data = preprocessing(adata, (112,112))
        if img_data.dtype != np.uint8:
            if img_data.max() <= 1.0:
                img_data = (img_data * 255).astype(np.uint8)
            else:
                img_data = img_data.astype(np.uint8)
        processed_patches = torch.stack([transform(patch) for patch in img_data])
        dataset = TensorDataset(processed_patches)
        data_loader = DataLoader(dataset, batch_size=1, shuffle=False)
        return extract_TP_faeture(data_loader, path = path)
    

def get_img_features_BiGAN(adata, path = ''):
    transform = transforms.Compose([
        transforms.ToPILImage(),
        transforms.ToTensor(),
        transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))
    ])
    img_data = preprocessing(adata, (36,36))
    processed_img_data = []
    for img in img_data:
        if img.dtype != np.uint8:
            img = (img * 255).astype(np.uint8)  
        processed_img_data.append(img)

    dataset = [transform(img) for img in processed_img_data]
    # Create a DataLoader
    batch_size = 8
    train_data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    

    test_loader = DataLoader(dataset, batch_size=64, shuffle=False, num_workers=4)
    
    
    return extract_BiGAN_faeture(train_data_loader, test_loader, path = path, sample_id  = list(adata.uns['spatial'].keys())[0] )

def preprocessing(adata,patch_size):
    # Convert images
    sample_id = list(adata.uns['spatial'].keys())[0] 

    img_hires = adata.uns['spatial'][sample_id]['images']['hires']
    print(img_hires.shape)
    # img_lowres = adata.uns['spatial'][sample_id]['images']['lowres']
    scale_factor = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']
    new_hir = img_hires.copy()
    coords = adata.obsm['spatial'] * scale_factor
    patches = []
    for i, spot in enumerate(coords):
        spot_x, spot_y = spot 
        spot_x = int(spot_x)
        spot_y = int(spot_y)
        patch = img_hires[spot_y:spot_y + patch_size[1], spot_x:spot_x + patch_size[0], :]
        
        if patch.shape[0] < patch_size[0] or patch.shape[1] < patch_size[1]:
            patch = cv2.copyMakeBorder(patch, 0, patch_size[0] - patch.shape[0], 0, patch_size[1] - patch.shape[1], cv2.BORDER_CONSTANT, value=0)
        patches.append(patch)
    
    # for i, entry in enumerate(patches[:10]):  # Display first 10 patches
    #     plt.subplot(2, 5, i + 1)
    #     plt.imshow(entry)
    #     plt.axis('off')
    # plt.show()
    img_data=np.array(patches)
    

    return img_data


def extract_UNI_feature(data_loader, is_model_exist = False, path = ''):
    all_features = []
    
    local_dir = path
    
    model = timm.create_model(
            "vit_large_patch16_224", img_size=224, patch_size=16, init_values=1e-5, num_classes=0, dynamic_img_size=True
        )
    if not is_model_exist:
        hf_hub_download("MahmoodLab/UNI", filename="pytorch_model.bin", local_dir=local_dir, force_download=True)
        model.load_state_dict(torch.load(os.path.join(local_dir, "pytorch_model.bin"), map_location="cpu"), strict=True)
        
    else:
        model.load_state_dict(torch.load(os.path.join(local_dir, "pytorch_model.bin"), map_location="cpu"), strict=True)    
    device = torch.device("cuda")
    model.to(device)
    model.eval()
    
    with torch.no_grad():
        for idx, batch in enumerate(data_loader):
            inputs = batch[0].to(device)  
            features = model(inputs)     
            features = features.cpu().numpy() 
            all_features.append(features)   
    
            if (idx + 1) % 1000 == 0:
                print(f"{idx + 1}")

    return np.squeeze(np.array(all_features), axis=1)

def extract_TP_faeture(data_loader, is_model_exist = False, path = ''):
    all_features = []
    
    model = ctranspath()
    model.head = nn.Identity()
    td = torch.load(path+'/'+'ctranspath.pth')
    model.load_state_dict(td['model'], strict=True)
    model.eval()
    device = torch.device("cuda")
    model = model.to(device)
    with torch.no_grad():
        for idx, batch in enumerate(data_loader):
            # 从 batch 中提取出张量
            inputs = batch[0].to(device)  # TensorDataset 返回的是一个元组，提取第一个元素
            inputs = inputs.to(device)
            features = model(inputs)       # 提取特征
            features = features.cpu().numpy()  # 移动到 CPU 并转换为 NumPy 数组
            all_features.append(features)      # 添加到列表中
    
            if (idx + 1) % 1000 == 0:
                print(f"{idx + 1}")
    return np.squeeze(np.array(all_features), axis=1)


def extract_BiGAN_faeture(train_data_loader,test_data_loader, path = '', sample_id=''):
    # Device configuration
    patch_size = (36,36)

    encoder = BiGAN_train(train_data_loader, patch_size,is_save = 'False', path = path, sample_id = sample_id)
    return BiGAN_encode(test_data_loader, encoder= encoder)

