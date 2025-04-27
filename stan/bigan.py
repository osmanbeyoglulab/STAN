# Author: Haoyu Wang

import os, torch
import torch.nn as nn
import torch.optim as optim
import torchvision.models as models
import torchvision.utils as vutils

from torch.nn.utils import spectral_norm
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms

# Weight initialization function
def weights_init(m):
    classname = m.__class__.__name__
    if classname.find('Conv') != -1 or classname.find('Linear') != -1:
        nn.init.kaiming_normal_(m.weight, a=0.2)
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)
    elif classname.find('BatchNorm') != -1:
        nn.init.constant_(m.weight, 1)
        nn.init.constant_(m.bias, 0)

# Encoder (E)
class Encoder(nn.Module):
    def __init__(self, latent_dim=128):
        super(Encoder, self).__init__()
        self.model = nn.Sequential(
            # Input: 3 x 26 x 26
            nn.Conv2d(3, 32, 4, 2, 1),    # Output: 32 x 13 x 13
            nn.BatchNorm2d(32),
            nn.LeakyReLU(0.2),
            nn.Conv2d(32, 64, 4, 2, 1),   # Output: 64 x 6 x 6
            nn.BatchNorm2d(64),
            nn.LeakyReLU(0.2),
            nn.Conv2d(64, 128, 4, 2, 1),  # Output: 128 x 3 x 3
            nn.BatchNorm2d(128),
            nn.LeakyReLU(0.2),
            nn.Flatten(),
            nn.Linear(128 * 3 * 3, latent_dim)
        )

    def forward(self, x):
        return self.model(x)

# Generator (G)
class Generator(nn.Module):
    def __init__(self, latent_dim=128):
        super(Generator, self).__init__()

        self.model = nn.Sequential(
            nn.Linear(latent_dim, 128 * 3 * 3),
            nn.LeakyReLU(0.2),
            nn.Unflatten(1, (128, 3, 3)),
            nn.BatchNorm2d(128),
            # First deconvolution
            nn.ConvTranspose2d(128, 64, 4, 2, 1),  # Output: 64 x 6 x 6
            nn.BatchNorm2d(64),
            nn.LeakyReLU(0.2),
            # Second deconvolution
            nn.ConvTranspose2d(64, 32, 4, 2, 1),   # Output: 32 x 12 x 12
            nn.BatchNorm2d(32),
            nn.LeakyReLU(0.2),
            # Third deconvolution
            nn.ConvTranspose2d(32, 16, 4, 2, 1),   # Output: 16 x 24 x 24
            nn.BatchNorm2d(16),
            nn.LeakyReLU(0.2),
            nn.Upsample(size=(36, 36), mode='bilinear', align_corners=False),
            nn.Conv2d(16, 3, 3, 1, 1),             # Output: 3 x 36 x 36
            nn.Tanh()
        )

    def forward(self, z):
        return self.model(z)

# Discriminator (D)
class Discriminator(nn.Module):
    def __init__(self):
        super(Discriminator, self).__init__()
        self.net = nn.Sequential(
            spectral_norm(nn.Conv2d(3, 32, 4, 2, 1)),
            nn.BatchNorm2d(32),
            nn.LeakyReLU(0.2, inplace=True),
            
            spectral_norm(nn.Conv2d(32, 64, 4, 2, 1)),     # Output: 64 x 9 x 9
            nn.BatchNorm2d(64),
            nn.LeakyReLU(0.2, inplace=True),
            
            spectral_norm(nn.Conv2d(64, 128, 4, 2, 1)),    # Output: 128 x 4 x 4
            nn.BatchNorm2d(128),
            nn.LeakyReLU(0.2, inplace=True),
            
            spectral_norm(nn.Conv2d(128, 256, 4, 1, 0)),   # Output: 256 x 1 x 1
            nn.BatchNorm2d(256),
            nn.LeakyReLU(0.2, inplace=True),
            
            nn.Flatten(),
            nn.Linear(256, 1)            
        )
    
    def forward(self, x):
        validity = self.net(x)
        return validity
        
class PerceptualLoss(nn.Module):
    def __init__(self):
        super(PerceptualLoss, self).__init__()
        self.vgg = models.vgg16(pretrained=True).features[:16].eval()
        for param in self.vgg.parameters():
            param.requires_grad = False

    def forward(self, x, y):
        return nn.functional.l1_loss(self.vgg(x), self.vgg(y))

def BiGAN_train(train_data_loader, patch_size, is_save = 'False', path = '',sample_id ='' ):
    # Device configuration
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # Initialize models
    latent_dim = 128
    encoder = Encoder(latent_dim).to(device)
    generator = Generator(latent_dim).to(device)
    discriminator = Discriminator().to(device)
    perceptual_loss = PerceptualLoss().to(device)
    
    # Apply weight initialization
    encoder.apply(weights_init)
    generator.apply(weights_init)
    discriminator.apply(weights_init)
    adversarial_loss = nn.BCEWithLogitsLoss()
    reconstruction_loss_function = nn.L1Loss()
    lambda_recon = 20  # Recon_loss_weight
    
    optimizer_D = optim.Adam(discriminator.parameters(), lr=1e-4, betas=(0.5, 0.8))
    optimizer_G = optim.Adam(generator.parameters(), lr=1e-3, betas=(0.5, 0.85))
    optimizer_E = optim.Adam(encoder.parameters(), lr=1e-3, betas=(0.5, 0.85))
    
    num_epochs = 25
    real_label = 0.9
    fake_label = 0.1
    
    os.makedirs('real_fake', exist_ok=True)
    
    for epoch in range(num_epochs):
        encoder.train()
        generator.train()
        discriminator.train()
        for batch_idx, batch in enumerate(train_data_loader):
            batch = batch.to(device)
            current_batch_size = batch.size(0)
    
            label_real = torch.full(
                (current_batch_size, 1), real_label, device=device, dtype=torch.float
            )
            label_fake = torch.full(
                (current_batch_size, 1), fake_label, device=device, dtype=torch.float
            )
    
            center_crop = batch[:, :, 5:31, 5:31]  #36-26
    
            optimizer_D.zero_grad()
    
            # Real
            validity_real = discriminator(batch)
            loss_D_real = adversarial_loss(validity_real, label_real)
    
            # Generate
            with torch.no_grad():
                z_fake = encoder(center_crop)
                x_fake = generator(z_fake)
            validity_fake = discriminator(x_fake.detach())
            loss_D_fake = adversarial_loss(validity_fake, label_fake)
    
            # Loss for D
            loss_D = (loss_D_real + loss_D_fake) / 2
            loss_D.backward()
    
            torch.nn.utils.clip_grad_norm_(discriminator.parameters(), max_norm=1.0)
    
            optimizer_D.step()
    
            # Train E and G
            optimizer_E.zero_grad()
            optimizer_G.zero_grad()
    
            z_fake = encoder(center_crop)
            x_fake = generator(z_fake)
            validity = discriminator(x_fake)
            adversarial_loss_value = adversarial_loss(validity, label_real)
    
            # recon_loss
            reconstruction_loss_value = reconstruction_loss_function(x_fake, batch)
    
            # total loss
            loss_GE = adversarial_loss_value + lambda_recon * reconstruction_loss_value + perceptual_loss(x_fake, batch)
    
            loss_GE.backward()
    
            torch.nn.utils.clip_grad_norm_(generator.parameters(), max_norm=3)
            torch.nn.utils.clip_grad_norm_(encoder.parameters(), max_norm=3)
    
            optimizer_E.step()
            optimizer_G.step()
    
            if batch_idx % 100 == 0:
                with torch.no_grad():
                    batch = batch.to(device)
    
                    # emb
                    z_fake = encoder(center_crop).detach()
                    # img
                    fake_images = generator(z_fake).cpu()
                    real_images = batch.cpu()
    
                    # create grid
                    real_grid = vutils.make_grid(real_images, nrow=8, normalize=True, padding=2)
                    fake_grid = vutils.make_grid(fake_images, nrow=8, normalize=True, padding=2)
    
                    # cat
                    combined_grid = torch.cat((real_grid, fake_grid), dim=1)  # dim=1 
    
                    # save
                    vutils.save_image(combined_grid, f'real_fake/real_fake_epoch_{epoch+1}_batch_{batch_idx}.png')
    
        
                    print(
                    f'Epoch [{epoch+1}/{num_epochs}], Batch [{batch_idx}/{len(train_data_loader)}], '
                    f'Loss D: {loss_D.item():.4f}, Loss GE: {loss_GE.item():.4f}'
                )
                with torch.no_grad():
                    print(
                        f'Validity Real Mean: {validity_real.mean().item():.4f}, '
                        f'Std: {validity_real.std().item():.4f}'
                    )
                    print(
                        f'Validity Fake Mean: {validity_fake.mean().item():.4f}, '
                        f'Std: {validity_fake.std().item():.4f}'
                    )
                    print(
                        f'Validity GE Mean: {validity.mean().item():.4f}, '
                        f'Std: {validity.std().item():.4f}'
                    )
                    print(
                        f'Reconstruction Loss: {reconstruction_loss_value.item():.4f}'
                    )
    if is_save:
        if path =='':
            os.makedirs('./results', exist_ok=True)
            torch.save(encoder.state_dict(), f'./results/{sample_id}_{patch_size[0]}_{num_epochs}l_BIGAN_encoder_new.pth')
            torch.save(generator.state_dict(), f'./results/{sample_id}_{patch_size[0]}_{num_epochs}l_BIGAN_generator.pth')
            torch.save(discriminator.state_dict(), f'./results/{sample_id}_{patch_size[0]}_{num_epochs}l_BIGAN_discriminator.pth')
        else:
            torch.save(encoder.state_dict(), f'{path}/{sample_id}_{patch_size[0]}_{num_epochs}l_BIGAN_encoder_new.pth')
            torch.save(generator.state_dict(), f'{path}/{sample_id}_{patch_size[0]}_{num_epochs}l_BIGAN_generator.pth')
            torch.save(discriminator.state_dict(), f'{path}/{sample_id}_{patch_size[0]}_{num_epochs}l_BIGAN_discriminator.pth')

    return encoder

def BiGAN_encode(test_data_loader,encoder):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    load_encoder = encoder

    load_encoder.eval()
    
    results = []
    os.makedirs('encode_results', exist_ok=True)
    
    with torch.no_grad():
        for batch_idx, batch in enumerate(test_data_loader):
            batch = batch.to(device)
            current_batch_size = batch.size(0)
            center_crop = batch[:, :, 5:31, 5:31]  # 36-26
            
            z_fake = encoder(center_crop)
            results.append(z_fake.cpu())
            
            # # x_fake = generator(z_fake).cpu()
            
            # # recon_loss, to check the result
            # reconstruction_loss = reconstruction_loss_function(x_fake, batch.cpu())
            # print(f'Test Batch [{batch_idx}], Reconstruction Loss: {reconstruction_loss.item():.4f}')
    
            # real_images = batch.cpu()
            # fake_images = x_fake
      
            # real_grid = vutils.make_grid(real_images, nrow=8, normalize=True, padding=2)
            # fake_grid = vutils.make_grid(fake_images, nrow=8, normalize=True, padding=2)
            
            # combined_grid = torch.cat((real_grid, fake_grid), dim=1)  # dim=1 
    
            # vutils.save_image(combined_grid, f'encode_results/test_real_fake_batch_{batch_idx}.png')
            
            # print(f'Encode Batch [{batch_idx}] processed.')
    
    print("Encoding Completed.")
    return torch.cat(results, dim=0).detach().cpu().numpy()
        
        
        
