#!/projects/opt/centos8/x86_64/miniconda3/py39_4.12.0/bin/python

import torch
import torchvision
from torchvision.io import read_image
from skimage.io import imread
import piq

def image_assessment(x, y):

    x = torch.tensor(imread(x)).permute(2, 0, 1)[None, ...]/255.
    y = torch.tensor(imread(y)).permute(2, 0, 1)[None, ...]/255.

    psnr_index = piq.psnr(x, y, data_range=1., reduction='none')
    print(f"PSNR index: {psnr_index.item():0.4f}")

    ssim_index = piq.ssim(x, y, data_range=1.)
    ssim_loss: torch.Tensor = piq.SSIMLoss(data_range=1.)(x, y)
    print(f"SSIM index: {ssim_index.item():0.4f}, loss: {ssim_loss.item():0.4f}")

    vsi_index: torch.Tensor = piq.vsi(x, y, data_range=1.)
    vsi_loss: torch.Tensor = piq.VSILoss(data_range=1.)(x, y)
    print(f"VSI index: {vsi_index.item():0.4f}, loss: {vsi_loss.item():0.4f}")

def main():
    img_err = ["1E-3", "3E-3", "5E-3", "7E-3", "9E-3"]
    
    for err in img_err:
        img_orig = "figures/NVB_C009_l10n512_S12345T692_z54.png"
        img_name = "figures/sz_rel__{}__NVB_C009_l10n512_S12345T692_z54.png".format(err)
        print(img_name)
        image_assessment(img_orig, img_name)

if __name__ == "__main__":
    main()
