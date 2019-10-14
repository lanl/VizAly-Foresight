"""
Train latent space of CAE for temporal modeling
"""
import numpy as np
import os, copy
import h5py
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import pywt
import argparse


def minmaxscaler(data):
    """ scale large turbulence dataset by channel"""
    nsnaps = data.shape[0]
    dim = data.shape[1]
    nch = data.shape[4]
    
    #scale per channel
    data_scaled = []
    rescale_coeffs = []
    for i in range(nch):
        data_ch = data[:,:,:,:,i]
        minval = data_ch.min(axis=0)
        maxval = data_ch.max(axis=0)
        temp = (data_ch - minval)/(maxval - minval)
        data_scaled.append(temp)
        rescale_coeffs.append((minval,maxval))
    data_scaled = np.stack(data_scaled, axis=4)
    np.save('rescale_coeffs_3DHIT', rescale_coeffs)
    return data_scaled

def inverse_minmaxscaler(data,filename):
    """ Invert scaling using previously saved minmax coefficients """
    rescale_coeffs = np.load(filename)
    nsnaps = data.shape[0]
    dim = data.shape[1]
    nch = data.shape[4]
    
    #scale per channel
    data_orig = []
    for i in range(nch):
        data_ch = data[:,:,:,:,i]
        (minval, maxval) = rescale_coeffs[i]
        temp = data_ch*(maxval - minval) + minval
        data_orig.append(temp)
    data_orig = np.stack(data_orig, axis=4)
    return data_orig

def convert_to_torchchannel(data):
    """ converts from  [snaps,dim1,dim2,dim3,nch] ndarray to [snaps,nch,dim1,dim2,dim3] torch tensor"""
    nsnaps = data.shape[0]
    dim1, dim2, dim3 = data.shape[1], data.shape[2], data.shape[3] 
    nch = data.shape[-1] #nch is last dimension in numpy input
    torch_permuted = np.zeros((nsnaps, nch, dim1, dim2, dim3))
    for i in range(nch):
        torch_permuted[:,i,:,:,:] = data[:,:,:,:,i]
    torch_permuted = torch.from_numpy(torch_permuted)
    return torch_permuted


def convert_to_numpychannel_fromtorch(tensor):
    """ converts from [snaps,nch,dim1,dim2,dim3] torch tensor to [snaps,dim1,dim2,dim3,nch] ndarray """
    nsnaps = tensor.size(0)
    dim1, dim2, dim3 = tensor.size(2), tensor.size(3), tensor.size(4)
    nch = tensor.size(1)
    numpy_permuted = torch.zeros(nsnaps, dim1, dim2, dim3, nch)
    for i in range(nch):
        numpy_permuted[:,:,:,:,i] = tensor[:,i,:,:,:]
    numpy_permuted = numpy_permuted.numpy()
    return numpy_permuted    


class conv3d_autoencoder(nn.Module):
    """ Conv. autoencoder for turbulent datasets, to input to ConvLSTM """
    def __init__(self):
        super(conv3d_autoencoder, self).__init__()
        self.nfilters = 25
        self.nch = 5
        self.encoder = nn.Sequential(
            nn.Conv3d(out_channels=self.nfilters, in_channels= self.nch, kernel_size=3, stride=2),
            nn.ReLU(),
            nn.Conv3d(out_channels=self.nfilters, in_channels=self.nfilters, kernel_size=3, stride=2),
            nn.ReLU(),
            nn.Conv3d(out_channels=self.nfilters, in_channels=self.nfilters, kernel_size=3, stride=2),
            nn.ReLU()        
        )
        self.decoder = nn.Sequential(
            nn.ConvTranspose3d(out_channels=self.nfilters ,in_channels=self.nfilters, kernel_size=3, stride=2),
            nn.ReLU(),
            nn.ConvTranspose3d(out_channels=self.nfilters ,in_channels=self.nfilters, kernel_size=3, stride=2),
            nn.ReLU(),
            nn.ConvTranspose3d(out_channels=self.nch ,in_channels=self.nfilters, kernel_size=3, stride=2, output_padding=1),
            nn.ReLU()
        )

    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        return x
    
    def encode(self, x):
        x = self.encoder(x)
        return x
        
    def decode(self, x):
        x = self.decoder(x)
        return x
    

def compress_data(data, seq_len, model):
    """ Pass data thru pre-trained conv. autoencoder"""
    nfilters  = 25
    latent_size = 15
    print('converting entire dataset into latent space...')
    datalen = data.size(0)
    #data = torch.from_numpy(data)
    nbatches = int(datalen/seq_len)
    lat_space = torch.zeros(nbatches,seq_len,nfilters,latent_size,latent_size,latent_size)
    for i in tqdm(range(nbatches)):
        cur_input = data[i:i+seq_len,::]
        #cur_input = cur_input.permute(0,4,1,2,3)
        #print(cur_input.size())
        lat_space[i,::] = model.encode(cur_input.type(torch.FloatTensor))

    return lat_space
    

def decompress_data(data, seq_len, model):
    nfilters = 5
    real_size = 128
    nbatches = data.size(0)
    real_space = torch.zeros(nbatches, seq_len, nfilters, real_size, real_size, real_size)
    
    for i in tqdm(range(nbatches)):
        real_space[i,::] = model.decode(data[i,::].type(torch.FloatTensor))
    
    return real_space

    
# Define dataset
#########################################################

path = '/home/arvindm/datasets/ScalarHIT/128cube/scalarHIT_fields100.h5'
f = h5py.File(path,'r')
fields = f['fields']
nch = 5
seq_len = 3
data = fields[:50,:,:,:,:nch] #all 5 fields including passive scalars
data = minmaxscaler(data) # scale data 0 - 1
predictFlag = True #load trained model and predict
resumeTrainFlag = True

# Create Latent space representation
##########################################################
print('Loading checkpoint to GPU...')
model = conv3d_autoencoder()
modelpath = '/home/arvindm/MELT/JoTruns/CCLSTM/CAE/scalarHIT_CAE_k3/l3/checkpoints'
checkpoint = torch.load(modelpath, map_location="cuda:0")
model.load_state_dict(checkpoint['model_state_dict'])
data2 = convert_to_torchchannel(data)

# COMPRESS
lat_space = compress_data(data2, seq_len, model)
lat_space = lat_space.detach() # very important, input dataset should not have autograd history!
print('Latent space full data size is',lat_space.size()) #size is (nbatches,seq_len,nfilters,dim,dim,dim) where seq_len in CSLTM = batchsize in CAE
print('Writing compressed data to file...')
f = h5py.File('compressedData.h5','w')
f.create_dataset('latentSpace',data=lat_space)
f.close()

#DECOMPRESS
real_space = decompress_data(lat_space, seq_len, model)
print('real space size is', real_space.size())

print('Writing decompressed data to file...')
g = h5py.File('decompressedData.h5','w')
g.create_dataset('decodedSpace',data=real_space.detach())
g.close()
