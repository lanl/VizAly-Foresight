# python3 ./inference.py inputData modelpath outputPath 
# python3 ./inference.py /bigData/Turbulence/scalarHIT_fields100.h5 checkpoints /home/pascal/projects/VizAly-Foresight/Analysis/pat/turbulence/data/test1

#####################################################
# LOAD TRAINED MODEL AND MAKE PREDICTIONS
#####################################################


import numpy as np
import torch
from torch import nn
from tqdm import tqdm
import h5py
import torch.nn.functional as F
from torch.autograd import Variable
from torch.utils import data
from torch.utils.data import DataLoader
from turb_funcs import diagnostics_np
#from pytorch_utilsCFDNS import CFDNSdataset, PreprocessingCFDNS
import os
import sys


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
    

path = sys.argv[1]
print("path: ", path)

modelpath = sys.argv[2]
print("modelpath path: ", modelpath)

outputpath = sys.argv[3]
print("output path: ", outputpath)
os.mkdir(outputpath)



f = h5py.File(path,'r')
fields = f['fields']
nch = 5
data = fields[:50,:,:,:,:nch]
data = minmaxscaler(data) # scale data 0 - 1
seq_len = 3
epochs = 500000
batch_size = seq_len


# reshape data for pytorch input and conver to torch tensor
inp_tensor = convert_to_torchchannel(data)
print(inp_tensor.size())

model = conv3d_autoencoder()
optimizer =  torch.optim.Adam(model.parameters(), lr=1e-04,
                             weight_decay=1e-5)

# Encode data
print('Loading checkpoint to GPU...')
#modelpath = 'checkpoints'
checkpoint = torch.load(modelpath, map_location="cuda:0")
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
epoch = checkpoint['epoch']
loss = checkpoint['loss']
#print('loading model %s for evaluation, trained for %5d epochs with %8.8f final loss' % ( str(blkID), epoch, loss))
print('loading model')
model.cuda()
model.eval()


# save 2 diagnostics datasets

# 1
randidx = np.random.randint(0,10)
test_input = inp_tensor[randidx:randidx+seq_len,::].cuda()
print(test_input.size())
lat_space = model.encode(test_input.type(torch.cuda.FloatTensor))
print('lat space size', lat_space.size())
real_space = model.decode(lat_space.type(torch.cuda.FloatTensor))
print(real_space.size())

ytest = convert_to_numpychannel_fromtorch(test_input.detach())

predt = convert_to_numpychannel_fromtorch(real_space.detach())

#rescale to physical limits
dns = inverse_minmaxscaler(ytest,'rescale_coeffs_3DHIT.npy')
mod = inverse_minmaxscaler(predt,'rescale_coeffs_3DHIT.npy')

dns = dns[1,:,:,:,:3]
mod = mod[1,:,:,:,:3]
print(dns.shape)

#diagnostics_np(mod,dns,save_dir='diagnosticsCAE/', iteration=1, pos=[0,0,1], dx=[0.049, 0.049, 0.049], diagnostics=['spectrum', 'intermittency', 'structure_functions','QR'])
diagnostics_np(mod,dns,save_dir=(outputpath+'/'), iteration=1, pos=[0,0,1], dx=[0.049, 0.049, 0.049], diagnostics=['spectrum', 'intermittency', 'structure_functions','QR'])