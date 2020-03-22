"""
Debugging weird error in Pascal converted Raw data
"""

import numpy as np
import os
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
import os, sys
import pdb


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


# Gather inputs
filename1 = sys.argv[1]
filename2 = sys.argv[2]
outputdir = sys.argv[3]
snapshots = 50
numScalars = 5
dimX = 128
dimY = 128
dimZ = 128

print('Load first snap')
f1 = h5py.File(filename1,'r')
fields1 = f1['fields']
data1 = fields1[:snapshots,:,:,:,:numScalars]
#pdb.set_trace()

print('Load second snap')
f2 = h5py.File(filename2,'r')
fields2 = f2['fields']
data2 = fields2[:snapshots,:,:,:,:numScalars]

mod = data1[0,:,:,:,:3]
dns = data2[0,:,:,:,:3]

outputPath = outputdir + '/'
#os.mkdir("outputPath")
diagnostics_np(mod,dns,save_dir=outputPath, iteration=1, pos=[0,0,1], dx=[0.049, 0.049, 0.049], diagnostics=['spectrum', 'intermittency', 'structure_functions','QR'])

print('done')