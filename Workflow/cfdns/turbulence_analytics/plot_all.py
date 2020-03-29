
# python3 turbulence_analytics/plot_all.py /home/pascal/projects/VizAly-Foresight/Analysis/pat/turbulence/data/test2 original original_plot
import matplotlib
matplotlib.use('Agg')

import numpy as np
import tensorflow as tf

import h5py
import sys
import os

import matplotlib.pyplot as plt
#import seaborn as sns
plt.style.use('ggplot')
import glob

from turb_funcs import plot_contour

#import skimage.measure
#from scipy.stats import kde

total_iterations = [0]
#total_iterations = [0, 4, 8]
#total_iterations = [10, 40, 70]
num_iters = len(total_iterations)

plot_size_x = 2.5
plot_size_y = 7.5

latnet_color = 'blue'
true_color = 'green'


dataPath = sys.argv[1]
print("dataPath: ",dataPath)

outputPrefix = sys.argv[2]
print("outputPrefix: ", outputPrefix)
#diagnosticsCAE = "diagnosticsCAE"

outputPath = sys.argv[3]
print("outputPath: ", outputPath)
os.makedirs(outputPath, exist_ok = True)

###############################################
## SPECTRUM PLOTS
###############################################

def read_latnet_spectrum(iteration=[0]):
  print(iteration)
  files = [] 
  #for it in iteration:
    #files = files + glob.glob('diagnosticsCAE/iter_' + str(it).zfill(4) + '*energy_spectrum.npz')

  for file in glob.glob( dataPath + "/" + "*energy_spectrum.npz"):
    files.append(file)
  print("energy_spectrum: ", files)

  latnet_e = []
  true_e = []
  for f in files:
    latnet_spectrum = np.load(f) 
    latnet_e.append(latnet_spectrum['net_e'])
    true_e.append(latnet_spectrum['true_e'])
  latnet_e = np.mean(latnet_e, axis=0) 
  true_e = np.mean(true_e, axis=0) 
  return latnet_e, true_e

fig, ax = plt.subplots(num_iters+1,1, figsize=(plot_size_x, plot_size_y), sharex=True, sharey=True)
for i in range(num_iters+1):
  axi = ax[i]
  if i < num_iters:
    axi.set_title("Iter " + str(total_iterations[i]))
    latnet_e, true_e = read_latnet_spectrum(iteration=total_iterations[i:i+1])
  else:
    axi.set_title("Averaged")
    latnet_e, true_e = read_latnet_spectrum(iteration=total_iterations)
  x = np.arange(1,latnet_e.shape[0])
  axi.loglog(x, latnet_e[1:],      label=outputPrefix, color=latnet_color)
  axi.loglog(x, true_e[1:],     label='DNS', color=true_color)
  axi.loglog(x, np.power(x, -(5.0/3.0)), label='k^(-5/3)')
  #plt.loglog(lstm_x, lstm_e[1:],      label='LSTM 2D flow')
  #plt.loglog(lstm_x, lstm_true_e[1:],      label='true 2D flow')

  axi.set_ylabel("E(k)")

  if i == num_iters:
    axi.set_xlabel("Wavenumber K")

  axi.set_xlim(1, 10e1)
  axi.set_ylim(10e-6, 10e1)

plt.legend(loc=0)
plt.savefig(outputPath + '/' + outputPrefix + '_energy_spectrum.png')
#plt.show()
plt.close()

###############################################
## INTERMITTENCY PLOTS
###############################################

def read_latnet_intermittency(iteration=[0]):
  files = [] 

  #for it in iteration:
  #  files = files + glob.glob('diagnosticsCLSTM/iter_' + str(it).zfill(4) + '*intermittency.npz')

  for file in glob.glob(dataPath + "/iter_" + "*intermittency.npz"):
    files.append(file)
  print("intermittency: ", files)

  latnet_ebins     = []
  latnet_zdns_plot = []
  true_ebins       = []
  true_zdns_plot   = []
  for f in files:
    latnet_intermittency = np.load(f) 
    latnet_ebins.append(latnet_intermittency['net_ebins'])
    latnet_zdns_plot.append(latnet_intermittency['net_zdns_plot'])
    true_ebins.append(latnet_intermittency['true_ebins'])
    true_zdns_plot.append(latnet_intermittency['true_zdns_plot'])

  latnet_ebins = np.mean(latnet_ebins, axis=0) 
  latnet_zdns_plot = np.mean(latnet_zdns_plot, axis=0) 
  true_ebins = np.mean(true_ebins, axis=0) 
  true_zdns_plot = np.mean(true_zdns_plot, axis=0) 
  #latnet_ebins = latnet_ebins[1] 
  #latnet_zdns_plot = latnet_zdns_plot[1]
  #true_ebins = true_ebins[1]
  #true_zdns_plot = true_zdns_plot[1]
  return latnet_ebins, latnet_zdns_plot, true_ebins, true_zdns_plot


fig, ax = plt.subplots(num_iters+1, 1, figsize=(plot_size_x, plot_size_y), sharex=True, sharey=True)
for i in range(num_iters+1):
  axi = ax[i]
  if i < num_iters:
    axi.set_title("Iter " + str(total_iterations[i]))
    latnet_ebins, latnet_zdns_plot, true_ebins, true_zdns_plot = read_latnet_intermittency(iteration=total_iterations[i:i+1])
  else:
    axi.set_title("Averaged")
    latnet_ebins, latnet_zdns_plot, true_ebins, true_zdns_plot = read_latnet_intermittency(iteration=total_iterations)

  axi.semilogy(latnet_ebins, latnet_zdns_plot, label=outputPrefix, color=latnet_color)
  axi.semilogy(true_ebins, true_zdns_plot, label='DNS', color=true_color)

  axi.set_xlim(-7, 7)
  axi.set_ylim(10e-6, 1.0)

  axi.set_ylabel("p(Z)")
  if i == num_iters:
    axi.set_xlabel("Z")

plt.legend(loc=0)
plt.savefig(outputPath + '/' + outputPrefix + '_intermittency.png')
#plt.show()
plt.close()

###############################################
## STRUCTURE PLOTS
###############################################

# def read_latnet_structure(iteration=[0]):
#   files = [] 
#   for it in iteration:
#     files = files + glob.glob('diagnosticsCLSTM/iter_' + str(it).zfill(4) + '*structure.npz')
#   latnet_zeta     = []
#   true_zeta   = []
#   for f in files:
#     latnet_structure = np.load(f) 
#     latnet_zeta.append(latnet_structure['net_zeta'])
#     true_zeta.append(latnet_structure['true_zeta'])

#   latnet_zeta = np.mean(latnet_zeta, axis=0)
#   true_zeta = np.mean(true_zeta, axis=0)
#   return latnet_zeta, true_zeta

# fig, ax = plt.subplots(num_iters+1, 1, figsize=(plot_size_x, plot_size_y), sharex=True, sharey=True)
# for i in range(num_iters+1):
#   axi = ax[i]
#   if i < num_iters:
#     axi.set_title("Iter " + str(total_iterations[i]))
#     latnet_zeta, true_zeta = read_latnet_structure(iteration=total_iterations[i:i+1])
#   else:
#     axi.set_title("Averaged")
#     latnet_zeta, true_zeta = read_latnet_structure(iteration=total_iterations)

#   orders = np.arange(1, 11)
#   struct_sim = np.zeros(len(orders))
#   mu = 0.25 
#   for n in range(len(orders)):
#     struct_sim[n]= 1./3.*orders[n]*(1.0-(1.0/6.0)*mu*(orders[n]-3.0));
#   axi.scatter(orders[1:], latnet_zeta[1:], label="Neural Network", color=latnet_color)
#   axi.scatter(orders[1:], true_zeta[1:], label="DNS", color=true_color)
#   axi.plot(orders, struct_sim, label='K62')
#   axi.plot(orders, orders/3, label='K41')

#   axi.set_xlim(1, 10)
#   axi.set_ylim(0.0, 3.5)

#   axi.set_ylabel("Cn")
#   if i == num_iters:
#     axi.set_xlabel("n")

# plt.legend(loc=0)
# plt.savefig('structure.pdf')
# plt.show()
# plt.close()


###############################################
## QR PLOTS
###############################################
coarse_graining = [0,8,32]

def read_latnet_qr(iteration=[0], coarse_graining=[0, 8, 32]):
  files = [] 
  #for it in iteration:
  #  files = files + glob.glob('diagnosticsCLSTM/iter_' + str(it).zfill(4) + '*qr_data.npz')

  for file in glob.glob(dataPath + "/iter_" + "*qr_data.npz"):
    files.append(file)
  print("qr_data: ", files)

  latnet_q = {}
  latnet_r = {} 
  true_q   = {}
  true_r   = {}
  for f in files:
    latnet_qr = np.load(f)
    for cg in coarse_graining:
      if cg not in latnet_q.keys():
        latnet_q[cg] = []
      if cg not in latnet_r.keys():
        latnet_r[cg] = []
      if cg not in true_q.keys():
        true_q[cg] = []
      if cg not in true_r.keys():
        true_r[cg] = []
      latnet_q[cg].append(latnet_qr['net_q_' + str(cg)])
      latnet_r[cg].append(latnet_qr['net_r_' + str(cg)])
      true_q[cg].append(latnet_qr['true_q_' + str(cg)])
      true_r[cg].append(latnet_qr['true_r_' + str(cg)])

  for cg in coarse_graining:
    latnet_q[cg] = np.concatenate(latnet_q[cg], axis=0) 
    latnet_r[cg] = np.concatenate(latnet_r[cg], axis=0) 
    true_q[cg] = np.concatenate(true_q[cg], axis=0) 
    true_r[cg] = np.concatenate(true_r[cg], axis=0) 
  return latnet_q, latnet_r, true_q, true_r

fig, ax = plt.subplots(num_iters+1,3, figsize=(2*plot_size_x, plot_size_y), sharex=True, sharey=True)
for i in range(num_iters+1):
  if i < num_iters:
    latnet_q, latnet_r, true_q, true_r = read_latnet_qr(iteration=total_iterations[i:i+1], coarse_graining=coarse_graining)
  else:
    latnet_q, latnet_r, true_q, true_r = read_latnet_qr(iteration=total_iterations, coarse_graining=coarse_graining)

  for j in range(len(coarse_graining)):
    axi = ax[i, j]
    if i < num_iters:
      axi.set_title("Iter " + str(total_iterations[i]) + ' r = ' + str(coarse_graining[j]))
    else:
      axi.set_title("Ave r = " + str(coarse_graining[j]))

    h1 = plot_contour(axi, latnet_r[coarse_graining[j]].flatten(), latnet_q[coarse_graining[j]].flatten(), label=outputPrefix, color=latnet_color)
    h2 = plot_contour(axi, true_r[coarse_graining[j]].flatten(), true_q[coarse_graining[j]].flatten(), label="DNS", color=true_color)

    axi.set_xlim(-10, 10)
    axi.set_ylim(-10, 10)


    if i == 1:
      axi.set_ylabel("Q/(Qw)")

    if j == 0:
      axi.set_ylabel("R/(Qw)^(3/2)")

    if i == num_iters and j == len(coarse_graining)-1:
      axi.legend([h1[0], h2[0]], [outputPrefix, 'DNS'])

 
plt.savefig(outputPath + '/' + outputPrefix + '_qr_plot.png')
#plt.show()
plt.close()


###############################################
## Image PLOTS
###############################################

def read_latnet_image(iteration=0):
  files = [] 
  #files = files + glob.glob('diagnosticsCLSTM/iter_' + str(iteration).zfill(4) + '*image.npz')
  for file in glob.glob(dataPath + "/iter_" + "*image.npz"):
    files.append(file)
  print("image: ", files)

  latnet_image_stream = np.load(files[0])
  latnet_image = latnet_image_stream['net_norm']
  true_image = latnet_image_stream['true_norm']

  return latnet_image, true_image

fig, ax = plt.subplots(num_iters ,2, figsize=(2*plot_size_x, plot_size_y), sharex=True, sharey=True)
for i in range(num_iters):
  latnet_image,  true_image  = read_latnet_image(total_iterations[i])
  for j in range(2):
    axi = ax[i, j]
    if j == 0:
      axi.set_title("Lat-Net Iter " + str(i*15))
      axi.imshow(latnet_image,cmap='jet')
    elif j == 1:
      axi.set_title("True Iter " + str(i*15))
      axi.imshow(true_image,cmap='jet')

plt.savefig(outputPath + '/' + outputPrefix + '_image_plot.png')
#plt.show()
plt.close()


