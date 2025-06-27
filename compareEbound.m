clear
clc
close all

addpath('stability/dt_comparison/full');

snapshot_data = 'results/shear_layer_1.000e+03_70x70_4/matlab_data.mat';
snapshots = load(snapshot_data,'uh_total');

M = [16,32,64,128,256];
