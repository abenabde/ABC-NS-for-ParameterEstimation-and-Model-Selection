%% Model Selection and Parameter Estimation using ABC-NS algorithm: Linear model vs. cubic model
clc;clear;close all;
load('training_data.mat');
load('excitation.mat');
% true_value = [0.05 50 1000];
tol = 100; % initial tolerance threshold 
accuracy=0.01; % change the value of accuracy as desired
theta_ABC = abc_ms(um_1,sig,tol,accuracy) 
