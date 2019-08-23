% Erivelton Gualter 
%
% Data created: 10/21/2018
%
% Homework 4
% ESC794: Selected Topics in Engineering Science Model Predictive Control

clear all
close all
clc

% Samples of Matrix
A = [-1 1; -0.5 3]; B = [1.5; 0.5];

step = 200; % It defines the resolution presicion. 

% Question 1-b
getBoundaries(A, B, step);

% Question 1-c
getBoundariesHorizon(A, B, step, 3);

% Question 2

% System Information
A = [1 0 1; 0 0 -1; 1 2 1]; 
B = [2 0; -1 0; 0 1]; 

% LQR Parameters
Q = eye(3); 
R = eye(2);
Qf = 10*eye(3);

X0 = [ 1 1 1 ]; % Initial Condition
tf = 7;         % Final time [s]
dt = 1e-1;      % Sampling time

FiniteHorizon_LQR(A, B, X0, Q, R, Qf, dt, tf)