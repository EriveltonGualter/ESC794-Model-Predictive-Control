
close all

%% Question 2

[tf, xf, uf] = hw5_MPCdlqrNoTerminal;
[tc, xc, uc] = MPCdlqrNoTerminal;

%% Question 3

hw5_CTdoubleIntMPCNoTerm
CTdoubleIntMPCeqNoTerm

%% Question 4
question4

%% Question 5
hw5_acadoDTocp

%% Coments
clc
display('Question 2');
display('For the Discrete MPC approach, the difference between the state');
display('response and control input using terminal Constraints are relativetly');
display('the same. Inclusing the time consuming to perform both operations');
display(' ');

display('Question 3');
display('FOr the Continuous approach, the advantages to use the terminal constraints');
display('is really clear. It boiwls down to the fact, the time the system reachs');
display('the equilibrium point with the terminal constraints is much faster');
display('Note that for the same number of interations, the control with terminal');
display('constraints reach less then 20 interaction, while with the free it');
display('does not even reach the equilibrium points');
