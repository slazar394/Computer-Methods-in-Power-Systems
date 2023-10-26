% Script to compare the power flow results obtained using Shirmohammadi's
% and DistFlow methods

clear variables
close all
clc

% Load the system data

load('13 bus system.mat','System') ; 

% Visualize the test system

G = graph(System.Branches.From_Bus,System.Branches.To_Bus) ; 
figure
plot(G,'Layout','Layered')

%  Perform the power flow analysis using Shirmohammadi's method

[V_sh,Iteration_sh] = shirmohammadi(System) ; 

% Perform the power flow analysis using DistFlow

[V_df,Iteration_df] = distflow(System) ; 

% Compare the voltage magnitudes

dV = abs(abs(V_sh)-V_df) ; 
disp(strcat("The mean voltage deviation is: ", num2str(mean(dV))))

% Compare the number of iterations

disp("The number of iterations to convergence:")
disp(strcat("Shirmohammadi - ",num2str(Iteration_sh)))
disp(strcat("DistFlow - ",num2str(Iteration_df)))