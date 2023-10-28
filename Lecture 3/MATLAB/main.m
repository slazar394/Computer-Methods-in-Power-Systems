% Script to compare the power flow results obtained using Gauss-Seidel's
% and Newton-Raphson's methods

clear variables
close all
clc

% Load the system data

load('9 bus system.mat','System') ; 

% Run the power flow using Gauss-Seidel's method

[V_gs,Theta_gs,Iteration_gs] = gauss_seidel(System) ;

% Run the power flow using Gauss-Seidel's method

[V_nr,Theta_nr,~,~,Iteration_nr] = newton_raphson(System) ;

% Compare the voltage magnitudes

dV = abs(abs(V_nr)-V_gs) ; 
disp(strcat("The mean voltage deviation is: ", num2str(mean(dV))))

% Compare the number of iterations

disp("The number of iterations to convergence:")
disp(strcat("Gauss-Seidel: ",num2str(Iteration_gs)))
disp(strcat("Newton-Rapshon: ",num2str(Iteration_nr)))