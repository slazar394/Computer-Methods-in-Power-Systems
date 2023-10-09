% Script to compare the methods of singular and non-singular transformation
% for the creation of the bus-admittance matrix

clear variables
close all
clc

% Load the system data

load("9 bus system.mat","System") ; 

% Create the bus-addmitance matrix using both methods

Y_b_1 = singular_transformation(System) ; 
Y_b_2 = nonsingular_transformation(System) ;

% Check the difference between the two matrices

dY = sum(sum(abs(Y_b_1-Y_b_2))) ;
disp(dY)