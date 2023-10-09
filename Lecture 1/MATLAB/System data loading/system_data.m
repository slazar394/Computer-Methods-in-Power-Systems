% Script to extract the system data from the .xlsx file and convert it to a
% .mat file

clear variables
close all
clc

% Load the bus and branch data

Buses = table2array(readtable("9 bus system.xlsx","Sheet","Bus data")) ; 
Branches = table2array(readtable("9 bus system.xlsx","Sheet","Branch data")) ;

% Store the bus and branch data in a System structure and save the data

System.Buses = Buses ;
System.Branches = Branches ; 
save("9 bus system.mat","System") ; 