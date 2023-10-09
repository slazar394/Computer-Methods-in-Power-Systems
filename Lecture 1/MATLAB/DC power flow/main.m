% Script to test the DC power flow method

clear varibles
close all
clc

% Load the system data

load('9 bus system.mat','System') ; 

% Perform the DC power flow

[Theta,P_branch] = DC_power_flow(System) ;