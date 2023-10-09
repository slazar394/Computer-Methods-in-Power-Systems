% Function to perform the power flow using the linear DC method

function [Theta,P_branch] = DC_power_flow(System)

% Determine the number of buses and the number of branches in the system

Number_of_Buses = size(System.Buses,1) ;
Number_of_Branches = size(System.Branches,1) ;

% Extract the system data

P_gen = System.Buses(:,7) ; 
P_load = System.Buses(:,5) ;
From_Bus = System.Branches(:,1) ;
To_Bus = System.Branches(:,2) ;
X = System.Branches(:,4) ;

% Create the bus-susceptance matrix

B = zeros(Number_of_Buses) ; 

for i = 1 : Number_of_Branches

    % Off-diagonal elements

    B(From_Bus(i),To_Bus(i)) = B(From_Bus(i),To_Bus(i)) + 1/X(i);
    B(To_Bus(i),From_Bus(i)) = B(To_Bus(i),From_Bus(i)) + 1/X(i) ; 

    % Diagonal elements

    B(From_Bus(i),From_Bus(i)) = B(From_Bus(i),From_Bus(i)) - 1/X(i) ;
    B(To_Bus(i),To_Bus(i)) = B(To_Bus(i),To_Bus(i)) - 1/X(i) ; 

end

% Remove the row and column associated with the slack bus to create a
% reduced bus-susceptance matrix

B_r = B ;
B_r(1,:) = [] ;
B_r(:,1) = [] ;

% Create the bus injection vector

P = P_gen - P_load ;

% Remove the row associated with the slack bus to create a reduced
% bus-injection vector

P_r = P(2:end) ; 

% Determine the voltage phase angles

Theta = zeros(Number_of_Buses,1) ; 
Theta(2:end) = -inv(B_r)*P_r ; 

% Calculate the branch active power flows

P_branch = zeros(Number_of_Branches,1) ; 

for i = 1 : Number_of_Branches

    P_branch(i) = (Theta(From_Bus(i))-Theta(To_Bus(i)))/X(i) ; 

end

end