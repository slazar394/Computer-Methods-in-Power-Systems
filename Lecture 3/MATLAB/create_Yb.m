% Function to create the bus-admittance matrix using the non-singular
% transformation method

function Yb = create_Yb(System)

% Determine the number of independent buses and the number of branches in
% the system

Number_of_Buses = size(System.Buses,1) ;
Number_of_Branches = size(System.Branches,1) ;

% Extract the branch data

From_Bus = System.Branches.From_Bus ;
To_Bus = System.Branches.To_Bus ;
R = System.Branches.R ;
X = System.Branches.X ;
B = System.Branches.B ;

% Initialize the bus-admittance matrix

Yb = zeros(Number_of_Buses) ;

% Form the bus-admittance matrix

for i = 1 : Number_of_Branches
    
    % Off-diagonal elements
    
    Yb(From_Bus(i),To_Bus(i)) = Yb(From_Bus(i),To_Bus(i)) - ...
        1/(R(i)+sqrt(-1)*X(i)) ;
    Yb(To_Bus(i),From_Bus(i)) = Yb(To_Bus(i),From_Bus(i)) - ...
        1/(R(i)+sqrt(-1)*X(i)) ;
    
    % Diagonal elements
    
    Yb(From_Bus(i),From_Bus(i)) = Yb(From_Bus(i),From_Bus(i)) + ...
        1/(R(i)+sqrt(-1)*X(i)) + sqrt(-1)*B(i)/2 ;
    Yb(To_Bus(i),To_Bus(i)) = Yb(To_Bus(i),To_Bus(i)) + ...
        1/(R(i)+sqrt(-1)*X(i)) + sqrt(-1)*B(i)/2 ;
    
end

end