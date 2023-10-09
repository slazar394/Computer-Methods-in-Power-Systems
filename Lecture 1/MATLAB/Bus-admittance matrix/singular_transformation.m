% Function to create the bus-admittance matrix using the singular
% transformation method

function Yb = singular_transformation(System)

% Determine the number of independent buses and the number of branches in
% the system

Number_of_Buses = size(System.Buses,1) ; 
Number_of_Branches = size(System.Branches,1) ; 

% Extract the branch data

From_Bus = System.Branches(:,1) ;
To_Bus = System.Branches(:,2) ; 
R = System.Branches(:,3) ; 
X = System.Branches(:,4) ; 
B = System.Branches(:,5) ;

% Form the branch-bus incidence matrix A and the branch admittance matrix Y
% assuming that each element is represented by a PI-section equivalent
% circuit consisting of a series impedance with two parallel admittances

A = zeros(Number_of_Buses,Number_of_Branches) ; 
Y = zeros(Number_of_Branches,1) ; 

for i = 1 : Number_of_Branches
    
    % Series elements between two independent buses
    
    A(From_Bus(i),i) = 1 ; 
    A(To_Bus(i),i) = -1 ;
    Y(i) = 1/(R(i)+sqrt(-1)*X(i)) ; 
    
    % Parallel elements between an independent bus and a referent bus: if
    % branch admittance is nonzero, two additional branches arise
    
    if B(i) ~= 0
        
       A = [A zeros(Number_of_Buses,2)] ;
       A(From_Bus(i),end-1) = 1 ;
       A(To_Bus(i),end) = 1 ;
       Y = [Y ; zeros(2,1)] ;
       Y(end-1) = sqrt(-1)*B(i)/2 ; 
       Y(end) = sqrt(-1)*B(i)/2 ; 
        
    end
    
end

% Since the branch admittance matrix is formed as a vector, it needs to be
% converted to an equivalent diagonal matrix

Y = diag(Y) ; 

% Form the bus admittance matrix

Yb = A*Y*A' ; 

end