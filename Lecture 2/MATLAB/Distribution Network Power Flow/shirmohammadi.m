function [V,Iteration] = shirmohammadi(System)

% Extract the bus data

Number_of_Buses = size(System.Buses,1) ;
P_load = System.Buses.P_load ;
Q_load = System.Buses.Q_load ;

% Extract the branch data

Number_of_Branches = size(System.Branches,1) ;
From_Bus = System.Branches.From_Bus ;
To_Bus = System.Branches.To_Bus ;
R = System.Branches.R ;
X = System.Branches.X ;

% Initialize the calculation variables: complex node power injections,
% branch impedances, complex node voltages, complex node current
% injections, complex branch currents, as well as iteration variables

S = - P_load - sqrt(-1)*Q_load ;
Z = R + sqrt(-1)*X ;
V = ones(Number_of_Buses,1) ;
I = zeros(Number_of_Buses,1) ;
J = zeros(Number_of_Branches,1) ;
V_old = 0*V ;
Tolerance = 1e-6 ;
Iteration = 0 ;

% Main loop

while any(abs(V-V_old)>Tolerance)

    % Update the iteration variables

    V_old = V ;
    Iteration = Iteration + 1 ;

    % Calculate the complex node current injections

    for i = 1 : Number_of_Buses

        I(i) = conj(S(i)/V(i)) ;

    end

    % Backward sweep: calculate the complex branch currents

    for i = Number_of_Branches : -1 : 1

        J(i) = - I(To_Bus(i)) ;

        for j = 1 : Number_of_Branches

            if To_Bus(i) == From_Bus(j)

                J(i) = J(i) + J(j) ;

            end

        end

    end

    % Forward sweep: calculate the complex node voltages

    for i = 1 : Number_of_Branches

        V(To_Bus(i)) = V(From_Bus(i)) - Z(i)*J(i) ;

    end

end

end