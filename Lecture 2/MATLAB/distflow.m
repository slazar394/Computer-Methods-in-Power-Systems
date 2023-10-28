% Function to perform power flow calculation using the DistFlow method

function [V,Iteration] = distflow(System)

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

% Initialize the calculation variables: arrays of newest and old voltage
% magnitudes, active and reactive power flows at the sending and receiving 
% ends of the branch, as well as iteration variables

V = ones(Number_of_Buses,1) ;
V_old = zeros(Number_of_Buses,1) ;
P = zeros(Number_of_Branches,1) ;
Q = zeros(Number_of_Branches,1) ;
P_rec = zeros(Number_of_Branches,1) ;
Q_rec = zeros(Number_of_Branches,1) ;
Iteration = 0 ;
Tolerance = 1e-6 ;

% Main loop

while any(abs(V-V_old)>=Tolerance)

    % Update the iteration variables

    V_old = V ;
    Iteration = Iteration + 1 ;

    % Backward sweep: calculate active and reactive powers at the sending
    % and receiving ends of each branch

    for i = Number_of_Branches : -1 : 1

        P_rec(i) = P_load(To_Bus(i)) ;
        Q_rec(i) = Q_load(To_Bus(i)) ;

        for j = 1 : Number_of_Branches

            % Check if branch "j" starts from the ending node of branch "i"

            if To_Bus(i) == From_Bus(j)

                P_rec(i) = P_rec(i) + P(j) ;
                Q_rec(i) = Q_rec(i) + Q(j) ;

            end

        end

        % Power flows at the sending end of the branch

        P(i) = P_rec(i) + R(i)*(P_rec(i)^2 + Q_rec(i)^2)/V(To_Bus(i)) ;
        Q(i) = Q_rec(i) + X(i)*(P_rec(i)^2 + Q_rec(i)^2)/V(To_Bus(i)) ;

    end

    % Forward sweep: calculate node voltages

    for i = 1 : Number_of_Branches

        V(To_Bus(i)) = V(From_Bus(i)) - 2*(P(i)*R(i) + Q(i)*X(i)) + ...
            (R(i)^2+X(i)^2)*(P(i)^2+Q(i)^2)/V(From_Bus(i)) ;

    end

end

% Perform voltage correction (for easier programming, up to this point, V
% represents an array of squared voltage magnitudes)

V = sqrt(V) ; 

end