function [V,Theta,Iteration] = gauss_seidel(System)

% Extract the system data

Number_of_Buses = size(System.Buses, 1) ;
Bus_Type = System.Buses.Bus_Type ;
V_specified = System.Buses.V ;
Theta = System.Buses.Theta ;
P_load = System.Buses.P_load ;
Q_load = System.Buses.Q_load ;
P_gen = System.Buses.P_gen ;
Q_gen = System.Buses.Q_gen ;
Q_gen_min = System.Buses.Q_gen_min ;
Q_gen_max = System.Buses.Q_gen_max ;

% Create the bus-admittance matrix

Y_b = create_Yb(System) ;

% Initialize the calculation variables

P = P_gen - P_load ;
Q = Q_gen - Q_load ;
V = V_specified.*exp(sqrt(-1)*Theta) ;
V_old = zeros(size(V)) ;
Iteration = 0 ; 
Tolerance = 1e-6 ;

% Precompute reciprocal of diagonal elements of Y_b

Y_b_diag_inv = 1./diag(Y_b) ;

% Main loop

while any(abs(V - V_old) > Tolerance)

    % Update the iteration variables
    
    V_old = V ;
    Iteration = Iteration + 1 ; 

    % Apply the Gauss-Seidel power flow equations to each node except the slack node
    
    for i = 2 : Number_of_Buses

        % Calculate the sum of node current contributions

        S = Y_b(i,:)*V ;  
        S = S - Y_b(i,i)*V(i) ; 

        if Bus_Type(i) == 2

            % If the current node represents a PV node, determine the
            % reactive power injection
            
            Q(i) = - imag(conj(V(i))*(S+Y_b(i,i)*V(i))) ;
            Q_gen(i) = Q(i) + Q_load(i) ;

            % Check the generator reactive power limits
            
            if Q_gen(i) < Q_gen_min(i) || Q_gen(i) > Q_gen_max(i)
                
                % If reactive power limits are not satisfied, limit the
                % reactive power value and update the voltage magnitude and
                % phase angle

                Q_gen(i) = max(min(Q_gen(i),Q_gen_max(i)),Q_gen_min(i)) ;
                Q(i) = Q_gen(i) - Q_load(i) ;
                V(i) = Y_b_diag_inv(i)*((P(i)-sqrt(-1)*Q(i))/conj(V(i))-S) ;
            
            else

                % If reactive power limits are satisfied, update only the
                % voltage phase angle

                V(i) = Y_b_diag_inv(i)*((P(i)-sqrt(-1)*Q(i))/conj(V(i))-S) ;
                V(i) = abs(V_specified(i))*exp(sqrt(-1)*angle(V(i))) ;
            
            end
        
        else
            
            % If the current node represents a PQ node, update the complex
            % node voltage
            
            V(i) = Y_b_diag_inv(i)*((P(i)-sqrt(-1)*Q(i))/conj(V(i))-S) ;
        
        end
    
    end

end

% Evaluate the voltage magnitudes and phase angles

Theta = angle(V)*180/pi ;
V = abs(V) ;

end