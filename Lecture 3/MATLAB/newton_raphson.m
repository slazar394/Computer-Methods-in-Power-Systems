% This function performs power flow analysis using Newton-Raphson's method
% and it neglects the generator reactive power limits

function [V,Theta,P,Q,Iteration] = newton_raphson(System)

% Extract the bus data and determine the IDs of PV and PQ nodes

Number_of_Buses = size(System.Buses, 1) ;
Bus_Type = System.Buses.Bus_Type ;
PQ_Buses = find(Bus_Type==3) ;

% Extract the specified voltage data

V_specified = System.Buses.V ;
Theta = System.Buses.Theta ;

% Extract the load and generation data

P_load = System.Buses.P_load ;
Q_load = System.Buses.Q_load ;
P_gen = System.Buses.P_gen ;
Q_gen = System.Buses.Q_gen ;

% Calculate the specified active and reactive power injections

P_specified = P_gen - P_load ;
Q_specified = Q_gen - Q_load ;

% Create the bus-admittance matrix and evaluate its magnitude and phase
% angle

Y_b = create_Yb(System) ;
G = real(Y_b) ; 
B = imag(Y_b) ; 

% Initialize the calculation variables

V = V_specified ;
V_old = zeros(size(V)) ;
Iteration = 0 ;
Tolerance = 1e-6 ;

% Main loop

while any(abs(V.*exp(sqrt(-1)*Theta)-V_old)>Tolerance)

    % Update the iteration variables

    V_old = V.*exp(sqrt(-1)*Theta) ;
    Iteration = Iteration + 1 ;

    % Determine active and reactive power injections

    P = zeros(Number_of_Buses,1) ;
    Q = zeros(Number_of_Buses,1) ;

    for i = 1 : Number_of_Buses

        for j = 1 : Number_of_Buses

            P(i) = P(i) + V(i)*V(j)*(G(i,j)*(cos(Theta(i)-Theta(j)))+B(i,j)*sin(Theta(i)-Theta(j))) ;
            Q(i) = Q(i) + V(i)*V(j)*(G(i,j)*(sin(Theta(i)-Theta(j)))-B(i,j)*cos(Theta(i)-Theta(j))) ;

        end

    end

    % Determine the active and reactive power injection deviations

    dP = P_specified - P ;
    dP = dP(2:Number_of_Buses) ;
    dQ = Q_specified - Q ;
    dQ = dQ(PQ_Buses) ;

    % Create the Jacobian submatrix J1 = dP/dTheta

    J1 = zeros(Number_of_Buses) ;

    for i = 1 : Number_of_Buses

        for j = 1 : Number_of_Buses

            if i == j

                for k = 1 : Number_of_Buses

                    J1(i,i) = J1(i,i) + V(i)*V(k)*(-G(i,k)*sin(Theta(i)-Theta(k))+B(i,k)*cos(Theta(i)-Theta(k))) ; 
         
                end

                J1(i,i) = J1(i,i) - V(i)^2*B(i,i) ;

            else

                J1(i,j) = V(i)*V(j)*(G(i,j)*sin(Theta(i)-Theta(j))-B(i,j)*cos(Theta(i)-Theta(j))) ; 

            end

        end

    end

    J1 = J1(2:end,2:end) ;

    % Create the Jacobian submatrix J2 = dP/dV

    J2 = zeros(Number_of_Buses) ;

    for i = 1 : Number_of_Buses

        for j = 1 : Number_of_Buses

            if i == j

                for k = 1 : Number_of_Buses

                    J2(i,i) = J2(i,i) + V(k)*(G(i,k)*cos(Theta(i)-Theta(k))+B(i,k)*sin(Theta(i)-Theta(k))) ;

                end

                J2(i,i) = J2(i,i) + G(i,i)*V(i) ; 

            else

                J2(i,j) = V(i)*(G(i,j)*cos(Theta(i)-Theta(j))+B(i,j)*sin(Theta(i)-Theta(j))) ;

            end

        end

    end

    J2 = J2(2:end,PQ_Buses) ;

    % Create the Jacobian submatrix J3 = dQ/dTheta

    J3 = zeros(Number_of_Buses) ;

    for i = 1 : Number_of_Buses

        for j = 1 : Number_of_Buses

            if i == j

                for k = 1 : Number_of_Buses

                    J3(i,i) = J3(i,i) + V(i)*V(k)*(G(i,k)*cos(Theta(i)-Theta(k))+B(i,k)*sin(Theta(i)-Theta(k))) ;

                end

                J3(i,i) = J3(i,i) - G(i,i)*V(i)^2 ; 

            else

                J3(i,j) = -V(i)*V(j)*(G(i,j)*cos(Theta(i)-Theta(j))+B(i,j)*sin(Theta(i)-Theta(j))) ;

            end

        end

    end

    J3 = J3(PQ_Buses,2:end) ;

    % Create the Jacobian submatrix J4 = dQ/dTheta

    J4 = zeros(Number_of_Buses) ;

    for i = 1 : Number_of_Buses

        for j = 1 : Number_of_Buses

            if i == j

                for k = 1 : Number_of_Buses

                    J4(i,i) = J4(i,i) + V(k)*(G(i,k)*sin(Theta(i)-Theta(k))-B(i,k)*cos(Theta(i)-Theta(k))) ;

                end

                J4(i,i) = J4(i,i) - B(i,i)*V(i) ; 

            else

                J4(i,j) = V(i)*(G(i,j)*sin(Theta(i)-Theta(j))-B(i,j)*cos(Theta(i)-Theta(j))) ;

            end

        end

    end

    J4 = J4(PQ_Buses,PQ_Buses) ;

    % Create the complete Jacobian matrix

    J = [J1 J2 ; J3 J4] ;

    % Solve the linear system of equations

    dX = J\[dP ; dQ] ;

    % Update the voltage magnitudes and voltage phase angles

    dTh = dX(1:Number_of_Buses-1) ;
    dV = dX(Number_of_Buses:end) ;
    Theta(2:Number_of_Buses) = Theta(2:Number_of_Buses) + dTh ;
    V(PQ_Buses) = V(PQ_Buses) + dV ;

end

Theta = Theta*180/pi ;

end