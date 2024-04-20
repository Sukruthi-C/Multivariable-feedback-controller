%% Final Project Part 3:  Reactive Ion Etching with MIMO Control

%% Part 3(A) -- DC Analysis With Oxygen Sensor

% Model
% Inputs: [Power; Throttle; %O2]
% Outputs: [|F|; Vbias; Pressure]
P1 = [zpk(-0.067,[-0.095 -19.69],0.49); ...
    zpk(-0.27,[-0.19 -62.42],12.23); ...
    zpk(0.006,[-0.19 -2.33],-0.011)];

P2 = [zpk(0.73,[-0.11; -39.76],4.85); tf(1.65,[1 0.16]); ...
      zpk([],[-0.18; -3],-0.97)];
P2.InputDelay = 0.42;

P3 = [tf(0.33,[1 0.17]); tf(0.25,[1 0.41]); tf(0.024,[1 0.4])];
P3.InputDelay = 0.77;

Pox = [P1 P2 P3];

% Use second-order Pade for plant
Pox = pade(Pox, 2);

% Input and output scalings based on equilibrium values

% Normalize Plant

% Condition number of DC gain

%% Part 3(B) -- DC Analysis With Oxygen Sensor

% State-space data for scaled plant
[As,Bs,Cs] = ssdata(PoxN);
[nx,nu] = size(Bs);
ny = size(Cs,1);


% Weighting matrices (Q,R,V,W)
% Assume Q of the form Q = blkdiag(alpha*Cs'*Cs, Qw)

% Augmented Plant with integrators

% Compute state feedback and observer gains

% Construct controller: 
% This includes the observer, integrators, and feedback gains.
%   Inputs: [|F|Ref; Vbias Ref; Press Ref; |F| Meas; Vbias Meas; Press Meas]
%   Outputs: [Power; Throttle; %O2]



% Form Closed-Loop
%   Inputs: [|F|Ref; Vbias Ref; Press Ref; |F| Noise; Vbias Noise; Press Noise]
%   Outputs: [|F|; Vbias; Press; Power; Throttle; %O2]



% Verify that closed-loop eigenvalues are the union of the observer
% and state-feedback eigenvalues. This is a useful debugging step
% to verify that you have correctly formed the closed-loop.




% Time vector
Tf = 40;
Nt = 400;
t = linspace(0,Tf,Nt);

% Step Responses (Without Noise)
