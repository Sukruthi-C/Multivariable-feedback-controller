%% Final Project Part 3:  Reactive Ion Etching with MIMO Control
clc;
clear;
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

DI = diag([1000, 12.5, 5]);
DO = diag([16.52, 340, 17.83]);

% Normalize Plant
PN = inv(DO)*Pox*DI;

%Get DC Gain:
PN_0 = evalfr(PN, 0);

% Condition number of DC gain
cond_no = cond(PN_0)    %far from inf, so far from singular = feasible.

%% Part 3(B) -- DC Analysis With Oxygen Sensor

% State-space data for scaled plant
PoxN = PN;
[As,Bs,Cs] = ssdata(PoxN);
[nx,nu] = size(Bs);
ny = size(Cs,1);

% Weighting matrices (Q,R,V,W)
% Assume Q of the form Q = blkdiag(alpha*Cs'*Cs, Qw)

alpha = 1;
Qw=diag([1, 1, 1]);
Q = blkdiag(alpha*Cs'*Cs, Qw);
Raug = Qw;

% Augmented Plant with integrators
Aaug=[As, zeros(26,3);
      Cs, zeros(3,3)];
Baug=[Bs; zeros(3,3)];
Caug=[Cs, zeros(3,3)];
D=0;

% Compute state feedback and observer gains
K = lqr(Aaug, Baug, Q, Raug);
Ki = K(:, 27:29);
K1 = K(:, 1:26);

V = Bs*Bs';
W = eye(3);

L = lqr(As', Cs', V, W);
L = L';

% Construct controller: 
% This includes the observer, integrators, and feedback gains.
%   Inputs: [|F|Ref; Vbias Ref; Press Ref; |F| Meas; Vbias Meas; Press Meas]
%   Outputs: [Power; Throttle; %O2]
A_cont = [As, -Bs*K1, -Bs*Ki;
         zeros(26,26), As-Bs*K1-L*Cs, -Bs*Ki;
         zeros(3,55)];
B_cont = [zeros(26,6);
          zeros(26,3), L;
          -eye(3),     eye(3)];
C_cont = [zeros(3,26), -K1, -Ki];
D_cont = 0;
Tr2u_obs = ss(A_cont, B_cont, C_cont, D_cont);


% Form Closed-Loop
%   Inputs: [|F|Ref; Vbias Ref; Press Ref; |F| Noise; Vbias Noise; Press Noise]
%   Outputs: [|F|; Vbias; Press; Power; Throttle; %O2]
A_obs_cl = [As, -Bs*K1, -Bs*Ki;
            L*Cs, As-Bs*K1-L*Cs, -Bs*Ki;
            Cs, zeros(3,26), zeros(3,3)];
B_obs_cl = B_cont;
C_obs_cl = [Cs, zeros(3,29);
            zeros(3,26), -K1, -Ki];
D_obs_cl = 0;
new_ss_obs = ss(A_obs_cl, B_obs_cl, C_obs_cl, D_obs_cl);


% Verify that closed-loop eigenvalues are the union of the observer
% and state-feedback eigenvalues. This is a useful debugging step
% to verify that you have correctly formed the closed-loop.
% eig(new_ss_obs);
% eig()


% Time vector
Tf = 40;
Nt = 400;
t = linspace(0,Tf,Nt);

% Step Responses (Without Noise)
figure(1)         %FIG 1
subplot(2,3,1)    %%top left
step(new_ss_obs(1,1),Tf);
hold on
step(new_ss_obs(2,1), Tf)
step(new_ss_obs(3,1), Tf)
hold off
title('Step to [F]')
ylim([-0.2, 1.2])

subplot(2,3,4)  %%bottom left
step(new_ss_obs(4,1), Tf)
hold on
step(new_ss_obs(5,1), Tf)
step(new_ss_obs(6,1), Tf)
hold off
title("")
ylim([-0.5, 2.5])

subplot(2,3,2)    %%top middle
step(new_ss_obs(1,2),Tf);     %treat as a matrix
hold on
step(new_ss_obs(2,2), Tf)
step(new_ss_obs(3,2), Tf)
hold off
title('$\textbf{Step in $V_{\textbf{bias}}$}$', 'Interpreter', 'latex')
ylim([-0.5, 1.5])

subplot(2,3,5)  %%bottom middle
step(new_ss_obs(4,2), Tf)
hold on
step(new_ss_obs(5,2), Tf)
step(new_ss_obs(6,2), Tf)
hold off
title('')
ylim([-4, 2])

subplot(2,3,3)    %%top right
step(new_ss_obs(1,3),Tf);     %treat as a matrix
hold on
step(new_ss_obs(2,3), Tf)
step(new_ss_obs(3,3), Tf)
hold off
title('Step in Pressure')
legend('[F]', '$V\_{bias}$', 'Pressure', 'Interpreter', 'latex')
ylim([-0.2, 1.2])

subplot(2,3,6)  %%bottom right
step(new_ss_obs(4,3), Tf)
hold on
step(new_ss_obs(5,3), Tf)
step(new_ss_obs(6,3), Tf)
hold off
title('')
legend('Power', 'Throttle', '%O2')
ylim([-2, 0.5])