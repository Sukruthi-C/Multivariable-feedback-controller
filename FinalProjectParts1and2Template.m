%% Final Project Parts 1 and 2:  Reactive Ion Etching with MIMO Control

% Plant Model from [Power; Throttle] to [|F|; Vbias]
P1 = [tf([0.17 0.7],[1 15 26.7]); tf(0.28,[1 0.97])];
P2 = [tf(-0.17,[1 0.24]); tf([2.41 9.75],[1 4 0.7])];
P2.InputDelay = 0.5;
P = [P1 P2];

% Normalized System
DO = diag([30 350]);
DI = diag([1000 12.5]);
PN = inv(DO)*P*DI;
PN = ss(PN);

% Use second-order Pade approximation for input delay
PN = pade(PN,2);
PN.InputName = 'u';
PN.OutputName = 'y';

% State-space matrices and dimensions
[AP,BP,CP,DP] = ssdata(PN);
[nx,nu] = size(BP);
ny = size(CP,1);

%% Part 1(A): Linear Quadratic Regulator with Integrators

% Augment state equations so that you can do integral control

% LQR Weighting Matrices

% LQ state feedback gain

% Closed loop state equations with state feedback and integrators.

% Verify that closed-loop is stable (Check to verify no bugs in code)

% Time vector
Tf = 50;
Nt = 500;
t = linspace(0,Tf,Nt);

% Step responses

%% Part 1(B): Linear Quadratic Regulator with Integrators

% Covariance matrices for loop transfer recovery observer

% Observer gain
% Note: We only need to estimate the plant states. We do not need the
% observer to construct an estimate of the integrator states.

% Verify that observer error is stable (Check to verify no bugs in code)

% Construct controller: 
% This includes the observer, integrators, and feedback gains.
%   Inputs: [|F| Ref.; Vbias Ref; |F| measurement; Vbias measurement]
%   Outputs: [Power; Throttle]

% Form Closed-Loop
%   Inputs are [|F| Ref.; Vbias Ref; |F| noise; Vbias noise]
%   Outputs: [|F|; Vbias; Power; Throttle]

% Verify that closed-loop eigenvalues are the union of the observer
% and state-feedback eigenvalues. This is a useful debugging step
% to verify that you have correctly formed the closed-loop.

% Step responses with noise on |F|



% Bode magnitude from |F| noise to [|F|; Vbias; Power; Throttle]

% Sigma magnitude from [|F| Ref.; Vbias Ref] to [Power; Throttle]

% Sigma magnitude from [|F| Noise; Vbias Noise] to [Power; Throttle]


%% Part 1(D): Stability Margins and Comparison With State Feedback

% Loop-at-a-time margins at the plant input


% Unstructured (fully-coupled) stability margin (USM) at the plant input

% Input loop transfer function: Compare Lsf to LI

% Input sensitivity: Compare Ssf to SI

% Input complementary sensitivity: Compare Tsf to TI

%% Part 2(B): Equivalent Controller

% Equivalent controller
%   Ceq = inv[ I+K1 inv(sI-A) B] (KI/s)



%% Part 2(C): Decentralized Approximation of Equivalent Controller




%% Part 2: Comment on Plant Transformation
M = [1 1; -1 1]/sqrt(2);
MP = M*PN;

figure(12)
subplot(2,1,1)
bodemag(PN(1,1),'b',PN(1,2),'r--',PN(2,1),'m-.',PN(2,2),'g-.',{1e-2,1e2});
legend('PN(1,1)','PN(1,2)','PN(2,1)','PN(2,2)','Location','Southwest');
grid on;
if exist('garyfyFigure','file'), garyfyFigure, end

subplot(2,1,2)
bodemag(MP(1,1),'b',MP(1,2),'r--',MP(2,1),'m-.',MP(2,2),'g-.',{1e-2,1e2});
legend('MP(1,1)','MP(1,2)','MP(2,1)','MP(2,2)','Location','Southwest');
grid on;
if exist('garyfyFigure','file'), garyfyFigure, end

