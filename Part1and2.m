%% Final Project Parts 1 and 2:  Reactive Ion Etching with MIMO Control
clear;
clc;

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
Aaug=[AP,zeros(8,2);
    CP, zeros(2,2)];
Baug=[BP;zeros(2,2)];
% Baug = [BP, zeros(8,2); zeros(2, 2), -eye(2,2)];
Caug=[CP,zeros(2,2)];
D=0;

% % Q and R
Q1=[1/0.2^2, 0;
    0, 1/0.2^2];
Q=CP'*Q1*CP;
QI=[1/1.5^2, 0;
    0, 1/1.5^2];
Qaug=[Q, zeros(8,2);
    zeros(2,8), QI];
% Raug=[0.5, 0;
%     0, 1/100^2];
Raug = QI;

% LQ state feedback gain
K = lqr(Aaug,Baug,Qaug,Raug);

% Closed loop state equations with state feedback and integrators.
K1=K(1,:);
K2=K(2,:);
K11=[K1(1:8);K2(1:8)];
Ki=[K1(9:10);K2(9:10)];

Acl=[AP-BP*K11, -BP*Ki;
    CP,zeros(2,2)];
Bcl=[zeros(8,2);
    -1,0;
    0,-1];
Ccl=[CP,zeros(2,2)];
Dcl=0;
new_ss=ss(Acl,Bcl,Ccl,Dcl);

% Verify that closed-loop is stable (Check to verify no bugs in code)
eig(Acl); %should be all negatives

% Time vector
Tf = 50;
Nt = 500;
t = linspace(0,Tf,Nt);

Tr2u=ss(Acl,Bcl,-K,0);

figure(1)         %FIG 1
subplot(2,2,1)    %%top left
step(new_ss(1,1),Tf);
hold on
step(new_ss(2,1), Tf)
hold off
title('Step to [F], SF')
legend('[F]', '$V\_{bias}$', 'Interpreter', 'latex')

subplot(2,2,3)  %%bottom left
step(new_ss(1,2), Tf)
hold on
step(new_ss(2,2), Tf)
hold off
title('$\textbf{Step to $V_{\textbf{bias}}$, SF}$', 'Interpreter', 'latex')
legend('[F]', '$V\_{bias}$', 'Interpreter', 'latex')

subplot(2,2,2)    %%top right
step(Tr2u(1,1),Tf);     %treat as a matrix
hold on
step(Tr2u(2,1), Tf)
hold off
title('Step to [F], SF')
legend('RF Power', 'Throttle', 'Interpreter', 'latex')

subplot(2,2,4)  %%bottom rightI;
step(Tr2u(1,2), Tf)
hold on
step(Tr2u(2,2), Tf)
hold off
title('$\textbf{Step to $V_{\textbf{bias}}$, SF}$', 'Interpreter', 'latex')
legend('RF Power', 'Throttle', 'Interpreter', 'latex')

% plot(tunit,yunit);

%% Part 1(B): Linear Quadratic Regulator with Integrators

% Covariance matrices for loop transfer recovery observer
q = 2.5;
V = q^2*BP*BP';
W = eye(2);

% Observer gain
% Note: We only need to estimate the plant states. We do not need the
% observer to construct an estimate of the integrator states.


% L = lqe(AP, BP, CP, W, V)     %????
L = lqr(AP', CP', V, W);
L = L';

% Construct controller: 
% This includes the observer, integrators, and feedback gains.
%   Inputs: [|F| Ref.; Vbias Ref; |F| measurement; Vbias measurement]
%   Outputs: [Power; Throttle]
%%%controller
A_cont = [AP, -BP*K11, -BP*Ki;
         zeros(8,8), AP-BP*K11-L*CP, -BP*Ki;
         zeros(2,18)];
B_cont = [zeros(8,4);
          zeros(8,2), L;
          -eye(2), eye(2)];
C_cont = [zeros(2,8), -K11, -Ki];
D_cont = 0;
Tr2u_obs = ss(A_cont, B_cont, C_cont, D_cont);

% Form Closed-Loop
%   Inputs are [|F| Ref.; Vbias Ref; |F| noise; Vbias noise]
%   Outputs: [|F|; Vbias; Power; Throttle]

A_obs_cl = [AP, -BP*K11, -BP*Ki;
            L*CP, AP-BP*K11-L*CP, -BP*Ki;
            CP, zeros(2,8), zeros(2,2)];
B_obs_cl = B_cont;
C_obs_cl = [CP, zeros(2,10);
            zeros(2,8), -K11, -Ki];
D_obs_cl = 0;
new_ss_obs = ss(A_obs_cl, B_obs_cl, C_obs_cl, D_obs_cl);

% Verify that closed-loop eigenvalues are the union of the observer
% and state-feedback eigenvalues. This is a useful debugging step
% to verify that you have correctly formed the closed-loop.
% (check eig values)

% Step responses with noise on |F|
t = linspace(0,Tf,Nt);
Fnoise =0.01*randn(Nt,1);       %INCLUDE IN REPORT - 1% Noise

figure(2)         %FIG 2

subplot(2,2,1)    %%top left
y1n = lsim(-new_ss_obs(1,3),Fnoise,t);
y2n = lsim(-new_ss_obs(2,3),Fnoise,t);
y1 = step(new_ss_obs(1,1),t);
y1_noisy=y1+y1n;
plot(t,y1_noisy,'b'); %
hold on
y2 = step(new_ss_obs(2,1),t);
y2_noisy = y2+y2n;          %noise? added
plot(t,y2_noisy,'r'); %
hold off
xlim([0,50])
title('Step to [F], MIMO')
legend('[F]', '$V\_{bias}$', 'Interpreter', 'latex')

subplot(2,2,3)    %%bottom left
% Bottom left
Fnoise =0.01*randn(Nt,1);
% y1n =lsim(-new_ss_obs(1,4),Fnoise,t);
% y2n = lsim(-new_ss_obs(2,4),Fnoise,t);
y1n =lsim(-new_ss_obs(1,3),Fnoise,t);   %changed to (1,3)
y2n = lsim(-new_ss_obs(2,3),Fnoise,t);  %... (2,3)
y1=step(new_ss_obs(1,2),t);
y1_noisy=y1+y1n;
plot(t,y1_noisy,'b'); %
hold on
y2=step(new_ss_obs(2,2),t);
y2_noisy=y2+y2n;        %noise? added
plot(t,y2_noisy,'r'); %
hold off
xlim([0,50])
title('$\textbf{Step to $V_{\textbf{bias}}$, MIMO}$', 'Interpreter', 'latex')
legend('[F]', '$V\_{bias}$', 'Interpreter', 'latex')

subplot(2,2,2)    %%top right
Fnoise =0.01*randn(Nt,1);
y1n =lsim(-new_ss_obs(3,3),Fnoise,t);
y2n = lsim(-new_ss_obs(4,3),Fnoise,t);
u1=step(new_ss_obs(3,1),t);
u1_noisy=u1+y1n;
plot(t,u1_noisy,'b'); %
hold on
u2=step(new_ss_obs(4,1),t);
u2_noisy=u2+y2n;        %noise? add +y2n?
plot(t,u2_noisy,'r'); %Tr2u
hold off
xlim([0,50])
title('Step to [F], MIMO')
legend('RF Power', 'Throttle', 'Interpreter', 'latex')

subplot(2,2,4)    %%bottom right
Fnoise =0.01*randn(Nt,1);
% y1n = lsim(-new_ss_obs(3,4),Fnoise,t);
% y2n = lsim(-new_ss_obs(4,4),Fnoise,t);
y1n = lsim(-new_ss_obs(3,3),Fnoise,t);  %changed to (3,3)
y2n = lsim(-new_ss_obs(4,3),Fnoise,t);  %changed to (4,3)
u1=step(new_ss_obs(3,2),t);
u1_noisy=u1+y1n;
plot(t,u1_noisy,'b'); %
hold on
u2=step(new_ss_obs(4,2),t);
u2_noisy=u2+y2n;        %noise?
plot(t,u2_noisy,'r'); %
hold off
xlim([0,50])
title('$\textbf{Step to $V_{\textbf{bias}}$, MIMO}$', 'Interpreter', 'latex')
legend('RF Power', 'Throttle', 'Interpreter', 'latex')

% Bode magnitude from |F| noise to [|F|; Vbias; Power; Throttle]
%FIG 3 in pdf
figure(3)
bode(new_ss_obs(1,3))
hold on
bodemag(new_ss_obs(2,3))
bodemag(new_ss_obs(3,3))
bodemag(new_ss_obs(4,3))
legend("|F|", "Vbias", "RF Power", "Throttle")
title('Closed loop Bode plots of response to F noise')

%check mag at w=1
freq = 1;
[mag1, ~, ~] = bode(new_ss_obs(1,3), freq);
mag1;    %dB
[mag2, ~, ~] = bode(new_ss_obs(2,3), freq);
mag2;    %dB
[mag3, ~, ~] = bode(new_ss_obs(3,3), freq);
mag3;    %dB
[mag4, ~, ~] = bode(new_ss_obs(4,3), freq);
mag4;    %dB
%adjust q until all mags are below 0 dB


%Below: Done in Part (d)Tr2u
% Sigma magnitude from [|F| Ref.; Vbias Ref] to [Power; Throttle]
% Sigma magnitude from [|F| Noise; Vbias Noise] to [Power; Throttle]


%% Part 1(D): Stability Margins and Comparison With State Feedback

% Loop-at-a-time margins at the plant input

Gaug = ss(Aaug,Baug,Caug,0); % Plant without Observer 
Kaug = Tr2u;
% Margin_O = allmargin(Gaug*Kaug)
% Margin_I = allmargin(Kaug*Gaug)
% Augmented Matrices for G_obs 
A_obs_aug = [AP, zeros(8,8), zeros(8,2);
             AP + L*CP, -L*CP zeros(8,2) ;
             CP, zeros(2,8), zeros(2,2)];

B_obs_aug = [ BP; BP; zeros(2,2)];

C_obs_aug = [ CP zeros(2,8) zeros(2,2)];

G_obs_aug = ss(A_obs_aug,B_obs_aug,C_obs_aug,0);
K_obs_aug = Tr2u_obs(:,1:2);
% Margin_O_obs = allmargin(G_obs_aug*K_obs_aug)
% Margin_I_obs = allmargin(K_obs_aug*G_obs_aug)

%Multi-loop margins
% MMIO = diskmargin(Gaug,Kaug);
% dm_c = MMIO.DiskMargin;

% MMIO_obs = diskmargin(G_obs_aug,K_obs_aug);
% dm_c_obs = MMIO_obs.DiskMargin;

% Unstructured (fully-coupled) stability margin (USM) at the plant input
T_obs = feedback(K_obs_aug*G_obs_aug,eye(2));
m_obs=1/norm(T_obs,inf);

T = feedback(Kaug*Gaug,eye(2));
m=1/norm(T,inf);

% Input loop transfer function: Compare Lsf to LI
s = tf('s');
Lsf = K*inv(s*eye(10) - Aaug)*Baug;
K_new = [K11 zeros(2,8) Ki];
% Lobs = K_new*inv(s*eye(18) - A_obs_aug)*B_obs_aug;%Tr2u
Lobs = G_obs_aug*K_obs_aug;

% Input complementary sensitivity: Compare Tsf to TI
Tsf = feedback(Lsf,eye(2));
Tobs = feedback(Lobs,eye(2));

% Input sensitivity: Compare Ssf to SI
Ssf = feedback(eye(2),Lsf);
Sobs = feedback(eye(2),Lobs);

figure(5)
sigma(Lsf, Lobs)
xlim([0.01,100])
legend('Lsf', 'Lobs')

figure(6)
sigma(Tsf, Tobs)
xlim([0.01,100])
legend('Tsf', 'Tobs')

figure(7)
sigma(Ssf, Sobs)
xlim([0.01,100])
legend('Ssf', 'Sobs')


%% Part 2(B): Equivalent Controller

% Equivalent controller
% Q1=[1/0.1^2, 0;
%     0, 1/0.1^2];
% 
% Q=CP'*Q1*CP;
% 
% QI=[1/1.2^2, 0;
%     0, 1/1.2^2];
% 
% Qaug=[Q, zeros(8,2);
%     zeros(2,8), QI];
% 
% Raug = [0.01, 0; 0, 0.01];
% 
% K = lqr(Aaug,Baug,Qaug,Raug);
inside = tf(ss(AP, BP, K11, 0));
Ceq = inv(eye(2)+ (K11 *inv(s*eye(8)-AP)*BP))* Ki/s ;
% Ceq = (eye(2) + inside)\(Ki.*tf(1,[1,0]));
figure(8)
%bode(Ceq(1,1), Ceq(1,2), Ceq(2,1), Ceq(2,2))
%legend('$C\{eq}(1,1)$', '$C\{eq}(1,2)$', '$C\{eq}(2,1)$', '$C\{eq}(2,2)$', 'Interpreter', 'latex')

T1 = feedback(Ceq*PN,eye(2));
step(T1);
%% Part 2(C): Decentralized Approximation of Equivalent Controller

wl11=2; % to be tuned
wl12=9;
s=tf('s');
beta_l11=1.5; 
beta_l12 = 5.5;
l11=(beta_l11*s+wl11)/(s+beta_l11*wl11);
l12 =(beta_l12*s+wl12)/(s+beta_l12*wl12);
Cd1 = (l11+l12)*0.1/s;


wl2 = 0.5;
beta_l2=2.5;
lead1=(beta_l2*s+wl2)/(s+beta_l2*wl2);
Cd2=0.4/s * lead1;

% figure(9)
% bode(Kl1)
% figure(10)
% bode(Kl2)

%Stitching Ceq and Cd's
figure(9)
bode(Ceq(1,1),Cd1,Ceq(2,2),Cd2)
legend("Ceq(1,1)","Cd1","Ceq(2,2)","Cd2")

Ceq_h = [Cd1 Cd1;
        -Cd2 Cd2;];
figure(11)
% P=CP*inv(s*eye(8)-AP)*BP;
step(feedback(Ceq*PN,eye(2)))
hold on;
step(feedback(Ceq_h*PN,eye(2)))



% Teq=inv(eye(2)+P*Ceq)*P*Ceq
% figure(10)
% step(Teq,Tf)


%% Part 2: Comment on Plant Transformation
% M = [1 1; -1 1]/sqrt(2);
% MP = M*PN;
% 
% figure(12)
% subplot(2,1,1)
% bodemag(PN(1,1),'b',PN(1,2),'r--',PN(2,1),'m-.',PN(2,2),'g-.',{1e-2,1e2});
% legend('PN(1,1)','PN(1,2)','PN(2,1)','PN(2,2)','Location','Southwest');
% grid on;
% if exist('garyfyFigure','file'), garyfyFigure, end
% 
% subplot(2,1,2)
% bodemag(MP(1,1),'b',MP(1,2),'r--',MP(2,1),'m-.',MP(2,2),'g-.',{1e-2,1e2});
% legend('MP(1,1)','MP(1,2)','MP(2,1)','MP(2,2)','Location','Southwest');
% grid on;
% if exist('garyfyFigure','file'), garyfyFigure, end

