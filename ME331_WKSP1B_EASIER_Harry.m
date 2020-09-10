% Wksp1B Harry Wei
% code computes Determine, as a function of engine speed, the energy rate 
% quantities (power, rate of heat transfer, eXergy destruction rate), the 
% total energy quantities (total work required, the total heat transferred
% , and the total eXergy destroyed).

clc;clear all;close all;
%% Input parameters
Vd = 700*10^-6; % m3 Displacement volume
r = 10; % compression ratio
B_S = 1.1; % bore to stroke ratio
L_a = 3; % connecting rod to crank radius ratio
P1 = 0.7*101325; % Pa initial pressure
T1 = 330; % K initial temperature
Twalls = 298; % K fixed inner wall temperature
T0 = 298; % K Dead state
hconv = 15; % W/m2/K convection coefficient
R = 287; % J/kg/K Gas constant for air
steps = 10000; % Number of steps for simulator
N = 2000/60; %RPS
theta0 = pi; %180 deg, or pi rad
%% Derived Parameters
Vc = Vd/(r-1);
B = ((Vd*4*B_S)/pi)^(1/3);
S = B/ B_S;
a = S/2;
L = L_a*a;
M = P1*(Vd+Vc)/(R*T1);
delta_t = 1/(2*N)/steps; % dt for one revolution and 1000 steps
%% Known Vectors
t = linspace(0,delta_t*steps,steps);
theta = theta0 - 2*pi*N*t;
gamma = asin(a/L*sin(theta));
V = Vc+Vd/2*((1-cos(theta)+L/a*(1-cos(gamma))));
dV_dt = -pi*Vd*N*(sin(theta)+tan(gamma).*cos(theta));
A = pi*B^2/2+4*(Vc/B)+pi*B*(a*(1-cos(theta))+L*(1-cos(gamma)));

%% Calorically Perfect Model
T_perfect = zeros(size(t));
T_perfect(1) = T1;
P_perfect = zeros(size(t));
P_perfect(1) = P1;

[cp,~] = cpair(T1);
cv = cp-R;
k = cp/cv;
for i = 2:length(T_perfect)
    T_perfect(i) = T_perfect(1)*(V(1)/V(i))^(k-1);
    P_perfect(i) = P_perfect(1)*(V(1)/V(i))^k;
end
fprintf("T1_perfect = %.2f K, Tf_perfect = %.2f K\n",T_perfect(1),T_perfect(end));
fprintf("P1_perfect = %.2f Pa, Pf_perfect = %.2f Pa\n",P_perfect(1),P_perfect(end));

Q_dot_in_perfect = hconv*A.*(Twalls-T_perfect);
W_dot_out_perfect = P_perfect.* dV_dt;
Xd_dot_perfect = T0*Q_dot_in_perfect.*(1./T_perfect - 1/Twalls);

%% Simulation
T = zeros(size(t));
T(1) = T1;
P = zeros(size(t));
P(1) = P1;
for i = 2:length(T)
    [cp_t,cv_t] = cpair(T(i-1));
    T(i) = T(i-1)+delta_t/(M*cv_t)*(hconv*A(i)*(Twalls - T(i-1))-P(i-1)*dV_dt(i-1));
    P(i) = M*R*T(i)/V(i);
end
fprintf("T1 = %.2f K, Tf = %.2f K\n",T_perfect(1),T_perfect(end));
fprintf("P1 = %.2f Pa, Pf = %.2f Pa\n",P_perfect(1),P_perfect(end));

Q_dot_in = hconv*A.*(Twalls-T);
W_dot_out = P.* dV_dt;
Xd_dot = T0*Q_dot_in.*(1./T - 1/Twalls);

% Not printing correct values yet
fprintf("Ideal Heat Out(J) = %.2f, Heat Out(J) = %.2f\n",-Q_dot_in_perfect(1),-Q_dot_in(1));
fprintf("Ideal Work In(J) = %.2f, Work In(J) = %.2f\n",max(-W_dot_out_perfect),max(-W_dot_out));
fprintf("Ideal Exergy Destroyed(J) = %.2f, Exergy Destroyed(J) = %.2f\n",Xd_dot_perfect(1),Xd_dot(1));

%% Results
figure();
subplot(1,2,1);
plot(theta,T_perfect,theta,T);
legend("Calorically Perfect","Thermally Perfect");
title("Temperature")
grid on;
xlabel("Theta [rad]");
ylabel("Temperature [K]")

subplot(1,2,2);
plot(theta,P_perfect,theta,P)
legend("Calorically Perfect","Thermally Perfect");
title("Pressure")
grid on;
xlabel("Theta [rad]");
ylabel("Pressure [Pa]")

figure();
subplot(2,1,1);
plot(theta,Q_dot_in_perfect,theta,W_dot_out_perfect,theta,Xd_dot_perfect);
legend("Q_{in}","Work_{out}","eXergy destruction");
title("Rate of quantities for calorically perfect model")
grid on;
xlabel("Theta [rad]");
ylabel("[J/s]");
ylim([-6*10^4,2*10^4]);

subplot(2,1,2);
plot(theta,Q_dot_in,theta,W_dot_out,theta,Xd_dot);
legend("Q_{in}","Work_{out}","eXergy destruction");
title("Thermally perfect model")
grid on;
xlabel("Theta [rad]");
ylabel("[J/s]");
ylim([-6*10^4,2*10^4]);
