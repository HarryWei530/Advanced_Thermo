% Wksp1A Harry Wei
% code computes cpair for temperature 300 - 3000K, plots cp, cv vs temp,
% and K = cp/cv vs temp

T = linspace(300,3000,1000);
[cp,cv] = cpair(T);
K = cp./cv;

figure(1);
subplot(2,1,1);
plot(T,cp,T,cv);
legend("Cp_{air}","Cv_{air}");
title("Specific heat vs Temperature")
xlabel("Temperature [K]")
ylabel("Specific heat constants [KJ/Kg/K]")
grid on;
subplot(2,1,2);
plot(T,K);
title("Ratio of specific heat vs Temperature")
xlabel("Temperature [K]")
ylabel("Specific heat ratio")
grid on;

