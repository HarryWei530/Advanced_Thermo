function [cp,cv] = cpair(T)
% [cp,cv] = cpair(temperature), returns cp and cv in strict SI unit, J/Kg/K
    M_air = 28.84; % Kg/Kmol
    R_air = 288; % J/Kg/K
    cp = 1000*(27.204 + 0.0050562*T + 2.4187*10^(-6)*T.^2 - 2.0043*10^(-9)*T.^3+3.4233*10^(-13)*T.^4)/M_air; %Cp in J/Kg/K
    cv = cp - R_air;
end
