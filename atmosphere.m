function [T, rho, mu, a] = atmosphere(h)
h = max(0, min(h, 18000));
if h >= 0 && h <= 11000 
    T = 288.16 - 0.0065*h;
    rho = 1.225 * (T/288.16)^4.259;
elseif h > 11000 && h <= 18000
        T = 216.66;
    rho = 2.0813 * exp(-1.5776e-4*h);
    
end
mu = 1.716e-5 * (T/288.15)^(3/2) * (288.15 + 110.4) / (T + 110.4); %sutherland viskozitesi
gamma = 1.4;
R = 287;
a = sqrt(gamma*R*T);
end