function [CD0,Re] = cd0_hull_wing(rho, V, mu, S, cmac, tc, M)

V = abs(V);              % hız negatif olsa bile Reynolds pozitif hesapla

Re = rho * V * cmac / mu;

if Re <= 1
    CD0 = NaN;
    return
end

Cf = 0.455 / ((log10(Re))^2.58);
CF = (1 + 0.2*M^2)^(-0.467);
IF = 1.20;
FF = 1 + 1.6*tc + 100*tc^4;
Swet = 2*S;

f = 1.1 * (Cf * CF * IF * FF * Swet);
CD0 = f / S;

end