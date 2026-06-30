clc; clear; close all;

%%
%initial guessler
%%
g = 9.81;
m = 3;                      %kg
w = m*g;

hinitial = 15000;           %başlangıç irtifası metre
xinitial = 0;               %başlangıç konumu
vinitial = 0.1;               % m/s
gama_deg = -5;              % deg
gamainitial = deg2rad(gama_deg);

tinitial = 0;               %ilk zaman
dt = 0.05;                  %zaman adımı
tmax = 150000;               %sonsuz döngüye karşı güvenlik sınırı

v_switch = 5;               %bu hızdan sonra normal glide denklemlerine geçilecek

b = 10;                     %metre mecburi
AR = 10;
S = b^2 / AR;               %area planform
chord = S / b;              % rectangular wing chord
e = 1 - 0.045*(AR^0.68);    %oswald efficincy tapered kullanırsan değişecek
kapa = 1/(pi*e*AR);         %CDi kapası kcl^2

alpha0_deg = 4.7;             % sabit alpha
alpha0_rad = deg2rad(alpha0_deg);

%%
%kanat noktaları
%%
coord = load('naca0012.txt');
coord = flipud(coord);

x = coord(:,1);
y = coord(:,2);

%burdan (t/c) oranı geldi AİRFOİL için
tc = (max(y) - min(y)) / (max(x) - min(x));

xmom = 0.25;
alpha_panel_deg = [-4 -2 0 2 4 6];
Cl_panel = zeros(size(alpha_panel_deg));

for i = 1:length(alpha_panel_deg)
    alpha_tmp = alpha_panel_deg(i);
    [~,~,~,Cl,~,~] = panel_ORJ(x, y, alpha_tmp, xmom);
    Cl_panel(i) = Cl;
end

pfit = polyfit(alpha_panel_deg, Cl_panel, 1);
clalpha_incomp_deg = pfit(1);
clalpha_incomp_rad = clalpha_incomp_deg * 180/pi;
alphaL0_deg = -pfit(2)/pfit(1);
alphaL0_rad = deg2rad(alphaL0_deg);

clalpha_theory = 6.28 + 4.7*tc;
kappa = clalpha_incomp_rad / clalpha_theory;

fprintf('clalpha_incomp_deg = %.6f\n', clalpha_incomp_deg);
fprintf('clalpha_incomp_rad = %.6f\n', clalpha_incomp_rad);
fprintf('alphaL0_deg = %.6f\n', alphaL0_deg);

phase2_started = false;

%%
% KAYDEDİLECEK ALAN
%%
xforgraph = [];
hforgraph = [];
vforgraph = [];
gamaforgraph = [];
tforgraph = [];
ClLforgraph = [];
cDforgraph = [];
Lforgraph = [];
Dforgraph = [];
phaseforgraph = [];
Mforgraph = [];


%%
%burda bir yerde döngüye girmesi lazım
%%
iter = 0;
itermax = ceil(tmax/dt);

while hinitial > 0 && iter < itermax
    iter = iter + 1;

    % atmosphere
    [T, rho, mu, a] = atmosphere(hinitial);

    if vinitial <= 0
        warning('Velocity became zero or negative. Solution stopped.');
        break
    end

    M = vinitial / a;
    if M >= 1
        warning('Mach >= 1. Solution stopped.');
        break
    end

    % wing lift slope
    radicand = 1 + (AR/(2*kappa))^2 * (1 - M^2);
    if radicand <= 0 || ~isfinite(radicand)
        warning('a0 radicand invalid. Solution stopped.');
        break
    end

    a0 = (pi*AR) / (1 + sqrt(radicand));

    % aerodynamics
    clL = a0 * (alpha0_rad - alphaL0_rad);

    [cd0, Re] = cd0_hull_wing(rho, max(vinitial,1e-6), mu, S, chord, tc, M);
    if ~isfinite(cd0) || cd0 <= 0
        warning('CD0 invalid. Solution stopped.');
        break
    end

    cD = cd0 + kapa*clL^2;

    % forces
    L = 0.5 * rho * vinitial^2 * S * clL;
    D = 0.5 * rho * vinitial^2 * S * cD;

    if ~isfinite(L) || ~isfinite(D) || ~isfinite(clL) || ~isfinite(cD)
        warning('Aerodynamic values invalid. Solution stopped.');
        break
    end

    
    % RELEASE / FALLING PHASE
   
    if ~phase2_started
    if vinitial < v_switch
        x_ave = vinitial * cos(gamainitial);
        h_ave = vinitial * sin(gamainitial);
        v_ave = -D/m - g*sin(gamainitial);
        gama_ave = 0;
        phaseforgraph(end+1) = 1;
    else
        phase2_started = true;

        x_ave = vinitial * cos(gamainitial);
        h_ave = vinitial * sin(gamainitial);
        v_ave = -D/m - g*sin(gamainitial);
        gama_ave = (L - w*cos(gamainitial)) / (m*vinitial);
        phaseforgraph(end+1) = 2;
    end
else
    % NORMAL GLIDE PHASE
    % Bu fazda artık flight mechanics denklemleri kullanılıyor.
    % Gamma denklemi yok edilmiyor.
    x_ave = vinitial * cos(gamainitial);
    h_ave = vinitial * sin(gamainitial);
    v_ave = -D/m - g*sin(gamainitial);
    gama_ave = (L - w*cos(gamainitial)) / (m*vinitial);
    phaseforgraph(end+1) = 2;
end
    %euler method
    xinitialnew = xinitial + x_ave * dt;
    hinitialnew = hinitial + h_ave * dt;
    vinitialnew = vinitial + v_ave * dt;
    gamainitialnew = gamainitial + gama_ave * dt;
    tinitialnew = tinitial + dt;

    % history
    xforgraph(end+1) = xinitialnew;
    hforgraph(end+1) = hinitialnew;
    vforgraph(end+1) = vinitialnew;
    gamaforgraph(end+1) = gamainitialnew;
    tforgraph(end+1) = tinitialnew;
    ClLforgraph(end+1) = clL;
    cDforgraph(end+1) = cD;
    Lforgraph(end+1) = L;
    Dforgraph(end+1) = D;
    Mforgraph(end+1) = M;

    %yeni değerleri değiştir bir sonrakinin eskisi yap
    xinitial = xinitialnew;
    hinitial = hinitialnew;
    vinitial = vinitialnew;
    gamainitial = gamainitialnew;
    tinitial = tinitialnew;

   % gamma pozitif olmasın
if gamainitial > deg2rad(-0.1)
    gamainitial = deg2rad(-0.1);
end

% gamma aşırı dik dalışa gitmesin
if gamainitial < deg2rad(-85)
    gamainitial = deg2rad(-85);
end

end

%%
% sonuç
%%
range_km = xinitial / 1000;
time_hr  = tinitial / 3600;

fprintf('\nTotal range = %.3f km\n', range_km);

%% Plots
figure;
plot(xforgraph/1000, hforgraph/1000, 'LineWidth', 1.6);
grid on;
xlabel('Range [km]');
ylabel('Altitude [km]');
title('Glide Trajectory');

figure;
plot(tforgraph, vforgraph, 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Velocity History');

figure;
plot(tforgraph, rad2deg(gamaforgraph), 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('\gamma [deg]');
title('Flight Path Angle History');

figure;
plot(tforgraph, ClLforgraph, 'LineWidth', 1.6); hold on;
plot(tforgraph, cDforgraph, 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('Coefficient');
legend('C_L','C_D','Location','best');
title('Aerodynamic Coefficients');

figure;
plot(tforgraph, Lforgraph, 'LineWidth', 1.6); hold on;
plot(tforgraph, Dforgraph, 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('Force [N]');
legend('L','D','Location','best');
title('Aerodynamic Forces');

figure;
stairs(tforgraph(1:length(phaseforgraph)), phaseforgraph, 'LineWidth', 1.6);
grid on;
xlabel('Time [s]');
ylabel('Phase');
title('1 = Low-speed phase, 2 = Full dynamics');
ylim([0.5 2.5]);

%% Plots

figure
plot(xforgraph/1000, hforgraph/1000, 'LineWidth', 1.5)
xlabel('Range [km]')
ylabel('Altitude [km]')
grid on
title('Glide Trajectory')

figure
plot(tforgraph, vforgraph, 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on
title('Velocity vs Time')

figure
plot(tforgraph, rad2deg(gamaforgraph), 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('\gamma [deg]')
grid on
title('Flight Path Angle vs Time')

figure
plot(tforgraph, ClLforgraph, 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('C_L')
grid on
title('Lift Coefficient vs Time')

figure
plot(tforgraph, cDforgraph, 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('C_D')
grid on
title('Drag Coefficient vs Time')

figure
plot(tforgraph, Lforgraph./w, 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('L/W')
grid on
title('Lift to Weight Ratio vs Time')

figure
plot(tforgraph, Dforgraph./w, 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('D/W')
grid on
title('Drag to Weight Ratio vs Time')

figure
plot(tforgraph, Mforgraph, 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('Mach Number')
grid on
title('Mach Number vs Time')