clear;
clc;

%% =========================
% Timeline
% ==========================
t_start = 0;              % 07:00
t_lunch = 5*3600;         % 12:00
t_move = 10*3600 + 50*60; % 17:50
t_dinner = 11*3600;       % 18:00
t_stop = 15*3600;         % 22:00

cooking_duration = 3600;

simulation_time_outdoor = [t_start t_move];
simulation_time_indoor = [t_move t_stop];

%% =========================
% Geometry
% ==========================
r_in_pot = 0.15;
h_in_pot = 0.08;
wall_in_pot = 2e-3;

r_out_pot = 0.21;
h_out_pot = 0.10;
wall_out_pot = 2e-3;

ins_thickness = 0.26;

r_insulation = r_out_pot + ins_thickness;
h_insulation = h_out_pot + ins_thickness;

%% =========================
% Physical parameters (food, pot)
% ==========================
m_rice = 0.075;
m_water = m_rice * 2;
m_mix = m_water + m_rice;
m_pot = 2;

pot_absorptivity = 0.95; % rough/dull aluminum

c_water = 4148;
c_rice = 1700;
c_pot = 500;

c_mix = (m_water * c_water + m_rice * c_rice) / (m_water + m_rice);

k_water = 0.6096;
k_rice = 0.1472;
k_aluminum = 237;
k_st_steel = 15;
k_pot = k_aluminum;

k_mix = (m_water * k_water + m_rice * k_rice) / (m_water + m_rice);

%% =========================
% Other values
% ==========================
t_cooking = 3600; % cooking time (1 hour)

% Convection coefficients
h_free_outdoor = 5;
h_for_outdoor = 0;
h_free_indoor = 2;
h_for_indoor = 0;

T_sat = 100 + 273.15;
T_sat_c = 2;
T_sat_use = T_sat - T_sat_c;

% Solar data
load('raytracing_output.mat', 'hours', 'captured_power');
start_hour = 7;
solar_time = (hours - start_hour) * 3600;

sb_constant = 5.67e-8;
pot_eps = 0.95;           
T_amb = 20 + 273.15;

%% =========================
% Insulation material
% ==========================
m_ins = 3;
k_wool = 0.035;
c_wool = 1300;
k_ins = k_wool;
c_ins = c_wool;

%% =========================
% PCM properties (Erythritol from the paper)
% ==========================
m_pcm = 2;
rho_pcm = 1480;
T_m_pcm = 119 + 273.15;
L_pcm = 340e3;
delta_Tm_pcm = 10;               % mushy zone width (K)

% Specific heat: solid and liquid
c_pcm_solid = 1380;              % J/(kg·K)
c_pcm_liquid = 2760;             % (doubled as per user)

% Thermal conductivity: solid and liquid (pure conduction)
k_pcm_solid = 0.733;             % W/(m·K)
k_pcm_liquid = 0.326;

% Natural convection parameters (liquid phase)
pcm_beta = 0.00067;              % thermal expansion coefficient (1/K)
pcm_nu = 4e-6;                   % kinematic viscosity (m²/s)
pcm_alpha = 1.2e-7;              % thermal diffusivity (m²/s)
g = 9.81;

%% =========================
% Layer division
% 1 = rice+water, 2-6 = PCM, 7-11 = insulation
% ==========================
n_layers = 11;
n_pcm_layers = 5;
n_ins_layers = 5;

layer_temp = ones(n_layers,1) * T_amb;
T0 = layer_temp;

%% =========================
% Helper functions
% ==========================
function result = compute_layer_area(radius, height)
    result = 2 * pi * radius * height + pi * radius^2;
end

function [layer_area, layer_radius, layer_height, layer_thickness] = generate_layer_params(...
    r_in_pot, r_out_pot, r_insulation, h_in_pot, h_out_pot, h_insulation, ...
    n_pcm_layers, n_ins_layers, wall_in_pot)

    n_layers = 1 + n_pcm_layers + n_ins_layers;
    layer_area = zeros(n_layers,1);
    layer_radius = zeros(n_layers,1);
    layer_height = zeros(n_layers,1);
    layer_thickness = zeros(n_layers,1);

    % Layer 1: food
    layer_radius(1) = r_in_pot;
    layer_height(1) = h_in_pot;
    layer_area(1) = compute_layer_area(r_in_pot, h_in_pot);
    layer_thickness(1) = r_in_pot/2 + wall_in_pot;

    % PCM layers
    pcm_boundary = linspace(r_in_pot, r_out_pot, n_pcm_layers+1);
    pcm_heights = linspace(h_in_pot, h_out_pot, n_pcm_layers);
    for j = 1:n_pcm_layers
        i = 1 + j;
        r_inner = pcm_boundary(j);
        r_outer = pcm_boundary(j+1);
        r_mid = (r_inner + r_outer)/2;
        layer_radius(i) = r_mid;
        layer_height(i) = pcm_heights(j);
        layer_area(i) = compute_layer_area(r_mid, layer_height(i));
        layer_thickness(i) = r_outer - r_inner;
    end

    % Insulation layers
    ins_boundary = linspace(r_out_pot, r_insulation, n_ins_layers+1);
    ins_heights = linspace(h_out_pot, h_insulation, n_ins_layers);
    for j = 1:n_ins_layers
        i = 1 + n_pcm_layers + j;
        r_inner = ins_boundary(j);
        r_outer = ins_boundary(j+1);
        r_mid = (r_inner + r_outer)/2;
        layer_radius(i) = r_mid;
        layer_height(i) = ins_heights(j);
        layer_area(i) = compute_layer_area(r_mid, layer_height(i));
        layer_thickness(i) = r_outer - r_inner;
    end
end

function [layer_k_wall_incl, layer_thickness_wall_incl] = include_pot_walls(...
    layer_k, layer_thickness, n_pcm_layers, k_pot, wall_in_pot, wall_out_pot)

    layer_k_wall_incl = layer_k;
    layer_thickness_wall_incl = layer_thickness;

    % Inner wall to first PCM layer (layer 2)
    R_inner = wall_in_pot/k_pot + layer_thickness(2)/layer_k(2);
    layer_thickness_wall_incl(2) = wall_in_pot + layer_thickness(2);
    layer_k_wall_incl(2) = layer_thickness_wall_incl(2) / R_inner;

    % Outer wall to last PCM layer
    pcm_e = 1 + n_pcm_layers;
    R_outer = wall_out_pot/k_pot + layer_thickness(pcm_e)/layer_k(pcm_e);
    layer_thickness_wall_incl(pcm_e) = wall_out_pot + layer_thickness(pcm_e);
    layer_k_wall_incl(pcm_e) = layer_thickness_wall_incl(pcm_e) / R_outer;
end

function k_eff = pcm_thermal_conductivity(T, T_m, delta_Tm, k_solid, k_liquid, ...
                    layer_height, layer_radius, layer_thickness, ...
                    g, beta, nu, alpha)
    % Returns effective k (W/m·K) considering natural convection in liquid
    T_lower = T_m - delta_Tm/2;
    T_upper = T_m + delta_Tm/2;

    if T <= T_lower
        k_eff = k_solid;
    elseif T >= T_upper
        % Liquid phase: natural convection enhancement
        L_char = layer_height;                 % characteristic length (vertical)
        delta_T = T - T_m;                     % superheat (K)
        if delta_T <= 0
            k_eff = k_liquid;
        else
            Ra = g * beta * delta_T * L_char^3 / (nu * alpha);
            % Simple Nusselt correlation for vertical cylinder / annulus
            Nu = 1 + 0.18 * Ra^(0.29);
            Nu = min(Nu, 10);                  % limit
            k_eff = k_liquid * Nu;
        end
    else
        % Mushy zone: linear interpolation
        f = (T - T_lower) / delta_Tm;
        k_eff = k_solid * (1-f) + k_liquid * f;
    end
end

function layer_c_ef = pcm_heat_capacity_change(layer_temp, layer_c_solid, layer_c_liquid, ...
                        n_pcm_layers, T_m_pcm, delta_Tm_pcm, L_pcm)
    % Effective heat capacity: solid/liquid + latent spike
    layer_c_ef = layer_c_solid;   % base, then overwrite PCM layers
    T_lower = T_m_pcm - delta_Tm_pcm/2;
    T_upper = T_m_pcm + delta_Tm_pcm/2;

    for i = 2:(1+n_pcm_layers)
        T = layer_temp(i);
        if T <= T_lower
            layer_c_ef(i) = layer_c_solid(i);
        elseif T >= T_upper
            layer_c_ef(i) = layer_c_liquid(i);
        else
            c_avg = (layer_c_solid(i) + layer_c_liquid(i)) / 2;
            layer_c_ef(i) = c_avg + L_pcm / delta_Tm_pcm;
        end
    end
end

function dTidt = layer_equation(i, i_e, Qs, Qr, Qc, Qco, Qv, layer_temp, ...
                layer_area, layer_k, layer_thickness, layer_mass, layer_c)

    if i == 1
        Q_cond = layer_area(i) * layer_k(i) / layer_thickness(i) * (layer_temp(i+1) - layer_temp(i));
        dTidt = (Q_cond - (Qco + Qv)) / (layer_mass(i) * layer_c(i));
    elseif i == i_e
        Q_cond = layer_area(i) * layer_k(i) / layer_thickness(i) * (layer_temp(i-1) - layer_temp(i));
        dTidt = (Qs - (Qr + Qc) + Q_cond) / (layer_mass(i) * layer_c(i));
    else
        Q_cond_left  = layer_area(i) * layer_k(i) / layer_thickness(i) * (layer_temp(i-1) - layer_temp(i));
        Q_cond_right = layer_area(i) * layer_k(i) / layer_thickness(i) * (layer_temp(i+1) - layer_temp(i));
        dTidt = (Q_cond_left + Q_cond_right) / (layer_mass(i) * layer_c(i));
    end
end

function Qr = radiation_loss(i_e, layer_temp, layer_area, T_amb, pot_eps, sb_constant)
    Qr = layer_area(i_e) * pot_eps * sb_constant * (layer_temp(i_e)^4 - T_amb^4);
end

function Qc = convection_loss(i_e, layer_temp, layer_area, T_amb, h_free, h_for)
    h_total = (h_free^3 + h_for^3)^(1/3);
    Qc = layer_area(i_e) * h_total * (layer_temp(i_e) - T_amb);
end

function Qco = cooking_loss(layer_temp, T_amb, c_mix, m_mix, t_cooking)
    Qco = (c_mix * m_mix / t_cooking) * (layer_temp(1) - T_amb);
end

function Qv = evaporation_loss(layer_temp, T_sat_use)
    Kv = 0.1;
    s = 2;
    if layer_temp(1) >= T_sat_use
        Qv = Kv * (layer_temp(1) - T_sat_use)^s;
    else
        Qv = 0;
    end
end

function dTdt = full_model(t, layer_temp, i_e, solar_time, captured_power, ...
                pot_absorptivity, layer_area, layer_k_base, layer_thickness, ...
                layer_mass, layer_c_solid, layer_c_liquid, T_amb, pot_eps, ...
                sb_constant, h_free, h_for, c_mix, m_mix, t_cooking, T_sat_use, ...
                n_pcm_layers, T_m_pcm, delta_Tm_pcm, L_pcm, cooking_time_points, ...
                cooking_duration, pcm_params)

    % pcm_params = struct with: g, beta, nu, alpha, k_solid, k_liquid, layer_heights, layer_radii
    dTdt = zeros(size(layer_temp));

    % Solar input
    Qs_captured = interp1(solar_time, captured_power, t, 'linear', 0);
    Qs = pot_absorptivity * Qs_captured;

    % Outer boundary losses (radiation + convection)
    Qr = radiation_loss(i_e, layer_temp, layer_area, T_amb, pot_eps, sb_constant);
    Qc = convection_loss(i_e, layer_temp, layer_area, T_amb, h_free, h_for);

    % Cooking active?
    is_cooking = false;
    for j = 1:length(cooking_time_points)
        if t >= cooking_time_points(j) && t <= cooking_time_points(j) + cooking_duration
            is_cooking = true;
            break;
        end
    end
    if is_cooking
        Qco = cooking_loss(layer_temp, T_amb, c_mix, m_mix, t_cooking);
    else
        Qco = 0;
    end

    Qv = evaporation_loss(layer_temp, T_sat_use);

    % Effective heat capacity (includes latent spike)
    layer_c_ef = pcm_heat_capacity_change(layer_temp, layer_c_solid, layer_c_liquid, ...
                    n_pcm_layers, T_m_pcm, delta_Tm_pcm, L_pcm);

    % Build effective conductivity array (temperature‑dependent for PCM)
    layer_k_eff = layer_k_base;
    for i = 2:(1+n_pcm_layers)
        k_eff = pcm_thermal_conductivity(layer_temp(i), T_m_pcm, delta_Tm_pcm, ...
                pcm_params.k_solid, pcm_params.k_liquid, ...
                pcm_params.layer_height(i), pcm_params.layer_radius(i), ...
                layer_thickness(i), pcm_params.g, pcm_params.beta, ...
                pcm_params.nu, pcm_params.alpha);
        layer_k_eff(i) = k_eff;
    end

    % Compute temperature derivatives
    for i = 1:i_e
        dTdt(i) = layer_equation(i, i_e, Qs, Qr, Qc, Qco, Qv, layer_temp, ...
                    layer_area, layer_k_eff, layer_thickness, layer_mass, layer_c_ef);
    end
end

%% =========================
% Generate geometry and layer properties
% ==========================
[layer_area, layer_radius, layer_height, layer_thickness] = generate_layer_params(...
    r_in_pot, r_out_pot, r_insulation, h_in_pot, h_out_pot, h_insulation, ...
    n_pcm_layers, n_ins_layers, wall_in_pot);

% Adjust food layer thickness
layer_thickness(1) = wall_in_pot + r_in_pot/10;

%% =========================
% Base thermal conductivities (will be updated during simulation)
% ==========================
layer_k_base = zeros(n_layers,1);
layer_k_base(1) = k_mix;
layer_k_base(2:(1+n_pcm_layers)) = k_pcm_solid;   % initial solid value
layer_k_base((2+n_pcm_layers):n_layers) = k_ins;

% Include pot walls (merges into PCM layers 2 and 6)
[layer_k_base, layer_thickness] = include_pot_walls(layer_k_base, layer_thickness, ...
    n_pcm_layers, k_pot, wall_in_pot, wall_out_pot);

%% =========================
% Layer masses
% ==========================
layer_mass = zeros(n_layers,1);
layer_mass(1) = m_water;
layer_mass(2:(1+n_pcm_layers)) = m_pcm / n_pcm_layers;
layer_mass((2+n_pcm_layers):n_layers) = m_ins / n_ins_layers;

%% =========================
% Specific heat arrays (solid and liquid for PCM)
% ==========================
layer_c_solid = zeros(n_layers,1);
layer_c_liquid = zeros(n_layers,1);

layer_c_solid(1) = c_water;
layer_c_liquid(1) = c_water;   % same for food

layer_c_solid(2:(1+n_pcm_layers)) = c_pcm_solid;
layer_c_liquid(2:(1+n_pcm_layers)) = c_pcm_liquid;

layer_c_solid((2+n_pcm_layers):n_layers) = c_ins;
layer_c_liquid((2+n_pcm_layers):n_layers) = c_ins;

%% =========================
% PCM natural convection parameters (passed to ODE)
% ==========================
pcm_params = struct();
pcm_params.g = g;
pcm_params.beta = pcm_beta;
pcm_params.nu = pcm_nu;
pcm_params.alpha = pcm_alpha;
pcm_params.k_solid = k_pcm_solid;
pcm_params.k_liquid = k_pcm_liquid;
pcm_params.layer_height = layer_height;
pcm_params.layer_radius = layer_radius;

%% =========================
% Outdoor simulation (07:00 - 17:50)
% ==========================
i_e_outdoor = 1 + n_pcm_layers;
cooking_time_points_outdoor = [t_lunch];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', 60);

[t_out, T_out] = ode45(@(t,T) full_model(t, T, i_e_outdoor, solar_time, captured_power, ...
    pot_absorptivity, layer_area, layer_k_base, layer_thickness, layer_mass, ...
    layer_c_solid, layer_c_liquid, T_amb, pot_eps, sb_constant, h_free_outdoor, ...
    h_for_outdoor, c_mix, m_mix, t_cooking, T_sat_use, n_pcm_layers, T_m_pcm, ...
    delta_Tm_pcm, L_pcm, cooking_time_points_outdoor, cooking_duration, pcm_params), ...
    simulation_time_outdoor, T0, options);

%% =========================
% CORRECTED PCM ENERGY STORED (outdoor phase)
% ==========================
% Function to compute effective c for a given PCM temperature (used in ODE)
c_eff_pcm = @(T_K) (T_K <= T_m_pcm - delta_Tm_pcm/2) * c_pcm_solid + ...
                   (T_K >= T_m_pcm + delta_Tm_pcm/2) * c_pcm_liquid + ...
                   ((T_K > T_m_pcm - delta_Tm_pcm/2) & (T_K < T_m_pcm + delta_Tm_pcm/2)) * ...
                   ( (c_pcm_solid + c_pcm_liquid)/2 + L_pcm/delta_Tm_pcm );

% Compute stored energy by integrating m * c_eff * dT for each PCM layer
E_pcm_stored = 0;
for i = 2:(1+n_pcm_layers)
    T_history = T_out(:, i);
    dTdt_history = gradient(T_history, t_out);  % approximate derivative
    c_eff_history = arrayfun(c_eff_pcm, T_history);
    power_history = layer_mass(i) * c_eff_history .* dTdt_history;
    % Only integrate positive power (energy entering the layer)
    power_history = max(power_history, 0);
    E_layer = trapz(t_out, power_history);
    E_pcm_stored = E_pcm_stored + E_layer;
end
E_pcm_stored_MJ = E_pcm_stored / 1e6;

%% =========================
% Indoor simulation (17:50 - 22:00)
% ==========================
i_e_indoor = n_layers;
solar_time_indoor = [t_move t_stop];
captured_power_indoor = [0 0];
cooking_time_points_indoor = [t_dinner];

T0_indoor = T_amb * ones(n_layers,1);
T0_indoor(1:i_e_outdoor) = T_out(end, 1:i_e_outdoor)';

[t_in, T_in] = ode45(@(t,T) full_model(t, T, i_e_indoor, solar_time_indoor, captured_power_indoor, ...
    pot_absorptivity, layer_area, layer_k_base, layer_thickness, layer_mass, ...
    layer_c_solid, layer_c_liquid, T_amb, pot_eps, sb_constant, h_free_indoor, ...
    h_for_indoor, c_mix, m_mix, t_cooking, T_sat_use, n_pcm_layers, T_m_pcm, ...
    delta_Tm_pcm, L_pcm, cooking_time_points_indoor, cooking_duration, pcm_params), ...
    simulation_time_indoor, T0_indoor, options);

%% =========================
% CORRECTED PCM ENERGY RELEASED (indoor phase)
% ==========================
E_pcm_released = 0;
for i = 2:(1+n_pcm_layers)
    T_history = T_in(:, i);
    dTdt_history = gradient(T_history, t_in);
    c_eff_history = arrayfun(c_eff_pcm, T_history);
    power_history = layer_mass(i) * c_eff_history .* dTdt_history;
    % Power leaving the layer is negative; integrate -power (positive)
    power_history = -min(power_history, 0);
    E_layer = trapz(t_in, power_history);
    E_pcm_released = E_pcm_released + E_layer;
end
E_pcm_released_MJ = E_pcm_released / 1e6;

%% =========================
% Plotting – original selected layers
% ==========================
t_all = [t_out; t_in];
T_all = [T_out; T_in];

selected_layers_full = [1,2,6,7,11];
figure;
plot(t_all/3600+7, T_all(:,selected_layers_full)-273.15, 'LineWidth',0.8);
hold on;
xlabel('Time of day [h]');
ylabel('Temperature [°C]');
title('Full-day heat transfer simulation (selected layers)');
legend('Layer1: food','Layer2: PCM inner','Layer6: PCM outer','Layer7: insulation inner','Layer11: insulation outer');
grid on;
xline(7+t_lunch/3600,'--','Lunch','HandleVisibility','off');
xline(7+t_move/3600,'--','Move indoors','HandleVisibility','off');
xline(7+t_dinner/3600,'--','Dinner','HandleVisibility','off');
yl = ylim;
patch([12 13 13 12],[yl(1) yl(1) yl(2) yl(2)],[0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
patch([7+t_dinner/3600, 7+(t_dinner+cooking_duration)/3600, 7+(t_dinner+cooking_duration)/3600, 7+t_dinner/3600], ...
    [yl(1) yl(1) yl(2) yl(2)],[0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
uistack(findobj(gca,'Type','line'),'top');
hold off;

figure;
plot((t_in-t_move)/60, T_in(:,selected_layers_full)-273.15, 'LineWidth',0.8);
hold on;
xlabel('Time after moving indoors [min]');
ylabel('Temperature [°C]');
title('Indoor simulation (selected layers)');
legend('Layer1: food','Layer2: PCM inner','Layer6: PCM outer','Layer7: insulation inner','Layer11: insulation outer');
grid on;
xline((t_dinner-t_move)/60,'--','Dinner','HandleVisibility','off');
yl = ylim;
patch([(t_dinner-t_move)/60, (t_dinner+cooking_duration-t_move)/60, ...
       (t_dinner+cooking_duration-t_move)/60, (t_dinner-t_move)/60], ...
      [yl(1) yl(1) yl(2) yl(2)],[0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
uistack(findobj(gca,'Type','line'),'top');
hold off;

%% =========================
% NEW: Plot all layers with color map
% ==========================
figure;
cmap = jet(n_layers);
for i = 1:n_layers
    plot(t_all/3600+7, T_all(:,i)-273.15, 'Color', cmap(i,:), 'LineWidth', 0.8);
    hold on;
end
xlabel('Time of day [h]');
ylabel('Temperature [°C]');
title('Full-day simulation – all 11 layers');
colormap(cmap);
c = colorbar;
c.Label.String = 'Layer number';
c.Ticks = 1:n_layers;
grid on;
% Add cooking and move lines
xline(7+t_lunch/3600,'--','Lunch','HandleVisibility','off');
xline(7+t_move/3600,'--','Move indoors','HandleVisibility','off');
xline(7+t_dinner/3600,'--','Dinner','HandleVisibility','off');
yl = ylim;
patch([12 13 13 12],[yl(1) yl(1) yl(2) yl(2)],[0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
patch([7+t_dinner/3600, 7+(t_dinner+cooking_duration)/3600, 7+(t_dinner+cooking_duration)/3600, 7+t_dinner/3600], ...
    [yl(1) yl(1) yl(2) yl(2)],[0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
hold off;

figure;
for i = 1:n_layers
    plot((t_in-t_move)/60, T_in(:,i)-273.15, 'Color', cmap(i,:), 'LineWidth', 0.8);
    hold on;
end
xlabel('Time after moving indoors [min]');
ylabel('Temperature [°C]');
title('Indoor simulation – all 11 layers');
colormap(cmap);
c = colorbar;
c.Label.String = 'Layer number';
c.Ticks = 1:n_layers;
grid on;
xline((t_dinner-t_move)/60,'--','Dinner','HandleVisibility','off');
yl = ylim;
patch([(t_dinner-t_move)/60, (t_dinner+cooking_duration-t_move)/60, ...
       (t_dinner+cooking_duration-t_move)/60, (t_dinner-t_move)/60], ...
      [yl(1) yl(1) yl(2) yl(2)],[0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
hold off;

%% =========================
% Energy and power analysis
% ==========================
fprintf('\n========== THERMAL POWER AND ENERGY ANALYSIS ==========\n');

% Outdoor phase
Qs_outdoor = pot_absorptivity * interp1(solar_time, captured_power, t_out, 'linear', 0);
Qr_outdoor = layer_area(i_e_outdoor) * pot_eps * sb_constant * (T_out(:,i_e_outdoor).^4 - T_amb^4);
h_total_outdoor = (h_free_outdoor^3 + h_for_outdoor^3)^(1/3);
Qc_outdoor = layer_area(i_e_outdoor) * h_total_outdoor * (T_out(:,i_e_outdoor) - T_amb);
Qloss_outdoor = Qr_outdoor + Qc_outdoor;
Quseful_outdoor = max(Qs_outdoor - Qloss_outdoor, 0);

E_solar_absorbed_outdoor = trapz(t_out, Qs_outdoor)/1e6;
E_loss_outdoor = trapz(t_out, Qloss_outdoor)/1e6;
E_useful_outdoor = trapz(t_out, Quseful_outdoor)/1e6;

fprintf('\n--- OUTDOOR PHASE (07:00 to 17:50) ---\n');
fprintf('Peak solar absorbed power:          %.1f W\n', max(Qs_outdoor));
fprintf('Mean solar absorbed power:          %.1f W\n', mean(Qs_outdoor(Qs_outdoor>0)));
fprintf('Peak thermal loss (Qr+Qc):          %.1f W\n', max(Qloss_outdoor));
fprintf('Mean thermal loss:                  %.1f W\n', mean(Qloss_outdoor));
fprintf('Peak useful power (Qs-losses):      %.1f W\n', max(Quseful_outdoor));
fprintf('Mean useful power:                  %.1f W\n', mean(Quseful_outdoor(Quseful_outdoor>0)));
fprintf('Total solar absorbed:               %.4f MJ\n', E_solar_absorbed_outdoor);
fprintf('Total loss energy:                  %.4f MJ\n', E_loss_outdoor);
fprintf('Useful energy (absorbed-loss):      %.4f MJ\n', E_useful_outdoor);
fprintf('Energy stored in PCM (correct):     %.4f MJ\n', E_pcm_stored_MJ);

% Indoor phase
Qs_indoor = zeros(size(t_in));
Qr_indoor = layer_area(i_e_indoor) * pot_eps * sb_constant * (T_in(:,i_e_indoor).^4 - T_amb^4);
h_total_indoor = (h_free_indoor^3 + h_for_indoor^3)^(1/3);
Qc_indoor = layer_area(i_e_indoor) * h_total_indoor * (T_in(:,i_e_indoor) - T_amb);
Qloss_indoor = Qr_indoor + Qc_indoor;

Qconduction_to_food = layer_area(1) * layer_k_base(1) / layer_thickness(1) .* (T_in(:,2) - T_in(:,1));
Qconduction_to_food = max(Qconduction_to_food, 0);

E_loss_indoor = trapz(t_in, Qloss_indoor)/1e6;
E_to_food_indoor = trapz(t_in, Qconduction_to_food)/1e6;

T_food_indoor_C = T_in(:,1) - 273.15;
dt = diff(t_in);
above = T_food_indoor_C(1:end-1) >= 100;
time_above_100 = sum(dt(above)) / 60;

fprintf('\n--- INDOOR PHASE (17:50 to 22:00) ---\n');
fprintf('Peak thermal loss (Qr+Qc):          %.1f W\n', max(Qloss_indoor));
fprintf('Mean thermal loss:                  %.1f W\n', mean(Qloss_indoor));
fprintf('Peak conduction to food:            %.1f W\n', max(Qconduction_to_food));
fprintf('Mean conduction to food:            %.1f W\n', mean(Qconduction_to_food(Qconduction_to_food>0)));
fprintf('Total loss energy (indoors):        %.4f MJ\n', E_loss_indoor);
fprintf('Energy delivered to food:           %.4f MJ\n', E_to_food_indoor);
fprintf('PCM energy released (correct):      %.4f MJ\n', E_pcm_released_MJ);
fprintf('Time food above 100°C:              %.1f min\n', time_above_100);

fprintf('\n--- FULL DAY SUMMARY ---\n');
fprintf('Total solar absorbed:               %.4f MJ\n', E_solar_absorbed_outdoor);
fprintf('Total loss (outdoor+indoor):        %.4f MJ\n', E_loss_outdoor+E_loss_indoor);
fprintf('PCM stored energy:                  %.4f MJ\n', E_pcm_stored_MJ);
fprintf('PCM released energy:                %.4f MJ\n', E_pcm_released_MJ);
fprintf('PCM efficiency (released/stored):   %.1f%%\n', 100*E_pcm_released_MJ/max(E_pcm_stored_MJ,1e-9));
fprintf('Food above 100°C:                   %.1f min\n', time_above_100);
fprintf('========================================================\n');