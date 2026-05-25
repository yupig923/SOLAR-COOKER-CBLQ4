%% =========================================================
%  SOLAR COOKER - RAY TRACING WITH PVGIS IRRADIANCE
%  Group 2, 4CBLW025
%
%   - Sun-tracking dish (fixed vertical dish gives zero output)
%   - Imperfect tracking: manual repositioning every N minutes
%   - Surface slope errors: Mylar wrinkle
%   - Receiver shading: pot blocks central rays
%   - Sensitivity study across tracking interval + slope error
%
%  DATA: PVGIS-SARAH3, Eindhoven 2023
%  Columns: time, Gb(i), Gd(i), Gr(i), H_sun, T2m, WS10m, Int
% =========================================================

clc; clear; close all;

CSV_PATH    = 'pvgis_hourly.csv';
DAY_OF_YEAR = 152;              % 152 = June 1

% Dish
f                    = 0.5;    % focal length [m]
radius               = 1.0;    % dish half-width [m]
reflectivity         = 0.95;   % Mylar reflectivity [-]
receiver_radius      = 0.105;   % receiver (pot base) radius [m]
N                    = 5000;   % Monte Carlo rays
sun_spread_deg       = 0.27;   % solar disc half-angle [deg]

% Improvement parameters
reposition_interval  = 15;     % manual tracking interval [min]
slope_error_deg      = 0.5;    % Mylar surface slope error [deg]

[hours, Gbi, elev_deg] = load_pvgis(CSV_PATH, DAY_OF_YEAR);
fprintf('Loaded %d points. Peak Gb(i): %.0f W/m2\n', ...
    numel(hours), max(Gbi));

%% Compute the dish geometry
x_plot = linspace(-radius, radius, 500);
y_plot = x_plot.^2 / (4*f);

hit_x = linspace(-radius, radius, N);
hit_y = hit_x.^2 / (4*f);

% Ideal surface normals
dy_dx         = hit_x ./ (2*f);
normals_ideal = [-dy_dx; ones(size(dy_dx))];
normals_ideal = normals_ideal ./ ...
    sqrt(normals_ideal(1,:).^2 + normals_ideal(2,:).^2);

receiver_center = [0; f];
dish_area       = pi * radius^2;
power_per_ray   = dish_area / N;

% Receiver shading: rays blocked by pot shadow
blocked = abs(hit_x) < receiver_radius;

% Pre-compute tracking parameters
deg_per_interval = 15 * (reposition_interval / 60);
slope_error_rad  = deg2rad(slope_error_deg);
spread_rad       = deg2rad(sun_spread_deg);

%% Time loop
n_steps            = numel(hours);
captured_power     = zeros(1, n_steps);
incident_power     = zeros(1, n_steps);
optical_efficiency = zeros(1, n_steps);
DNI_vec            = zeros(1, n_steps);

for t_idx = 1:n_steps
    elev    = elev_deg(t_idx);
    Gbi_now = Gbi(t_idx);

    if elev <= 5 || Gbi_now < 10
        continue;
    end

    elev_rad = deg2rad(elev);

    %Imperfect tracking, snap to last repositioning step
    elev_snapped       = floor(elev / deg_per_interval) * deg_per_interval;
    tracking_error_rad = deg2rad(elev - elev_snapped);

    % Implement surface slope errors, redraw normals each step
    phi          = slope_error_rad * randn(1, N);
    cos_phi      = cos(phi);
    sin_phi      = sin(phi);
    normals      = zeros(2, N);
    normals(1,:) = cos_phi .* normals_ideal(1,:) ...
                 - sin_phi .* normals_ideal(2,:);
    normals(2,:) = sin_phi .* normals_ideal(1,:) ...
                 + cos_phi .* normals_ideal(2,:);

    % DNI from horizontal beam irradiance
    DNI_now = Gbi_now / sin(elev_rad);
    DNI_vec(t_idx) = DNI_now;

    % Incident power (corrected for receiver shadow)
    shadow_fraction          = (2 * receiver_radius) / (2 * radius);
    P_incident               = DNI_now * dish_area * (1 - shadow_fraction^2);
    incident_power(t_idx)    = P_incident;

    % Ray directions: tracking error + solar disc spread
    theta = spread_rad * randn(1, N) + tracking_error_rad;
    dir   = [sin(theta); -cos(theta)];

    % Energy per ray
    ray_energy = (DNI_now * power_per_ray) * ones(1, N);

    % Receiver shading
    ray_energy(blocked) = 0;

    % Reflect and check capture
    captured = 0;
    for i = 1:N
        if ray_energy(i) == 0, continue; end

        d = dir(:, i);
        n = normals(:, i);

        r = d - 2 * dot(d, n) * n;
        ray_energy(i) = ray_energy(i) * reflectivity;

        t_param = linspace(0, 2, 80);
        rx      = hit_x(i) + r(1) * t_param;
        ry      = hit_y(i) + r(2) * t_param;
        dist    = sqrt((rx - receiver_center(1)).^2 + ...
                       (ry - receiver_center(2)).^2);

        if any(dist < receiver_radius)
            captured = captured + ray_energy(i);
        end
    end

    captured_power(t_idx) = captured;
    if P_incident > 0
        optical_efficiency(t_idx) = captured / P_incident;
    end
end

%% Results
total_energy_MJ = trapz(hours * 3600, captured_power) / 1e6;
daytime         = optical_efficiency > 0;

fprintf('\n===== RESULTS =====\n');
fprintf('Day of year:          %d\n',   DAY_OF_YEAR);
fprintf('Reposition interval:  %d min\n', reposition_interval);
fprintf('Surface slope error:  %.2f deg\n', slope_error_deg);
fprintf('Peak DNI:             %.1f W/m2\n', max(DNI_vec));
fprintf('Peak captured power:  %.1f W\n',    max(captured_power));
fprintf('Total energy (day):   %.3f MJ\n',   total_energy_MJ);
fprintf('Mean optical eff:     %.1f%%\n', ...
    100 * mean(optical_efficiency(daytime)));

save('raytracing_output.mat', 'hours', 'captured_power');

%% Main plots
figure;
yyaxis left
bar(hours, Gbi, 0.6, 'FaceColor',[0.8 0.85 0.95],'EdgeColor','none');
ylabel('G_b(i) horizontal beam irradiance (W/m^2)');
yyaxis right
plot(hours, captured_power, 'r-o', 'LineWidth', 2, 'MarkerSize', 5);
ylabel('Captured power at receiver (W)');
xlabel('Hour of day (solar time)');
title(sprintf('PVGIS irradiance and captured power -- DoY %d, Eindhoven', ...
    DAY_OF_YEAR));
legend('Gb(i) from PVGIS', 'Captured power');
grid on;

figure;
plot(hours(daytime), 100 * optical_efficiency(daytime), ...
    'g-^', 'LineWidth', 2, 'MarkerSize', 6);
yline(95, '--k', 'Reflectivity limit (95%)');
xlabel('Hour of day');
ylabel('Optical efficiency (%)');
title(sprintf('Optical efficiency -- reposition=%d min, slope=%.1f deg', ...
    reposition_interval, slope_error_deg));
ylim([0 105]); grid on;

%% Sensitivity study
intervals    = [5, 10, 15, 20, 30];
slope_errors = [0, 0.25, 0.5, 1.0, 2.0];

fprintf('\nRunning sensitivity study...\n');
mean_eff_grid = zeros(numel(intervals), numel(slope_errors));

for ki = 1:numel(intervals)
    for ks = 1:numel(slope_errors)
        mean_eff_grid(ki,ks) = run_day(hours, Gbi, elev_deg, ...
            intervals(ki), slope_errors(ks), ...
            f, radius, reflectivity, receiver_radius, ...
            N, sun_spread_deg, normals_ideal, hit_x, hit_y);
    end
    fprintf('  interval=%d min done\n', intervals(ki));
end

figure;
hold on;
markers = {'-o','-s','-^','-d','-v'};
for ks = 1:numel(slope_errors)
    plot(intervals, 100 * mean_eff_grid(:,ks), markers{ks}, ...
        'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', sprintf('\\sigma_s = %.2f\\circ', slope_errors(ks)));
end
xlabel('Repositioning interval (min)');
ylabel('Mean daily optical efficiency (%)');
title('Sensitivity: tracking interval vs surface slope error');
legend('Location','southwest');
grid on;

% Tracking only (slope=0)
figure;
plot(intervals, 100 * mean_eff_grid(:,1), 'b-o', 'LineWidth', 2, ...
    'MarkerSize', 8);
xlabel('Repositioning interval (min)');
ylabel('Mean daily optical efficiency (%)');
title('Effect of repositioning interval on efficiency (perfect surface)');
grid on;

% Slope only (tracking = 15 min)
figure;
plot(slope_errors, 100 * mean_eff_grid(1,:), 'r-s', 'LineWidth', 2, ...
    'MarkerSize', 8);
xlabel('Surface slope error \sigma_s (deg)');
ylabel('Mean daily optical efficiency (%)');
title('Effect of Mylar surface quality on efficiency (5 min tracking)');
grid on;

%%  Function: run_day
%  Runs one full day simulation and returns mean efficiency.
%  Called by the sensitivity study loop.

function mean_eff = run_day(hours, Gbi, elev_deg, ...
        reposition_interval, slope_error_deg, ...
        f, radius, reflectivity, receiver_radius, ...
        N, sun_spread_deg, normals_ideal, hit_x, hit_y)

    deg_per_interval = 15 * (reposition_interval / 60);
    slope_error_rad  = deg2rad(slope_error_deg);
    spread_rad       = deg2rad(sun_spread_deg);
    dish_area        = pi * radius^2;
    power_per_ray    = dish_area / N;
    receiver_center  = [0; f];
    shadow_fraction  = (2 * receiver_radius) / (2 * radius);
    blocked          = abs(hit_x) < receiver_radius;

    n_steps            = numel(hours);
    captured_power     = zeros(1, n_steps);
    incident_power_vec = zeros(1, n_steps);
    optical_efficiency = zeros(1, n_steps);

    for t_idx = 1:n_steps
        elev    = elev_deg(t_idx);
        Gbi_now = Gbi(t_idx);
        if elev <= 5 || Gbi_now < 10, continue; end

        elev_rad           = deg2rad(elev);
        elev_snapped       = floor(elev / deg_per_interval) * deg_per_interval;
        tracking_error_rad = deg2rad(elev - elev_snapped);

        % Perturb normals
        phi          = slope_error_rad * randn(1, N);
        cos_phi      = cos(phi);
        sin_phi      = sin(phi);
        normals      = zeros(2, N);
        normals(1,:) = cos_phi .* normals_ideal(1,:) ...
                     - sin_phi .* normals_ideal(2,:);
        normals(2,:) = sin_phi .* normals_ideal(1,:) ...
                     + cos_phi .* normals_ideal(2,:);

        DNI_now    = Gbi_now / sin(elev_rad);
        P_incident = DNI_now * dish_area * (1 - shadow_fraction^2);
        incident_power_vec(t_idx) = P_incident;

        theta = spread_rad * randn(1, N) + tracking_error_rad;
        dir   = [sin(theta); -cos(theta)];

        ray_energy          = (DNI_now * power_per_ray) * ones(1, N);
        ray_energy(blocked) = 0;

        captured = 0;
        for i = 1:N
            if ray_energy(i) == 0, continue; end
            d = dir(:,i);
            n = normals(:,i);
            r = d - 2 * dot(d,n) * n;
            ray_energy(i) = ray_energy(i) * reflectivity;

            t_param = linspace(0, 2, 80);
            rx      = hit_x(i) + r(1) * t_param;
            ry      = hit_y(i) + r(2) * t_param;
            dist    = sqrt((rx - receiver_center(1)).^2 + ...
                           (ry - receiver_center(2)).^2);
            if any(dist < receiver_radius)
                captured = captured + ray_energy(i);
            end
        end

        captured_power(t_idx) = captured;
        if P_incident > 0
            optical_efficiency(t_idx) = captured / P_incident;
        end
    end

    daytime  = optical_efficiency > 0;
    if any(daytime)
        mean_eff = mean(optical_efficiency(daytime));
    else
        mean_eff = 0;
    end
end

%%  Function: load_pvgis
%  Parses PVGIS hourly CSV.
%  Columns: time, Gb(i), Gd(i), Gr(i), H_sun, T2m, WS10m, Int

function [hours_out, Gbi_out, elev_out] = load_pvgis(filepath, doy)

    fid = fopen(filepath, 'r');
    if fid < 0
        error('Cannot open: %s', filepath);
    end

    rows = {};
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ischar(line) && numel(line) >= 8 && ...
                all(isstrprop(line(1:8), 'digit'))
            rows{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);

    if isempty(rows)
        error('No data rows found. Check file format.');
    end

    n        = numel(rows);
    doy_all  = zeros(n,1);
    hour_all = zeros(n,1);
    Gbi_all  = zeros(n,1);
    Hsun_all = zeros(n,1);

    for k = 1:n
        p   = strsplit(rows{k}, ',');
        ts  = strtrim(p{1});
        yr  = str2double(ts(1:4));
        mo  = str2double(ts(5:6));
        dy  = str2double(ts(7:8));
        hh  = str2double(ts(10:11));
        mm  = str2double(ts(12:13));

        doy_all(k)  = doy_from_date(yr, mo, dy);
        hour_all(k) = hh + mm/60;
        Gbi_all(k)  = max(0, str2double(strtrim(p{2})));
        Hsun_all(k) = str2double(strtrim(p{5}));
    end

    mask = (doy_all == doy);
    if ~any(mask)
        error('Day %d not found in CSV.', doy);
    end

    hours_out = hour_all(mask);
    Gbi_out   = Gbi_all(mask);
    elev_out  = Hsun_all(mask);

    [hours_out, idx] = sort(hours_out);
    Gbi_out  = Gbi_out(idx);
    elev_out = elev_out(idx);
end

%%  Function: doy_from_date

function doy = doy_from_date(yr, mo, dy)
    dpm = [31 28 31 30 31 30 31 31 30 31 30 31];
    if (mod(yr,4)==0 && mod(yr,100)~=0) || mod(yr,400)==0
        dpm(2) = 29;
    end
    doy = sum(dpm(1:mo-1)) + dy;
end