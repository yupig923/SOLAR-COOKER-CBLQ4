%New 3D Heat Transfer Model (Partial differetial equations toolbox needed)


%Clear and clc
clear;
clc;
global solar_time captured_power pot_eps;

addpath("functions");
savepath;
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
%simulation_time_outdoor_lunch = [t_lunch t_lunch+3600];
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
V_food = pi * (0.15)^2 * 0.08; % Volume of core food chamber (m³)


pot_absorptivity = 0.95; % rough/dull aluminum

c_water = 4148;
c_rice = 1700;
c_pot = 500;

c_mix = (m_water * c_water + m_rice * c_rice) / (m_water + m_rice);

%

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




%Initialize the 3D model
thermalModel_indoors = createpde('thermal', 'transient');
thermalModel_outdoors = createpde('thermal', 'transient');
%Multicylinder 3D model (Will be replaced by actual pot model)
g_food = multicylinder(r_in_pot, h_in_pot);
g_pcm  = multicylinder([r_in_pot, r_out_pot], h_out_pot);
g_ins  = multicylinder([r_out_pot, r_insulation], h_insulation);

gm_indoors=multicylinder([r_in_pot, r_out_pot, r_insulation], h_insulation);
thermalModel_indoors.Geometry = gm_indoors;

gm_outdoors=g_pcm;
thermalModel_outdoors.Geometry = gm_outdoors;

generateMesh(thermalModel_outdoors, 'Hmax', 0.08);

%For new pot model (Guus will make it)
%gm = importGeometry(thermalModel, 'cooking_pot.stl', 'VectorLengthScale', 0.001);


%Check the unit of the dimensions if incorrect change vector lengthscale



% Plot the 3D geometry
figure;
pdegplot(thermalModel_outdoors, 'FaceAlpha', 0.5, 'CellLabels', 'on');
title('Verified 3D Concentric Pot Geometry Without Insulation');

figure;
pdegplot(thermalModel_indoors, 'FaceAlpha', 0.5, 'CellLabels', 'on');
title('Verified 3D Concentric Pot Geometry With Insulation');

%Layers replaced by 3D domains


%Assign thermal properties


%For outdoors model

% Food Subdomain (Cell 1)
thermalProperties(thermalModel_outdoors, 'Cell', 1, 'ThermalConductivity', k_mix, ...
                                          'MassDensity', 1000, ...
                                          'SpecificHeat', c_mix);

% PCM Subdomain (Cell 2) 
thermalProperties(thermalModel_outdoors, 'Cell', 2, 'ThermalConductivity', @pcm_3d_k, ...
                                          'MassDensity', 1480, ...
                                          'SpecificHeat', @pcm_3d_cp);



%For indoors model
% Food Subdomain (Cell 1)
thermalProperties(thermalModel_indoors, 'Cell', 1, 'ThermalConductivity', k_mix, ...
                                          'MassDensity', 1000, ...
                                          'SpecificHeat', c_mix);

% PCM Subdomain (Cell 2) 
thermalProperties(thermalModel_indoors, 'Cell', 2, 'ThermalConductivity', @pcm_3d_k, ...
                                           'MassDensity', 1480, ...
                                           'SpecificHeat', @pcm_3d_cp);

% Insulation Subdomain (Cell 3)
thermalProperties(thermalModel_indoors, 'Cell', 3, 'ThermalConductivity', k_ins, ...
                                          'MassDensity', 100, ... 
                                          'SpecificHeat', c_ins);




%% =========================
% CORRECTED PCM ENERGY STORED (outdoor phase)
%No cooking From 0700-1300 1400-1750
% ==========================
%
T_init = T_amb;
thermalIC(thermalModel_outdoors, T_init);


figure;
pdegplot(thermalModel_outdoors, 'FaceLabels', 'on', 'FaceAlpha', 0.4);
title('Boundary Face Labels');
view(30,20);

outerFaces = [4];
innerFaces=[2];
topFace    = [3];
bottomFace = [1];
% The faces may be incoreect double check
allExternalFaces = unique([outerFaces topFace bottomFace]);
allSolarFaces=unique([outerFaces topFace]);

A_abs = sum(arrayfun(@(f) faceArea(gm_outdoors,f), allSolarFaces));
%BoundaryConditions For Outdoor Charging Phase
%if not cooking 
thermalBC(thermalModel_outdoors, 'Face', allSolarFaces, ...
    'HeatFlux', @(region,state)solar_flux_3d(region,state,A_abs), ...
    'ConvectionCoefficient', h_free_outdoor, ...
    'AmbientTemperature', T_amb);


cooking_time_points = [t_lunch, t_dinner];

internalHeatSource(thermalModel_outdoors, @(region, state) food_volumetric_losses(region, state, ...
    V_food, c_mix, m_mix, t_cooking, T_amb, T_sat_use, cooking_time_points), 'Cell', 1);

%if cooking 
% thermalBC(thermalModel_outdoors, 'Face', allExternalFaces, ...
%     'HeatFlux', @solar_flux_3d, ...
%     'ConvectionCoefficient', h_free_outdoor, ...
%     'AmbientTemperature', T_amb);


%% =========================
% Outdoor Simulation
%% =========================

fprintf('Running outdoor simulation...\n');

n_steps_outdoor = 400;
tlist_outdoor = linspace(simulation_time_outdoor(1), ...
                         simulation_time_outdoor(2), ...
                         n_steps_outdoor);

result_outdoor = solve(thermalModel_outdoors, tlist_outdoor);
fprintf('Outdoor simulation finished...\n');

%% =========================
% Extract Temperatures
%% =========================

T_outdoor1 = result_outdoor.Temperature;
time_outdoor1 = result_outdoor.SolutionTimes;

% Average temperatures by cell
food_nodes = findNodes(result_outdoor.Mesh,'region','Cell',1);
pcm_nodes  = findNodes(result_outdoor.Mesh,'region','Cell',2);

T_food_avg = zeros(size(time_outdoor1));
T_pcm_avg  = zeros(size(time_outdoor1));


for i = 1:length(time_outdoor1)
    T_food_avg(i) = mean(T_outdoor1(food_nodes,i));
    T_pcm_avg(i)  = mean(T_outdoor1(pcm_nodes,i));

end


%% =========================
% PCM Stored Energy Calculation
%% =========================

Q_pcm = zeros(size(T_pcm_avg));

for i = 1:length(T_pcm_avg)

    Tpcm = T_pcm_avg(i);

    if Tpcm < (T_m_pcm - delta_Tm_pcm/2)

        % Fully solid
        Q_pcm(i) = m_pcm * c_pcm_solid * (Tpcm - T_init);

    elseif Tpcm <= (T_m_pcm + delta_Tm_pcm/2)

        % Mushy region
        sensible_solid = m_pcm * c_pcm_solid * ...
            ((T_m_pcm - delta_Tm_pcm/2) - T_init);

        melt_fraction = (Tpcm - (T_m_pcm - delta_Tm_pcm/2)) ...
            / delta_Tm_pcm;

        latent = melt_fraction * m_pcm * L_pcm;

        mushy_cp = m_pcm * c_pcm_liquid * ...
            (Tpcm - (T_m_pcm - delta_Tm_pcm/2));

        Q_pcm(i) = sensible_solid + latent + mushy_cp;

    else

        % Fully liquid
        sensible_solid = m_pcm * c_pcm_solid * ...
            ((T_m_pcm - delta_Tm_pcm/2) - T_init);

        latent = m_pcm * L_pcm;

        sensible_liquid = m_pcm * c_pcm_liquid * ...
            (Tpcm - (T_m_pcm + delta_Tm_pcm/2));

        Q_pcm(i) = sensible_solid + latent + sensible_liquid;

    end
end
pdeplot3D(thermalModel_outdoors,'ColorMapData',result_outdoor.Temperature(:,end))

% Set time range for indoor simulation
indoorStartHour = 17 + 50/60; % 17:50 in hours
indoorCookingStartHour=18;
indorCookingEndHour=19;
indoorEndHour = 22;            % 22:00 in hours
indoorTime = (hours - indoorStartHour) * 3600; % Convert to seconds


%% =========================
% Indoor simulation (17:50 - 22:00)
% ==========================


thermalIC(thermalModel_indoors, T_init)



