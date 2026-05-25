function Qv_vol = evaporation_loss_3d(~, state)
    T_food = state.u;
    T_sat_use = (100 - 2) + 273.15; % 371.15 K
    Kv = 0.1;
    s = 2;
    
    V_food = pi * (0.15)^2 * 0.08;
    Qv_vol = zeros(size(T_food));
    
    % Evaporation kicks in if local node elements pass the threshold
    boiling_idx = (T_food >= T_sat_use);
    Q_watts = Kv * (T_food(boiling_idx) - T_sat_use).^s;
    
    % Return as negative volumetric heat sink
    Qv_vol(boiling_idx) = -Q_watts / V_food;
end