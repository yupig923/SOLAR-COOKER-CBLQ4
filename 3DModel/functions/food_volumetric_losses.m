function Q_vol = food_volumetric_losses(~, state, V_food, c_mix, m_mix, t_cooking, T_amb, T_sat_use, cooking_time_points)
    t = state.time;
    T = state.u;
    
    is_cooking = false;
    for j = 1:length(cooking_time_points)
        if t >= cooking_time_points(j) && t <= (cooking_time_points(j) + 3600)
            is_cooking = true;
            break;
        end
    end
    
    if is_cooking
        Qco = (c_mix * m_mix / t_cooking) * (T - T_amb);
    else
        Qco = zeros(size(T));
    end
    
    Kv = 0.1;
    Qv = zeros(size(T));
    boiling_idx = (T >= T_sat_use);
    if any(boiling_idx)
        Qv(boiling_idx) = Kv * (T(boiling_idx) - T_sat_use).^2;
    end
    
    Q_vol = -(Qco + Qv) / V_food;
end