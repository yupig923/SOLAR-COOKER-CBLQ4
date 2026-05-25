function Qco_vol = cooking_loss_3d(~, state)
    T_food = state.u;
    T_amb = 20 + 273.15;
    
    % Parameters from your file
    m_mix = 0.075 + (0.075 * 2); % 0.225 kg
    c_mix = 3332.1;              % Calculated blend Cp
    t_cooking = 3600;            % 1 hour duration
    
    % Volume of food domain (r = 0.15m, h = 0.08m)
    V_food = pi * (0.15)^2 * 0.08; 
    
    % Total dynamic Watt loss per unit volume (W/m^3)
    Q_total_watts = (c_mix * m_mix / t_cooking) * (T_food - T_amb);
    
    % Return as a negative internal source (sink)
    Qco_vol = -Q_total_watts / V_food;
end