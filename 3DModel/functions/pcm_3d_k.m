function k_eff = pcm_3d_k(region, state)

T = state.u;

T_m = 119 + 273.15;
dT  = 10;

k_s = 0.733;
k_l = 0.326;

k_eff = zeros(size(T));

for k = 1:length(T)

    if T(k) < (T_m - dT/2)

        k_eff(k) = k_s;

    elseif T(k) <= (T_m + dT/2)

        frac = (T(k) - (T_m - dT/2)) / dT;

        k_eff(k) = k_s + frac * (k_l - k_s);

    else

        k_eff(k) = k_l;

    end
end

end