function cp = pcm_3d_cp(region, state)

T = state.u;

T_m = 119 + 273.15;
L   = 340e3;
dT  = 10;

cp_s = 1380;
cp_l = 2760;

cp = zeros(size(T));

for k = 1:length(T)

    if T(k) < (T_m - dT/2)

        cp(k) = cp_s;

    elseif T(k) <= (T_m + dT/2)

        latent_cp = L / dT;

        frac = (T(k) - (T_m - dT/2)) / dT;

        cp_base = cp_s + frac * (cp_l - cp_s);

        cp(k) = cp_base + latent_cp;

    else

        cp(k) = cp_l;

    end
end

end