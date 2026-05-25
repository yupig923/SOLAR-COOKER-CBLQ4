


function q = solar_flux_3d(region, state,A_abs)

global solar_time captured_power;

if isempty(state.time)
    fprintf('empty time')
    q = zeros(size(region.x));
    return;
end
Psolar = interp1(solar_time, captured_power, state.time, 'linear', 0);


%Area calculations (Needs to be replaced by actual area from region)
%r_abs=0.21;

qsolar = Psolar / A_abs;

% PDE Toolbox convention: positive = into domain
q = qsolar * ones(size(region.x));



end