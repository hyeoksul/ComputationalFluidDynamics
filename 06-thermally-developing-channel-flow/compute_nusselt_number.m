function heat_transfer = compute_nusselt_number(model, temperature)
%COMPUTE_NUSSELT_NUMBER Calculate bulk temperature and local heat transfer.

bulk_temperature = zeros(1,model.nx);
local_nusselt = nan(1,model.nx);
correlation_nusselt = nan(1,model.nx);
velocity_integral = trapz(model.y,model.velocity);

for i = 1:model.nx
    bulk_temperature(i) = trapz(model.y, ...
        model.velocity.*temperature(:,i))/velocity_integral;
    if i == 1
        continue;
    end

    lower_gradient = (-3*model.wall_temperature ...
        +4*temperature(2,i)-temperature(3,i))/(2*model.dy);
    upper_gradient = (3*model.wall_temperature ...
        -4*temperature(end-1,i)+temperature(end-2,i))/(2*model.dy);
    heat_flux = 0.5*model.k*(-lower_gradient+upper_gradient);
    coefficient = heat_flux/(model.wall_temperature-bulk_temperature(i));
    local_nusselt(i) = coefficient*model.hydraulic_diameter/model.k;

    x_star = model.x(i)/(model.reynolds_hydraulic*model.prandtl ...
        *model.hydraulic_diameter);
    if x_star < 0.001
        correlation_nusselt(i) = 1.233*x_star^(-1/3)+0.4;
    else
        correlation_nusselt(i) = 7.541 ...
            +6.874*(1000*x_star)^(-0.488)*exp(-245*x_star);
    end
end

developed_index = find(abs(local_nusselt-7.541)/7.541 < 0.01,1,'first');
if isempty(developed_index)
    fully_developed_x = NaN;
else
    fully_developed_x = model.x(developed_index);
end

heat_transfer = struct('bulk_temperature',bulk_temperature, ...
    'local_nusselt',local_nusselt, ...
    'correlation_nusselt',correlation_nusselt, ...
    'fully_developed_x',fully_developed_x);
end
