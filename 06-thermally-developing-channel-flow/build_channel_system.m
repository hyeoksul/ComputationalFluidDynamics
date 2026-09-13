function model = build_channel_system(cell_size, scheme)
%BUILD_CHANNEL_SYSTEM Assemble a finite-volume convection-diffusion system.

arguments
    cell_size (1,1) double {mustBePositive}
    scheme (1,1) string
end

half_gap = 1.0;
channel_length = 50*half_gap;
channel_height = 2*half_gap;
rho = 1.177;
mu = 1.846e-5;
cp = 1007;
k = 0.02624;
reynolds_gap = 50;
inlet_temperature = 300;
wall_temperature = 400;

nx_intervals = round(channel_length/cell_size);
ny_intervals = round(channel_height/cell_size);
assert(abs(nx_intervals*cell_size-channel_length) < 1e-12);
assert(abs(ny_intervals*cell_size-channel_height) < 1e-12);

x = linspace(0, channel_length, nx_intervals+1);
y = linspace(0, channel_height, ny_intervals+1).';
nx = numel(x);
ny = numel(y);
dx = x(2)-x(1);
dy = y(2)-y(1);

mean_velocity = reynolds_gap*mu/(rho*channel_height);
maximum_velocity = 1.5*mean_velocity;
velocity = maximum_velocity*(1-((y-half_gap)/half_gap).^2);

switch lower(scheme)
    case "central"
        weighting = @(peclet) 1-0.5*abs(peclet);
    case "upwind"
        weighting = @(peclet) ones(size(peclet));
    case "hybrid"
        weighting = @(peclet) max(0, 1-0.5*abs(peclet));
    case "powerlaw"
        weighting = @(peclet) max(0, (1-0.1*abs(peclet))).^5;
    otherwise
        error('Unknown convection scheme: %s', scheme);
end

number_of_unknowns = nx*ny;
matrix = spalloc(number_of_unknowns, number_of_unknowns, ...
    5*number_of_unknowns);
rhs = zeros(number_of_unknowns, 1);
boundary_indices = false(number_of_unknowns, 1);
streamwise_cv_width = dx*ones(1, nx);
streamwise_cv_width([1,end]) = dx/2;

for i = 1:nx
    for j = 1:ny
        row = sub2ind([ny,nx], j, i);
        if j == 1 || j == ny
            matrix(row,row) = 1;
            rhs(row) = wall_temperature;
            boundary_indices(row) = true;
            continue;
        elseif i == 1
            matrix(row,row) = 1;
            rhs(row) = inlet_temperature;
            boundary_indices(row) = true;
            continue;
        end

        east_west_area = dy;
        north_south_area = streamwise_cv_width(i);
        flow = rho*cp*velocity(j)*east_west_area;

        west_diffusion = k*east_west_area/dx;
        west_peclet = flow/west_diffusion;
        a_west = west_diffusion*weighting(west_peclet) + flow;

        if i < nx
            east_diffusion = k*east_west_area/dx;
            east_peclet = flow/east_diffusion;
            a_east = east_diffusion*weighting(east_peclet);
        else
            a_east = 0;
        end

        a_north = k*north_south_area/dy;
        a_south = a_north;
        a_point = a_east+a_west+a_north+a_south;

        matrix(row,row) = a_point;
        matrix(row,sub2ind([ny,nx],j,i-1)) = -a_west;
        if i < nx
            matrix(row,sub2ind([ny,nx],j,i+1)) = -a_east;
        end
        matrix(row,sub2ind([ny,nx],j+1,i)) = -a_north;
        matrix(row,sub2ind([ny,nx],j-1,i)) = -a_south;
    end
end

model = struct();
model.A = matrix;
model.b = rhs;
model.boundary_indices = boundary_indices;
model.x = x;
model.y = y;
model.nx = nx;
model.ny = ny;
model.dx = dx;
model.dy = dy;
model.velocity = velocity;
model.rho = rho;
model.mu = mu;
model.cp = cp;
model.k = k;
model.half_gap = half_gap;
model.hydraulic_diameter = 4*half_gap;
model.reynolds_gap = reynolds_gap;
model.reynolds_hydraulic = rho*mean_velocity*model.hydraulic_diameter/mu;
model.prandtl = cp*mu/k;
model.mean_velocity = mean_velocity;
model.inlet_temperature = inlet_temperature;
model.wall_temperature = wall_temperature;
model.scheme = scheme;
end
