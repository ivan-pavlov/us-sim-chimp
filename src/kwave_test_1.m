clearvars;

% simulation settings
DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 20;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]

% set total number of grid points not including the PML
Nx = 256 - 2 * pml_x_size;      % [grid points]
Ny = 128 - 2 * pml_y_size;      % [grid points]
Nz = 128 - 2 * pml_z_size;      % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                      % [m]

% calculate the spacing between the grid points
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 1.5e6;        % [Hz]
tone_burst_cycles = 4;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 32;  	% total number of transducer elements
transducer.element_width = 2;       % width of each element [grid points]
transducer.element_length = 24;  	% length of each element [grid points]
transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = 20e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Hanning';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
number_active_elements = 32;
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

% define a large image size to move across
number_scan_lines = 96;
Nx_tot = Nx;
Ny_tot = Ny + number_scan_lines * transducer.element_width;
Nz_tot = Nz;

% define a random distribution of scatterers for the medium
background_map_mean = 1;
background_map_std = 0.008;
background_map = background_map_mean + background_map_std * randn([Nx_tot, Ny_tot, Nz_tot]);

% define a random distribution of scatterers for the highly scattering
% region
scattering_map = randn([Nx_tot, Ny_tot, Nz_tot]);
scattering_c0 = c0 + 25 + 75 * scattering_map;
scattering_c0(scattering_c0 > 1600) = 1600;
scattering_c0(scattering_c0 < 1400) = 1400;
scattering_rho0 = scattering_c0 / 1.5;

% define properties
sound_speed_map = c0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;

% define a sphere for a highly scattering region
radius = 6e-3;      % [m]
x_pos = 27.5e-3;    % [m]
y_pos = 20.5e-3;    % [m]
scattering_region1 = makeBall(Nx_tot, Ny_tot, Nz_tot, round(x_pos/dx), round(y_pos/dx), Nz_tot/2, round(radius/dx));

% assign region
sound_speed_map(scattering_region1 == 1) = scattering_c0(scattering_region1 == 1);
density_map(scattering_region1 == 1) = scattering_rho0(scattering_region1 == 1);

% define a sphere for a highly scattering region
radius = 5e-3;      % [m]
x_pos = 30.5e-3;    % [m]
y_pos = 37e-3;      % [m]
scattering_region2 = makeBall(Nx_tot, Ny_tot, Nz_tot, round(x_pos/dx), round(y_pos/dx), Nz_tot/2, round(radius/dx));

% assign region
sound_speed_map(scattering_region2 == 1) = scattering_c0(scattering_region2 == 1);
density_map(scattering_region2 == 1) = scattering_rho0(scattering_region2 == 1);

% define a sphere for a highly scattering region
radius = 4.5e-3;    % [m]
x_pos = 15.5e-3;    % [m]
y_pos = 30.5e-3;    % [m]
scattering_region3 = makeBall(Nx_tot, Ny_tot, Nz_tot, round(x_pos/dx), round(y_pos/dx), Nz_tot/2, round(radius/dx));

% assign region
sound_speed_map(scattering_region3 == 1) = scattering_c0(scattering_region3 == 1);
density_map(scattering_region3 == 1) = scattering_rho0(scattering_region3 == 1);

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% preallocate the storage
scan_lines = zeros(number_scan_lines, kgrid.Nt);

% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

number_scan_lines = 1;

% run the simulation if set to true, otherwise, load previous results from
% disk
if RUN_SIMULATION

    % set medium position
    medium_position = 1;
    
    % loop through the scan lines
    for scan_line_index = 1:number_scan_lines
        
        % update the command line status
        disp('');
        disp(['Computing scan line ' num2str(scan_line_index) ' of ' num2str(number_scan_lines)]);

        % load the current section of the medium
        medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny - 1, :);
        medium.density = density_map(:, medium_position:medium_position + Ny - 1, :);
        
        % run the simulation
        sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});

        % extract the scan line from the sensor data
        scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);
        
        % update medium position
        medium_position = medium_position + transducer.element_width;

    end

    % save the scan lines to disk
    save example_us_bmode_scan_lines scan_lines;
    
else
    
    % load the scan lines from disk
    load example_us_bmode_scan_lines;
    
end



