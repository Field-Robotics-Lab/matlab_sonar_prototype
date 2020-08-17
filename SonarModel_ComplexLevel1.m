%% Single Beam Sonar
%%  Sonar Point-Scatter Model
% Contributors: Andreina Rascon, Derek Olson, Woeng-Sug Choi

clear;
close all;
clc;

%% Input Parameters
% The BlueView P900-45 is used as an example sonar for the purposes of
% demonstrating the point-scatter model


% Sonar Parameters
sonarFreq = 900E3; % Sonar frquency
bandwidth = 29.5e4; % [Hz]
freqResolution = 100e2;
nRays = 3;
beamWidth = 0.1; % radians


% Enviromental Parameters
soundSpeed = 1500.0; % [m/s]

% Target information (which should be obtained from ray gazebo plugin)
% 3 Targets, One dimension, no incident angle
ray.num = nRays;
ray.distance = [2, 5, 10];
ray.alpha = [0, 0, 0];
ray.azimuthAngleWidth = [0.1, 0.1, 0.1];
ray.elevationAngleWidth = [0.1, 0.1, 0.1];

% Surface Reflectivity
mu = 10e-4;


%% Frequency Domain
% 

% Frequency spectrum
fmin = sonarFreq - bandwidth/2*4;
fmax = sonarFreq + bandwidth/2*4;

% Sampling periods
max_T = max(ray.distance)*2/soundSpeed;
delta_f = 1/max_T;
delta_t = 1/(fmax - fmin);

nFreq = round((fmax - fmin) / delta_f);
freq = linspace(fmin,fmax,nFreq);
time = linspace(0,max_T,nFreq);

% Transmit spectrum, frequency domain
S_f = exp(-(freq-sonarFreq).^2.*pi^2./bandwidth^2); 


%% Point Scattering model
% Calculate echo level using the point scatter model

P_f = zeros(ray.num,length(S_f));   % Initialize vector
for k = 1:ray.num    
    xi_z = rand;   % generate a random number, (Gaussian noise)
    xi_y = rand;   % generate another random number, (Gaussian noise)
    
    % Assume every ray hit targets
    Nscatters = ray.num;  % This should be the number of surfaces a ray intersects with. For now, using constant equal to ray number for simplicity
    
    for j = 1:Nscatters
        % scattering amplitude, includes overall Gaussian noise and target strength
        % alpha needs to be the angle between incident ray and normal vector
        alpha = ray.alpha(j);
        distance = ray.distance(j);
        azimuthAngleWidth = ray.azimuthAngleWidth(j);
        elevationAngleWidth = ray.elevationAngleWidth(j);
        amplitude(j) = ((xi_z+1j*xi_y)/sqrt(2))*(sqrt(mu*cos(alpha)^2*distance^2*azimuthAngleWidth*elevationAngleWidth));
        
        % travel time
        travelTime = distance*2/soundSpeed;
        
        % define wave vector
        k = 2*pi*freq/soundSpeed;
        
        % Summation of Echo returned from a signal
        P_f(j,:) = P_f(j,:) + S_f.*amplitude(j).*exp(1i.*k.*distance*2)./(distance^2);
    end
end

%% Plotting

% Plot inverse fast fourier transform (frequency domain to time domain)
p_t = ifft(P_f,[],2);

figure (1)
plot(time,p_t)

% Plot Sound Pressure Level of Echo Level
SPL = 20*log(abs(p_t)); % sound pressure level, [dB]

figure (2)
plot(time,SPL)

%% Beam Forming

% Beam Design Parameters
NB = nRays; % dimensionless number of beams
NT = nFreq; % dimensionless number of time samples

% th = absolute look angle
th = beamWidth.*(-NB/2+1:NB/2); % assuming the sonar is centered at zero theta

% th0 = look angle of one of the beams (constant)
th0 = 0.02;

sinc = @(x) sin(x)/x;
B_th = sinc((th - th0)/beamWidth); % beam pressure response

p_t_modified = zeros(nRays,NT);

for ii = 1:NB
        thetaCurrent = th(ii);
        beamPatternCurrent = sinc((th - thetaCurrent)./beamWidth);
        p_t_modified(ii,:) = p_t(ii,:).*beamPatternCurrent; %p_t is NBxNT, and beamPatternCurrent is NBx1
end

%% Plotting Beam Formed Data

% Plot Sound Pressure Level of Echo Level
SPL_modified = 20*log(abs(p_t_modified)); % sound pressure level, [dB]


figure (3)
plot(time,p_t_modified)

figure (4)
plot(time,SPL_modified)




