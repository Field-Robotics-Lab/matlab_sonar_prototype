%% Single Beam Sonar
%%  Sonar Point-Scatter Model
% Contributors: Andreina Rascon, Derek Olson, Woeng-Sug Choi

clear;
close all;
clc;

% Input Parameters
% The BlueView P900-45 is used as an example sonar for the purposes of
% demonstrating the point-scatter model


%% Sonar Parameters
sonarFreq = 900E3; % Sonar frquency
bandwidth = 29.5e4; % [Hz]
freqResolution = 100e2;
beamWidth = 0.1; % radians
nRays = 3;
nBeams = 1;

% Beam Parameters
NB = nRays; % dimensionless number of beams

% th = absolute look angle
th = beamWidth.*(-NB/2+1:NB/2); % assuming the sonar is centered at zero theta
% th0 = look angle of one of the beams (constant)
th0 = 0.0; % only one beam right now

sinc = @(x) sin(x)/x;
B_th = sinc((th - th0)/beamWidth); % beam pressure response

%% Input Parameters
soundSpeed = 1500.0; % [m/s]

% % Incidence angle calculation
% ray.vector = [z_xi z_yi z_zi];  % directions, magnitude is 1
% norm.vector = [n_xi n_yi n_zi]; % directions, magnitude is 1
% psi = acos(dot(ray.vector,norm.vector));    % [rad]
% ray.alpha = pi-psi;

% Target information (which should be obtained from ray gazebo plugin)
% 3 Targets, One dimension, no incident angle

ray.distance = [2, 5, 10];
ray.alpha = [0, 0, 0];
ray.azimuthAngleWidth = [0.1, 0.1, 0.1];
ray.elevationAngleWidth = [0.1, 0.1, 0.1];

% Surface Reflectivity
mu = 10e-4;


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
NT = nFreq; % dimensionless number of time samples

% Transmit spectrum, frequency domain
S_f = exp(-(freq-sonarFreq).^2.*pi^2./bandwidth^2);

% define wave vector
kw = 2*pi*freq/soundSpeed;
        
% add attenuation
absorption = 0.0354; % [dB/m]
attenuation = absorption*log(10)/20;
K = kw + 1j*attenuation;

%% Point Scattering model
% Calculate echo level using the point scatter model

P_f = zeros(nBeams,nRays,nFreq);    % pre-allocate P(f)
P_f_ray = zeros(nRays,nFreq);       % pre-allocate ray P(f)

% Beamforming Loop

for k = 1:nBeams
    
    % Ray Loop
    for i = 1:nRays
        xi_z = rand;   % generate a random number, (Gaussian noise)
        xi_y = rand;   % generate another random number, (Gaussian noise)
        
        alpha = ray.alpha(i);
        distance = ray.distance(i);
        azimuthAngleWidth = ray.azimuthAngleWidth(i);
        elevationAngleWidth = ray.elevationAngleWidth(i);
        amplitude(i) = ((xi_z+1j*xi_y)/sqrt(2))*(sqrt(mu*cos(alpha)^2*distance^2*azimuthAngleWidth*elevationAngleWidth));
        
        % travel time
        travelTime = distance*2/soundSpeed;
        
        
        % Summation of Echo returned from a signal
        P_f_ray(i,:) = (P_f_ray(i,:) + S_f.*amplitude(i).*exp(1i.*K.*distance*2)./(distance^2));
        
    end
        P_f(k,i,:) = P_f_ray(i,:); 

%     thetaCurrent = th(k);
%     beamPatternCurrent = sinc((th - thetaCurrent)./beamWidth);
%     p_t_modified(k,,i,:) = p_t(k,,i,:).*beamPatternCurrent; %p_t is NBxNT, and beamPatternCurrent is NBx1
    
end

P_beam = squeeze(sum(P_f,3));

P_t = ifft(P_beam,[],2);

% summation = pulsint(P_f,'noncoherent');
%% Plotting

% Plot inverse fast fourier transform (frequency domain to time domain)
% p_t = ifft(P_f,[],2);

% figure (1)
% plot(time,P_t)

% % Plot Sound Pressure Level of Echo Level
% SPL = 20*log(abs(p_t)); % sound pressure level, [dB]
% 
% 
% % Plotting Beam Formed Data
% 
% % Plot Sound Pressure Level of Echo Level
% SPL_modified = 20*log(abs(p_t_modified)); % sound pressure level, [dB]
% 
% 
% figure (3)
% plot(time,p_t_modified)
% 
% figure (4)
% plot(time,SPL_modified)




