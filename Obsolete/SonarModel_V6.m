%% Single Beam Sonar
%% Sonar Point-Scatter Model
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


% Beam Parameters
nBeams = 7;
th = linspace(-beamWidth,beamWidth,nBeams);

%% Input Parameters
soundSpeed = 1500.0; % [m/s]

% % Incidence angle calculation
% ray.vector = [z_xi z_yi z_zi];  % directions, magnitude is 1
% norm.vector = [n_xi n_yi n_zi]; % directions, magnitude is 1
% psi = acos(dot(ray.vector,norm.vector));    % [rad]
% ray.alpha = pi-psi;

% Target information (which should be obtained from ray gazebo plugin)
% 3 Targets, One dimension, no incident angle

ray.distance = [1, 5, 1];
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
S_f = 1e11*exp(-(freq-sonarFreq).^2.*pi^2./bandwidth^2);

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
P_t_modified = zeros(nBeams,nFreq);
% beamPatternCurrent = zeros(nBeams);


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
        amplitude = ((xi_z+1j*xi_y)/sqrt(2))*(sqrt(mu*cos(alpha)^2*distance^2*azimuthAngleWidth*elevationAngleWidth));
        
        % travel time
        travelTime = distance*2/soundSpeed;
        
        
        % Summation of Echo returned from a signal
        P_f_ray(i,:) = (P_f_ray(i,:) + S_f.*amplitude.*exp(1i.*K.*distance*2)./(distance^2));
        
    end
        P_f(k,i,:) = P_f_ray(i,:); 
     
        P_beam = squeeze(sum(P_f,2));
        
        beamPatternCurrent = (sinc((pi*0.88./(pi*.88./th))./beamWidth)).^2;
        
        % frequency domain to time domain
        P_t = ifft(P_beam,[],2);
     
        P_t_modified = beamPatternCurrent'.*P_t; %p_t is nFreqxnBeams, and beamPatternCurrent is nBeamsx1
    
end



%% Plotting

% Plot inverse fast fourier transform 


% % figure (1)
plot(time,P_t)
xlabel('Time, [s]')
ylabel('Pressure, [Pa]')

% % % Plot Sound Pressure Level of Echo Level
SPL = 20*log(abs(P_t)); % sound pressure level, [dB]

figure (2)
plot(time,SPL)
xlabel('Time, [s]')
ylabel('Sound Pressure Level, [dB]')


% % Plotting Beam Formed Data
 
SPL_modified = 20*log(abs(P_t_modified)); % sound pressure level, [dB]

figure (3)
plot(time,SPL_modified)
xlabel('Time, [s]')
ylabel('Sound Pressure Level, [dB]')

% Plot Image
% use imagesc, zoom in

% imagesc(SPL_modified')

function y = sinc(x)
     y = abs(sin(x)./x);
     y(x==0) = 1;
end

