%% Sonar Point-Scatter Model
%% 

clear;clc;close all;

bandwidth = 29.5e4; % [Hz]
freqResolution = 100e2;
nRays = 3;

%-------------------------------------
% Given Input parameters
%-------------------------------------
soundSpeed = 1500.0; % [m/s]

% Sonar properties
% BlueView P900-90 Modified
sonarFreq = 900E3; % Sonar frquency

% Target information (which should be obtained from ray gazebo plugin)
% 3 Targets, One dimension, no incident angle
ray.num = nRays;
ray.distance = [2, 5, 10];
ray.alpha = [0, 0, 0];
ray.azimuthAngleWidth = [0.1, 0.1, 0.1];
ray.elevationAngleWidth = [0.1, 0.1, 0.1];



% Surface Reflectivity
mu = 10e-4;

%-------------------------------------
% Computational Variables
%-------------------------------------
% Transmit Spectrum
fmin = sonarFreq - bandwidth/2*4;
fmax = sonarFreq + bandwidth/2*4;
% freq = linspace(fmin,fmax,freqResolution); % Frequency vector

max_T = max(ray.distance)*2/soundSpeed;
delta_f = 1/max_T;
delta_t = 1/(fmax - fmin);

nFreq = round((fmax - fmin) / delta_f);
freq = linspace(fmin,fmax,nFreq);
time = linspace(0,max_T,nFreq);

S_f = exp(-(freq-sonarFreq).^2.*pi^2./bandwidth^2); % Transmit spectrum, frequency domain


%-------------------------------------
%------  Point Scattering model ------
%-------------------------------------
% Calculate echo level with TS considered 
%-------------------------------------
P_f = zeros(ray.num,length(S_f));
for k = 1:ray.num    
    xi_z = rand;   % generate a random number, part of Gaussian noise
    xi_y = rand;   % generate another random number, part of Gaussian noise
    
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

% figure;plot(freq,real(P_f));
% figure;plot(freq,phase(P_f));

p_t = ifft(P_f,[],2);

figure (1)
plot(time,p_t)

SPL = 20*log(abs(p_t));

figure (2)
plot(time,SPL)

