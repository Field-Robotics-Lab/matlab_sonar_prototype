% point scatterer model for sphere
f0 = 900e3; % Center Frequency: Hz
c = 1500; % sound velocity: m/s
lambda = c/f0; % wavelength
% attenuation - right now it's the imaginary part of the wavnenumber, but
% we can use the Francois-Garrison model to compute this using S,T,P,pH
% TODO: add in physical model of attenuation
kpp = 0.00001;
bw = 29.9e3; % bandwidth
tau = 1/bw; % pulse length,
tauX = tau*c/2; %pulse resolution, meters;

sourceLevel = 220; % db re 1 muPa;
pref = 1e-6; %1 micro pascal (muPa)
s0 = sqrt(10.^(220/10))*pref; % source term

% define signal parameters

maxRange = 5;
maxT = maxRange*2/c;
delta_f = 1/maxT;
delta_t = 1/bw;
Nfreq = ceil(bw/delta_f);
delta_f = bw/Nfreq;
Ntime = Nfreq; % #becausefft
t_vector = (0:Ntime-1).*delta_t;
range_vector = t_vector.*c./2;
windowFreq = hamming(Nfreq); % this can be user specified.
% normalize the window so it has unit (incoherent) gain.
windowFreq = windowFreq./sqrt(sum(windowFreq.^2));
windowFreq = windowFreq(:).';% make a row vector so we can use broadcasting.
if mod(Nfreq,2) == 0
    f_vec = delta_f.*(-Nfreq/2+1 : Nfreq/2);
else
    f_vec = delta_f.*(-(Nfreq-1)/2+1 : (Nfreq-1)/2);
end

k_vec = f_vec.*2.*pi./c;

k_vec_complex = k_vec + 1i*kpp;



% define sonar geometry.
% position is defined in an E-N-U coordinate system.
% local UTM for E-N, and z coordinate is referenced to the sea surface
sonarPos = [0, 0, 0];
targetPos = [2.1,0 0];
% orientation is defined as heading-pitch-roll in radians
% could also define it in spherical coordinates
% heading is angle from true north - rotation around the z axis
% pitch is angle from horizontal, rotation around the y-axis.
% roll is angle from vertical - rotation around the x axis
sonarOrientation = [20 , 0 , 0]*pi/180;
% number of beams. 
Nbeams = 128;
% width of each beam in radians
sonarBeamwidth = 2*pi/180;
if mod(Nbeams,2) == 0
    sonarBeams =  (-Nbeams/2+1:Nbeams/2)*sonarBeamwidth;
else
    sonarBeams=  (-(Nbeams-1)/2+1:(Nbeams-1)/2)*sonarBeamwidth;
end
sonarBeams = sonarBeams(:);
R_heading = [cos(sonarOrientation(1)) -sin(sonarOrientation(1)) 0;
             sin(sonarOrientation(1)) cos(sonarOrientation(1))  0;
             0                        0                         1];
         
R_pitch = [cos(sonarOrientation(2)) 0 sin(sonarOrientation(2));
           0                        1 0                       ;
           -sin(sonarOrientation(2)) 0 cos(sonarOrientation(2))];
       
      
R_roll = [1 0                        0                        ;
          0 cos(sonarOrientation(3)) -sin(sonarOrientation(3));
          0 sin(sonarOrientation(3)) cos(sonarOrientation(3))];
R_total = R_heading * R_pitch * R_roll;
% the beams are oreinted so that the center beam (0 degrees) points north.
% We need to rotate the beams
beamNormals = [cos(sonarBeams) sin(sonarBeams) zeros(size(sonarBeams))];
beamNormalsRotated = beamNormals*R_total;
% we also need to see where (1 0 0), (0 1 0) vector and (0 0 1) vectors get rotated
% to, since those define the horizontal and vertical beam patterns
xHatPrime = [1 0 0]*R_total;
yHatPrime = [0 1 0]*R_total;
zHatPrime = [0 0 1]*R_total;
% for the beampattern calcuations, we also need a way to compute 
% make sure the sonar is rotating the correct way
plotBeams = 0;
if plotBeams
figure(1)
clf
quiver3(zeros(size(beamNormals(:,1))),zeros(size(beamNormals(:,2))), zeros(size(beamNormals(:,3))),beamNormals(:,1),beamNormals(:,2),beamNormals(:,3))
hold on
quiver3(zeros(size(beamNormals(:,1))),zeros(size(beamNormals(:,2))), zeros(size(beamNormals(:,3))),beamNormalsRotated(:,1),beamNormalsRotated(:,2),beamNormalsRotated(:,3),'--')
axis equal
axis tight
view(3)
end
% define a beam pattern function. It's a sinc right ntaow, but it could
% theoretically be anything.
% matlab sinc is sin(pi*x)/(pi*x). have to divide by pi on the input to

% sinc - this is done last to differntiate it from the other factors of pi
beamPatternHorizontal = @(phi) sinc((2*pi*0.88/sonarBeamwidth)./2.*sin(phi)./pi);


%% uncomment to create a sphere using the distmesh software
% (http://persson.berkeley.edu/distmesh/)
% Need to run  "mex trisurfup.cpp" in the directory to use distmeshsurface.m
% 
sphereRadius =tauX*20;
sphere_mesh_size = tauX/3;
% 
% fd=@(p) dsphere(p,0,0,0,sphereRadius);
% [spherePts,sphereTriangles]=distmeshsurface(fd,@huniform,sphere_mesh_size,sphereRadius*1.1*[-1,-1,-1;1,1,1]);
% 
% save(fullfile(dataDir,'sphereModel'),'spherePts','sphereTriangles')


%% set up the geometry

% get surface normals.
% load(fullfile(sphereModel,'sphereModel'))
load('sphereModel.mat')
targetFaces = sphereTriangles;
targetPoints = spherePts;
% translate target to target position;
targetPoints = targetPoints + targetPos;
targetMu = 1e-3;
%



FV.faces = targetFaces;
FV.vertices = targetPoints;
centroids = 1/3 * (targetPoints(targetFaces(:,1),:) + ...
             targetPoints(targetFaces(:,2),:) + ...
             targetPoints(targetFaces(:,3),:));
% normals is used for the scattering calculation from the target.
[targetNormals,targetAreas] = patchnormals(FV);
nTargetPatches = size(targetFaces,1);
%% verify that normals are outward-faceing
figure(1)
clf
% plot mesh
patch(FV,'FaceNormals',targetNormals,'FaceColor','Red')
hold on
% plot normals on the mesh
quiver3(centroids(:,1),centroids(:,2),centroids(:,3),targetNormals(:,1),targetNormals(:,2),targetNormals(:,3))

%% geometry that is used in all point scattering models.
% in a GPU or C++ implementaiton, these should be calcuated within the for
% loop over beams and scatterers. Anything that depends on the target point
% scatterer position or orientation should be calcuated on the fly on the
% GPU.

% define differences in position between source and target
deltaX = centroids(:,1) - sonarPos(:,1);
deltaY = centroids(:,2) - sonarPos(:,2);
deltaZ = centroids(:,3) - sonarPos(:,3);
% distance to each face
R = sqrt(deltaX.^2 + deltaY.^2 + deltaZ.^2);
Rhat = [deltaX deltaY deltaZ]./R;
% compute dot product between normal and rRhat


% target coordinates in the local sonar xyz frame
xPrime = xHatPrime(1,1)*centroids(:,1) + xHatPrime(1,2)*centroids(:,2) + xHatPrime(1,3)*centroids(:,3,:);
yPrime = yHatPrime(1,1)*centroids(:,1) + yHatPrime(1,2)*centroids(:,2) + yHatPrime(1,3)*centroids(:,3,:);
zPrime = zHatPrime(1,1)*centroids(:,1) + zHatPrime(1,2)*centroids(:,2) + zHatPrime(1,3)*centroids(:,3,:);

thetaPrime = atan(zPrime ./ xPrime);
phiPrime = atan(yPrime./xPrime);
% thetaPrime = atan2(xPrime, zPrime);
% phiPrime = atan2(xPrime , yPrime);
rPrime = [xPrime yPrime zPrime];
randAmps = (randn([nTargetPatches,1]) + 1i*randn([nTargetPatches,1]))./sqrt(2);
%% Point scattering model 1: Most accurate
% All points on the target contribute to all the beams.
% THis is the most accurate way of doing things (for the point scattering
% model

p_1 = zeros(Nbeams , Nfreq);


% This loop uses vectorization / broadcasting to perform all calculations
% per scatterer at once.
% parallel for loop for each beam.
% for a GPU, we can use thread and blocks and grids to organize the data.
tic
for ii = 1 : Nbeams
    currentPhi = sonarBeams(ii);
    currentTheta = 0;
    currentNormal = beamNormalsRotated(ii,:);
    cosAngle = -(currentNormal(1)*targetNormals(:,1) ...
        + currentNormal(2)*targetNormals(:,2) ...
        + currentNormal(3)*targetNormals(:,3));
    
    lambertTerm = targetMu.*cosAngle.^2;
    % if the angle is negative, then it's zero - this is fake shadowing.
    % This is already done by occlusion checking, I thik
    lambertTerm = step(cosAngle).*lambertTerm; 
%     lambertTerm = ones(size(lambertTerm));
    propagationTerm = 1./R.^2.*exp(-2.*kpp.*R); %two-way loss of pressure (not intensity). 
    phiBeam = phiPrime - currentPhi;
    thetaBeam = thetaPrime - currentTheta;
    directivityTerm = beamPatternHorizontal(phiBeam);
    
    amplitude = randAmps ...
        *s0.*propagationTerm.*directivityTerm.* ...
        sqrt(lambertTerm.*targetAreas); % sqrt is because the lambert model is for intensity, and we are working with pressure.

    pFreq = sum(exp(1i.*2*R.*k_vec).*amplitude.*windowFreq,1); %add up all terms and time delays; can take a ton of memory...
    % put this in the order expected by matlab ifft
%     pFreq = [pFreq(Nfreq/2:end) pFreq(1:Nfreq/2-1)];
    pTime = fft(pFreq,[],2).*delta_f;
%     pTime = fftshift(pTime,2);
    p_1(ii,:) = pTime;
    disp([ii Nbeams])
end
% tt = exp(1i.*2*R.*k_vec).*amplitude.*windowFreq;
% ttt = fft(tt,[],2);
% [maxval,maxind] = max(abs(ttt),[] , 2);
ctime1 = toc;
fprintf('Full Method: %g s.\n',ctime1);

%% Point scattering model 2: Less accurate
% beams only see within the beamwidth
% point spread function is added in after the fact by adding beams and
% weighting using the beam pattern. This is theoretically identical to the
% first method, but we don't need to calculate so many point scatterers per
% beam.
% This is the method Andi Rascon used in 
p_2 = zeros(Nbeams , Nfreq);
tic
for ii = 1 : Nbeams
% for ii = Nbeams/2
    currentPhi = sonarBeams(ii);
    currentTheta = 0;
    currentNormal = beamNormalsRotated(ii,:);
    cosAngle = -(currentNormal(1)*targetNormals(:,1) ...
        + currentNormal(2)*targetNormals(:,2) ...
        + currentNormal(3)*targetNormals(:,3));
    
    lambertTerm = targetMu.*cosAngle.^2;
    % if the angle is negative, then it's zero - this is fake shadowing.
    % This is already done by occlusion checking, I thik
    lambertTerm = step(cosAngle).*lambertTerm; 
%     lambertTerm = ones(size(lambertTerm));
    propagationTerm = 1./R.^2.*exp(-2.*kpp.*R); %two-way loss of pressure (not intensity). 
    phiBeam = phiPrime - currentPhi;
    validInds = abs(phiBeam) <= sonarBeamwidth;
    % how many scatterers do we ahve
    nValid = sum(validInds);
    thetaBeam = thetaPrime - currentTheta;
    directivityTerm = beamPatternHorizontal(phiBeam(validInds));
    
    amplitude =randAmps(validInds)*s0.*propagationTerm(validInds).*directivityTerm.* ...
        sqrt(lambertTerm(validInds).*targetAreas(validInds)); % sqrt is because the lambert model is for intensity, and we are working with pressure.

    pFreq = sum(exp(1i.*2*R(validInds).*k_vec).*amplitude.*windowFreq,1); %add up all terms and time delays; can take a ton of memory...
    % put this in the order expected by matlab ifft
%     pFreq = [pFreq(Nfreq/2:end) pFreq(1:Nfreq/2-1)];
    pTime = fft(pFreq,[],2).*delta_f;
%     pTime = fftshift(pTime,2);
    p_2(ii,:) = pTime;
%     disp([ii Nbeams])
end
ctime2 = toc;
p_2_corrected = zeros(size(p_2));
tic
for ii = 1 : Nbeams
    weights = beamPatternHorizontal(sonarBeams - sonarBeams(ii));
    p_2_corrected(ii,:) = sum(weights.*p_2./sqrt(sum(weights.^2)),1);
end
ctime2a = toc;
fprintf('Beam culling: %g s.\n',ctime2);
fprintf('Correction to Beam culling: %g s.\n',ctime2a);

%% Method 3, incoherently average scatterers that fall in the same time sample
% This could also be called a form of rasterization for the rays.
% Right now, it's not faster than method 2, but this is probably because of
% the small problem size (number of scatterers).
%
% Also, matlab is quite slow at for-loops, so the inner loop below
% (for iii = 1 : Ntime)
% may be much faster in C/C++
%
% This method does not use the random amplitudes generated above, since it
% computes an incoherent average, which destroys the phase of each
% scattere. Instead, I create a new one on the fly for each beam and each
% sample in time. This is not a great way to do it, but it may generate
% realistic images.
%
% Derek: TODO: the amplitudes of this problem don't line up very well
% might need a fix to correct the rms pressure later
%
p_3 = zeros(Nbeams , Nfreq);
tic
for ii = 1 : Nbeams
% for ii = Nbeams/2
    currentPhi = sonarBeams(ii);
    currentTheta = 0;
    currentNormal = beamNormalsRotated(ii,:);
    cosAngle = -(currentNormal(1)*targetNormals(:,1) ...
        + currentNormal(2)*targetNormals(:,2) ...
        + currentNormal(3)*targetNormals(:,3));
    
    lambertTerm = targetMu.*cosAngle.^2;
    % if the angle is negative, then it's zero - this is fake shadowing.
    % This is already done by occlusion checking, I thik
    lambertTerm = step(cosAngle).*lambertTerm; 
%     lambertTerm = ones(size(lambertTerm));
    propagationTerm = 1./R.^2.*exp(-2.*kpp.*R); %two-way loss of pressure (not intensity). 
    phiBeam = phiPrime - currentPhi;
    validInds = abs(phiBeam) <= sonarBeamwidth;
    % how many scatterers do we ahve
    nValid = sum(validInds);
    thetaBeam = thetaPrime - currentTheta;
    directivityTerm = beamPatternHorizontal(phiBeam(validInds));
    
    amplitude =s0.*propagationTerm(validInds).*directivityTerm.* ...
        sqrt(lambertTerm(validInds).*targetAreas(validInds)); % sqrt is because the lambert model is for intensity, and we are working with pressure.
    timeDelay = 2.*R(validInds)./c;
    for iii = 1 : Ntime
        timeInds = timeDelay > (iii-1)*delta_t & timeDelay <= iii*delta_t;
        rmsP = sum(amplitude(timeInds));
%         numInds(ii,iii) = sum(timeInds);
        p_3(ii,iii) = rmsP.*(randn + 1i*randn)./sqrt(2);
    end
    
    p_3(ii,:) = ifft(fft(p_3(ii,:),[],2).*windowFreq,[],2);
%     pFreq = sum(exp(1i.*2*R(validInds).*k_vec).*amplitude.*windowFreq,1); %add up all terms and time delays; can take a ton of memory...
    % put this in the order expected by matlab ifft
%     pFreq = [pFreq(Nfreq/2:end) pFreq(1:Nfreq/2-1)];
%     pTime = fft(pFreq,[],2).*delta_f;
%     pTime = fftshift(pTime,2);
%     p_2(ii,:) = pTime;
%     disp([ii Nbeams])
end
ctime3 = toc;
p_3_corrected = zeros(size(p_3));
tic
for ii = 1 : Nbeams
    weights = beamPatternHorizontal(sonarBeams - sonarBeams(ii));
    p_3_corrected(ii,:) = sum(weights.*p_3./sqrt(sum(weights.^2)),1);
end
ctime2a = toc;
fprintf('Beam culling & raster : %g s.\n',ctime2);
fprintf('Correction to raster: %g s.\n',ctime2a);

%% Plot the three methods
figure(1)
clf
scatterPointSize = 10;
maxP = max(abs(p_1(:)));
% plot 60 dB of dynamic range on a common scale
clims = [-60 0] + 20*log10(maxP);
subplot(2,4,1)
pcolor(t_vector*c/2 , sonarBeams*180/pi, 10*log10(abs((p_1)).^2))
caxis(clims)
colormap(hot)
shading flat
xlabel('Range [m]')
ylabel('Beam Angle [Deg]')
title('Method 1: \phi - R')
h = colorbar;
ylabel(h,'Echo Level')
% colormap(viridis)

x = range_vector.*cos(sonarBeams);
y = range_vector.*sin(sonarBeams);

subplot(2,4,5)
scatter(x(:),y(:),scatterPointSize,20*log10(abs(p_1(:))),'filled')
caxis(clims)

colorbar
title('Method 1: x-y')  
xlabel('X [m]')
ylabel('Y [Deg]')
xlim(1.02*[-maxRange maxRange])
ylim(1.02*[-maxRange maxRange])
h = colorbar;
ylabel(h,'Echo Level')
set(gca,'Color','k')
axis equal
axis tight

subplot(2,4,2)
pcolor(t_vector*c/2 , sonarBeams*180/pi, 10*log10(abs((p_2)).^2))
% caxis(clims)
shading flat
xlabel('Range [m]')
ylabel('Beam Angle [Deg]')
title('Method 2: \phi - R')
h = colorbar;
ylabel(h,'Echo Level')
colormap(hot)

x = range_vector.*cos(sonarBeams);
y = range_vector.*sin(sonarBeams);

subplot(2,4,6)
scatter(x(:),y(:),scatterPointSize,20*log10(abs(p_2(:))),'filled')
% caxis(clims)

colorbar
title('Method 2: x-y')
xlabel('X [m]')
ylabel('Y [Deg]')
xlim(1.02*[-maxRange maxRange])
ylim(1.02*[-maxRange maxRange])
h = colorbar;
ylabel(h,'Echo Level')
set(gca,'Color','k')
axis equal
axis tight


subplot(2,4,3)
pcolor(t_vector*c/2 , sonarBeams*180/pi, 10*log10(abs((p_2_corrected)).^2))
caxis(clims)
shading flat
xlabel('Range [m]')
ylabel('Beam Angle [Deg]')
title('Method 2 (corrected): \phi - R')
h = colorbar;
ylabel(h,'Echo Level')
colormap(hot)




subplot(2,4,7)
scatter(x(:),y(:),scatterPointSize,20*log10(abs(p_2_corrected(:))),'filled')
caxis(clims)

colorbar
title('Method 2 (corrected): x-y')
xlabel('X [m]')
ylabel('Y [Deg]')
xlim(1.02*[-maxRange maxRange])
ylim(1.02*[-maxRange maxRange])
h = colorbar;
ylabel(h,'Echo Level')
set(gca,'Color','k')
axis equal
axis tight

subplot(2,4,4)
pcolor(t_vector*c/2 , sonarBeams*180/pi, 10*log10(abs((p_3_corrected)).^2))
caxis([-60 0] + 20*log10(max(abs(p_3_corrected(:)))));
colormap(hot)
shading flat
xlabel('Range [m]')
ylabel('Beam Angle [Deg]')
title('Method 3 Raster in time: \phi - R')
h = colorbar;
ylabel(h,'Echo Level')
% colormap(viridis)

x = range_vector.*cos(sonarBeams);
y = range_vector.*sin(sonarBeams);

subplot(2,4,8)
scatter(x(:),y(:),scatterPointSize,20*log10(abs(p_3_corrected(:))),'filled')
% caxis(clims)
caxis([-60 0] + 20*log10(max(abs(p_3_corrected(:)))));
colorbar
title('Method 3 Raster in time: x-y')
xlabel('X [m]')
ylabel('Y [Deg]')
xlim(1.02*[-maxRange maxRange])
ylim(1.02*[-maxRange maxRange])
h = colorbar;
ylabel(h,'Echo Level')
set(gca,'Color','k')
axis equal
axis tight

% print('-dpng',['ScatteringModelComparisons_Dave' date])
