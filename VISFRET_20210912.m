clc, clear all, close all
fprintf('\n');
global P_flu_emit QY1 absorbedincident P_exc_abs detFraction P_in
%% Description

for iter=1:1
%% G
model = initializeMCmatlabModel();
lambda=456;
model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication

model.G.nx                = 200; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = .8; % [cm] x size of simulation cuboid
model.G.Ly                = .4; % [cm] y size of simulation cuboid
model.G.Lz                = .4; % [cm] z size of simulation cuboid

model.G.mediaPropParams   = lambda; % Cell array containing any additional parameters to be passed to the getMediaProperties function
model.G.mediaPropertiesFunc = @mediaProp_FRET1; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @GeometryDefinition_FRET1; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

% Execution, do not modify the next line:
model = defineGeometry(model);
plotMCmatlabGeom(model);

%% MC
model = clearMCmatlabModel(model,'MC'); % Only necessary if you want to run this section repeatedly, re-using previous G data
model.MC.useGPU = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
model.MC.nPhotonsRequested        = 1e7; % # of photons to launch

model.MC.silentMode               = false; % (Default: false) Disables command window text and progress indication
% model.MC.useAllCPUs               = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
% model.MC.calcNFR                  = true; % (Default: true) If true, the 3D fluence rate output array NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
% model.MC.calcNFRdet               = false; % (Default: false) If true, the 3D fluence rate output array NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
model.MC.nExamplePaths            = 0; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes
% model.MC.farFieldRes              = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.

model.MC.matchedInterfaces        = true; % (Default: true) If true, assumes all refractive indices are 1. If false, uses the refractive indices defined in getMediaProperties
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
model.MC.wavelength               = 450; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.beam.beamType            = 5; % 0: Pencil beam, 1: Isotropically emitting point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.beam.NF.XDistr           = 0; % X near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.XWidth           = .2; % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.NF.YDistr           = 0; % Y near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
model.MC.beam.NF.YWidth           = .2; % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
model.MC.beam.FF.XDistr           = 2; % X far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.XWidth           = pi/8; % [rad] X far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.FF.YDistr           = 2; % Y far field distribution - 0: Top-hat, 1: Gaussian, 2: Cosine (Lambertian), Array: Custom. Doesn't need to be normalized.
model.MC.beam.FF.YWidth           = pi/8; % [rad] Y far field 1/e^2 half-angle if top-hat or Gaussian or half-angle of the full distribution if custom
model.MC.beam.psi                 = 0; % [rad] (Default: 0) Axial rotation angle of beam, relevant only for XY distributed beams

model.MC.beam.theta               = 0; % [rad] Polar angle of beam center axis
model.MC.beam.phi                 = 0; % [rad] Azimuthal angle of beam center axis



model.MC.useLightCollector      = true;
model.MC.LC.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.LC.y         = 0; % [cm] y position
model.MC.LC.z         = 0.00; % [cm] z position

model.MC.LC.theta     = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.LC.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.LC.f         = Inf; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.LC.diam      = .17; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.LC.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.LC.NA        = 0.4; % [-] Fiber NA. Only used for infinite f.

model.MC.LC.res       = 1; %50 X and Y resolution of light collector in pixels, only used for finite f

% model.MC.LC.tStart    = -1e-13; % [s] Start of the detection time interval
% model.MC.LC.tEnd      = 5e-12; % [s] End of the detection time interval
% model.MC.LC.nTimeBins = 30; % (Default: 0) Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.

% Execution, do not modify the next line:
model.MC.useGPU                   = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
model = runMonteCarlo(model);
plotMCmatlab(model);


cmax=80;
for c=1:cmax

lambda=round(456+(c/cmax*254)); %456

        %% Donor
        %% G
        
        
        model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication
        
        model.G.nx                = 200; % Number of bins in the x direction
        model.G.ny                = 100; % Number of bins in the y direction
        model.G.nz                = 100; % Number of bins in the z direction
        model.G.Lx                = .8; % [cm] x size of simulation cuboid
        model.G.Ly                = .4; % [cm] y size of simulation cuboid
        model.G.Lz                = .4; % [cm] z size of simulation cuboid
        
        model.G.mediaPropParams   = lambda; % Cell array containing any additional parameters to be passed to the getMediaProperties function
        model.G.mediaPropertiesFunc = @mediaProp_FRET1; % Media properties defined as a function at the end of this file
        model.G.geomFunc          = @GeometryDefinition_FRET1; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
        model.G.geomFuncParams    = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
        
        % Execution, do not modify the next line:
        model = defineGeometry(model);
        plotMCmatlabGeom(model);
        %% FMC
        global P_flu_emit QY1 absorbedincident P_exc_abs detFraction P_in
        model = clearMCmatlabModel(model,'FMC');
        
        model.FMC.useGPU = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
        
        
        model.FMC.nPhotonsRequested        = 1e7; % # of photons to launch
        
        model.FMC.silentMode               = false; % (Default: false) Disables command window text and progress indication
        % model.FMC.useAllCPUs             = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
        % model.FMC.calcNFR                = true; % (Default: true) If true, the 3D fluence rate output array NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
        % model.FMC.calcNFRdet             = false; % (Default: false) If true, the 3D fluence rate output array NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
        % model.FMC.nExamplePaths          = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes
        % model.FMC.farFieldRes            = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.
        
        model.FMC.matchedInterfaces        = true; % (Default: true) If true, assumes all refractive indices are 1. If false, uses the refractive indices defined in getMediaProperties
        model.FMC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
        model.FMC.wavelength               = lambda; % [nm] Excitation wavelength, used for determination of optical properties for excitation light. Should be emission?
        
        model.FMC.useLightCollector      = true;
        model.FMC.LC.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
        model.FMC.LC.y         = 0; % [cm] y position
        model.FMC.LC.z         = 0.00; % [cm] z position
        
        model.FMC.LC.theta     = 0; % [rad] Polar angle of direction the light collector is facing
        model.FMC.LC.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing
        
        model.FMC.LC.f         = Inf; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
        model.FMC.LC.diam      = .17; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
        model.FMC.LC.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
        model.FMC.LC.NA        = 0.4; % [-] Fiber NA. Only used for infinite f.
        
        model.FMC.LC.res       = 1; % X and Y resolution of light collector in pixels, only used for finite f
        
        % Execution, do not modify the next line:
        
        model = runMonteCarlo(model,'fluorescence');
        plotMCmatlab(model,'fluorescence');
        if c>2
            model.MC.mediaProperties_funcHandles(4).Y =model.FMC.mediaProperties_funcHandles(4).Y;
            %model.MC.mediaProperties_funcHandles =model.FMC.mediaProperties_funcHandles;
        end
        
        detectord1(c,1)=lambda;
        detectord1(c,2)=absorbedincident*P_flu_emit/(P_exc_abs)*(detFraction/P_in)/100;
        datad(c,1)=lambda;
        datad(c,2)=absorbedincident;
        datad(c,3)=P_flu_emit/(P_exc_abs);
        datad(c,4)=(detFraction/P_in)/100;
        datad(c,5)=model.MC.mediaProperties(4).Y;
        datad(c,6)=model.FMC.mediaProperties(2).mua;
        datad(c,7)=model.FMC.mediaProperties(1).mua;
    end
    %% Acceptor

for c=1:cmax
lambda=round(503+(c/cmax*197));
        model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication
        
        model.G.nx                = 200; % Number of bins in the x direction
        model.G.ny                = 100; % Number of bins in the y direction
        model.G.nz                = 100; % Number of bins in the z direction
        model.G.Lx                = .8; % [cm] x size of simulation cuboid
        model.G.Ly                = .4; % [cm] y size of simulation cuboid
        model.G.Lz                = .4; % [cm] z size of simulation cuboid
        
        model.G.mediaPropParams   = lambda; % Cell array containing any additional parameters to be passed to the getMediaProperties function
        model.G.mediaPropertiesFunc = @mediaProp_FRET2; % Media properties defined as a function at the end of this file
        model.G.geomFunc          = @GeometryDefinition_FRET2; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
        model.G.geomFuncParams    = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
        
        % Execution, do not modify the next line:
        model = defineGeometry(model);
        plotMCmatlabGeom(model);
        
        %% FMC
        global P_flu_emit QY1 absorbedincident P_exc_abs detFraction P_in
        model = clearMCmatlabModel(model,'FMC');
        
        model.FMC.useGPU = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
        
        
        model.FMC.nPhotonsRequested        = 1e7; % # of photons to launch
        
        model.FMC.silentMode               = false; % (Default: false) Disables command window text and progress indication
        % model.FMC.useAllCPUs             = true; % (Default: false) If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
        % model.FMC.calcNFR                = true; % (Default: true) If true, the 3D fluence rate output array NFR will be calculated. Set to false if you have a light collector and you're only interested in the image output.
        % model.FMC.calcNFRdet             = false; % (Default: false) If true, the 3D fluence rate output array NFRdet will be calculated. Only photons that end up on the light collector are counted in NFRdet.
        % model.FMC.nExamplePaths          = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes
        % model.FMC.farFieldRes            = 50; % (Default: 0) If nonzero, photons that "escape" will have their energies tracked in a 2D angle distribution (theta,phi) array with theta and phi resolutions equal to this number. An "escaping" photon is one that hits the top cuboid boundary (if boundaryType == 2) or any cuboid boundary (if boundaryType == 1) where the medium has refractive index 1.
        
        model.FMC.matchedInterfaces        = true; % (Default: true) If true, assumes all refractive indices are 1. If false, uses the refractive indices defined in getMediaProperties
        model.FMC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping
        model.FMC.wavelength               = lambda; % [nm] Excitation wavelength, used for determination of optical properties for excitation light. Should be emission?
        
        model.FMC.useLightCollector      = true;
        model.FMC.LC.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
        model.FMC.LC.y         = 0; % [cm] y position
        model.FMC.LC.z         = 0.00; % [cm] z position
        
        model.FMC.LC.theta     = 0; % [rad] Polar angle of direction the light collector is facing
        model.FMC.LC.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing
        
        model.FMC.LC.f         = Inf; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
        model.FMC.LC.diam      = .17; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
        model.FMC.LC.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
        model.FMC.LC.NA        = 0.4; % [-] Fiber NA. Only used for infinite f.
        
        model.FMC.LC.res       = 1; % X and Y resolution of light collector in pixels, only used for finite f
        
        % Execution, do not modify the next line:
        
        model = runMonteCarlo(model,'fluorescence');
        plotMCmatlab(model,'fluorescence');
        if c>1
            model.MC.mediaProperties_funcHandles(4).Y =model.FMC.mediaProperties_funcHandles(4).Y;
            %model.MC.mediaProperties_funcHandles =model.FMC.mediaProperties_funcHandles;
        end
        detectora1(c,2)=absorbedincident*P_flu_emit/(P_exc_abs)*(detFraction/P_in)/100;
        detectora1(c,1)=lambda;
        dataa(c,1)=lambda;
        dataa(c,2)=absorbedincident;
        dataa(c,3)=P_flu_emit/(P_exc_abs);
        dataa(c,4)=(detFraction/P_in)/100;
        dataa(c,5)=model.MC.mediaProperties(4).Y;       
        dataa(c,6)=model.FMC.mediaProperties(2).mua;
        dataa(c,7)=model.FMC.mediaProperties(1).mua;
    end



%%

%% Post-processing
spectrumoverlap=.62; %.257 for new values; .2868 for old, %.486 for new new
r=3.5; %nm, 3.5 old value
r0= 5; %nm, forster radius
%E=(r0^6/(r0^6+r^6));
E=.803; %Fret efficiency
% Determine Intensities for donor and acceptor emissions reaching detector
clearvars FMCinput Ginput Goutput MCinput MCoutput
tic


for i=1:size(detectord1,1)
    Detectord(i,2)=detectord1(i,2)*(1-spectrumoverlap*E*.023*.932*3/2); %spectrum overlap, E=distance effect, _=% of bound Con-A, _=% of bound APTS-MT, 3/2=k^2 orientation factor
    Detectord(i,1)=detectord1(i,1);
   
end


for i=1:size(detectora1,1)
    Detectora(i,2)=detectora1(i,2)*(spectrumoverlap*E*58.275*.023*3/2);
    Detectora(i,1)=detectora1(i,1);
end
close all
% Processing
Detectoraa=LinearInterpolator(Detectora,0);
Detectordd=LinearInterpolator(Detectord,0);
Detector=[Detectoraa;Detectordd];
Detector1=sortrows(Detector,1,'ascend');
c=1;
%
%
for i=2:size(Detector1,1)-1
    if Detector1(i,1)==Detector1(i+1,1)
        Detector2(c,1)=Detector1(i,1);
        Detector2(c,2)=Detector1(i,2)+Detector1(i+1,2);
        c=c+1;
    elseif Detector1(i,1)==Detector(i-1,1)
        continue
    else Detector3(i,:)=Detector1(i,:);
    end
end

Detector4=[Detector2; Detector3(1:57,:)];
Detector4=sortrows(Detector4,1,'ascend');
Detector4( ~any(Detector4,2), : ) = [];
plot(Detector4(:,1),Detector4(:,2))
Detector4=unique(Detector4,'rows');
Detector4(1:9,:)=[];
Detector5(:,:,iter)=Detector4;

clearvars -except data1 data2 detectord1 detectora1 Detector4 spectrumoverlap E zztop Detector5

end

%% Geometry function(s)
function M = GeometryDefinition_FRET1(X,Y,Z,parameters)
M=2*ones(size(X));
M(Z<.01)=1;
M(Z>.16)=6;
M(X>-.325 & X<.325 & Y>-.035 & Y<.035 & Z>.155 & Z<.245)=3;
M(X>-.325+.2 & X<-.225+.2 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=4;
M(X>-.325+.5 & X<-.225+.5 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=4;

M(X>-.325+.35 & X<-.225+.35 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=5;
M(X>-.325+.05 & X<-.225+.05 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=5;
end

function M = GeometryDefinition_FRET2(X,Y,Z,parameters)
M=2*ones(size(X));
M(Z<.01)=1;
M(Z>.16)=6;
M(X>-.325 & X<.325 & Y>-.035 & Y<.035 & Z>.155 & Z<.245)=3;
M(X>-.325+.2 & X<-.225+.2 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=4;
M(X>-.325+.5 & X<-.225+.5 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=4;

M(X>-.325+.35 & X<-.225+.35 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=5;
M(X>-.325+.05 & X<-.225+.05 & Y>-.025 & Y<.025 & Z>.175 & Z<.225)=5;
end

%% Media Properties function
function mediaProperties = mediaProp_FRET1(wavelength,parameters)
load spectralLIB.mat
MU(:,1) = interp1(nmLIB,muaoxy,wavelength);
MU(:,2) = interp1(nmLIB,muadeoxy,wavelength);
MU(:,3) = interp1(nmLIB,muawater,wavelength);
MU(:,4) = interp1(nmLIB,muamel,wavelength);
Muafat = interp1(nmLIB,muafat,wavelength);
muabase = (7.84*10^8).*wavelength.^(-3.255);

j=1;
mediaProperties(j).name  = 'Epidermis'; %(thickness: https://www.ncbi.nlm.nih.gov/pubmed/14690333)
B = 0;
S = 0.75;
W = 0.75;
Me=.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = 0.1*(MU*X+muabase);
mediaProperties(j).mus =  musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 2;%1.3;
mediaProperties(j).VHC = 3391*1.109e-3;
mediaProperties(j).TC  = 0.37e-2;

j=2;
mediaProperties(j).name = 'Dermis';
B = 0.002;
S = 0.67;
W = 0.65;
Me = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua =  MU*X+muabase;
mediaProperties(j).mus =  musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;

j=3;
mediaProperties(j).name  = 'PEGDA Hydrogel';
mediaProperties(j).mua   = 1e-6; %10
mediaProperties(j).mus   = .0005;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.4;

j=4;
a2=xlsread('APTS.xlsx');
b1=a2(:,1)-parameters;
c1=find(b1(:,1)==0);
%Measuring the FRET response from the FRET LED- calculating FRET
mediaProperties(j).name  = 'APTS';%'FRET';
mediaProperties(j).mua   = .001716; 
mediaProperties(j).mus   = 1e-6;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.3;
mediaProperties(j).Y     =a2(c1,3)*.2;%

j=5;
mediaProperties(j).name  = 'Phosphorescence Lifetime Chemistry';%'PdTPTBP';
mediaProperties(j).mua   = 168;%.173; %(nm (M-1 cm-1))630nm
mediaProperties(j).mus   = 1e-6;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.4;
mediaProperties(j).Y     =0;%0.21*a1(c,3);%*a1(c,3);%.21= quantum yield

j=6;
mediaProperties(j).name  = 'Subcutaneous Fat';
B = 0.05;
S = 0.60;
W = 0.7;
Me = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
gg      = 0.90;
musp = 1.08*1e8*wavelength.^(-2.525)+157.494*wavelength.^(-0.345);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = MU*X;
mediaProperties(j).mus =  musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
end

function mediaProperties = mediaProp_FRET2(wavelength,parameters)
load spectralLIB.mat
MU(:,1) = interp1(nmLIB,muaoxy,wavelength);
MU(:,2) = interp1(nmLIB,muadeoxy,wavelength);
MU(:,3) = interp1(nmLIB,muawater,wavelength);
MU(:,4) = interp1(nmLIB,muamel,wavelength);
Muafat = interp1(nmLIB,muafat,wavelength);
muabase = (7.84*10^8).*wavelength.^(-3.255);

j=1;
mediaProperties(j).name  = 'Epidermis'; %(thickness: https://www.ncbi.nlm.nih.gov/pubmed/14690333)
B = 0;
S = 0.75;
W = 0.75;
Me=.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = 0.1*(MU*X+muabase);
mediaProperties(j).mus = musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
mediaProperties(j).VHC = 3391*1.109e-3;
mediaProperties(j).TC  = 0.37e-2;

j=2;
mediaProperties(j).name = 'Dermis';
B = 0.002;
S = 0.67;
W = 0.65;
Me = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(wavelength/500).^-4 + (1-fray)*(wavelength/500).^-bmie);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = MU*X+muabase;
mediaProperties(j).mus = musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
mediaProperties(j).VHC = 3391*1.109e-3;
mediaProperties(j).TC  = 0.37e-2;

j=3;

mediaProperties(j).name  = 'PEGDA Hydrogel';
mediaProperties(j).mua   = 1e-6; %10
mediaProperties(j).mus   = .0005;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.4;

j=4;
a3=xlsread('TRITC.xlsx');
b2=a3(:,1)-(parameters);
c2=find(b2(:,1)==0);
%Measuring the FRET response from the FRET LED- calculating FRET
mediaProperties(j).name  = 'TRITC';%'FRET';
mediaProperties(j).mua   = .001716; 
mediaProperties(j).mus   = 1e-6;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.3;
mediaProperties(j).Y     =a3(c2,3)*.1;%

j=5;
mediaProperties(j).name  = 'PdTCPP';%'PdTPTBP';
mediaProperties(j).mua   = 168;%.173; %(nm (M-1 cm-1))630nm
mediaProperties(j).mus   = 1e-6;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.4;
mediaProperties(j).Y     =0;%.21*a1(c,3);%*a1(c,3);%.21= quantum yield

j=6;
mediaProperties(j).name  = 'Fat';
B = 0.05;
S = 0.60;
W = 0.7;
Me = 0;
musp500 = 20;
fray    = 0.2;
bmie    = 1.0;
gg      = 0.90;
musp = 1.08*1e8*wavelength.^(-2.525)+157.494*wavelength.^(-0.345);
X = [B*S B*(1-S) W Me]';
mediaProperties(j).mua = MU*X;
mediaProperties(j).mus = musp/(1-gg);
mediaProperties(j).g   = gg;
mediaProperties(j).n   = 1.3;
end
