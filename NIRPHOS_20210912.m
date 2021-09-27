clc, clear all, close all
fprintf('\n');
global P_flu_emit QY1 absorbedincident P_exc_abs detFraction P_in
%% Description
set(0,'DefaultFigureVisible','off');

for r=1:3
    tic
    %% Geometry definition
    model = initializeMCmatlabModel();
    
    model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication
    
    model.G.nx                = 100; % Number of bins in the x direction
    model.G.ny                = 100; % Number of bins in the y direction
    model.G.nz                = 100; % Number of bins in the z direction
    model.G.Lx                = .8; % [cm] x size of simulation cuboid
    model.G.Ly                = .4; % [cm] y size of simulation cuboid
    model.G.Lz                = .4; % [cm] z size of simulation cuboid
    
    model.G.mediaPropParams   = 700; % Cell array containing any additional parameters to be passed to the getMediaProperties function
    model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
    model.G.geomFunc          = @GeometryDefinition_PhosBarcode2; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
    model.G.geomFuncParams    = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
    
    % Execution, do not modify the next line:
    model = defineGeometry(model);
    plotMCmatlabGeom(model);
    
    % Monte Carlo simulation
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
    model.MC.wavelength               = 680; % [nm] Excitation wavelength, used for determination of optical properties for excitation light
    
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
    
    
    %%
    
    cmax=80;
    for c=1:cmax
        lambda=round(700+(150/cmax*c))-1; %342
        model = clearMCmatlabModel(model,'G');
        model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication
        
        model.G.nx                = 100; % Number of bins in the x direction
        model.G.ny                = 100; % Number of bins in the y direction
        model.G.nz                = 100; % Number of bins in the z direction
        model.G.Lx                = .8; % [cm] x size of simulation cuboid
        model.G.Ly                = .4; % [cm] y size of simulation cuboid
        model.G.Lz                = .4; % [cm] z size of simulation cuboid
        
        model.G.mediaPropParams   = lambda; % Cell array containing any additional parameters to be passed to the getMediaProperties function
        model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
        model.G.geomFunc          = @GeometryDefinition_PhosBarcode2; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
        model.G.geomFuncParams    = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
        
        % Execution, do not modify the next line:
        model = defineGeometry(model);
        plotMCmatlabGeom(model);
        
        %
        
        model.FMC.useGPU = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
        %model.FMC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
        model.FMC.nPhotonsRequested        = 1e7; % # of photons to launch
        model.FMC.silentMode               = false; % (Default: false) Disables command window text and progress indication
        
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
        if c>1
            model.MC.mediaProperties_funcHandles(5).Y =model.FMC.mediaProperties_funcHandles(5).Y;
        end
        model = runMonteCarlo(model,'fluorescence');
        plotMCmatlab(model,'fluorescence');
        
        data(c,2,r)=absorbedincident*P_flu_emit/(P_exc_abs)*(detFraction/P_in)/100;
        data(c,1,r)=lambda;
        
        
    end
    
end
toc

%% Geometry function(s)
% A geometry function takes as input X,Y,Z matrices as returned by the
% "ndgrid" MATLAB function as well as any parameters the user may have
% provided in the definition of model.G. It returns the media matrix M,
% containing numerical values indicating the media type (as defined in
% mediaPropertiesFunc) at each voxel location.
function M = GeometryDefinition_PhosBarcode2(X,Y,Z,parameters)
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
% The media properties function defines all the optical and thermal
% properties of the media involved by constructing and returning a
% "mediaProperties" struct with various fields. As its input, the function
% takes the wavelength as well as any other parameters you might specify
% above in the model file, for example parameters that you might loop over
% in a for loop. Dependence on excitation fluence rate FR, temperature T or
% fractional heat damage FD can be specified as in examples 12-15.
function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
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
%guesstimated values **DO NOT USE**
mediaProperties(j).name  = 'PEGDA Hydrogel';
mediaProperties(j).mua   = 1e-6; %10
mediaProperties(j).mus   = .0005;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.4;

j=4;
mediaProperties(j).name  = 'FRET Sensing Chemistry- Phosphorescence';
mediaProperties(j).mua   = 250; %.0002.*wavelength.^2-.3449.*wavelength+175.36;
mediaProperties(j).mus   = 1e-6;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.4;
mediaProperties(j).Y     =0;

j=5;
a1=xlsread('PDBMAP3_Aug14.xlsx');
b=a1(:,1)-(parameters);
c=find(b(:,1)==0);
mediaProperties(j).name  = 'PdBP';
mediaProperties(j).mua   = 2.817*4.1764;
mediaProperties(j).mus   = 1e-6;
mediaProperties(j).g     = .75;
mediaProperties(j).n     = 1.4;
mediaProperties(j).Y     =.21*a1(c,3);%*a1(c,3);%.21= quantum yield

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