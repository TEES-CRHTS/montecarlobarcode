function [ppix] = Autofluorescence_450nm_ppix_Feb21(nphoton,resolution,lambda1,lambda2)
min=583;
max=699;
fprintf('\n');
global P_flu_emit QY1 absorbedincident P_exc_abs detFraction P_in
%% Description
lambda=lambda1;
for r=1:1
    
    
    %% Geometry definition
    model = initializeMCmatlabModel();
    
    model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication
    
    model.G.nx                = 200; % Number of bins in the x direction
    model.G.ny                = 100; % Number of bins in the y direction
    model.G.nz                = 100; % Number of bins in the z direction
    model.G.Lx                = .8; % [cm] x size of simulation cuboid
    model.G.Ly                = .4; % [cm] y size of simulation cuboid
    model.G.Lz                = .4; % [cm] z size of simulation cuboid
    
    model.G.mediaPropParams   = 585; % Cell array containing any additional parameters to be passed to the getMediaProperties function
    model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
    model.G.geomFunc          = @GeometryDefinition_PhosBarcode2; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
    model.G.geomFuncParams    = {0.03}; % Cell array containing any additional parameters to pass into the geometry function, such as media depths, inhomogeneity positions, radii etc.
    
    % Execution, do not modify the next line:
    model = defineGeometry(model);
    plotMCmatlabGeom(model);
    
    % Monte Carlo simulation
    model = clearMCmatlabModel(model,'MC'); % Only necessary if you want to run this section repeatedly, re-using previous G data
    model.MC.useGPU = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
    model.MC.nPhotonsRequested        = nphoton; % # of photons to launch
    
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
    model.MC.beam.NF.XWidth           = .02; % [cm] X near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
    model.MC.beam.NF.YDistr           = 0; % Y near field distribution - 0: Top-hat, 1: Gaussian, Array: Custom. Doesn't need to be normalized.
    model.MC.beam.NF.YWidth           = .01; % [cm] Y near field 1/e^2 radius if top-hat or Gaussian or half-width of the full distribution if custom
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
    model.MC.LC.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
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
    counter=1;
    while lambda <lambda2
        if counter>1
            lambda=lambda+resolution;
        end
        if lambda<min | lambda>max
            data(counter,1)=lambda;
            data(counter,2)=0;
            data(counter,3)=0;
            data(counter,4)=(0);
            data(counter,5)=0;
            data(counter,6)=0;
            ppix=data;
            counter=counter+1;
            continue
        end
        model.G.silentMode        = false; % (Default: false) Disables command window text and progress indication
        
        model.G.nx                = 200; % Number of bins in the x direction
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
        model = clearMCmatlabModel(model,'FMC');
        
        model.FMC.useGPU = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
        

        model.FMC.nPhotonsRequested        = nphoton; % # of photons to launch
        
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
        model.FMC.LC.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
        model.FMC.LC.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
        model.FMC.LC.NA        = 0.4; % [-] Fiber NA. Only used for infinite f.
        
        model.FMC.LC.res       = 1; % X and Y resolution of light collector in pixels, only used for finite f
        
        % Execution, do not modify the next line:
        model = runMonteCarlo(model,'fluorescence');
        plotMCmatlab(model,'fluorescence');
        
        if counter>1
            model.MC.mediaProperties_funcHandles(4).Y =model.FMC.mediaProperties_funcHandles(4).Y;
            %model.MC.mediaProperties_funcHandles =model.FMC.mediaProperties_funcHandles;
        end
        
        data(counter,1)=lambda;
        data(counter,2)=absorbedincident;
        data(counter,3)=P_flu_emit;
        data(counter,4)=(detFraction)/100;
        data(counter,5)=model.MC.mediaProperties_funcHandles(4).Y;
        data(counter,6)=(absorbedincident)*(P_flu_emit/P_exc_abs)*(detFraction/P_in)/100;
        ppix=data;
        counter=counter+1;
        %% Post-processing
        
    end
end
%% Geometry function(s)
    function M = GeometryDefinition_PhosBarcode2(X,Y,Z,parameters)
        M=2*ones(size(X));
        M(Z<.01)=1;
        M(Z>.16)=3;
        M(Z<.01 & Z>.002)=4;
    end

%% Media Properties function
    function mediaProperties = mediaPropertiesFunc(wavelength,parameters)
        test=load('spectralLIB.mat');
        MU(:,1) = interp1(test.nmLIB,test.muaoxy,wavelength);
        MU(:,2) = interp1(test.nmLIB,test.muaoxy,wavelength);
        MU(:,3) = interp1(test.nmLIB,test.muaoxy,wavelength);
        MU(:,4) = interp1(test.nmLIB,test.muaoxy,wavelength);
        Muafat = interp1(test.nmLIB,test.muaoxy,wavelength);
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
        % mediaProperties(j).mua = MU*X;
        mediaProperties(j).mua = Muafat;
        mediaProperties(j).mus = musp/(1-gg);
        mediaProperties(j).g   = gg;
        mediaProperties(j).n   = 1.3;
        
        j=4;
        mediaProperties(j).name = 'PPIX';
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
        a1=xlsread('ppix_feb21.xlsx');
        b=a1(:,1)-parameters;
        c=find(b(:,1)==0);
        mediaProperties(j).mua  = .046*2.3; %MU*X; .08%
        mediaProperties(j).mus  = musp/(1-gg);
        mediaProperties(j).g    = gg;
        mediaProperties(j).n    = 1.3;
        mediaProperties(j).Y    = .085*a1(c,3); %.06
    end
end
