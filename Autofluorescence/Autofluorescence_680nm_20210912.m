%clc, clear all, close all
set(0,'DefaultFigureWindowStyle','docked')

nphoton=1e7; %photons
resolution=2; %nm spacing between output wavelengths
lambda1=690;%nm
lambda2= 740;%nm
for r=1:1
    Melanin=Autofluorescence_680nm_Melanin_Aug20(nphoton,resolution,lambda1,lambda2);
    
 
    
    Mel(:,1)=Melanin(:,1);
    Mel(:,2)=Melanin(:,6);
    
    %%

    Mel=LinearInterpolator(Mel,0);

    %%
    autofluorescence(:,1,r)=Mel(:,1);
    autofluorescence(:,2,r)=Mel(:,2);

    
    clearvars FAD Col Mel ppix VitA

        
    end
