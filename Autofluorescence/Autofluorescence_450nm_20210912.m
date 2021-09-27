%clc, clear all, close all
set(0,'DefaultFigureWindowStyle','docked')

nphoton=1e7; %photons
resolution=5; %nm spacing between output wavelengths
lambda1=465;%nm
lambda2= 680;%nm
for r=1:2
    FAD1=Autofluorescence_450nm_FAD_May4(nphoton,resolution,lambda1,lambda2);
    Collagen=Autofluorescence_450nm_Collagen2_May4(nphoton,resolution,lambda1,lambda2);
    PPIX=Autofluorescence_450nm_PPIX_May4(nphoton,resolution,lambda1,lambda2);
    VitaminA=Autofluorescence_450nm_VitaminA_May4(nphoton,resolution,lambda1,lambda2);
    Melanin=Autofluorescence_450nm_Melanin_May4(nphoton,resolution,lambda1,lambda2);
    
    %%
    
    FAD(:,1)=FAD1(:,1);
    FAD(:,2)=FAD1(:,6);
    
    Col(:,1)=Collagen(:,1);
    Col(:,2)=Collagen(:,6);
    
    Mel(:,1)=Melanin(:,1);
    Mel(:,2)=Melanin(:,6);
    
    ppix(:,1)=PPIX(:,1);
    ppix(:,2)=PPIX(:,6);
    
    VitA(:,1)=VitaminA(:,1);
    VitA(:,2)=VitaminA(:,6);
    %%
    FAD=LinearInterpolator(FAD,0);
    Mel=LinearInterpolator(Mel,0);
    ppix=LinearInterpolator(ppix,0);
    VitA=LinearInterpolator(VitA,0);
    Col=LinearInterpolator(Col,0);
    %%
    autofluorescence(:,1,r)=linspace(lambda1,lambda2,(lambda2-lambda1)+1);
    autofluorescence(:,2,r)=FAD(:,2);
    autofluorescence(:,3,r)=Mel(:,2);
    autofluorescence(:,4,r)=Col(:,2);
    autofluorescence(:,5,r)=VitA(:,2);
    autofluorescence(:,6,r)=ppix(:,2);
    
    for i=1:size(autofluorescence,1)
        autofluorescence(i,7,r)=sum(autofluorescence(i,2:6,r));
    end
    
    clearvars FAD Col Mel ppix VitA
    if r==1
        beep
    end
end