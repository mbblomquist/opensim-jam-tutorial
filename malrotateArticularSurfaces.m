function [done,femMalrot,tibMalrot] = malrotateArticularSurfaces( Params , iSim )
% -------------------------------------------------------------------------
% this script is used to vary V-V and I-E positions of articular surfaces
% 
% Overview:
% 
% 1. Load OBJ file of components
% 2. Translate point clouds to desired origin for transformation
% 3. Transform point clouds
% 4. Translate point clouds back to original origin
% 5. Write new OBJ file of components
% 
% 
% WRITTEN BY JOSH ROTH
% v1 (JDR): Initial Release
% v2 (MBB): Compatible with OpenSim
% v3 (MBB): Updated to fit into OpenSim_JAM pipeline
% -------------------------------------------------------------------------
%% Parse inputs
%==============

filePathOut = Params.fileOut ;
femName = Params.femImplant ;
tibName = Params.tibImplant ;
simulationProg = Params.simulationProg ;
fem = Params.Fem ; % femur properties
tib = Params.Tib ; % tibia properties
vvFemMalrot = Params.vvFemMalrot ;
ieFemMalrot = Params.ieFemMalrot ;
vvTibMalrot = Params.vvTibMalrot ;
ieTibMalrot = Params.ieTibMalrot ;

switch simulationProg
    case 'OpenSim'

        femPts = fem.pts;
        femTri = fem.tri;
        femNorm = fem.norm;

        tibPts = tib.pts;
        tibTri = tib.tri;
        tibNorm = tib.norm;
        
end

% Plot points
plotComponents = 0; % "1" to plot components at initial and final steps and "0" to not plot anything
plotTroubleshooting = 0; % "1" to plot components at all steps and "0" to not plot anything
    
%% APPLY OFFSET FOR TRANSFORMATION
%   translate *.stl points such that origin for tranformation is located at
%   the geometric center in AP and ML and the most distal in PD

% Plot stl points
if plotComponents == 1
    figure()
    plot3(femPts(:,1),femPts(:,2), femPts(:,3),'g.')
    hold on
    plot3(tibPts(:,1),tibPts(:,2), tibPts(:,3),'g.')
    axis equal tight
end


% Define offsets for femoral articular surface based on axis of malrotation
%-------------------------------------------------------------------------
if vvFemMalrot(iSim) > 0
    mlAxisFem = 'LatDist' ;
elseif vvFemMalrot(iSim) < 0
    mlAxisFem = 'MedDist' ;
else
    mlAxisFem = 'Center' ;
end

pdOffsetFem = min( femPts(:,2) ) ; % find most distal point

switch mlAxisFem
    case 'LatDist'
        latCompartmentIndicesFem = femPts(:,3)<0 ;
        latCompartmentFem = femPts(latCompartmentIndicesFem,:);
        [ ~ , latMostDistIndexFem ] = min(latCompartmentFem(:,2));
        mlOffsetFem = latCompartmentFem(latMostDistIndexFem,3); % ml offset = ml value of most distal point
        %     mlOffsetFem = (max(fem(:,3)) - min(fem(:,3)))/4; %rotation about center of lateral condyle
    case 'Center'
        mlOffsetFem = -(max(femPts(:,3)) + min(femPts(:,3)))/2; %rotation about center of component
    case 'MedDist'
        medCompartmentIndicesFem = femPts(:,3)>0 ;
        medCompartmentFem = femPts(medCompartmentIndicesFem,:);
        [ ~ , medMostDistIndexFem ] = min(medCompartmentFem(:,2));
        mlOffsetFem = medCompartmentFem(medMostDistIndexFem,3); % ml offset = ml value of most distal point
        %     mlOffsetFem = -(max(fem(:,3)) - min(fem(:,3)))/4; %rotation about center of medial condyle
end

% apOffsetFem = mean(femPts(:,1)); for rotating about AP center of femoral component
apOffsetFem = -min(femPts(:,1));

% Origin now max posterior, max distal, medial/laterally at max distal point
femCentered = [femPts(:,1)+apOffsetFem,femPts(:,2)-pdOffsetFem,femPts(:,3)+mlOffsetFem];

% Define offsets for tibial articular surface based on axis of malrotation
%-------------------------------------------------------------------------
if vvTibMalrot(iSim)>0
    mlAxisTib = 'LatDist';
%     mlAxis = 'LatEdge';
elseif vvTibMalrot(iSim)<0
    mlAxisTib = 'MedDist';
%     mlAxis = 'MedEdge';
else 
    mlAxisTib = 'Center';   
end

pdOffsetTib = min(tibPts(:,2)); % find most distal point

switch mlAxisTib
    case 'LatDist'
        latCompartmentIndicesTib = tibPts(:,3)<0;
        latCompartmentTib = tibPts(latCompartmentIndicesTib,:);
        [ ~ , latMostDistIndexTib] = min(latCompartmentTib(:,2));
        mlOffsetTib = latCompartmentTib(latMostDistIndexTib,3);
        %     mlOffsetTib = (max(tibPts(:,3)) - min(tibPts(:,3)))/4; %rotation about center of lateral condyle
    case 'LatEdge'
        latCompartmentIndicesTib = tibPts(:,3)>0;
        latCompartmentTib = tibPts(latCompartmentIndicesTib,:);
        [ ~ , latMostLatIndexTib] = max(latCompartmentTib(:,3));
        mlOffsetTib = latCompartmentTib(latMostLatIndexTib,3);
    case 'Center'
        mlOffsetTib = -(max(tibPts(:,3)) + min(tibPts(:,3)))/2;
    case 'MedDist'
        medCompartmentIndicesTib = tibPts(:,3)<0;
        medCompartmentTib = tibPts(medCompartmentIndicesTib,:);
        [ ~ , medMostDistIndexTib] = min(medCompartmentTib(:,2));
        mlOffsetTib = -medCompartmentTib(medMostDistIndexTib,3);
        %     mlOffsetTib = -(max(tibPts(:,3)) - min(tibPts(:,3)))/4; %rotation about center of medial condyle
    case 'MedEdge'
        medCompartmentIndicesTib = tibPts(:,3)<0;
        medCompartmentTib = tibPts(medCompartmentIndicesTib,:);
        [ ~ , medMostMedIndexTib] = min(medCompartmentTib(:,3));
        mlOffsetTib = medCompartmentTib(medMostMedIndexTib,3);
        %     mlOffsetTib = -(max(tibPts(:,3)) - min(tibPts(:,3)))/4; %rotation about center of medial condyle
end

apOffsetTib = -min(tibPts(:,1));

tibCentered = [tibPts(:,1)+apOffsetTib,tibPts(:,2)-pdOffsetTib,tibPts(:,3)-mlOffsetTib];

if plotComponents == 1 && plotTroubleshooting == 1
    plot3(tibCentered(:,1),tibCentered(:,2), tibCentered(:,3),'cx')
    axis equal tight
%     plot3(femCentered(:,1),femCentered(:,2), femCentered(:,3),'cx')
%     axis equal tight
end


%% ALIGN TO GS COORDINATE SYSTEM

%re-aligned fem with x = FE, y = VV, z = IE per Grood-Suntay Coordinate
%System

thetaY = deg2rad(-90);

Ry = [cos(thetaY)   0   sin(thetaY);
      0             1   0;
      -sin(thetaY)  0   cos(thetaY)];

femGS = femCentered * Ry;
femNormGS = femNorm * Ry;

tibGS = tibCentered * Ry;
tibNormGS = tibNorm * Ry;

thetaX = deg2rad(-90);
Rx = [  1   0               0;
        0   cos(thetaX)     -sin(thetaX);
        0   sin(thetaX)     cos(thetaX)];

femGS = femGS * Rx;
femNormGS = femNormGS * Rx;

tibGS = tibGS * Rx;
tibNormGS = tibNormGS * Rx;

if plotComponents == 1 && plotTroubleshooting == 1
%     plot3(femGS(:,1),femGS(:,2), femGS(:,3),'ms')
    plot3(tibGS(:,1),tibGS(:,2), tibGS(:,3),'ms')
    axis equal tight
end


%% TRANSFORM ARTICULAR SURFACES USING RANDOMLY GENERATED TRANSFORMATION MATRIX
%X = A-P, Y = P-D, Z = M-L
%
%Rx = [ 1   0   0   0;
%       0   c   -s  0;
%       0   s   c   0;
%       0   0   0   1];
%
%Ry = [ c   0   s   0;
%       0   1   0   0;
%       -s  0   c   0;
%       0   0   0   1];
%
%Rz = [ c   -s  0   0;
%       s   c   0   0;
%       0   0   1   0;
%       0   0   0   1];


feFem   = deg2rad(0); %flexion-extension rotation in rad
vvFem   = deg2rad(vvFemMalrot(iSim)); %varus-valgus rotation in rad +Valgus -Varus
ieFem   = deg2rad(ieFemMalrot(iSim)); %internal-external rotation in rad +External -Internal

%transformation applied using Cardan sequence fe-vv-ie
transFem = [cos(vvFem)*cos(ieFem)                                       -cos(vvFem)*sin(ieFem)                                  sin(vvFem)              ;
            cos(feFem)*sin(ieFem) + cos(ieFem)*sin(feFem)*sin(vvFem)    cos(feFem)*cos(ieFem)-sin(feFem)*sin(vvFem)*sin(ieFem)  -cos(vvFem)*sin(feFem)  ;
            sin(feFem)*sin(ieFem) - cos(feFem)*cos(ieFem)*sin(vvFem)    cos(ieFem)*sin(feFem)+cos(feFem)*sin(vvFem)*sin(ieFem)  cos(feFem)*cos(vvFem)   ];
        
% femCenteredTrans = femCentered * transFem;
femCenteredTrans = femGS * transFem;
femNormTrans = femNormGS * transFem;

feTib   = deg2rad(0); %flexion-extension rotation in rad
vvTib   = deg2rad(vvTibMalrot(iSim)); %varus-valgus rotation in rad
ieTib   = deg2rad(ieTibMalrot(iSim)); %internal-external rotation in rad

%transformation applied using Cardan sequence fe-vv-ie
transTib = [cos(vvTib)*cos(ieTib)                                       -cos(vvTib)*sin(ieTib)                                  sin(vvTib)              ;
            cos(feTib)*sin(ieTib) + cos(ieTib)*sin(feTib)*sin(vvTib)    cos(feTib)*cos(ieTib)-sin(feTib)*sin(vvTib)*sin(ieTib)  -cos(vvTib)*sin(feTib)  ;
            sin(feTib)*sin(ieTib) - cos(feTib)*cos(ieTib)*sin(vvTib)    cos(ieTib)*sin(feTib)+cos(feTib)*sin(vvTib)*sin(ieTib)  cos(feTib)*cos(vvTib)   ];
        
tibCenteredTrans = tibGS * transTib;
tibNormTrans = tibNormGS * transTib;

%% ALIGN BACK TO IMPORTED COORDINATE SYSTEM

% Re-aligned fem with x = FE, y = VV, z = IE per Grood-Suntay Coordinate System

thetaX = deg2rad(90);
Rx = [  1   0               0;
        0   cos(thetaX)     -sin(thetaX);
        0   sin(thetaX)     cos(thetaX)];

femCenteredTrans = femCenteredTrans * Rx;
femNormTrans = femNormTrans * Rx;

tibCenteredTrans = tibCenteredTrans * Rx;
tibNormTrans = tibNormTrans * Rx;

thetaY = deg2rad(90);

Ry = [cos(thetaY)   0   sin(thetaY);
      0             1   0;
      -sin(thetaY)  0   cos(thetaY)];

% thetaZ = deg2rad(-90);
% 
% Rz = [ cos(thetaZ)   -sin(thetaZ)  0;
%       sin(thetaZ)   cos(thetaZ)   0;
%       0             0             1];

femCenteredTrans = femCenteredTrans * Ry;
femNormTrans = femNormTrans * Ry;

tibCenteredTrans = tibCenteredTrans * Ry;
tibNormTrans = tibNormTrans * Ry;

if plotComponents == 1 && plotTroubleshooting == 1
    plot3(femCenteredTrans(:,1), femCenteredTrans(:,2), femCenteredTrans(:,3),'b.')
    plot3(tibCenteredTrans(:,1), tibCenteredTrans(:,2), tibCenteredTrans(:,3),'b.')
    axis equal tight
end

%% REMOVE OFFSET FOR TRANSFORMATION TO BRING BACK TO ORIGINAL MODEL COORDINATE SYSTEM

femTrans = [femCenteredTrans(:,1) - apOffsetFem,femCenteredTrans(:,2) + pdOffsetFem,femCenteredTrans(:,3) - mlOffsetFem];
tibTrans = [tibCenteredTrans(:,1) - apOffsetTib,tibCenteredTrans(:,2) + pdOffsetTib,tibCenteredTrans(:,3) + mlOffsetTib];

if plotComponents == 1
    figure()
%     plot3(femTrans(:,1), femTrans(:,2), femTrans(:,3), 'k.')
    hold on
    plot3(tibTrans(:,1), tibTrans(:,2), tibTrans(:,3), 'k.')
    axis equal tight
    
%     plot3(femPts(:,1),femPts(:,2), femPts(:,3),'g.')
%     hold on
    plot3(tibPts(:,1),tibPts(:,2), tibPts(:,3),'g.')

end

%assemble output structures
%--------------------------
femMalrot.pts = femTrans;
femMalrot.norm = femNormTrans;
femMalrot.tri = femTri;

tibMalrot.pts = tibTrans;
tibMalrot.norm = tibNormTrans;
tibMalrot.tri = tibTri;

%% WRITE TRANSFORMED ART SURFACE FILES

newFemName = [ femName(1:end-4) , '_' num2str( iSim ) , '.stl' ] ;
newTibName = [ tibName(1:end-4) , '_' num2str( iSim ) , '.stl' ] ;

switch simulationProg
    case 'OpenSim'

        femTR = triangulation( femMalrot.tri , femMalrot.pts ) ;
        tibTR = triangulation( tibMalrot.tri , tibMalrot.pts ) ;
        
        stlwrite( femTR , fullfile( filePathOut , 'Geometry' , newFemName ) ) ;
        stlwrite( tibTR , fullfile( filePathOut , 'Geometry' , newTibName ) ) ;
        
end

done = iSim;