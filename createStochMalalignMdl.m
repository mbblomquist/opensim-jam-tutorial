function doneMessage = createStochMalalignMdl( Params )
%
% -------------------------------------------------------------------------
% this script is used to
% run monte carlo simulations on vv alignment of femoral component,
% ie alignment of femoral component, and vv alignment of tibial component
%
% Overview:
%
% 1. Import surface model
% 2. Define stochastic alignment parameters
% 3. Create malaligned component files
%
%
% WRITTEN BY JOSH ROTH
% v1 (JDR): Initial Release
% v2 (MBB): Compatible with OpenSim
% v3 (MBB): Updated to fit into OpenSim_JAM pipeline
% -------------------------------------------------------------------------
%% Parse inputs
%==============

filePathOut = Params.fileOut ; % file to save stoch implant files
filePathIn = Params.fileIn ; % implant file locations
femName = Params.femImplant ; % femur implant name
tibName = Params.tibImplant ; % tibia implant name
distType = Params.distType ; % distribution type
numMdls = Params.numModels ; % number of models to create
femRot = Params.femRot ; % distribution parameters for femur
tibRot = Params.tibRot ; % distribution parameters for tibia

%% Load articular surface models
%===============================
% Simulations can either be run using OpenSim or SIMM
% Each program requires different files type
%
% OpenSim requires OBJ files
% SIMM requires ASC files

if ~exist( filePathOut, 'dir' )
    mkdir( filePathOut ) ;
    copyfile( filePathIn , fullfile( filePathOut , 'Geometry' ) )
end

% Check extension name of femur part
[ ~, ~, fileExt ] = fileparts( femName );

switch fileExt
    case '.stl'
        simulationProg = 'OpenSim';
        % Compute points, triangulation, and normals
        [ Fem.pts, Fem.tri, Fem.norm ] = STL_Import( fullfile( filePathIn, femName ), 1 ) ;
        [ Tib.pts, Tib.tri, Tib.norm ] = STL_Import( fullfile( filePathIn, tibName ), 1 );
    case '.asc'
        simulationProg = 'SIMM';
        % load asc's into structure to pass to malrotation function
        [ Fem.pts, Fem.norm, Fem.tri ] = load_asc( fullfile( filePathIn, femName ) );
        [ Tib.pts, Tib.norm, Tib.tri ] = load_asc( fullfile( filePathIn, tibName ) );
    otherwise
        disp('File type not supported');
        return
end

% Insert into Params structure for malrotateArticularSurfaces
Params.Fem = Fem ;
Params.Tib = Tib ;

%% Define stochastic alignment parameters
%==============================

%!!!TODO!!! add all DOFs

switch distType
    case 'normal'
        % Femural component
        % Fem vv
        if isfield( femRot , 'vv' )
            mu = femRot.vv( 1 );
            sd = femRot.vv( 2 );
            vvFemMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            vvFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Fem ie
        if isfield( femRot , 'ie' )
            mu = femRot.ie( 1 ) ;
            sd = femRot.ie( 2 ) ;
            ieFemMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            ieFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Tibial component
        % Tib vv
        if isfield( tibRot , 'vv' )
            mu = tibRot.vv( 1 );
            sd = tibRot.vv( 2 );
            vvTibMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            vvTibMalrot = zeros( numMdls, 1 ) ;
        end
        % Tib ie
        if isfield( tibRot , 'ie' )
            mu = tibRot.ie( 1 );
            sd = tibRot.ie( 2 );
            ieTibMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            ieTibMalrot = zeros( numMdls, 1 ) ;
        end
    case 'uniform'
        % Femoral component
        % Fem vv
        if isfield( femRot , 'vv' )
            lowLim = femRot.vv( 1 );
            upLim = femRot.vv( 2 );
            vvFemMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            vvFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Fem ie
        if isfield( femRot , 'ie' )
            lowLim = femRot.ie( 1 );
            upLim = femRot.ie( 2 );
            ieFemMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            ieFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Tibial component
        % Tib vv
        if isfield( tibRot , 'vv' )
            lowLim = tibRot.vv( 1 );
            upLim = tibRot.vv( 2 );
            vvTibMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            vvTibMalrot = zeros( numMdls, 1 ) ;
        end
        % Tib ie
        if isfield( tibRot , 'ie' )
            lowLim = tibRot.ie( 1 );
            upLim = tibRot.ie( 2 );
            ieTibMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            ieTibMalrot = zeros( numMdls, 1 ) ;
        end
end

%% Create malaligned component files
%===================================

%Save MAT-file containing all simulation parameters
%--------------------------------------------------
vvFem = round( vvFemMalrot , 2 ) ;
ieFem = round( ieFemMalrot , 2 ) ;
vvTib = round( vvTibMalrot , 2 ) ;
ieTib = round( ieTibMalrot , 2 ) ;

% Store into Params struct for malrotateArticularSurfaces.m
Params.simulationProg = simulationProg ;
Params.vvFemMalrot = vvFemMalrot ;
Params.ieFemMalrot = ieFemMalrot ;
Params.vvTibMalrot = vvTibMalrot ;
Params.ieTibMalrot = ieTibMalrot ;

mdlNums = 1 : numMdls ;

malrotLog.hdr = { 'Sim Num', 'vvFem', 'ieFem', 'vvTib', 'ieTib' };
malrotLog.data = [ mdlNums' , vvFem , ieFem , vvTib , ieTib ];

save( [ filePathOut , '\malAlignImplantData.mat' ] , 'malrotLog' )

%Generate stochastic malrotation models
%--------------------------------------
for iMdl = 1 : numMdls
    [ ~ , ~ , ~ ] = malrotateArticularSurfaces( Params , iMdl ) ;
end

doneMessage = 'Stoch articular surfaces created!';
