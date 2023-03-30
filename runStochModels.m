%% Run Stoch Models
% =========================================================================
% Author: Matthew Blomquist
%
% Revised code from Colin Smith:
% https://github.com/opensim-jam-org/jam-resources/tree/main/matlab
%
% Purpose: To create and run models with stochastic parameters
%
% Output: Creates 1) executables to run forward simulations and 
%   2) stochastic models with altered parameters
%
% Other .m files required:
%   createStochModels.m
%   createStochMalalignMdl.m
%
% Revision history:
%   v1      11-27-2022      First commit (MBB)
%   v2      03-28-2023      Update code for local or HT (MBB)
%
%==========================================================================
clear ; close all ; clc ;

import org.opensim.modeling.*
% Logger.setLevelString( 'Info' ) ;

%% ======================= Specify Settings ===========================
% Look through this entire section to change parameters to whatever you
% wish to run. You shouldn't have to change anything other than this
% section. This entire section is creating a structure full of parameters
% to create the models and the files to run simulations
% =====================================================================

% ------------------------------------------------------------------------
% -------------------------- SPECIFY LOCAL VS HT -------------------------
% ------------------------------------------------------------------------

% Specify whether to run locally or whether you will run the models on the
% high-throughput grid
% Options: 'local' or 'HT'
Params.localOrHT = 'HT' ;

% Set base name of output folder where models, exectuables, and inputs
% should be created
switch Params.localOrHT
    case 'local'
        % You don't need to change this. The output directory is in the
        % current directory (pwd = print working directory)
        Params.baseOutDir = pwd ;
    case 'HT'
        % IF RUNNING ON CHTC, CHANGE THIS TO PLACE YOU WANT FILES TO GO
        % I would do it outside of this folder because it's too many files
        % for git to track. I usually create them in a folder on my
        % desktop, but another folder in documents works as well
        Params.baseOutDir = 'C:\Users\mbb201\Desktop\htcTKArelease\testNewCode2' ;
        % Also specify which study ID for BAM lab work (not too important,
        % but this is what some files will have for a prefix in their name)
        Params.studyId = 'bam014' ;
end

% ------------------------------------------------------------------------
% ----------------------- SPECIFY MODEL PARAMETERS -----------------------
% ------------------------------------------------------------------------

% Number of models to create and run
Params.numModels = 10 ;

% Base model to use. Options are in lenhart2015 folder
%   Current options =
%       'lenhart2015' (intact model)
%       'lenhart2015_implant' (TKA model - implants and no ACL or MCLd)
Params.baseMdl = 'lenhart2015' ;

% Names of ligaments to change
%   Options: 'allLigs' to change all the ligaments in the model
%                           OR
%     [cell array with each ligament you want to change]
%       'MCLd' , 'MCLs', 'MCLp', 'ACLpl' , 'ACLam' , 'LCL', 'ITB', 'PFL',
%       'pCAP', 'PCLpm', 'PCLal', 'PT', 'lPFL', 'mPFL'
Params.ligNamesToChange = 'allLigs' ;

% Ligament properties to change [cell array]
%   Options: 'linear_stiffness', 'slack_length'
Params.ligPropsToChange = { 'linear_stiffness' , 'slack_length' } ;

% Probability distribution type [cell array for each ligPropsToChange]
%   Options: 'normal' , 'uniform'
Params.probDistType = { 'normal' , 'normal' } ;

% Probability distribution reference [cell array for each ligPropsToChange]
%   Options: 'relativePercent' , 'relativeAbs' , 'absolute'
Params.probDistRef = { 'relativePercent' , 'relativePercent' } ;

% Probability distribution parameters (in percent change from baseline model)
%  [cell array for each ligPropsToChange]
%   For 'normal': [ <mean> , <std> ]
%       Example: [ 0, 0.3 ] = mean of 0% change (same as baseline model)
%       with a standard deviation of 30% from the baseline model
%   For 'uniform': [ <lower_limit> , <upper_limit> ]
%       Example: [ -0.2, 0.2 ] = limits of distribution are -20 to 20% of
%       the baseline model value
Params.probDistParams = { [ 0, 0 ] , [ 0 , 0 ] } ;

% ------------------------------------------------------------------------
% --------------------- SPECIFY IMPLANT PARAMETERS -----------------------
% ------------------------------------------------------------------------

% Only need to specify implant parameters if the model is
% lenhart2015_implant
switch Params.baseMdl
    case 'lenhart2015_implant'

        % Femur and tibia implant names
        Params.femImplant = 'femur_component_surface_gc_katka-lenhart_updated.stl' ;
        Params.tibImplant = 'tibial_insert_surface_gc_katka-lenhart_updated.stl' ;

        % Directory with implant files
        Params.implantDir = fullfile( pwd , 'lenhart2015' , 'Geometry' ) ;

        % Probability distribution type for the implants
        %   Options: 'normal' , 'uniform'
        Params.distType = 'uniform' ;

        % Specify lower and upper limits (if uniform) or mean and std (if
        % normal) for femur and tibia Varus-valgus and internal-external
        % rotation of implants. Comment out the sections if you don't want
        % to malalign in that degree of freedom. NOTE: All values in deg
        Params.femRot.vv = [ -2 , 2 ] ;
        % Params.femRot.ie = [ -2 , 0 ] ;
        Params.tibRot.vv = [ -2 , 2 ] ;
        % Params.tibRot.ie = [ 0 2 ] ;
end

% ------------------------------------------------------------------------
% -------------------- SPECIFY SIMULATION PARAMETERS ---------------------
% ------------------------------------------------------------------------

% Specify forward simulation test(s) to run [cell array]
%   For Laxity Tests, options are:
%       Anterior: 'ant'
%       Posterior: 'post'
%       Varus: 'var'
%       Valgus: 'val'
%       Internal Rotation: 'ir'
%       External Rotation: 'er'
%       Compression: 'comp'
%       Distraction: 'dist'
%   For Passive Flexion, options are:
%       Passive: 'flex'
Params.testDOFs = { 'var' , 'val' } ;

% Specify flexion angle(s) of knee during each simulation [cell array]
%   For passive flexion, you can leave blank
%   Each testDOF will be run at each flexion angle (so the total number of
%   simulations will be length(testDOFs) * length( kneeFlexAngles )
Params.kneeFlexAngles = { 0 , 20 } ;

% Specify external load(s) applied, one for each testDOFs [cell array]
%   Put 0 if passive flexion test
%   Keep this number positive
Params.externalLoads = { 10 , 10 } ;

%% ============ Checks to make sure Params is set up correctly ============
% Throw an error before running code if something in Params is not set up
% correctly
% =========================================================================

if length( Params.ligPropsToChange ) ~= length( Params.probDistType )
    error( 'Error: ligPropsToChange needs to be the same length as probDistType' )
elseif length( Params.ligPropsToChange ) ~= length( Params.probDistRef )
    error( 'Error: ligPropsToChange needs to be the same length as probDistRef' )
elseif length( Params.ligPropsToChange ) ~= length( Params.probDistParams )
    error( 'Error: ligPropsToChange needs to be the same length as probDistParams' )
elseif length( Params.testDOFs ) ~= length( Params.externalLoads )
    error( 'Error: testDOFs needs to be the same length as externalLoads' )
end

%% ======================== Compute Trial Name ===========================
% Computes the trial name(s) that will be used for creating executables and
% input files. Puts it in a cell array
% ========================================================================

% For passive flexion:
%   1) flexion, 2) passive, 3) start flexion, 4) end flexion
%   Ex: 'flex_passive_0_90'
% For laxity tests:
%   1) laxity, 2) degree of freedom, 3) force applied, 4) flexion angle
%   Ex: 'lax_var_frc10_25'

Params.trialNames = { } ; % initialize
trialCounter = 1 ; % counter for loop
for iDOF = 1 : length( Params.testDOFs )
    if ~isequal( Params.testDOFs{iDOF} , 'flex' ) % if laxity test
        for iAng = 1 : length( Params.kneeFlexAngles ) % loop through flexion angles
            Params.trialNames{ trialCounter } = ...
                [ 'lax_' , Params.testDOFs{iDOF} , '_frc' , num2str( Params.externalLoads{iDOF} ) , '_' , num2str( Params.kneeFlexAngles{iAng} ) ] ;
            trialCounter = trialCounter + 1 ;
        end
    elseif isequal( Params.testDOFs{iDOF} , 'flex' ) % if passive flexion test
        Params.trialNames{ trialCounter } = 'flex_passive_0_90' ;
        trialCounter = trialCounter + 1 ;
    end
end

Params.numTrials = length( Params.trialNames ) ;

clear trialCounter iDOF iAng

%% ======================== Create Stoch Models ==========================
% Sets up folders to put models and executables, then create stochastic
% models (and implants if applicable) and put them in desired output folder
% ========================================================================
Params.baseMdlFile = fullfile( pwd , 'lenhart2015' , [ Params.baseMdl , '.osim' ] ) ;

switch Params.localOrHT
    case 'local'
        % Check if inputs, stochModels, and results directories exist. If
        % they don't, create them. If they do, delete contents to start
        % from scratch

        % inputs
        if ~exist( 'inputs' , 'dir' )
            mkdir( 'inputs' )
        else
            xmlFiles = dir( fullfile( 'inputs' , '*.xml' ) ) ;
            for iFile = 1 : length( xmlFiles )
                delete( fullfile( 'inputs' , xmlFiles(iFile).name ) )
            end
            stoFiles = dir( fullfile( 'inputs' , '*.sto' ) ) ;
            for iFile = 1 : length( stoFiles )
                delete( fullfile( 'inputs' , stoFiles(iFile).name ) )
            end
        end

        % stochModels
        if ~exist( 'stochModels' , 'dir' )
            mkdir( 'stochModels' )
            copyfile( 'lenhart2015\Geometry' , 'stochModels\Geometry' )
        else
            osimFiles = dir( fullfile( 'stochModels' , '*.osim' ) ) ;
            for iFile = 1 : length( osimFiles )
                delete( fullfile( 'stochModels' , osimFiles(iFile).name ) )
            end
            matFiles = dir( fullfile( 'stochModels' , '*.mat' ) ) ;
            for iFile = 1 : length( matFiles )
                delete( fullfile( 'stochModels' , matFiles(iFile).name ) )
            end
        end

        % results
        if ~exist( 'results' , 'dir' )
            mkdir( 'results' )
        else
            stoFiles = dir( fullfile( 'results' , '*.sto' ) ) ;
            for iFile = 1 : length( stoFiles )
                delete( fullfile( 'results' , stoFiles(iFile).name ) )
            end
        end

    case 'HT'
        % Check to see if Params.baseOutDir exists. If it doesn't, then
        % create subfolders to be able to use on the high throughput grid
        if ~exist( Params.baseOutDir , 'dir' )
            mkdir( Params.baseOutDir )
        end
        if ~exist( fullfile( Params.baseOutDir , 'input' ) , 'dir' )
            mkdir( fullfile( Params.baseOutDir , 'input' ) )
            for iMdl = 1 : Params.numModels
                mkdir( fullfile( Params.baseOutDir , 'input' , num2str(iMdl-1) ) )
            end
        end
        for iDOF = 1 : length( Params.testDOFs )
            if ~exist( fullfile( Params.baseOutDir , Params.testDOFs{iDOF} ) , 'dir' )
                mkdir( fullfile( Params.baseOutDir , Params.testDOFs{iDOF} ) )
                copyfile( 'sharedFilesForHT\' , fullfile( Params.baseOutDir , Params.testDOFs{iDOF} , 'shared' ) )
                mkdir( fullfile( Params.baseOutDir , Params.testDOFs{iDOF} ) , 'results' )
                for iMdl = 1 : Params.numModels
                    mkdir( fullfile( Params.baseOutDir , Params.testDOFs{iDOF} , 'results' , num2str(iMdl-1) ) )
                end
            end
        end

end

%% ======================== Create Stoch Models ==========================
% Create stochastic models and (if applicable) malaligned implants
% ========================================================================

% Run createStochMalalignMdl code if running implant code
switch Params.baseMdl
    case 'lenhart2015_implant'
        doneMsg = createStochMalalignMdl( Params ) ;
end

% Create Stochastic models
StochMdlParams = createStochModels( Params ) ;

%% =================== Define Simulation Time Points =====================
% Define the duration of each portion of the simulation (flexion, settle,
% load, etc)
% ========================================================================

% Simulation consists of three to five phases:
% All simulations consist of:
%   Settle : allow knee to settle into equilbrium
%   Flex   : period of knee flexion
%   Settle : allow knee to settle into equilbrium
% If laxity test, add these two:
%   Force  : ramp up the desired external force
%   Settle : hold force constant and allow knee to settle into equilbrium

% Time step for simulation
timeStep = 0.01;

% Initialize
numPhases = zeros( Params.numTrials , 1 ) ;
time = cell( Params.numTrials , 1 ) ;
timePoints = cell( Params.numTrials , 1 ) ;
numSteps = zeros( Params.numTrials , 1 ) ;

% Loop through each test
for iTrial = 1 : Params.numTrials

    % Define temporary trial name to extract temporary DOF
    tempTrialName = Params.trialNames{ iTrial } ;
    splitName = split( tempTrialName , '_' ) ;
    if isequal( splitName{1} , 'flex' )
        tempDOF = 'flex' ;
    else
        tempDOF = splitName{2} ;
    end

    if strcmp( tempDOF , 'flex' ) % passive flexion
        % Specify duration of each phase
        phaseDurations = [ ...
            0.5 ... % Settle 1 Duration
            1.0 ... % Flex Duration
            0.5 ... % Settle 2 Duration
            ] ;
    else % laxity test
        % Specify duration of each phase
        phaseDurations = [ ...
            0.5 ... % Settle 1 Duration
            1.0 ... % Flex Duration
            0.5 ... % Settle 2 Duration
            1.0 ... % External Force Duration
            0.5 ... % Settle 3 Duration
            ] ;
    end

    % Number of phases
    numPhases(iTrial) = length( phaseDurations ) ;

    % Total duration of simulation
    simDuration = sum( phaseDurations ) ;

    % Time array of simulation
    time{iTrial} = 0 : timeStep : simDuration ;

    % Compute starting time points for each phase
    timePoints{iTrial} = zeros( 1 , numPhases(iTrial) + 1 ) ; % initialize array before for loop
    for iTimePt = 1 : numPhases(iTrial)
        timePoints{iTrial}( 1 , iTimePt+1 ) = sum( phaseDurations( 1 : iTimePt ) ) ;
    end

    % Compute the number of steps in the simulation
    numSteps(iTrial) = length( time{iTrial} ) ;

end

%% ================ Create Prescribed Coordinates File ===================
% The Prescribed Coordinates file sets the hip, knee, and pelvis flexion
% angles to their prescibed values
% ========================================================================

% For passive flexion:
%   1) flexion, 2) passive, 3) start flexion, 4) end flexion
%   Ex: 'flex_passive_0_90'
% For laxity tests:
%   1) laxity, 2) degree of freedom, 3) force applied, 4) flexion angle
%   Ex: 'lax_var_frc10_25'

% Hip flexion angle
hipFlexAngle = 0 ;

% Pelvis Tilt
pelvisTilt = 90 ; % 0 = standing, 90 = supine

% Loop through each test
for iTrial = 1 : Params.numTrials

    % Define temporary trial name to extract temporary DOF and angle
    tempTrialName = Params.trialNames{ iTrial } ;
    splitName = split( tempTrialName , '_' ) ;
    if isequal( splitName{1} , 'flex' )
        tempDOF = 'flex' ;
    else
        tempDOF = splitName{2} ;
    end
    tempAng = str2double( splitName{4} ) ;

    % .sto File Name
    prescribedCoordinatesFileName = ...
        [ 'prescribed_coordinates_' , tempTrialName , '.sto' ] ;

    coord_data.time = time{iTrial} ;

    % Knee data
    kneeFlexData = [ zeros( 1  , 2 ) , ones( 1 , numPhases(iTrial)-1 ) * tempAng ] ;
    smoothKneeFlexData = interp1( timePoints{iTrial} , kneeFlexData , time{iTrial} , 'pchip' ) ;
    coord_data.knee_flex_r = smoothKneeFlexData' ;

    % Hip data
    hipFlexData = [ zeros( 1  , 2 ) , ones( 1 , numPhases(iTrial)-1 ) * hipFlexAngle ] ;
    smoothHipFlexData = interp1( timePoints{iTrial} , hipFlexData , time{iTrial} , 'pchip' ) ;
    coord_data.hip_flex_r = smoothHipFlexData' ;

    % Pelvis data
    coord_data.pelvis_tilt = ones( length(time{iTrial}) , 1 ) * pelvisTilt ;

    % Function distributed in OpenSim Resources\Code\Matlab\Utilities
    coord_table = osimTableFromStruct( coord_data ) ;

    switch Params.localOrHT
        case 'local'
            outDir = fullfile( pwd , 'inputs' ) ;
        case 'HT'
            outDir = fullfile( Params.baseOutDir , tempDOF , 'shared' ) ;
    end

    % Write file
    STOFileAdapter.write( coord_table , fullfile( outDir , prescribedCoordinatesFileName ) ) ;

end


%% ==================== Create External Load Files =======================
% The external load files are to specify where, how much, and in which
% direction external loads should be applied during the simulation.
% ========================================================================

% Loop through each test
for iTrial = 1 : Params.numTrials

    % Define temporary trial name to extract temporary DOF
    tempTrialName = Params.trialNames{ iTrial } ;
    splitName = split( tempTrialName , '_' ) ;
    if isequal( splitName{1} , 'flex' )
        tempDOF = 'flex' ;
    else
        tempDOF = splitName{2} ;
    end

    % Create external loads files if it is a laxity test (i.e., not a
    %   passive flexion test)
    if ~strcmp( tempDOF , 'flex' ) % if not passive flexion test

        % External load magnitude
        tempLoad = str2double( splitName{3}(4:end) ) ;

        % .sto and .xml File Name
        externalLoadsSto = [ 'external_loads_' , tempTrialName , '.sto' ] ;
        externalLoadsXml = [ 'external_loads_' , tempTrialName , '.xml' ] ;

        % Define positive and negative directions
        switch tempDOF
            case { 'ant', 'ir', 'val', 'comp' }
                loadSign = 1 ;
            case { 'post', 'er', 'var', 'dist' }
                loadSign = -1 ;
        end

        % Define location and magnitude of load based on testDof and externalLoad
        switch tempDOF
            case { 'ant' , 'post' }
                loadPointHeight = -0.1 ; % Apply at the tibial tuberosity height similar to KT-1000 test
                loadMagnitude = tempLoad * loadSign ;
            case { 'var' , 'val' }
                loadPointHeight = -0.3 ; % Apply near ankle similiar to coronal laxity test
                loadMagnitude = tempLoad / abs(loadPointHeight) * loadSign ; % Moment, so account for moment arm
            case { 'ir' , 'er' , 'dist' , 'comp' }
                loadMagnitude = tempLoad * loadSign ; % Apply at location = 0
        end

        % Create external load array
        loadArray = [ zeros( 1 , 4 ) , loadMagnitude , loadMagnitude ] ;
        smoothLoadArray = interp1( timePoints{iTrial} , loadArray , time{iTrial} , 'pchip' );

        % Construct arrays for sto file based on
        loadData.time = time{iTrial} ;
        tempNumSteps = numSteps(iTrial) ;
        switch tempDOF
            case { 'ant', 'post' }
                % Applied load in x-direction distal to knee
                loadData.tibia_proximal_r_force_vx = smoothLoadArray' ;
                loadData.tibia_proximal_r_force_vy = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_vz = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_px = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_py = ones( tempNumSteps , 1 ) * loadPointHeight ;
                loadData.tibia_proximal_r_force_pz = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_x = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_y = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_z = zeros( tempNumSteps , 1 ) ;
            case { 'var', 'val' }
                % Applied load in z-direction at ankle
                loadData.tibia_proximal_r_force_vx = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_vy = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_vz = smoothLoadArray' ;
                loadData.tibia_proximal_r_force_px = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_py = ones( tempNumSteps , 1 ) * loadPointHeight ;
                loadData.tibia_proximal_r_force_pz = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_x = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_y = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_z = zeros( tempNumSteps , 1 ) ;
            case { 'ir', 'er' }
                % Applied torque about y-direction
                loadData.tibia_proximal_r_force_vx = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_vy = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_vz = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_px = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_py = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_pz = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_x = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_y = smoothLoadArray' ;
                loadData.tibia_proximal_r_torque_z = zeros( tempNumSteps , 1 ) ;
            case { 'comp', 'dist' }
                % Applied force about y-direction
                loadData.tibia_proximal_r_force_vx = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_vy = smoothLoadArray' ;
                loadData.tibia_proximal_r_force_vz = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_px = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_py = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_force_pz = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_x = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_y = zeros( tempNumSteps , 1 ) ;
                loadData.tibia_proximal_r_torque_z = zeros( tempNumSteps , 1 ) ;
        end

        % Function distributed in OpenSim Resources\Code\Matlab\Utilities
        loadTable = osimTableFromStruct( loadData ) ;
        loadTable.addTableMetaDataString( 'header' , [ tempTrialName , ' External Load' ] )

        switch Params.localOrHT
            case 'local'
                outDir = fullfile( pwd , 'inputs' ) ;
            case 'HT'
                outDir = fullfile( Params.baseOutDir , tempDOF , 'shared' ) ;
        end
        STOFileAdapter.write( loadTable , fullfile( outDir , externalLoadsSto ) );

        % Construct External Force
        extForce = ExternalForce() ;
        extForce.setName( [ tempTrialName , '_load' ] );
        extForce.set_applied_to_body( 'tibia_proximal_r' );
        extForce.set_force_expressed_in_body( 'tibia_proximal_r' );
        extForce.set_point_expressed_in_body( 'tibia_proximal_r' );
        extForce.set_force_identifier( 'tibia_proximal_r_force_v' );
        extForce.set_point_identifier( 'tibia_proximal_r_force_p' );
        extForce.set_torque_identifier( 'tibia_proximal_r_torque_' );

        % Construct External Loads
        extLoads = ExternalLoads() ;
        extLoads.setDataFileName( externalLoadsSto  ) ;
        extLoads.cloneAndAppend( extForce ) ;
        extLoads.print( fullfile( outDir , externalLoadsXml ) ) ;
    end
end
%% ====== Create Forsim Settings Files And Run/Create Linux Files ========
% Create the forward simulation settings file(s) that specifies how the
% forward simulation should be run, which files it calls on, and where the
% outputs should go.
% Then, for local simulations, run each model under the specified forward
% simulations. For HT simulations, create .sh and .sub files to run on the
% high throughput grid.
% ========================================================================

% Loop through each trial
for iTrial = 1 : Params.numTrials

    % Define temporary trial name to extract temporary DOF
    tempTrialName = Params.trialNames{ iTrial } ;
    splitName = split( tempTrialName , '_' ) ;
    if isequal( splitName{1} , 'flex' )
        tempDOF = 'flex' ;
    else
        tempDOF = splitName{2} ;
    end

    switch Params.localOrHT

        % ----------------------------------------------------------------
        % ---------------------------- LOCAL -----------------------------
        % ----------------------------------------------------------------

        case 'local'

            % Loop through each Model to create forsim_settings and run the model
            for iMdl = 1 : Params.numModels

                % If implant model, then switch out stoch implants so that
                % lenhart model reads the correct one
                switch Params.baseMdl
                    case 'lenhart2015_implant'
                        tempStochFemName = [ Params.femImplant(1:end-4) , '_' , num2str( iMdl ) , '.stl' ] ;
                        tempStochTibName = [ Params.tibImplant(1:end-4) , '_' , num2str( iMdl ) , '.stl' ] ;

                        movefile( fullfile( pwd , 'stochModels' , 'Geometry' , tempStochFemName ) , ...
                            fullfile( pwd , 'stochModels' , 'Geometry' , Params.femImplant ) )
                        movefile( fullfile( pwd , 'stochModels' , 'Geometry' , tempStochTibName ) , ...
                            fullfile( pwd , 'stochModels' , 'Geometry' , Params.tibImplant ) )
                end

                % Specify settings
                forsimSettingsFileName = [ 'forsim_settings_' , tempTrialName , '.xml' ] ;
                modelFile = [ './stochModels/lenhart2015_stoch' , num2str(iMdl) , '.osim' ] ;
                forsimResultDir = './results';
                resultsBasename = [ tempTrialName , '_' , num2str( iMdl ) ] ;

                % Set the integrator accuracy
                %   Choose value between speed (1e-2) vs accuracy (1e-8)
                integratorAccuracy = 1e-2 ;

                % Create ForsimTool
                forsim = ForsimTool() ;

                % Create AnalysisSet
                analysisSet = AnalysisSet() ;

                % Create ForceReporter
                frcReporter = ForceReporter();

                % Set settings
                frcReporter.setName( 'ForceReporter' ) ;
                analysisSet.cloneAndAppend( frcReporter ) ;
                forsim.set_AnalysisSet( analysisSet ) ;
                forsim.set_model_file( modelFile ) ; % location and name of model
                forsim.set_results_directory( forsimResultDir ) ; % location to put results files
                forsim.set_results_file_basename( resultsBasename ) ; % basename of results files
                forsim.set_start_time( -1 ) ; % set to -1 to use data from input files
                forsim.set_stop_time( -1 ) ; % set to -1 to use data from input files
                forsim.set_integrator_accuracy( integratorAccuracy ) ; % accuracy of the solver
                forsim.set_constant_muscle_control( 0.001 ) ; % 0.001 to represent passive state
                forsim.set_use_activation_dynamics( true ) ; % use activation dynamics
                forsim.set_use_tendon_compliance( true ) ; % use tendon compliance
                forsim.set_use_muscle_physiology( true ) ; % use muscle physiology
                % Set all coordinates to be unconstrained except flexion-extension
                forsim.set_unconstrained_coordinates( 0 , '/jointset/knee_r/knee_add_r' ) ;
                forsim.set_unconstrained_coordinates( 1 , '/jointset/knee_r/knee_rot_r' ) ;
                forsim.set_unconstrained_coordinates( 2 , '/jointset/knee_r/knee_tx_r' ) ;
                forsim.set_unconstrained_coordinates( 3 , '/jointset/knee_r/knee_ty_r' ) ;
                forsim.set_unconstrained_coordinates( 4 , '/jointset/knee_r/knee_tz_r' ) ;
                forsim.set_unconstrained_coordinates( 5 , '/jointset/pf_r/pf_flex_r' ) ;
                forsim.set_unconstrained_coordinates( 6 , '/jointset/pf_r/pf_rot_r' ) ;
                forsim.set_unconstrained_coordinates( 7 , '/jointset/pf_r/pf_tilt_r' ) ;
                forsim.set_unconstrained_coordinates( 8 , '/jointset/pf_r/pf_tx_r' ) ;
                forsim.set_unconstrained_coordinates( 9 , '/jointset/pf_r/pf_ty_r' ) ;
                forsim.set_unconstrained_coordinates( 10 , '/jointset/pf_r/pf_tz_r' ) ;
                forsim.set_prescribed_coordinates_file( [ 'prescribed_coordinates_' , tempTrialName , '.sto' ] ) ;
                if strcmp( tempDOF , 'flex' ) % passive flexion
                    forsim.set_external_loads_file( '' ) ;
                else % laxity test
                    forsim.set_external_loads_file( [ 'external_loads_' , tempTrialName , '.xml' ] ) ;
                end
                forsim.set_use_visualizer( false ) ; % use visualizer while running (true or false)
                forsim.print( fullfile( pwd , 'inputs' , forsimSettingsFileName ) ) ;

                tic % compute time that simulation runs
                disp( [ 'Running Forsim Tool, Model '  , num2str( iMdl ) ] )
                forsim.run() ;
                toc

                % If implant model, then switch back names of implants
                switch Params.baseMdl
                    case 'lenhart2015_implant'
                        movefile( fullfile( pwd , 'stochModels' , 'Geometry' , Params.femImplant ) , ...
                            fullfile( pwd , 'stochModels' , 'Geometry' , tempStochFemName ) )
                        movefile( fullfile( pwd , 'stochModels' , 'Geometry' , Params.tibImplant ) , ...
                            fullfile( pwd , 'stochModels' , 'Geometry' , tempStochTibName ) )
                end

            end

            % ----------------------------------------------------------------
            % ------------------------------ HT ------------------------------
            % ----------------------------------------------------------------

        case 'HT'

            % Specify settings
            forsimSettingsFileName = [ 'forsim_settings_' , tempTrialName , '.xml' ] ;
            modelFile = './lenhart2015_stoch.osim' ;
            forsimResultDir = './' ;
            resultsBasename = tempTrialName ;

            % Set the integrator accuracy
            %   Choose value between speed (1e-2) vs accuracy (1e-8)
            integratorAccuracy = 1e-2 ;

            % Create ForsimTool
            forsim = ForsimTool() ;

            % Create AnalysisSet
            analysisSet = AnalysisSet() ;

            % Create ForceReporter
            frcReporter = ForceReporter();

            % Set settings
            frcReporter.setName( 'ForceReporter' ) ;
            analysisSet.cloneAndAppend( frcReporter ) ;
            forsim.set_AnalysisSet( analysisSet ) ;
            forsim.set_model_file( modelFile ) ; % location and name of model
            forsim.set_results_directory( forsimResultDir ) ; % location to put results files
            forsim.set_results_file_basename( resultsBasename ) ; % basename of results files
            forsim.set_start_time( -1 ) ; % set to -1 to use data from input files
            forsim.set_stop_time( -1 ) ; % set to -1 to use data from input files
            forsim.set_integrator_accuracy( integratorAccuracy ) ; % accuracy of the solver
            forsim.set_constant_muscle_control( 0.001 ) ; % 0.001 to represent passive state
            forsim.set_use_activation_dynamics( true ) ; % use activation dynamics
            forsim.set_use_tendon_compliance( true ) ; % use tendon compliance
            forsim.set_use_muscle_physiology( true ) ; % use muscle physiology
            % Set all coordinates to be unconstrained except flexion-extension
            forsim.set_unconstrained_coordinates( 0 , '/jointset/knee_r/knee_add_r' ) ;
            forsim.set_unconstrained_coordinates( 1 , '/jointset/knee_r/knee_rot_r' ) ;
            forsim.set_unconstrained_coordinates( 2 , '/jointset/knee_r/knee_tx_r' ) ;
            forsim.set_unconstrained_coordinates( 3 , '/jointset/knee_r/knee_ty_r' ) ;
            forsim.set_unconstrained_coordinates( 4 , '/jointset/knee_r/knee_tz_r' ) ;
            forsim.set_unconstrained_coordinates( 5 , '/jointset/pf_r/pf_flex_r' ) ;
            forsim.set_unconstrained_coordinates( 6 , '/jointset/pf_r/pf_rot_r' ) ;
            forsim.set_unconstrained_coordinates( 7 , '/jointset/pf_r/pf_tilt_r' ) ;
            forsim.set_unconstrained_coordinates( 8 , '/jointset/pf_r/pf_tx_r' ) ;
            forsim.set_unconstrained_coordinates( 9 , '/jointset/pf_r/pf_ty_r' ) ;
            forsim.set_unconstrained_coordinates( 10 , '/jointset/pf_r/pf_tz_r' ) ;
            forsim.set_prescribed_coordinates_file( [ 'prescribed_coordinates_' , tempTrialName , '.sto' ] ) ;
            if strcmp( tempDOF , 'flex' ) % passive flexion
                forsim.set_external_loads_file( '' ) ;
            else % laxity test
                forsim.set_external_loads_file( [ 'external_loads_' , tempTrialName , '.xml' ] ) ;
            end
            forsim.set_use_visualizer( false ) ; % use visualizer while running (true or false)
            forsim.print( fullfile( Params.baseOutDir , tempDOF , 'shared' , forsimSettingsFileName ) ) ;

    end

end

%% ========== Create Sh and Sub Files if Running on the Grid =============
% Running on linux requires a couple additional files to run called a
% shared file (.sh) and a submit file (.sub). So, if running on the grid,
% create these extra files
% ========================================================================

switch Params.localOrHT
    case 'HT'
        for iDOF = 1 : length( Params.testDOFs )

            Params.tempTestDOF = Params.testDOFs{iDOF} ;

            % ----------------
            % Tar shared files
            % ----------------
            Params.sharedTarName = 'shared.tar.gz' ;
            Params.sharedDir = fullfile( Params.baseOutDir , Params.tempTestDOF , 'shared' ) ;

            disp( 'Tarring files...' )
            tar( fullfile( Params.sharedDir , Params.sharedTarName ) , '.' , Params.sharedDir )
            disp( [ 'Tar file created in ' , Params.sharedDir ] )

            % Delete rest of contents in shared folder
            sharedDirContents = dir( Params.sharedDir ) ;
            for iFile = 1 : length( sharedDirContents )
                if ~isequal( sharedDirContents(iFile).name , '.' ) && ...
                        ~isequal( sharedDirContents(iFile).name , '..' ) && ...
                        ~isequal( sharedDirContents(iFile).name , Params.sharedTarName )
                    % If not the tar file, the current dir, or upper dir,
                    % then delete the file
                    if sharedDirContents(iFile).isdir
                        rmdir( fullfile( Params.sharedDir , sharedDirContents(iFile).name ) , 's' )
                    else
                        delete( fullfile( Params.sharedDir , sharedDirContents(iFile).name ) )
                    end
                end
            end

            % ----------------
            % Create .sub file
            % ----------------
            Params.date = char( datetime( 'today', 'format', 'yyyy-MM-dd' ) ) ;
            Params.trialDir = fullfile( Params.baseOutDir , Params.tempTestDOF ) ;
            Params.shName = 'runJob.sh' ;
            outMsg = writeHtcSubFile( Params ) ;
            disp( outMsg )

            % ---------------
            % Create .sh file
            % ---------------
            outMsg = writeHtcShFile( Params ) ;
            disp( outMsg )

        end
end