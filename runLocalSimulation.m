%% Example Anterior Laxity
%==========================================================================
% Author: Colin Smith
%   https://github.com/opensim-jam-org/jam-resources/tree/main/matlab
%
% Revised by: Matthew Blomquist
%==========================================================================
clear ; close all ; clc ;

import org.opensim.modeling.*
% Logger.setLevelString( 'Info' ) ;

%% Set Simulation Parameters

% Specify starting directory and location of model file
startDir = 'C:\Users\mbb201\Documents\MATLAB\Research\opensimModeling\opensim-jam-code' ;

% Laxity Test or Passive Flexion Test
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
testDof = 'flex' ;

% Specify flexion angle of knee during simulation
%   For passive flexion, put angle to flex to (e.g., 90 for 90deg)
kneeFlexAngle = 90 ;

% Specify external load applied
%   Put 0 if passive flexion test
%   Keep this number positive
externalLoad = 0 ;

%% Compute trial name

if strcmp( testDof , 'flex' ) % passive flexion
    trialName = [ 'flex_passive_0_' , num2str( kneeFlexAngle ) ] ;
else % laxity test
    trialName = [ 'lax_' , testDof , '_frc' , num2str( externalLoad ) , '_' , num2str( kneeFlexAngle ) ] ;
end

%% Define Simulation Time Points
%-------------------------------
% Simulation consists of three to five phases:
% All simulations:
%   Settle : allow knee to settle into equilbrium
%   Flex   : hip and knee flexion
%   Settle : allow knee to settle into equilbrium
% If laxity test, add these two:
%   Force  : ramp up the desired external force
%   Settle : hold force constant and allow knee to settle into equilbrium

% Time step for simulation
timeStep = 0.01;

if strcmp( testDof , 'flex' ) % passive flexion
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
numPhases = length( phaseDurations ) ;

% Total duration of simulation
simDuration = sum( phaseDurations ) ;

% Time array of simulation
time = 0 : timeStep : simDuration ;

% Compute starting time points for each phase
timePoints = zeros( 1 , numPhases + 1 ) ; % initialize array before for loop
for iTimePt = 1 : numPhases
    timePoints( 1 , iTimePt+1 ) = sum( phaseDurations( 1 : iTimePt ) ) ;
end

% Compute the number of steps in the simulation
numSteps = length( time ) ;

%% Create Prescribed Coordinates File
%------------------------------------

% .sto File Name
prescribedCoordinatesFileName = [ 'prescribed_coordinates_' , trialName , '.sto' ] ;

coord_data.time = time ;

% Knee data
kneeFlexData = [ zeros( 1  , 2 ) , ones( 1 , numPhases-1 ) * kneeFlexAngle ] ;
smoothKneeFlexData = interp1( timePoints , kneeFlexData , time , 'pchip' ) ;
coord_data.knee_flex_r = smoothKneeFlexData' ;

% Hip data
hipFlexAngle = 0 ;
hipFlexData = [ zeros( 1  , 2 ) , ones( 1 , numPhases-1 ) * hipFlexAngle ] ;
smoothHipFlexData = interp1( timePoints , hipFlexData , time , 'pchip' ) ;
coord_data.hip_flex_r = smoothHipFlexData' ;

% Pelvis data
pelvisTilt = 90 ; % 0 = standing, 90 = supine
coord_data.pelvis_tilt = ones( length(time) , 1 ) * pelvisTilt ;

% Function distributed in OpenSim Resources\Code\Matlab\Utilities
coord_table = osimTableFromStruct( coord_data ) ;
STOFileAdapter.write( coord_table , fullfile( startDir , 'inputs' , prescribedCoordinatesFileName ) ) ;

%% Create External Loads Files
%------------------------------

if ~strcmp( testDof , 'flex' ) % if not passive flexion test

    % .sto and .xml File Name
    externalLoadsSto = [ 'external_loads_' , trialName , '.sto' ] ;
    externalLoadsXml = [ 'external_loads_' , trialName , '.xml' ] ;

    % Define positive and negative directions
    switch testDof
        case { 'ant', 'ir', 'val', 'comp' }
            loadSign = 1 ;
        case { 'post', 'er', 'var', 'dist' }
            loadSign = -1 ;
    end

    % Define location and magnitude of load based on testDof and externalLoad
    switch testDof
        case { 'ant' , 'post' }
            loadPointHeight = -0.1 ; % Apply at the tibial tuberosity height similar to KT-1000 test
            loadMagnitude = externalLoad * loadSign ;
        case { 'var' , 'val' }
            loadPointHeight = -0.3 ; % Apply near ankle similiar to coronal laxity test
            loadMagnitude = externalLoad / abs(loadPointHeight) * loadSign ; % Moment, so account for moment arm
        case { 'ir' , 'er' , 'dist' , 'comp' }
            loadMagnitude = externalLoad * loadSign ; % Apply at location = 0
    end

    % Create external load array
    loadArray = [ 0 , 0 , 0 , 0 , loadMagnitude , loadMagnitude ] ;
    smoothLoadArray = interp1( timePoints , loadArray , time , 'pchip' );

    % Construct arrays for sto file based on
    loadData.time = time ;
    switch testDof
        case { 'ant', 'post' }
            % Applied load in x-direction distal to knee
            loadData.tibia_proximal_r_force_vx = smoothLoadArray' ;
            loadData.tibia_proximal_r_force_vy = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_vz = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_px = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_py = ones(numSteps,1) * loadPointHeight ;
            loadData.tibia_proximal_r_force_pz = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_x = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_y = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_z = zeros( numSteps , 1 ) ;
        case { 'var', 'val' }
            % Applied load in z-direction at ankle
            loadData.tibia_proximal_r_force_vx = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_vy = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_vz = smoothLoadArray' ;
            loadData.tibia_proximal_r_force_px = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_py = ones(numSteps,1) * loadPointHeight ;
            loadData.tibia_proximal_r_force_pz = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_x = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_y = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_z = zeros( numSteps , 1 ) ;
        case { 'ir', 'er' }
            % Applied torque about y-direction
            loadData.tibia_proximal_r_force_vx = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_vy = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_vz = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_px = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_py = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_pz = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_x = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_y = smoothLoadArray' ;
            loadData.tibia_proximal_r_torque_z = zeros( numSteps , 1 ) ;
        case { 'comp', 'dist' }
            % Applied force about y-direction
            loadData.tibia_proximal_r_force_vx = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_vy = smoothLoadArray' ;
            loadData.tibia_proximal_r_force_vz = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_px = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_py = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_force_pz = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_x = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_y = zeros( numSteps , 1 ) ;
            loadData.tibia_proximal_r_torque_z = zeros( numSteps , 1 ) ;
    end

    % Function distributed in OpenSim Resources\Code\Matlab\Utilities
    loadTable = osimTableFromStruct( loadData ) ;

    loadTable.addTableMetaDataString( 'header' , [ trialName , ' External Load' ] )

    STOFileAdapter.write( loadTable , fullfile( startDir , 'inputs' , externalLoadsSto ) );

    % Construct External Force
    extForce = ExternalForce() ;
    extForce.setName( [ trialName , '_load' ] );
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
    extLoads.print( fullfile( startDir , 'inputs' , externalLoadsXml ) ) ;
end

%% Perform Simulation with ForsimTool

% Specify settings
forsimSettingsFileName = [ 'forsim_settings_' , trialName , '.xml' ] ;
modelFile = './lenhart2015/lenhart2015.osim';
forsimResultDir = './results';
resultsBasename = trialName ;
integratorAccuracy = 1e-2 ; % Speed (1e-2) vs accuracy (1e-8)

% Create ForsimTool
forsim = ForsimTool() ;

% Create AnalysisSet
analysisSet = AnalysisSet() ;

% Create Force Reporter
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
forsim.set_prescribed_coordinates_file( [ './inputs/' , prescribedCoordinatesFileName ] ) ;
if strcmp( testDof , 'flex' ) % passive flexion
    forsim.set_external_loads_file( '' ) ;
else % laxity test
    forsim.set_external_loads_file( [ './inputs/' , externalLoadsXml ] ) ;
end
forsim.set_use_visualizer( false ) ; % use visualizer while running (true or false)
forsim.print( fullfile( startDir , 'inputs' , forsimSettingsFileName ) ) ;

disp( 'Running Forsim Tool...' )
forsim.run() ;
disp( 'Simulation Complete.' )
