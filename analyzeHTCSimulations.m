%% Analyze High throughput Batch Simulation(s)
%==========================================================================
% Author: Matthew Blomquist
%
% Purpose: To process and analyze the data from batch simulations
%
% Output:
%   simData - structure with simulation data based on specified parameters
%
% Other .m files required:
%   batchProcessHTCSims.m
%
% Revision history:
%   v1      03-04-2023      First commit (MBB)
%
%==========================================================================

clear ; close all ;

%% Specify parameters (change these to desired parameters and file locations)
%----------------------------------------------------------------------------

% Results directory
Params.resultsDir = 'C:\Users\mbb201\Desktop\htcTKArelease\cohort2std_MBB' ;

% Name of results tar file
Params.resultsTarFile = 'job_results.tar.gz' ;

% Number of models that were run
Params.numModels = 200 ;

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
testDof = 'var' ;

% Specify flexion angle of knee during simulation
kneeFlexAngle = 30 ;

% Specify external load applied
externalLoad = 10 ;

% Choose model:
%   options: lenhart2015
Params.model = 'lenhart2015' ;

% Specify which data to pull:
%----------------------------

% Joint kinematics (6 degree-of-freedom)
%   Options: there are a lot, but the two common ones are 'knee_r' (right
%   tibiofemoral joint) and 'pf_r' (right patellofemoral joint)
Params.jointKinematics = { 'knee_r' , 'pf_r' } ;

% Muscle names
%   Options: there are a lot, but we usually only worry about the
%   quadriceps { 'recfem_r' , 'vasint_r' , 'vaslat_r' , 'vasmed_r' },
%   hamstrings { 'bflh_r' , 'bfsh_r' , 'semimem_r' , 'semiten_r' }, and
%   calf { 'gaslat_r' , 'gasmed_r' , 'soleus_r' } muscles
Params.muscleNames = { 'bflh_r' , 'bfsh_r' , 'semimem_r' , 'semiten_r' , ...
    'gaslat_r' , 'gasmed_r' , 'soleus_r' , ...
    'recfem_r' , 'vasint_r' , 'vaslat_r' , 'vasmed_r' } ;

% Muscle properties
%   Options: { 'force' , 'fiber_length' }
Params.muscleProperties = { 'force' } ;

% Ligament names
%   Options: { 'MCLd' , 'MCLs' , 'MCLp' , 'ACLpl' , 'ACLam' , 'PCLal' , ...
%   'PCLpm' , 'LCL' , 'PFL' , 'ITB' , 'PT' , 'lPFL' , 'mPFL' ,  'pCAP' }
Params.ligamentNames = { 'MCLd' , 'MCLs' , 'MCLp' , 'ACLpl' , 'ACLam' , 'PCLal' , ...
    'PCLpm' , 'LCL' , 'PFL' , 'ITB' , 'PT' , 'lPFL' , 'mPFL' ,  'pCAP' } ;

% Ligament properties
%   Options: { 'force_spring' , 'force_damping' , 'force_total' , 'length' ,
%   'lengthening_speed' , 'strain' , 'strain_rate' }
Params.ligamentProperties = { 'force_total' , 'length' , 'strain' } ;

% Compartment names for contact data
%   Common options are 'tf_contact' (tibiofemoral joint) or 'pf_contact'
%   (patellofemoral joint)
Params.contactCompartmentNames = { 'tf_contact' , 'pf_contact' } ;

% Contact data properties
%   Many options, most common are 'mean_pressure' , 'max_pressure' ,
%   'center_of_pressure_x' , 'center_of_pressure_y' , 'center_of_pressure_z' ,
%   'contact_force_x' , 'contact_force_y' , 'contact_force_z' ,
%   'contact_moment_x' , 'contact_moment_y' , 'contact_moment_z'
%   NOTE: x = anterior(+)-posterior, y = superior(+)-inferior, z =
%   lateral(+)-medial
Params.contactForces = { 'contact_force_x' , 'contact_force_y' , ...
    'contact_force_z' , 'contact_moment_x' , 'contact_moment_y' , ...
    'contact_moment_z' } ;

%% Set up for post-processing
%----------------------------

if strcmp( testDof , 'flex' ) % passive flexion
    trialNames = { [ 'flex_passive_0_' , num2str( kneeFlexAngle ) ] } ;
else % laxity tests
    trialNames = { [ 'lax_' , testDof , '_frc' , num2str( externalLoad ) , '_' , num2str( kneeFlexAngle ) ] } ;
end
Params.trialNames = trialNames ;

%% Run processLocalSims
%----------------------

if ~exist( fullfile( Params.resultsDir , trialNames{1} , 'simData.mat' ) , 'file' )
    simData = batchProcessHTCSims( Params ) ;
    save( fullfile( Params.resultsDir , trialNames{1} , 'simData.mat' ) , 'simData' ) % Save data
else
    load( fullfile( Params.resultsDir , trialNames{1} , 'simData.mat' ) ) ;
end

%% Plot Data
plotData = 1 ;

if plotData
    load( 'TI.mat' ) % load tolerance intervals
    figure() ; hold on ; grid on ; % Create figure

    % Plot tolerance intervals
    switch testDof
        case { 'var' , 'val' }
            xVals = [ 2.5 , 5 , 7.5 , 10 , 10 , 7.5 , 5 , 2.5 ] ;
            dofAxis = 'vv' ;
    end
    yVals = [ TI.( testDof ).( [ 'ang' , num2str( kneeFlexAngle ) ] ).min , fliplr( TI.( testDof ).( [ 'ang' , num2str( kneeFlexAngle ) ] ).max ) ] ;
    patch( xVals , yVals , 'k' , 'FaceAlpha' , 0.3 )

    % Plot simulation data
    loadIdcs = 201:301 ; % external load indices

    for iMdl = 1 : Params.numModels
        plot( simData.( trialNames{1} ).extLoad(loadIdcs,iMdl) , simData.( trialNames{1} ).kine.tf.(dofAxis)(loadIdcs,iMdl) , 'r' )
    end

    title( [ testDof , ' ' , num2str( kneeFlexAngle ) , 'deg' ] ) ; xlabel( 'Load (N)' ) ; ylabel( 'Translation (mm)' ) ;
    legend( 'Tolerance Interval' , 'Simulations' )
end