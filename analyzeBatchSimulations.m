%% Analyze Batch Simulation(s)
%==========================================================================
% Author: Matthew Blomquist
%
% Purpose: To process and analyze the data from batch simulations
%
% Output:
%   simData - structure with simulation data based on specified parameters
%
% Other .m files required:
%   processLocalSims.m
%
% Revision history:
%   v1      11-17-2022      First commit (MBB)
%
%==========================================================================

clear ; close all ;

%% Specify parameters (change these to desired parameters and file locations)
%----------------------------------------------------------------------------

% Results directory
Params.resultsDir = 'C:\Users\mbb201\Documents\MATLAB\Research\opensimModeling\opensim-jam-code\results' ;

% Number of models that were run
Params.numModels = 10 ;

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
kneeFlexAngle = 25 ;

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
Params.muscleProperties = { 'force' , 'fiber_length' } ;

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
simData = batchProcessSims( Params ) ;

%% Load Stoch Data File
%-----------------------

load( 'stochModels\stochMdlParams.mat' )

%% Analyze Data

% Compute varus laxity values for each model
%   Laxity = kinematics change from end of load to start of load
varLaxity = simData.lax_var_frc10_25.kine.tf.vv( end , : ) - ...
    simData.lax_var_frc10_25.kine.tf.vv( 201 , : ) ;

% Compute LCL slack length changes from baseline
%   Multiply by 100 to change to %
LCLchanges = StochMdlParams.stoch.LCL.slack_length * 100 ;

figure() ; grid on ; hold on ;
scatter( varLaxity , LCLchanges , 'filled' ) 
xlabel( 'Varus Laxity (deg)' ) ; ylabel( 'LCL Slack Length Change (%)' ) ;