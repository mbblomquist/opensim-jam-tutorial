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
% v1    11-17-2022  First commit (MBB)
% v2    03-30-2023  Updated for local or HT (MBB)
%
%==========================================================================
clc ; clear ; close all ;

%% ======================= Specify Settings ===========================
% Look through this entire section to change parameters to whatever you
% wish to run. This entire section is creating a structure full of
% parameters to analyze the simulations and extract parameters that you
% want
% =====================================================================

% ------------------------------------------------------------------------
% -------------------------- SPECIFY LOCAL VS HT -------------------------
% ------------------------------------------------------------------------

% Specify whether to run locally or whether you will run the models on the
% high-throughput grid
% Options: 'local' or 'HT'
Params.localOrHT = 'local' ;

% Set base name of output folder where models, exectuables, and inputs
% should be created
Params.baseOutDir = 'C:\Users\mbb201\Desktop\htcTKArelease\localTest' ;

% Name of results file (no need to change this) - This is for if you ran on
% the HT grid
Params.resultsTarFile = 'results.tar.gz' ;

% ------------------------------------------------------------------------
% ----------------------- SPECIFY MODEL PARAMETERS -----------------------
% ------------------------------------------------------------------------
% Specif the model parameters that you want to extract

% Number of models that were run
Params.numModels = 1 ;

% Base model used. Options are in lenhart2015 folder
%   Current options =
%       'lenhart2015' (intact model)
%       'lenhart2015_implant' (TKA model - implants and no ACL or MCLd)
%       'lenhart2015_BCRTKA' (BCR-TKA model - implants with ACL and MCLd)
%       'lenhart2015_SarahISTA_PCL' (for Sarah's ISTA abstract)
Params.baseMdl = 'lenhart2015' ;

% Joint kinematics (6 degree-of-freedom)
%   Options: there are a lot, but the two common ones are 'knee_r' (right
%   tibiofemoral joint) and 'pf_r' (right patellofemoral joint)
Params.jointKinematics = { 'knee_r' , 'pf_r' } ;

% Muscle names
%   Options: there are a lot, but we usually only worry about the
%   quadriceps { 'recfem_r' , 'vasint_r' , 'vaslat_r' , 'vasmed_r' },
%   hamstrings { 'bflh_r' , 'bfsh_r' , 'semimem_r' , 'semiten_r' }, and
%   calf { 'gaslat_r' , 'gasmed_r' , 'soleus_r' } muscles
Params.muscleNames = { } ;

% Muscle properties
%   Options: { 'force' , 'fiber_length' }
Params.muscleProperties = { 'force' } ;

% Ligament names
%   Options: 'allLigs' to extract data on all of the ligaments in the model
%                           OR
%     [cell array with each ligament you want to extract]
%       'MCLd' , 'MCLs', 'MCLp', 'ACLpl' , 'ACLam' , 'LCL', 'ITB', 'PFL',
%       'pCAP', 'PCLpm', 'PCLal', 'PT', 'lPFL', 'mPFL'
Params.ligamentNames = 'allLigs' ;

% Ligament properties
%   Options: { 'force_spring' , 'force_damping' , 'force_total' , 'length' ,
%   'lengthening_speed' , 'strain' , 'strain_rate' }
Params.ligamentProperties = { 'force_total' , 'length' , 'strain' } ;

% Compartment names for contact data
%   Common options are 'tf_contact' (tibiofemoral joint) or 'pf_contact'
%   (patellofemoral joint)
Params.contactCompartmentNames = { 'tf_contact' , 'pf_contact' } ;
% Params.contactCompartmentNames = { 'tf_contact_medial' , 'tf_contact_lateral' , 'pf_contact' } ;

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
Params.testDOFs = { 'var' } ;

% Specify flexion angle(s) of knee during each simulation [cell array]
%   For passive flexion, specify the end flexion angle (starts at 0)
%   Each testDOF will be run at each flexion angle (so the total number of
%   simulations will be length(testDOFs) * length( kneeFlexAngles )
Params.kneeFlexAngles = { 30 } ;

% Specify external load(s) applied, one for each testDOFs [cell array]
%   Put 0 if passive flexion test
%   Keep this number positive
Params.externalLoads = { 10 } ;

%% ======================== Compute Trial Name ===========================
% Computes the trial name(s) that will be used for running through the
% results files
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
    if isequal( Params.testDOFs{iDOF} , 'flex' ) % if passive flexion test
        for iAng = 1 : length( Params.kneeFlexAngles ) % loop through flexion angles
            Params.trialNames{ trialCounter } = ...
                [ 'flex_passive_0_' , num2str( Params.kneeFlexAngles{iAng} ) ] ;
            trialCounter = trialCounter + 1 ;
        end
    elseif contains( Params.testDOFs{iDOF} , '-' ) % if combined loading
        for iAng = 1 : length( Params.kneeFlexAngles ) % loop through flexion angles
            Params.trialNames{ trialCounter } = ...
                [ 'lax_' , Params.testDOFs{iDOF} , '_frc' , num2str( Params.externalLoads{iDOF}(1) ) , '-' , num2str( Params.externalLoads{iDOF}(2) ) , '_' , num2str( Params.kneeFlexAngles{iAng} ) ] ;
            trialCounter = trialCounter + 1 ;
        end
    else % if laxity test
        for iAng = 1 : length( Params.kneeFlexAngles ) % loop through flexion angles
            Params.trialNames{ trialCounter } = ...
                [ 'lax_' , Params.testDOFs{iDOF} , '_frc' , num2str( Params.externalLoads{iDOF} ) , '_' , num2str( Params.kneeFlexAngles{iAng} ) ] ;
            trialCounter = trialCounter + 1 ;
        end
    
    end
end

Params.numTrials = length( Params.trialNames ) ;

clear trialCounter iDOF iAng

%% Run processLocalSims
%----------------------
simData = batchProcessSims( Params ) ;
