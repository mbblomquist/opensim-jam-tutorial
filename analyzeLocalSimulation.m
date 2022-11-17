%% Analyze Local Simulation

clear ; close all ;

%% Specify parameters (change these to desired parameters and file locations)
%============================================================================

% Choose file from pop-up window
[ file , path ] = uigetfile( '.\results\*.sto' , 'MultiSelect' , 'on' ) ;

% Choose model:
%   options: lenhart2015
Params.model = 'lenhart2015' ;

% Specify which data to pull:
%----------------------------

% Joint kinematics (6 degree-of-freedom)
%   Options: there are a lot, but the two common ones are 'knee_r' (right
%   tibiofemoral joint) and 'pf_r' (right patellofemoral joint)
Params.jointKinematics = { 'knee_r' , 'pf_r' } ;

% Muscle forces
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



%% Run post processing (no need to change anything in this portion)
%==================================================================

% Specify trialnames to run
Params.resultsDir = path ;

if iscell(file)
    trialNames = cell( length( file ) , 1 ) ;
    splitFile = split( file' , '_' , 1 ) ;
else
    trialNames = cell( 1 , 1 ) ;
    splitFile = split( file , '_' ) ;
end
for iTrial = 1 : size( splitFile , 2 )
    trialNames{iTrial} = [ splitFile{1,iTrial} , '_' , splitFile{2,iTrial} , '_' , splitFile{3,iTrial} , '_' , splitFile{4,iTrial} ] ;
end
Params.trialNames = trialNames ;

% Run processLocalSims
simData = processLocalSims( Params ) ;
