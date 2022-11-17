%% Analyze Local Simulation(s)
%==========================================================================
% Author: Matthew Blomquist
%
% Purpose: To analyze the results of forward simulations and extract the
%   parameters of interest
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

% Pull results directory
Params.resultsDir = path ;

% Specify trialnames to run
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

%% Run processLocalSims
%----------------------
simData = processLocalSims( Params ) ;

%% ADD CODE HERE TO EXTRACT DATA OF INTEREST
%--------------------------------------------

% Switch this between 1 and 0. If you want to plot the results for example
% 1 (below), then set this to 1. If you don't want to plot the results for
% example 1, then set this to 0.
plotExample1 = 0 ;

if plotExample1
    % Example: If you wanted to plot LCL tension vs flexion angle
    figure() ; hold on ; grid on ;

    flexionAngles = simData.flex_passive_0_90.kine.tf.fe ; % 202 x 1 double of flexion angles
    LCLtension = simData.flex_passive_0_90.lig.LCL.allStrands.force_total ; % 202 x 1 double of LCL tensions

    plot( flexionAngles , LCLtension )

    xlabel( 'Flexion Angle (^o)' ) ; ylabel( 'LCL Tension (N)' ) ;
end

% Switch this between 1 and 0. If you want to plot the results for example
% 1 (below), then set this to 1. If you don't want to plot the results for
% example 1, then set this to 0.
plotExample2 = 0 ;

if plotExample2
    % Example: If you wanted to plot tension of each MCLs strand vs flexion angle
    figure() ; hold on ; grid on ;

    flexionAngles = simData.flex_passive_0_90.kine.tf.fe ; % 202 x 1 double of flexion angles

    % There are 6 strands for the MCLs. The best way to plot would be to
    % loop through them all using a for loop, but if you aren't as
    % comfortable with MATLAB, then this way should make more sense:
    MCLs1_tension = simData.flex_passive_0_90.lig.MCLs.strand1.force_total ; % 202 x 1 double
    MCLs2_tension = simData.flex_passive_0_90.lig.MCLs.strand2.force_total ; % 202 x 1 double
    MCLs3_tension = simData.flex_passive_0_90.lig.MCLs.strand3.force_total ; % 202 x 1 double
    MCLs4_tension = simData.flex_passive_0_90.lig.MCLs.strand4.force_total ; % 202 x 1 double
    MCLs5_tension = simData.flex_passive_0_90.lig.MCLs.strand5.force_total ; % 202 x 1 double
    MCLs6_tension = simData.flex_passive_0_90.lig.MCLs.strand6.force_total ; % 202 x 1 double

    % Plot each strand's tension vs flexion angle. They will all stay on
    % the same plot since we added "hold on" on line 92.
    plot( flexionAngles , MCLs1_tension )
    plot( flexionAngles , MCLs2_tension )
    plot( flexionAngles , MCLs3_tension )
    plot( flexionAngles , MCLs4_tension )
    plot( flexionAngles , MCLs5_tension )
    plot( flexionAngles , MCLs6_tension )

    xlabel( 'Flexion Angle (^o)' ) ; ylabel( 'MCLs Tension (N)' ) ;
    legend( { 'Strand 1' , 'Strand 2' , 'Strand 3' , 'Strand 4' , 'Strand 5' , 'Strand 6' } )

end

% Make various plots here to get more comfortable with the simData
% structure and plotting in MATLAB! Some things you could try (PLEASE reach
% out to me if you are unable to get any of these things, I'm happy to
% help!)

% 1) Tibiofemoral superior-inferior (S-I) contact force vs flexion angle
% (S-I is the "y" direction in OpenSim)

% 2) Anterior-Posterior, Compression-Distraction, and Medial-Lateral
% kinematics on the same plot vs flexion angle

% 3) LCL length vs tension for each strand (there are 4 strands for the
% LCL)

% 4) Come up with something on your own and let me know and we can double
% check that we have the same results!

% 5) Now, run a laxity test! Run an anterior laxity test, 90deg flexion, 
% 100N load. Plot external anterior load vs anterior-posterior kinematics