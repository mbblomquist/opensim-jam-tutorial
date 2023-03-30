function outMsg = writeHtcSubFile( Params )
%==========================================================================
% Purpose:  Create *.sub file to launch jobs on the UW CHTC network
% 
% Inputs:
%   - Params: Structure of parameters to create .sub file
% 
% Outputs:
%   - Creates a .sub file in folder specified in Params
% 
% Written By: Joshua D. Roth (2019-07-01)
%
%--------------------------------------------------------------------------
% Revision history:
% -----------------
% v1    2019-07-01(JDR)     inital release
% v2    2023-03-29(MBB)     Update code for new OpenSim-JAM version
% 
%==========================================================================

%% ============== Define params ============== %%

fld = Params.trialDir ;
jobName = [ Params.studyId, '_', Params.date ] ;


%% ============== Write *.sub file ============== %%

subFileId = fopen( fullfile( fld, [ jobName, '.sub' ] ), 'w+' );

fprintf( subFileId, '%s\n', '# HTCondor Submit File' ) ;
fprintf( subFileId, '%s\n\n', '#==============================================================================' ) ;

fprintf( subFileId, '%s\n\n', [ 'job_name = ', jobName ] ) ;

fprintf( subFileId, '%s\n', 'universe = vanilla' ) ;
fprintf( subFileId, '%s\n', 'log = results/$(Process)/$(job_name)_$(Process)_$(Cluster).log' ) ;
fprintf( subFileId, '%s\n', 'error = results/$(Process)/$(job_name)_$(Process).err' ) ;
fprintf( subFileId, '%s\n\n', 'output = results/$(Process)/$(job_name)_$(Process).out' ) ;

fprintf( subFileId, '%s\n', [ 'executable = ./shared/', Params.shName ] ) ;
fprintf( subFileId, '%s\n\n', 'arguments = $(Process)' ) ;

fprintf( subFileId, '%s\n', 'should_transfer_files = YES' ) ;
fprintf( subFileId, '%s\n', 'transfer_input_files = ../input/$(Process)/, shared/' ) ;
fprintf( subFileId, '%s\n', 'when_to_transfer_output = ON_EXIT' ) ;
fprintf( subFileId, '%s\n\n', 'transfer_output_remaps = "results.tar.gz = results/$(Process)/results.tar.gz"' ) ;

fprintf( subFileId, '%s\n', 'request_cpus = 1' ) ;
fprintf( subFileId, '%s\n', 'request_memory = 4GB' ) ;
fprintf( subFileId, '%s\n\n', 'request_disk = 4GB' ) ;

fprintf( subFileId, '%s\n\n', 'Requirements = (HasCHTCSoftware == true) && (OpSysMajorVer == 7)' ) ;

fprintf( subFileId, '%s\n', '+WantFlocking = true' ) ;
fprintf( subFileId, '%s\n\n', '+WantGlideIn = true' ) ;

fprintf( subFileId, '%s\n\n', [ 'queue ', num2str( Params.numModels ) ] ) ;

%% ============== Close *.sub file ============== %%

outMsg = [ 'Created ' , fullfile( fld, [ jobName, '.sub' ] ) , ' file' ] ;
fclose( subFileId ) ;

end