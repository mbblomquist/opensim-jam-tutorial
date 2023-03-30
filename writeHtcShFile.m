function outMsg = writeHtcShFile( Params )
%==========================================================================
% Purpose:  Write *.sh files for running jobs on HTC grid
%
% Inputs:
%   - Params: Structure of parameters to create .sh file
% 
% Outputs:
%   - Creates a .sh file in folder specified in Params
% 
% Written By: Joshua D. Roth (2019-07-01)
%
%--------------------------------------------------------------------------
% Revision history:
% -----------------
% v1    2019-01-03(JDR)     inital release
% v2    2021-01-14(MBB)     added ant, post, flex for osim
% v3    2023-03-29(MBB)     Update code for new OpenSim-JAM version
%
%==========================================================================

%% Unpack Params structure
%=========================

testDOF = Params.tempTestDOF ;
fld = Params.sharedDir ;

%% Write start of *.sh file
% =========================

if exist( fullfile( fld, Params.shName ), 'file' )
    overwriteShFileYesNo = input( 'Do you want to overwrite the existing sh file [y/n]?\n\n', 's' );
    if strcmp( overwriteShFileYesNo, 'n' )
        error( 'Action stopped to prevent overwriting existing *.sh file' )
    end
end

shFileId = fopen( fullfile( fld, Params.shName ), 'w' ) ;


% Standard info for all sh files
% -------------------------------
fprintf( shFileId, '%s\n\n' , '#!/bin/sh' ) ;

fprintf( shFileId, '%s\n' , '#=======================' ) ;
fprintf( shFileId, '%s\n' , '# Run Forward Simulation' ) ;
fprintf( shFileId, '%s\n\n' , '#=======================' ) ;

fprintf( shFileId, '%s\n\n' , 'JAM_DIR=$PWD/opensim-jam' ) ;

fprintf( shFileId, '%s\n' , '# Add path to OpenSim lib for executables to run' ) ;
fprintf( shFileId, '%s\n' , 'export LD_LIBRARY_PATH=$JAM_DIR:$LD_LIBRARY_PATH' ) ;
fprintf( shFileId, '%s\n\n' , 'export PATH=$JAM_DIR:$PATH' ) ;

fprintf( shFileId, '%s\n' , '# Extract Tar Files' ) ;
fprintf( shFileId, '%s%s\n' , 'tar -xzf ' , Params.sharedTarName ) ;
fprintf( shFileId, '%s%s\n\n' , 'tar -xzf ' , Params.opensimLibTarName ) ;

fprintf( shFileId, '%s\n' , '# Give permission to execute opensim-cmd' ) ;
fprintf( shFileId, '%s\n\n' , 'chmod +x opensim-jam/opensim-cmd' ) ;

%% Write call for laxity forward simulation
% =========================================

if isequal( testDOF , 'flex' ) % passive flexion
    fprintf( shFileId , '%s\n' , '# Run Passive Flexion Test' ) ;
    % Call for forward simulation executable
    fprintf( shFileId , '%s%s%s\n\n' , ...
        'opensim-cmd ' , 'run-tool ' , 'forsim_settings_flex_passive_0_90.xml' ) ;
else % laxity test
    % Indices of trials with this DOF in it
    trialIdc = contains( Params.trialNames , Params.tempTestDOF ) ;
    % Names of trials
    trialsToRun = Params.trialNames( trialIdc ~= 0 ) ;
    for iTrial = 1 : length( trialsToRun )
        tempTrialName = trialsToRun{iTrial} ;
        fprintf( shFileId, '%s\n', [ '# Run ' , tempTrialName ] ) ;
        % Call for forward simulation executable
        fprintf( shFileId , '%s%s%s\n\n' , ...
            'opensim-cmd ' , 'run-tool ' , [ 'forsim_settings_' , tempTrialName , '.xml' ] ) ;
    end
end

%% Write call to compress results files
% =====================================

fprintf( shFileId , '%s\n' , '# Compress Results Files' ) ;
resultsString = cell(1,1) ;
if isequal( testDOF , 'flex' ) % passive flexion
    resultsString{1} = 'flex_passive_0_90_ForceReporter_forces.sto flex_passive_0_90_states.sto' ;
else % laxity test
    % Indices of trials with this DOF in it
    trialIdc = contains( Params.trialNames , Params.tempTestDOF ) ;
    % Names of trials
    trialsToRun = Params.trialNames( trialIdc ~= 0 ) ;
    for iTrial = 1 : length( trialsToRun )
        tempTrialName = trialsToRun{iTrial} ;
        resultsString{1} = [ resultsString{1} , ...
            tempTrialName , '_ForceReporter_forces.sto ' , ...
            tempTrialName , '_states.sto ' ] ;
    end
end

% Call to tar results files
fprintf( shFileId, '%s%s\n', ...
    'tar czf results.tar.gz ', resultsString{1} ) ;


%% Set output message
%====================

fclose( shFileId ) ;

outMsg = [ 'Wrote ', Params.shName, ' to ', fld, '.' ] ;

end

%% Extra stuff to use for COMAK simulations...?

% % =================
% % comak simulation
% % =================
% case 'comak_simm'
%     %for troubleshooting purposes only
%     disp( 'Write comak sh file' )
%
%     %standard info for all sh files
%     %------------------------------
%     fprintf( shFileId, '%s\n', '#!/bin/sh' ) ;
%     fprintf( shFileId, '%s\n\n', '#Run IK AND COMAK Simulation' ) ;
%     fprintf( shFileId, '%s\n\n', '#===========================' ) ;
%
%     fprintf( shFileId, '%s\n', '#Extract Files' ) ;
%     fprintf( shFileId, '%s%s\n\n', 'tar -xzf ', sharedInputsName ) ;
%
%     fprintf( shFileId, '%s\n', '#add present working directory (pwd) to path so c-libraries can be found by exes' ) ;
%     fprintf( shFileId, '%s\n', 'LD_LIBRARY_PATH="$(pwd)"' ) ;
%     fprintf( shFileId, '%s\n', 'export LD_LIBRARY_PATH' ) ;
%     fprintf( shFileId, '%s\n\n', 'export PATH=$PATH:$(pwd)' ) ;
%
%     %write call for inverse kinematics
%     %---------------------------------
%     if any( strcmp( Params.simSteps, 'ik' ) )
%         fprintf( shFileId, '%s\n', '#inverse kinematics' ) ;
%         fprintf( shFileId, '%s\n', '#==================' ) ;
%         fprintf( shFileId, '%s\n', '#exe params scaledModel moCapData outputMot' ) ;
%         %call for inverse kinematics exe
%         % exe params model moCapData resultsFile
%         fprintf( shFileId, '%s%s%s%s%s%s%s%s\n\n', ...
%             './invkin ', Params.invKinParams, ' ', mdlName, ' ', Params.inputTrc, ' ', [ resultsName, 'Ik' ] ) ;
%     end
%
%     %write call for comak and/or contact
%     %-----------------------------------
%     if any( strcmp( Params.simSteps, 'comak' ) )
%         fprintf( shFileId, '%s\n', '#comak' ) ;
%         fprintf( shFileId, '%s\n', '#=====' ) ;
%         %call for forward simulation exe
%         % exe params scaledModel invKinResultMot outputFile startTime endTime grfSegmentFp1 grfSegmentFp2 grfSegmentFp3
%         switch Params.activity
%             case 'walk'
%                 fprintf( shFileId, '%s\n', '#exe params scaledModel invKinResultMot outputFile startTime endTime grfSegmentFp1 grfSegmentFp2 grfSegmentFp3' ) ;
%                 fprintf( shFileId, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n\n', ...
%                     './mfd ', Params.comakParams, ' ', mdlName, ' ', [ resultsName, 'Ik' ], ' ', ...
%                     resultsName, ' ', Params.startTime, ' ', Params.endTime, ' ', ...
%                     Params.grfSegFrcPlate1, ' ', Params.grfSegFrcPlate2, ' ', Params.grfSegFrcPlate3 ) ;
%             case 'stairs'
%                 fprintf( shFileId, '%s\n', '#exe params scaledModel invKinResultMot outputFile startTime endTime grfSegmentFp1 grfSegmentFp2' ) ;
%                 fprintf( shFileId, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n\n', ...
%                     './mfd ', Params.comakParams, ' ', mdlName, ' ', [ resultsName, 'Ik' ], ' ', ...
%                     resultsName, ' ', Params.startTime, ' ', Params.endTime, ' ', ...
%                     Params.grfSegFrcPlate1, ' ', Params.grfSegFrcPlate2 ) ;
%         end
%     end
%
%     if any( strcmp( Params.simSteps, 'contact' ) )
%         fprintf( shFileId, '%s\n', '#contact' ) ;
%         fprintf( shFileId, '%s\n', '#=======' ) ;
%         fprintf( shFileId, '%s\n', '#exe params scaledModel motFile' ) ;
%         %call for contact exe
%         % exe params modelresultsFile
%         fprintf( shFileId, '%s%s%s%s%s%s\n\n', ...
%             './contact ', cntParams,' ', mdlName, ' ', resultsName ) ;
%     end
%
%     %write call for compress output files
%     %------------------------------------
%     fprintf( shFileId, '%s\n', '#Compress Results Files' ) ;
%     if any( strcmp( Params.simSteps, 'contact' ) )
%         fprintf( shFileId, '%s%s%s%s%s%s%s\n', ...
%             'tar czf job_results.tar.gz ', resultsName, 'Ik.mot ', resultsName, '.mot ', resultsName, '.h5' ) ;
%     else
%         fprintf( shFileId, '%s%s%s%s%s\n', ...
%             'tar czf job_results.tar.gz ', resultsName, 'Ik.mot ', resultsName, '.mot ') ;
%     end