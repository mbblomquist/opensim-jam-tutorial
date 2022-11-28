%% Analyze Local Simulation(s)
%==========================================================================
% Written by: Josh Roth
% Revised by: Matthew Blomquist
%
% Purpose: To create a large structure by processing local forward
%   simulations and saving kinematics and kinetics of interest
%
% Output:
%   StoData - structure with simulation data based on specified parameters
%
% Other .m files required:
%   none (nested files at bottom)
%
% Revision history:
%   v1      11-17-2022      First commit (MBB)
%
%==========================================================================

function StoData = batchProcessSims( Params )

% Results Directory
resultsDir = Params.resultsDir ;

%% Read results
%---------------
trialNames = Params.trialNames ;

for iMdl = 1 : Params.numModels

    % Read sto files and organize results
    for iTrial = 1 : length( trialNames )

        % Names of results files
        laxResultFile = [ trialNames{ iTrial } , '_' , num2str( iMdl ) , '_states.sto' ] ;
        frcResultsFile = [ trialNames{ iTrial } , '_' , num2str( iMdl ) , '_ForceReporter_forces.sto' ] ;

        % Use read_opensim_mot function to read .sto files
        [ tempData.kine.data, tempData.kine.labels, tempData.kine.header ] = ...
            read_opensim_mot( fullfile( resultsDir, laxResultFile ) ) ; % Kinematics files
        [ tempData.frc.data, tempData.frc.labels, tempData.frc.header ] = ...
            read_opensim_mot( fullfile( resultsDir, frcResultsFile ) ) ; % Force files

        % Find indices for each trial
        StoData.idx = parseOpenSimSto( tempData, Params ) ;

        %============%
        % Kinematics %
        %============%
        if isfield( StoData.idx , 'kine' ) % Make sure kinematics is a field

            % Loop through joints
            joints = fieldnames( StoData.idx.kine ) ;
            for iJoint = 1 : length( joints )

                % Switch extracting data based on joint
                switch joints{ iJoint }
                    case 'tf'

                        % Fieldnames of the tibiofemoral degrees of freedom
                        kineDofs = fieldnames( StoData.idx.kine.tf ) ;

                        % Loop through degrees of freedom
                        for iDof = 1 : length( kineDofs )
                            tempDof = kineDofs{ iDof } ;

                            % Convert data separately based on DOF
                            switch tempDof
                                case { 'vv' , 'ie' , 'fe' } % rotation values
                                    scaleFactor = 180 / pi() ; % radians to degrees
                                    StoData.( trialNames{ iTrial } ).kine.tf.( tempDof )( : , iMdl ) = ...
                                        tempData.kine.data( : , StoData.idx.kine.tf.(tempDof) ) * scaleFactor ;

                                case { 'ap' , 'pd' , 'ml' }
                                    scaleFactor = 1000 ; % m to mm
                                    % For translations, we need to account for which
                                    %   flexion angle the knee is at and apply
                                    %   a rotation matrix to it
                                    flexAngles = StoData.( trialNames{ iTrial } ).kine.tf.fe( : , iMdl ) ;
                                    numFlexAngles = length( flexAngles ) ;

                                    switch tempDof
                                        case 'ap'
                                            axisNum = 1 ; % x axis
                                        case 'pd'
                                            axisNum = 2 ; % y axis
                                        case 'ml'
                                            axisNum = 3 ; % z axis
                                    end

                                    newTransValues = zeros( numFlexAngles , 1 ) ; % initialize

                                    % Loop through flexion angles
                                    for iAngle = 1 : numFlexAngles
                                        % Create rotation matrix based on flexion angle
                                        rotMatrix = ...
                                            [   cosd( iAngle )  -sind( iAngle ) 0 ; ...
                                            sind( iAngle )  cosd( iAngle )  0 ; ...
                                            0               0               1 ] ;

                                        % Create vector of translations from OpenSim
                                        transValues = [ tempData.kine.data( iAngle , StoData.idx.kine.tf.ap ) ; ...
                                            tempData.kine.data( iAngle , StoData.idx.kine.tf.pd ) ; ...
                                            tempData.kine.data( iAngle , StoData.idx.kine.tf.ml ) ] ;

                                        % Compute new translation values
                                        tempValues = rotMatrix * transValues ;

                                        % Extract DOF of interest
                                        newTransValues( iAngle , 1 ) = tempValues( axisNum ) ;
                                    end % iAngle

                                    % Add to structure
                                    StoData.( trialNames{ iTrial } ).kine.tf.( tempDof )( : , iMdl ) = ...
                                        scaleFactor * ( newTransValues ) ;

                            end % switch tempDofName
                        end % for tempDof

                    case 'pf'
                        % ===========================================
                        % MIGHT NEED TO CHANGE THIS
                        %   Also, need to add offsets of translations
                        % ===========================================

                        % Fieldnames of the patellofemoral degrees of freedom
                        kineDofs = fieldnames( StoData.idx.kine.pf ) ;

                        for iDof = 1 : length( kineDofs )
                            tempDof = kineDofs{ iDof } ;

                            StoData.( trialNames{ iTrial } ).kine.pf.( tempDof )( : , iMdl ) = ...
                                tempData.kine.data( : , StoData.idx.kine.pf.(tempDof) ) ;

                        end % for iDof

                end % switch joints{ iJoint }
            end % for iJoint
        end % if isField

        %=============%
        % Muscle data %
        %=============%

        if isfield( StoData.idx , 'msl' ) % Make sure muscles are a field
            mslNames = fieldnames( StoData.idx.msl ) ;
            mslProp = fieldnames( StoData.idx.msl.( mslNames{1} ) ) ;

            % Loop through muscles
            for iMsl = 1 : length( mslNames )
                tempMslName = mslNames{ iMsl } ;

                % Loop through muscle properties
                for iProp = 1 : length( mslProp )

                    switch mslProp{ iProp }
                        case 'force' % extract from tempData.frc
                            StoData.( trialNames{ iTrial } ).msl.( tempMslName ).frc( : , iMdl ) = ...
                                tempData.frc.data( : , StoData.idx.msl.( tempMslName ).force ) ; % muscle force
                        case 'fiber_length' % extract from tempData.kine
                            StoData.( trialNames{ iTrial } ).msl.( tempMslName ).length( : , iMdl ) = ...
                                tempData.kine.data( : , StoData.idx.msl.( tempMslName ).fiber_length ) ; % muscle fiber length
                    end % switch mslProp{iProp}

                end % for iProp
            end % iMsl
        end % isfield

        %===============%
        % Ligament data %
        %===============%

        if isfield( StoData.idx , 'lig' ) % Make sure ligaments are a field
            ligNames = fieldnames( StoData.idx.lig ) ;
            ligProp = fieldnames( StoData.idx.lig.( ligNames{1} ) ) ;

            % Loop through ligaments
            for iLig = 1 : length( ligNames )
                tempLigName = ligNames{ iLig } ;
                tempNumStrands = size( StoData.idx.lig.( ligNames{ iLig } ).( ligProp{ 1 } ) , 1 ) ;

                % Loop through ligament properties
                for iProp = 1 : length( ligProp )
                    tempProp = ligProp{ iProp } ;
                    for iStrand = 1 : tempNumStrands
                        StoData.( trialNames{ iTrial } ).lig.( tempLigName ).( [ 'strand' , num2str( iStrand ) ] ).( tempProp )( : , iMdl ) = ...
                            tempData.frc.data( : , StoData.idx.lig.( ligNames{ iLig } ).( tempProp )( iStrand ) ) ;

                        % Add data for addition structure to compute allStrands value
                        tempLigData.( tempProp )( : , iStrand ) = ...
                            StoData.( trialNames{ iTrial } ).lig.( tempLigName ).( [ 'strand' , num2str( iStrand ) ] ).( tempProp )( : , iMdl ) ;
                    end

                    % Sum or average (depending on property) to get one value for ligament (allStrands)
                    if contains( tempProp , 'force' ) % Sum the force_total, force_damping, or force_spring
                        StoData.( trialNames{ iTrial } ).lig.( tempLigName ).allStrands.( tempProp )( : , iMdl ) = ...
                            sum( tempLigData.( tempProp ) , 2 ) ;
                    else % Average the length, lengthening_speed, strain, strain_rate
                        StoData.( trialNames{ iTrial } ).lig.( tempLigName ).allStrands.( tempProp )( : , iMdl ) = ...
                            mean( tempLigData.( tempProp ) , 2 ) ;
                    end % if contains

                end % for iProp

                clear tempLigData

            end % for iLig
        end % isfield

        %====================%
        % Contact force data %
        %====================%

        if isfield( StoData.idx , 'cnt' ) % Make sure cnt are a field
            compNames = fieldnames( StoData.idx.cnt ) ;
            cntLoadNames = fieldnames( StoData.idx.cnt.( compNames{1} ) ) ;

            % Loop through contact compartments
            for iCntComp = 1 : length( compNames )

                % Loop through contact loads
                for iCntLoad = 1 : length( cntLoadNames )
                    StoData.( trialNames{ iTrial } ).cnt.( compNames{ iCntComp } ).( cntLoadNames{ iCntLoad } )( : , iMdl ) = ...
                        tempData.frc.data( : , StoData.idx.cnt.( compNames{ iCntComp } ).( cntLoadNames{ iCntLoad } ) ) ;

                end % for iCntLoad
            end % for iCntComp
        end % isfield

        %====================%
        % External load data %
        %====================%

        if isfield( StoData.idx , 'extLoad' )
            extLoads = fieldnames( StoData.idx.extLoad ) ;

            % Loop through external load variables (forces, moments, force points)
            for iExtLoad = 1 : length( extLoads )
                tempExtLoadData.( extLoads{iExtLoad} ) = ...
                    tempData.frc.data( : , StoData.idx.extLoad.( extLoads{iExtLoad} ) ) ;

            end % for iExtLoad

            % Anterior-Posterior = forces
            if contains( trialNames{ iTrial } , 'ant' ) || contains( trialNames{ iTrial } , 'post' )
                StoData.( trialNames{ iTrial } ).extLoad( : , iMdl ) = ...
                    sqrt( tempExtLoadData.load_Fx.^2 + tempExtLoadData.load_Fy.^2 + tempExtLoadData.load_Fz.^2 ) ;

                % Varus-Valgus = forces * point load
            elseif contains( trialNames{ iTrial } , 'var' ) || contains( trialNames{ iTrial } , 'val' )
                StoData.( trialNames{ iTrial } ).extLoad( : , iMdl ) = ...
                    sqrt( tempExtLoadData.load_Fx.^2 + tempExtLoadData.load_Fy.^2 + tempExtLoadData.load_Fz.^2 ) .* -tempExtLoadData.load_py ;

                % Internal-External Rotation = Moments
            elseif contains( trialNames{ iTrial } , 'ir' ) || contains( trialNames{ iTrial } , 'er' )
                StoData.( trialNames{ iTrial } ).extLoad( : , iMdl ) = ...
                    sqrt( tempExtLoadData.load_Tx.^2 + tempExtLoadData.load_Ty.^2 + tempExtLoadData.load_Tz.^2 ) ;

                % Distraction-Compression = Forces
            elseif contains( trialNames{ iTrial } , 'dist' ) || contains( trialNames{ iTrial } , 'comp' )
                StoData.( trialNames{ iTrial } ).extLoad( : , iMdl ) = ...
                    sqrt( tempExtLoadData.load_Fx.^2 + tempExtLoadData.load_Fy.^2 + tempExtLoadData.load_Fz.^2 ) ;

                % Passive flexion = 0
            elseif contains( trialNames{ iTrial } , 'flex' )
                StoData.( trialNames{ iTrial } ).extLoad( : , iMdl ) = 0 ;
            end

        end % isfield

    end % for iTrial

end % for iMdl

end % function


%% Nested Functions
%==================

function Idx = parseOpenSimSto( StoDataStruct, Params )
%=======================================================================
% Purpose: Parse sto-file outputs from OpenSim-JAM simulations using
% simpleKnee_lab04.osim model (has 1 strand for MCL, LCL, ACL, and PCL)
%
% Inputs:
% 	StoDataStruct: data and labels from OpenSim-JAM (structure)
% 	Params: parameters for simulation (structure)
%
% Outputs:
% 	Idx: index values for data columns (structure)
%
% Example:
% 	Idx = parseOpenSimSto( StoDataStruct, Params )
% 	- Important parameter is Params.simType, which determines which data is
% 	pulled
%
% Other m-files required:
%	None
%
% Written By: Josh Roth
% Revised By: Matthew Blomquist
%
% Project: OpenSim-JAM
%
% --------------------------------------------
%
% Revision History:
% v1	2020-11-11	initial release (JDR)
% v2    2022-11-17  updated code for new JAM release (MBB)
%=======================================================================

%% Parse sto file
%===============

switch Params.model
    case 'lenhart2015'
        kineLabels = StoDataStruct.kine.labels ;
        frcLabels = StoDataStruct.frc.labels ;

        % FIND INDICES:
        %==============

        % Time
        %------
        Idx.time = find( strcmp( kineLabels, 'time' ) ) ;

        % Kinematics
        %------------
        %   Tibiofemoral Joint
        if any( contains( Params.jointKinematics , 'knee_r' ) )
            Idx.kine.tf.fe = find( strcmp( kineLabels, '/jointset/knee_r/knee_flex_r/value' ) ) ;
            Idx.kine.tf.vv = find( strcmp( kineLabels, '/jointset/knee_r/knee_add_r/value' ) ) ;
            Idx.kine.tf.ie = find( strcmp( kineLabels, '/jointset/knee_r/knee_rot_r/value' ) ) ;
            Idx.kine.tf.ap = find( strcmp( kineLabels, '/jointset/knee_r/knee_tx_r/value' ) ) ;
            Idx.kine.tf.pd = find( strcmp( kineLabels, '/jointset/knee_r/knee_ty_r/value' ) ) ;
            Idx.kine.tf.ml = find( strcmp( kineLabels, '/jointset/knee_r/knee_tz_r/value' ) ) ;
        end
        %   Patellofemoral Joint
        if any( contains( Params.jointKinematics , 'pf_r' ) )
            Idx.kine.pf.fe = find( strcmp( kineLabels, '/jointset/pf_r/pf_flex_r/value' ) ) ;
            Idx.kine.pf.rot = find( strcmp( kineLabels, '/jointset/pf_r/pf_rot_r/value' ) ) ;
            Idx.kine.pf.tilt = find( strcmp( kineLabels, '/jointset/pf_r/pf_tilt_r/value' ) ) ;
            Idx.kine.pf.ap = find( strcmp( kineLabels, '/jointset/pf_r/pf_tx_r/value' ) ) ;
            Idx.kine.pf.pd = find( strcmp( kineLabels, '/jointset/pf_r/pf_ty_r/value' ) ) ;
            Idx.kine.pf.ml = find( strcmp( kineLabels, '/jointset/pf_r/pf_tz_r/value' ) ) ;
        end

        % Muscles
        %---------
        mslNames = Params.muscleNames ;
        numMsls = length( mslNames ) ;
        mslProp = Params.muscleProperties ;
        numProps = length( mslProp ) ;

        % Go through muscles
        for iMsl = 1 : numMsls
            mslId = mslNames{ iMsl } ;

            % Go through muscle properties
            for iProp = 1 : numProps
                switch mslProp{ iProp }
                    case 'force'
                        Idx.msl.( mslId ).( mslProp{ iProp } ) = ...
                            find( strcmp( frcLabels, mslId ) ) ; % Look in frcLabels
                    case 'fiber_length'
                        Idx.msl.( mslId ).( mslProp{ iProp } ) = ...
                            find( contains( kineLabels, mslId ) & endsWith( kineLabels , mslProp{iProp} ) ) ; % Look in kineLabels
                end
            end

        end

        % Ligaments
        %-----------
        ligNames = Params.ligamentNames ;
        numLigs = length( ligNames ) ;
        ligProp = Params.ligamentProperties ;
        numProps = length( ligProp ) ;

        % Go through ligaments
        for iLig = 1 : numLigs
            ligId = ligNames{ iLig } ;

            % Go through ligament properties
            for iProp = 1 : numProps
                Idx.lig.( ligNames{ iLig } ).( ligProp{iProp} ) = ...
                    find( strncmpi( frcLabels, ligId, length( ligId ) ) & endsWith( frcLabels , ligProp{iProp} ) ) ;
            end

        end


        % Contact Forces
        %----------------
        cntCompNames = Params.contactCompartmentNames ;
        numCntComp = length( cntCompNames ) ;
        contactLoadNames = Params.contactForces ;
        numContLoads = length( contactLoadNames ) ;

        % Loop through compartment names
        for iCntComp = 1 : numCntComp

            % Specify which cartilage surface to find data from
            switch cntCompNames{ iCntComp }
                case 'tf_contact'
                    cart = 'tibia_cartilage' ;
                case 'pf_contact'
                    cart = 'patella_cartilage' ;
            end

            % Loop through contact loads
            for iCntLoad = 1 : numContLoads
                Idx.cnt.( cntCompNames{ iCntComp } ).( contactLoadNames{ iCntLoad } ) = ...
                    find( strncmpi( frcLabels, cntCompNames{ iCntComp }, length( cntCompNames{ iCntComp } ) ) & ... % contains compartment name
                    contains( frcLabels, contactLoadNames{ iCntLoad } ) & ... % contains contact load name
                    contains( frcLabels , cart ) ) ; % contains cartilage surface
            end

        end

        % External Loads
        %---------------

        % F = forces, p = load locations, T = torques
        extLoads = { 'load_Fx' , 'load_Fy' , 'load_Fz' , ...
            'load_px' , 'load_py' , 'load_pz' , ...
            'load_Tx' , 'load_Ty' , 'load_Tz' } ;
        numExtLoad = length( extLoads ) ;

        % Loop through external loads
        for iExtLoad = 1 : numExtLoad
            Idx.extLoad.( extLoads{iExtLoad} ) = ...
                find( endsWith( frcLabels , extLoads{iExtLoad} ) ) ;
        end


end

end


function [data, labels, header] = read_opensim_mot(file)
%%=========================================================================
%READ_OPENSIM_MOT
%--------------------------------------------------------------------------
%Author(s): Colin Smith
%Date: 5/14/2018

%The kneemos_matlab toolkit is a collection of code for developing and
%analyzing musculoskeletal simulations in SIMM and OpenSIM. The developers
%are based at the University of Wisconsin-Madison and ETH Zurich. Please
%see the README.md file for more details. It is your responsibility to
%ensure this code works correctly for your use cases.
%
%Licensed under the Apache License, Version 2.0 (the "License"); you may
%not use this file except in compliance with the License. You may obtain a
%copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
%Unless required by applicable law or agreed to in writing, software
%distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
%WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
%License for the specific language governing permissions and limitations
%under the License.
%
%--------------------------------------------------------------------------
%file : str
%
%
%data : [nFrames x nLabels] matrix
%
%
%labels : {1 x nLabels} cell of strings
%
%
%header : struct
%
%
%==========================================================================

if ~exist('file','var')
    [infile, inpath]=uigetfile('*.mot','Select input file');
    file=[inpath infile];
end

fid=fopen(file,'r');

if fid <0
    mot=[];labels=[];
    disp('File Not Found:\n');
    disp([file '\n']);
    return
end

disp(['Loading file...' file] );


%read the file name line
header.filename=strtrim(fgetl(fid));


% Read Header
line = fgetl(fid);
while ~strncmpi(line,'endheader',length('endheader'))


    if (line == -1)
        disp('ERROR: Reached EOF before "endheader"')
        return
    end
    line_space = strrep(line,'=',' ');
    split_line = strsplit(line_space);

    if (length(split_line)==2)
        var = split_line{1};
        value = split_line{2};

        if strcmpi(var,'version')
            header.version = str2double(value);
        elseif strcmpi(var,'nRows') || strcmpi(var,'datarows')
            nr = str2double(value);
            header.nRows = nr;
        elseif strcmpi(var,'nColumns') || strcmpi(var,'datacolumns')
            nc = str2double(value);
            header.nColumns = nc;
        elseif strcmpi(var,'indegrees')
            header.indegrees = strtrim(value);
        end
    end

    line = fgetl(fid);
end


%Load Column Headers
line=fgetl(fid);

labels=cell(nc,1);

j=1;
jl=length(line);
for i=1:nc
    name=sscanf(line(j:jl),'%s',1);
    ii = findstr(line(j:jl), name);
    j=j+ii(1)+length(name);
    labels(i,1)=cellstr(name);
end

% Now load the data
data=zeros(nr,nc);
i=0;
while ((feof(fid)==0)&&(i<nr))
    i=i+1;
    line=fgetl(fid);
    try
        data(i,:)=sscanf(line,'%f');
    catch
        line = strrep( line, 'nan(ind)', '0' ) ;
        data(i,:)=sscanf(line,'%f');
        disp( ' ' )
        disp( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
        disp( '!Warning! Data contained nan-values that were set to 0!' )
        disp( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
        disp( ' ' )
    end
end
if (i<nr)
    disp(['Number of rows (',num2str(i),') is less than that specified in header (',num2str(nr),')']);
    data=data(1:i,:);
end
fclose(fid);



if (nargout>1)
    % return all data in a single array
    mot=data;
elseif (nargout==1)
    % return all information in a single structure
    mot.data=data;
    mot.hdr=labels;
end


    function [t,q]=load_exp(file);
        global NQ NM;
        rawdata=load_motionfile(file);
        [t,data]=extractcolumns(rawdata,1);
        [q,data]=extractcolumns(data,NQ);
    end


    function [x,outdata]=extractcolumns(data,nc);
        x=data(:,1:nc);
        [m,n]=size(data);
        outdata=data(:,(nc+1):n);
    end

end
