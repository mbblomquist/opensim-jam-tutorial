%% Create Stoch Models
%==========================================================================
% Purpose:  Generate stochastic OpenSim knee models for use with JAM.
%
% Inputs: Params (parameter variables)
%
% Outputs: StochMdlParams (stochastic values and model values for
% ligaments and properties to change
%
% Written By: Joshua D. Roth (2020-08-17)
% Revised By: Matthew B. Blomquist (2022-11-27)
%
% See Also: <Enter other functions used>
%--------------------------------------------------------------------------
% Revision history:
% -----------------
% v1    2020-08-17(JDR)     inital release
% v2    2022-11-27(MBB)     updated to be a function
% v3    2023-03-28(MBB)     Update code for local or HT (MBB)
%
%==========================================================================
function StochMdlParams = createStochModels( Params )

%% Extract Params structure
%--------------------------

ligProps = Params.ligPropsToChange ;
probDistType = Params.probDistType ;
probDistRef = Params.probDistRef ;
probDistParams = Params.probDistParams ;

if isequal( Params.changeLigCoords, 1)

    ligCoord = Params.ligNamesCoord ;
    ligPropCoord = Params.ligPropertiesCoord ;
    probDistTypeCoord = Params.probDistTypeCoord ;
    probDistRefCoord = Params.probDistRefCoord ;
    probDistParamsCoord = Params.probDistParamsCoord ;

end

% Use nested function to create ligNames
ligNames = findLigNames( Params ) ;

%% Generate stochastic model parameters
%--------------------------------------

% Loop through each model, ligament, and property to create stochastic
% model values
for iMdl = 1 : Params.numModels
    for iLig = 1 : length( ligNames )
        tempLig = ligNames{ iLig } ;
        if isequal( Params.changeIdvStrandsProp , 1 )
            numStrands = findNumStrands( tempLig ) ;
        else
            numStrands = 1 ;
        end
        for iStrand = 1 : numStrands
            for iProp = 1 : length( ligProps )

                % Switch depending on if distribution is uniform or normal
                switch probDistType{ iProp }
                    case 'uniform'

                        StochMdlParams.stoch.( tempLig ).( ligProps{ iProp } )( iMdl, iStrand ) = ...
                            probDistParams{ iProp }( 1 ) + ...
                            rand() * ( probDistParams{ iProp }( 2 ) - probDistParams{ iProp }( 1 ) ) ;

                    case 'normal'

                        StochMdlParams.stoch.( ligNames{ iLig } ).( ligProps{ iProp } )( iMdl, iStrand ) = ...
                            normrnd( probDistParams{ iProp }( 1 ) , probDistParams{ iProp }( 2 ) )  ;

                end % switch probDistType
            end % for iProp
        end % for iStrand
    end % for iLig
end % for iMdl

if isequal(Params.changeLigCoords, 1)
    for iMdl = 1 : Params.numModels
        for iLig = 1 : length( ligCoord )
            for iProp = 1 : length( ligPropCoord )
                switch probDistTypeCoord{ iProp }
                    case 'uniform'

                        StochMdlParams.stoch.( ligCoord{ iLig } ).( ligPropCoord{ iProp } )( iMdl, 1 ) = ...
                            probDistParamsCoord{ iProp }( 1 ) + ...
                            rand() * ( probDistParamsCoord{ iProp }( 2 ) - probDistParamsCoord{ iProp }( 1 ) ) ;

                end %switch probDistType and probDistRef
            end %for iProp
        end %for iLig
    end %for iMdl
end

%% Setup OpenSim API
%--------------------

import org.opensim.modeling.*

%% Identify ligament(s) in model to change
%-------------------------

% Load base model
model = Model( Params.baseMdlFile ) ;

% Identify ligament(s) in model
ligStrandCnt = zeros( length( ligNames ) , 1 ) ;
for iFrc = 0 : model.getForceSet.getSize( ) - 1
    force = model.getForceSet.get( iFrc ) ; % get force

    for iLig = 1 : length( ligNames )
        tempLigName = ligNames{ iLig } ;
        % If force matches ligament name, add to structure
        if contains( char( force.getName( ) ) , tempLigName )
            ligChar = char( force.getName( ) ) ;

            if ligChar(1) == tempLigName(1) % Need this for PFL and mPFL/lPFL
                ligStrandCnt( iLig ) = ligStrandCnt( iLig ) + 1 ;
                MdlLigNames.( tempLigName ){ ligStrandCnt( iLig ) } = force.getName(  ) ;
            end % if ligChar(1)

        end % if contains
    end % for iLig
    clear force
end % for iFrc

%% Loop through models
%---------------------

% Go through each model to set new parameters
for iMdl = 1 : Params.numModels
    stochModel = Model( model ) ; % Creates a copy of template model

    % % Loop through ligament and strand to get force
    % for iLig = 1 : length( ligNames )
    %     for iStrand = 1 : ligStrandCnt( iLig )
    %         force = stochModel.getForceSet.get( MdlLigNames.( ligNames{ iLig } ){ iStrand } ) ;
    % 
    %         initLen = 1 ;
    % 
    %         % Loop through properties to change values
    %         for iProp = 1 : length( ligPropCoord )
    %             propName = ligPropCoord{ iProp } ;
    % 
    %             switch probDistRefCoord{ iProp }
    %                 case 'relativeAbs'
    %                     abs_prop = force.getPropertyByName( 'GeometryPath' ) ;
    %                     lig_obj = abs_prop.updValueAsObject() ;
    %                     geomPath = GeometryPath.safeDownCast( lig_obj ) ;
    %                     attachPts = geomPath.updPathPointSet() ;
    % 
    %                     if isequal( initLen , 1 )
    %                         femAttachPts = PathPoint.safeDownCast( attachPts.get( 0 ) ) ;
    %                         femAttachPtsArray = osimVec3ToArray( femAttachPts.get_location() ) ;
    %                         tibAttachPts = PathPoint.safeDownCast( attachPts.get( 1 ) ) ;
    %                         tibAttachPtsArray = osimVec3ToArray( tibAttachPts.get_location() ) ;
    %                         initialLength = norm( femAttachPtsArray - tibAttachPtsArray ) ;
    %                         initLen = 0 ;
    %                     end
    % 
    %                     if contains( propName , 'Fem' ) % femur
    %                         femAttachPts = PathPoint.safeDownCast( attachPts.get( 0 ) ) ;
    %                         femAttachPtsArray = osimVec3ToArray( femAttachPts.get_location() ) ;
    %                         if contains( propName, 'x' ) % x coord
    %                             newCoord = femAttachPtsArray(1) + ...
    %                                 StochMdlParams.stoch.( ligCoord{ iLig } ).( propName )( iMdl, 1 ) ;
    %                             StochMdlParams.mdl.( ligNames{ iLig } ).( propName )( iMdl, iStrand ) = newCoord ;
    %                             femAttachPts.set_location( Vec3( newCoord, femAttachPtsArray(2), femAttachPtsArray(3) ) )
    %                         elseif contains( propName , 'y' ) % y coord
    %                             newCoord = femAttachPtsArray(2) + ...
    %                                 StochMdlParams.stoch.( ligCoord{ iLig } ).( propName )( iMdl, 1 ) ;
    %                             StochMdlParams.mdl.( ligNames{ iLig } ).( propName )( iMdl, iStrand ) = newCoord ;
    %                             femAttachPts.set_location( Vec3( femAttachPtsArray(1), newCoord, femAttachPtsArray(3) ) )
    %                         end
    %                     end
    % 
    %             end % switch PropDistRefCoord
    %         end % for iProp
    % 
    %         femAttachPts = PathPoint.safeDownCast( attachPts.get( 0 ) ) ;
    %         femAttachPtsArray = osimVec3ToArray( femAttachPts.get_location() ) ;
    %         tibAttachPts = PathPoint.safeDownCast( attachPts.get( 1 ) ) ;
    %         tibAttachPtsArray = osimVec3ToArray( tibAttachPts.get_location() ) ;
    %         currentLength = norm( femAttachPtsArray - tibAttachPtsArray ) ;
    % 
    %         % Loop through properties to change values
    %         for iProp = 1 : length( ligProps )
    %             propName = ligProps{ iProp } ;
    % 
    %             ph = PropertyHelper();
    %             currentValue = ph.getValueDouble( force.getPropertyByName( propName ) ) ;
    %             newValue = (currentLength/initialLength) * currentValue * ...
    %                 ( 1 + StochMdlParams.stoch.( ligNames{ iLig } ).( propName )( iMdl, iStrand ) ) ;
    %             StochMdlParams.mdl.( ligNames{ iLig } ).( propName )( iMdl, iStrand ) = newValue ;
    %             ph.setValueDouble( newValue, force.getPropertyByName( propName ) ) ;
    % 
    %         end
    % 
    %         clear force ph currentValue newValue
    % 
    %     end % for iStrand
    % end % for iLig

    % Loop through ligament and strand to get force
    for iLig = 1 : length( ligNames )
        for iStrand = 1 : ligStrandCnt( iLig )
            force = stochModel.getForceSet.get( MdlLigNames.( ligNames{ iLig } ){ iStrand } ) ;

            % Loop through properties to change values
            for iProp = 1 : length( ligProps )
                propName = ligProps{ iProp } ;

                switch probDistRef{ iProp }
                    case 'relativePercent'
                        ph = PropertyHelper();
                        currentValue = ph.getValueDouble( force.getPropertyByName( propName ) ) ;
                        newValue = currentValue * ...
                            ( 1 + StochMdlParams.stoch.( ligNames{ iLig } ).( propName )( iMdl, 1 ) ) ;
                        StochMdlParams.mdl.( ligNames{ iLig } ).( propName )( iMdl, iStrand ) = newValue ;
                        ph.setValueDouble( newValue, force.getPropertyByName( propName ) ) ;
                    case 'relativeAbs'
                        ph = PropertyHelper();
                        currentValue = ph.getValueDouble( force.getPropertyByName( propName ) ) ;
                        newValue = currentValue + ...
                            StochMdlParams.stoch.( ligNames{ iLig } ).( propName )( iMdl, 1 ) ;
                        StochMdlParams.mdl.( ligNames{ iLig } ).( propName )( iMdl, iStrand ) = newValue ;
                        ph.setValueDouble( newValue, force.getPropertyByName( propName ) ) ;
                    case 'absolute'
                        ph = PropertyHelper();
                        newValue = StochMdlParams.stoch.( ligNames{ iLig } ).( propName )( iMdl, 1 ) ;
                        StochMdlParams.mdl.( ligNames{ iLig } ).( propName )( iMdl, iStrand ) = newValue ;
                        ph.setValueDouble( newValue, force.getPropertyByName( propName ) ) ;
                end % switch probDistRef

            end % for iProp
            clear force ph currentValue newValue

        end % for iStrand
    end % for iLig

    stochModel.initSystem();

    % Set up output files for local or HT running
    switch Params.localOrHT
        case 'local'
            outDir = fullfile( Params.baseOutDir , 'stochModels' ) ;
            stochModelName = [ 'lenhart2015_stoch' , num2str(iMdl) , '.osim' ] ;
            saveFileDir = outDir ;
        case 'HT'
            outDir = fullfile( Params.baseOutDir , 'input' , num2str(iMdl-1) ) ; % Specify output directory
            stochModelName = 'lenhart2015_stoch.osim' ;
            saveFileDir = Params.baseOutDir ;
    end

    stochModelFile = fullfile( outDir , stochModelName ) ;
    stochModel.print( stochModelName ) ;
    movefile( stochModelName , stochModelFile )
    clear stochModel
    disp( [ 'Wrote model to: ', stochModelFile ] )
end

save( fullfile( saveFileDir , 'stochMdlParams.mat' ) , 'StochMdlParams' )

doneMsg = 'Done creating models' ;
disp( doneMsg )

end

%% ========================== NESTED FUNCTIONS ===========================
% ========================================================================

function ligNames = findLigNames( Params )

if contains( Params.ligNamesToChange , 'allLigs' )
    switch Params.baseMdl
        case { 'lenhart2015' , 'lenhart2015_OA' , 'lenhart2015_meniscus' }
            ligNames = { 'MCLd' , 'MCLs', 'MCLp', 'ACLpl' , 'ACLam' , ...
                'LCL', 'ITB', 'PFL', 'pCAP', 'PCLpm', 'PCLal', 'PT', ...
                'lPFL', 'mPFL' } ;
        case { 'lenhart2015_implant' }
            ligNames = { 'MCLs', 'MCLp', 'LCL', 'ITB', 'PFL', ...
                'pCAP', 'PCLpm', 'PCLal', 'PT', 'lPFL', 'mPFL' } ;
    end
else
    ligNames = Params.ligNamesToChange ;
end

end

function numStrands = findNumStrands( ligamentName )
switch ligamentName
    case 'ITB'
        numStrands = 1 ;
    case { 'LCL' }
        numStrands = 4 ;
    case { 'MCLd' , 'PCLal' , 'PCLpm' , 'MCLp' , 'PFL' }
        numStrands = 5 ;
    case { 'MCLs' , 'ACLpl' , 'ACLam' , 'PT' , 'mPFL' }
        numStrands = 6 ;
    case { 'lPFL' , 'pCAP' }
        numStrands = 8 ;
end
end