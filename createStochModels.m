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

% Use nested function to create ligNames
ligNames = findLigNames( Params ) ;

%% Generate stochastic model parameters
%--------------------------------------

% Loop through each model, ligament, and property to create stochastic
% model values
for iMdl = 1 : Params.numModels
    for iLig = 1 : length( ligNames )
        for iProp = 1 : length( ligProps )

            % Switch depending on if distribution is uniform or normal
            switch probDistType{ iProp }
                case 'uniform'

                    StochMdlParams.stoch.( ligNames{ iLig } ).( ligProps{ iProp } )( iMdl, 1 ) = ...
                        probDistParams{ iProp }( 1 ) + ...
                        rand() * ( probDistParams{ iProp }( 2 ) - probDistParams{ iProp }( 1 ) ) ;

                case 'normal'

                    StochMdlParams.stoch.( ligNames{ iLig } ).( ligProps{ iProp } )( iMdl, 1 ) = ...
                        normrnd( probDistParams{ iProp }( 1 ) , probDistParams{ iProp }( 2 ) )  ;

            end % switch probDistType
        end % for iProp
    end % for iLig
end % for iMdl

%% Setup OpenSim API
%--------------------

import org.opensim.modeling.*

%% Build stochastic models
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

% Go through each model to set new parameters
for iMdl = 1 : Params.numModels
    stochModel = Model( model ) ; % Creates a copy of template model

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
            case { 'lenhart2015' , 'lenhart2015_BCRTKA' }
                ligNames = { 'MCLd' , 'MCLs', 'MCLp', 'ACLpl' , 'ACLam' , ...
                    'LCL', 'ITB', 'PFL', 'pCAP', 'PCLpm', 'PCLal', 'PT', ...
                    'lPFL', 'mPFL' } ;
            case { 'lenhart2015_implant' , 'lenhart2015_SarahISTA_PCL' }
                ligNames = { 'MCLs', 'MCLp', 'LCL', 'ITB', 'PFL', ...
                    'pCAP', 'PCLpm', 'PCLal', 'PT', 'lPFL', 'mPFL' } ;
            case 'lenhart2015_SarahISTA_noPCL'
                ligNames = { 'MCLs', 'MCLp', ...
                    'LCL', 'ITB', 'PFL', 'pCAP', 'PT', ...
                    'lPFL', 'mPFL' } ;
            case 'lenhart2015_UKA'
                ligNames = { 'MCLd' , 'MCLs', 'MCLp', 'ACLpl' , 'ACLam' , ...
                    'LCL', 'ITB', 'PFL', 'pCAP', 'PCLpm', 'PCLal', 'PT', ...
                    'lPFL', 'mPFL' } ;
        end
    else
        ligNames = Params.ligNamesToChange ;
    end

end