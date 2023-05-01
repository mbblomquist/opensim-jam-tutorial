function doneMessage = createStochMalalignMdl( Params )
%
% -------------------------------------------------------------------------
% this script is used to
% run monte carlo simulations on vv alignment of femoral component,
% ie alignment of femoral component, and vv alignment of tibial component
%
% Overview:
%
% 1. Import surface model
% 2. Define stochastic alignment parameters
% 3. Create malaligned component files
%
%
% WRITTEN BY JOSH ROTH
% v1 (JDR): Initial Release
% v2 (MBB): Compatible with OpenSim
% v3 (MBB): Updated to fit into OpenSim_JAM pipeline
% -------------------------------------------------------------------------
%% Parse inputs
%==============

numMdls = Params.numModels ; % number of models to create

femName = Params.femImplant ; % femur implant name
tibName = Params.tibImplant ; % tibia implant name
distType = Params.distType ; % distribution type
if isfield( Params , 'femRot' )
    femRot = Params.femRot ; % distribution parameters for femur
else
    femRot = '' ;
end
if isfield( Params , 'tibRot' )
    tibRot = Params.tibRot ; % distribution parameters for tibia
else
    tibRot = '' ;
end


%% Load articular surface models
%===============================

% Compute points, triangulation, and normals
[ Fem.pts, Fem.tri, Fem.norm ] = ...
    STL_Import( fullfile( Params.implantDir, femName ), 1 ) ;
[ Tib.pts, Tib.tri, Tib.norm ] = ...
    STL_Import( fullfile( Params.implantDir, tibName ), 1 ) ;

% Insert into Params structure for malrotateArticularSurfaces
Params.Fem = Fem ;
Params.Tib = Tib ;

%% Define stochastic alignment parameters
%==============================

%!!!TODO!!! add all DOFs

switch distType
    case 'normal'
        % Femural component
        % Fem vv
        if isfield( femRot , 'vv' )
            mu = femRot.vv( 1 );
            sd = femRot.vv( 2 );
            vvFemMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            vvFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Fem ie
        if isfield( femRot , 'ie' )
            mu = femRot.ie( 1 ) ;
            sd = femRot.ie( 2 ) ;
            ieFemMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            ieFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Tibial component
        % Tib vv
        if isfield( tibRot , 'vv' )
            mu = tibRot.vv( 1 );
            sd = tibRot.vv( 2 );
            vvTibMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            vvTibMalrot = zeros( numMdls, 1 ) ;
        end
        % Tib ie
        if isfield( tibRot , 'ie' )
            mu = tibRot.ie( 1 );
            sd = tibRot.ie( 2 );
            ieTibMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            ieTibMalrot = zeros( numMdls, 1 ) ;
        end
        % Tib fe
        if isfield( tibRot , 'fe' )
            mu = tibRot.fe( 1 );
            sd = tibRot.fe( 2 );
            feTibMalrot = normrnd( mu , sd , numMdls , 1 ) ;
        else
            feTibMalrot = zeros( numMdls, 1 ) ;
        end
    case 'uniform'
        % Femoral component
        % Fem vv
        if isfield( femRot , 'vv' )
            lowLim = femRot.vv( 1 );
            upLim = femRot.vv( 2 );
            vvFemMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            vvFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Fem ie
        if isfield( femRot , 'ie' )
            lowLim = femRot.ie( 1 );
            upLim = femRot.ie( 2 );
            ieFemMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            ieFemMalrot = zeros( numMdls, 1 ) ;
        end
        % Tibial component
        % Tib vv
        if isfield( tibRot , 'vv' )
            lowLim = tibRot.vv( 1 );
            upLim = tibRot.vv( 2 );
            vvTibMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            vvTibMalrot = zeros( numMdls, 1 ) ;
        end
        % Tib ie
        if isfield( tibRot , 'ie' )
            lowLim = tibRot.ie( 1 );
            upLim = tibRot.ie( 2 );
            ieTibMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            ieTibMalrot = zeros( numMdls, 1 ) ;
        end
        % Tib fe
        if isfield( tibRot , 'fe' )
            lowLim = tibRot.fe( 1 );
            upLim = tibRot.fe( 2 );
            feTibMalrot = lowLim + ( upLim - lowLim ) .* rand( numMdls , 1 ) ;
        else
            feTibMalrot = zeros( numMdls, 1 ) ;
        end
end

%% Create malaligned component files
%===================================

%Save MAT-file containing all simulation parameters
%--------------------------------------------------
vvFem = round( vvFemMalrot , 2 ) ;
ieFem = round( ieFemMalrot , 2 ) ;
vvTib = round( vvTibMalrot , 2 ) ;
ieTib = round( ieTibMalrot , 2 ) ;
feTib = round( feTibMalrot , 2 ) ;

% Store into Params struct for malrotateArticularSurfaces.m
Params.vvFemMalrot = vvFemMalrot ;
Params.ieFemMalrot = ieFemMalrot ;
Params.vvTibMalrot = vvTibMalrot ;
Params.ieTibMalrot = ieTibMalrot ;
Params.feTibMalrot = feTibMalrot ;

mdlNums = 1 : numMdls ;

malrotLog.hdr = { 'Sim Num', 'vvFem', 'ieFem', 'vvTib', 'ieTib' , 'feTib' };
malrotLog.data = [ mdlNums' , vvFem , ieFem , vvTib , ieTib , feTib ];

% Save stochastic values
switch Params.localOrHT
    case 'local'
        outDir = fullfile( Params.baseOutDir , 'stochModels' ) ;
    case 'HT'
        outDir = Params.baseOutDir ; % Specify output directory
end
save( fullfile( outDir , 'malAlignImplantData.mat' ) , 'malrotLog' )

%Generate stochastic malrotation models
%--------------------------------------
for iMdl = 1 : numMdls
    [ ~ , ~ , ~ ] = malrotateArticularSurfaces( Params , iMdl ) ;
end

doneMessage = 'Stoch articular surfaces created!';

end

%% ========================== NESTED FUNCTIONS ===========================
% ========================================================================

%% malrotateArticularSurfaces
function [done,femMalrot,tibMalrot] = malrotateArticularSurfaces( Params , iSim )
% -------------------------------------------------------------------------
% this script is used to vary positions of articular surfaces
% 
% Overview:
% 
% 1. Load STL file of components
% 2. Translate point clouds to desired origin for transformation
% 3. Transform point clouds
% 4. Translate point clouds back to original origin
% 5. Write new STL file of components
% 
% 
% WRITTEN BY JOSH ROTH
% v1 (JDR): Initial Release
% v2 (MBB): Compatible with OpenSim
% v3 (MBB): Updated to fit into OpenSim_JAM pipeline
% v4 (MBB): Add in ability to change tibial slope
% -------------------------------------------------------------------------

%==============
% Parse inputs
%==============

femName = Params.femImplant ;
tibName = Params.tibImplant ;
fem = Params.Fem ; % femur properties
tib = Params.Tib ; % tibia properties
vvFemMalrot = Params.vvFemMalrot ;
ieFemMalrot = Params.ieFemMalrot ;
vvTibMalrot = Params.vvTibMalrot ;
ieTibMalrot = Params.ieTibMalrot ;
feTibMalrot = Params.feTibMalrot ;

femPts = fem.pts;
femTri = fem.tri;
femNorm = fem.norm;

tibPts = tib.pts;
tibTri = tib.tri;
tibNorm = tib.norm;

% Plot points for troubleshooting help
plotTroubleshooting = 0 ; % "1" to plot components at all steps and "0" to not plot anything
    
% APPLY OFFSET FOR TRANSFORMATION
%   translate *.stl points such that origin for tranformation is located at
%   - the most medial or lateral point for V-V rotation
%   - the A-P and M-L center for I-E rotation
%   - the most anterior point for F-E rotation

%--------------------------------------------------------------------------
% Define offsets for femoral articular surface based on axis of malrotation
%--------------------------------------------------------------------------

% =============================
% NEED TO FIX FEMORAL ROTATIONS
% =============================

if vvFemMalrot(iSim) > 0
    mlAxisFem = 'LatDist' ;
elseif vvFemMalrot(iSim) < 0
    mlAxisFem = 'MedDist' ;
else
    mlAxisFem = 'Center' ;
end

pdOffsetFem = min( femPts(:,2) ) ; % find most distal point

switch mlAxisFem
    case 'LatDist'
        latCompartmentIndicesFem = femPts(:,3)<0 ;
        latCompartmentFem = femPts(latCompartmentIndicesFem,:);
        [ ~ , latMostDistIndexFem ] = min(latCompartmentFem(:,2));
        mlOffsetFem = latCompartmentFem(latMostDistIndexFem,3); % ml offset = ml value of most distal point
        %     mlOffsetFem = (max(fem(:,3)) - min(fem(:,3)))/4; %rotation about center of lateral condyle
    case 'Center'
        mlOffsetFem = -(max(femPts(:,3)) + min(femPts(:,3)))/2; %rotation about center of component
    case 'MedDist'
        medCompartmentIndicesFem = femPts(:,3)>0 ;
        medCompartmentFem = femPts(medCompartmentIndicesFem,:);
        [ ~ , medMostDistIndexFem ] = min(medCompartmentFem(:,2));
        mlOffsetFem = medCompartmentFem(medMostDistIndexFem,3); % ml offset = ml value of most distal point
        %     mlOffsetFem = -(max(fem(:,3)) - min(fem(:,3)))/4; %rotation about center of medial condyle
end

% apOffsetFem = mean(femPts(:,1)); for rotating about AP center of femoral component
apOffsetFem = -min(femPts(:,1));

% Origin now max posterior, max distal, medial/laterally at max distal point
femCentered = [femPts(:,1)+apOffsetFem,femPts(:,2)-pdOffsetFem,femPts(:,3)+mlOffsetFem];

%-------------------------------------------------------------------------
% Define offsets for tibial articular surface based on axis of malrotation
%-------------------------------------------------------------------------

if vvTibMalrot(iSim) ~= 0

    % Find most medial or lateral point
    if vvTibMalrot(iSim) < 0
        [ ~ , rotPt ] = max( tibPts(:,3) );
    elseif vvTibMalrot(iSim) > 0
        [ ~ , rotPt ] = min( tibPts(:,3) ) ;        
    end

    % Extract A-P, P-D, and M-L coordinates of medial or lateral point
    apOffsetTib = -tibPts(rotPt,1) ;
    pdOffsetTib = -tibPts(rotPt,2) ;
    mlOffsetTib = -tibPts(rotPt,3) ;

elseif ieTibMalrot ~= 0

    midAP = mean( [ max( tibPts(:,1) ) , min( tibPts(:,1) ) ] ) ;
    midML = mean( [ max( tibPts(:,3) ) , min( tibPts(:,3) ) ] ) ;

    % Move to be midpoint of A-P and M-L (P-D doesnt matter)
    apOffsetTib = -midAP ;
    pdOffsetTib = 0 ;
    mlOffsetTib = -midML ;

elseif feTibMalrot ~= 0

    % Find most anterior point
    [ ~ , rotPt ] = max( tibPts(:,1) );

    % Extract A-P, P-D, and M-L coordinates of anterior point
    apOffsetTib = -tibPts(rotPt,1) ;
    pdOffsetTib = -tibPts(rotPt,2) ;
    mlOffsetTib = -tibPts(rotPt,3) ;

end

% Translate point
tibCentered = [ tibPts(:,1) + apOffsetTib , tibPts(:,2) + pdOffsetTib , tibPts(:,3) + mlOffsetTib ] ;

% TRANSFORM ARTICULAR SURFACES USING RANDOMLY GENERATED TRANSFORMATION MATRIX
%X = A-P, Y = P-D, Z = M-L
%
%Rx = [ 1   0   0   0;
%       0   c   -s  0;
%       0   s   c   0;
%       0   0   0   1];
%
%Ry = [ c   0   s   0;
%       0   1   0   0;
%       -s  0   c   0;
%       0   0   0   1];
%
%Rz = [ c   -s  0   0;
%       s   c   0   0;
%       0   0   1   0;
%       0   0   0   1];

feFem   = deg2rad(0); %flexion-extension rotation in rad
vvFem   = deg2rad(vvFemMalrot(iSim)); %varus-valgus rotation in rad +Valgus -Varus
ieFem   = deg2rad(ieFemMalrot(iSim)); %internal-external rotation in rad +External -Internal

%transformation applied using Cardan sequence fe-vv-ie
transFem = [cos(vvFem)*cos(ieFem)                                       -cos(vvFem)*sin(ieFem)                                  sin(vvFem)              ;
            cos(feFem)*sin(ieFem) + cos(ieFem)*sin(feFem)*sin(vvFem)    cos(feFem)*cos(ieFem)-sin(feFem)*sin(vvFem)*sin(ieFem)  -cos(vvFem)*sin(feFem)  ;
            sin(feFem)*sin(ieFem) - cos(feFem)*cos(ieFem)*sin(vvFem)    cos(ieFem)*sin(feFem)+cos(feFem)*sin(vvFem)*sin(ieFem)  cos(feFem)*cos(vvFem)   ];
        
femCenteredTrans = femCentered * transFem;
femNormTrans = femNorm * transFem;

vv2Tib   = deg2rad(vvTibMalrot(iSim)); %flexion-extension rotation in rad
ie2Tib   = deg2rad(ieTibMalrot(iSim)); %varus-valgus rotation in rad
fe2Tib   = deg2rad(feTibMalrot(iSim)); %internal-external rotation in rad

%transformation applied using Cardan sequence fe-vv-ie
transTib = [cos(ie2Tib)*cos(fe2Tib)                                       -cos(ie2Tib)*sin(fe2Tib)                                  sin(ie2Tib)              ;
            cos(vv2Tib)*sin(fe2Tib) + cos(fe2Tib)*sin(vv2Tib)*sin(ie2Tib)    cos(vv2Tib)*cos(fe2Tib)-sin(vv2Tib)*sin(ie2Tib)*sin(fe2Tib)  -cos(ie2Tib)*sin(vv2Tib)  ;
            sin(vv2Tib)*sin(fe2Tib) - cos(vv2Tib)*cos(fe2Tib)*sin(ie2Tib)    cos(fe2Tib)*sin(vv2Tib)+cos(vv2Tib)*sin(ie2Tib)*sin(fe2Tib)  cos(vv2Tib)*cos(ie2Tib)   ];
        
tibCenteredTrans = tibCentered * transTib;
tibNormTrans = tibNorm * transTib;

% ALIGN BACK TO IMPORTED COORDINATE SYSTEM

femPtsNew = [femCenteredTrans(:,1) - apOffsetFem,femCenteredTrans(:,2) + pdOffsetFem,femCenteredTrans(:,3) - mlOffsetFem];
tibPtsNew = [tibCenteredTrans(:,1) - apOffsetTib,tibCenteredTrans(:,2) - pdOffsetTib,tibCenteredTrans(:,3) - mlOffsetTib];

if isequal( plotTroubleshooting , 1 )
    figure() ; hold on ;
    plot3(tibPts(:,1),tibPts(:,2), tibPts(:,3),'g.')
    plot3(tibCentered(:,1),tibCentered(:,2), tibCentered(:,3),'c.')
    plot3(tibCenteredTrans(:,1),tibCenteredTrans(:,2), tibCenteredTrans(:,3),'b.')
    plot3(tibPtsNew(:,1), tibPtsNew(:,2), tibPtsNew(:,3), 'k.')
    xlabel('x'); ylabel('y') ; zlabel('z') ;
    axis equal tight
    legend( 'Old' , 'Transform 1' , 'Transform 2' , 'New' )
end

%assemble output structures
%--------------------------
femMalrot.pts = femPtsNew;
femMalrot.norm = femNormTrans;
femMalrot.tri = femTri;

tibMalrot.pts = tibPtsNew;
tibMalrot.norm = tibNormTrans;
tibMalrot.tri = tibTri;

% WRITE TRANSFORMED ART SURFACE FILES

femTR = triangulation( femMalrot.tri , femMalrot.pts ) ;
tibTR = triangulation( tibMalrot.tri , tibMalrot.pts ) ;

% Save stoch implant STLs
switch Params.localOrHT
    case 'local'
        newFemName = [ femName(1:end-4) , '_' num2str( iSim ) , '.stl' ] ;
        newTibName = [ tibName(1:end-4) , '_' num2str( iSim ) , '.stl' ] ;

        outDir = fullfile( Params.baseOutDir , 'stochModels' ) ;
        stlwrite( femTR , fullfile( outDir , 'Geometry' , newFemName ) ) ;
        stlwrite( tibTR , fullfile( outDir , 'Geometry' , newTibName ) ) ;
    case 'HT'
        outDir = fullfile( Params.baseOutDir , 'input' , num2str(iSim-1) )  ;
        stlwrite( femTR , fullfile( outDir , femName ) ) ;
        stlwrite( tibTR , fullfile( outDir , tibName ) ) ;
end

done = iSim;
end

%% STL_Import

function varargout=STL_Import(filename,mode)
% STL_Import is a tool designed to import into MATLAB both binary and ASCII STL files.
% 
% This scprit is mainly a collage betwwen file axchange fileid 22409 and 3642, plus
% some other features that can be considered new on FEX.
% 
% SYNOPSIS:
% 
% 
% %mode 1 (default)
% [p,t,tnorm]=STL_Import(filename,mode)
% 
% %mode 2
% [v,tnorm])=STL_Import(filename,mode)
% 
% 
% INPUT:
% 
% filename: string representing the name fo the file
% 
% mode:
% 
% 
%  mode=1 (if omitted is automatically set to one)
%     
%   set the the output to:
%     
%     output=[p,t,tnorm]
%     
%     where        
%     
%     p=points (unique) of the model nx3 array

%     t=triangles indexes of the model

%     tnorm= normals of triangles
%     

%  mode=2
%       
%   set the the output to:
%     
%     output=[v,tnorm]
%     
%     where        
%     
%     v=  vertex of the model(not unique points) of the model nx3 array. Each
%         trhee points we have a triagnle in consecutive order.

%     tnorm= normals of triangles
% 
% EXAMPLES:    
%    
%    [p,t,tnorm]=STL_Import('link1.stl',1); 
%    [pv,tnorm]=STL_Import('link1.stl',2);
% 
%
% Visit:
% 
%  http://giaccariluigi.altervista.org/blog/
%  
% Author: Giaccari Luigi  (giaccariluigi@msn.com)



if nargin<2
    mode=1;%default value
end


if ~(mode==1 || mode==2)
    error('invalid mode')
end

if nargout<3 && mode==1
    error('invalid input number /mode setting')
end
if nargout>2 && mode==2
    error('invalid input number /mode setting')
end


%open file
fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
if fid == -1
    error('File could not be opened, check name or path.')
end


 M = fread(fid,inf,'uint8=>uint8');
    fclose(fid);    
    
if( isbinary(M) )
     [v,tnorm]=ImportSTL_binary(M);
    
else
    clear M;    
    [v,tnorm]=ImportSTL_ASCII(filename);
   
end

clear M

varargout = cell(1,nargout);
switch mode
    case 1
        [p,t]=fv2pt(v,length(v)/3);%gets points and triangles
        
        varargout{1} = p;
        varargout{2} = t;
        varargout{3} = tnorm;
    case 2
        varargout{1} = v;
        varargout{2} = tnorm;
end
end



function [v,tnorm]=ImportSTL_ASCII(filename)

%counting the number of vertex
vnum=0;
fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
while feof(fid) == 0                    % test for end of file, if not then do stuff
    tline = fgetl(fid);                 % reads a line of data from file.
    fword = sscanf(tline, '%s ');       % make the line a character string
    if strncmpi(fword, 'v',1)   % Checking if a "V"ertex line, as "V" is 1st char.
        vnum = vnum + 1;                % If a V we count the # of V's
    end
end
numt=ceil(vnum/3);%triangles number equals vertex number/3

tnorm=zeros(numt,3);%preallocate for normals
v=zeros(vnum,3);%not unique vertex

c=0;%vertex counter
fnum=0;
fid=fopen(filename, 'r'); %REOpen the file
while feof(fid) == 0                    % test for end of file, if not then do stuff
    tline = fgetl(fid);                 % reads a line of data from file.
    fword = sscanf(tline, '%s ');       % make the line a character string
    
    % Check vertex
    if strncmpi(fword, 'v',1)    % Checking if a "V"ertex line, as "V" is 1st char.
        c = c + 1;                % If a V we count the # of V's
        v(c,:) = sscanf(tline, '%*s %f %f %f'); % & if a V, get the XYZ data of it.
        
        % Check facet normal
    elseif strncmpi(fword, 'f',1)   % Checking if a "V"ertex line, as "V" is 1st char.
        fnum =fnum + 1;                % If a V we count the # of V's
        tnorm(fnum,:) = sscanf(tline, '%*s %*s %f %f %f'); % & if a V, get the XYZ data of it.
        
        % %% Check for color
        %     elseif strncmpi(fword, 'c',1) ;    % Checking if a "C"olor line, as "C" is 1st char.
        %         VColor = sscanf(tline, '%*s %f %f %f'); % & if a C, get the RGB color data of the face.
        %    % Keep this color, until the next color is used.
        
    end
    
end
end


function [p,t]=fv2pt(v,fnum)

%gets points and triangle indexes given vertex and facet number

c=size(v,1);

%triangles with vertex id data
t=zeros(3,fnum);
t(:)=1:c;


%now we have to keep unique points fro vertex
[p,~,j]=unique(v,'rows'); %now v=p(j) p(i)=v;
t(:)=j(t(:));
t=t';

end

% 


function tf = isbinary(A)
% ISBINARY determines if an STL file is binary or ASCII.

    % Look for the string 'endsolid' near the end of the file
    if isempty(A) || length(A) < 16
        error('MATLAB:stlread:incorrectFormat', ...
              'File does not appear to be an ASCII or binary STL file.');
    end
    
    % Read final 16 characters of M
    i2  = length(A);
    i1  = i2 - 100;%100 empirical value
    str = char( A(i1:i2)' );
    
    k = strfind(lower(str), 'endsolid');
    if ~isempty(k)
        tf = false; % ASCII
    else
        tf = true;  % Binary
    end
end


function [V,N]=ImportSTL_binary(M)

 
    
    if length(M) < 84
        error('MATLAB:stlread:incorrectFormat', ...
              'Incomplete header information in binary STL file.');
    end
    
    % Bytes 81-84 are an unsigned 32-bit integer specifying the number of faces
    % that follow.
    numFaces = typecast(M(81:84),'uint32');
    %numFaces = double(numFaces);
    if numFaces == 0
        warning('MATLAB:stlread:nodata','No data in STL file.');
        return
    end
    
    T = M(85:end);

    V = NaN(3*numFaces,3);
    N = NaN(numFaces,3);
    
    numRead = 0;
    while numRead < numFaces
        % Each facet is 50 bytes
        %  - Three single precision values specifying the face normal vector
        %  - Three single precision values specifying the first vertex (XYZ)
        %  - Three single precision values specifying the second vertex (XYZ)
        %  - Three single precision values specifying the third vertex (XYZ)
        %  - Two unused bytes
        i1    = 50 * numRead + 1;
        i2    = i1 + 50 - 1;
        facet = T(i1:i2)';
        
        n  = typecast(facet(1:12),'single');
        v1 = typecast(facet(13:24),'single');
        v2 = typecast(facet(25:36),'single');
        v3 = typecast(facet(37:48),'single');
        
        n = double(n);
        v = double([v1; v2; v3]);
        
        % Figure out where to fit these new vertices, and the face, in the
        % larger F and V collections.        
        fInd  = numRead + 1;        
        vInd1 = 3 * (fInd - 1) + 1;
        vInd2 = vInd1 + 3 - 1;
        
        V(vInd1:vInd2,:) = v;
        N(fInd,:)        = n;
        
        numRead = numRead + 1;
    end
    
end