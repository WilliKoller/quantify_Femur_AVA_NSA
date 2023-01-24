%-------------------------------------------------------------------------%
% Copyright (c) 2022 Koller W.                                            %
%    Author:   Willi Koller,  2022                                        %
%    email:    willi.koller@univie.ac.at                                  %
% ----------------------------------------------------------------------- %

% calcFrom3Dgeometry (c) by Willi Koller, University of Vienna
%
% calcFrom3Dgeometry is licensed under a
% Creative Commons Attribution-NonCommercial 4.0 International License.
%
% You should have received a copy of the license along with this
% work. If not, see <http://creativecommons.org/licenses/by-nc/4.0/>.

clear;
[file, path] = uigetfile('*.stl');
filenameNoExt = strrep(file, '.stl', '');
side = questdlg('Which side is this femur?', 'Side?', 'left', 'right', 'left');
if strcmp(side, 'left')
    side = 'L';
else
    side = 'R';
end

%%

triGeom_set = createTriGeomSet({filenameNoExt}, path);
triGeom_set.(['femur_' lower(side)]) = triGeom_set.(filenameNoExt);
triGeom_set = rmfield(triGeom_set, filenameNoExt);

% code from GIBOC (file "femur_guess_CS")
Femur = triGeom_set.(['femur_' lower(side)]);

connList = triGeom_set.(['femur_' lower(side)]).ConnectivityList;
points = triGeom_set.(['femur_' lower(side)]).Points;

[ V_all, ~ ] = TriInertiaPpties( Femur );
Z0 = V_all(:,1);

%% Get the central 60% of the bone -> The femur diaphysis
LengthBone = max(Femur.Points*Z0) - min(Femur.Points*Z0);
L_ratio = 0.20;
% First remove the top 20% percent
alt_top = max(Femur.Points*Z0) - L_ratio* LengthBone;
ElmtsTmp1 = find(Femur.incenter*Z0<alt_top);
TrTmp1 = TriReduceMesh( Femur, ElmtsTmp1);
TrTmp1 = TriFillPlanarHoles( TrTmp1 );
% Then remove the bottom 20% percent
alt_bottom = min(Femur.Points*Z0) + L_ratio* LengthBone;
ElmtsTmp2 = find(TrTmp1.incenter*Z0>alt_bottom);
TrTmp2 = TriReduceMesh( TrTmp1, ElmtsTmp2);
FemurDiaphysis = TriFillPlanarHoles( TrTmp2 );

%% Get the principal inertia axis of the diaphysis (potentially wrongly orientated)
[ V_all, ~ ] = TriInertiaPpties( FemurDiaphysis );
diaphyseVector = V_all(:,1);

try
    [~, JCS, ~, ~, CS] = GIBOC_femur(triGeom_set.(['femur_' lower(side)]), side, 'cylinder');
catch
    [~, JCS, ~, ~, CS] = GIBOC_femur(triGeom_set.(['femur_' lower(side)]), side, 'ellipsoids');
end

hipOrigin = JCS.(['hip_' lower(side)]).Origin';
rotMatrixToOpenSim = JCS.(['hip_' lower(side)]).V;
kneeOrigin = (JCS.(['knee_' lower(side)]).Origin' - hipOrigin) * rotMatrixToOpenSim;
fem_head_radius = CS.RadiusFH_Renault;

diaphyseVector = diaphyseVector' * rotMatrixToOpenSim;
if diaphyseVector(2) < 0
    diaphyseVector = diaphyseVector * (-1);
end

%% import the STL
[stlStruct] = import_STL(fullfile(path, file));
vertexMAT = (stlStruct.solidVertices{1, 1} - hipOrigin) * rotMatrixToOpenSim;

points = (points - hipOrigin) * rotMatrixToOpenSim;
triGeom_setOpenSim = triangulation(connList, points);

D = pdist2(vertexMAT, [0 0 0]);

% get neck axis cutoff
f = figure;
hold on;
rotate3d on;
b = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23], 'value',30, 'min',15, 'max',50);

temp = vertexMAT;
global neckAxisCutOff;
neckAxisCutOff = 30;
temp(D > neckAxisCutOff, :) = [];
global h1;
global h2;
h1 = plotNodes(temp, 'r');
temp = vertexMAT;
temp(or(D < neckAxisCutOff, D > 100), :) = [];
h2 = plotNodes(temp, 'b');

daspect([1 1 1]);
if strcmp(side, 'R')
    view([-100, -20]);
else
    view([-100, 20]);
end

b.Callback = @(es,ed) updateNeckAxisLengthPlot(es.Value, vertexMAT, D);

drawnow;
figure(f);

answer = MFquestdlg([0.4, 0.4], 'Press yes when the selection for determining the neck axis is okay.', 'Confirm Selection', 'Yes', 'No', 'Yes');
if ~strcmp(answer, 'Yes')
    return;
end
close(f);

%
neckNodes = vertexMAT;
neckNodes(D > neckAxisCutOff, :) = [];
allNeckNodes = neckNodes;

if strcmp(side, 'R')
    neckAxisGuess = [-0.3 -0.5 0.3];
else
    neckAxisGuess = [-0.3 -0.5 -0.3];
end

p2 = neckAxisGuess * (-1) * fem_head_radius;
D = pdist2(neckNodes, p2);

neckNodes(D < fem_head_radius, :) = [];

% use only ~ 5000 nodes
nodesUsedForCylinderFit = neckNodes(floor(linspace(1, size(neckNodes, 1), 5000)), :);

[~, neckAxisGuess] = lscylinder(nodesUsedForCylinderFit, [0 0 0]', neckAxisGuess', 15, 0.1, 0.1);

latestGuessNeckAxis = neckAxisGuess;
iterationCounter = 1;
while(1)
    neckNodes = allNeckNodes;
    p2 = latestGuessNeckAxis * (-1) * fem_head_radius;
    D = pdist2(neckNodes, p2');
    neckNodes(D < fem_head_radius, :) = [];
    
    % use only ~ 5000 nodes
    nodeIndizes = unique(floor(linspace(1, size(neckNodes, 1), 5000)));
    nodesUsedForCylinderFit = neckNodes(nodeIndizes, :);
    [~, neckAxis] = lscylinder(nodesUsedForCylinderFit, [0 0 0]', latestGuessNeckAxis', 15, 0.1, 0.1);
    
    % calc angle between neckAxis and latestGuessNeckAxis
    angleChange = atan2d(vecnorm(cross(neckAxis,latestGuessNeckAxis)),dot(neckAxis,latestGuessNeckAxis));
    
    if angleChange < 0.2
        break;
    else
        latestGuessNeckAxis = neckAxis;
    end

    iterationCounter = iterationCounter + 1;
    if iterationCounter > 100
        disp('fail safe - finding neck axis did not converge after 100 iterations, changes were bigger than 0.2° each iteration');
        break;
    end
end

figure;
hold on;
daspect([1 1 1]);
rotate3d on;
plotTriangLight(triGeom_setOpenSim, JCS.(['hip_' lower(side)]), 0, 0.5);
plotNodes(nodesUsedForCylinderFit, 'r');
quiver3(neckAxis(1)*(-60),neckAxis(2)*(-60),neckAxis(3)*(-60), neckAxis(1), neckAxis(2), neckAxis(3), 120, 'LineWidth', 2);
quiver3(kneeOrigin(1),kneeOrigin(2),kneeOrigin(3), diaphyseVector(1), diaphyseVector(2), diaphyseVector(3), 350, 'LineWidth', 2);
quiver3(kneeOrigin(1),kneeOrigin(2),kneeOrigin(3), 0, 0, -1, 70, 'LineWidth', 2);



%% calculate angles %%%% 
%NSA in 3D
NSA = atan2d(vecnorm(cross(neckAxis,diaphyseVector)),dot(neckAxis,diaphyseVector));
disp(['Neck-Shaft Angle = ' num2str(NSA, 4) '°']);
% AVA from transverse plane --> set y coordinate to 0;  %%%%%%% needs change! femur might be wrongly oriented!
u = neckAxis([1 3]);
if strcmpi(side, 'r') 
    v = [0 1];
else
    v = [0 -1];
end
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
AVA_inTransverse = real(acosd(CosTheta));
disp(['Anteversion in transverse plane = ' num2str(AVA_inTransverse, 3) '°']);



%% function to update plot when choosing neck axis cutoff
function updateNeckAxisLengthPlot(distance, nodes, D)
global h1;
global h2;
set(h1,'Visible','off');
set(h2,'Visible','off');
temp = nodes;
temp(D > distance, :) = [];
oldView = get(gca, 'View');
h1 = plotNodes(temp, 'r');
temp = nodes;
temp(or(D < distance, D > 100), :) = [];
h2 = plotNodes(temp, 'b');
view(oldView);
daspect([1 1 1]);
global neckAxisCutOff;
neckAxisCutOff = distance;
end