%-------------------------------------------------------------------------%
% Copyright (c) 2022 Koller W.                                            %
%    Author:   Willi Koller,  2022                                        %
%    email:    willi.koller@univie.ac.at                                  %
% ----------------------------------------------------------------------- %

% Method_Sangeux_2015 (c) by Willi Koller, University of Vienna
%
% Method_Sangeux_2015 is licensed under a
% Creative Commons Attribution-NonCommercial 4.0 International License.
%
% You should have received a copy of the license along with this
% work. If not, see <http://creativecommons.org/licenses/by-nc/4.0/>.

% this script uses GIBBON- and MSK-STAPLE toolboxes
% https://github.com/modenaxe/msk-STAPLE
% https://github.com/gibbonCode/

clear;
addpath(genpath('./lib'));
side = 'L';
STLfileName = [];

data.markups.controlPoints(1).label = 'HJC';
data.markups.controlPoints(1).position = [-0.056276 0.85151 0.07726]';
data.markups.controlPoints(2).label = 'GT';
data.markups.controlPoints(2).position = [-0.0788038 0.80126 0.13722]';
data.markups.controlPoints(3).label = 'PS';
data.markups.controlPoints(3).position = [-0.0662013 0.780346 0.116023]';
data.markups.controlPoints(4).label = 'DS';
data.markups.controlPoints(4).position = [-0.0467056 0.439209 0.07726]';
data.markups.controlPoints(5).label = 'MC';
data.markups.controlPoints(5).position = [-0.0790706 0.439788 0.0532338]';
data.markups.controlPoints(6).label = 'LC';
data.markups.controlPoints(6).position = [-0.0790706 0.439788 0.100893]';

for i = 1 : length(data.markups.controlPoints)
    switch data.markups.controlPoints(i).label
        case 'HJC'
            hjcI = i;
        case 'GT'
            gtI = i;
        case 'PS'
            psI = i;
        case 'DS'
            dsI = i;
        case 'LC'
            lcI = i;
        case 'MC'
            mcI = i;
    end
end

neckAxis = data.markups.controlPoints(hjcI).position - data.markups.controlPoints(gtI).position;
shaftAxis = data.markups.controlPoints(psI).position - data.markups.controlPoints(dsI).position;
kneeAxis = data.markups.controlPoints(lcI).position - data.markups.controlPoints(mcI).position;

if ~isempty(STLfileName)
    filenameNoExt = strrep(STLfileName, '.stl', '');
    [stlStruct] = import_STL(fullfile(folder, STLfileName));
    
    triGeom_set = createTriGeomSet({filenameNoExt}, folder);
    triGeom_set.(['femur_' lower(side)]) = triGeom_set.(filenameNoExt);
    triGeom_set = rmfield(triGeom_set, filenameNoExt);
    
    try
        [~, JCS, ~, ~, CS] = GIBOC_femur(triGeom_set.(['femur_' lower(side)]), side, 'cylinder', 0);
    catch
        [~, JCS, ~, ~, CS] = GIBOC_femur(triGeom_set.(['femur_' lower(side)]), side, 'ellipsoids', 0);
    end
    
    neckAxisOpenSim = neckAxis' * JCS.(['hip_' lower(side)]).V ;
    kneeAxisOpenSim = kneeAxis' * JCS.(['hip_' lower(side)]).V;
    shaftAxisOpenSim = shaftAxis' * JCS.(['hip_' lower(side)]).V;
    
    vertexMAT = stlStruct.solidVertices{1, 1};
    vertexMAT = vertexMAT * JCS.(['hip_' lower(side)]).V;

    connList = triGeom_set.(['femur_' lower(side)]).ConnectivityList;
    points = triGeom_set.(['femur_' lower(side)]).Points;
    points = points * JCS.(['hip_' lower(side)]).V;
    triGeom_setOpenSim = triangulation(connList, points);
end

% NSA in 3D
NSA = atan2d(vecnorm(cross(neckAxis*(-1),shaftAxis)),dot(neckAxis*(-1),shaftAxis));
disp(['Neck-Shaft Angle = ' num2str(NSA, 4) '°']);

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0; ...
              norm(cross(A,B)) dot(A,B)  0; ...
              0              0           1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

UU = @(Fi,G) Fi*G*inv(Fi);
a=shaftAxis; 
a = a /norm(a);
b=[0 0 1]';
U = UU(FFi(a,b), GG(a,b));
norm(U);
norm(b-U*a);
U;

neckAxisRotated = U * neckAxis;
kneeAxisRotated = U * kneeAxis;
shaftAxisRotated = U * shaftAxis;

u = neckAxisRotated([1 2]);
v = kneeAxisRotated([1 2]) *(-1);
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
AVA_inPlanePerpendicularToShaft = real(acosd(CosTheta));
disp(['Anteversion in plane perpendicular to shaft = ' num2str(AVA_inPlanePerpendicularToShaft, 3) '°']);

figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.5]);
tiledlayout(1, 2);
nexttile;
hold on;
rotate3d on;
if ~isempty(STLfileName)
    plotTriangLight(triGeom_set.(['femur_' lower(side)]), JCS.(['hip_' lower(side)]), 0, 0.3);
end
daspect([1 1 1]);
quiver3(data.markups.controlPoints(gtI).position(1) ,data.markups.controlPoints(gtI).position(2) ,data.markups.controlPoints(gtI).position(3) , neckAxis(1), neckAxis(2), neckAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(dsI).position(1) ,data.markups.controlPoints(dsI).position(2) ,data.markups.controlPoints(dsI).position(3) , shaftAxis(1), shaftAxis(2), shaftAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(mcI).position(1) ,data.markups.controlPoints(mcI).position(2) ,data.markups.controlPoints(mcI).position(3) , kneeAxis(1), kneeAxis(2), kneeAxis(3), 'LineWidth', 2);
leg = legend({'Neck Axis', 'Shaft Axis', 'Knee Axis'});
leg.Layout.Tile = "north";
title('View on plane which is perpendicular to shaft axis')
view(shaftAxis);

nexttile;
hold on;
rotate3d on;
daspect([1 1 1]);
if ~isempty(STLfileName)
    plotTriangLight(triGeom_set.(['femur_' lower(side)]), JCS.(['hip_' lower(side)]), 0, 0.3);
end
quiver3(data.markups.controlPoints(gtI).position(1) ,data.markups.controlPoints(gtI).position(2) ,data.markups.controlPoints(gtI).position(3) , neckAxis(1), neckAxis(2), neckAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(dsI).position(1) ,data.markups.controlPoints(dsI).position(2) ,data.markups.controlPoints(dsI).position(3) , shaftAxis(1), shaftAxis(2), shaftAxis(3), 'LineWidth', 2);
quiver3(data.markups.controlPoints(mcI).position(1) ,data.markups.controlPoints(mcI).position(2) ,data.markups.controlPoints(mcI).position(3) , kneeAxis(1), kneeAxis(2), kneeAxis(3), 'LineWidth', 2);
title('View in transverse plane (if in OpenSim coordinate system)')
view(0, 0);

%% run only if geometry (stl) is also available

if ~isempty(STLfileName)
    figure;
    hold on;rotate3d on;
    plotTriangLight(triGeom_setOpenSim, JCS.(['hip_' lower(side)]), 0, 0.2);
    daspect([1 1 1]);
    gtRotated = data.markups.controlPoints(gtI).position' * JCS.(['hip_' lower(side)]).V;
    dsRotated = data.markups.controlPoints(dsI).position' * JCS.(['hip_' lower(side)]).V;
    mcRotated = data.markups.controlPoints(mcI).position' * JCS.(['hip_' lower(side)]).V;
    quiver3(gtRotated(1) ,gtRotated(2) ,gtRotated(3) , neckAxisOpenSim(1), neckAxisOpenSim(2), neckAxisOpenSim(3), 'LineWidth', 2);
    quiver3(dsRotated(1) ,dsRotated(2) ,dsRotated(3) , shaftAxisOpenSim(1), shaftAxisOpenSim(2), shaftAxisOpenSim(3), 'LineWidth', 2);
    quiver3(mcRotated(1) ,mcRotated(2) ,mcRotated(3) , kneeAxisOpenSim(1), kneeAxisOpenSim(2), kneeAxisOpenSim(3), 'LineWidth', 2);
    
    u = neckAxisOpenSim([1 3]);
    v = kneeAxisOpenSim([1 3]) *(-1);
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    AVA_inTransverse = real(acosd(CosTheta));
    disp(['Anteversion in transverse plane (OpenSim coordinate system) = ' num2str(AVA_inTransverse, 3) '°']);
end
