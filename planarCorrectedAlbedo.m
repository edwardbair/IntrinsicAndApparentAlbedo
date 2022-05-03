function [S] = planarCorrectedAlbedo(Z,R,T,usePC)
% [S] = planarCorrectedAlbedo(Z,R,T,usePC)
% calculate albedo corrected for planar representation of surface
%
%Input
% Z - 2D or 3D elevations, 3rd dimension corresponding to dates in T
% R - raster reference
% T - table of radiation values
% usePC - if true, use MATLAB point-cloud routines, if false, use fit
%
%Output
% structure with uncorrected and planar corrected albedos

% lat-lon of surface, same for all layers of Z
[x,y] = worldGrid(R);
[lat,lon] = projinv(R.ProjectedCRS,x,y);
meanLat = mean(lat,'all','omitnan');
meanLon = mean(lon,'all','omitnan');

% check that Z and T correspond
assert(size(Z,3)==height(T),'3rd dimension of Z must equal height of T')

% solar geometry
[declin,~,omega] = EarthEphemeris(T.TIMESTAMP);
[mu0,phi0] = sunang(meanLat,meanLon,declin,omega);

% max values of likely reflectance
bbmax = 0.9;
nirmax = 0.8;

% radiation, diffuse fraction, albedo from T
tbadbb = abs(T.UpLookingPlatformClear-T.UpLookingBunkerClear)>40;
tbadnir = abs(T.UpLookingPlatformRed-T.UpLookingBunkerRed)>35;
irradBroadband = T.UpLookingPlatformClear;
diffuseFraction = T.SPN1Diffuse./T.SPN1Global;
irradNIR = T.UpLookingPlatformRed;
reflectedBroadband = T.DownLookingBunkerClear;
reflectedBroadband(tbadbb) = NaN;
reflectedNIR = T.DownLookingBunkerRed;
reflectedNIR(tbadnir) = NaN;
apparentAlbedo.Broadband = filterReflectance(T.TIMESTAMP,...
    reflectedBroadband./irradBroadband,bbmax);
apparentAlbedo.nIR = filterReflectance(T.TIMESTAMP,...
    reflectedNIR./irradNIR,nirmax);

% planar correction, by day
M = zeros(size(irradBroadband));
Mcos = zeros(size(M));
% broadcast variables
passVals = {x,y,R,T};
W = parallel.pool.Constant(passVals);
% warning('off','curvefit:fit:iterationLimitReached')
warning('off')
% weird point cloud 30 Jan 2021
badPCdate = datetime('30-Jan-2021 10:45:00','TimeZone','Etc/GMT+8');
parfor k=1:size(Z,3)
    % broadcast variables
    pv = W.Value;
    Tx = pv{4};
    % properties of the plane
    % option to use the MATLAB point cloud routines, vs using fit
    Zq = fitPlane(usePC || Tx.TIMESTAMP(k)==badPCdate,...
        pv{1},pv{2},Z(:,:,k),pv{3});
    % illumination angle on the plane
    [Sp,Ap] = topographicSlope(Zq,pv{3});
    muS = sunslope(mu0(k),phi0(k),mean(Sp(:),'omitnan'),mean(Ap(:),'omitnan'));
    M(k) = muS/mu0(k);
    Mcos(k) = mean(1./cosd(Sp),'all','omitnan');
end
% warning('on','curvefit:fit:iterationLimitReached')
warning('on')

% adjustment of diffuse parameters
Xadj = albedoAdjustments(mu0,mean(Z,'all','omitnan'),true);

% planar corrected albedos
bbDiffuse = diffuseFraction.*Xadj.D.Broadband(mu0);
nirDiffuse = diffuseFraction.*Xadj.D.nearIR(mu0);
planeAlbedo.Broadband = filterReflectance(T.TIMESTAMP,reflectedBroadband./...
    (M.*(1-bbDiffuse).*irradBroadband+bbDiffuse.*irradBroadband),bbmax*1.03);
planeAlbedo.nIR = filterReflectance(T.TIMESTAMP,reflectedNIR./...
    (M.*(1-nirDiffuse).*irradNIR+nirDiffuse.*irradNIR),nirmax*1.03);
S.Uncorrected = apparentAlbedo;
S.Plane = planeAlbedo;
S.M = M;
S.Mcos = Mcos;
end

function Zhat = fitPlane(usePointCloud,x,y,z,R)
nt = ~isnan(z);
if usePointCloud
    pc = pointCloud([x(nt) y(nt) z(nt)]);
    % set max # trials to 5x pixel spacing, otherwise got warning
    % message
    plane = pcfitplane(pc,5*mean([R.CellExtentInWorldX R.CellExtentInWorldY]));
    A=plane.Parameters(1);
    B=plane.Parameters(2);
    C=plane.Parameters(3);
    D=plane.Parameters(4);
    % equation for plane Ax + By + Cz + D = 0
    Zhat = -(D+A*x+B*y)/C; % Solve for z data for the plane
else
    [Fp,~,O] = fit([x(nt),y(nt)],z(nt),'poly11','robust','bisquare','normalize','on');
    if O.exitflag<1
        Fp = fit([x(nt),y(nt)],z(nt),'poly11','normalize','on');
    end
    Zhat = Fp(x,y);
end
end