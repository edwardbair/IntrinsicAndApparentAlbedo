function [S] = terrainCorrectedAlbedo(illumCosine,viewF,Z,R,T,whichRefl)
% [S] = terrainCorrectedAlbedo(illumCosine,viewF,Z,R,T,whichRefl)
%terrain adjusted reflectance based on intrinsic direct and diffuse
%Input
% illumCosine - 2D or 3D spatial values of cosine of illumination
% viewF - 2D or 3D spatial values of view Factors
% Z - 2D or 3D spatial values of elevations
% R - raster reference
% T - table of values
% whichRefl - 'bband','vis', or 'nir'
%
%Output
% structure with Uncorrected and Spatially Corrected albedos, along with
% intrinsic values

% lat-lon of surface, same for all layers of illumCosine and viewF
[x,y] = worldGrid(R);
[lat,lon] = projinv(R.ProjectedCRS,x,y);
meanLat = mean(lat,'all','omitnan');
meanLon = mean(lon,'all','omitnan');

% check that illumCosine, viewF, and T correspond
assert(isequal(size(illumCosine,3),size(viewF,3),height(T)),...
    '3rd dimension of illumCosine and viewF must equal height of T')

% solar geometry
[declin,~,omega] = EarthEphemeris(T.TIMESTAMP);
[mu0,~] = sunang(meanLat,meanLon,declin,omega);

% functions to adjustdiffuse fractions for incoming radiation and "bounce" for
% radiation scattered in ablation hollows (for bband, vis, or nir regions)
nmu = ceil((max(mu0)-min(mu0))/.05);
Xadj = albedoAdjustments(linspace(min(mu0),max(mu0),nmu),...
    mean(Z,'all','omitnan'));

% max values of likely reflectance
bbmax = 0.9;
nirmax = 0.8;
vismax = 0.95;

% radiation, diffuse fraction, total, vis, nIR
tbadbb = abs(T.UpLookingPlatformClear-T.UpLookingBunkerClear)>40;
tbadnir = abs(T.UpLookingPlatformRed-T.UpLookingBunkerRed)>35;
totalIrrad = T.UpLookingPlatformClear;
diffuseFraction = T.SPN1Diffuse./T.SPN1Global;
totalReflected = T.DownLookingBunkerClear;
nirIrrad = T.UpLookingPlatformRed;
nirReflected = T.DownLookingBunkerRed;
visIrrad = totalIrrad-nirIrrad;
visReflected = T.DownLookingBunkerClear-T.DownLookingBunkerRed;
totalReflected(tbadbb) = NaN;
nirReflected(tbadnir) = NaN;
visReflected(tbadbb|tbadnir) = NaN;

switch lower(whichRefl)
    case 'bband'
        apparentAlbedo = filterReflectance(T.TIMESTAMP,totalReflected./totalIrrad,bbmax);
        diffuseFraction = diffuseFraction.*Xadj.D.Broadband(mu0);
        bounceFunction = Xadj.E.Broadband;
        maxValue = bbmax;
    case 'vis'
        apparentAlbedo = filterReflectance(T.TIMESTAMP,visReflected./visIrrad,vismax);
        diffuseFraction = diffuseFraction.*Xadj.D.visible(mu0);
        bounceFunction = Xadj.E.visible;
        maxValue = vismax;
    case 'nir'
        apparentAlbedo = filterReflectance(T.TIMESTAMP,nirReflected./nirIrrad,nirmax);
        diffuseFraction = diffuseFraction.*Xadj.D.nearIR(mu0);
        bounceFunction = Xadj.E.nearIR;
        maxValue = nirmax;
    otherwise
        error('whichFlag must be ''bband'', ''vis'', or ''nir''')
end

%% terrain effect by day
gridAlbedo = nan(size(illumCosine));
avgModel = zeros(size(apparentAlbedo));
intrinsicRefl = zeros(size(apparentAlbedo));
noReRefl = zeros(size(apparentAlbedo));
firstRefl = zeros(size(apparentAlbedo));
trapped = zeros(size(apparentAlbedo));
M = zeros(size(apparentAlbedo));
Mcos = zeros(size(apparentAlbedo));
for k=1:size(illumCosine,3)
    [avgModel(k),gridAlbedo(:,:,k),intrinsicRefl(k),noReRefl(k),firstRefl(k),trapped(k)] =...
        spatialEstimate(illumCosine(:,:,k),viewF(:,:,k),mu0(k),...
        diffuseFraction(k),apparentAlbedo(k),bounceFunction);
    [slope,~] = topographicSlope(Z(:,:,k),R);
    Mcos(k) = mean(1./cosd(slope),'all','omitnan');
    M(k) = mean(illumCosine(:,:,k),'all','omitnan')/mu0(k);
end
S.measuredAlbedo = apparentAlbedo;
S.gridAlbedo = gridAlbedo;
S.modelAlbedo = avgModel;
S.intrinsicAlbedo = filterReflectance(T.TIMESTAMP,intrinsicRefl,maxValue*1.03);
S.firstRefl = filterReflectance(T.TIMESTAMP,firstRefl,maxValue);
S.withoutReRefl = noReRefl;
S.trapped = trapped;
S.M = M;
S.Mcos = Mcos;
end

function [avgM,modA,intrinsicR,noReRefl,firstRefl,trapped] =...
    spatialEstimate(muS,vF,mu0,dF,measA,bounceFunction)
% solve for directRefl and diffuseRefl that best match the measured
% apparent albedo
M = muS/mu0;
noReRefl = mean(measA./(M*(1-dF)+dF*vF),'all','omitnan');
x0 = [max(0.01,measA/2) min(.99,measA*2)];
[x,fval,exitflag,output] = fzero(@measMinusModel,x0); %#ok<ASGLU>
[avgM,modA,firstRefl,trapped] = spatialCorrection(M,vF,dF,bounceFunction,x);
intrinsicR = x;
    function xdiff = measMinusModel(x)
        averageModel = spatialCorrection(M,vF,dF,bounceFunction,x);
        xdiff = measA-averageModel;
    end
end
function [avgM,varargout] = spatialCorrection(M,viewF,df,bounceF,intrinsicR)
% model to match measured albedo based on topography correction of intrinsic
% albedo
% optional output separates initially reflected albedo from those that go
% through multiple boundes
origSize = size(viewF);
% change to column vectors for calculation
notNaN = ~isnan(viewF) & ~isnan(M);
vF = viewF(notNaN);
M = M(notNaN);
% direct and diffuse reflectance the same, hook here to change if necessary
directR = intrinsicR;
diffuseR = intrinsicR;
% all reflectance at top, goes in either direction
topIn = (1-df)*directR.*M+df*diffuseR.*vF+(1-vF).*diffuseR.^2;
% initial reflections that escape and don't escape
N = 20; % max iterations, loop always stops earlier
aResc = zeros(length(vF),N);
aRint = zeros(size(aResc));
aResc(:,1) = vF.*topIn;
aRint(:,1) = (1-vF).*topIn;
% re-reflected or escaping internally
tol = 1e-6;
thisRefl = bounceF(repmat(diffuseR,1,N-1),1:(N-1));
for n=2:N
    aResc(:,n) = thisRefl(n-1)*aRint(:,n-1).*vF;
    aRint(:,n) = thisRefl(n-1)*aRint(:,n-1).*(1-vF);
    if mean(aResc(:,n))<=tol
        break;
    end
end
apparentRefl = sum(aResc,2);
avgM = mean(apparentRefl);
if nargout>1
    reflMap = nan(origSize);
    reflMap(notNaN) = apparentRefl;
    varargout{1} = reflMap;
    if nargout>2
        varargout{2} = mean(aResc(:,1));
        if nargout>3
            varargout{3} = mean(aRint(:,1)-sum(aResc(:,2:end),2));
        end
    end
end
end