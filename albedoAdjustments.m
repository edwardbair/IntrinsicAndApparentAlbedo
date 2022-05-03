function [X] = albedoAdjustments(mu0v,Z,varargin)
% [X] = albedoAdjustments(mu0v,Z,varargin)
%albedoAdjustments use SMARTS + SPIReS_fwd to identify statistical adjustment to
%intrinsic albedo for re-reflected radiation in broadband and near-infrared regions
% input
% mu0v - vector of cosZ values to consider
% Z - elevation
%
% output
%   X - structure with values

% variable argument, set to true to return diffuse factor only
if nargin>2
    if varargin{1}
        dfOnly = true;
    end
else
    dfOnly = false;
end

% SMARTS atmosphere
atmosphere = 'mlw';
% consider 4 snow types: clean fine-grain, clean coarse, dirty fine, dirty
% coarse (radius in um, San Juan dust)
nvals = 8;
snowGrainRadius = [50 1000];
log10Dust = [-8 -3];
snowGrainRadius = linspace(sqrt(snowGrainRadius(1)),sqrt(snowGrainRadius(2)),nvals).^2;
dustConcentration = 10.^(log10Dust(1):log10Dust(2));

% spectral ranges of PSP, red band, SPN1
bbandLimits = [285 2800]; % nm
nearIRLimits = [700 2800];
diffuseFractionLimits = [400 2700];
visLimits = [285 700];

% typical SMARTS spectra under these conditions
T = generateSMARTS(atmosphere,Z,mu0v);

% diffuse fraction in bband, nearIR, and vis compared to that in diffuseFraction
D = diffuseFactor(T,mu0v,bbandLimits,nearIRLimits,visLimits,diffuseFractionLimits);
X.D = D;
if dfOnly
    return
end

% enhancement factor in bband and nearIR as function of number of
% scattering events, could add vis
A = successiveReflection(D.Fdif,D.Fglo,T.waveL,snowGrainRadius,dustConcentration,...
    mu0v,cat(1,bbandLimits,nearIRLimits,visLimits),{'Broadband','nearIR','visible'});

% create interpolating functions for
fn = {'Broadband','nearIR','visible'};
for k=1:length(fn)
    E.(fn{k}) = interpFunction(A.(fn{k}).albedoTbl.albedo);
end
%output (X.D already exists)
X.E = E;
end

function X = diffuseFactor(T,mu0v,bbandLimits,nearIRLimits,visLimits,diffuseFractionLimits)
bbandDiffuseFactor = zeros(length(mu0v),1);
nirDiffuseFactor = zeros(size(bbandDiffuseFactor));
visDiffuseFactor = zeros(size(bbandDiffuseFactor));
for k=1:length(mu0v)
    Fdif = fit(T.waveL,T.HorzDiffuse(:,k),'pchip');
    Fglo = fit(T.waveL,T.HorzGlobal(:,k),'pchip');
    bbandDiffuseFraction = integrate(Fdif,bbandLimits(2),bbandLimits(1))/...
        integrate(Fglo,bbandLimits(2),bbandLimits(1));
    SPNdiffuseFraction = integrate(Fdif,diffuseFractionLimits(2),diffuseFractionLimits(1))/...
        integrate(Fglo,diffuseFractionLimits(2),diffuseFractionLimits(1));
    nearIRdiffuseFraction = integrate(Fdif,nearIRLimits(2),nearIRLimits(1))/...
        integrate(Fglo,nearIRLimits(2),nearIRLimits(1));
    visDiffuseFraction = integrate(Fdif,visLimits(2),visLimits(1))/...
        integrate(Fglo,visLimits(2),visLimits(1));
    bbandDiffuseFactor(k) = bbandDiffuseFraction/SPNdiffuseFraction;
    nirDiffuseFactor(k) = nearIRdiffuseFraction/SPNdiffuseFraction;
    visDiffuseFactor(k) = visDiffuseFraction/SPNdiffuseFraction;
end
whichFactor = {bbandDiffuseFactor,nirDiffuseFactor,visDiffuseFactor};
nameFactor = {'Broadband','nearIR','visible'};
for k=1:3
    y = whichFactor{k};
    F = fit(mu0v(:),y(:),'poly2','robust','bisquare');
    X.(nameFactor{k}) = F; % y = p1*x^2+p2*x+p3
end
X.Fdif = Fdif;
X.Fglo = Fglo;
end

function T = generateSMARTS(atmosphere,Z,mu0)
SetSMARTSversion(295);
for k=1:length(mu0)
    P = defaultSMARTSinput(atmosphere,'altit',Z/1000,'cosZ',mu0(k));
    S = SMARTSMain(P);
    if k==1
        waveL = S.spectralTbl.waveL;
        HorzDiffuse = S.spectralTbl.HorzDiffuse;
        HorzGlobal = S.spectralTbl.HorzGlobal;
    else
        HorzDiffuse = cat(2,HorzDiffuse,S.spectralTbl.HorzDiffuse);
        HorzGlobal = cat(2,HorzGlobal,S.spectralTbl.HorzGlobal);
    end
end
T = table(waveL,HorzDiffuse,HorzGlobal);
end

function F = interpFunction(aMatrix)
% create interpolation function for albedo values = f(aMatrix(:,1),n),
% where n is the number of bounces

% first column has the intial albedos, responses bx are in the other
% columns
bx = aMatrix(:,2:end);
allN = 1:size(bx,2); % column indices, = number of reflections
% create scattered interpolant as f(aMatrix(:,1),n)
[x,nn] = ndgrid(aMatrix(:,1),allN);
G = scatteredInterpolant(x(:),nn(:),bx(:),'linear','nearest');
% create gridded Interpolant from scattered interpolant, evenly spaced
% in albedo
minA = max(0,round(min(aMatrix(:,1)),2)-.01);
maxA = min(1,round(max(aMatrix(:,1)),2)+.01);
[y,ni] = ndgrid(minA:.01:maxA,allN);
z = G(y,ni);
z(z>.999) = .999; % (shouldn't be any)
% options for final fit, use gridded interpolant to constrain range but SLM
% tools to ensure monotonic in albedo and # of bounces
X = cat(2,y(:,1),z);
Xm = X;
for k=1:size(X,2)
    slm = slmengine(1:size(X,1),X(:,k),'increasing','on','knots',-6);
    Xm(:,k) = slmeval(1:size(X,1),slm);
end
% check to make sure monotonic in both directions
[dx,dy] = gradient(Xm);
if any(dx<0,'all') || any(dy<0,'all')
    warning('albedo growth grid is not monotonic in one or both directions, check code')
end
F = griddedInterpolant(y,ni,Xm(:,2:end),'linear','nearest');
% figure for paper, comment out to remove for ordinary operation
[xa,xb] = ndgrid(y(:,1),1:10);
% surf(xa,xb,F(xa,xb),'EdgeColor','none');

%alternate method of plotting
delta_alpha=F(xa,xb)-xa;
figure;imagesc(flipud(delta_alpha'));
xlabel('\alpha^{(0)}')
ylabel('n, number of reflections');
c=colorbar;
xt=1:4:size(xb,1);
xtl=linspace(min(xa(:,1)),max(xa(:,1)),length(xt));
xtl=num2str(xtl',"%0.2f");
yt=1:1:size(xb,2);
ytl=linspace(max(yt),min(yt),length(yt));
c.Label.String={'\alpha^{(n)}-\alpha^{(0)}'};
set(gca,'XTick',xt,'XTickLabel',xtl,...
    'YTick',yt,'YTickLabel',ytl);


end

function S = successiveReflection(Fdif,Fglo,waveL,radius,dust,cosZ,bandLimits,bandNames)
% adjustment factors for snow reflectance as f(number of scattering events),
% recognizing that the re-reflected radiation in the ablation hollows is stronger in the
% wavelengths where the snow itself is reflective and that the diffuse
% irradiance is biased toward the wavelengths where the snow is reflective
%
%input
% Fdif - diffuse irradiance interpolant
% Fglo - global irradiance interpolant
% waveL - wavelengths
% radius, dust - vectors of equal length for grain size and dust concentration
% cosZ - vector of cosines of illumination angles
% bandLimits (N x 2 matrix) - ranges of spectral bands to consider
% bandNames - cell vector of names of bands
%
%output
% S - structure with tables for statistical analysis

% all combinations of radius, dust, and cosZ
[r,d,mu] = ndgrid(radius,dust,cosZ);
r = r(:);
d = d(:);
mu = mu(:);

N = 20; % max number of bounces
for k=1:length(r)
    % direct and diffuse albedo of snow for all wavelengths
    P = setPrescription('snow','cosZ',mu(k),'radius',r(k),'LAP','dust',...
        'wavelength',waveL,'waveunit','nm','LAPfraction',d(k),'lookup',true);
    directRefl = SPIReS_fwd(P);
    P.cosZ = [];
    diffuseRefl = SPIReS_fwd(P);
    % interpolations
    Sdir = fit(waveL,directRefl,'pchip');
    Sdif = fit(waveL,diffuseRefl,'pchip');

    % initial values include both direct and diffuse irradiance but all
    % reflection is diffuse
    reflectedInit = Sdir(waveL).*(Fglo(waveL)-Fdif(waveL))+Sdif(waveL).*Fdif(waveL);

    % each wavelength band in loop
    for m=1:size(bandLimits,1)
        waveRange = bandLimits(m,:);
        % matrices to hold output, radiation at each iteration, and wavelength
        % integrated albedo
        twave = waveL>=waveRange(1) & waveL<=waveRange(2);
        w = waveL(twave);
        reflected = zeros(length(w),N);
        albedo = zeros(1,N);
        totalReflected = zeros(size(albedo));
        for n=1:N
            if n==1
                reflected(:,n) = reflectedInit(twave);
                totalReflected(n) = trapz(w,reflected(:,n));
                albedo(n) = totalReflected(n)/integrate(Fglo,waveRange(2),waveRange(1));
            else
                reflected(:,n) = Sdif(w).*reflected(:,n-1);
                totalReflected(n) = trapz(w,reflected(:,n));
                albedo(n) = totalReflected(n)/totalReflected(n-1);
            end
        end
        % results in table
        thisTbl = table(r(k),d(k),mu(k),albedo,'VariableNames',{'radius','dust','cosZ','albedo'});
        if k==1
            S.(bandNames{m}).albedoTbl = thisTbl;
        else
            S.(bandNames{m}).albedoTbl = [S.(bandNames{m}).albedoTbl; thisTbl];
        end
    end
end
end