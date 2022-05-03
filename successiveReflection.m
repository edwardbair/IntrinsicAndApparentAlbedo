function S = successiveReflection(radius,dust,cosZ)
%increase in broadband snow albedo with successive reflections
%

% spectral ranges for typical instruments
bbandLimits = [285 2800]; % nm
nearIRLimits = [700 2800];

% typical SMARTS spectrum under these conditions
atmosphere = 'mlw';
Z = 3000;
SetSMARTSversion(298);
P = defaultSMARTSinput(atmosphere,'altit',Z,'cosZ',cosZ);
X = SMARTSMain(P);
T = X.spectralTbl(:,{'waveL','HorzDiffuse','HorzGlobal'});
% interpolations
Fglo = fit(T.waveL,T.HorzGlobal,'pchip');
Fdif = fit(T.waveL,T.HorzDiffuse,'pchip');

% direct and diffuse albedo of snow for all wavelengths
P = setPrescription('snow','cosZ',cosZ,'radius',radius,'LAP','dust',...
    'wavelength',T.waveL,'waveunit','nm','LAPfraction',dust,'lookup',false);
directRefl = SPIReS_fwd(P);
P.cosZ = [];
diffuseRefl = SPIReS_fwd(P);
% interpolations
Sdir = fit(T.waveL,directRefl,'pchip');
Sdif = fit(T.waveL,diffuseRefl,'pchip');

% initial values include both direct and diffuse irradiance but all
% reflection is diffuse
reflectedInit = Sdir(T.waveL).*(Fglo(T.waveL)-Fdif(T.waveL))+Sdif(T.waveL).*Fdif(T.waveL);

N = 40;
% broadband and near-IR in loop
for k=1:2
    switch k
        case 1
            waveRange = bbandLimits;
        case 2
            waveRange = nearIRLimits;
    end
    % matrices to hold output, radiation at each iteration, and wavelength
    % integrated albedo
    twave = T.waveL>=waveRange(1) & T.waveL<=waveRange(2);
    w = T.waveL(twave);
    reflected = zeros(length(w),N);
    albedo = zeros(N,1);
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
    switch k
        case 1
            S.Broadband.reflected = reflected;
            S.Broadband.albedo = albedo;
            S.Broadband.totalReflected = totalReflected;
        case 2
            S.nearIR.reflected = reflected;
            S.nearIR.albedo = albedo;
            S.nearIR.totalReflected = totalReflected;
    end
end
end