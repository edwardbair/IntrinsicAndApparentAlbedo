
wl=[0.28:0.01:2.8]';

el=3e3;
units='um';
fSCA=[1:-0.1:0]';

R0=table(wl,zeros(length(wl),1),'VariableNames',...
    {'wavelength','reflectance'});
R0.Properties.VariableUnits={units,''};
fshade=1-fSCA;
rg=100:200:1000;
cosZ=0.1:0.2:1;
lapconc=0:200:1000;
a_diff=zeros(length(fSCA),length(cosZ),length(rg),length(lapconc));

for i=1:length(fSCA)
    for j=1:length(cosZ)
        for k=1:length(rg)
            parfor kk=1:length(lapconc)
        P = setPrescription('snow','radius',rg(k),'wavelength',wl,...
            'r0',R0,'cosZ',cosZ(j),'waveUnit',units,'LAPradius',3,'LAP','dust',...
            'LAPfraction',lapconc(kk)*1e-6,'elevation',el,'lookup',true,'fSCA',[fSCA(i) fshade(i)]);
        
        R = SPIReS_fwd(P);
        S=SolarScale('units',units,'elevation',el,'cosZ',cosZ(j),'wavelength',wl);
        
        a_apparent=trapz(wl,S.Global.*R)/trapz(wl,S.Global);
        P2=P;
        P2.snow.fSCA=[1 0];
        R2 = SPIReS_fwd(P2);
        a_intrinsic=trapz(wl,S.Global.*R2)/trapz(wl,S.Global);
        a_diff(i,j,k,kk)=a_intrinsic-a_apparent;
            end
        end
    end
end

% X=repmat(fshade',[1 length(cosZ) length(rg) length(lapconc)]);
figure;

X=[fshade;flipud(fshade)];
Y=zeros(length(fshade),2);
for i=1:length(fshade)
    yy=a_diff(i,:,:,:);
    yy=yy(:);
    Y(i,1)=min(yy);
    Y(i,2)=max(yy);
end
Y=[Y(:,1);flipud(Y(:,2))];

fill(X,Y,[0.6 0.6 0.6],'EdgeColor','none');
xlim([0 1]);ylim([0 1]);
hold on;
plot([0 1],[0 1],'--k','LineWidth',2);
xlabel('fshade');ylabel('\alpha_{intrinsic}-\alpha_{uncorrected}');
axis square
set(gca,'FontSize',20)