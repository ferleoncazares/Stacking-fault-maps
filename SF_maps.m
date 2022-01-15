%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the stress orientation dependence for the appearance of 
% the different stacking fault (SF) and superlattice stacking fault (SSF)
% configurations in Ni-based superalloys strengthened by cuboidal
% precipitates. The derivation and further details of maps developed can be
% found in:

  % F.D. León-Cázares, F. Monni, C.M.F. Rae. Stress orientation dependence
  % for the propagation of stacking faults and superlattice stacking faults
  % in nickel-based superalloys. Acta Materialia, 199 (2020) 209-224.
  
% Coded by F.D. León-Cázares

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Material parameters (elastic, microstructure, fault energies) can be 
% selected in the first section. The two following sections produce plots:
% - Regions of stacking fault configurations: full stress map with the 
%   different regions.
% - Uniaxial loading inverse pole figure: plots the stress orientation maps
%   for all possible uniaxial loading directions in tension or compression.
%   This section needs mtex (tested with version 5.1.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Material inputs
nu = 0.33;                      % Poisson's ratio [-]
mu = 58.6 * 1e9;                % Shear modulus [Pa]
b = 0.253 * 1e-9;               % Burgers vector [m]
theta = 90;                     % Dislocation character angle [degree]: 0 -> screw, 90 -> edge 
gamma.ISF = 0.011;              % Intrinsic stacking fault energy [J/m^2]         - RR1000 @700C [Galindo SF]
gamma.ESF = gamma.ISF;          % Extrinsic stacking fault energy [J/m^2]                - Assumed to be the same as gammaISF
gamma.gTB = gamma.ISF/2;        % Matrix twin boundary energy [J/m^2]                    - Assumed as half of gammaISF [Chandran2011]
gamma.SISF = 0.062;             % Superlattice intrinsic stacking fault energy [J/m^2] - [Vorontsov 2011 - Ni3Al]
gamma.SESF = 0.074;             % Superlattice extrinsic stacking fault energy [J/m^2] - [Vorontsov 2011 - Ni3Al]
gamma.gpTB = gamma.SESF/2;      % Twin boundary fault energy [J/m2]
gamma.APB = 0.276;              % Antiphase boundary energy [J/m^2]     - [Vorontsov 2011 - Ni3Al]
alpha = 1;                      % Line tension parameter [-]
H = 20*1e-9;                    % Channel width [m]
Ff = 0*100e6*b/3;               % Friction forces on perfect dislocations [N/m]
Ff1 = 0.5*Ff;                   % Friction forces of each partial [N/m]
Ff2 = Ff-Ff1;


%% Regions of stacking fault configurations
%

%%% Input - plotting range
SS = [0,1.5e9];                 % Schmid stress [Pa]
SE = [-1.5e9,1.5e9];            % Escaig stress [Pa]
dS = 1e7;                       % Stress step size [Pa]

%%% Calculations
[R,Rt,~] = sf_channel(mu,b,nu,alpha,theta,H); % Precipitate bypassing resistance
[tS,tE] = meshgrid(SS(1):dS:SS(2),SE(1):dS:SE(2));
[~,~,~,MRt,Clrs,Lbls] = stacking_faults_region(tS,tE,gamma,b,R,Rt,Ff1,Ff2); % Regions

%%% Plot
figure
scatter(reshape(tS',[],1)/1e9,reshape(tE',[],1)/1e9,20,reshape(MRt',[],1),'filled')
colormap(Clrs(min(min(MRt)):max(max(MRt)),:))
ga = gca;
gf = gcf;
xlim([SS(1)/1e9,SS(2)/1e9])
xlabel('Schmid stress [GPa]')
ylabel('Escaig stress [GPa]')
set(gca,'fontsize',12)
cb = colorbar;
cb.Ticks = 1.5+(0:(length(Lbls)-1))*(length(Lbls)-2)/(length(Lbls)-1);
cb.TickLabels = Lbls;
ga.Position(3) = ga.Position(3)/2;
gf.Position(3) = gf.Position(3)*0.9;

%}

%% Uniaxial loading inverse pole figure (IPF)
%

% Make sure mtex has been initialised!
% startup_mtex

%%% Inputs
S = 800e6;                      % Applied stress [Pa]
mode = 1;                       % Loading mode: 1 -> Tension, -1 -> Compression
res = 0.4;                      % Resolution for the IPF [degree]
ms = 3;                         % Marker size for the IPF plot

%%% Calculations
CS = crystalSymmetry('m-3m');
sR = fundamentalSector(CS);                         % The fundamental sector only
r = plotS2Grid(sR,'resolution',res*degree);         % Map of the points to analyze
r = Miller(r,CS);

[R,Rt,A] = sf_channel(mu,b,nu,alpha,theta,H); % Precipitate bypassing resistance

Fdirs = r(:);                           % Analysing multiple loading orientations
Fdirs = [Fdirs.h,Fdirs.k,Fdirs.l];
[mS,mE] = deal(zeros(size(Fdirs,1),12));    
for i = 1:size(Fdirs,1)              
    [Tr] = m24(mode,Fdirs(i,:),12,0);   % Schmid and Escaig factors
    mS(i,:)= Tr.mS;
    mE(i,:) = Tr.mE;
end
tS = S*mS(:,1);                         % Schmid and Escaig stresses
tE = S*mE(:,1);
[~,~,~,MRt,Clrs,Lbls] = stacking_faults_region(abs(tS),tE,gamma,b,R,Rt,Ff1,Ff2); % Regions

%%% Plots

figure          % IPF mapped onto the SOM
scatter(reshape(tS',[],1)/1e9,reshape(tE',[],1)/1e9,20,reshape(MRt',[],1),'filled')
colormap(Clrs(min(min(MRt)):max(max(MRt)),:))
ga = gca;
xlabel('Schmid stress [GPa]')
ylabel('Escaig stress [GPa]')
set(gca,'fontsize',12)
cb = colorbar;
cb.Ticks = 1.5+(0:(length(Lbls)-1))*(length(Lbls)-2)/(length(Lbls)-1);
cb.TickLabels = Lbls;
axis equal
xlim([min(min(tS))/1e9-0.05,max(max(tS))/1e9+0.05])
ylim([min(min(tE))/1e9-0.05,max(max(tE))/1e9+0.05])
ga.Position(1) = 0;
box on
grid on

figure          % IPF
plot(r,reshape(MRt',[],1),'MarkerSize',ms)
ga = gca;
gf = gcf;
set(gca,'fontsize',12)
colormap(Clrs(min(min(MRt)):max(max(MRt)),:))
cb = colorbar;
cb.Ticks = 1.5+(0:(length(Lbls)-1))*(length(Lbls)-2)/(length(Lbls)-1);
cb.TickLabels = Lbls;
gf.Position(3) = gf.Position(3)*1.5;

%}

