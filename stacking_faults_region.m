%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                      SFs regions                       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script determines the SFs and SSFs that will form for a given
% material and stress state.
%
% Inputs:
% tS                Schmid stresses [Pa]
% tE                Escaig stresses [Pa]
% gamma             Array with the fault energies [J/m2]
%   - gamma.ISF     Intrinsic stacking fault energy
%   - gamma.ESF     Extrinsic stacking fault energy 
%   - gamma.gTB     Matrix twin boundary energy
%   - gamma.SISF    Superlattice intrinsic stacking fault energy 
%   - gamma.SESF    Superlattice extrinsic stacking fault energy 
%   - gamma.gpTB    Twin boundary fault energy 
%   - gamma.APB     Antiphase boundary energy
% b                 Burgers vector [m]
% R                 Athermal obstacle bypassing forces [N/m]
% Rt                Thermally activated bypassing mechanisms [N/m]
% Ff1               Friction force on the leading partial [N/m]
% Ff2               Friction force on the trailing partial [N/m]
%
% Outputs:
% C                 Conditions
% M                 Regions for individual configurations
% MRa               Regions for athermal SF formation
% MRt               Regions for athermal SF and thermal SSf formation
% Clrs              Colour codes for the existing regions
% Lbls              Labels for the existing regions

function [C,M,MRa,MRt,Clrs,Lbls] = stacking_faults_region(tS,tE,gamma,b,R,Rt,Ff1,Ff2)

% Conditions
C.C1 = b/2*tS - b/(2*sqrt(3))*tE - gamma.ISF - Ff1 >= R.R1;             % The leading partial will always bypass the obstacles
C.C2 = b/2*tS + b/(2*sqrt(3))*tE + gamma.ISF + Ff2 >= -R.R1;            % The trailing partial will always bypass the obstacles
C.C3 = b/2*tS + b/(2*sqrt(3))*tE + gamma.ISF - Ff2 < R.R2;              % The trailing partial will never bypass the obstacles
C.C4 = tS > 1/b*(R.R1 + R.R2 + Ff1 + Ff2);                              % The trailing partial will become mobile before the leading one gets blocked
C.C5 = tE <= -sqrt(3)/b*(2*gamma.ISF + R.R1 - R.R2 + Ff1 - Ff2);        % The force on the leading partial will always be higher than that on the trailing one

C.C1I = 2*(b/2*tS + b/(2*sqrt(3))*tE) - gamma.SISF - 2*Ff1 - Ff2 > Rt.Rt1I; % The leading SISF partials will bypass the precipitate     
C.C2I = b/2*tS + b/(2*sqrt(3))*tE + gamma.ISF - Ff2 > R.R1;                 % The trailing SISF partial will bypass the precipitate
C.C1E = 2*(b/2*tS - b/(2*sqrt(3))*tE) - gamma.SESF - 2*Ff1 > Rt.Rt1E;       % The leading SESF partials will bypass the precipitate
C.C2E = b/2*tS - b/(2*sqrt(3))*tE + gamma.ESF - 2*Ff2 > R.R1;               % Both trailing SESF partials will bypass the precipitate
C.C3E = b/2*tS + b/(2*sqrt(3))*tE + gamma.ESF - gamma.ISF - Ff2 > R.R1;     % Only one trailing SESF partial will bypass the precipitate

% Regions
M.i = ~C.C1 & ~C.C4;                                        % Immobile considering only athermal
M.inf = ~C.C2;
M.sinf = C.C2 & C.C1 & C.C3;
M.md = ~C.C3 & C.C4 & C.C5;
M.mc = C.C4 & ~C.C5;

M.ti = ~C.C1 & ~C.C4 & ~C.C1I & ~C.C2I & ~C.C1E & ~C.C2E;   % Immobile considering also thermal
M.I1 = ~C.C1 & ~C.C4 & C.C1I & ~C.C2I;
M.I2 = ~C.C1 & ~C.C4 & C.C2I;
M.E1 = ~C.C1 & ~C.C4 & C.C1E & ~C.C2E & ~C.C3E;
M.E2 = ~C.C1 & ~C.C4 & C.C2E;
M.E3 = ~C.C1 & ~C.C4 & C.C1E & ~C.C2E & C.C3E;

% Color-coded regions: athermal & thermal
MRa =  M.i + 2*M.inf + 3*M.sinf + 4*M.md + 5*M.mc;
MRt = M.ti + 2*M.inf + 3*M.sinf + 4*M.md + 5*M.mc + 6*M.I1 + 7*M.I2 + 8*M.E1 +10*M.E2 + 12*M.E3;

% Selecting existing indices - labels and colors
v1 = [1:8,10,12,14:19];      % Current indices
v2 = 1:length(v1);           % New indices
for i = 1:length(v1)
   if v1(i) ~= v2(i)
      MRt(MRt == v1(i)) = v2(i); % Changing MRt to new indices
   end
end

MRti = unique(MRt)';        % Existing indices in MRTi
Lbls = {'Immobile dislocations';'Infinitely long ISF';'Semi-infinite ISF';...
        'Glide of decorrelated partials';'Glide of correlated partials';'Extended SISF+ISF';...
        'Isolated SISF';'Extended SESF+ESF';'Isolated SESF';...
        'Extended SESF+ISF';'Ext. SISF+ISF & Ext. SESF+ESF';'Iso. SISF & Ext. SESF+ESF';...
        'Ext. SISF+ISF & Iso. SESF';'Iso. SISF & Iso. SESF';'Ext. SISF+ISF & Ext. SESF+ISF';...
        'Iso. SISF & Ext. SESF+ISF'};
Lbls = Lbls(MRti);
Clrs = [217,217,217;... % Colors
        152,145, 38;...
        237,224, 34;...
         88,179,  0;...
          0,113, 51;...
          0, 69,140;...
        100,137,197;...
        157,  4, 21;...
        187, 92,160;...
        242,159, 13;...
       105,  2, 14;...
        52,  1,  7;...
        125, 61,107;...
         93, 45, 79;...
        110, 72,  6;...
        220,145, 12]/255;
Clrs = Clrs(MRti,:);

% Final transformation of MRti
v3 = unique(MRt)';                     % Existing indices
v4 = 1:length(unique(MRt));           % Final indices
for i = 1:length(v3)
   if v3(i) ~= v4(i)
      MRt(MRt == v3(i)) = v4(i); % Changing MRt to new indices
   end
end

end