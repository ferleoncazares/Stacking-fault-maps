%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---------------   m24 factors   ----------------%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Slip planes analysed in the next order: (111), (-1-11), (1-1-1), (-11-1)
% Order of partials to ensure ISF formation and use of the same line vector
%
% 3 coordinate systems considered (CSsample -R1-> CScrystal -R2-> CSss): 
%   - CSsample: CS so that a uniaxial force is applied along [1,0,0]
%   - CScrystal: Miller indices
%   - CSss: CS aligned for each slip system so that
%           b = [1,0,0], cross(b,n) = [0,1,0], n = [0,0,1]
%           (b = slip direction, n = slip plane normal)
%
% Inputs:
% mode          Loading mode:
%               mode = 1 -> Tension; mode = -1 -> compression;
%               mode = 2 -> Pure shear of a (111) plane; mode = -2 -> Same (opposite direction)
%               size(mode)=[3,3] -> mode = Scrystal normalised
% Fdirection    Direction of the force or shear appplied
% m             Display format (12 or 24 slip systems)
% dispT         Print options (1 -> Tr,  2 -> Tr sorted by mS)
%
% Outputs:
% Tr            Table with the results
%  - id         Slip system identifier
%  - n          Slip normal
%  - b          Slip direction (perfect a/2<110>{111} dislocations)
%  - mS         Schmid factor
%  - mE         Escaig factor
%  - mN         Normal factor
%  - bp1        Slip direction of leading partials (a/6<112>{111} dislocations)
%  - bp2        Slip direction of trailing partials (a/6<112>{111} dislocations)
%  - mSp1       Schmid factor of bp1 leading partial
%  - mSp2       Schmid factor of bp2 trailing partial
%
% taup          Shear stresses on each (111) plane
% phi           Angle [-30deg,30deg] of the shear stress
%               phi = -30 -> pointing towards a leading partial
%               phi = 0   -> pointing towards a perfect dislocation
%               phi = 30  -> pointing towards a trailing partial
% mN            Normal factor for each (111) plane
% Scrystal      Stress tensor given in crystal reference system

function [Tr,taup,phi,mN,Scrystal] = m24(mode,Fdirection,m,dispT)

if ~((m == 12) || (m == 24))
    error('m must be 12 or 24')
end
if ~(all(size(mode) == [3,3]) || (abs(mode) == 1) || (abs(mode) == 2))
    error('mode must be -1, 1, -2, 2 or a 3-by-3 matrix')
end

%%%%%%%%%%%%% Calculation
%%% fcc
% Slip system ID
id = (1:24)';

% Glide Planes
n =  [ 1, 1, 1;     % d
       1, 1, 1;
       1, 1, 1;
      -1,-1, 1;     % c
      -1,-1, 1;
      -1,-1, 1;
       1,-1,-1;     % b
       1,-1,-1;
       1,-1,-1;
      -1, 1,-1;     % a
      -1, 1,-1;
      -1, 1,-1];
n = [n;-n];         % By artifitially inverting nsp rather than nsd -> (mS13 = -mS1 & mE13 = mE1) 
% Glide Directions
b =  [-1, 0, 1;     % CB
       1,-1, 0;     % BA
       0, 1,-1;     % AC
       1, 0, 1;     % DA
      -1, 1, 0;     % AB
       0,-1,-1;     % BD        
       1, 1, 0;     % DC
       0,-1, 1;     % CA
      -1, 0,-1;     % AD      
       0, 1, 1;     % DB
       1, 0,-1;     % BC
      -1,-1, 0];    % CD
b = [b;b];

% Slip direction of the partials: 1 = leading, 2 = trailing
% This analysis considers the partials 1=right and 2=left, 
% with mS>0 inducing a force in the dislocation towards the right.
bp1i= [-2, 1, 1;    % dB
        1,-2, 1;    % dA
        1, 1,-2;    % dC
        2,-1, 1;    % gA
       -1, 2, 1;    % gB
       -1,-1,-2;    % gD        
        1, 2,-1;    % bC
        1,-1, 2;    % bA
       -2,-1,-1;    % bD         
       -1, 1, 2;    % aB
        2, 1,-1;    % aC
       -1,-2,-1];   % aD
bp2i= [-1,-1, 2;    % Cd
        2,-1,-1;    % Bd
       -1, 2,-1;    % Ad
        1, 1, 2;    % Dg
       -2, 1,-1;    % Ag
        1,-2,-1;    % Bg         
        2, 1, 1;    % Db
       -1,-2, 1;    % Cb
       -1, 1,-2;    % Ab
        1, 2, 1;    % Da
        1,-1,-2;    % Ba
       -2,-1, 1];   % Ca   
bp1 = [bp1i;bp2i];
bp2 = [bp2i;bp1i];


% Determining Scrystal
if all(size(mode) == [3,3])     % Stress tensor given by the user
    Scrystal = mode;
elseif abs(mode) == 1           % Tension/compression
    if Fdirection(2) == 0                           % To avoid division by 0
        Fdirection(2) = max(Fdirection)/1e10;
    end
    F = Fdirection'/norm(Fdirection);               % Normalising load direction
    Ssample = [mode,0,0;0,0,0;0,0,0];               % Uniaxial stress state in CSsample

    G = zeros(3,1);
    G(1) = F(2)/sqrt(F(1)^2+F(2)^2);
    G(2) = -F(1)*G(1)/F(2);
    H = cross(F,G);
    R1 = F*[1,0,0] + G*[0,1,0] + H*[0,0,1];         % Rotation matrix: CSsample -R1-> CScrystal
    Scrystal = R1*Ssample*R1';                      % Stress state in CScrystal
elseif abs(mode) == 2           % Shear
    if Fdirection(2) == 0                           % To avoid division by 0
        Fdirection(2) = max(Fdirection)/1e10;
    end
    F = Fdirection'/norm(Fdirection);               % Normalising shear direction
    P = [1;1;1]/norm([1,1,1]);                      % Normal to the slip plane
    G = cross(P,F);
    if sum(F.*P) ~= 0
        error('Pure shear -> F must be perpendicular to [1,1,1]')
    end
    Ssample = [0,0,mode;0,0,0;mode,0,0]/2;          % Pure shear in the CSsample
    R1 = F*[1,0,0] + G*[0,1,0] + P*[0,0,1];         % Rotation matrix: CSsample -R1-> CScrystal
    Scrystal = R1*Ssample*R1';                      % Stress state in CScrystal
end

E = sort(eig(Scrystal));
Scrystal = Scrystal/(max(E)-min(E));

% Schmid and Escaig factors
mS = zeros(24,1);
mE = zeros(24,1);
mN = zeros(24,1);
mSp1 = zeros(24,1);
mSp2 = zeros(24,1);
for i = 1:24  
    % Rotation matrix: CScrystal -R2-> CSss
    w = cross(n(i,:),b(i,:));
    R2 = [1;0;0]*b(i,:)/sqrt(2) + [0;1;0]*w/sqrt(6) + [0;0;1]*n(i,:)/sqrt(3);
    Sss = R2*Scrystal*R2';                                              % Stress state in CSss
    mS(i) = Sss(1,3);                                                   % Schmid factor
    mE(i) = Sss(2,3);                                                   % Escaig factor
    if mS(i) < 0                                                        % Reassigning bp1 and bp2 vectors in case mS<0 -> dislocations i and i+12 will now have the same mSp1
        [bp1(i,:),bp2(i,:)] = deal(-bp2(i,:),-bp1(i,:));
    end
    mSp1(i) = dot(Scrystal*n(i,:)'/sqrt(3),bp1(i,:)'/sqrt(6));          % Schmid factor for each partial dislocation: mSp1 + mSp2 = mS
    mSp2(i) = dot(Scrystal*n(i,:)'/sqrt(3),bp2(i,:)'/sqrt(6));
    mN(i) = Sss(3,3);                                                   % Stress normal to the slip plane
end

% mS(abs(mS)<1e-15) = 0;          % To ignore small numerical deviations from 0
% mE(abs(mE)<1e-15) = 0;
% mS2(i) = (sum(F'.*n(i,:),2)/sqrt(3) .* sum(F'.*b(i,:),2)/sqrt(2))./sum(F'.^2);    % Gives the same result, directly from the calculations of the individual cosines
% mS3(i) = dot(Scrystal*n(i,:)'/sqrt(3),b(i,:)'/sqrt(2));                           % ",                              from the defintion of the Schmid factor m=Sn.b 

% Summary
Tr = table(id,n,b,mS,mE,mN,bp1,bp2,mSp1,mSp2);              % Table with all the results 
if m == 12
    Tr(13:end,:) = [];
end

[taup,phi,mN] = deal(zeros(4,1));
np = reshape(1:12,3,4)';        % Indicates which slip systems belong to each slip plane
for i=1:4
    taup(i) = sqrt(Tr.mS(np(i,1)).^2+Tr.mE(np(i,1)).^2);
    [~,I] = max(abs(Tr.mS(np(i,:))));
    phi(i) = atan2d(Tr.mE(np(i,I)),abs(Tr.mS(np(i,I))));
    mN(i) = Tr.mN(np(i,1));
end  
phi(abs(phi)>30.0000001) = NaN;

% Display options
if dispT == 1           % Tr   
    disp(Tr)
elseif dispT == 2       % Tr sorted by mS
    disp(sortrows(Tr,{'mS','mSCS','mE','mECS'},'descend'))
end


end