%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Stacking fault obstacle: Cuboidal Precipitates         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script calculates the forces needed to bypass an array of cuboidal
% precipitates, both athermally via glide through the channels and
% thermally by the formation of superlattice stacking faults.
%
% Nomenclature:
% pd            perfect dislocation
% p1            leading partial
% p2            trailing partial
% 
% Inputs:
% mu            Shear modulus [Pa]
% b             Burgers vector [m]
% nu            Poisson's ratio [-]
% thetaT        Angle between b and the channel orientation [deg]
% H             Channel width [m]
%
% Outputs:
% R             Athermal obstacle bypassing forces
% Rt            Thermally activated bypassing mechanisms
% A             Elastic parameter

function [R,Rt,A] = sf_channel(mu,b,nu,alpha,theta,H)

A = (2+nu-4*nu*cosd(theta)^2)/(24*pi*(1-nu)) * mu*b^2;  % Isotropic elasticity term [N] 
theta2 = theta+90;                                      % Dislocation character in the direction perpendicular to the channel
dc = H.*(4*pi*(1-nu)-2*nu.*cosd(2*theta2)+nu.*sind(30+2*theta2)-sqrt(16*pi^2*(1-nu).^2 + ...
    2*(2-nu).^2-8*pi*(2-3*nu+nu.^2)-4*(2-nu).*nu.*cosd(2*theta2)+4*nu.^2.*cosd(2*theta2).^2- ...
    (4-8*pi*(1-nu)-2*nu).*nu.*sind(30+2*theta2)+nu.^2.*sind(30+2*theta2).^2)) ./ ...
    (-4+8*pi*(1-nu)+2*nu+2*nu.*sind(30+2*theta2));      % Separation of the partials pushing against the precipitates [m]

% Athermal
R1 = alpha*mu*b^2./(3*H);                               % Force [N/m] for p1 gliding through the channel
R2 = alpha*mu*b^2./(3*(H-2*dc));                        % Force [N/m] for p2 gliding through the channel
Rd = alpha*mu*b^2./H;                                   % Force [N/m] for the perfect dislocation precipitate bypassing

% Thermal
Rt1I = 0;                                               % (1T) Force [N/m] for thermally activated SISF formation
Rt1E = 0;                                               % (1T) Force [N/m] for thermally activated SESF formation
Rt2Il = R1;
Rt2El = R1;
Rt2Ell = R1;
Rt2ET = 0;

%%% Summary
R.R1 = R1;                      % Obstacle athermal forces [N]   
R.R2 = R2;
R.Rd = Rd;
R.def1 = 'l';
R.def2 = 'll';
R.R2s = Inf;

Rt.Rt1I = Rt1I;                 % Obstacle thermal forces [N]
Rt.Rt1E = Rt1E;
Rt.Rt2Il = Rt2Il;
Rt.Rt2El = Rt2El;
Rt.Rt2Ell = Rt2Ell;
Rt.Rt2ET = Rt2ET;

end