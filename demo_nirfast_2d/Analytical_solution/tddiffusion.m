function Phi = tddiffusion(mua, musp, v, n, srcpos,detpos, t)       %Replaced Reff with n since using different boundary conditions
%
%  Phi = tddiffusion(mua, musp, v, Reff, srcpos,detpos, t)
%    
%  semi-infinite medium analytical solution to diffusion model
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        mua:   the absorption coefficients in 1/mm
%        musp:  the reduced scattering coefficients in 1/mm
%        v:     the speed of light
%        Reff:  the effective reflection coeff.
%        srcpos:array for the source positions (x,y,z)
%        detpos:array for the detector positions (x,y,z)
%        t:     a list of time in s at which to evaluate the 
%               analytical diffusion solution
%
%    output:
%        Phi:  the output fluence for all time points
%
%    this file is part of Monte Carlo eXtreme (MCX)
%    License: GPLv3, see http://mcx.sf.net for details
%    see Boas2002
%

D = 1/(3*(mua+musp));
% zb = (1+Reff)/(1-Reff)*2*D;       %EBC, no mismatch between boundaries (?)
%zb=0;
                                    %Robin
R0 = (n - 1)^2/(n + 1)^2;           
% critical angle (total internal reflection (from tissue to air))
theta_incidence = asin(1/n);
A = (2/(1-R0) - 1 + abs(cos(theta_incidence))^3)/(1 - abs(cos(theta_incidence))^2);
zb = 2*D*A;

z0 = 1/(musp+mua);
%3D
% r=getdistance([srcpos(:,1:2) srcpos(:,3)+z0],detpos);
% r2=getdistance([srcpos(:,1:2) srcpos(:,3)-z0-2*zb],detpos);
%2D
r=getdistance([srcpos(:,1) srcpos(:,2)+z0],detpos);        
r2=getdistance([srcpos(:,1) srcpos(:,2)-z0-2*zb],detpos);

s=4*D*v*t;
% size(s);
% r;
% size(r);
% size(r.^2);

% unit of phi:  1/(mm^2*s)
Phi =v./((s*pi).^(3/2)).*exp(-mua*v*t).*(exp(-(r.^2)./s) - exp(-(r2.^2)./s));
