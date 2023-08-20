% based on maketissueRaman.m
% maketissueRaman2.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcraman2.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcraman.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcraman.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueListRaman2.m
%
%   To use, 
%       1. Prepare makeTissueListRaman2.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissueRaman.m.
%
%   Note: mcraman2.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%  Dayle Kotturi, adapted for Raman Feb, 2020.  
%  Dayle Kotturi, change wavelength after inelastic scattering Jan, 2023.  

clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

myname      = 'raman10min';% name for files: myname_T.bin, myname_H.mci  
time_min    = 10;      	% time duration of the simulation [min] <----- run time -----
% nm          = 785;   	% desired wavelength of simulation
% 1/29/23 kdk: nm will become an array of all wavelengths of interest (excitation, peak1, peak2, ...)
nm          = [785 857 884 896 906];
Nbins       = 200;    	% # of bins in each dimension of cube 
binsize     = 0.01; 	% size of each bin in [cm]. Increased by 20x for SERS-sensor scale

% Set Monte Carlo launch flags
mcflag      = 0;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
launchflag  = 0;        % 0 = let mcraman.c calculate launch trajectory
                        % 1 = manually set launch vector.
boundaryflag = 1;       % 0 = no boundaries, 1 = escape at boundaries
                        % 2 = escape at surface only. No x, y, bottom z
                        % boundaries

% Sets position of source
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 0.0101;  	% z of source

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 or 3 (not 2=isotropic pt.)
radius      = 0.0300;   % 1/e radius of beam at tissue surface
waist       = 0.0300;  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0.7;      % trajectory projected onto x axis
uy0         = 0.4;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% 
% Prepare Monte Carlo 
%%%

% Create tissue properties
Nw = length(nm)
for j=1:Nw
    switch j
        case 1
            tissue_excite = makeTissueListRaman2(nm(j)); % also --> global tissue(1:Nt).s
            % 1/29/23 kdk: the above gets called for all values in nm
            Nt = length(tissue_excite);
            % 1/29/23 kdk: why do we ned to copy these when we have them in tissue?
            for i=1:Nt
                muav(j,i)  = tissue_excite(i).mua;
                musv(j,i)  = tissue_excite(i).mus;
                gv(j,i)    = tissue_excite(i).g;
            end
        case 2
            tissue_refpeak1 = makeTissueListRaman2(nm(j)); % also --> global tissue(1:Nt).s
            % 1/29/23 kdk: the above gets called for all values in nm
            Nt = length(tissue_refpeak1);
            % 1/29/23 kdk: why do we ned to copy these when we have them in tissue?
            for i=1:Nt
                muav(j,i)  = tissue_refpeak1(i).mua;
                musv(j,i)  = tissue_refpeak1(i).mus;
                gv(j,i)    = tissue_refpeak1(i).g;
            end
        case 3
            tissue_pHpeak1 = makeTissueListRaman2(nm(j)); % also --> global tissue(1:Nt).s
            % 1/29/23 kdk: the above gets called for all values in nm
            Nt = length(tissue_pHpeak1);
            % 1/29/23 kdk: why do we ned to copy these when we have them in tissue?
            for i=1:Nt
                muav(j,i)  = tissue_pHpeak1(i).mua;
                musv(j,i)  = tissue_pHpeak1(i).mus;
                gv(j,i)    = tissue_pHpeak1(i).g;
            end
        case 4
            tissue_refpeak2 = makeTissueListRaman2(nm(j)); % also --> global tissue(1:Nt).s
            % 1/29/23 kdk: the above gets called for all values in nm
            Nt = length(tissue_refpeak2);
            % 1/29/23 kdk: why do we ned to copy these when we have them in tissue?
            for i=1:Nt
                muav(j,i)  = tissue_refpeak2(i).mua;
                musv(j,i)  = tissue_refpeak2(i).mus;
                gv(j,i)    = tissue_refpeak2(i).g;
            end
        case 5
            tissue_pHpeak2 = makeTissueListRaman2(nm(j)); % also --> global tissue(1:Nt).s
            % 1/29/23 kdk: the above gets called for all values in nm
            Nt = length(tissue_pHpeak2);
            % 1/29/23 kdk: why do we ned to copy these when we have them in tissue?
            for i=1:Nt
                muav(j,i)  = tissue_pHpeak2(i).mua;
                musv(j,i)  = tissue_pHpeak2(i).mus;
                gv(j,i)    = tissue_pHpeak2(i).g;
            end
    end
    
end

% Specify Monte Carlo parameters    
Nx = Nbins;
Ny = Nbins;
Nz = Nbins;
dx = binsize;
dy = binsize;
dz = binsize;
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

if isinf(zfocus), zfocus = 1e12; end

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.
% 1/29/23 kdk: we only need one of these 

T = double(zeros(Ny,Nx,Nz)); 

% 1. Start with making everything skin and then adjust where 
T = T + 4;      % fill background with skin (dermis)

zsurf = 0.0100;  % position of air/skin surface

for iz=1:Nz % for every depth z(iz)

    % air
    if iz<=round(zsurf/dz)
        T(:,:,iz) = 2; 
    end

    % epidermis (60 um thick)
    if iz>round(zsurf/dz) & iz<=round((zsurf+0.0060)/dz)
        T(:,:,iz) = 5; 
    end
    
%     % Cylindrical SERS-active hydrogel @ xc, zc (2mm) , radius, oriented along y axis
%     xc      = 0;            % [cm], center of cylindrical sensor
%     zc      = 0.2;      	% [cm], center of cylindrical sensor 
%     sensorradius  = 0.0500; % [cm], cylindrical sensor radius
%     for ix=1:Nx
%             xd = x(ix) - xc;	% x distance from sensor center
%             zd = z(iz) - zc;   	% vessel, z distance from sensor center                
%             r  = sqrt(xd^2 + zd^2);	% r from vessel center
%             if r<=sensorradius     	% if r is within vessel
%                 for iy=1:Ny
%                     if iy > Ny/4 & iy < 3*Ny/4 % sensor is centered and
%                         % spans half the structure, which is 1cm out of 2cm
%                         % from y=-0.5 to y=+0.5 cm
%                         T(iy,ix,iz) = 10; % SERS-active hydrogel tissue type
%                         %fprintf('sensor at %d,%d,%d\n', ix,iy,iz);
%                     end 
%                 end %iy
%             end
%     end %ix 
% end % iz

% Disc SERS-active hydrogel @ xc, zc (2mm) , radius, oriented along y axis
    xc      = 0;        % [cm], center of cylindrical disc sensor
    yc      = 0;      	% [cm], center of cylindrical disc sensor 
    sensorradius  = 0.300; % [cm], cylindrical sensor radius
    
    izDepth = iz * 2/Nz;
    if izDepth >= 0.2 && izDepth <= 0.3 % sensor is centered and   
        for ix=1:Nx
            xd = x(ix) - xc;	% x distance from sensor center
            for iy= 1:Ny
                yd = y(iy) - yc;   	% y distance from sensor center                
                r  = sqrt(xd^2 + yd^2);	% r from vessel center
                if r<=sensorradius     	% if r is within vessel
                    % spans half the structure, which is 1cm out of 2cm
                    % from y=-0.5 to y=+0.5 cm
                    T(iy,ix,iz) = 10; % SERS-active hydrogel tissue type
                    if iz == 25
                        fprintf('sensor at %d,%d,%d\n', ix,iy,iz);
                    end
                end
            end % iy
        end %ix
    end % if izDepth
end % iz


%%
if SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',time_min);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,launchflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',xfocus);
        fprintf(fid,'%0.4f\n',yfocus);
        fprintf(fid,'%0.4f\n',zfocus);
        fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',waist);
        
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        fprintf(fid,'%d\n',Nw); % 1/29/23 kdk: number of wavelengths
        % 1/29/23 kdk: write out the wavelengths
        for j=1:Nw
            fprintf(fid,'%d\n', nm(j));
        end
        for j=1:Nw
            for i=1:Nt
                fprintf(fid,'%0.4f\n',muav(j,i));
                fprintf(fid,'%0.4f\n',musv(j,i));
                fprintf(fid,'%0.4f\n',gv(j,i));
            end
        end
    fclose(fid);

    %% write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

    toc
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
%% This shows cross-sectional view of cylindrical SERS sensor
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx

%% Look at structure of Tzy at ix=Nx/2
%% This shows side view of cylindrical SERS sensor
Tzyx = shiftdim(T,2); %Tyxz --> Tzyx
Tzy = Tzyx(:,:,Nx/2); % Tzy

%% Look at structure of Tyx at ix=Nx/2
%% This shows top view of cylindrical SERS sensor
Tyx = T(:,:,25); % Tyx

%%
figure(1); clf
sz = 12;  fz = 10; 
imagesc(x,z,Tzx,[1 Nt]);
hold on
set(gca,'fontsize',sz);
xlabel('x [cm]');
ylabel('z [cm]');
title('Tissue structure: front view','FontSize',18);
c = addLegend(Nt, Nz, dz, x, xmin, xmax, zmax, zmin, tissue_excite, mcflag, radius, zs, z);
%%
figure(2); clf
sz = 18;  
imagesc(y,z,Tzy,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('y [cm]')
ylabel('z [cm]')
title('Tissue structure: side view','FontSize',18);
c = addLegend(Nt, Nz, dz, x, xmin, xmax, zmax, zmin, tissue_excite, mcflag, radius, zs, z);
%%
figure(3); clf
sz = 18;  
imagesc(x,y,Tyx,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('y [cm]')
title('Tissue structure: top view at depth z=0.25 cm','FontSize',18);
c = addLegend(Nt, Ny, dy, x, xmin, xmax, ymax, ymin, tissue_excite, 4, radius, zs, z);



function c = addLegend(Nt, Nz, dz, x, xmin, xmax, zmax, zmin, tissue, mcflag, radius, zs, z)   
    fz = 10; 
    colorbar
    cmap = makecmapRaman(Nt);
    colormap(cmap)
    set(colorbar,'fontsize',1)
    % label colorbar
    zdiff = zmax-zmin;
    %%%
    if mcflag == 4
        for i=1:Nt
            yy = ((Nt-i)/(Nt-1)-0.5)*Nz*dz;
            text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
        end  
    else
        for i=1:Nt
            yy = (Nt-i)/(Nt-1)*Nz*dz;
            text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
        end 
    end
    
    text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
    axis equal image
    axis([xmin xmax zmin zmax])

    %%% draw launch
    N = 20; % # of beam rays drawn
    switch mcflag
        case 0 % uniform
              for i=0:N
                  % plot((-radius + 2*radius*i/N)*[1 1],[zs max(z)],'r-')
                  % Stop it at the sensor
                  % plot((-radius + 2*radius*i/N)*[1 1],[zs 0.2],'r-')
                  % Make it at 45 degrees
                  % factor(i) = 2*radius*i/N
                  % x = (-0.75 + factor(i))*[0.1 01];
                  x = [-0.2 0.07];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.202 0.06];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.204 0.05];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.206 0.04];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.208 0.03];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.21 0.02];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.212 0.01];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.214 0.0];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
                  x = [-0.216 -0.01];
                  z = [-0.1 0.2];
                  plot(x,z,'r-');
              end

        case 1 % Gaussian
            for i=0:N
                plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'r-')
            end

        case 2 % iso-point
            for i=1:N
                th = (i-1)/19*2*pi;
                xx = Nx/2*cos(th) + xs;
                zz = Nx/2*sin(th) + zs;
                plot([xs xx],[zs zz],'r-')
            end

        case 3 % rectangle
            zz = max(z);
            for i=1:N
                xx = -radius + 2*radius*i/20;
                plot([xx xx],[zs zz],'r-')
            end
        case 4 % do nothing
            fprintf("draw no excitation\n");
    end
    c = 1;
end


