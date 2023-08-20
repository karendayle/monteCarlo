% based on maketissue_18apr17.m
% maketissueRaman.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcraman.c.
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
%       makeTissueListRamanDisc.m
%
%   To use, 
%       1. Prepare makeTissueListRaman.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissueRaman.m.
%
%   Note: mcraman.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%  Dayle Kotturi, adapted for Raman Feb, 2020.  
%  Dayle Kotturi, adapted for a 6 mm diameter, 1 mm thick implant,
%                 centered at (0,0) and with its top surface at z=0.2 cm.
%                 Change laser beam to unfocused, collimated 4 mm diam
%                 May 23, 2022

clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

myname      = 'skindisc3min';% name for files: myname_T.bin, myname_H.mci  
time_min    = 3;      	% time duration of the simulation [min] <----- run time -----
nm          = 785;   	% desired wavelength of simulation
Nbins       = 200;    	% # of bins in each dimension of cube 
binsize     = 0.01; 	% size of each bin in [cm]. Increased by 20x for SERS-sensor scale

% Set Monte Carlo launch flags
mcflag      = 0;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                        % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
launchflag  = 1;        % 0 = let mcraman.c calculate launch trajectory
                        % 1 = manually set launch vector. 5/25/22 try this
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
radius      = 0.0200;   % 1/e radius of beam at tissue surface
waist       = 0.0200;  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
% ux0         = 0.7;      % trajectory projected onto x axis
% uy0         = 0.4;      % trajectory projected onto y axis
ux0         = 0.7;      % trajectory projected onto x axis 5/25/22 try this
uy0         = 0.;      % trajectory projected onto y axis 5/25/22 try this
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% 
% Prepare Monte Carlo 
%%%

% Create tissue properties
tissue = makeTissueListRaman(nm); % also --> global tissue(1:Nt).s
Nt = length(tissue);
for i=1:Nt
    muav(i)  = tissue(i).mua;
    musv(i)  = tissue(i).mus;
    gv(i)    = tissue(i).g;
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

T = double(zeros(Ny,Nx,Nz)); 

% 1. Start with making everything skin and then adjust where 
T = T + 4;      % fill background with skin (dermis)

zsurf = 0.1;  % position of air/skin surface

for iz=1:Nz % for every depth z(iz)

    % air
    if iz<=round(zsurf/dz)
        T(:,:,iz) = 1; 
    end

    % epidermis (60 um thick) -> update 5/23/22 to 400 um taking avg from
    % ref = https://www.umb.edu.pl/photo/pliki/progress-file/current_issue/progres_8.2/doi/223-228__niczyporuk_.pdf
    if iz>round(zsurf/dz) && iz<=round((zsurf+0.0400)/dz)
        T(:,:,iz) = 5; 
    end
    
    % Cylindrical SERS-active hydrogel @ xc, zc (2mm) , radius, oriented along y axis
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
                    T(ix,iy,iz) = 10; % SERS-active hydrogel tissue type
                    %fprintf('sensor at %d,%d,%d\n', ix,iy,iz);
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
        for i=1:Nt
            fprintf(fid,'%0.4f\n',muav(i));
            fprintf(fid,'%0.4f\n',musv(i));
            fprintf(fid,'%0.4f\n',gv(i));
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

% 5/23/22 changed setting of T from Tyxz to Txyz, for the sake of figure3,
% It means that Txzy and Tzyx are flipped now, about x=y
%% Look at structure of Tzx at iy=Ny/2
%% This shows cross-sectional view of cylindrical SERS sensor
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx

%% Look at structure of Tzy at ix=Nx/2
%% This shows side view of cylindrical SERS sensor
Tzyx = shiftdim(T,2); %Tyxz --> Tzyx
Tzy = Tzyx(:,:,Nx/2); % Tzy

%% Look at structure of Txy at iz=Nz/10
%% This shows top view of cylindrical SERS sensor
%Txyz = shiftdim(T,3);
Txy = T(:,:,Nz/10);

%%
figure(1); clf
sz = 12;  fz = 10; 
imagesc(x,z,Tzx,[1 Nt]);
hold on
set(gca,'fontsize',sz);
xlabel('x [cm]');
ylabel('z [cm]');
title('Tissue structure: front view','FontSize',18);
c = addFrontAndSideLegend(Nt, Nz, dz, x, xmin, xmax, zmax, zmin, tissue, mcflag, radius, zs, z);
%%
figure(2); clf
sz = 12;  
imagesc(y,z,Tzy,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('y [cm]')
ylabel('z [cm]')
title('Tissue structure: side view','FontSize',18);
c = addFrontAndSideLegend(Nt, Nz, dz, x, xmin, xmax, zmax, zmin, tissue, mcflag, radius, zs, z);
%%
figure(3); clf
sz = 12;  
imagesc(x,y,Txy,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('y [cm]')
title('Tissue structure: top view at z=0.25 cm','FontSize',18);
c = addTopLegend(Nt, Nz, dz, x, xmin, xmax, -1*ymax, ymax, tissue, mcflag, radius, zs, z);
hold on
drawGrid();

function c = addFrontAndSideLegend(Nt, Nz, dz, x, xmin, xmax, zmax, zmin, tissue, mcflag, radius, zs, z)   
    fz = 10; 
    colorbar
    cmap = makecmapRaman(Nt);
    colormap(cmap)
    set(colorbar,'fontsize',1)
    % label colorbar
    zdiff = zmax-zmin;
    %%%
    for i=1:Nt
        yy = (Nt-i)/(Nt-1)*Nz*dz;
        text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
    end
    
    text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
    axis equal image
    axis([xmin xmax zmin zmax])

    %%% draw launch
    N = 20; % # of beam rays drawn
    switch mcflag
        case 0 % uniform
            for i=0:N
                zDepth = max(z)/N + 0.005;
                plot((-radius + 2*radius*i/N)*[1 1],[zs zDepth],'r')
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
    end
    c = 1;
end

function c = addTopLegend(Nt, Nz, dz, x, xmin, xmax, ymin, ymax, tissue, mcflag, radius, zs, z)   
    fz = 10; 
    colorbar
    cmap = makecmapRaman(Nt);
    colormap(cmap)
    set(colorbar,'fontsize',1)
    % label colorbar
    ydiff = ymax-ymin;
    %%%
    for i=1:Nt
        yy = (Nt-i)/(Nt-1)*Nz*dz + ymin;
        text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
    end
    
    text(xmax,ymin - ydiff*0.06, 'Tissue types','fontsize',fz)
    axis equal image
    axis([xmin xmax ymin ymax])

    c = 1;
end

function d = drawGrid() 
    % vertical lines
    x = 0.;
    y = -0.5;
    for i = 1:11
        for j = 1:11
            plot([x x], [y y+1], 'Color', 'k');
            hold on      
        end
        x = x + 0.1;
    end

    % horizontal lines
    x = 0.;
    y = -0.5;
    for i = 1:11
        for j = 1:11
            plot([x x+1], [y y], 'Color', 'k');
            hold on     
        end
        y = y + 0.1;
    end
end
