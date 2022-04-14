%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo_nirfast_2d.m
% demo for nirfast simulation
% created on 2022.03.02 Jingjing Jiang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add paths for required packages
% define nirfast paths
% Nirfaster: GPU-fascilitated model
% Nirfast8: old CPU version
pathNirfaster = '/home/letizia/Documents/Nirfaster';
pathNirfast8 = '/home/letizia/Documents/nirfast8';
addpath(genpath(pathNirfast8));
pathPioneerIR = '/home/letizia/Documents/PioneerImageReconstruction';
addpath(genpath(pathPioneerIR));
addpath(genpath('/home/letizia/2D-simulation/demo_nirfast_2d/Analytical_solution'));
%% load mesh  / create mesh
mesh = load_mesh('/home/letizia/PioneerImageReconstruction/demo_nirfast_2d/mesh/Rectangle-stnd-mesh');
%mesh = load_mesh('/home/laninil/Documents/2D-simulation/demo_nirfast_2d/mesh2/Rectangle-stnd-mesh');
mesh_homo = mesh;
%% remove short source-detector distances
nsrc = length(mesh_homo.source.num);
ndet = length(mesh_homo.meas.num);
psrc = mesh_homo.source.coord;
pdet = mesh_homo.meas.coord;
link_all_selected = ones(ndet*nsrc,1);
thresh_d = 10;
for isrc = 1:nsrc
    dist_mesh(isrc).value = ...
        sqrt(sum((repmat(psrc(isrc,:),ndet,1) - pdet).^2,2));
    idx_delete = find(dist_mesh(isrc).value<thresh_d);
    link_all_selected((isrc-1)*ndet+idx_delete) = 0;
end
id_sel = find(link_all_selected);
mesh_homo.link = mesh.link(id_sel,:);

%% visualize mesh

h_mesh = figure()
% plot boundary nodes
mesh_1 = mesh_homo; 
plot(mesh_1.nodes(logical(mesh_1.bndvtx),1),mesh_1.nodes(logical(mesh_1.bndvtx),2),'c.');
hold on
% plot sources
plot(mesh_1.source.coord(:,1),mesh_1.source.coord(:,2),'rx','LineWidth',2,'MarkerSize',8);
plot(mesh_1.source.coord(1,1),mesh_1.source.coord(1,2),'mx','LineWidth',2,'MarkerSize',12);
% plot detectors
plot(mesh_1.meas.coord(:,1), mesh_1.meas.coord(:,2),'bo','LineWidth',2,'MarkerSize',8);
plot(mesh_1.meas.coord(1,1), mesh_1.meas.coord(1,2),'co','LineWidth',2,'MarkerSize',10);

%% add optical properties
wavList = {'690', '830'};
 
muas_bulk = [0.008 0.01];% Target values
mus_r_bulk = [0.3 0.3]; % Target values
%  mus_r_bulk - 0.54;
muas_si = muas_bulk .*5;
mus_r_si = mus_r_bulk;
% mus_r_bulk = 0.54;
iwav = 1;
val.mua=muas_bulk(iwav); %692
val.mus=mus_r_bulk(iwav);
 
val.ri=1.43;
mesh_homo = set_mesh(mesh_homo,0,val);

%% add blob
clear blob
blob.x= 6;
blob.y= 5;
blob.r= 3;
blob.mua=muas_si(iwav);
blob.mus= mus_r_si(iwav);
blob.region=2;
mesh_anom = add_blob(mesh_homo,blob);

blob.x= -6;
blob.y= 5;
blob.r= 3;
blob.mua=muas_si(iwav);
blob.mus= mus_r_si(iwav);
blob.region=2;
mesh_anom = add_blob(mesh_anom,blob);

id_inc = find(mesh_anom.mua>0.008);
length(mesh_anom.nodes(id_inc,:))
plotimage(mesh_anom, mesh_anom.mua);

%% reduce mesh node size
mesh_homo.nodes = mesh_homo.nodes(:,1:2);
mesh_anom.nodes = mesh_anom.nodes(:,1:2);
%% calculate moments
addpath(genpath(pathNirfaster))
% test up to 2nd moment order (0th, 1st, and 2nd)
max_order = 2;

% increase the absolute tolerance for moments as the values will be low.
% The relative tolerance stays the same.
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;

% GPU
tic;
[data_Moments] = femdata_stnd_TR_moments(mesh_homo, ...
    max_order,'field',[],solver_name_GPU,OPTIONS);
[data_Moments_anom] = femdata_stnd_TR_moments(mesh_anom, ...
    max_order,'field',[],solver_name_GPU,OPTIONS);
figure,
semilogy(data_Moments.moments)
hold on
semilogy(data_Moments_anom.moments)
%% reconstruction: moments
% reconstruction
tic
REG_THIKONOV = 1;
MESH_COARSE = [50 30]*0.5;
MAX_ITERATIONS = 20;
[meshRec_moments, pj_error_moments] = reconstruct_stnd_moments(mesh_homo, max_order, ...
     data_Moments_anom, 'mua',[], MESH_COARSE, MAX_ITERATIONS, REG_THIKONOV);
toc
plotimage(meshRec_moments, meshRec_moments.mua)

%% simulate frequency domain data 
% The unit of frequency in Nirfaster is Hz
% The unit of frequency in Nirfast8 is MHz
addpath(genpath(pathNirfaster))
frequencies = 100e6; % Hz
tic; [data_FD_NIRFAST_anom] = femdata_stnd_FD(mesh_anom,frequencies); 
figure,
subplot(121)
plot(data_FD_NIRFAST_anom.amplitude)
title('amplitude')
subplot(122)
plot(data_FD_NIRFAST_anom.phase)
title('phase')
%% reconstruction: frequency domain 
% data_FD_NIRFAST_anom = add_noise(data_FD_NIRFAST_anom, 2,2);
REG_THIKONOV =1;
tic
[meshRec_1F, pj_error_1F] = reconstruct_stnd_FD(mesh_homo,frequencies(1), ...
                             data_FD_NIRFAST_anom, 'mua',...
                            [],MESH_COARSE, MAX_ITERATIONS, REG_THIKONOV) 
toc
meshRec_1F
plotimage(mesh_anom, meshRec_1F.mua)

%% calculate temporal data
addpath(genpath(pathNirfaster))

% mesh_homo.link = mesh.link(find(link_all_selected),:);
% mesh_anom.link = mesh_homo.link;
cfg.tstart = 0;
cfg.tstep = 0.0488*1e-9/2;  
cfg.num_gates = 20*2;
cfg.tend = cfg.tstep* cfg.num_gates ;
cfg.timeGates = cfg.tstart+cfg.tstep:cfg.tstep:cfg.tend;

t = (1:steps)*dt - dt/2

tic; [data_TR_anom] = femdata_stnd_TR(mesh_anom,cfg.tend,cfg.tstep,'field', 'BiCGStab_GPU'); toc;   
tic; [data_TR_homo] = femdata_stnd_TR(mesh_homo,cfg.tend,cfg.tstep,'field', 'BiCGStab_GPU'); toc;  
% tic; [data_TR_anom] = femdata_stnd_TR(mesh_anom,cfg.tend,cfg.tstep,'field', 'BiCGStab_CPU'); toc;   
% tic; [data_TR_homo] = femdata_stnd_TR(mesh_homo,cfg.tend,cfg.tstep,'field', 'BiCGStab_CPU'); toc; 
%%  plot tpsfs 

% Phi=zeros(length(mesh_1.source.coord), length(mesh_1.meas.coord), length(cfg.timeGates));
% for i = 1:length(mesh_1.source.coord)
%     for j=1:length(mesh_1.meas.coord)
%         Phi(i,j,:) = tddiffusion(muas_si(iwav), ...
%             mus_r_si(iwav), 299792458000/val.ri, ...
%             0.493, ...
%             mesh_1.source.coord(i,:), mesh_1.meas.coord(j,:), cfg.timeGates);
%     end
% end
% Phi = sum(Phi, [1 2])                             %Sum over all source-detector pairs for every time gate, and plot that sum
% Phi = reshape(Phi, [1 length(cfg.timeGates)])

Phi=zeros(length(mesh_homo.link(:,2)), 1, length(cfg.timeGates));    %Sum only over the source-detector pairs with enough distance between them (those in mesh_homo.link)
for i = 1:length(mesh_homo.link(:,2))
    Phi(i,1,:) = tddiffusion(val.mua, ...
        val.mus, 299792458000/val.ri, ...
        val.ri, ...
        mesh_homo.source.coord(mesh_homo.link(i,1),:), mesh_homo.meas.coord(mesh_homo.link(i,2),:), ...
        cfg.timeGates);
end
Phi = sum(Phi, 1);
Phi = reshape(Phi, [1 length(cfg.timeGates)]);
% figure()
% plot(cfg.timeGates, Phi)

%Analytical solution from Nirfaster package
Phi2 = zeros(length(mesh_homo.link(:,2)), 1, length(cfg.timeGates));
for i = 1:length(mesh_homo.link(:,2))
    rsd2 = getdistance(mesh_homo.source.coord(mesh_homo.link(i,1),:), mesh_homo.meas.coord(mesh_homo.link(i,2),:));
%   rsd2 = sqrt(sum((mesh.source.coord(mesh_homo.link(i,1),:)-mesh.meas.coord(mesh_homo.link(i,2),:)).^2,2));     %gives same result! From femdata_stnd_TR line 251
    Phix = semi_infinite_TR(val.mua, val.mus, val.ri, rsd2, cfg.tend, cfg.tstep, 'EBC-Robin');
    Phi2(i,1,:) = Phix(:,2);
end
Phi2 = sum(Phi2, 1);
Phi2 = reshape(Phi2, [1, length(cfg.timeGates)]);

% figure()
% plot(cfg.timeGates, Phi2, 'DisplayName', 'getdistance')
% legend()

% Phi_single = semi_infinite_TR(val.mua, val.mus, val.ri, rsd(10), cfg.tend, cfg.tstep, 'EBC-Robin');

%mixed function from tddiffusion and boundary condition from semi_infinite_TD
Phi3 = zeros(length(mesh_homo.link(:,2)), 1, length(cfg.timeGates));        
for i = 1:length(mesh_homo.link(:,2))
    source = mesh_homo.source.coord(mesh_homo.link(i,1),:);
    detector = mesh_homo.meas.coord(mesh_homo.link(i,2),:);
    Phixx = semi_infinite_TR_tddiffusion(val.mua, val.mus, val.ri, source, detector, cfg.tend, cfg.tstep, 'EBC-Robin');
    Phi3(i,1,:) = Phixx(:,2);
end
Phi3 = sum(Phi3, 1);
Phi3 = reshape(Phi3, [1, length(cfg.timeGates)]);

itpsf = 6;  %was 6      %-> to select a signle source-detector pair!
%otherwise sum over all pairs for each time step:
Phi_homo = sum(data_TR_homo.tpsf, 1);
Phi_homo = reshape(Phi_homo, [1, length(cfg.timeGates)]);

Phi_anom = sum(data_TR_anom.tpsf, 1);
Phi_anom = reshape(Phi_anom, [1, length(cfg.timeGates)]);

max_Phi = max(Phi);
max_Phi2 = max(Phi2);
max_Phi3 = max(Phi3);
max_homo = max(data_TR_homo.tpsf(itpsf,:));
max_homo_all = max(Phi_homo);
max_anom = max(data_TR_anom.tpsf(itpsf,:));
max_anom_all = max(Phi_anom);
%max_single = max(Phi_single(:,2));

h_psf=figure()
% plot(cfg.timeGates*1e9, data_TR_homo.tpsf(itpsf,:)/max_homo, 'DisplayName', 'Homogeneous')
% hold on
plot(cfg.timeGates*1e9, Phi_homo/max_homo_all, 'DisplayName', 'Homogeneous (sum)')
hold on
% plot(cfg.timeGates*1e9, data_TR_anom.tpsf(itpsf,:)/max_anom,'DisplayName', 'Heterogeneous')
% hold on
plot(cfg.timeGates*1e9, Phi_anom/max_anom_all, 'DisplayName', 'Heterogeneous (sum)')
hold on
plot(cfg.timeGates*1e9, Phi/max_Phi, 'DisplayName', 'Analytical solution (sum)')
hold on
plot(cfg.timeGates*1e9, Phi2/max_Phi2, 'DisplayName', 'Analytical solution Nirfaster (sum)')
hold on
plot(cfg.timeGates*1e9, Phi3/max_Phi3, 'DisplayName', 'Analytical solution mix (sum)')
%hold on
%plot(cfg.timeGates, Phi_single(:,2)/max_single, 'DisplayName', 'Analytical solution Nirfaster')
legend();
xlabel('time [ps]')

%try normalization
figure()
% plot(cfg.timeGates*1e9, data_TR_homo.tpsf(itpsf,:), 'DisplayName', 'Homogeneous')
% hold on
plot(cfg.timeGates*1e9, Phi_homo, 'DisplayName', 'Homogeneous (sum)')
hold on
% plot(cfg.timeGates*1e9, data_TR_anom.tpsf(itpsf,:),'DisplayName', 'Heterogeneous')
% hold on
plot(cfg.timeGates*1e9, Phi_anom, 'DisplayName', 'Heterogeneous (sum)')
hold on
plot(cfg.timeGates*1e9, Phi*cfg.tstep, 'DisplayName', 'Analytical solution (sum)')
hold on
plot(cfg.timeGates*1e9, Phi2, 'DisplayName', 'Analytical solution Nirfaster (sum)')
hold on
plot(cfg.timeGates*1e9, Phi3, 'DisplayName', 'Analytical solution mix (sum)')
legend();
xlabel('time [ps]')

%% compare different source-detector pairs

pair_set = [mesh_homo.link(10,1:2); mesh_homo.link(15,1:2); mesh_homo.link(18,1:2); mesh_homo.link(21,1:2); mesh_homo.link(27,1:2); mesh_homo.link(22,1:2); mesh_homo.link(25,1:2); mesh_homo.link(40,1:2); mesh_homo.link(45,1:2); mesh_homo.link(1,1:2)];

%in time domain
figure()
for i = 1:length(pair_set)
    source = mesh_homo.source.coord(pair_set(i,1),:);
    detector = mesh_homo.meas.coord(pair_set(i,2),:);
    Phi_single = semi_infinite_TR_tddiffusion(val.mua, val.mus, val.ri, source, detector, cfg.tend, cfg.tstep, 'EBC-Robin');
    Phi_single_max = max(Phi_single(:,2));
    distance = getdistance(mesh_homo.source.coord(pair_set(i,1),:), mesh_homo.meas.coord(pair_set(i,2),:));
    %plot(cfg.timeGates, Phi_single(:,2), 'DisplayName', sprintf('Pair %d', pair_set(i,1), pair_set(i,2)));
    %plot(cfg.timeGates, Phi_single(:,2)/Phi_single_max, 'DisplayName', sprintf('%d', pair_set(i,:)));
    plot(cfg.timeGates, Phi_single(:,2), 'DisplayName', sprintf('%d', distance));
    %plot(cfg.timeGates, Phi_single(:,2)/Phi_single_max, 'DisplayName', sprintf('%d', distance));
    hold on;
    legend('-DynamicLegend');
end
%plot(cfg.timeGates, data_TR_homo.tpsf(itpsf,:)/max_homo, 'DisplayName', 'Homogeneous', 'color', 'r')
hold off;

%compare analytical and measured for pair 6-2
source = mesh_homo.source.coord(6,:);
detector = mesh_homo.meas.coord(2,:);
Phi_single = semi_infinite_TR_tddiffusion(val.mua, val.mus, val.ri, source, detector, cfg.tend, cfg.tstep, 'EBC-Robin');

figure()
plot(cfg.timeGates*1e9, data_TR_homo.tpsf(30,:), 'DisplayName', 'Homogeneous')
hold on
plot(cfg.timeGates*1e9, data_TR_anom.tpsf(30,:), 'DisplayName', 'Heterogeneous')
hold on
plot(cfg.timeGates*1e9, Phi_single(:,2)*max(data_TR_homo.tpsf(30,:))/max(Phi_single(:,2)), 'DisplayName', 'Analytical')
hold off
legend()

%in frequency domain
rsd = sqrt(sum((mesh.source.coord(mesh.link(:,1),:) - mesh.meas.coord(mesh.link(:,2),:)).^2,2));
Phi_FD = semi_infinite_FD(val.mua, val.mus, val.ri, frequencies, rsd)
figure()
plot(Phi_FD)

for i = 1:length(pair_set)
    rsd = getdistance(mesh_homo.source.coord(mesh_homo.link(i,1),:), mesh_homo.meas.coord(mesh_homo.link(i,2),:));
%   rsd = sqrt(sum((mesh.source.coord(mesh_homo.link(i,1),:)-mesh.meas.coord(mesh_homo.link(i,2),:)).^2,2));     %gives same result! From femdata_stnd_TR line 251
    Phi_FD = semi_infinite_FD(val.mua, val.mus, val.ri, frequencies, rsd)
end
% figure()
% plot(Phi_FD)

% mesh_pair = mesh_homo;
% for i = 1:length(pair_set)
%     mesh_pair.source = 
%     [data_FD_NIRFAST_homo] = femdata_stnd_FD(mesh_homo,frequencies);
%     
%     figure()
%     subplot(121)
%     plot(data_FD_NIRFAST_homo.amplitude)
%     title('amplitude')
%     subplot(122)
%     plot(data_FD_NIRFAST_homo.phase)
%     title('phase')
% end

%% image reconstruction: Time domain
t_cfg.range = cfg.tend;
t_cfg.step = cfg.tstep;
t_cfg.num_gates = cfg.num_gates;
tic
[meshRec_TD, pj_error_TG] = reconstruct_stnd_TD(mesh_homo,t_cfg, ...
        data_TR_anom, 'mua',[], MESH_COARSE, MAX_ITERATIONS, REG_THIKONOV)
     toc
    plotimage(meshRec_TD, meshRec_TD.mua)


