clearvars
close all

%% set path

addpath(genpath('...\fcn'))
addpath(genpath('...\ext'))

%% set directories 

data_dir = '...\mat';
out_dir = '...\out';
comm_dir = '...\out\optimized_comm';

%% load data

load(fullfile(data_dir, 'grandaverage_HCP'))
load(fullfile(data_dir, 'yeo7_200'));
load(fullfile(out_dir, 'TSE_HCP_curve'))
load(fullfile(comm_dir, 'anneal_comm_100iter.mat'))
load(fullfile(mrcc_dir, 'MRCC_out.mat'))

FCr = FC(randorder, randorder);
yeo7r = yeo7(randorder, 1);
[~, back2ord] = sort(randorder);

nc = size(C_out_all, 3);

yeo_lab = {'Vis', 'SM', 'DAN', 'VAN', 'Lim', 'FP', 'DMN'};

%% set up parc plot

[surfStruct, annotMap, annotName] = setup_parcplot_schaefer(200);

%% schematic of the algorithm

nc = 5;
it = 1;

seed = randi([1 nc],N,1);
part = squeeze(C_out_all(back2ord,it,nc));

mod_size = histcounts(part, [0.5 : 1 : nc+0.5]);
mod_size_seed = histcounts(seed, [0.5 : 1 : nc+0.5]);

for m=1:nc
    ind = find(part==m);
    TCmod(m) = calcI_alt(FC(ind,ind));
    ind_seed = find(seed==m);
    TCmod_seed(m) = calcI_alt(FC(ind_seed,ind_seed));
end

ff = figure;
l1 = plot(1:200, Ilvl, 'Color', rgb('SteelBlue'), 'LineWidth', 1.1);
hold on 
l2 = plot(1:200, Ilvl_max, 'Color', rgb('midnightblue'), 'LineWidth', 1.1);
s1 = scatter(mod_size, TCmod, 30, 'MarkerFaceColor', rgb('Crimson'), 'MarkerEdgeColor', 'none');
s2 = scatter(mod_size_seed, TCmod_seed, 15, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0]);
axis square; grid on;
xlabel('subset size', 'FontSize', 15)
ylabel('TC (normalized)', 'FontSize', 15)
set(gca, 'FontSize', 15)
legend([l1, l2, s1, s2], 'avg. integration', 'max. integration', 'TC opt. modules', 'TC seed modules',...
    'Location','northwest')

% zooom in
ff = figure;
hold on 
l1 = plot(1:200, Ilvl, 'Color', rgb('SteelBlue'), 'LineWidth', 1.1);
l2 = plot(1:200, Ilvl_max, 'Color', rgb('midnightblue'), 'LineWidth', 1.1);
xlim([27 53])
axis square; grid on; box on
xlabel('subset size', 'FontSize', 15)
ylabel('TC (normalized)', 'FontSize', 15)
set(gca, 'FontSize', 15, 'XTick', 30:5:50) 
s1 = scatter(mod_size, TCmod, 60, 'MarkerFaceColor', rgb('Crimson'), 'MarkerEdgeColor', 'none');
for ii=1:nc
    line([mod_size(ii) mod_size(ii)], [Ilvl(mod_size(ii)) TCmod(ii)],...
        'LineStyle', ':', 'Color', rgb('Crimson'), 'LineWidth', 1.3)
    line([mod_size_seed(ii) mod_size_seed(ii)], [Ilvl(mod_size_seed(ii)) TCmod_seed(ii)],...
        'LineStyle', ':', 'Color', [0 0 0], 'LineWidth', 1.3)
end
s2 = scatter(mod_size_seed, TCmod_seed, 25, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0]);

% plot on cortex
col_comm = cbrewer('qual', 'Set3', 7, 'pchip');
for i=1:N
    col_nodes(i,:) = col_comm(part(i), :);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)

for i=1:N
    col_nodes(i,:) = col_comm(seed(i), :);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)

% similarity between output partitions
for ii=2:size(C_out_all,3)
    jj = 1;
    for it= 1:iter
        for it2= it+1:iter
            ami(ii,jj) = AMI(C_out_all(:,it,ii), C_out_all(:,it2,ii));
            [vi(ii,jj), nmi(ii,jj)] = partition_distance(C_out_all(:,it,ii), C_out_all(:,it2,ii));
            jj = jj+1;
        end
    end
end
vi_comp = 1-vi;


colors = cbrewer('qual', 'Paired', 6, 'pchip');
col2 = colors(2,:);
ff= figure;
hold on
line([1 13], [0 0], 'LineWidth', 1.1, 'LineStyle', '--', 'Color', [0.1 0.1 0.1])
line([1 13], [1 1], 'LineWidth', 1.1, 'LineStyle', '--', 'Color', [0.4 0.4 0.4])
for ii= 2:size(sim,1)
    r1 = 0.9 + ii-1;
    r2 = 1.1 + ii-1;
    jit = r1 + (r2-r1).*rand(1,size(sim,2));
    s1 = scatter(jit, sim(ii,:), 25, 'MarkerEdgeColor', [0.7 0.7 0.7],...
        'MarkerFaceColor', [0.9 0.9 0.9]);
    hold on
    x1 = 0.75 + ii-1;
    y1 = prctile(sim(ii,:), 25);
    y2 = prctile(sim(ii,:), 75);
    rectangle('Position', [x1 y1 0.5 y2-y1], 'LineWidth', 1.5,...
        'EdgeColor', col2)
    ly1 = prctile(sim(ii,:), 5);
    ly2 = prctile(sim(ii,:), 95);
    line([ii ii], [ly1 y1], 'Color', col2, 'LineWidth', 1.5)
    line([ii ii], [y2 ly2], 'Color', col2, 'LineWidth', 1.5)
    ym = mean(sim(ii,:));
    line([x1 x1+0.5], [ym ym], 'Color', col2, 'LineWidth', 1.5)
end
xlim([1 13])
ylim([-0.02 1.02])
axis square; box on; grid on;
set(gca, 'XTick', 2:12, 'XTickLabel', {'2','','4','','6','','8','','10','','12'},...
    'XTickLabelRotation', 0, 'FontSize', 15)
xlabel('Number of modules')
ylabel('VI')


%% relation with Q


% compute modularity on the partitions obtained with TCmods_anneal
for jj= 1:nc
    for it= 1:iter
        q_anneal(it,jj) = modularity_signed(FCr, C_out_all(:,it,jj), 1);
    end
end
clear jj

% compute TCscore on MRCC partitions
for jj= 1:size(ciall,2)
    tmp_nc = length(unique(ciall(:,jj)));
    TCmod = zeros(1, tmp_nc);
    mod_size = histcounts(ciall(:,jj), [0.5:1:tmp_nc+0.5]);   
    for m= 1:tmp_nc
        ff = find(ciall(:,jj)==m);
        TCmod(m) = calcI_alt(FC(ff,ff));
    end
    TCscore_mrcc(jj) = mean(TCmod - Ilvl(mod_size));
    q_mrcc(jj) = modularity_signed(FC, ciall(:,jj), gamma(jj));
end
clear TCmod

% plot qmrcc vs tcscoremrcc

comm_lab = unique(max(ciall));
nc_mrcc = max(ciall);
colors = crameri('batlow', length(comm_lab));
ff = figure;
hold on
for i= 1:length(comm_lab)
    ind = find(nc_mrcc == comm_lab(i));
    scatter(q_mrcc(ind), TCscore_mrcc(ind), 30,...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor',colors(i,:))
%     plot(q_mrcc(ind), TCscore_mrcc(ind), '.', 'Color', colors(i,:));
end
box on; grid on; axis square
xlabel('Q')
ylabel('TC_{score}')
set(gca, 'xlim', [0.1 0.95], 'ylim', [0 7], 'XTick', 0.2:0.2:0.9,...
    'YTick', 1:1:7, 'FontSize', 15, 'LineWidth', 0.7)
[rho_qtc, pval_qtc] = corr(q_mrcc', TCscore_mrcc')
text(0.21, 6.3, 'r = 0.91, p_{val}<10^{-15}', 'FontSize', 14, 'Color', [0.4 0.4 0.4])

ff= figure;
imagesc(1:49)
colormap(colors)
set(gcf, 'Position', [360 198 240 420])
set(gca, 'DataAspectRatio', [1 0.3 1], 'FontSize', 17, 'XTick', 10:10:40,...
    'YTick', [])
title('number of modules', 'FontSize', 15)


% plot qmrcc vs tcmrcc for every resolution
colscatt = rgb('steelblue'); 
for i= 22 %1:12 %length(comm_lab)
    ff = figure;
    ind = find(nc_mrcc == comm_lab(i));
    scatter(q_mrcc(ind), TCscore_mrcc(ind), 60, 'MarkerFaceColor', colscatt,...
        'MarkerEdgeColor', [0 0 0], 'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha', 0.95);
    axis square; box on
    xlabel('Q'); ylabel('TC_{score}')
    set(gca, 'FontSize', 20, 'LineWidth', 0.7)
    [rrr, ppp] = corr(q_mrcc(ind)', TCscore_mrcc(ind)');
    SW = [min(xlim) min(ylim)]+[diff(xlim) diff(ylim)]*0.05;
    text(SW(1), SW(2), sprintf('r = %g', round(rrr, 2)), 'FontSize', 16, 'Color', [0.4 0.4 0.4],...
        'VerticalAlignment','bottom', 'HorizontalAlignment','left')
    title(sprintf('M = %g', comm_lab(i)), 'FontSize',20)
end
clear rrr ppp ind ff

% plot correlation for every resolution
col_sig = rgb('darkorange');
col_nonsig = rgb('lightslategray');
ff = figure;
hold on
for i= 1:length(comm_lab)
    ind = find(nc_mrcc == comm_lab(i));
    [rrr(i), ppp(i)] = corr(q_mrcc(ind)', TCscore_mrcc(ind)');
    if ppp(i)<0.05
        scatter(comm_lab(i), rrr(i), 60, 'MarkerFaceColor', col_sig, 'MarkerEdgeColor', 'none');
    else
        scatter(comm_lab(i), rrr(i), 60, 'MarkerFaceColor', col_nonsig, 'MarkerEdgeColor', 'none');
    end
end
line([0 50], [0 0], 'Color', [0 0 0], 'LineWidth', 1.2)
xlim([0 50])
ylim([-1.1 0.2])
box on; grid on; axis square
set(gca, 'FontSize', 15)
xlabel('Number of modules')
ylabel('corr ( Q ; TC_{score} )')


% plot 2 extreme brains for nc=7

col_comm = cbrewer('qual', 'Paired', 10, 'pchip');
col_comm = col_comm([4 10 2 1 3 6 7], :);

ind7 = find(nc_mrcc == 7);

[maxTC7, pos_maxTC7] = max(TCscore_mrcc(ind7));
max_part = ciall(:,ind7(pos_maxTC7));
for i=1:N
    col_nodes(i,:) = col_comm(max_part(i,1),:);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)

[minTC7, pos_minTC7] = min(TCscore_mrcc(ind7));
min_part = ciall(:,ind7(pos_minTC7));
min_part_al = get_AlignedPartitions(max_part, min_part);
for i=1:N
    col_nodes(i,:) = col_comm(min_part_al(i,1),:);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)

% plot 2 extreme brains for nc=5

col_comm = cbrewer('qual', 'Paired', 10, 'pchip');
col_comm = col_comm([1 7 3 10 6], :);

ind5 = find(nc_mrcc == 5);

[maxTC5, pos_maxTC5] = max(TCscore_mrcc(ind5));
max_part = ciall(:,ind5(pos_maxTC5));
for i=1:N
    col_nodes(i,:) = col_comm(max_part(i,1),:);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)

[minTC5, pos_minTC5] = min(TCscore_mrcc(ind5));
min_part = ciall(:,ind5(pos_minTC5));
min_part_al = get_AlignedPartitions(max_part, min_part);
for i=1:N
    col_nodes(i,:) = col_comm(min_part_al(i,1),:);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)

% symmetry of partitions and dimension of modules 
hem1 = 1:100; hem2 = 101:200;
for i= 1:length(comm_lab)
    ind = find(nc_mrcc == comm_lab(i));
    [~, tmp_pos_maxTC] = max(TCscore_mrcc(ind)); 
    [~, tmp_pos_minTC] = min(TCscore_mrcc(ind));
    max_ciall = ciall(:,ind(tmp_pos_maxTC)); 
    min_ciall = ciall(:,ind(tmp_pos_minTC));
    for jj= 1:comm_lab(i)
        pos_max = find(max_ciall==jj);
        unb_max(jj) =...
            abs(length(find(ismember(pos_max, hem1))) - length(find(ismember(pos_max, hem2))))/length(pos_max);
        pos_min = find(min_ciall==jj);
        unb_min(jj) =...
            abs(length(find(ismember(pos_min, hem1))) - length(find(ismember(pos_min, hem2))))/length(pos_min);
    end
    simm_max(i) = mean(unb_max); clear unb_max
    simm_min(i) =  mean(unb_min); clear unb_min
    mod_size_max(i) = std(histcounts(max_ciall, [0.5:1:comm_lab(i)+0.5]));
    mod_size_min(i) = std(histcounts(min_ciall, [0.5:1:comm_lab(i)+0.5]));
end

colors = cbrewer('qual', 'Paired', 10, 'pchip');
col_max = colors(2,:); col_min = colors(1,:);

ff = figure; hold on;
for i= 1:length(comm_lab)
    line([comm_lab(i) comm_lab(i)], [mod_size_min(i) mod_size_max(i)],...
        'LineStyle', '--', 'Color', [0.3 0.3 0.3])
end
s1 = scatter(comm_lab, mod_size_min, 50, 'MarkerFaceColor', col_min, 'MarkerEdgeColor', 'none');
s2 = scatter(comm_lab, mod_size_max, 50, 'MarkerFaceColor', col_max, 'MarkerEdgeColor', 'none');
xlim([0 50])
box on; grid on; axis square
set(gca, 'FontSize', 15)
xlabel('Number of modules')
ylabel('Variance of modules size')
legend([s2, s1], 'max TCscore', 'min  TCscore')

ff = figure; hold on;
for i= 1:length(comm_lab)
    line([comm_lab(i) comm_lab(i)], [simm_min(i) simm_max(i)],...
        'LineStyle', '--', 'Color', [0.3 0.3 0.3])
end
s1 = scatter(comm_lab, simm_min, 50, 'MarkerFaceColor', col_min, 'MarkerEdgeColor', 'none');
s2 = scatter(comm_lab, simm_max, 50, 'MarkerFaceColor', col_max, 'MarkerEdgeColor', 'none');
xlim([0 50])
box on; grid on; axis square
set(gca, 'FontSize', 15)
xlabel('Number of modules')
ylabel('Symmetry of modules')
legend([s2, s1], 'max TCscore', 'min  TCscore', 'Location', 'southeast')


% count negative edges within modules
for i= 1:length(comm_lab)
    ind = find(nc_mrcc == comm_lab(i));
    [~, tmp_pos_maxTC] = max(TCscore_mrcc(ind)); 
    [~, tmp_pos_minTC] = min(TCscore_mrcc(ind));
    max_ciall = ciall(:,ind(tmp_pos_maxTC)); 
    min_ciall = ciall(:,ind(tmp_pos_minTC));
    for jj= 1:comm_lab(i)
        posmin = find(min_ciall==jj);
        count_min(jj) = length(find(FC(posmin,posmin)<0));
        posmax = find(max_ciall==jj);
        count_max(jj) = length(find(FC(posmax,posmax)<0));
    end
    neg_min(i) = mean(count_min);
    neg_max(i) = mean(count_max);
    clear count_min count_max
end

ff = figure; hold on;
for i= 1:length(comm_lab)
    line([comm_lab(i) comm_lab(i)], [neg_min(i) neg_max(i)],...
        'LineStyle', '--', 'Color', [0.3 0.3 0.3])
end
s1 = scatter(comm_lab, neg_min, 50, 'MarkerFaceColor', col_min, 'MarkerEdgeColor', 'none');
s2 = scatter(comm_lab, neg_max, 50, 'MarkerFaceColor', col_max, 'MarkerEdgeColor', 'none');
xlim([0 50])
box on; grid on; axis square
set(gca, 'FontSize', 15)
xlabel('Number of modules')
ylabel('Negative edges within module')
legend([s2, s1], 'max TCscore', 'min  TCscore')

%% relation with YEO systems

for ii= 2:nc
    [vi_yeo(ii,:), nmi_yeo(ii,:)] = partition_distance(C_out_all(:,:,ii), yeo7r);
    for it= 1:iter
        ami_yeo(ii,it) = AMI(C_out_all(:,it,ii), yeo7r);
    end
    [r_yeo(ii), pval_yeo(ii)] = corr(cost_out_all(:,ii), ami_yeo(ii,:)');
end

% plot similarity to yeo

colors = cbrewer('qual', 'Paired', 6, 'pchip');
col2 = colors(2,:);
ff= figure;
hold on
for ii= 2:nc
    r1 = 0.9 + ii-1;
    r2 = 1.1 + ii-1;
    jit = r1 + (r2-r1).*rand(1,iter);
    s1 = scatter(jit, ami_yeo(ii,:), 25, 'MarkerEdgeColor', [0.7 0.7 0.7],...
        'MarkerFaceColor', [0.9 0.9 0.9]);
    hold on
    x1 = 0.75 + ii-1;
    y1 = prctile(ami_yeo(ii,:), 25);
    y2 = prctile(ami_yeo(ii,:), 75);
    rectangle('Position', [x1 y1 0.5 y2-y1], 'LineWidth', 1.5,...
        'EdgeColor', col2)
    ly1 = prctile(ami_yeo(ii,:), 5);
    ly2 = prctile(ami_yeo(ii,:), 95);
    line([ii ii], [ly1 y1], 'Color', col2, 'LineWidth', 1.5)
    line([ii ii], [y2 ly2], 'Color', col2, 'LineWidth', 1.5)
    ym = mean(ami_yeo(ii,:));
    line([x1 x1+0.5], [ym ym], 'Color', col2, 'LineWidth', 1.5)
end
xlim([1 13])
ylim([0.15 0.7])
axis square
set(gca, 'XTick', 2:12, 'XTickLabelRotation', 0, 'FontSize', 15)
box on; grid on
xlabel('Number of modules')
ylabel('AMI (opt.part. ; Yeo systems)')

% plot YEO - consensus at 7 clusters - switching nodes

col_comm = cbrewer('qual', 'Set3', 10, 'pchip');

freq7 = find_FrequentPartition(squeeze(C_out_all(:,:,7)));
freq7 = freq7(back2ord);
freq7 = get_AlignedPartitions(yeo7, freq7);
for i=1:N
    col_yeo(i,:) = col_comm(yeo7(i,1),:);
    col_tcpart(i,:) = col_comm(freq7(i,1),:);
end

ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_yeo)

ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_tcpart)

col_switching = cbrewer('seq', 'Greys', 100, 'pchip'); col_switching = col_switching(1:75, :);
diff_part = zeros(N,1);
for it= 1:iter
    part7 = get_AlignedPartitions(yeo7, C_out_all(back2ord,it,7));
    diff_part = diff_part + logical(yeo7-part7);
end
col_switch = get_proportionalColors(diff_part(logical(diff_part)), col_switching);
col_zeros = repmat(col_switching(1,:), length(find(diff_part==0)), 1);
col2plot(logical(diff_part==0),:) = col_zeros;
col2plot(logical(diff_part),:) = col_switch;
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col2plot)

ff= figure;
imagesc(1:75)
colormap(col_switching)
set(gcf, 'Position', [360 198 240 420])
set(gca, 'DataAspectRatio', [1 0.15 1], 'FontSize', 17, 'YTick', [],...
    'XTick', [1 75], 'XTickLabel', [min(diff_part) max(diff_part)])


for ii=1:7
    diff_sys(ii) = sum(diff_part(find(yeo7==ii),1));
end
diff_sys = diff_sys./max(diff_sys);
col_bar = cbrewer('seq', 'Greys', 10, 'pchip'); col_bar = col_bar(1:7,:);
col_bar = get_proportionalColors(diff_sys', col_bar);
ff = figure;
bb = bar(diff_sys);
bb.FaceColor = 'flat';
bb.CData = col_bar;
set(gca, 'XTick', 1:7, 'XTickLabel', yeo_lab, 'Fontsize', 17)
ylim([0 1.1])
ylabel('freq.')

% make matrix to see how the modules split - use centroid
colscatt = col_comm(1:7,:);
colscatt(2,:) = [1 1 102/255]; %rgb('gold');
ff= figure;
hold on
for ii= 1:7
    pos = find(yeo7==ii);
    xc = 0.8*rand(1,length(pos)) + ii-1;
    laby = unique(freq7(pos));
    yc = [];
    for jj= 1:length(laby)
        yc = cat(2, yc, (0.8*rand(1,length(find(freq7(pos)==laby(jj)))) + laby(jj)-1));
    end
    scatter(xc, yc, 30, 'MarkerEdgeColor', rgb('darkslategray'),...
        'MarkerFaceColor', colscatt(ii,:))
    if ii<7
        line([ii ii], [0 7], 'LineWidth', 0.4, 'Color', rgb('lightgray'))
        line([0 7], [ii ii], 'LineWidth', 0.4, 'Color', rgb('lightgray'))
    end
end
axis square
box on
xlim([0 7])
ylim([0 7])
set(gca, 'XTick', 0.5:1:6.5, 'XTickLabel', yeo_lab, 'YTick', 0.5:1:6.5,...
    'YTickLabel', {'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7'},...
    'FontSize', 15)

% make matrix to see how the modules split - use all iterations
clear pos
for ii= 1:7
    pos{ii} = find(yeo7==ii);
    dimyeo(ii) = length(pos{ii});
end
for it= 1:iter
    part7 = get_AlignedPartitions(yeo7, C_out_all(back2ord,it,7));
    for jj= 1:7
        labp7 = unique(part7(pos{jj}));
        for kk=1:length(labp7)
            mat_iter(jj,labp7(kk),it) = length(find(part7(pos{jj})==labp7(kk)))/dimyeo(jj);
        end
    end
end
mat_iter_av = mean(mat_iter,3);
colors = cbrewer('seq', 'Greys', 100, 'pchip');
colprop = get_proportionalColors(mat_iter_av, colors);

ff= figure;
hold on
for ii= 1:7
    for jj= 1:7
        rectangle('Position', [ii-1 jj-1 1 1], 'FaceColor', colprop(ii,:,jj),...
            'LineWidth', 0.5, 'EdgeColor', rgb('silver'))
    end
end
axis square;
xlim([0 7])
ylim([0 7])
set(gca, 'XTick', 0.5:1:6.5, 'XTickLabel', yeo_lab, 'YTick', 0.5:1:6.5,...
    'YTickLabel', {'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7'},...
    'FontSize', 15)

ff= figure;
imagesc([1:100]')
colormap(flipud(colors))
set(gcf, 'Position', [360 198 420 240])
set(gca, 'DataAspectRatio', [0.1 1 1], 'FontSize', 17, 'YTick', [2 99],...
    'YTickLabel', [1 0], 'XTick', [], 'YAxisLocation','right')


%% TC scores vs number of clusters

% plot TC across number of modules

colors = cbrewer('qual', 'Paired', 10, 'pchip');
col1 = colors(1,:);
col2 = colors(2,:);% rgb('Navy');%

nc = 12;
num_comm = 2:nc;

ff= figure;
hold on
for ii= num_comm
    r1 = 0.9 + ii-1;
    r2 = 1.1 + ii-1;
    jit = r1 + (r2-r1).*rand(1,iter);
    s1 = scatter(jit, cost_out_all(:,ii), 25, 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', col1);
    hold on
    x1 = 0.75 + ii-1;
    y1 = prctile(cost_out_all(:,ii), 25);
    y2 = prctile(cost_out_all(:,ii), 75);
    rectangle('Position', [x1 y1 0.5 y2-y1], 'LineWidth', 1.5,...
        'EdgeColor', col2)
    ly1 = prctile(cost_out_all(:,ii), 5);
    ly2 = prctile(cost_out_all(:,ii), 95);
    line([ii ii], [ly1 y1], 'Color', col2, 'LineWidth', 1.5)
    line([ii ii], [y2 ly2], 'Color', col2, 'LineWidth', 1.5)
    ym = mean(cost_out_all(:,ii));
    line([x1 x1+0.5], [ym ym], 'Color', col2, 'LineWidth', 1.5)
end
xlim([1 13])
ylim([3 8.5])
set(gca, 'XTick', 2:12, 'XTickLabel', {'2','3','4','5','6','7','8','9','10','11','12'},...
    'YTick', 3.5:0.5:8, 'YTickLabel', {'','4','','5','','6','','7','','8'}, 'FontSize', 15)
box on; grid on
xlabel('Number of modules')
ylabel('TC_{score}')

% compute oinfo rand
for ii= num_comm
    for it= 1:iter
        for s=1:100
            nod_rand = randperm(N, ii);
            [Orand(s), ~, ~] = calcO_logdet(FCr(nod_rand,nod_rand));
        end
        Obc_rand(it,ii) = mean(Orand);
    end
    clear Orand
end

% plot oinfo between
col1_obs = colors(9,:); col2_obs = colors(10,:);
col1_null = rgb('lightgrey'); col2_null = rgb('darkslategrey');
ff= figure;
hold on
for bb= num_comm(2:end)
    % null
    data2plot = Obc_rand(:,bb);
    r1 = 0.9 + bb-1;
    r2 = 1.1 + bb-1;
    jit = r1 + (r2-r1).*rand(1,100);
    s2 = scatter(jit, data2plot, 25, 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', col1_null);
    hold on
    x1 = 0.75 + bb-1;
    y1 = prctile(data2plot, 25);
    y2 = prctile(data2plot, 75);
    rectangle('Position', [x1 y1 0.5 y2-y1], 'LineWidth', 1.5,...
        'EdgeColor', col2_null);
    ly1 = prctile(data2plot, 5);
    ly2 = prctile(data2plot, 95);
    p2 = line([bb bb], [ly1 y1], 'Color', col2_null, 'LineWidth', 1.5);
    line([bb bb], [y2 ly2], 'Color', col2_null, 'LineWidth', 1.5)
    ym = mean(data2plot);
    line([x1 x1+0.5], [ym ym], 'Color', col2_null, 'LineWidth', 1.5)

    % observed
    data2plot = Obc_all(:,bb);
    r1 = 0.9 + bb-1;
    r2 = 1.1 + bb-1;
    jit = r1 + (r2-r1).*rand(1,100);
    s1 = scatter(jit, data2plot, 25, 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', col1_obs);
    hold on
    x1 = 0.75 + bb-1;
    y1 = prctile(data2plot, 25);
    y2 = prctile(data2plot, 75);
    rectangle('Position', [x1 y1 0.5 y2-y1], 'LineWidth', 1.5,...
        'EdgeColor', col2_obs);
    ly1 = prctile(data2plot, 5);
    ly2 = prctile(data2plot, 95);
    p1 = line([bb bb], [ly1 y1], 'Color', col2_obs, 'LineWidth', 1.5);
    line([bb bb], [y2 ly2], 'Color', col2_obs, 'LineWidth', 1.5)
    ym = mean(data2plot);
    line([x1 x1+0.5], [ym ym], 'Color', col2_obs, 'LineWidth', 1.5)
        
end
xlim([2 13])
ylim([-0.02 0.45])
set(gca, 'XTick', 2:12, 'XTickLabel', {'2','3','4','5','6','7','8','9','10','11','12'},...
    'YTick', 0:0.05:0.4, 'YTickLabel', {'0','','0.1','','0.2','','0.3','','0.4'}, 'FontSize', 15)
box on; grid on;
xlabel('Number of modules')
ylabel('O-information (between modules)')
legend([p1 p2], {'observed', 'null'}, 'Location', 'northwest')

% ttest between obs and rand
for ii= num_comm(2:end)
    [h(ii), pval(ii)] = ttest2(Obc_all(:,ii), Obc_rand(:,ii));
    effect(ii) = meanEffectSize(Obc_all(:,ii), Obc_rand(:,ii), 'Effect','meandiff').Effect;
end
col1 = rgb('orange');
col2 = rgb('lightslategray');
ff = figure; hold on
line([2 13], [0.01 0.01], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.1, 'LineStyle', '--')
bb = bar(abs(effect));
bb.FaceColor = 'flat';
axis square; box on
xlim([2 13])
bb.CData(1:4,:) = repmat(col2, 4,1);
bb.CData(5:end,:) = repmat(col1, 8,1);
bb.FaceAlpha = 0.9;
set(gca, 'XTick', 3:12, 'YTick', 0.01:0.03:0.14, 'FontSize', 15)
xlabel('Number of modules')
ylabel('effect size (Oinfo_{obs} ; Oinfo_{rand})')

% correlation between oinfo and tcscore
for ii= num_comm
    [rho_obs(ii), pval_obs(ii)] = corr(cost_out_all(:,ii), Obc_all(:,ii)); 
    [rho_null(ii), pval_null(ii)] = corr(cost_out_all(:,ii), Obc_rand(:,ii)); 
end
ff= figure; % plot for 5 modules
scatter(cost_out_all(:,5), Obc_all(:,5), 60, 'MarkerFaceColor', rgb('SteelBlue'),...
        'MarkerEdgeColor', [0 0 0]);
axis square; box on; grid on
xlim([6.45 6.85])
ylim([0.006 0.026])
set(gca,'XTick', 6.5:0.1:6.8, 'YTick', 0.008:0.004:0.024, 'FontSize', 15)
xlabel('TC_{score} (within modules)')
ylabel('O-information (between modules)')
SW = [min(xlim) min(ylim)]+[diff(xlim) diff(ylim)]*0.05;
text(SW(1), 0.0235, sprintf('r = %g ; p_{val} < 0.005', round(rho_obs(5),2)),...
    'FontSize', 14, 'Color', [0.4 0.4 0.4], 'VerticalAlignment','bottom',...
    'HorizontalAlignment','left')
title('5 MODULES')
exportgraphics(ff, fullfile(fig_dir,'Obm_vs_TC_5mod.png'), 'Resolution', 300)
col1 = rgb('orange');
col2 = rgb('lightslategray');
sigrho = find(pval_obs(3:end)<0.05);
nonsigrho = find(pval_obs(3:end)>=0.05);
ff = figure; % plot rho
bb = bar(rho_obs(3:end));
bb.FaceColor = 'flat';
xlim([0 11]); ylim([-0.4 0.1])
axis square
bb.CData(nonsigrho,:) = repmat(col2, length(nonsigrho),1);
bb.CData(sigrho,:) = repmat(col1, length(sigrho),1);
bb.FaceAlpha = 0.9;
set(gca, 'XTick', 1:10, 'XTickLabel', {'3','4','5','6','7','8','9','10','11','12'},...
    'YTick', -0.4:0.1:0.1, 'FontSize', 15)
xlabel('Number of modules')
ylabel('corr. (TC_{score} ; O-information)')

% plot partition at 5 modules
freq5 = find_FrequentPartition(C_out_all(back2ord,:,5));
ciu5al = get_AlignedPartitions(freq5, ciu(back2ord,5));
col_comm = cbrewer('qual', 'Set3', 5, 'pchip');
for i=1:N
    col_nodes(i,:) = col_comm(freq5(i,1),:);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)


% try to find an optimal scale

tc_norm = mean(cost_out_all(:,3:end),1)/(max(mean(cost_out_all(:,3:end),1)));
o_norm = mean(Obc_all(:,3:end),1)/(max(mean(Obc_all(:,3:end),1)));
eff = abs(effect(3:end));
eff_norm = eff./(max(eff));
index = tc_norm.*((1-o_norm).*eff);

colors = cbrewer('qual', 'Paired', 10, 'pchip');

ff = figure;
plot(3:12, index, 'Color', colors(5,:))
hold on
scatter(3:12, index, 50, 'MarkerFaceColor', colors(6,:), 'MarkerEdgeColor', 'none')
xlim([2 13])
ylim([0 0.025])
axis square; grid on
set(gca, 'FontSize', 15, 'XTick', 3:12)
xlabel('Number of modules')
ylabel('B (segregation-integration balance)')

% plot partition at 9 modules
freq9 = find_FrequentPartition(C_out_all(back2ord,:,9));
ciu9al = get_AlignedPartitions(freq9, ciu(back2ord,9));
col_comm = cbrewer('qual', 'Set3', 9, 'pchip');
for i=1:N
    col_nodes(i,:) = col_comm(freq9(i,1),:);
end
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)


%% relative integration coefficient


N = 200; allN = 1:N;

nc = 7;

comm = squeeze(C_out_all(:,:,nc));
comm = comm(back2ord,:);

% work on centroid
cent = find_FrequentPartition(comm);

mod_size = histcounts(cent, [0.5 : 1 : nc+0.5]);

for m=1:nc
    ind = find(cent==m);
    TCmod(m) = calcI_alt(FC(ind,ind));
    TCmod_norm(m) = calcI_alt(FC(ind,ind))/length(ind);
end
TCall = calcI_alt(FC);
TCall_norm = calcI_alt(FC)/N;

for ii= 1:N
    lab = cent(ii);
    ind = find(cent==lab);
    newind = ind(not(ismember(ind,ii)));
    newallN = allN(not(ismember(allN,ii)));
    dec_within(ii) = TCmod(lab) - calcI_alt(FC(newind,newind));
    dec_all(ii) = TCall - calcI_alt(FC(newallN,newallN));
end

decwn = dec_within./(max(dec_within));
decan = dec_all./(max(dec_all));

ric = decwn./decan;

connector = find(ric<prctile(ric,10));
provincial = find(ric>1);

colors = cbrewer('div', 'PiYG', 15, 'pchip');
col_nodes = repmat(colors(8,:), N,1);
col_nodes(provincial,:) = repmat(colors(14,:), length(provincial), 1);
col_nodes(connector,:) = repmat(colors(2,:), length(connector), 1);
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_nodes)

for ii= 1:length(yeo_lab)
    pos = find(yeo7==ii);
    ricyeo(1,ii) = length(find(ismember(pos,connector)));
    ricyeo(2,ii) = length(find(ismember(pos,provincial)));
end
ff= figure;
by = bar(1:7, ricyeo, 'FaceColor', 'flat', 'Horizontal','on', 'BarWidth', 1,...
    'FaceAlpha',0.8);
by(1).CData = colors(2,:);
by(2).CData = colors(14,:);
axis square; box on
xlim([0 17])
set(gca, 'XTick', [0 5:5:15], 'YTick', 1:7, 'YTickLabel', yeo_lab,...
    'FontSize', 15, 'XGrid','on')
xlabel('nodes count')

ff= figure; hold on
rectangle('Position', [min(ric)-0.05 0 prctile(ric,10)-(min(ric)-0.05) 60],...
            'FaceColor', colors(6,:), 'EdgeColor', 'none')
rectangle('Position', [1 0 (max(ric)+0.05)-1 60],...
            'FaceColor', colors(10,:), 'EdgeColor', 'none')
histogram(ric, 'FaceColor', rgb('lightsteelblue'), 'FaceAlpha', 1)
line([1 1], [0 60], 'Color', rgb('darkslategray'), 'LineStyle', '--', 'LineWidth', 1.5)
xlim([min(ric)-0.051 max(ric)+0.051])
ylim([-0.1 60.1])
axis square; box on;
set(gca, 'FontSize', 15)
xlabel('Relative Integration Coefficient')
ylabel('Number of Nodes')


% heatmap considering all the iterations

ric_heat = zeros(N,iter);
connector_heat = zeros(N,1); provincial_heat = zeros(N,1);
for it= 1:iter
    
    tmp = squeeze(comm(:,it));

    mod_size = histcounts(tmp, [0.5 : 1 : nc+0.5]);
    
    for m=1:nc
        ind = find(tmp==m);
        TCmod(m) = calcI_alt(FC(ind,ind));
    end
    TCall = calcI_alt(FC);
    
    for ii= 1:N
        lab = tmp(ii);
        ind = find(tmp==lab);
        newind = ind(not(ismember(ind,ii)));
        newallN = allN(not(ismember(allN,ii)));
        dec_within(ii) = TCmod(lab) - calcI_alt(FC(newind,newind));
        dec_all(ii) = TCall - calcI_alt(FC(newallN,newallN));
    end
    
    decwn = dec_within./(max(dec_within));
    decan = dec_all./(max(dec_all));
    
    ric_heat(:,it) = decwn'./decan';
    
    ind_conn = find(ric_heat(:,it)<prctile(ric_heat(:,it),10));
    ind_prov = find(ric_heat(:,it)>1);
    
    connector_heat(ind_conn,1) = connector_heat(ind_conn,1)+1;
    provincial_heat(ind_prov,1) = provincial_heat(ind_prov)+1;

end

colors_heat = cbrewer('div', 'PiYG', 201, 'pchip');

col_prov_heat = get_proportionalColors(provincial_heat, colors_heat(100:201,:));
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_prov_heat)

col_conn_heat = get_proportionalColors(connector_heat, flipud(colors_heat(1:101,:)));
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_conn_heat)


% colorbar
ff= figure;
imagesc(1:100)
colormap(colors_heat(100:201,:))
set(gcf, 'Position', [360 198 240 420])
set(gca, 'DataAspectRatio', [1 0.09 1], 'FontSize', 17, 'XTick', [1 100],...
    'YTick', [], 'XTickLabel', [0 100])
exportgraphics(ff, fullfile(fig_dir,...
    sprintf('RIC_colorbarPROV_heatiter_%g.png', nc)), 'Resolution', 300)
ff= figure;
imagesc(1:100)
colormap(flipud(colors_heat(1:101,:)))
set(gcf, 'Position', [360 198 240 420])
set(gca, 'DataAspectRatio', [1 0.09 1], 'FontSize', 17, 'XTick', [1 100],...
    'YTick', [], 'XTickLabel', [0 100])


% heat across number of clusters

ric_heat2 = zeros(N,10);
connector_heat2 = zeros(N,1); provincial_heat2 = zeros(N,1);
for nn= 3:12
    clear TCmod
    comm_heat = squeeze(C_out_all(:,:,nn));
    comm_heat = comm_heat(back2ord,:);

    cent_heat = find_FrequentPartition(comm_heat);

    mod_size = histcounts(cent_heat, [0.5 : 1 : nn+0.5]);

    for m=1:nn
        ind = find(cent_heat==m);
        TCmod(m) = calcI_alt(FC(ind,ind));
    end
    TCall = calcI_alt(FC);

    for ii= 1:N
        lab = cent_heat(ii);
        ind = find(cent_heat==lab);
        newind = ind(not(ismember(ind,ii)));
        newallN = allN(not(ismember(allN,ii)));
        dec_within(ii) = TCmod(lab) - calcI_alt(FC(newind,newind));
        dec_all(ii) = TCall - calcI_alt(FC(newallN,newallN));
    end
    

    decwn = dec_within./(max(dec_within));
    decan = dec_all./(max(dec_all));

    ric_heat2(:,nn) = decwn./decan;
    
    ind_conn = find(ric_heat2(:,nn)<prctile(ric_heat2(:,nn),10));
    ind_prov = find(ric_heat2(:,nn)>1);
    
    connector_heat2(ind_conn,1) = connector_heat2(ind_conn,1)+1;
    provincial_heat2(ind_prov,1) = provincial_heat2(ind_prov)+1;
    
end

colors_heat = cbrewer('div', 'PiYG', 201, 'pchip');

col_prov_heat = get_proportionalColors(provincial_heat2, colors_heat(100:201,:));
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_prov_heat)

col_conn_heat = get_proportionalColors(connector_heat2, flipud(colors_heat(1:101,:)));
ff = figure;
parc_plot(surfStruct,annotMap,annotName,1:200,'cMap',col_conn_heat)

ff= figure;
imagesc(1:100)
colormap(colors_heat(100:201,:))
set(gcf, 'Position', [360 198 240 420])
set(gca, 'DataAspectRatio', [1 0.09 1], 'FontSize', 17, 'XTick', [1 100],...
    'XTickLabel', [0 10], 'YTick', [])

ff= figure;
imagesc(1:100)
colormap(flipud(colors_heat(1:101,:)))
set(gcf, 'Position', [360 198 240 420])
set(gca, 'DataAspectRatio', [1 0.09 1], 'FontSize', 17, 'XTick', [1 100],...
    'XTickLabel', [0 10], 'YTick', [])

colors = cbrewer('div', 'PiYG', 15, 'pchip');
yl = [120 100 60 70 60 50 60 90 80 70];
for nn= [3 4 5 6 8 10 11 12]
    figure('Visible', 'off'); histogram(ric_heat2(:,nn))
    xl = get(gca, 'XLim');

    ff= figure; hold on
    rectangle('Position', [xl(1) 0 prctile(ric_heat2(:,nn),10)-xl(1) yl(nn-2)],...
        'FaceColor', colors(6,:), 'EdgeColor', 'none')
    rectangle('Position', [1 0 xl(2)-1 yl(nn-2)],...
        'FaceColor', colors(10,:), 'EdgeColor', 'none')
    histogram(ric_heat2(:,nn), 'FaceColor', rgb('lightsteelblue'), 'FaceAlpha', 1)
    line([1 1], [0 60], 'Color', rgb('darkslategray'), 'LineStyle', '--', 'LineWidth', 1.5)
    xlim([xl(1) xl(2)])
    ylim([-0.1 yl(nn-2)+0.1])
    axis square; box on;
    set(gca, 'FontSize', 15)
    xlabel('Relative Integration Coefficient')
    ylabel('Number of Nodes')
    title(sprintf('%g MODULES', nn))

end