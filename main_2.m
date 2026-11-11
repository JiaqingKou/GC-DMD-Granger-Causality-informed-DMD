%% 0. Settings and input paths — Figure 5: three Reynolds numbers
% Only path initialization and the plotting entry point differ from the source project.
clear; clc; close all;
root=fileparts(mfilename('fullpath')); addpath(fullfile(root,'methods'),fullfile(root,'utils'),fullfile(root,'figures')); % Shared project paths
cfg=cavity_plot_settings();
cfg.output_dir=prepare_results_dir(root,2);
missing={};
for k=1:numel(cfg.data_files)
    if ~isfile(cfg.data_files{k}), missing{end+1}=cfg.data_files{k}; end %#ok<SAGROW>
end
if ~isempty(missing)
    error('CavityPlot:MissingData', ...
        'Missing input MAT files:\n%s\nCopy files into data/ or edit cavity_plot_settings.m.',strjoin(missing,newline));
end
%% 1. Load snapshots; 2. Mean subtraction and snapshot SVD/POD
% Figure 5 compares reduced representations across Reynolds numbers;
% this script does not fit DMD or GC-DMD evolution operators.
Results=cell(1,numel(cfg.reynolds));
for k=1:numel(Results)
    fprintf('\nProcessing Re=%d ...\n',cfg.reynolds(k));
    D=load_cavity_snapshots(cfg.data_files{k},cfg.struct_name,cfg.field_name,cfg.n_snapshots,0);
    Xraw=D.Xraw;
    mu=mean(Xraw,2);
    Xc=Xraw-mu;
    Results{k}=compute_cavity_pod(Xraw,Xc,mu,D.N,D.time(:),cfg.reynolds(k),cfg);
end
%% 3. Plot manuscript Figure 5 (two original output panels)
plot_fig5(Results,cfg);
fprintf('\nFigures and numerical results: %s\n',cfg.output_dir);
