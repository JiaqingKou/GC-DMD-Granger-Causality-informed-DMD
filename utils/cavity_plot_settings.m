function cfg=cavity_plot_settings()
% Central settings. Only input files must be provided by the user.
root=fileparts(fileparts(mfilename('fullpath'))); % Project root; utils/ is one level below
cfg.root=root;
cfg.data_files={fullfile(root,'data','Cavity13k.mat'), ...
    fullfile(root,'data','Cavity19k.mat'),fullfile(root,'data','Cavity30k.mat')};
cfg.reynolds=[13000 19000 30000];
cfg.struct_name='VelocityField'; cfg.field_name='Psi';
cfg.n_snapshots=800; cfg.snapshot_indices=[1 400 800];
cfg.pod_energy_target=0.99999; % 99.999%, NOT a prescribed integer rank
cfg.normalization='maxabs'; % q1/max(abs(q1)), independently for each case
cfg.field_source='truth'; % 'truth' (raw snapshots) or 'pod' (rank-r reconstruction)
cfg.interpolation_points=401; % original toolbox: -1:0.005:1
cfg.contour_levels=100; % original toolbox contour convention
cfg.field_color_limit=[]; % optional scalar or [limit13k limit19k limit30k]
% t labels indicate snapshot/time-step indices, not VelocityField.Time.
%  % []: percentile per case across its three fields
cfg.field_color_percentile=95; % stronger interior contrast, display only
cfg.field_colormap='smooth_redblue'; % smooth_redblue or original
cfg.show_colorbars=false;
cfg.coefficient_line_styles={'-','--','-.'};
cfg.slice_color=[0.93 0.69 0.08]; % golden yellow, readable on white
cfg.slice_line_width=1.8;
cfg.font_name='Times New Roman'; cfg.font_size=28; cfg.label_size=30;
cfg.legend_size=24; cfg.annotation_size=22; cfg.panel_label_size=22;
cfg.nondimensional_dt=0.1;
cfg.line_width=2.0; cfg.grid=false;
cfg.coefficient_position=[80 60 1500 800];
cfg.field_position=[100 40 1200 1455];
cfg.show_panel_labels=true; % case/time labels, no long figure titles
cfg.export_pdf=false; cfg.export_png=true; cfg.export_fig=false; % lightweight defaults
cfg.export_individual_fields=false; % optional nine separate field files
cfg.png_resolution=300; % set to 600 for final high-resolution raster output
cfg.output_dir=fullfile(root,'results');
end
