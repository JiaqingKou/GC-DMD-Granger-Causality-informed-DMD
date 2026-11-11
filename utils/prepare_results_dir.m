function output_dir = prepare_results_dir(root, main_id)
%PREPARE_RESULTS_DIR Prepare a shared results folder without removing other mains' figures.
% main_id = 1: VAR(1), figs 1-4; 2: POD, fig 5; 3: cavity 13k, figs 6-7.
% Only the figures owned by the current main are deleted on rerun.
if nargin < 2
    error('prepare_results_dir:MissingMainID', 'Specify main_id (1, 2 or 3).');
end
output_dir = fullfile(root,'results');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
switch main_id
    case 1
        ids = 1:4;
    case 2
        ids = 5;
    case 3
        ids = 6:7;
    case 4
        ids = 8:12;
    case 7
        ids = [];
        for name = {'figA1a.eps','figA1b.eps','figA1c.eps','figA1d.eps','figA2.eps'}
            target=fullfile(output_dir,name{1});
            if exist(target,'file'), delete(target); end
        end
    case 6
        ids = 13:16;
    case 5
        ids = []; % Appendix B uses B1/B2 filenames rather than main-text figure numbers.
        for name = {'B1.eps','B2a.png','B2b.png','B2c.png','B2d.png'}
            target=fullfile(output_dir,name{1});
            if exist(target,'file'), delete(target); end
        end
    otherwise
        error('prepare_results_dir:UnknownMainID','Unknown main_id: %d',main_id);
end
files = dir(output_dir);
for k = 1:numel(files)
    if files(k).isdir, continue; end
    name = files(k).name;
    % Match only a complete figure number, optionally followed by a subpanel
    % letter, and only EPS/PNG outputs. Thus fig1 never matches fig19.
    tok = regexp(name, '^fig([0-9]+)([a-zA-Z]?)\.(eps|png)$', 'tokens', 'once', 'ignorecase');
    if isempty(tok), continue; end
    number = str2double(tok{1});
    if ismember(number,ids)
        delete(fullfile(output_dir,name));
    end
end
end
