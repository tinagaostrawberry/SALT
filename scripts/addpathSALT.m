function root = addpathSALT()
%ADDPATHSALT  Put the SALT tree (and MATPOWER) on the MATLAB path.
%
%   Only the source folders are added. After adding them every SALT
%   function is checked to make sure it actually resolves to the SALT copy:
%   several file names are shared with the legacy solver folders, and
%   MATLAB gives the current folder priority over the path, so a stale
%   working directory would otherwise silently substitute the old models.
%
%   Every source folder has to be listed in srcDirs below. A folder that
%   is missing is silently skipped, and the first thing that needs it
%   then fails with an undefined-function error far from the cause.

root = fileparts(fileparts(mfilename('fullpath')));

srcDirs = ["solver","scripts","scripts_contingency","utils","validation", ...
           fullfile("solver","indices"),   fullfile("solver","ss_models"), ...
           fullfile("solver","td_models"), fullfile("solver","init"), ...
           fullfile("solver","scripts")];
for f = srcDirs
    d = fullfile(root,f);
    if isfolder(d), addpath(d); end
end

% MATPOWER supplies psse2mpc / runpf / loadcase / define_constants
mp = 'C:\Users\tinagao\Desktop\Work\Software\matpower8.0';
if isfolder(mp) && isempty(which('runpf'))
    addpath(genpath(mp));
end

assertNoShadowing(root,srcDirs);
end

function assertNoShadowing(root,srcDirs)
shadowed = strings(0,1);
shadowedBy = strings(0,1);
for f = srcDirs
    d = fullfile(root,f);
    if ~isfolder(d), continue, end
    files = dir(fullfile(d,'*.m'));
    for k = 1:numel(files)
        [~,name] = fileparts(files(k).name);
        resolved = which(name);
        if isempty(resolved), continue, end
        if ~strcmpi(fileparts(resolved),d)
            shadowed(end+1,1)   = string(name);      %#ok<AGROW>
            shadowedBy(end+1,1) = string(fileparts(resolved)); %#ok<AGROW>
        end
    end
end
if isempty(shadowed), return, end
error('SALT:shadowedFunctions', ...
    ['%d SALT function(s) are shadowed by another folder, so the wrong ' ...
     'model code would run:\n  %s\nShadowing folder(s):\n  %s\n' ...
     'Change the current folder away from the legacy solver directory, ' ...
     'or remove it from the path.'], ...
    numel(shadowed), strjoin(shadowed.',', '), ...
    strjoin(unique(shadowedBy).',sprintf('\n  ')));
end
