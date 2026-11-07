% GET_SUBJECT_IDS List subject subdirectories in deterministic name order.
% Written by Junrui Zhang, 2025
%
% Input: root_dir, a directory containing only subject folders.
% Output: subject_ids, a column cell array of subject-folder names.
% Dependencies: MATLAB.

function subject_ids = get_subject_ids(root_dir)
assert(isfolder(root_dir), 'The subject-input directory does not exist.');
entries = dir(root_dir);
is_subject = [entries.isdir] & ~ismember({entries.name}, {'.', '..'});
names = {entries(is_subject).name};
subject_ids = sort(names(:));
assert(~isempty(subject_ids), 'The input directory contains no subject folders.');
end
