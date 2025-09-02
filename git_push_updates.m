function git_push_updates(commit_message)
% GIT_PUSH_UPDATES - Quick script to commit and push research updates
% 
% Usage: git_push_updates('Your commit message here')
%        git_push_updates()  % Will prompt for message
%
% Author: d.s.fraser@bham.ac.uk

if nargin < 1
    commit_message = input('Commit message: ', 's');
    if isempty(commit_message)
        commit_message = 'Update research files';
    end
end

fprintf('=== Pushing Research Updates ===\n');

% Check repository status
fprintf('Checking repository status...\n');
system('git status');

% Add all changes
fprintf('\nAdding changes...\n');
system('git add .');

% Commit changes
fprintf('\nCommitting changes...\n');
[status, result] = system(sprintf('git commit -m "%s"', commit_message));
if status ~= 0
    if contains(result, 'nothing to commit')
        fprintf('No changes to commit.\n');
        return;
    else
        error('Commit failed: %s', result);
    end
end

% Push to GitHub
fprintf('\nPushing to GitHub...\n');
[status, result] = system('git push origin main');
if status == 0
    fprintf('✓ Successfully pushed to GitHub!\n');
else
    fprintf('Push failed: %s\n', result);
end

end
