function safeStr = escapeUnderscores(str)
    % Helper function to escape underscores in strings for MATLAB titles/labels
    safeStr = strrep(str, '_', '\_');
end