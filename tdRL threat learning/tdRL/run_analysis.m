function run_analysis(inputFile, outputDir)
%RUN_ANALYSIS Fit both models and perform the animal-level comparison.
%   RUN_ANALYSIS(INPUTFILE, OUTPUTDIR) accepts an .xlsx, .xls, .csv, or .txt
%   file in the format described in DATA_FORMAT.md. When arguments are
%   omitted, Example data.xlsx and the local results folder are used.

    repositoryDir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(inputFile)
        inputFile = fullfile(repositoryDir, 'Example data.xlsx');
    end
    if nargin < 2 || isempty(outputDir)
        outputDir = fullfile(repositoryDir, 'results');
    end

    inputFile = char(inputFile);
    outputDir = char(outputDir);
    addpath(repositoryDir, fullfile(repositoryDir, 'function'));
    check_dependencies();

    fprintf('\nInput:  %s\n', inputFile);
    fprintf('Output: %s\n', outputDir);
    fprintf('\n[1/3] Fitting the R-W model...\n');
    run(fullfile(repositoryDir, 'RW_G.m'));
    fprintf('\n[2/3] Fitting the tdRL model...\n');
    run(fullfile(repositoryDir, 'tdRL_G.m'));
    fprintf('\n[3/3] Comparing models...\n');
    run(fullfile(repositoryDir, 'compare_RW_vs_tdRL_G.m'));
    fprintf('\nAnalysis complete. Results are in: %s\n', outputDir);
end
