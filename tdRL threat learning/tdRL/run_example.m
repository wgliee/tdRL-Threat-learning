function run_example()
%RUN_EXAMPLE Reproduce the complete workflow using Example data.xlsx.

    repositoryDir = fileparts(mfilename('fullpath'));
    run_analysis(fullfile(repositoryDir, 'Example data.xlsx'), ...
        fullfile(repositoryDir, 'results'));
end
