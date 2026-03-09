function jobFile_synthesis(data_dim, type_spl, spl_file, input_file, ... 
    output_file, gS_input_path, gS_output_path, gS_job_path, gS_synth_path, type_discr)


% Function to automatically create a geoSplinter_synthesis job file from 
% the given input variables.
% 
% INPUT VARIABLES:
% Dimension of dataset (1D/2D): data_dim      
% Type of splines (linear/cubic): type_spl  
% Type of delta discretization (1d/2d): type_discr
% Spline file: spl_file 
% Input file: input_file                
% Output file: output_file
% Path to geoSplinter input folder: gS_input_path
% Path to geoSplinter output folder: gS_output_path
% Path to geoSplinter job folder: gS_job_path
% Path to geoSplinter synthesis folder: gS_synth_path
%
%
% Roberto Monti
% Politecnico di Milano
% Last update: February 2025


% Check on the number of input variables
if nargin == 9

    if type_spl ~= 2
        error('The number of input variables is not compatible with the chosen spline type')
    end

elseif nargin == 10

    if type_spl ~= 1
        error('The number of input variables is not compatible with the chosen spline type')
    end

end

% Create the job file
job_file = fullfile('.', gS_job_path, [output_file, '.job']);
fid = fopen(job_file, 'w');

if type_spl == 1   % linear splines

    fprintf(fid, '\n%s\n', num2str(data_dim));
    fprintf(fid, '%s\n', num2str(type_spl));
    fprintf(fid, '%s\n', num2str(type_discr));
    fprintf(fid, '%s\n', strcat(gS_output_path, filesep, spl_file));
    fprintf(fid, '%s\n', 'y');
    fprintf(fid, '%s\n', strcat(gS_input_path, filesep, input_file));
    fprintf(fid, '%s\n', strcat(gS_synth_path, filesep, output_file, '.txt'));

elseif type_spl == 2   % cubic splines
 
    fprintf(fid, '\n%s\n', num2str(data_dim));
    fprintf(fid, '%s\n', num2str(type_spl));
    fprintf(fid, '%s\n', strcat(gS_output_path, filesep, spl_file));
    fprintf(fid, '%s\n', 'y');
    fprintf(fid, '%s\n', strcat(gS_input_path, filesep, input_file));
    fprintf(fid, '%s\n', strcat(gS_synth_path, filesep, output_file, '.txt'));

end

fclose(fid);
