function jobFile_analysis(data_dim, type_spl, input_file, output_file, n_obs, ...
    n_spl, t_0, t_end, lambda, n_sig, gS_input_path, gS_output_path, gS_job_path, type_discr)

% Function to automatically create a geoSplinter_analysis job file from the
% given input variables.
% 
% INPUT VARIABLES:
% Dimension of dataset (1D/2D): data_dim      
% Type of splines (linear/cubic): type_spl  
% Type of delta discretization (1d/2d): type_discr
% Input file: input_file                
% Output file: output_file               
% Number of observations: n_obs            
% Number of nodes (= splines): n_spl        
% First abscissa (= time): t_0   
% Last abscissa (= time): t_end         
% Regularization parameter (λ): lambda        
% Number of significant digits: n_sig
% Path to geoSplinter input folder: gS_input_path
% Path to geoSplinter output folder: gS_output_path
% Path to geoSplinter job folder: gS_job_path
%
%
% Roberto Monti
% Politecnico di Milano
% Last update: May 2024


% Check on the number of input variables
if nargin == 13

    if type_spl ~= 2
        error('The number of input variables is not compatible with the chosen spline type')
    end

elseif nargin == 14

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
    fprintf(fid, '%s\n', strcat('.', filesep, gS_input_path, filesep, input_file));
    fprintf(fid, '%s\n', strcat('.', filesep, gS_output_path, filesep, output_file));
    fprintf(fid, '%s\n', num2str(n_obs));
    fprintf(fid, '%s\n', num2str(n_spl));
    fprintf(fid, '%s\n', num2str(t_0));
    fprintf(fid, '%s\n', num2str(t_end));
    fprintf(fid, '%s\n', num2str(lambda));
    fprintf(fid, '%s\n', num2str(n_sig));
    fprintf(fid, '%s\n', '.');

elseif type_spl == 2   % cubic splines

    % Manipulation of splines to add one at both ends
    % take the same time gap and add one spline before first epoch and one
    % after the last one.
    delta_t = (t_end - t_0) / (n_spl - 1);
    n_spl = n_spl + 2;
    t_0 = t_0 - delta_t;
    t_end = t_end + delta_t;    

    fprintf(fid, '\n%s\n', num2str(data_dim));
    fprintf(fid, '%s\n', num2str(type_spl));
    fprintf(fid, '%s\n', strcat('.', filesep, gS_input_path, filesep, input_file));
    fprintf(fid, '%s\n', strcat('.', filesep, gS_output_path, filesep, output_file));
    fprintf(fid, '%s\n', num2str(n_obs));
    fprintf(fid, '%s\n', num2str(n_spl));
    fprintf(fid, '%s\n', num2str(t_0));
    fprintf(fid, '%s\n', num2str(t_end));
    fprintf(fid, '%s\n', num2str(lambda));
    fprintf(fid, '%s\n', num2str(n_sig));
    fprintf(fid, '%s\n', '.');

end