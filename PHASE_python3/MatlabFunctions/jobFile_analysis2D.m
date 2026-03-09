function jobFile_analysis2D(data_dim, type_spl, input_file, output_file, n_obs, ...
    n_row, n_col, x_min, x_max, y_min, y_max, lambda, n_sig, gS_input_path, gS_output_path, gS_job_path, type_discr)

% Function to automatically create a geoSplinter_analysis job file from the
% given input variables.
% 
% INPUT VARIABLES:
% Dimension of dataset (1D/2D): data_dim      
% Type of splines (bilinear/bicubic): type_spl  
% Type of delta discretization (1d/2d): type_discr
% Input file: input_file                
% Output file: output_file               
% Number of observations: n_obs            
% Number of rows: n_row
% Number of columns: n_col 
% First X coordinate: x_min   
% Last X coordinate: x_max
% First Y coordinate: y_min   
% Last Y coordinate: y_max
% Regularization parameter (λ): lambda        
% Number of significant digits: n_sig
% Path to geoSplinter input folder: gS_input_path
% Path to geoSplinter output folder: gS_output_path
% Path to geoSplinter job folder: gS_job_path
%
%
% Roberto Monti
% Politecnico di Milano
% Last update: March 2025


% Check on the number of input variables
if nargin == 16

    if type_spl ~= 2
        error('The number of input variables is not compatible with the chosen spline type')
    end

elseif nargin == 17

    if type_spl ~= 1
        error('The number of input variables is not compatible with the chosen spline type')
    end

end

% Create the job file
job_file = fullfile('.', gS_job_path, [output_file, '.job']);
fid = fopen(job_file, 'w');

if data_dim == 2   % 2D splines

    if type_spl == 1   % linear splines
    
        fprintf(fid, '\n%s\n', num2str(data_dim));
        fprintf(fid, '%s\n', num2str(type_spl));
        fprintf(fid, '%s\n', num2str(type_discr));
        fprintf(fid, '%s\n', strcat('.', filesep, gS_input_path, filesep, input_file));
        fprintf(fid, '%s\n', strcat('.', filesep, gS_output_path, filesep, output_file));
        fprintf(fid, '%s\n', num2str(n_obs));
        fprintf(fid, '%s\n', num2str(n_row));
        fprintf(fid, '%s\n', num2str(n_col));
        fprintf(fid, '%s\n', num2str(x_min));
        fprintf(fid, '%s\n', num2str(x_max));
        fprintf(fid, '%s\n', num2str(y_min));
        fprintf(fid, '%s\n', num2str(y_max));
        fprintf(fid, '%s\n', num2str(lambda));
        fprintf(fid, '%s\n', num2str(n_sig));
        fprintf(fid, '%s\n', '.');
    
    elseif type_spl == 2   % cubic splines
    
        % Manipulation of splines to add one at both ends
        % take the same time gap and add one spline before first epoch and one
        % after the last one.
        delta_x = (x_max - x_min) / (n_col - 1);
        n_col = n_col + 2;
        x_min = x_min - delta_x;
        x_max = x_max + delta_x; 

        delta_y = (y_max - y_min) / (n_row - 1);
        n_row = n_row + 2;
        y_min = y_min - delta_y;
        y_max = y_max + delta_y; 
    
        fprintf(fid, '\n%s\n', num2str(data_dim));
        fprintf(fid, '%s\n', num2str(type_spl));
        fprintf(fid, '%s\n', strcat('.', filesep, gS_input_path, filesep, input_file));
        fprintf(fid, '%s\n', strcat('.', filesep, gS_output_path, filesep, output_file));
        fprintf(fid, '%s\n', num2str(n_obs));
        fprintf(fid, '%s\n', num2str(n_row));
        fprintf(fid, '%s\n', num2str(n_col));
        fprintf(fid, '%s\n', num2str(x_min));
        fprintf(fid, '%s\n', num2str(x_max));
        fprintf(fid, '%s\n', num2str(y_min));
        fprintf(fid, '%s\n', num2str(y_max));
        fprintf(fid, '%s\n', num2str(lambda));
        fprintf(fid, '%s\n', num2str(n_sig));
        fprintf(fid, '%s\n', '.');
    
    end

elseif data_dim == 1

    error('For 1D splines use jobFile_analysis.m')  

end

fclose(fid);