data_dir = '/nfs/data3/gergo/Mijung/Figure4/Output_data_highres';
data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat';
code_dir = '/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code';

Model =  'NSFR';
% Model = 'PLDS';

for k =3 
    for runnum =7 
        cd(code_dir);
        output_folder = [data_dir filesep 'k' num2str(k) '_run' num2str(runnum)];
        if ~exist(output_folder, 'dir'), mkdir(output_folder); end
        
        % Call runAlexsdata_prediction_func(data_file, output_folder, code_dir, k, Model)
        
        %Create a SLURM script file in the corresponding folder with the
        %correct parameters, then call srun via system
        script_file_path = create_slurm_script( data_file, output_folder, code_dir, k, Model );
        
        
        %Call via command line of linux
%         matlab_cmd = ['/opt/matlab-R2013a/bin/matlab -nodesktop -nosplash -singleCompThread -logfile ' output_folder filesep 'logfile.txt'];
        func_call = [' -r "cd Figure4;' 'runAlexsdata_prediction_func ' data_file ' ' output_folder ' ' code_dir ' ' num2str(k) ' ' Model ';" '];
        eval(func_call(6:end-2))
% %         
% %         system([matlab_cmd func_call]);

        
    end
end

% exit;

