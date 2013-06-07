% uncomment a location below -- adapt to local installation

% location = 'home_matteo';
% location = 'home_luca';
location = 'work_matteo';

restoredefaultpath

if strmatch(location,'home_matteo','exact')    
    dir1='C:\dynare\4.3.1\matlab';
    dir2='C:\E\occbin_beta\toolkit_files_20130531';
    dir3='C:\E\occbin\occbin_20130531\toolkit_files';
    rmpath 'C:\Program Files (x86)\MATLAB\R2009b\toolbox\ident\idobsolete'
elseif strmatch(location,'home_luca','exact')
    dir1='/users/jason/documents/matlab/dynare/4.1.0/matlab';
    dir2='../toolkit_files';
elseif strmatch(location,'work_matteo','exact')  
    dir1='C:\E\dynare\4.3.1\matlab';
    dir2='C:\E\occbin\occbin_20130531\toolkit_files_private';
    dir3='C:\E\occbin\occbin_20130531\toolkit_files';
else 
    error('Specify path to Dynare installation')
end

path(dir1,path);
path(dir2,path);
path(dir3,path);

dynare_config

