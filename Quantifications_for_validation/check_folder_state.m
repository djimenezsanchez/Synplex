function [] = check_folder_state(folderName)
    if ~exist(folderName, 'dir')
       mkdir(folderName)
    end
end