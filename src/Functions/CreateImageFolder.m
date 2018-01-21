function [ folder ] = CreateImageFolder( folder_name )
    images_path = '../gif';
    mkdir(images_path, folder_name);
    folder = strcat(images_path,'/',folder_name);
end