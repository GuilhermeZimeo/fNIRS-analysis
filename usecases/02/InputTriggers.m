function InputTriggers(rest_dur, task_dur, pathname, overwrt)
%Function to define trigger markers in the *.nirs file
%
%this will only work for systems with 10Hz sampl. freq.
%for now, only for block design with a single condition
%i.e. alternating rest and task blocks (in this order)
%
%Last Update: 2018-03-03
%
%By: Guilherme A. Zimeo Morais
%Contact: gazmorais@gmail.com

if nargin < 3
    pathname = cd; %default pathname = current directory
end

if nargin < 4
    overwrt = 0; %do not overwrite original file
end

% [filename, pathname] = uigetfile('*.nirs', 'Pick *.nirs file');
% load([pathname filename],'-mat'); %load file in Homer format

%this follows code of loadDirectory to retrieve all files at once
files = rdir(fullfile(pathname,'**','*.nirs'));

%iterate over all *.nirs files found
for i=1:length(files)
    
    try 
        load(files(i).name,'-mat'); %load file in Homer format
    catch
        disp(['Error loading ' files(i).name]);
    end
    
    dpt = t; %get time (s) of each data point

    %find all time samples multiple of (rest+task)_dur
    rest_tr = find( ~mod(dpt, (rest_dur+task_dur)) );
    
    %find all time samples multiple of task_dur
    all_tr = find( ~mod(dpt, task_dur) );
    
    %get the time samples corresponding to task onset
    rest_idx = ismember(all_tr, rest_tr);
    task_tr = all_tr(~rest_idx);

    %first onset found should be 'task'    
    if(rest_tr(1) < task_tr(1))
        rest_tr(1) = [];
    end

    %last onset found should be 'rest'
    if(rest_tr(end) < task_tr(end))
        task_tr(end) = [];
    end

    s(task_tr) = 2; %define task trigger markers ("2")
    s(rest_tr) = 1; %trigger to compute task duration

    s = im2double(s); %convert int to double precision

    if overwrt
        filepath = files(i).name; %overwriting original file
    else
        filepath = [files(i).folder filesep 'new_file.nirs'];
    end
    
    %save current file with trigger markers
    save(filepath,'aux','brainsight','d','ml','s','SD','t');

end