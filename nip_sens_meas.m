function h = nip_sens_meas(y, labels, a_handle)
% h = nip_sens_meas(y, labels, a_handle)
% Shows the activity 'data' in the volume described by 'mesh'.
% Input:
%       y -> NcxNt. Contains the readings of Nc sensors over a time course
%       of Nt.
%       labels -> cell. Cell containing the labels of the channels
%       a_handle -> axes handle. This is where the plot is going to be
%               drawn.
% Output:
%       h -> plot handle.
%
% Additional comments :
% Example of labels: 
% eeg1010_labels={'Cz','Oz','PO3','PO4','O1','O2','P3','P4','P8', ...
%            'P7','C3','T7','C1','C2','C4','T8','Fz','F4','F8','F3','F7', ...
%            'AF3','AF4','Fp1','FC5','FC1','FC2','FC6','CP5','CP1','CP2', ...
%            'CP6','Fp2','Pz'};
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 26 Jan 2013


if nargin == 2
    if (iscell(labels))
        a_handle = gca;
    else
        a_handle = labels;
    end
end

const_offset = 1.5;
axes(a_handle);
offset = max(max(y));
t = linspace(0,1,size(y,2));
ytick_vector = zeros(size(y,1),1);

% Introduce an offset to each channel so they don't get overlapped
for i = 1:size(y,1)
    ytick_vector(i) = const_offset*offset*i;
    plot(t,y(i,:) + ytick_vector(i));
    hold on
end
xlabel('Time (secs)');
set(a_handle,'Ytick',ytick_vector);
if iscell(labels)  
    set(gca,'YTickLabel',labels);
end
ylim([ytick_vector(1)-const_offset*offset, ytick_vector(end) + const_offset*offset]);
h = gca;