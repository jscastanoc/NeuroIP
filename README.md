NeuroIP
=======

Matlab Toolbox to solve the EEG inverse problem.
It has some basic functions that can be useful to those that are new to the topic.
This is more like a "personal" repo, but if you find something useful, feel free to used it for whatever you want.


Required Toolboxes (put them in the 'external/' dir):
toolbox_graph -> http://www.mathworks.com/matlabcentral/fileexchange/5355

nway(support for PARAFAC)  -> http://www.mathworks.com/matlabcentral/fileexchange/1088-the-n-way-toolbox

ltfat (support for Time-Frequency transformations (in some test_scripts) -> http://ltfat.sourceforge.net/

spm			-> http://www.fil.ion.ucl.ac.uk/spm/

fieldtrip 	-> http://fieldtrip.fcdonders.nl/

DAL (for S-FLEX solver) -> http://www.ibis.t.u-tokyo.ac.jp/ryotat/dal/

#-------- LEAD FIELD DATA ------------#
The lead field Matrix was taken from: 
http://eeg.pl/Members/jarekz/lead-field-data-for-subjects-form-the-paper/view

And was calculated using the Fieltrip toolbox:
http://fieldtrip.fcdonders.nl/


#------- SAMPLE SCRIPTS -------------#
script_LORETA : Example to reconstruct using MinimumNorm-like approaches
