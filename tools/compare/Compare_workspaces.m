% This script is useful for comparing workspace files generated in the
% model mapping process, or when analysing the associated rule of thumb. It
% relies on the function comp_struct.m (written by Michael Arant and available at
% https://uk.mathworks.com/matlabcentral/fileexchange/22752-compare-structures?s_tid=mwa_osa_a)
% 
% Update: 10-10-2019

close all; % close all figures
clearvars;
Model_Mapping_workspace = 1; % If true, the file is a workspace created from Model_Mapping_code
% Otherwise it's a workspace created to study the rule of thumb
subfolder = 'GB'; % GB, SL, SA, assortativity
filename = 'R015_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050'; % Copy here the file name
varname = 'theta_A'; % Name of variable I want to compare
varlist = {'Rg','Rg_A','Rg_H','Rh_H','zAH','zAHsim','zA','zAsim','zH','zHsim',...
    'piAHsim','piAsim','piHsim','tAHsim','tAsim','tHsim',...
    'hfs_AH','hfs_H','hfs_by_ic','inVSout','inVSout_A',...
    'r_AH_best','r_AH_high','r_A','r_H_best','r_H_high',...
    'v_AH','vH_AH',...
    'NGM_G_stored','NGM_A_stored','afs_each_typeAH','NGM_At_stored'};

% Path stuff
current_dir = cd;
if ispc
    if Model_Mapping_workspace
        folder = '-workspaces\';
    else
        folder = '-rule-of-thumb\';
    end
    fo = [current_dir,'\output',folder,subfolder,'\',filename,'.mat'];
    fs = [current_dir,'\saved',folder,subfolder,'\',filename,'.mat'];
    ft = [current_dir,'\tools\'];
else
    base = [current_dir,'/'];
    if Model_Mapping_workspace
        folder = '-workspaces/';
    else
        folder = '-rule-of-thumb/';
    end
    fo = [current_dir,'/output',folder,subfolder,'/',filename,'.mat'];
    fs = [current_dir,'/saved',folder,subfolder,'/',filename,'.mat'];
    ft = [current_dir,'/tools/'];
end
cd(ft);
[output_workspace, saved_workspace, differences] = comp_struct(load(fo),load(fs))
cd(current_dir);
% matlab.io.saveVariablesToScript('output.m')
% matlab.io.saveVariablesToScript('saved.m')
% open('output.m')
% open('saved.m')

l1 = load(fo);
l2 = load(fs);
% Summary
lv = length(varlist);
for iv = 1:lv
    var = varlist{iv};
    var1 = eval(['l1.',var]);
    var2 = eval(['l2.',var]);
    diff = var1(10:end,15:16)-var2(10:end,15:16);
    disp(['Max discrepancy in ',var,':'])
    disp(max(abs(diff(:))));
    disp('');
end
% Specific variable
var1 = eval(['l1.',varname]);
var2 = eval(['l2.',varname]);
disp(var1(10:end,15:16)-var2(10:end,15:16))

