function [IO,fOK] = SX_Startup (IO)
%% HELP %%
% This function perform really basic operation needed when SVEC is started:
%   * if saving, check if filename selected is not used yet, eventually generates a new file name
%   * if loading, creates a new filename for new file eventually saved (to avoid overwriting)
%   * obtain the time
%   * create log file
% 
% INPUT 
% IO  : Struct fs save parameters
% 
% OUTPUT
% IO  : Updated struct fs save parameters
% fOK : flag that control if error occours during simulation. (1 = no errors <> 0= errors, stop simulation)

    %% IO INPUT CHECK %%
    fOK = 1;

    f1  = IO.fMODE ~= 0 && IO.fMODE ~= 1;
    f2  = IO.fDPLT ~= 0 && IO.fDPLT ~= 1;
    f3  = IO.fDRPT ~= 0 && IO.fDRPT ~= 1;
    f4  = IO.fSRPT ~= 0 && IO.fSRPT ~= 1;
    f5  = IO.fSPLT ~= 0 && IO.fSPLT ~= 1;

    if f1
        warning('S1_InputCheck:IO','Unknown value of IO.fMODE. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end

    if f2
        warning('S1_InputCheck:IO','Unknown value of IO.fDPLT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f3
        warning('S1_InputCheck:IO','Unknown value of IO.fDRPT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f4
        warning('S1_InputCheck:IO','Unknown value of IO.fSPLT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    if f5
        warning('S1_InputCheck:IO','Unknown value of IO.fSRPT. IO variables can either be 0 or 1');
        SX_Logfile ('e',{lastwarn});
        fOK = 0;
    end
    clear f1 f2 f3 f4 f5

    
    %% CREATE SAVE FOLDER %%
    mkdir (fullfile(IO.save_dir));
    
    %% CHECK FILE NAME %%
    if IO.fMODE % saving
        outFile         = SX_FileCheck([IO.filename '.mat'],IO.save_dir,3,1); % check if .mat file with same name exist and eventually obtain new name
        [~,shortFile,~] = fileparts(outFile);                                 % obtain file name with no extension
        IO.filename     = shortFile;                                          % new file name is setted
        IO.fullsavename = fullfile(IO.save_dir, IO.filename);                 % fullfile saving name
        IO.fullloadname = IO.fullsavename;                                    % fullfile loading name
        SX_Logfile ('c', [IO.fullsavename,' logfile.log']);                   % creation of logfile

    elseif ~IO.fMODE % loading
        outFile         = SX_FileCheck([IO.filename '.mat'],IO.save_dir,3,1);    % obtain file name to eventually save new outputs
        [~,shortFile,~] = fileparts(outFile);                                    % obtain file name with no extension
        IO.newname      = shortFile;                                             % new saving name setted
        IO.fullloadname = fullfile(IO.save_dir, IO.filename);                    % fullfile load name
        IO.fullsavename = fullfile(IO.save_dir, IO.newname);                     % fullfile saving name
        SX_Logfile ('c', [IO.fullsavename,' logfile.log']);                      % creation of logfile
        warning('SX_Startup:LOAD','Loading simulation file: %s.mat',IO.filename)
        SX_Logfile ('v', {lastwarn});
    end

end