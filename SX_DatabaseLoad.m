function [structLOAD,fOK] = SX_DatabaseLoad (type,selection,fOK)
%% HELP %%
% This function load the data of the selected item from the corresponding
% database % DATABASE TO BE COMPLETED !!!
% INPUT:
% type      : database type
% selection : name of object whose data are required
%
% OUTPUT
% structLOAD : struct of required data
% fOK        : error flag
%
% NOTE: - this function can be shortened by using 'eval' function
%       - 'eval' function can be pretty tricky,buggy and unstable. It is better avoid it
% * HOW TO DO IT WITH EVAL FUNCTION *
% function [structLOAD,fOK] = DatabaseLoad (type,selection,fOK)
% if fOK
%  load(['_',type,'DB.mat'])
%  tempDB = eval([upper('type'),'db']);
%  if isfield(tempDB,selection)
%   structLOAD=tempDB.(selection);
%  else
%   warning('CALL_SVECmodelMain:Selection','Unknown value of %s_selected',type);
%   structLOAD = [];   fOK = 0;
%   SX_Logfile ('w',{lastwarn});
%  end, clear tempDB
% end, end

    %% LOAD DATABASE %%
    if fOK
        switch type
            case 'Mach'
                load('DBmach.mat')
                if isfield(MACHdb,selection)
                    structLOAD = MACHdb.(selection);
                else
                    warning('DatabaseLoad:selection','Unknown value of Machine_selected');
                    structLOAD = [];   fOK = 0;
                    SX_Logfile ('w',{lastwarn});
                end, clear MATTEISTDdb

            case 'Gas'
                load ('DBgas.mat')
                if isfield(GASdb,selection)
                    structLOAD = GASdb.(selection);
                else
                    warning('DatabaseLoad:selection','Unknown value of Gas_selected');
                    SX_Logfile ('w',{lastwarn});
                    structLOAD = [];   fOK = 0;
                end, clear GASdb

            case 'Oil'
                load ('DBoil.mat');
                if isfield(OILdb,selection)
                    structLOAD = OILdb.(selection);
                else
                    warning('DatabaseLoad:selection','Unknown value of Liquid_selected');
                    SX_Logfile ('e',{lastwarn});
                    structLOAD = [];   fOK = 0;
                end, clear OILdb

            case 'Vane'
                load ('DBvane.mat')
                if isfield(VANEdb,selection)
                    structLOAD = VANEdb.(selection);
                else
                    warning('DatabaseLoad:selection','Unknown value of Vane_selected');
                    SX_Logfile ('e',{lastwarn});
                    structLOAD = [];   fOK = 0;
                end, clear VANEdb

            case 'Nozzle'
                load('DBnozzle.mat')
                if isfield(NOZZLEdb,selection)
                    structLOAD = NOZZLEdb.(selection);
                else
                    warning('DatabaseLoad:selection','Unknown type of nozzle selected in NOZZLES.name_nz');
                    SX_Logfile ('e',{lastwarn});
                    structLOAD = [];   fOK = 0;
                end, clear NOZZLESdb

            otherwise
                warning('DatabaseLoad:type','Database type unknown')
                SX_Logfile ('e',{lastwarn});
                structLOAD = [];   fOK = 0;
        end
    else
        structLOAD = [];
    end
end