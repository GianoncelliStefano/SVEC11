function SX_Logfile(identifier,varagin)
% This function create and update a log file
% INPUT:
% identifier  : select what the function needs to do
%                 'c' - create or change log file
%                 'v' - log event
%                 'd' - log debug
%                 'w' - log warning
%                 'e' - log error
%                 'a' - log abort
% varagin     : cell array containing the input data {'savedirectory','log file name'}
% 
% EXAMPLE: create a new log file:
% SX_Logfile ('c')                                * create a file named 'mylog.txt in current directory'
% SX_Logfile ('c',{'mylogfile.txt'})              * create log file named 'mylogfile' in current directory
% SX_Logfile ('c',{'mydirectory\mylogfile.txt'})  * create logfile named mylogfile.txt in userdefined directory
%
% EXAMPLE: add a log
% SX_Logfile ('v',{'Simulation started'})                                    * add an event
% SX_Logfile ('d',{'Values not corresponding. Error computed: %2.2f', err})  * add a debug log
% logfile ('w',{'Inlet gas temperature seems to high'})                   * add a warning log
%
% NOTE   : is possible to use the command 'lastwarn' or 'lasterr' to get the last warning/error 
%           complete message as a string (and eventually its message identifier), to easy perform
%           a display warning and log it to the file via this function
% EXAMPLE:
% warning('Noozle %i is missing a field',3)
% SX_Logfile('w',{lastwarn})
% ans = WARNING  Nozzle 3 is missing a field

    %% DEFINITIONS %%
    persistent logname;
    switch identifier
        case 'c'
            if nargin == 1
                outFile = SX_FileCheck('mylog.log',pwd,3,1);  % check if a file with the same name exist, and eventually generates a new file name
                logname = outFile;                                    % fullfile of logfile
            elseif nargin == 2
                logname = varagin;
            end  
            
            % creates header of logfile
            logID = fopen (logname,'at');
            if isempty (fgetl(logID))
                fprintf (logID,['===========================================================','\n']);
                fprintf (logID,['\t','TIME','\t\t','TYPE','\t\t','MESSAGE','\n']);
                fprintf (logID,['===========================================================','\n\n']);
                fclose(logID);
                SX_Logfile('v',{'LogFile created'})
            else
                fclose(logID);
            end
            return
            
        case 'v'
            type  = 'EVENT';
        case 'd'
            type  = 'DEBUG';
        case 'w'
            type  = 'WARNING';
        case 'e'
            type  = 'ERROR';
        case 'a'
            type  = 'ABORT';
    end
    
    %% WRITE ON LOG FILE %%
    date   = datestr(now,'dd/mm/yyyy HH:MM:SS');    % get current date
    Finfo  = dbstack;                               % get current identated functions names
    numfz  = size(Finfo,1);                         % get number of identated functions
    msg    = varagin{1};                            % get message to display
    numval = size(varagin,2)-1;                     % get number of variables to display\print
    values = NaN(1,numval);                         % preallocation
    for i=1:numval                                  % set to values the variables to be printed
        values(i) = varagin {numval+1};
    end
    logID = fopen (logname,'at');                                % open logfile
    fprintf (logID,[date,'\t',type,'\t ',msg,'\n'],values);      % print first line (date + log type + message)
    for i=2:numfz                                                % print identated function names
        fprintf (logID,['\t\t\t\t In ',Finfo(i).name,' (line ',num2str(Finfo(i).line),')\n']);
    end
    fprintf(logID,'\n');
    fclose(logID);                                               % close logfile
    
end
