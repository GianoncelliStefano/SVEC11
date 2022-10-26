function outFile = SX_FileCheck(inFile,inPath,numDigits,startCount)
% Function to check if a given file already exists and create a new
% filename that can be used to avoid overwriting. New file names are
% created by appending a sequence like '_001' to the file name. The number
% is hereby increased until an unused file name is found. 
%
% Usage: outFile = SX_FileCheck(inFile,inPath,numDigits)
%
% INPUT
% inFile [string]  : The name of the file that should be checked. Make sure to add filetype if required.         
% inPath [string]  : The path of the file. Can be omitted to only check in the Matlab path.
% numDigits [-]    : The number of digits used for enumration. Default is 2 digits, so filenames are set up as e.g. 'inFile_01.mat' and so on.
% startCount [-]   : Defines the starting number for enumeration. Standard is 0 so the first file is created as 'inFile_00.mat'. 
%               
% OUTPUT
% outFile: The filename that can be used to save a file without overwriting existing data.
%       
% EXAMPLE: 
% fName = 'Example.mat'; %example file
% outFile1 = SX_FileCheck(fName,pwd,3,1);
% save(outFile1); %save empty file
% disp(['First file name: ' outFile1]);
% outFile2 = SX_FileCheck(fName,pwd,3,1); %second file
% save(outFile2);
% disp(['Second file name: ' outFile2]);
% outFile3 = SX_FileCheck(fName,pwd,3,1); %third file
% save(outFile3);
% disp(['Third file name: ' outFile3]);
% delete(outFile1,outFile2,outFile3) %delete empty files again

    %% CHECK INPUT %%
    if ~exist('inPath','var')
        inPath = [];     %number of digits when adding numbers to a filename
    end

    if ~exist('numDigits','var')
        numDigits = 2;     %number of digits when adding numbers to a filename
    end
    numchar = num2str(numDigits);

    if ~exist('startCount','var')
        startCount = 0;     %first value for counter. This determines the first filename during enumeration.
    end

    if ~strcmp(inPath(end),filesep) %make sure path has a seperator at the end
        inPath = [inPath filesep];
    end

    %% CHECK IF FILE EXIST AND ENUMERATE %%
    if exist([inPath inFile],'file') == 2           % file exists already, check for alternative
        [~,shortFile,fileType] = fileparts(inFile); % exclude the file type
        checker = true;                             % check for alternate file names
        lgt      = length(shortFile);               % length of shortFile
        lastchar = shortFile ((lgt-numDigits+1):lgt);
        elem4    = shortFile(lgt-numDigits);
        
        if strcmp(elem4,'_') && isnumeric(str2num(lastchar))  % file already enumerated
            realname = shortFile (1:lgt-numDigits-1);         % gain clean name os filename
            Cnt      = str2num(lastchar)+1;                   % counter for file name
            
        else                                                  % file with no enumeration
            realname = shortFile;                             % gain clean name os filename
            Cnt = startCount;                                 % counter for file name
        end
        
        % first attempt filename
        
        
        %increase counter until a non-existing file name is found
        while checker
            testPath = [inPath realname '_' num2str(Cnt, ['%0' numchar 'i']) fileType];
            if exist(testPath,'file') == 2   
                Cnt = Cnt + 1;            
            else
                checker = false;
            end
        end
        
        outFile = [realname '_' num2str(Cnt, ['%0' numchar 'i']) fileType];

    else
        outFile = inFile;
    end
end