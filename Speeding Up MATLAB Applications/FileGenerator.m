%% Creating a large file
% To create a large file to be used for the textscan demo.
% Contains a array for different times for different
% models performances at each time step.

Nfiles=10;

for k=1:Nfiles

    %% Generating Model Names
    Nmodels= 100;   % Total Number of Models
    Ntimes= floor(10000*rand) ; % Total Number of Time Steps
    maxTime=1; % run for 1 hr
    Filename=strcat('ModelResults',num2str(k),'.txt');
    fid = fopen(Filename,'a+');

    %%
    for i=1:Nmodels
        modelname{i}=['Test  ' num2str(i)];
    end

    timecounter= ['TimeSteps  '  num2str(Ntimes)];
    
  %  modelnamechar=char(modelname{:});
  %  timecounterchar=char(timecounter);

    %% Generating Model Values

    timevector=linspace(0,maxTime,Ntimes);
    basevalues=timevector.^2 + (Ntimes^2)* sin(timevector/5).^2;
    modelResults=repmat(basevalues,Nmodels,1);
    modelResults=modelResults.*(rand(Nmodels,Ntimes).^2);

    %% Writing Model Values to a File
    
    fprintf(fid,timecounter);
    fprintf(fid,'\n');
    
     for i=1:Nmodels
         myModel = modelname{i};
         fprintf(fid,myModel);
         fprintf(fid,'\n');
         fprintf(fid,'%E %E \n',[timevector;modelResults(i,:)]);
     end

    %%
    fclose(fid);
end