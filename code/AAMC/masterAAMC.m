masterclock=clock;
hour=num2str(masterclock(4));
minu=num2str(masterclock(5));
diary(strcat('masterAAMC_',date,'_',hour,'_',minu,'.txt'));
diary on;
clear
clc
pwd
ls
AAMC(80,1,100)
clear
AAMC(80,2,100)
clear
AAMC(80,3,100)
%clear
%AAMC(250,1,200)
%clear
%AAMC(250,2,200)
%clear
%AAMC(250,3,200)
%clear
diary off
