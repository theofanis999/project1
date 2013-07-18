masterclock=clock;
hour=num2str(masterclock(4));
minu=num2str(masterclock(5));
diary(strcat('masterAAMC_',date,'_',hour,'_',minu,'.txt'));
diary on;
clear
clc
pwd
ls
AAMC(1,2)
clear
AAMC(2,2)
clear
AAMC(3,2)
clear
AAMC(1,100)
clear
AAMC(2,100)
clear
AAMC(3,100)
clear
diary off
