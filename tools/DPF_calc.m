function DPF = DPF_calc(din, age_all)
%This calculates the DPF based on subject's age and wavelengths used
%based on Scholkmann and Wolf 2013 - DOI: 10.1117/1.JBO.18.10.105004
%
%Last Updated: 2018-03-04
%
%By: Guilherme A. Zimeo Morais
%Contact: gazmorais@gmail.com

if nargin < 2 %age not provided as input
    
    %retrieve demographics from data provided
    demographics = nirs.createDemographicsTable(din);
    
    %attempt to retrieve 'age' category from demographics table
    idxctg = strcmpi(demographics.Properties.VariableNames,'age');
    
    if any(idxctg) %category found
        age_all = table2array(demographics(:,idxctg));
        
    else %if no age info, assume age = 25
        age_all = 25*ones(1,length(din)); 
        disp('Age info not available. Approximating as 25 years old.');
    end
    
end

%retrieves wavelengths used in the measurement
%note: assumes all data collected with same WLs
WLi = unique(din(1).probe.link.type)';

%initialize DPF matrix (subjs x WLs)
DPF = zeros(length(din),length(WLi));

%iterates over all sessions
for n=1:length(din)

    age = age_all(n);
              
    %check if age info is not reliable
    if ~isnumeric(age) || age == 0
        age = 25; %assume 25y old
    end

    %iterate over wavelengths
    for w=1:size(DPF,2)
        
        %retrieve wavelength
        WL = WLi(w);    
        
        %apply DPF equation as suggested in the manuscript
        DPF(n,w) = 223.3 + 0.05624*(age^0.8493) + (-5.723)*(10^-7)*(WL^3) + 0.001245*(WL^2) + (-0.9025)*WL;
        
    end

end

end