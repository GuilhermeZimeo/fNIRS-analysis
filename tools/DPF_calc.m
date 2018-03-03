function DPF = DPF_calc(din, age_all)
%This calculates the DPF based on subject's age and wavelengths used
%based on Scholkmann and Wolf 2013 - DOI: 10.1117/1.JBO.18.10.105004
%
%Last Updated: 2018-02-14
%
%By: Guilherme A. Zimeo Morais
%Contact: gazmorais@gmail.com

%retrieves wavelengths used in the measurement
%note: assumes all data collected with same WLs
WLi = unique(din(1).probe.link.type)';

%initialize DPF
DPF = zeros(length(din),length(WLi));

%iterates over all sessions
for n=1:length(din)

    if nargin < 2
        
        %attempt to retrieve Age from demographics table
        if din(n).demographics.values{2} > 0
            age = din(n).demographics.values{2};
            
        %if age not available (i.e. 0), initialize as 25    
        else
            age = 25;
        end
    
    %age array provided as input for the function   
    else
        age = age_all(n);
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