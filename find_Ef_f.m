function[Ef,NN,roEf]=find_Ef_f(Ec,E,ro,Ntot,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e=1.602176487E-19;              %% electron charge [C]
kB = 1.3806488E-23;             %% Boltzmann's constant [J/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%% Computes the Fermi level at any T %%%%%%%%%%%%%%%%%%%%%%%%%%
if T==0
   T=1e-10; 
end

Ef=Ec(1);
Fermi= 1./(1+exp((E-Ef)/(kB*T/e))) ;
NtotX=0;

roEf=[];

for i=1:length(Ec)
    roEf(:,i)=ro(:,i).*Fermi';
    NN(i)= trapz(E,roEf(:,i),1)  ;
    NtotX = NtotX + NN(i) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Now, it will try to get as close as posible to the real Ef with an
% error of 0.1% by dichotomy
ddE=0.01; % eV
Ef1=Ef;
Ef2=Ef1+ddE;

while  abs(NtotX - Ntot)/Ntot > 0.001  % find the Fermi level at any temperature

    NN(length(Ec))=0;
    if NtotX > Ntot
        Ef  = Ef - abs(Ef1-Ef2)/2 ;  
        Ef1 = Ef ;
    else
        Ef  = Ef + abs(Ef1-Ef2)/2 ;
        Ef2 = Ef ; 
    end
    
    Fermi = 1./(1+exp((E-Ef)/(kB*T/e))) ;  % Fermi Dirac distribution function
    
    NtotX=0;
    for i=1:length(Ec)
          roEf(:,i) = ro(:,i).*Fermi';
          NN(i)= trapz(E,roEf(:,i),1)  ;
          NtotX = NtotX + NN(i) ;
    end
end

end