%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% last update 10June2020, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program computes the Schrodinger-Poisson equation in heterostructures
% In order to keep the code fast but still usefull, the mass is kept constant
% all over the structure. It means that meff should be set at the value of the
% well. Obviously, the non-parabolicity of the band are also not considered in 
% the Schrodinger solver and the density of states.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the code doesn t converge:
% -> decrease the doping
% -> increase the resolution dz
% -> increase the temperature (T=0K is very bad while T=10K is already much better)
% -> increase the amount of loops, Nloops
% -> The Newton-Raphson algorithm is not really helping... It is slower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ref: Newton-Raphson algorithm
% "Newton-Raphson solution of Poisson's equation in a pn diode"
% R. A. Jabr, M. Hamad and Y. M. Mohanna
% International Journal of Electrical Engineering Education
%  Volume: 44 issue: 1, page(s): 23-33, Issue published: January 1, 2007
% https://doi.org/10.7227/IJEEE.44.1.3
% https://journals.sagepub.com/doi/10.7227/IJEEE.44.1.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h    = 6.62606896E-34;              %% Planck constant [J.s]
hbar = h/(2*pi);
e    = 1.602176487E-19;             %% electron charge [C]
m0   = 9.10938188E-31;              %% electron mass [kg]
Epsi0= 8.854187817620E-12;          %% Vaccum dielectric constant [F/m]
kB   = 1.3806488E-23;               %% Boltzmann's constant [J/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloops = 50;                  % number of loops
NNewRap= 100;                 % number of the loop at which starts the Newton-Raphson algorithm
n      = 5;                   % number of solution asked per model
ScF    = 0.1;                 % scaling factor to plot the wave function [Without Dimension]
dz     = 1e-10;               % resolution of the grid [m]
F0     = 0;%-2e7;%-6e7;       % Electric field [Volt/meter]
T      = 300;                 % Temperature [Kelvin], react on the Fermi function only

plot_density=1;               % Activate the plot 0 or 1
plot_convergence=0;           % Activate the plot 0 or 1
plot_field=0;                 % Activate the plot 0 or 1
plot_Vbending=0;              % Activate the plot 0 or 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt   = M(:,2)*1e-9;         % conversion of the length from Angstrom to meter
Dopt = M(:,3)*1e18*1e6;     % n doping conversion from cm-3 to m-3
CBOt = M(:,1);              % Conduction Band Offset [eV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z, the potential V0 and that values that are needed

z=0; V0=CBOt(1); Dop=Dopt(1);

for i=1:length(zt)
    t=zt(i);
    zv= (z(end)+dz): dz : (z(end)+dz)+t;
    z=[z zv];
    V0  = [ V0     ones(size(zv)) * CBOt(i)  ];
    Dop = [ Dop    ones(size(zv)) * Dopt(i)  ];
end

V0=V0-min(V0);             % Shift the band in order to get the bottom of the well at zero
V0=(F0*z)+V0;              % adding the electric field to the potential

Ntott=Dopt.*zt;
Ntot=sum(Ntott);   % total number of charges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Building operator matrix for Newton-Raphson Algorithm %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nz  = length(z);
DZ2 = (-2)*diag(ones(1,Nz)) + (1)*diag(ones(1,Nz-1),-1) + (1)*diag(ones(1,Nz-1),1);
DZ2 = DZ2/dz^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Starting of the Poisson s loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vs=zeros(size(z)); Vsold=Vs;
ntot=0;
nloop=1;

ErrVec=1;
sumVtotVec=1;

if Dopt==0
    Nloops=2;
end

while nloop<Nloops
    
    nloop
    x = 1;
    Vbending=Vs*x + Vsold*(1-x);
    Vtot=V0+Vbending;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% schrodinger solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    [Ec,psic] = Schroed1D_FEM_f(z,Vtot,meff,n);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Here, I re-define the energy grid in order optimize the meshing
    dE1=1e-4; dE2=1e-2;

    E1 = Ec(1):dE1:Ec(1)+0.1 ;
    E2 = E1(end):dE2:max(Vtot);
    En=sort([E1 E2]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ro=[];
    for i=1:length(Ec)
        ro( En>Ec(i),i) = e*meff*m0/(pi*hbar^2);
        ro( En<Ec(i),i) = 0;
    end

    [Ef,NN,roEf]=find_Ef_f(Ec,En,ro,Ntot,T);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ntot2 = repmat(NN,[length(z) 1]).*abs(psic).^2 ;
    
    ntot = sum(ntot2,2)' - Dop;  % remove the positive charges (ions)
    
    
    if (nloop>NNewRap-1)
        dVV=1e-5;
        dVtot=V0+Vbending+dVV;
           
        [dEc,dpsic] = Schroed1D_FEM_f(z,dVtot,meff,n);
        E11 = dEc(1):dE1:dEc(1)+0.1 ;
        E22 = E11(end):dE2:max(dVtot);
        dEn=sort([E11 E22]);
        dro=[];
        for i=1:length(dEc)
            dro( dEn>dEc(i),i) = e*meff*m0/(pi*hbar^2);
            dro( dEn<dEc(i),i) = 0;
        end
        [dEf,dNN,droEf]=find_Ef_f(dEc,dEn,dro,Ntot,T);
        
        ntot2 = repmat(dNN,[length(z) 1]).*abs(dpsic).^2 ;
        dntot = sum(ntot2,2)' - Dop - ntot;  % remove the positive charges (ions)
    
    end
    
    if nloop<NNewRap        % => Damping injection method
    
        %%%%%%%%%%%%%%%%%%%%%%%%%% Electrical Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F = e*cumtrapz(z,ntot)./(Epsi0*Epsi);
        MF  = trapz(z,F)/(z(end)-z(1));  % MF=mean(F) on a nonlinear grid z
        F = F-MF;
        %%%%%%%%%%%%%%%%%%%%%%%%%%% New Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Vsold = Vs;
        Vs    = -cumtrapz(z,F);
            
    elseif nloop>NNewRap-1  % => Newton Raphson algorithm
            
        FFF = -(DZ2*Vtot')*Epsi*Epsi0/e - ntot';  %% FFF is the fonctionnel that should converge to zero
        JJJ = sparse( -DZ2*Epsi*Epsi0/e - diag(dntot/dVV) ); %% JJJ is the Jacobian
        
        Vsold = Vs;
        %Vs = Vtot - (inv(JJJ)*FFF)';
        Vs = Vtot - (JJJ\FFF)';
        Vs=Vs-Vs(1);
        
        F=-e*cumtrapz(z,ntot)/Epsi/Epsi0; %% Computes the E-field for the plot
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Convergence analysis/plot %%%%%%%%%%%%%%%%%%%%%%%%%%

    Err = abs(  1 - sumVtotVec(end)/sum(Vs)  );
    sumVtotVec(nloop) = sum(Vs);
    ErrVec = [ErrVec Err];

    nloop=nloop+1;
    
    if Err<1e-10
       break 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Scaling and shifting the wavefunctions %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Ec)
    PSIc(:,i)=abs(psic(:,i)).^2/max(abs(psic(:,i)).^2)*ScF + Ec(i); % normalisation for the plotting
end

for ii=1:length(Ec)
    psic_M = repmat(psic(:,ii),[1,length(En)])';
    ROEf(:,:,ii) = repmat(roEf(:,ii),[1 length(z)]) .* abs(psic_M.^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 1000 700],'color','w');
subplot(1,1,1,'fontsize',15)
hold on;grid on;
col=colormap(jet);

if plot_density==1
    grid off
    pcolor(z*1e9,En,sum(ROEf,3)*1e-6 )
    set(gca,'color',col(1,:))
    shading flat
    hcb=colorbar;
    title(hcb,'\fontsize{8}cm-3')
    
    plot(z*1e9,V0,  'w--','linewidth',1)
    plot(z*1e9,Vtot,'w-' ,'linewidth',1)
    
elseif plot_density==0
    plot(z*1e9,V0,  'b--','linewidth',1)
    plot(z*1e9,Vtot,'b-' ,'linewidth',1)
end

plot([z(1) z(end)]*1e9,[1 1]*Ef,'g','linewidth',1)
text(z(end)*1e9*0.95,Ef+0.01,'\color{green}Fermi')

for i=1:length(Ec)
    plot(z*1e9,PSIc(:,i),'color','r','linewidth',1)
end

xlabel('z (nm)')
ylabel('Energy (eV)')

title(strcat('\fontsize{12}T=',num2str(T),'K; meff=',num2str(meff),'; Epsilon=',num2str(Epsi),'; dz=',num2str(dz*1e9),'nm; Ntot=',num2str(Ntot*1e-4,'%.1e'),'cm-2'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_convergence==1

    figure
    semilogy(1:nloop,ErrVec,'bo-')
    hold on; grid on;
    xlabel('Cycles')
    ylabel('Convergence (norm. units)')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_field==1
    
    figure
    hold on;grid on;
    [AX,H1,H2]=plotyy(z*1e9,F*1e-2*1e-3,z*1e9,Dop*1e-18*1e-6);
        
    set(H1,'color','r')
    set(H2,'color','b')
    
    xlabel('z (nm)')
    ylabel(AX(1),'E- field (kV/cm)','color','red')
    ylabel(AX(2),'Doping (1e18 cm-3)','color','blue')
    
    set(AX(1),'ycolor','red')
    set(AX(2),'ycolor','blue')
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Vbending==1

    figure
    hold on; grid on;

    [AX,H1,H2]=plotyy(z*1e9,Vbending,z*1e9,ntot*1e-18*1e-6);

    set(H1,'color','r')
    set(H2,'color','b')
    
    xlabel('z (nm)')
    ylabel(AX(1),'Vbending (eV)','color','red')
    ylabel(AX(2),'ntot (1e18 cm-3)','color','blue')
    
    set(AX(1),'ycolor','red')
    set(AX(2),'ycolor','blue')

end