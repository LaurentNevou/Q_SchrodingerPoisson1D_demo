%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% layers structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the conduction band offset in eV
% second column is the length of the layer in nm
% third column is the n doping volumique of that layer in 1e18cm-3 

% You have to put a resonable amount of doping! Otherwise, it will diverge 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InGaAs=0;
AlInAs=0.52;
AlInAs2=0.2;
meff   = 0.042;
Epsi   = 10;

M=[
AlInAs      5   0
AlInAs      1   2
AlInAs      5   0
InGaAs     6   0
AlInAs      5   0
InGaAs     7   0
AlInAs      5   0
AlInAs      1   5
AlInAs      5   0
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GaAs=0;
% AlGaAs40=0.36;
% meff   = 0.067;
% Epsi   = 10;
% 
% M=[
% AlGaAs40   5  0
% AlGaAs40   1  5
% AlGaAs40   5  0
% GaAs       8  0
% AlGaAs40   2  0
% GaAs       5  0
% AlGaAs40   5  0
% AlGaAs40   1  5
% AlGaAs40   5  0
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% InGaAs=0;
% GaAs=0.2;
% meff   = 0.05;
% Epsi   = 10;
% 
% M=[
% GaAs    10  0
% GaAs     1  5
% GaAs     5  0
% InGaAs   15  0
% GaAs     5  0
% GaAs     1  5
% GaAs    10  0
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GaN  = 0;
% AlN  = 1.8;
% meff = 0.22;
% Epsi = 10;
% DF   = 10;        % Electrical field discontinuity [MV/cm]
% 
% Lb = 3;           % barrier thickness [nm]
% Lw = 1;           % well thickness [nm]
% Ls = 0.2;         % doping spike thickness [nm] (in order to get the E-field)
% 
% Fb = +DF*(Lw+2*Ls)/(Lw+Lb+4*Ls)*1e6*1e2;   %[V/m]
% Fw = -DF*(Lb+2*Ls)/(Lw+Lb+4*Ls)*1e6*1e2;   %[V/m]
% 
% dopS = DF*1e6*1e2*Epsi*Epsi0/e;   % charge/m2 MUST BE added on the interface for GaN/AlN for Wurtzite
% dopV = dopS/(Ls*1e-9);            % charge/m3
% dopV = dopV*1e-6*1e-18;           % charge 1e18cm-3
% 
% M=[
% AlN     Ls   0
% AlN     Lb   0
% AlN     Ls   0
% 
% GaN     Ls  -dopV
% GaN     Lw   50
% GaN     Ls  +dopV
% 
% AlN     Ls   0
% AlN     Lb   0
% AlN     Ls   0
% 
% GaN     Ls  -dopV
% GaN     Lw   0
% GaN     Ls  +dopV
% ];
