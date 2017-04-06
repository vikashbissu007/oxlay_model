function Fc = oxlay (V, tl, w, alpha)

%% constants for AISI 1045
% values for S, K, Tr are to be replaced by actual values
A = 553.1; B = 600.8; C = 0.0134; ...
    n = 0.234; m = 1; Tm = 1460;...
    S = 1; K = 1; Tr = 1;

%% flow chart implementation
% exp is arbitary variable and to be replaced by expression
% for each variable
%
% first value for C0 is to be replaced here
%
C0 = 1; exp = 1;

for delta = 0.001:0.001:0.005
    while true
        for phi = 4:delta:45-delta
            Tab = Tr;
            while true
                E_ab = exp;
                Eab = exp;
                k_ab = exp;
                beta = exp;
                Fs = exp;
                NewTab = exp;
                if Tab == NewTab
                    break;
                else
                    Tab = NewTab;
                end
            end
            to_int = exp;
            sigmaN = exp;
            sigmaN_ = exp;
            E_int = exp;
            Eint = exp;
            Tint = exp;
            k_chip = exp;
            if k_chip == to_int
                break;
            end
        end
        if sigmaN == sigmaN_
            break;
        else
            C0 = exp;
        end
    end
    delta_final = exp;
    if delta == delta_final
        break;
    end
end

%% determine delta with min Fc
% no idea :(
