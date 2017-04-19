function opt_delta = oxlay (V, t1, w, Cp, K)
%% flow chart implementation

A = 553.1; B = 600.8;
alpha = -5;
C = 0.0134;
m = 1;
n = 1.45;
Tw = 313;
ro = 7580;

ita = 0:0.1:1;
psi = 0:0.1:1;
sigma_diff = 1e5;            % assumed init. diff. between sigmaN & N_
opt_C0 = 2;                  % assumed optimum C0

for delta = 0.005:0.005:0.2
    itr = 1;
    for C0 = 2:0.1:10
        for phi = 5:0.1:45
            to_int = zeros(size(phi));
            k_chip = zeros(size(phi));
            Fc = zeros(size(phi));
            
            l = t1/sind(phi);
            Vs = V*cosd(alpha)/cosd(phi-alpha);
            GammaAB = 0.5*cos(alpha)/sin(phi)/cos(phi-alpha);
            Gamma_AB = C0*Vs/l;
            Eab = GammaAB/sqrt(3);
            E_ab = Gamma_AB/sqrt(3);
            Tab = Tw;
            while true
                Cp = 420 + 0.504*Tab;
                K = 52.61 - 0.0281*Tab;
                k_ab = sqrt(1/3)*(A+B*Eab^n)...
                    *(1+C*log(E_ab/E_0)*(1-((Tab-Tw)/(tm-tw))^m));
                Fs = k_ab*l*w;
                Rt = ro*Cp*V*t1/K;
                
                if Rt*tan(phi)<=10 || Rt*tan(phi) >0.04
                    beta = 0.5 - 0.35*log10(Rt*tan(phi));
                elseif Rt*tan(phi)>10
                    beta = 0.3 - 0.15*log10(Rt*tan(phi));
                end
                m_chip = ro*V*t1*w;
                Tsz = (1 - beta)*Fs*Vs/m_chip/Cp;
                NewTab = Tw + mean(ita * Tsz);

                if abs(Tab - NewTab) > 0.1
                    Tab = NewTab;
                else
                    break;
                end
            end
            n_eq = n*B*Eab^n/(A+B*Eab^n); 
            theta = atan(1 + 2*(pi/4 -phi) +C0*n);
            lambda = theta + alpha - phi;
            R = Fs/cosd(theta);
            F = R*sind(lambda);
            N = R*cosd(lambda);
            Fc(itr) = R*cosd(lambda-alpha);
            Ft = R*sind(lambda-alpha);
            t2 = t1*cosd(phi-alpha)/sind(phi);
            Vc = V*sind(phi)/cosd(phi-alpha);
            h = t1*sind(theta)/cosd(lambda)/sin(phi)*...
                (1+ (C0*n_eq/( 3*(1+2*(pi/4-phi)-C0*n_eq) )) );
            to_int(itr) = F/h/w;
            GammaM = h/delta/t2;
            GammaInt = 2*GammaAB + 0.5*GammaM;
            Eint = GammaInt/sqrt(3);
            E_int = V/sqrt(3)/delta/t2;
            Tc = Tw + Tsz;
            while true
                Cp = 420 + 0.504*Tc;
                K = 52.61 - 0.0281*Tc;
                del_Tc = F*Vc/m_chip/Cp;
                NewTc = Tw + Tsz + del_Tc;
                if abs(NewTc - Tc) > 0.1
                    Tc = NewTc;
                else
                    break;
                end
            end
            Rt = ro*Cp*V*t1/K;
            del_Tm = delta_Tc*...
                10^(0.06-0.195*delta*sqrt(Rt*t2/h)+0.5*log10(Rt*t2/h));
            Tint = Tw + Tsz + mean(psi*del_Tm);
            k_chip(itr) = sqrt(1/3)*(A+B*Eint^n)...
                *(1+C*log(E_int/E_0)*(1-((Tint-Tw)/(tm-tw))^m));
        end
        Fc_min = [Fc_min min(Fc)];
        % coomputing optimum phi
        diff = abs(to_int - k_chip);
        phi_index = logical(diff==min(diff((5:0.1:40)*10-50)));             %not sure
        fi = 5:0.1:45;
        opt_phi = fi(phi_index);
        disp(opt_phi);
        
        sigmaN_ = k_ab*(1+pi/2-2*alpha-2*C0*n_eq);
        sigmaN = N/h/w; 
        
        %computing optimum C0
        if abs(sigmaN - sigmaN_) < sigma_diff
            opt_C0 = C0;
            sigma_diff = abs(sigmaN - sigmaN_);
        end
    end
    Fc_min_m = [Fc_min_m min(Fc_min)];
    disp(opt_C0);
    
    new_Fc = [newFc exp];                   %no idea
end
%computing optimum delta
del_index = logical(new_Fc==min(new_FC));
delta = 0.005:0.005:0.2;
opt_delta = delta(del_index);
disp(opt_delta);
