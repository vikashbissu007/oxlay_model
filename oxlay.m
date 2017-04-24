function [Fc, Delta] = oxlay (V, t1, DoC)
%% flow chart implementation



w = 1.402*DoC;

A = 553.1; B = 600.8;
alpha = -5;
C = 0.0134;
m = 1;
n = 0.234;
Tw = 313;
ro = 7580;
Tm = 1701;

ita = 0:0.1:1;
psi = 0:0.1:1;
sigma_diff = 1e5;            % assumed init. diff. between sigmaN & N_
opt_C0 = 2;                  % assumed optimum C0
diff = 1e5;
opt_phi = 5;

E_0 = 1;

for delta = 0.005:0.005:0.2
    Fc_min_m = [];
    for C0 = 2:0.1:10
        Fc_min = [];
        Fc = [];

        for phi = 5:0.1:45
            itr = 1;
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
                    *(1+C*log(E_ab/E_0)*(1-((Tab-Tw)/(Tm-Tw))^m));
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
            to_int = F/h/w;
            GammaM = h/delta/t2;
            GammaInt = 2*GammaAB + 0.5*GammaM;
            Eint = GammaInt/sqrt(3);
            E_int = V/sqrt(3)/delta/t2;
            Tc = Tw + Tsz;
            while true
                Cp = 420 + 0.504*Tc;
                K = 52.61 - 0.0281*Tc;
                delta_Tc = F*Vc/m_chip/Cp;
                NewTc = Tw + Tsz + delta_Tc;
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
            k_chip = sqrt(1/3)*(A+B*Eint^n)...
                *(1+C*log(E_int/E_0)*(1-((Tint-Tw)/(Tm-Tw))^m));
            
            % computing optimum phi
            if phi<=40 && abs(to_int - k_chip) < diff
                opt_phi = phi;
                diff = abs(to_int - k_chip);
            end
%             disp(opt_phi);
        end
        Fc_min = [Fc_min min(Fc)];
        
        sigmaN_ = k_ab*(1+pi/2-2*alpha-2*C0*n_eq);
        sigmaN = N/h/w; 
        
        %computing optimum C0
        if abs(sigmaN - sigmaN_) < sigma_diff
            opt_C0 = C0;
            sigma_diff = abs(sigmaN - sigmaN_);
        end
    end
    Fc_min_m = [Fc_min_m min(Fc_min)];
%     disp(opt_C0);
    
end
disp(min(Fc_min_m));
delta = 0.005:0.005:0.2;
Delta = delta(find(Fc_min_m == min(Fc_min_m)));
