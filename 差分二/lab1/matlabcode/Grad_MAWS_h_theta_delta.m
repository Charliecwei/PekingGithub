function [Fu,Grad_Fu] = Grad_MAWS_h_theta_delta(u,h,G)
%%To calculate the Grad_u of MAWS_h_theta_delta

L = size(G,1);

Nu = length(u);

for l = 1:L
e = G(l,:);
[Max_deltaeu,Gradu_Max_deltaeu] = GraduOfMax_deltaeu(u,h,e);
  
    if l == 1
           Grad_Fu = Gradu_Max_deltaeu;
           Fu = Max_deltaeu;
    else
        Part_xmin_delta = Partial_xmin_delta(Max_deltaeu,Fu);
        Part_ymin_delta = Partial_xmin_delta(Fu,Max_deltaeu);

        Grad_Fu = sparse(1:Nu,1:Nu,Part_xmin_delta,Nu,Nu)*Gradu_Max_deltaeu...
                       + sparse(1:Nu,1:Nu,Part_ymin_delta,Nu,Nu)*Grad_Fu;

        Fu = min_delta(Max_deltaeu,Fu);
    end



end
    
    
 
          
 

    
end




function [Max_deltaeu,Gradu_Max_deltaeu] = GraduOfMax_deltaeu(u,h,e)
%%Calculate the product plus and its gradient along 
%%the perpendicular direction of e and e£¬ Sparse matrix storage
Nu = length(u);
e_v = [-e(2),e(1)]; 
Dletaeu_delta_plus = U_plus_delta(Delta_e(u,e,h));
Dletae_vu_delta_plus = U_plus_delta(Delta_e(u,e_v,h));
Gradu_Deltaeu = Grad_uOfDelta_eu(u,e,h);
Gradu_Deltae_vu = Grad_uOfDelta_eu(u,e_v,h);

Partial_xDeltaeu_plux = Partial_xmax_delta(Dletaeu_delta_plus,0);
Partial_xDeltae_vu_plux = Partial_xmax_delta(Dletae_vu_delta_plus,0);

S_1  = sparse(1:Nu,1:Nu,Dletae_vu_delta_plus.*Partial_xDeltaeu_plux,Nu,Nu)*Gradu_Deltaeu;
S_2  = sparse(1:Nu,1:Nu,Dletaeu_delta_plus.*Partial_xDeltae_vu_plux,Nu,Nu)*Gradu_Deltae_vu;

Gradu_Max_deltaeu = S_1+S_2;

Max_deltaeu = Dletaeu_delta_plus.*Dletae_vu_delta_plus;


end
    
