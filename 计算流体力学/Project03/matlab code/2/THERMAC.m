function k = THERMAC(mu)
   %Thermal conductivity
   global c_p Pr;
    k = mu.*c_p./Pr;
end
