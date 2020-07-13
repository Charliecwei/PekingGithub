function phi = getphi3DP3N(lambda)

     la1 = lambda(:,1);
     la2 = lambda(:,2);
     la3 = lambda(:,3);
     la4 = lambda(:,4);





phi(:,1)=dian(la1);
phi(:,2)=dian(la2);
phi(:,3)=dian(la3);
phi(:,4)=dian(la4);

phi(:,5)=mian(la2,la3,la4);
phi(:,6)=mian(la3,la4,la1);
phi(:,7)=mian(la4,la1,la2);
phi(:,8)=mian(la1,la2,la3);

phi(:,[9,10])=bian(la1,la2);
phi(:,[11,12])=bian(la1,la3);
phi(:,[13,14])=bian(la1,la4);
phi(:,[15,16])=bian(la2,la3);
phi(:,[17,18])=bian(la2,la4);
phi(:,[19,20])=bian(la3,la4);

phi(:,21) = 9* (la1.^2+la2.^2+la3.^2+la4.^2)...
                   -11*(la1.^3+la2.^3+la3.^3+la4.^3)...
                   -72*(la1.*la2.*la3+la1.*la2.*la4+la1.*la3.*la4+la2.*la3.*la4);



   %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算顶点处值
function phis=dian(la)
phis=0.5*la.*(3*la-1).*(3*la-2);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算'1，2，3'为面上的基函数值
    function phis=mian(la1,la2,la3)
        phis=27*la1.*la2.*la3;
    end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算'1，2'边上基函数值
    function phis=bian(la1,la2)
        phis(:,1)=4.5*la1.*la2.*(3*la1-1);
        phis(:,2)=4.5*la2.*la1.*(3*la2-1);
    end

end