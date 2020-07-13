function Ds = Dyscheme(s,Dy,var)
%%Y direction difference for s
Ds = zeros(size(s));
Ds(:,1) = (s(:,2)-s(:,1))/Dy;
Ds(:,end) = (s(:,end)-s(:,end-1))/Dy;
switch var
    case 'forward'
        Ds(:,2:end-1) = (s(:,3:end)-s(:,2:end-1))/Dy;
    case 'backward'
        Ds(:,2:end-1) = (s(:,2:end-1)-s(:,1:end-2))/Dy;
    case 'center'
        Ds(:,2:end-1) = (s(:,3:end)-s(:,1:end-2))/(2*Dy);
end
end