function q_y = QY(k,T,Dy,state)

switch state
    case 'forecast'
        q_y =  -k.*Dyscheme(T,Dy,'backward');
    case 'correct'
        q_y =  -k.*Dyscheme(T,Dy,'forward');
end

end