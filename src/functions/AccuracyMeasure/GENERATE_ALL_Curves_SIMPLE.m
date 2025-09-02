function [xy,theta_curve]=GENERATE_ALL_Curves_SIMPLE(Curve_Type, nu_list) 

theta_curve=[];

switch char(Curve_Type)
    case 'Line'
        [xy]=DRAW_Line();
    case {'Circle','Circle2'}
        [xy]=DRAW_Circle();
    case 'Line-Circle'
        [xy]=DRAW_Line_and_Circle();
    case {'Sinusoidal','Spiral'}
        [xy,theta_curve]=DRAW_Sinusoidal(char(nu_list));
    otherwise
        xy=[];
end

if ~isempty(xy)
    xy=100*[xy(1,:)+2; -xy(2,:)];
    xy=reshape(repmat(xy,[2,1]),[2,2*length(xy)]);   xy(:,[1,end])=[];
end
if ~isempty(theta_curve)
    theta_curve=reshape(repmat(theta_curve,[2,1]),[1,2*length(theta_curve)]);   theta_curve(:,[1,end])=[];
end

end
%%
function [path,theta]=DRAW_Sinusoidal(nu_list)
    load('Sinusoidal_Curves')
    disp(char(nu_list))
    switch char(nu_list)
        case '0',   k=1;
        case '2/33',k=2;
        case '2/5', k=3;
        case '4/5', k=4;
        case '4/3', k=5;
        case '2',   k=6;
        case '3',   k=7;
        case '4',   k=8;
        case '6',   k=9;
        otherwise error()
    end
    path=path_all{k}; 
    theta=theta_curve_all{k};
end
%% DRAW_Line
function [xy]=DRAW_Line()
    xy=0.3937*9.5*[-1.5, 1.5; 0.65,0.65];  %     xy=0.9*0.3937*[real(z);imag(z)]*9.5/max(imag(z));

end

%% DRAW_Circle
function [xy]=DRAW_Circle()

    xy=0.3937*9.5*[cos(linspace(0,2*pi,301)); -0.2+sin(linspace(0,2*pi,301))];  %     xy=0.9*0.3937*[real(z);imag(z)]*9.5/max(imag(z));
    xy=0.8*xy;
end

%% DRAW_Line_and_Circle
function [xy]=DRAW_Line_and_Circle(subtype)

        xy{1}=0.3937*9.5*[-1.5, 1.5; 0.65,0.65];  %     xy=0.9*0.3937*[real(z);imag(z)]*9.5/max(imag(z));

        xy{2}=0.3937*9.5*[cos(linspace(0,2*pi,301)); -0.2+sin(linspace(0,2*pi,301))];  %     xy=0.9*0.3937*[real(z);imag(z)]*9.5/max(imag(z));
        xy{2}=0.8*xy{2};
end

