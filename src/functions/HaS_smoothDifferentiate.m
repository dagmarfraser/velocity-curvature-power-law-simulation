function [x0,x1,x2]=HaS_SmoothDifferentiate(x,N,m,samplerate)
% x: DATA (2nd dimension corresponds to time)
% N: Filter width (must be an odd integer. Larger N yields smoother
% results)
% m: Ahh..Accuracy (integer. larger m gives more accurate differentiation)
% samplerate: sampling rate (Hz). 
% x1: 1st order differentiated data (unit: 1/sec)
% x2: 2st order differentiated data (unit: 1/sec^2)
% Each differentiation decreases the length of the DATA by N-1.
% written by Ben Dongsung Huh 2011, UCSD, Salk Institute, 
% adopted by Jen Cook Lab - Dagmar Fraser July 2019

    if ~mod(N,2)
        error('N must be odd')
    else
        n = (N-1)/2;
    end

    h_d1=Diff_Filter1(n,m);
    % in line function Diff_Filter1
    
    x1=0;
    for k=1:N
        x1=x1+samplerate*h_d1(k)*x(:,k+(0:end-N));
    end
    
    if nargout==2
        x0=x(:,1+n:end-n);
    elseif nargout==3
        
%         N2=2*N-1;
%         n2=2*n; m2=2*m-1;
% %         N2=N;
% %         n2=n; m2=m;
%         h_d2=Diff_Filter2(n2,m2);
%         x2=0;
%         for k=1:N2
%             x2=x2+samplerate^2*h_d2(k)*x(:,k+(0:end-N2));
%         end
% %         x0=x(:,1+n:end-n);

        x2=0;
        for k=1:N
            x2=x2+samplerate*h_d1(k)*x1(:,k+(0:end-N));
        end
        x0=x(:,1+2*n:end-2*n);
        x1=x1(:,1+n:end-n);
    end
end

%% 1st Order Derivative Filter
function h_d1=Diff_Filter1(n,m)
    if m<1 || n-m<1, error('m is too small/large.'), end
    for k1=1:m
        M1(k1,:)=(1:n).^(2*k1-1);
    end
    for k2=1:n-m
        M1(m+k2,:)=(1:n).^(2*k2-1).*(-1).^(1:n);
    end
    h_d1=inv(M1)*[1/2; zeros(n-1,1)];
    h_d1=[-h_d1(end:-1:1); 0; h_d1];
end

%% 2nd Order Derivative Filter
function h_d2=Diff_Filter2(n,m)
    if m<1 || n-m<1, error('m is too small/large.'), end
    for k1=1:m
        M2(k1,:)=(1:n).^(2*k1);
    end
        M2(m+1,:)=(1:n).^(0).*((-1).^(1:n)-1);
    for k2=2:n-m
        M2(m+k2,:)=(1:n).^(2*k2-2).*(-1).^(1:n);
    end
    h_d2=inv(M2)*[1; zeros(n-1,1)];
    h_d2=[h_d2(end:-1:1); -2*sum(h_d2); h_d2];
end