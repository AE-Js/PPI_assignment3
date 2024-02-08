%% Description
% Function that calculates what modes are activated based on the rheology
% mode and the forcing mode. It calculates the activated modes based on
% some selection rules.
%%
function modes=next_coupling(mode0,modeR,order)
n0=mode0(1); 
m0=mode0(2);
ST=mode0(3); 
n1=modeR(1); 
m1=modeR(2);
kmode=1; 
num1 = 1/2*(n0+n1-abs(n0-n1));
num2 = (1/2*(n0+n1-abs(n0-n1))-1);
num_modes = num1 + num2 + 2;
modes_test = zeros(num_modes,4);
for j=0:num1
    m_new=m0+m1;
    n_new=abs(n0-n1)+2*j; 
    if ST==1
       ST_new=1; 
    else
       ST_new=2; 
    end
    if abs(m_new)<n_new+1
        modes_test(kmode,1)=n_new;    
        modes_test(kmode,2)=m_new;
        modes_test(kmode,3)=ST_new;
        modes_test(kmode,4)=order;
        kmode=kmode+1;
    end
end
for j=0:num2
    m_new=m0+m1;
    n_new=abs(n0-n1)+2*j+1;
    if ST==1
        ST_new=2; 
    else
        ST_new=1; 
    end
    if abs(m_new)<n_new+1 && abs(m_new)+abs(m0)+abs(m1)~=0
        modes_test(kmode,1)=n_new;    
        modes_test(kmode,2)=m_new;
        modes_test(kmode,3)=ST_new;
        modes_test(kmode,4)=order;
        kmode=kmode+1; 
    end
end
modes = modes_test(1:kmode-1,:);
end