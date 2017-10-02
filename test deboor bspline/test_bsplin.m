clear all
clc
close all
%%
% M_4=[-1 3 -3 1; ...
%     3 -6 3 0; ...
%     -3 0 3 0; ...
%     1 4 1 0];
% 
% k_4=1/6;
% 
% P=[0 10 7 0 15 15 0];
% 
% ti=[0 1 2 3 4 5 6];
% tstep=[0:0.01:1];
% 
% M1=k_4*M_4*[P(1) P(1) P(1) P(2)]';
% M2=k_4*M_4*[P(1) P(1) P(2) P(3)]';
% M3=k_4*M_4*P(1:4)';
% M4=k_4*M_4*P(2:5)';
% M5=k_4*M_4*P(3:6)';
% 
% p1=[tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M1;
% p2=[tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M2;
% p3=[tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M3;
% p4=[tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M4;
% p5=[tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M5;
% 
% figure;
% clf
% hold on
% plot(tstep-2,p1,'k');
% plot(tstep-1,p2,'r');
% plot(tstep,p3);
% plot(tstep+1,p4,'r');
% plot(tstep+2,p5,'k');
% hold off
% 
% figure;
% clf
% hold on
% plot(tstep-2,p1,'k');
% plot(tstep-1,p2,'r');
% plot(tstep,p3);
% plot(tstep*3+1,p4,'r');
% plot(tstep+4,p5,'k');
% hold off

%%
% M_6_=([
%   1  -5 10  -10 5   -1;
%   0  5  -20 30  -20 5;
%   0  10 -20 0   20  -10;
%   0  10 20  -60 20  10;
%   0  5  50  0   -50 -5;
%   0  1  26  66  26  1
%  ]);
% 
% M_6=M_6_(:,6:-1:1);
% 
% 
% k_6=1/120;
% 
% P=[0 10 7 0 15 15 0];
% 
% ti=[0 1 2 3 4 5 6];
% tstep=[0:0.01:1];
% 
% M1=k_6*M_6*[P(1) P(1) P(1) P(1) P(1) P(2)]';
% M2=k_6*M_6*[P(1) P(1) P(1) P(1) P(2) P(3)]';
% M3=k_6*M_6*[P(1) P(1) P(1) P(2) P(3) P(4)]';
% M4=k_6*M_6*[P(1) P(1) P(2) P(3) P(4) P(5)]';
% M5=k_6*M_6*[P(1) P(2) P(3) P(4) P(5) P(6)]';
% 
% p1=[tstep'.^5 tstep'.^4 tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M1;
% p2=[tstep'.^5 tstep'.^4 tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M2;
% p3=[tstep'.^5 tstep'.^4 tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M3;
% p4=[tstep'.^5 tstep'.^4 tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M4;
% p5=[tstep'.^5 tstep'.^4 tstep'.^3 tstep'.^2 tstep'.^1 tstep'.^0]*M5;
% 
% 
% figure;
% clf
% hold on
% plot(tstep-2,p1,'k');
% plot(tstep-1,p2,'r');
% plot(tstep,p3);
% plot(tstep+1,p4,'r');
% plot(tstep+2,p5,'k');
% hold off
%%
% clear all
% close all
% P=[0 10 7 0 15 15 0];
% ti=[0 1 2 3 4 5 6];
% tstep_1=[0:0.1:1-0.01]';
% tstep_2=[0:0.1:1-0.01]';
% tstep_3=[0:0.1:1-0.01]';
% 
% % N0=1;
% n=1;
% % N.t1=one;
% % N.ti(n-2)=[1 0 0];
% % N.ti(n-1)=[1 0 0];
% % ti=1;
% N=zeros(size(tstep_1,1)+size(tstep_2,1)+size(tstep_3,1),3);
% N(1:size(tstep_1,1),1)=ones(1:size(tstep_1,1),1);
% N(size(tstep_1,1)+1:size(tstep_1,1)+size(tstep_2,1),2)=ones(1:size(tstep_2,1),1);
% N(size(tstep_1,1)+size(tstep_2,1)+1:size(tstep_1,1)+size(tstep_2,1)+size(tstep_3,1),3)=ones(1:size(tstep_3,1),1);
% for i=1:n
%     N(:,1)=N(:,1).*(tstep_1-ti(1))/(ti(2)-ti(1))+ N(:,2).*(ti(1+n+1)-tstep_1)/(ti(1+n+1)-ti(1+1));
%     N(:,2)=N(:,2).*(tstep_2-ti(2))/(ti(3)-ti(2))+ N(:,3).*(ti(2+n+1)-tstep_2)/(ti(2+n+1)-ti(2+1));
%     N(:,3)=N(:,3).*(tstep_3-ti(3))/(ti(4)-ti(3));
% end

% % %%
% % P=[0 10 7 0 15 15 0 5 10];
% % 
% % ti=[0 1 2 3 4 5 6 7 8 9 10 11 12];
% % ti(6:end)=ti(6:end)+1;
% % 
% % n=3;
% % i=4;
% % 
% % tstep=[ti(i):0.01:ti(i+1)-0.01];
% % 
% % v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% % v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% % v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% % v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% % 
% % v=[v1' v2' v3' v4'];
% % 
% % % p=v*P(i:i+3)';
% % P(i-3:i)'
% % p=v*P(i-3:i)';
% % 
% % tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% % 
% % s=tm*p;
% % 
% % figure;
% % clf
% % hold on
% % plot(ti(1:9),P,'o')
% % plot(tstep,s)
% % hold off
% % %%
% % n=3;
% % i=5;
% % 
% % tstep=[ti(i):0.01:ti(i+1)-0.01];
% % 
% % v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% % v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% % v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% % v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% % v=[v1' v2' v3' v4'];
% % 
% % % p=v*P(i:i+3)';
% % p=v*P(i-3:i)';
% % 
% % tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% % 
% % s=tm*p;
% % hold on
% % plot(tstep,s,'r')
% % hold off
% % %%
% % n=3;
% % i=6;
% % 
% % tstep=[ti(i):0.01:ti(i+1)-0.01];
% % 
% % v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% % v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% % v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% % v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% % v=[v1' v2' v3' v4'];
% % 
% % % p=v*P(i:i+3)';
% % p=v*P(i-3:i)';
% % 
% % tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% % 
% % s=tm*p;
% % hold on
% % plot(tstep,s,'k')
% % hold off
% % % %%
% % % n=
% % % v=[-(ti(n)^3) (3*ti(n)^2) -(3*ti(n)) 1]/((-ti(n)+ti(n+1))*(-ti(n)+ti(n+2))*(-ti(n)+ti(n+2)));
% % %%
% % n=3;
% % i=3;
% % 
% % tstep=[ti(i):0.01:ti(i+1)-0.01];
% % 
% % v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% % v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
% %     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% % v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% % v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% % 
% % v=[v2' v3' v4'];
% % 
% % % p=v*P(i:i+3)';
% % p=v*P(i-2:i)';
% % 
% % tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% % 
% % s=tm*p;
% % hold on
% % plot(tstep,s,'k')
% % hold off
% % %%
% % n=3;
% % i=2;
% % 
% % tstep=[ti(i):0.01:ti(i+1)-0.01];
% % 
% % v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
% %     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% % v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% % 
% % v=[v3' v4'];
% % 
% % % p=v*P(i:i+3)';
% % p=v*P(i-1:i)';
% % 
% % tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% % 
% % s=tm*p;
% % hold on
% % plot(tstep,s,'r')
% % hold off
% % %%
% % n=3;
% % i=1;
% % 
% % tstep=[ti(i):0.01:ti(i+1)-0.01];
% % 
% % v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% % 
% % v=[v4'];
% % 
% % % p=v*P(i:i+3)';
% % p=v*P(i:i)';
% % 
% % tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% % 
% % s=tm*p;
% % hold on
% % plot(tstep,s,'b')
% % hold off


%%
clear all
P=[0 10 7 0 15 15 0 5 10];

ti=[0 1 2 3 4 5 6 7 8 9 10 11 12];
ti(6:end)=ti(6:end)+1;

n=3;
i=4;

tstep=[ti(i):0.01:ti(i+1)-0.01];
tstep=[ti(i):0.01:ti(i+1)];


v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
 -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));

v=[v1' v2' v3' v4'];

% p=v*P(i:i+3)';

tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

% syms a b c d
% [sol_a, sol_b, sol_c, sol_d]=solve([tm(1,:);tm_d(1,:);tm_dd(1,:);tm(end,:)]*v*[a;b;c;d]==[0;0;0;1]);

p=v*P(i-3:i)';
% p=v*[sol_a; sol_b; sol_c; sol_d];
s=tm*p;
% s=tm_d*p;
% s=tm_dd*p;

figure;
clf
hold on
plot(ti(1:9),P,'o')
plot(tstep,s)
hold off
%%%%
%%
n=3;
i=5;

tstep=[ti(i):0.01:ti(i+1)-0.01];

v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
 -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
v=[v1' v2' v3' v4'];

% p=v*P(i:i+3)';
p=v*P(i-3:i)';


tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

s=tm*p;
% s=tm_d*p;
% s=tm_dd*p;

hold on
plot(tstep,s,'r')
hold off
%%
n=3;
i=6;

tstep=[ti(i):0.01:ti(i+1)-0.01];

v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
 -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
v=[v1' v2' v3' v4'];

% p=v*P(i:i+3)';
p=v*P(i-3:i)';


tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

s=tm*p;
% s=tm_d*p;
% s=tm_dd*p;

hold on
plot(tstep,s,'k')
hold off

%%
p=v*[P(i-3:i-1)';-15];


tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

s=tm*p;
% s=tm_d*p;
% s=tm_dd*p;

hold on
plot(tstep,s,'m')
hold off
%%
n=3;
i=6;

tstep=[ti(i):0.01:ti(i+1)-0.01];

% ti=ti-ti(i);
v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
    1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
 -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
    -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
v=[v1' v2' v3' v4'];

% p=v*P(i:i+3)';
p=v*P(i-3:i)';


tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

syms a b c d e f
[sol_a, sol_b, sol_c]=solve([tm(1,:);tm_d(1,:);tm_dd(1,:)]*v*[a;b;c;-15]==[tm(1,:)*p;tm_d(1,:)*p;tm_dd(1,:)*p]);

% p=v*P(i-3:i)';
p=v*[sol_a; sol_b; sol_c; 15];

s=tm*p;
% s=tm_d*p;
% s=tm_dd*p;

hold on
plot(tstep,s,'*k')
hold off
%%
clear all
P=[0 10 7 0 15 15 0 5 10];

ti=[0 1 2 3 4 5 6 7 8 9 10 11 12];
% ti(6:end)=ti(6:end)+1;

n=5;
i=6;

tstep=[ti(i):0.01:ti(i+1)-0.01];

% v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% v1=[];
% v2=[];
% v3=[];
% v4=[];
for j=1:n+1
    run(['v' num2str(j) '_n5']);
end
v=[v1' v2' v3' v4' v5' v6'];

% p=v*P(i:i+3)';
tm=[];
for n=0:n
    tm=[tm tstep'.^(n)];
end
tm_d=[];
for n=0:n
    tm_d=[tm_d n*tstep'.^(n-1)];
end
tm_dd=[];
for n=0:n
    tm_dd=[tm_dd n*(n-1)*tstep'.^(n-2)];
end
% tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
% tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

% syms a b c d
% [sol_a, sol_b, sol_c, sol_d]=solve([tm(1,:);tm_d(1,:);tm_dd(1,:);tm(end,:)]*v*[a;b;c;d]==[0;0;0;1]);

p=v*P(i-n:i)';
% p=v*[sol_a; sol_b; sol_c; sol_d];
% s=tm*p;
% s=tm_d*p;
s=tm_dd*p;

figure;
clf
hold on
plot(ti(1:9),P,'o')
plot(tstep,s)
hold off
%%%%
%%
n=5;
i=7;

tstep=[ti(i):0.01:ti(i+1)-0.01];

% v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% v1=[];
% v2=[];
% v3=[];
% v4=[];
for j=1:n+1
    run(['v' num2str(j) '_n5']);
end
v=[v1' v2' v3' v4' v5' v6'];

% p=v*P(i:i+3)';
tm=[];
for n=0:n
    tm=[tm tstep'.^(n)];
end
tm_d=[];
for n=0:n
    tm_d=[tm_d n*tstep'.^(n-1)];
end
tm_dd=[];
for n=0:n
    tm_dd=[tm_dd n*(n-1)*tstep'.^(n-2)];
end
% tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
% tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

% syms a b c d
% [sol_a, sol_b, sol_c, sol_d]=solve([tm(1,:);tm_d(1,:);tm_dd(1,:);tm(end,:)]*v*[a;b;c;d]==[0;0;0;1]);

p=v*P(i-n:i)';
% p=v*[sol_a; sol_b; sol_c; sol_d];
% s=tm*p;
% s=tm_d*p;
s=tm_dd*p;

hold on

plot(tstep,s,'r')
hold off
%%%%
%%
n=5;
i=8;

tstep=[ti(i):0.01:ti(i+1)-0.01];

% v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% v1=[];
% v2=[];
% v3=[];
% v4=[];
for j=1:n+1
    run(['v' num2str(j) '_n5']);
end
v=[v1' v2' v3' v4' v5' v6'];

% p=v*P(i:i+3)';
tm=[];
for n=0:n
    tm=[tm tstep'.^(n)];
end
tm_d=[];
for n=0:n
    tm_d=[tm_d n*tstep'.^(n-1)];
end
tm_dd=[];
for n=0:n
    tm_dd=[tm_dd n*(n-1)*tstep'.^(n-2)];
end
% tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
% tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

% syms a b c d
% [sol_a, sol_b, sol_c, sol_d]=solve([tm(1,:);tm_d(1,:);tm_dd(1,:);tm(end,:)]*v*[a;b;c;d]==[0;0;0;1]);

p=v*P(i-n:i)';
% p=v*[sol_a; sol_b; sol_c; sol_d];
% s=tm*p;
% s=tm_d*p;
s=tm_dd*p;

hold on

plot(tstep,s,'k')
hold off
%%%%
%%
n=5;
i=8;

tstep=[ti(i):0.01:ti(i+1)-0.01];

% v1=[(ti(i+1)^3) -(3*ti(i+1)^2) (3*ti(i+1)) -(1)]/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)));
% v2=[-((ti(i-2)*ti(i+1)^2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(ti(i-1)*ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i)*ti(i+2)^2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     (2*ti(i-2)*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+ti(i+1)^2/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+(ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i+1)*ti(i+2))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(2*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)^2/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     -(ti(i-2)/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))))-(2*ti(i+1))/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))-ti(i-1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i+2)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(2*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2))) ...
%     1/((-ti(i-2)+ti(i+1))*(-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1)))+1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))];
% v3=[(ti(i-1)^2*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+(ti(i-1)*ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(ti(i)^2*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%  -(ti(i-1)^2/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-(2*ti(i-1)*ti(i+1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))-(ti(i-1)*ti(i))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i-1)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-(ti(i)*ti(i+2))/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-ti(i)^2/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))-(2*ti(i)*ti(i+3))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     (2*ti(i-1))/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i+1)/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2)))+ti(i-1)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+ti(i+2)/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))+(2*ti(i))/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))+ti(i+3)/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3))) ...
%     -(1/((-ti(i-1)+ti(i+1))*(-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))))-1/((-ti(i)+ti(i+1))*(-ti(i-1)+ti(i+2))*(-ti(i)+ti(i+2)))-1/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)))];
% v4=[-(ti(i)^3) (3*ti(i)^2) -((3*ti(i))) 1]/((-ti(i)+ti(i+1))*(-ti(i)+ti(i+2))*(-ti(i)+ti(i+3)));
% v1=[];
% v2=[];
% v3=[];
% v4=[];
for j=1:n+1
    run(['v' num2str(j) '_n5']);
end
v=[v1' v2' v3' v4' v5' v6'];

% p=v*P(i:i+3)';
tm=[];
for n=0:n
    tm=[tm tstep'.^(n)];
end
tm_d=[];
for n=0:n
    tm_d=[tm_d n*tstep'.^(n-1)];
end
tm_dd=[];
for n=0:n
    tm_dd=[tm_dd n*(n-1)*tstep'.^(n-2)];
end
% tm=[(tstep.^0)' (tstep.^1)' (tstep.^2)' (tstep.^3)'];
% tm_d=[0*(tstep.^0)' 1*(tstep.^0)' 2*(tstep.^1)' 3*(tstep.^2)'];
% tm_dd=[0*(tstep.^0)' 0*1*(tstep.^0)' 1*2*(tstep.^0)' 2*3*(tstep.^1)'];

syms a b c d
[sol_a, sol_b, sol_c]=solve([tm(1,:);tm_d(1,:);tm_dd(1,:);tm(end,:)]*v*[a;b;c;]==[0;0;0;1]);

% p=v*P(i-n:i)';7 0 15 15 0 5
% p=v*[P(i-n:i-1)';-100];
% p=v*[P(i-n:i-2)';100;-100];
p=v*[-100;P(i-n+1:i-1)';-100];



% p=v*[sol_a; sol_b; sol_c; sol_d];
% s=tm*p;
% s=tm_d*p;
s=tm_dd*p;

hold on

plot(tstep,s,'*k')
hold off
%%%%