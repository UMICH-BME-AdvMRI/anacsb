function [msig,m] = sprofile(rf_signal,gss,tau_rf,T1,T2,pos,df,dt,gammabar)


for jj = 1: length(df)
    M =[0;0;1];	                    % Starting magnetization.
    [A,B] = freeprecess(tau_rf/2,T1,T2,df(jj));
    for ii=1:length(pos)
        for k = 1:length(rf_signal)
	        M = A*M+B;
	        gssrot = 2*pi*gammabar*(pos(ii))*gss(k)*dt/2;
            Rz = [cos(gssrot) -sin(gssrot) 0;sin(gssrot) cos(gssrot) 0; 0 0 1];
	        M = Rz*M;
        
            M = throt(abs(rf_signal(k)),angle(rf_signal(k))) * M;	% RF Rotation.
        
	        M = A*M+B;
	        gssrot = 2*pi*gammabar*(pos(ii))*gss(k)*dt/2;
            Rz = [cos(gssrot) -sin(gssrot) 0;sin(gssrot) cos(gssrot) 0; 0 0 1];
    
	        M = Rz*M;
        end
    
        m(:,ii,jj) = M;
        msig(jj,ii) = M(1)+i*M(2);
    end
end
