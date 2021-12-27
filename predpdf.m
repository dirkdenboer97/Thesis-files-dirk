function pdf = predpdf(hsf,stanR,inno_win,win,binss,binsize)


        %state1
        %hs = hsf*min(stanR(1),iqr(inno_win(1,:))/1.34)*win^(-0.2);
        hs = hsf*max(1e-10,stanR(1))*win^(-0.2);
        xi=ones(win+1,1)*binss; 
        data2i=inno_win(1,:)'*ones(1,binsize); 
        ud=(xi-data2i)/hs; 
        preddist=(.75/hs)*(1-ud.^2).*(abs(ud)<=1);
        pdf = sum(preddist,1)./sum(sum(preddist,1));


end