function [mu,rhok] = isgauss(inno,Gaussvec)
        
    innon = inno./std(inno')';
    rhok = [kstest2(innon(1,:));kstest2(innon(2,:));kstest2(innon(3,:));kstest2(innon(4,:))];
    mu = mean(Gaussvec,2);

    end

    
    
    
    
    
    
    