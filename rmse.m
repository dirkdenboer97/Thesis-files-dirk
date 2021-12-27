function rmsetab = rmse(ckf_err,ukf_err,packf_err,pfckf_err,ckfapf_err)

RMSE_ckf = [    mean(sqrt((ckf_err(1,1:400).^2)+(ckf_err(2,1:400).^2)))...
                mean(sqrt((ckf_err(1,400:800).^2)+(ckf_err(2,400:800).^2))) ...
                mean(sqrt((ckf_err(1,800:1200).^2)+(ckf_err(2,800:1200).^2))) ...
                mean(sqrt((ckf_err(1,1200:1600).^2)+(ckf_err(2,1200:1600).^2))) ...
                mean(sqrt((ckf_err(1,1600:1700).^2)+(ckf_err(2,1600:1700).^2))) ...
                mean(sqrt((ckf_err(1,1700:2000).^2)+(ckf_err(2,1700:2000).^2))) ...
                mean(sqrt((ckf_err(1,2000:2400).^2)+(ckf_err(2,2000:2400).^2))) ...
                mean(sqrt((ckf_err(1,:).^2)+(ckf_err(2,:).^2)))];

RMSE_ukf = [    mean(sqrt((ukf_err(1,1:400).^2)+(ukf_err(2,1:400).^2)))...
                mean(sqrt((ukf_err(1,400:800).^2)+(ukf_err(2,400:800).^2))) ...
                mean(sqrt((ukf_err(1,800:1200).^2)+(ukf_err(2,800:1200).^2))) ...
                mean(sqrt((ukf_err(1,1200:2000).^2)+(ukf_err(2,1200:2000).^2))) ...
                mean(sqrt((ukf_err(1,1600:1700).^2)+(ukf_err(2,1600:1700).^2))) ...
                mean(sqrt((ukf_err(1,1700:2000).^2)+(ukf_err(2,1700:2000).^2))) ...
                mean(sqrt((ukf_err(1,2000:2400).^2)+(ukf_err(2,2000:2400).^2))) ...
                mean(sqrt((ukf_err(1,:).^2)+(ukf_err(2,:).^2)))];
            
RMSE_packf = [    mean(sqrt((packf_err(1,1:400).^2)+(packf_err(2,1:400).^2)))...
                mean(sqrt((packf_err(1,400:800).^2)+(packf_err(2,400:800).^2))) ...
                mean(sqrt((packf_err(1,800:1200).^2)+(packf_err(2,800:1200).^2))) ...
                mean(sqrt((packf_err(1,1200:2000).^2)+(packf_err(2,1200:2000).^2))) ...
                mean(sqrt((packf_err(1,1600:1700).^2)+(packf_err(2,1600:1700).^2))) ...
                mean(sqrt((packf_err(1,1700:2000).^2)+(packf_err(2,1700:2000).^2))) ...
                mean(sqrt((packf_err(1,2000:2400).^2)+(packf_err(2,2000:2400).^2))) ...
                mean(sqrt((packf_err(1,:).^2)+(packf_err(2,:).^2)))];
            
RMSE_pfckf = [    mean(sqrt((pfckf_err(1,1:400).^2)+(pfckf_err(2,1:400).^2)))...
                mean(sqrt((pfckf_err(1,400:800).^2)+(pfckf_err(2,400:800).^2))) ...
                mean(sqrt((pfckf_err(1,800:1200).^2)+(pfckf_err(2,800:1200).^2))) ...
                mean(sqrt((pfckf_err(1,1200:2000).^2)+(pfckf_err(2,1200:2000).^2))) ...
                mean(sqrt((pfckf_err(1,1600:1700).^2)+(pfckf_err(2,1600:1700).^2))) ...
                mean(sqrt((pfckf_err(1,1700:2000).^2)+(pfckf_err(2,1700:2000).^2))) ...
                mean(sqrt((pfckf_err(1,2000:2400).^2)+(pfckf_err(2,2000:2400).^2))) ...
                mean(sqrt((pfckf_err(1,:).^2)+(pfckf_err(2,:).^2)))];     
            
RMSE_ckfapf = [    mean(sqrt((ckfapf_err(1,1:400).^2)+(ckfapf_err(2,1:400).^2)))...
                mean(sqrt((ckfapf_err(1,400:800).^2)+(ckfapf_err(2,400:800).^2))) ...
                mean(sqrt((ckfapf_err(1,800:1200).^2)+(ckfapf_err(2,800:1200).^2))) ...
                mean(sqrt((ckfapf_err(1,1200:2000).^2)+(ckfapf_err(2,1200:2000).^2))) ...
                mean(sqrt((ckfapf_err(1,1600:1700).^2)+(ckfapf_err(2,1600:1700).^2))) ...
                mean(sqrt((ckfapf_err(1,1700:2000).^2)+(ckfapf_err(2,1700:2000).^2))) ...
                mean(sqrt((ckfapf_err(1,2000:2400).^2)+(ckfapf_err(2,2000:2400).^2))) ...
                mean(sqrt((ckfapf_err(1,:).^2)+(ckfapf_err(2,:).^2)))]; 

rmsetab = [RMSE_ckf ; RMSE_ukf ; RMSE_packf ; RMSE_pfckf ; RMSE_ckfapf];

end