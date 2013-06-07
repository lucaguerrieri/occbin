if exist('simC')==1
    c_s = simC';
    b_s = simB';
end

    y_p=exp(y_difference);
    y_l=exp(y_uncdifference);
    
    b_p=b_difference+M;
    blag_p=lg(b_difference,1)+M;
    blag_p(1)=b_ss+M*bstock(ibstock);
    b_l=b_uncdifference+M;
    blag_l=lg(b_uncdifference,1)+M;
    blag_l(1)=b_ss+M*bstock(ibstock);
    
    c_p = y_p + b_p - R*blag_p ;
    c_l = y_l + b_l - R*blag_l ;

