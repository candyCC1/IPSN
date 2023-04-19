%注意，X的元素数目为L+snap,前面L个留给tau_el.
function c = CumConjNQ(X,m,n,p,q,HasItem3)
    N_x = size(X,2);
    %计算cum(X(:,m),X(:,n),X(:,p),X(:,q))
    %注意，不省略第三项！！！因为有很多类型的信号不满足零均值等条件。
    s_temp1=0;
    s_temp2=0;
    s_temp3=0;
    s_temp4=0;
    s_temp5=0;
    s_temp6=0;
    s_temp7=0;
    for c=1:N_x
        % 第一项
        s_temp1=s_temp1+X(m,c)*conj(X(n,c))*X(p,c)*conj(X(q,c));
        % 第二项
        s_temp2=s_temp2+X(m,c)*conj(X(n,c));
        s_temp3=s_temp3+conj(X(p,c))*X(q,c);
        % 第三项
        s_temp4=s_temp4+X(m,c)*X(q,c);
        s_temp5=s_temp5+conj(X(n,c))*conj(X(p,c));
        % 第四项
        s_temp6=s_temp6+X(m,c)*conj(X(p,c));
        s_temp7=s_temp7+conj(X(n,c))*X(q,c);
    end
    s_temp1=s_temp1/N_x;
    s_temp2=s_temp2/N_x;
    s_temp3=s_temp3/N_x;
    s_temp4=s_temp4/N_x;
    s_temp5=s_temp5/N_x;
    s_temp6=s_temp6/N_x;
    s_temp7=s_temp7/N_x;
        % 总项
    if HasItem3
        c=s_temp1-s_temp2*s_temp3-s_temp4*s_temp5-s_temp6*s_temp7;
    else
        c=s_temp1-s_temp2*s_temp3-s_temp6*s_temp7;
    end
end