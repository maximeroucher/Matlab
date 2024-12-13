function A=matdif(d)
unit=ones(d-1,1);
%
B1 = diag(unit,-1);
B2 = diag(unit,1);
%
A = 2*eye(d) - B1 - B2;
end
