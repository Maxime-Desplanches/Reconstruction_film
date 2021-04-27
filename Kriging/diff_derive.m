function [D] = diff_derive(A,B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
H=1/16*[-1 -2 -1;-2 12 -2;-1 -2 -1];
A=imfilter(A,H);
B=imfilter(B,H);
nb=size(A,2);
for i=1:size(A,1)
    A(i,1)=0;
    A(i,nb)=0;
    B(i,1)=0;
    B(i,nb)=0;
end
D=A-B;
end

