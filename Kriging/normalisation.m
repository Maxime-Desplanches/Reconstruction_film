function [A] = normalisation(A,laquelle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
switch laquelle
    case 'inf'
        A=A/max(A(:));
    case 'normal'
        A=A/mean(A(:));
        
end

