function [phA,phB,phC] = calc3Ph(phA)
phB = phA * (-0.5-1j*0.5*sqrt(3)); % -120 deg
phC = phA * (-0.5+1j*0.5*sqrt(3)); % 120 deg
end