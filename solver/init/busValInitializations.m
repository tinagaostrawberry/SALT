function [bus_VComplex_3ph,bus_VComplex] = busValInitializations(bus_Vm,bus_Va)

bus_VComplex = bus_Vm.*exp(1j*bus_Va);
[phA,phB,phC] = calc3Ph(bus_VComplex);
bus_VComplex_3ph = [phA;phB;phC];

% Sanity check that the set is balanced
vm_realMag = sqrt((2/3)*(real(phA).^2 + real(phB).^2 + real(phC).^2));
vm_posSeq = abs(phA + exp(1j*2/3*pi)*phB + (exp(1j*2/3*pi))^2*phC)/3;
assert(all(abs(vm_realMag-vm_posSeq)<1e-12))
assert(all(abs(vm_realMag-abs(phA))<1e-12))


end