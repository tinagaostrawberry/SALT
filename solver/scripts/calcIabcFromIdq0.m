function [iabc,iabc_genrou,invT] = calcIabcFromIdq0(THETA,id,iq,genrouToBus_I)
invT = getInvParkMatrices(THETA);
iabc_genrou = invT*[id;iq;0];
iabc = iabc_genrou * genrouToBus_I;
end