function Kc = computeKc(G,A,Bw,Xi)
    Kc = G * inv(A) * Bw * Xi;
end