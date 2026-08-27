q = 4*u^2-6*u+3;
W1 = r*(q*r-(2*u-1));
W2 = q*r^2-(4*u^2-4*u+2)*r+(2*u-1);
W4 = (16*u^4-40*u^3+40*u^2-18*u+3)*r^2+(-16*u^3+28*u^2-18*u+4)*r+(4*u^2-4*u+1);
pr = W1*W2;
print("W1*W2 = ", pr);
print("factor(W1*W2) over Q(r): ", factor(pr));
\\ is it r * (something)^2 ?
print("W1*W2 / r  square test: ", issquare(pr/r));
print("factor(W1*W4): ", factor(W1*W4));
quit
