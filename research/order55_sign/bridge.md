# Exact bridge and sign moments (2026-09-11)

Convention: C(a)[i,j] = a[(j-i) mod n]. For x=2a-1 and k=sum a,
C(x)=2C(a)-J. The all-ones line is an eigenspace with eigenvalue 2k-n;
on its orthogonal complement J vanishes. Equivalently f_x(zeta^r)=2f_a(zeta^r)
at every nonzero frequency and f_x(1)=2k-n. Thus, including singular
circulants, k det C(x)=2^(n-1)(2k-n) det C(a).
For odd n and 0<k<n, conjugate frequency pairing gives
det C(a)=k product(q_r), q_r=|f_a(zeta^r)|^2 >=0.
Hence R(a)=|n-2k| product(q_r).

Complement a -> 1-a maps x -> -x. It preserves absolute determinant,
and the constant endpoints have rank one (n>1). Therefore
S_pm(55)=max_{1<=k<=27} ((55-2k)/k) D_circ(55,k).
This is a separate objective from the complementary binary coefficient 55-k.

Care with the sign convention: for x=2a-1 the SIGNED row sum is 2k-55.
The positive magnitude on the reduced weights is 55-2k.
Writing sum(x)=55-2k would reverse the chosen alphabet encoding.

Let c_s=sum a_i a_(i+s) and rho_s=sum x_i x_(i+s).
Expanding proves rho_s=55-4k+4c_s (also at s=0), so rho_s=3 mod 4.
rho_0=55 and 55+2 sum_(s=1)^27 rho_s=(2k-55)^2.

Write p_r=|f_x(zeta^r)|^2=4q_r for the 27 frequency pairs and s=2k-55.
Parseval applied first to x and then its autocorrelation gives
s^2+2 sum p_r=55^2,
s^4+2 sum p_r^2=55*(55^2+2 sum rho_s^2).
Equivalently sum q_r=k(55-k)/2 and
sum q_r^2=[55*(k^2+2 sum c_s^2)-k^4]/2.

The conversion between c and rho is an invertible affine map at a fixed
weight; its congruences encode exactly the binary correlation integrality.
There is no additional restriction solely from changing coordinates.
Integer balancing minimizes sum c_s^2 with sum c_s=k(k-1)/2.
Multiplying the rigorously bounded product of q_r by 55-2k yields all
sign-specific first, fourth and capped moment bounds.

The resultant factorization is det C(a)=k N5 N11 N55, because
z^55-1=(z-1)Phi5 Phi11 Phi55 and the nontrivial factors have even degree.
Consequently R(a)=(55-2k)N5 N11 N55 on reduced weights. N_m is the
resultant of f_a and Phi_m, a nonnegative product of conjugate pairs.

For mod-m folds r, T_m=(m sum r_i^2-k^2)/2 and
L_m=(m sum h_i^2-k^4)/2, h autocorrelation of r.
The primitive-55 first and fourth moments subtract T5,T11 and L5,L11
from the full moments. Its 20 q values give the fold-pair bound
(55-2k) N5 N11 (T55/20)^20, or the stronger moment product bound.

The exact small-weight 0/1 certificates transfer by the bridge, without
rerunning their enumerations. If the prior global binary theorem is assumed,
complement duality also gives D_circ(55,27)=27*D01(55)/28 and
S_pm(55,27)=D01(55)/28. The supplied checkout's independent full fiber-union
audit is NOT_COMPLETED; new global sign pruning does not rely on that
uncompleted independent audit or on old binary pass/fail flags.

All proof pruning must be strict below an exact certified incumbent.
Equality is retained to classify every maximizer.

