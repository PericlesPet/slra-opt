m_t  = (ell+1) * (m_in + p_out);
n_t  = (T-ell);
np = length(p);
phi  = eye(m_t);
pext = [0; p];
tts  = s2s(s,np);
vec_tts = tts(:); NP = 1:np;
bfs = vec_tts(:, ones(1, np)) == NP(ones(m_t * n_t, 1), :);
s0_2   = zeros(m_t, n_t);

myHank = p(tts);   % EQ. to ( phi * (s0_2 + pext(tts + 1)) )
%% TEST HOW u,y COMPARE TO HANKEL MATRIX (Placement?)
first3 = myHank(1:12,1); second3 = myHank(1:12,2);
u1_31 = u(1:3,1); u2_31 = u(1:3,2); y1_31 = y(1:3,1); y2_31 = y(1:3,2);
u1_32 = u(2:4,1); u2_32 = u(2:4,2); y1_32 = y(2:4,1); y2_32 = y(2:4,2);
[[first3  - [u1_31;u2_31;y1_31;y2_31]] [second3 - [u1_32;u2_32;y1_32;y2_32]] ] 
% Check R*Hank(X) CONDITIONS 
[norm(RmyOpt * myHank) norm(info.Rh * myHank) norm(info.Rh * myhHank) ... 
    norm(info3.Rh * myHank) norm(info3.Rh) norm(myHank)]

%% This uses the 720x1 p (input & output), whereas the next one only uses the output (y) part
R_test  = Rini;
th_test = R2th(R_test, phi * (s0_2 + pext(tts + 1)), psi, opt.R0);
myObj = @(th)Mslra_ext(th2R(th),tts,p,[],bfs,phi,s0) ... 
    +opt.g*norm(C(th),'fro')^2;
[M_val, next_p] = Mslra_ext(R_test, tts, p, [], bfs, phi, s0_2); 

%% Test G(R)*p = h(R)
% R_test = info3.Rh;
% np2 = length(p_new);
% g = reshape(R_test * reshape(bfs, m_t, n_t * np2), size(R_test, 1) * n_t, np2);
% h = g * p_new + vec(R_test * s0_2);
% g*(p_n(361:end)-ph(361:end)) - h
% norm(g*(p(361:end)-ph(361:end)) - h)
%% Test ALM Function
Y = eye(size(R_test * myHank));
ALMval = 1/2*norm(ph - p)^2 + norm(R_test * myHank + Y)^2 - 1/2*norm(Y)^2;
dFdx = ph - p; 




%% PH IS [U1 ... U_m ; Yh1 ... Yh_p] ,  [U1 ... U_m] DOESNT CHANGE CUZ WE KNOW IT
% THEREFORE MSLRA ONLY GIVES YH, NOT THE WHOLE THING
% CHECK THAT ALM ACTUALLY DECREASES WITH SLRA ITERATIONS -> IT SHOULD MAKE
% SENSE
% Also: find ph_ini, Rini

Ms              = [];
Phs             = [];
ALMs            = [];
ALM_component1s = [];
ALM_component2s = [];

mslra_handle = @(R) Mslra_ext(R, wtfdata.tts, wtfdata.p, [], ... 
    wtfdata.bfs, wtfdata.phi, wtfdata.s0);

for i = 1:100
    R_current = info.RhK(:,:,i);
    [M, ph_n] = mslra_handle(R_current);
    Ms = [Ms M];
    ph_actual = [p(isinf(s.w)); ph_n];
    Phs = [Phs ph_actual];
    hankel_current = ph_actual(tts);
    
    ALM = 1/2*norm(ph_actual - p)^2 + ... 
        1/2*norm(R_current * hankel_current + Y)^2 - 1/2*norm(Y)^2;
    ALM_component1 = 1/2*norm(ph_actual - p)^2 ;
    ALM_component2 = 1/2*norm(R_current * hankel_current + Y)^2 - 1/2*norm(Y)^2;

    ALMs = [ALMs ALM]; 
    ALM_component1s = [ALM_component1s ALM_component1]; 
    ALM_component2s = [ALM_component2s ALM_component2]; 
end
plot(1:length(ALMs), ALMs, '--', ... 
    1:length(ALMs), ALM_component1s, '-.', ... 
    1:length(ALMs) , ALM_component2s, ... 
    1:length(ALMs), Ms)
legend('ALM', 'ALM Comp1', 'ALM Comp2', 'M')
%% NEW HANKEL MATRIX FOR NEW Ph - ALSO CHECK LOW RANK THRU IMAGE OR SVD
myhHank = ph(tts);
[uh sh vd] = svd(myhHank);
image(sh*10)
sh(1:12,1:12)


%% CE(X) GRAD TEST
cex = @(X)reshape(... 
        reshape(X(np+1:end), size(Rini)) * X(tts), ... 
        size(Rini, 1) * (T-ell), 1 ...
        );

    
ncex = @(X) norm( ...
    reshape(...             % Reshaping to a vector actually improves speed by ~2x
    reshape(X(np+1:end), size(Rini)) * X(tts), ... 
    size(Rini, 1) * (T-ell), 1 ...
    ));

gstep = 1/1e7;
grad = FDGradient(ncex, X, gstep);
gradFcn = @(X) FDGradient(ncex, X, gstep);



%% ALM test
dce = problem.dce;
df = problem.df;
%%

pinv(dce(x0))*df(x0)
%%
size(pinv(dce(x0)))
size(dce(x0))
%%
size(df(x0))


