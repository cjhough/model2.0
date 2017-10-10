# Parameters that work u,g,z = run([.2;.2],1.,0.,2.3,10.,1000.,0.05,1.,50000,0)



function gain(x)
    sqrt(max(x,0))
end

function dynamics(u::Array{Float64,1},a::Array{Float64,1},z::Array{Float64,1},S::Array{Float64,1},beta::Float64,alpha::Float64,gamma::Float64,tau::Float64,taua::Float64,sigma::Float64,mu::Float64,dt::Float64)

    z[1] += -dt*z[1] + sqrt(dt)*sigma*randn()
    z[2] += -dt*z[2] + sqrt(dt)*sigma*randn()

    # println(gain(S[1] - beta*u[2] + alpha*u[1] - gamma*a[1] + z[1]))
    # println(S[2] - beta*u[1] + alpha*u[2] - gamma*a[2] + z[2])

    v = u[1] + dt*(gain(S[3]+S[1] - beta*a[2]*u[2] + alpha*u[1]  + z[1]*(1.+mu*u[1])) - u[1])/tau
    u[2] += dt*(gain(S[3]+S[2] - beta*a[1]*u[1] + alpha*u[2] + z[2]*(1.+mu*u[2])) - u[2])/tau
    u[1] = v

    a[1] += dt*(1 - a[1] - gamma*a[1]*u[1])/taua
    a[2] += dt*(1 - a[2] - gamma*a[2]*u[2])/taua

    return u,a,z
end


# function dynamics(u::Array{Float64,1},a::Array{Float64,1},z::Array{Float64,1},S::Array{Float64,1},beta::Float64,alpha::Float64,gamma::Float64,tau::Float64,taua::Float64,sigma::Float64,dt::Float64)
#
#     z[1] += -dt*z[1] + sqrt(dt)*sigma*randn()
#     z[2] += -dt*z[2] + sqrt(dt)*sigma*randn()
#
#     # println(gain(S[1] - beta*u[2] + alpha*u[1] - gamma*a[1] + z[1]))
#     # println(S[2] - beta*u[1] + alpha*u[2] - gamma*a[2] + z[2])
#
#     v = u[1] + dt*(gain(S[1] - beta*a[2]*u[2]*(1+0*z[1]) + alpha*u[1]  + z[1]) - u[1])/tau
#     u[2] += dt*(gain(S[2] - beta*a[1]*u[1]*(1+0*z[2]) + alpha*u[2] + z[2]) - u[2])/tau
#     u[1] = v
#
#     a[1] += dt*(1 - a[1] - gamma*a[1]*u[1])/taua
#     a[2] += dt*(1 - a[2] - gamma*a[2]*u[2])/taua
#
#     return u,a,z
# end

# function dynamics(u::Array{Float64,1},a::Array{Float64,1},z::Array{Float64,1},S::Array{Float64,1},beta::Float64,alpha::Float64,gamma::Float64,tau::Float64,taua::Float64,sigma::Float64,dt::Float64)
#
#     z[1] += -dt*z[1] + sqrt(dt)*sigma*randn()
#     z[2] += -dt*z[2] + sqrt(dt)*sigma*randn()
#
#     # println(gain(S[1] - beta*u[2] + alpha*u[1] - gamma*a[1] + z[1]))
#     # println(S[2] - beta*u[1] + alpha*u[2] - gamma*a[2] + z[2])
#
#     v = u[1] + dt*(gain(S[1] - beta*g[2]*u[2] + alpha*u[1] - gamma*a[1] + z[1]) - u[1])/tau
#     u[2] += dt*(gain(S[2] - beta*g[1]*u[1] + alpha*u[2] - gamma*a[2] + z[2]) - u[2])/tau
#     u[1] = v
#
#     a[1] += dt*(u[1] - a[1])/taua
#     a[2] += dt*(u[2] - a[2])/taua
#
#     g[1] += dt*(1 - g[1] - phi*g[1]*u[1])
#     g[2] += dt*(1 - g[2] - phi*g[2]*u[2])
#
#     return u,a,z
# end

function run(u::Array{Float64,1},a::Array{Float64,1},z::Array{Float64,1},S::Array{Float64,1},beta::Float64,alpha::Float64,gamma::Float64,tau::Float64,taua::Float64,sigma::Float64,mu::Float64,dt::Float64,total::Int)

    uout = Array{Float64,2}(total,2)
    aout = Array{Float64,2}(total,2)
    zout = Array{Float64,2}(total,2)

    for t = 2:total
        u,a,z = dynamics(u,a,z,S,beta,alpha,gamma,tau,taua,sigma,mu,dt)
        # println(u)
        uout[t,:] = u'
        aout[t,:] = a'
        zout[t,:] = z'
    end

    # println(u[59])

    return uout,aout,zout
end

function run(S::Array{Float64,1},beta::Float64,alpha::Float64,gamma::Float64,tau::Float64,taua::Float64,sigma::Float64,mu::Float64,dt::Float64,total::Int,transient::Int)

    # println(S,beta,alpha,gamma,tau,taua,sigma)

    uout = Array{Float64,2}(total,2)
    aout = Array{Float64,2}(total,2)
    zout = Array{Float64,2}(total,2)

    u = rand(2)
    a = ones(2)
    z = zeros(2)

    for t = 1:transient
        u,a,z = dynamics(u,a,z,[0;0;S[3]],beta,alpha,gamma,tau,taua,sigma,mu,dt)
    end

    # println(u)

    uout[1,:] = u'
    aout[1,:] = a'
    zout[1,:] = z'

    # println(aout[1,:])

    for t = 2:total
        u,a,z = dynamics(u,a,z,S,beta,alpha,gamma,tau,taua,sigma,mu,dt)
        # println(u)
        uout[t,:] = u'
        aout[t,:] = a'
        zout[t,:] = z'

    end
    # uout,aout,zout = run(u,a,z,S,beta,alpha,gamma,tau,taua,sigma,dt,total)

    # println(u[59])

    return uout,aout,zout
end

function domtime(u::Array{Float64,2},dt,thresh)

    dom1 = Array{Float64,1}(0);
    dom2 = Array{Float64,1}(0);
    on1 = Array{Float64,1}(0);
    on2 = Array{Float64,1}(0);

    # println(length(on1))

    difp = (u[1,1]-u[1,2])/(u[1,1]+u[1,2])
    for t = 2:size(u)[1]
        dif = (u[t,1]-u[t,2])/(u[t,1]+u[t,2])
        if dif >= thresh && difp < thresh
            push!(on1,t*dt)
        elseif dif < thresh && difp >= thresh
            if length(on1) > 0
                push!(dom1,t*dt - on1[end])
            end
        elseif dif <= -thresh && difp > -thresh
            push!(on2,t*dt)
        elseif dif > -thresh && difp <= -thresh
            if length(on2) > 0
                push!(dom2,t*dt - on2[end])
            end
        end
        difp = dif
    end
    if length(on1) < 1
        on1 = [100000.]
    end
    if length(on2) < 1
        on2 = [100000.]
    end

    return on1,dom1,on2,dom2

end

function cv(u,dt,thresh)
    out = domtime(u,dt,thresh)
    return std(out[2])/mean(out[2]), std(out[4])/mean(out[4])
end

function firstepoch(S::Array{Float64,1},beta::Float64,alpha::Float64,gamma::Float64,tau::Float64,taua::Float64,sigma::Float64,mu::Float64,dt::Float64,thresh::Float64)

    transient = round(Int,taua)
    total = round(Int,taua)

    C1 = 0
    C2 = 0
    nexps = 800
    for i = 1:nexps
        u,a,z = run(S,beta,alpha,gamma,tau,taua,sigma,mu,dt,total,transient)
        out = domtime(u,dt,thresh)
        # println(out[1][1]-out[3][1])
        if length(out[1]) > 0 && length(out[3]) > 0
            if out[3][1] > out[1][1]
                C1 += 1
            elseif out[1][1] > out[3][1]
                C2 += 1
            end
        end
    end

    if C1 == 0
        C1 = .01
    end
    if C2 == 0
        C2 = 0.01
    end

    return C1/(C1+C2), C2/(C1+C2)

end

function meanstd(S::Array{Float64,1},beta::Float64,alpha::Float64,gamma::Float64,tau::Float64,taua::Float64,sigma::Float64,mu::Float64,dt::Float64,thresh::Float64,total::Int,transient::Int)

    u,a,z = run(S,beta,alpha,gamma,tau,taua,sigma,mu,dt,total,transient)
    out = domtime(u,dt,thresh)

    return mean(out[2]),std(out[2]),mean(out[4]),std(out[4])
end

function outputs(p,dt,thresh,total,transient)

    for i = 1:5

        C = i-1

        beta = p[1]
        alpha = 0.
        gamma = p[2]
        tau = 10.
        taua = p[3]
        sigma = p[4]
        mu = p[8]
        S = [p[5] + p[6]*C;p[5] - p[6]*C;p[7]]

        # println(p)

        P1,P2 = firstepoch(S,beta,alpha,gamma,tau,taua,sigma,mu,dt,thresh)
        TD1,STD1,TD2,STD2 = meanstd(S,beta,alpha,gamma,tau,taua,sigma,mu,dt,thresh,total,transient)
        println(P1," ",TD1," ",STD1," ",TD2," ",STD2)

    end
end

function chisq(p,dt,thresh,total,transient)

    # DataP1 = [.5;.96;.99;1.;1.]
    # DataTD1 = [1800.;1900.;2200.;2200.;3500.]
    # DataSTD1 = [560.;590.;850.;820.;2400.]
    # DataTD2 = [1700.;1700.;1700.;1700.;1700.]
    # DataSTD2 = [400.;500.;600.;400.;500.]

    DataP1 = [.5;.95;.99;1.;1.]
    DataTD1 = [1800;1900;2200;2700;3200]
    DataSTD1 = DataTD1*.6
    DataTD2 = [1700;1700;1700;1700;1700]
    DataSTD2 = DataTD2*.6

    per = .25
    vP1 = (DataP1*.1).^2
    vTD1 = (DataTD1*per).^2
    vSTD1 = (DataSTD1*per).^2
    vTD2 = (DataTD2*per).^2
    vSTD2 = (DataTD2*per).^2

    # vP1 = [.1^2;.1^2;.1^2;.1^2;.1^2]
    # vTD1 = [400.^2;400.^2;400.^2;400.^2;400.^2]
    # vSTD1 = [200.^2;200.^2;200.^2;200.^2;200.^2]
    # vTD2 = [400.^2;400.^2;400.^2;400.^2;400.^2]
    # vSTD2 = [200.^2;200.^2;200.^2;200.^2;200.^2]

    chisq = Array{Float64,1}(5)

    for i = 1:5

        C = i-1
        beta = p[1]
        alpha = 0.
        gamma = p[2]
        tau = 10.
        taua = p[3]
        sigma = p[4]
        mu = p[8]
        S = [p[5] + 0*p[6]*C;p[5] - p[6]*C;p[7]]
        # println(p)

        P1,P2 = firstepoch(S,beta,alpha,gamma,tau,taua,sigma,mu,dt,thresh)
        TD1,STD1,TD2,STD2 = meanstd(S,beta,alpha,gamma,tau,taua,sigma,mu,dt,thresh,total,transient)

        if isnan(TD1)
            TD1 = 10000
        end
        if isnan(TD2)
            TD2 = 10000
        end
        if isnan(STD1)
            STD1 = 0
        end
        if isnan(STD2)
            STD2 = 0
        end

        chisq[i] = (P1-DataP1[i])^2/vP1[i] + (TD1-DataTD1[i])^2/vTD1[i] +(STD1-DataSTD1[i])^2/vSTD1[i] + (TD2-DataTD2[i])^2/vTD2[i] +(STD2-DataSTD2[i])^2/vSTD2[i]

    end
    return chisq
end

function mcmc(p,nsteps,dt=1.)

    # p = [1.;2.3;1000.;.05;.1;.1;.01]
    plength = length(p)
    pml = p
    pk = p

    thresh = .3
    total = 5000000
    transient = 3000

    chisqp = chisq(p,dt,thresh,total,transient)
    chisqk = sum(chisqp)
    chisqml = chisqk

    for step = 1:nsteps
        p = max(pk + .03*max(pk,.00001).*randn(plength),0.)
        chisqp = chisq(p,dt,thresh,total,transient)
        chisqt = sum(chisqp)
        if exp(-.5*(chisqt-chisqk)) > rand()
            chisqk = chisqt
            pk = p
            if chisqk < chisqml
                chisqml = chisqk
                pml = pk
            end
        end
        println(chisqk)
        println(pk)
    end
    return pml,chisqml
end
