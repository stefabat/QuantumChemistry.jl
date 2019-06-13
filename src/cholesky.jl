# DEVELOPMENT FILE




function cholesky(A)
    N = size(A,1)
    L = zeros(N,N)

    for k = 1:N
        L[k,k] = sqrt(A[k,k] - sum( L[k,1:k-1].^2 ) )
        for i = k+1:N
            L[i,k] = 1/L[k,k] * (A[i,k] - sum(L[i,1:k-1] .* L[k,1:k-1]) )
        end
    end

    return L
end


function cholesky2(A)
    N = size(A,1)
    L = zeros(N,N)
    Lt = zeros(N)
    D = diag(A)

    for k = 1:N

        if k > 1
            D[k:end] -= L[k:end,k-1].^2
        end

        L[k,k] = sqrt(D[k])

        for n = k+1:N

            if k > 1
                Lt[n] += L[n,k-1] * L[k,k-1]
            end

            # compute ith entry of Jth Cholesky vector
            L[n,k] = (A[n,k] - Lt[n]) / L[k,k]

        end
    end

    return L
end



function cholesky_ERI(basis, τ::Real = 1e-6)

    # initialize array of diagonal entries
    N = length(basis)
    D = zeros(N*N)
    Lt = zeros(N*N)
    L = zeros(N*N,N*N)

    # compute diagonal
    for (i,α) in enumerate(basis)
        for (j,β) in enumerate(basis)
            eri = Integrals.repulsion(α, β, α, β)
            k = (i-1)*N + j
            @inbounds D[k] = eri
        end
    end

    # max entry is simply the first element
    # Dmax = D[1]
    # println("Dmax = ", Dmax)

    # println("Prescreening (αβ|αβ)")
    # println("Discarding elements smaller than τ²/Dmax = ",τ^2/Dmax)

    # find smallest element
    # idxmin = findfirst(x -> x <= τ^2/Dmax, D[pid])

    # in case all elements are above the threshold, set it to N*N
    # if idxmin == 0
        # idxmin = N*N
    # end

    # println("Discarded ",N*N - idxmin," elements")

    # array containing permutation indices
    p = collect(1:N*N)

    # generate Cholesky vectors
    for k = 1:N*N

        # obtain D̃ at iter k by subtracting L² elementwise from iter k-1
        if k > 1
            D[p[k:end]] -= L[k:end,k-1].^2
        end

        println("Diagonal element ",D[p[k:end]])

        Dmax,Imax = findmax(D[p[k:end]])
        J = p[k:end][Imax]
        t = findin(p,J)
        p[t] = p[k]
        p[k] = J

        if Dmax <= 1e-10
            return L,p
        end

        L[k,k] = sqrt(D[p[k]])

        # find which Gaussian product |γδ) corresponds to index Imax
        γ = fld(J-1, N) + 1; δ = mod(J-1, N) + 1
        
        # generate entries of the Jth Cholesky vector 
        for n = k+1:N*N
            
            # retrieve index corresponding to Dmax 
            i = p[n]

            # find which Gaussian product (αβ| corresponds to index i
            α = fld(i-1, N) + 1; β = mod(i-1, N) + 1
            # println("Computing entry i = ",i," → |γδ) = |",γi," ",δi,")")

            # compute 2e-integral
            eri = Integrals.repulsion(basis[α],basis[β],basis[γ],basis[δ])
            println("computing (αβ|γδ) = (",α," ",β,"|",γ," ",δ,") = ",eri)

            if k > 1
                Lt[n] += L[n,k-1] * L[k,k-1]
            end

            # compute ith entry of Jth Cholesky vector
            L[n,k] = (eri - Lt[n]) / L[k,k]

        end
    end

    return L,p

end


function order(A,p)

    B = zeros(size(A))
    for i = 1:length(p)
        B[:,p[i]] = A[:,i]
        # B[p[i],:] = A[i,:]
    end

    return B
end
