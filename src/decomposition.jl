using LinearAlgebra
export decomposition, get_entropy

function decomposition(P::Array)
    D = ndims(P)
    r = []
    P = reshape(P, 2, :)
    F = svd(P)
    push!(r, rank(Diagonal(F.S),1e-13))
    P = Diagonal(F.S)*F.Vt
    d = 2*size(F.U,2)
    for i in 2:D-1
        P = reshape(P, d, :)
        F = svd(P)
        push!(r, rank(Diagonal(F.S), 1e-13))
        P = Diagonal(F.S)*F.Vt
        d = 2*size(F.U,2)
    end
    r
end

function get_entropy(P::Array)
    D = ndims(P)
    S = []
    P = reshape(P, 2, :)
    F = svd(P)
    push!(S, -sum(((F.S).^2).*log.((F.S).^2)))
    P = Diagonal(F.S)*F.Vt
    d = 2*size(F.U,2)
    for i in 2:D-1
        P = reshape(P, d, :)
        F = svd(P)
        push!(S, -sum(F.S.*log.(F.S)))
        P = Diagonal(F.S)*F.Vt
        d = 2*size(F.U,2)
    end
    S
end
