using Yao
using YaoExtensions
export get_P_tensor

function get_P_tensor(H::AbstractBlock, β::Real)
    H1 = Matrix(H)
    d = size(H1, 1)
    n = Int(log2(d))
    F = eigen(H1)
    val = sort(Real.(F.values))

    Psize = []
    for i in 1:n
        push!(Psize, 2)
    end
    P = zeros(Tuple(Psize))

    for i in 1:d
        P[i] = exp(-β*val[i])
    end
    P./sum(P)
end
