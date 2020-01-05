using Yao
using YaoExtensions
export get_P_tensor

function get_P_tensor(H::AbstractBlock, β::Real)
    H1 = Matrix(H)
    d = size(H1, 1)
    n = Int(log2(d))

    Psize = []
    for i in 1:n
        push!(Psize, 2)
    end
    P = zeros(Tuple(Psize))

    for i in 1:d
        ϕ = product_state(n, i-1)
        P[i] = exp(-β*Real(expect(H, ϕ)))
    end
    P
end
