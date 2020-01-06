push!(LOAD_PATH, "mps-decomposition\\src")
push!(LOAD_PATH, "thermalvqe-master\\src")
using Main.MPS
using Main.ThermalVQE:hamiltonian, TFIM
using Yao
using YaoExtensions
using LinearAlgebra
using PyCall
using PyPlot


function ising_decomposition(nx::Int, ny::Int,β::Real, Γ::Real)
    H = hamiltonian(TFIM(nx, ny; Γ=Γ, periodic=false))
    P = get_P_tensor(H, β)
    decomposition(P)
end

function scan_beta(nx::Int, ny::Int, Γ::Real)
    H = hamiltonian(TFIM(nx, ny; Γ=Γ, periodic=false))
    println("$nx x $ny TFIM for Γ = $Γ")

    nbit = nx*ny
    rang = 0:0.1:2
    l = length(rang)
    S = zeros(l, nbit-1)

    for i in 1:l
        β = rang[i]
        P = get_P_tensor(H, β)
        s = get_entropy(P)
        println("β = $β, entanglement entropy is $s ")
        S[i,:] = s[:]
    end
    pygui(:true)
    gcf()
    plot(rang, S)
    plt.title("$nx x $ny TFIM for Γ = $Γ")
    plt.xlabel("β")
    plt.ylabel("entropy")

end

function scan_gamma(nx::Int, ny::Int, β::Real)
    println("$nx x $ny TFIM for β = $β")
    for Γ in 0:1:5
        r = ising_decomposition(nx, ny, β, Γ)
        println("Γ = $Γ, bond dimonsion is $r")
    end
end
