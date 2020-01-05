push!(LOAD_PATH, "mps-decomposition\\src")
push!(LOAD_PATH, "thermalvqe-master\\src")
using Main.MPS
using Main.ThermalVQE:hamiltonian, TFIM
using Yao
using YaoExtensions

function ising_decomposition(nx::Int, ny::Int,β::Real, Γ::Real)
    H = hamiltonian(TFIM(nx, ny; Γ=Γ, periodic=false))
    P = get_P_tensor(H, β)
    decomposition(P)
end

function scan_beta(nx::Int, ny::Int, Γ::Real)
    H = hamiltonian(TFIM(nx, ny; Γ=Γ, periodic=false))
    println("$nx x $ny TFIM for Γ = $Γ")
    for β in 0:1:10
        P = get_P_tensor(H, β)
        S = get_entropy(P)
        println("β = $β, entanglement entropy is $S")
    end
end

function scan_gamma(nx::Int, ny::Int, β::Real)
    println("$nx x $ny TFIM for β = $β")
    for Γ in 0:1:5
        r = ising_decomposition(nx, ny, β, Γ)
        println("Γ = $Γ, bond dimonsion is $r")
    end
end
