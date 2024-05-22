using Pkg
cd(@__DIR__)
Pkg.activate(@__DIR__)
if !isfile("instalado")
    println("Instalando por primera vez, va a tardar un rato...")
    Pkg.instantiate()
    touch("instalado")
end
using Pluto
Pluto.run()