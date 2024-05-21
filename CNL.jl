using Pkg
cd(@__DIR__)
Pkg.activate(@__DIR__)
if !isfile("instalado")
    prinln("Instalando por primera vez, va a tardar un rato...")
    Pkg.instantiate()
end
using Pluto
Pluto.run()