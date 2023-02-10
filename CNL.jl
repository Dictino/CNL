using Pkg
cd(@__DIR__)
Pkg.activate(@__DIR__)
Pkg.instantiate()
using Pluto
if isfile("imagen.so")
    Pluto.run(;sysimage="imagen.so")
else
    Pluto.run()
end