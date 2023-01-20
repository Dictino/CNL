using Pkg
Pkg.activate(@__DIR__)
cd(@__DIR__)
using Pluto
if isfile("imagen.so")
    Pluto.run(;sysimage="imagen.so")
else
    Pluto.run()
end