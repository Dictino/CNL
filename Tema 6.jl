### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ ad471560-28e9-11eb-16b5-f7025684bc5f
using SymEngine

# ╔═╡ 58b659de-28fb-11eb-3ad9-53a988ea9e56
md"""# OJOOOOOOOOOOOOOOOOOO falta poner las ks y simular"""

# ╔═╡ 80213b30-28e7-11eb-388d-5501de9711e7
md"""# Vamos a empezar por el final
## Ejemplo final de la teoría, estabilizar:

$$\dot x_1=x_1^2-x_1^3+x_2$$
$$\dot x_2=x_3$$
$$\dot x_3 =u$$

# Plan de acción
* Con $u$ controlo $x_3$
* Con $x_3$ controlo $x_2$
* Con $x_2$ controlo $x_1$ 

Pero el plan se ejecuta al revés ☺

"""

# ╔═╡ b5d53800-28ea-11eb-1858-4360a7d9fafc
md"De nuevo usamos la función de tema 4"

# ╔═╡ bd1a0730-28ea-11eb-2355-c1d9285733c5
begin
	function derivada(funcion,variables,derivadas)
		d=0
		for i=1:length(variables)
			d=d+diff(funcion,variables[i])*derivadas[i]
		end
		return d
	end
	
	function derivada(f,x,dx,n) #multiple distpach
		d=f
		for i=1:n
			d=derivada(d,x,dx)
		end
		return d
	end
end

# ╔═╡ 2890d210-28ea-11eb-2f16-23af078db24a
x1,x2,x3,u=symbols("x1,x2,x3,u")

# ╔═╡ e0c749a0-28e9-11eb-2a9d-3b5f6566605a
begin
	dx1=x1^2-x1^3+x2
    dx2=x3
	dx3=u
end

# ╔═╡ 4cce2a0e-28ea-11eb-26c6-672aeb7bfdf9
md"## Paso 1, controlar $x_1$"

# ╔═╡ 646d6000-28ea-11eb-342d-b1c902896547
V1=x1^2/2

# ╔═╡ 224b2df0-28eb-11eb-2a45-c9eb740f149c
dV1=expand(derivada(V1,[x1],[dx1]))

# ╔═╡ 4684dc70-28eb-11eb-1fde-2575c48354e8
md"""Ahora si usásemos la variable $x_2$ como control poríamos eliminar el término **malo** vamos a jugar un poco"""

# ╔═╡ 5159e1d0-28ec-11eb-14f4-99de193e35a1
x2_deseado=-x1^2 - x1  #el ϕ_1 de la teoría

# ╔═╡ 220f19e0-28ec-11eb-057e-0b369e37930c
expand(subs(dV1,x2=>x2_deseado))

# ╔═╡ 883fd650-28ec-11eb-0ba3-014589882c54
md"Esto estaría guay, pero no puedo hacer que $x_2$ valga lo que yo quiera (todavía) así que defino un *error* y continuo"

# ╔═╡ 9a568450-28f2-11eb-250f-4b7af09c70d7
md"""## Ahora definimos 
$z_1=x_2-x_{2deseado}=x_2-ϕ_1$ 
y... ¡OJO NO, no lo hagáis *en el código*!  
### *Truco*, NO hagamos las sustituciones
Si le decimos al ordenador lo que es z lo va a sustituir inmetiatamente, no queremos eso queremos usar z como una variable, por eso lo definimos como tal y le decimos al ordenador cuál es su derivada """

# ╔═╡ 0a93d7f2-28f2-11eb-1037-0be3114dcc6b
begin
	ϕ1,dϕ1,z1=symbols("ϕ1,dϕ1,z1")
	#z1=x2-x2_deseado #OJO NO HACER ESTO 
	#Si lo hacemos el ordenador sustituirá todas las z1 y no queremos que lo haga
end

# ╔═╡ 97ad19c0-28f3-11eb-2b30-c96dbc9d8337
md"Y defino la derivada que es lo único que necesito para continuar"

# ╔═╡ 559eb340-28f3-11eb-3720-059fec2196ee
dz1=dx2+dϕ1

# ╔═╡ 42851430-28ed-11eb-2823-a5014cb3eda0
md"""## Paso 2, controlar $x_1$ y $z_1$
Hemos reducido el problema de controlar $x_1$ a hacer $z_1=0$
¿Qué hemos ganado?
Paciencia, la idea es que nos estamos acercando a $u$
## Ampliemos nuestro V"""

# ╔═╡ 984ebfb0-28ed-11eb-2cee-fdbf3d493de8
V2=V1+z1^2/2

# ╔═╡ a5582280-28f4-11eb-2c05-51acb08f478b
md"Y derivemos a ver que sale"

# ╔═╡ 7d9e7930-28ed-11eb-20fd-970954cf18b7
begin
	dV2=derivada(V2,[x1 z1],[dx1 dz1])
	dV2=subs(dV2,x2=>z1+x2_deseado) #La sustiucion es para hacer aparecer z1
	dV2=expand(dV2)
end

# ╔═╡ da519e60-28f6-11eb-3598-c12f14f24180
md"""queremos hacer que $x_1z_1 + x_3z_1 + z_1 \dot ϕ_1 = -z_1^2$, despejamos

$x_1 + x_3 + \dot ϕ_1 = -z_1$

$x_3 = -z_1 -x_1 - \dot ϕ_1$

**Nota:** esto se podría despejar usando un paquete CAS más potente pero no merece la pena, hacerlo a mano ayuda a *entender* lo que estamos haciendo. El CAS lo usaremos para verificar los cálculos
"""

# ╔═╡ 2e63aed0-28f7-11eb-27a6-cbf8cb83c09c
x3_deseado= -z1 -x1 - dϕ1

# ╔═╡ a58593c0-28f7-11eb-1674-51101547d596
md"""de nuevo lo mismo, ahora:
$z_2=x_3-x_{3deseado}=x_3-\phi_2$
y lo dejo *sin definir* para que no lo sustituya demasiado pronto como antes.
"""

# ╔═╡ 3094ca2e-28f8-11eb-2cc0-a1d2cebc1da6
ϕ2,dϕ2,z2=symbols("ϕ2,dϕ2,z2")

# ╔═╡ 9f114912-28f9-11eb-17f1-7fa461ff3371
dz2=dx3-dϕ2

# ╔═╡ 24332a00-28fa-11eb-2a88-8b5e6f124e6d
md"Aparece $u$ ¡Yupi!"

# ╔═╡ fb7102a0-28f8-11eb-131d-c1792c1aaa27
md"""
## Paso 3, controlar $x_1$, $z_1$ y $z_2$
Hemos llegado a dónde queríamos por fin va a salir $u$
## Ampliemos nuestro V otra vez"""

# ╔═╡ 28437790-28f9-11eb-18f1-11364b9029ee
V3=V2+z2^2/2

# ╔═╡ 85a6c900-28f9-11eb-229f-f78f9cd042e4
begin
	dV3=derivada(V3,[x1 z1 z2],[dx1 dz1 dz2])
	dV3=subs(dV3,x2=>z1+x2_deseado)
	dV3=subs(dV3,x3=>z2+x3_deseado) #esta es nueva
	dV3=expand(dV3)
end

# ╔═╡ c9b94770-28fa-11eb-25ba-2384e1b7672a
md"""queremos hacer que $uz_2 - z_2 \dot ϕ_2 + z_2z_1 = -z_2^2$

$u - \dot ϕ_2 + z_1 = -z_2$

$u =  \dot ϕ_2 - z_1 -z_2$

## Ya tenemos ley de control, ahora hay que implementarla

"""

# ╔═╡ c259a7d0-28fb-11eb-3cc6-ab3e44ea196c


# ╔═╡ Cell order:
# ╠═58b659de-28fb-11eb-3ad9-53a988ea9e56
# ╟─80213b30-28e7-11eb-388d-5501de9711e7
# ╠═ad471560-28e9-11eb-16b5-f7025684bc5f
# ╟─b5d53800-28ea-11eb-1858-4360a7d9fafc
# ╟─bd1a0730-28ea-11eb-2355-c1d9285733c5
# ╠═2890d210-28ea-11eb-2f16-23af078db24a
# ╠═e0c749a0-28e9-11eb-2a9d-3b5f6566605a
# ╠═4cce2a0e-28ea-11eb-26c6-672aeb7bfdf9
# ╠═646d6000-28ea-11eb-342d-b1c902896547
# ╠═224b2df0-28eb-11eb-2a45-c9eb740f149c
# ╟─4684dc70-28eb-11eb-1fde-2575c48354e8
# ╠═5159e1d0-28ec-11eb-14f4-99de193e35a1
# ╠═220f19e0-28ec-11eb-057e-0b369e37930c
# ╟─883fd650-28ec-11eb-0ba3-014589882c54
# ╟─9a568450-28f2-11eb-250f-4b7af09c70d7
# ╠═0a93d7f2-28f2-11eb-1037-0be3114dcc6b
# ╠═97ad19c0-28f3-11eb-2b30-c96dbc9d8337
# ╠═559eb340-28f3-11eb-3720-059fec2196ee
# ╟─42851430-28ed-11eb-2823-a5014cb3eda0
# ╠═984ebfb0-28ed-11eb-2cee-fdbf3d493de8
# ╟─a5582280-28f4-11eb-2c05-51acb08f478b
# ╠═7d9e7930-28ed-11eb-20fd-970954cf18b7
# ╟─da519e60-28f6-11eb-3598-c12f14f24180
# ╠═2e63aed0-28f7-11eb-27a6-cbf8cb83c09c
# ╟─a58593c0-28f7-11eb-1674-51101547d596
# ╠═3094ca2e-28f8-11eb-2cc0-a1d2cebc1da6
# ╠═9f114912-28f9-11eb-17f1-7fa461ff3371
# ╟─24332a00-28fa-11eb-2a88-8b5e6f124e6d
# ╟─fb7102a0-28f8-11eb-131d-c1792c1aaa27
# ╠═28437790-28f9-11eb-18f1-11364b9029ee
# ╠═85a6c900-28f9-11eb-229f-f78f9cd042e4
# ╠═c9b94770-28fa-11eb-25ba-2384e1b7672a
# ╠═c259a7d0-28fb-11eb-3cc6-ab3e44ea196c
