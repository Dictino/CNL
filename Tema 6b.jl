### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ ad471560-28e9-11eb-16b5-f7025684bc5f
using SymEngine

# ╔═╡ 3f007f22-4148-11eb-1c8b-9b6e06ccc541
using OrdinaryDiffEq

# ╔═╡ 42931bc0-4148-11eb-26b9-e5fee42c7580
using Plots

# ╔═╡ 470b0f50-4148-11eb-3e5a-810c88e0c029
using LaTeXStrings

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
md"""Ahora si usásemos la variable $x_2$ como control poríamos eliminar el término **malo**  bastaría con hacer que $x_2x_1 + x_1^3=-x_1^2 \to x_2=-x_1 - x_1^2$"""

# ╔═╡ 5159e1d0-28ec-11eb-14f4-99de193e35a1
ϕ1=-x1^2 - x1  #el  valor deseado de x2

# ╔═╡ bbd3f930-4181-11eb-1007-7fecf56f39aa
md"""Veamos que no nos hemos equivocado"""

# ╔═╡ 220f19e0-28ec-11eb-057e-0b369e37930c
expand(subs(dV1,x2=>ϕ1)) #si pudiera hacer esto realmente ya habría terminado

# ╔═╡ 883fd650-28ec-11eb-0ba3-014589882c54
md"Esto estaría guay, pero no puedo hacer que $x_2$ valga lo que yo quiera (todavía) así que defino un *error* $z_1=x_2-\phi_1$ y continuo, lo que tengo en realidad es que $x_2=\phi_1+z_1$"

# ╔═╡ c5680400-418b-11eb-32bb-ad1a2f2077ed
z1=symbols("z1")

# ╔═╡ 95969500-418d-11eb-0132-fd7f0ed24652
md"De modo que al derivar lo que realmente sale es"

# ╔═╡ 5f805e82-418b-11eb-05a4-dba1827925ec
expand(subs(dV1,x2=>ϕ1+z1)) 

# ╔═╡ 42851430-28ed-11eb-2823-a5014cb3eda0
md"""## Paso 2, controlar $x_1$ y $z_1$
Hemos reducido el problema de controlar $x_1$ a hacer $z_1=0$
¿Qué hemos ganado?
Paciencia, la idea es que nos estamos acercando a $u$
"¿Cómo me *desago* del término con $z_1$?
## Ampliemos nuestro V con $\frac{z_1^2}{2}$"""

# ╔═╡ 984ebfb0-28ed-11eb-2cee-fdbf3d493de8
V2=V1+z1^2/2

# ╔═╡ a5582280-28f4-11eb-2c05-51acb08f478b
md"Y derivemos a ver que sale"

# ╔═╡ 559eb340-28f3-11eb-3720-059fec2196ee
begin
	dz1=dx2-derivada(ϕ1,[x1],[dx1])
	dϕ1=symbols("dϕ1")
	dz1=dx2-dϕ1
end

# ╔═╡ 03ef90e0-4190-11eb-1f86-9f2bf2582d96
md"**Truco**: para no arrastrar *chorizos* podíamos crear una variable que sea la derivada"

# ╔═╡ 7d9e7930-28ed-11eb-20fd-970954cf18b7
let
	dV2=derivada(V2,[x1 z1],[dx1 dz1])
	dV2=subs(dV2,x2=>z1+ϕ1) #La sustiucion es para hacer aparecer z1
	dV2=expand(dV2)
end

# ╔═╡ da519e60-28f6-11eb-3598-c12f14f24180
md"""queremos hacer que $x_1z_1 + x_3z_1 - z_1 \dot ϕ_1 = -z_1^2$, despejamos

$x_1 + x_3 - \dot ϕ_1 = -z_1$

$x_3 = -z_1 -x_1 + \dot ϕ_1$

Y hacemos $ϕ_2=\dotϕ_1 - x_1 - z_1$

**Nota:** esto se podría despejar usando un paquete CAS más potente pero no merece la pena, hacerlo a mano ayuda a *entender* lo que estamos haciendo. El CAS lo usaremos para verificar los cálculos
"""

# ╔═╡ 2e63aed0-28f7-11eb-27a6-cbf8cb83c09c
ϕ2= -z1 -x1 + dϕ1

# ╔═╡ 0dc9ed20-4192-11eb-1039-93d7bfd5d09f
md"veamos que funciona"

# ╔═╡ ebea590e-4182-11eb-2a0b-a9161adc7f97
let
	dV2=derivada(V2,[x1 z1],[dx1 dz1])
	dV2=subs(dV2,x2=>z1+ϕ1) 
	dV2=subs(dV2,x3=>ϕ2) #lo nuevo (si pudiese hacerlo claro)
	dV2=expand(dV2)
end

# ╔═╡ a58593c0-28f7-11eb-1674-51101547d596
md"""de nuevo lo mismo **no puedo** hacer eso así que 
ahora:
$z_2=x_3-\phi_2$"""

# ╔═╡ 3094ca2e-28f8-11eb-2cc0-a1d2cebc1da6
z2=symbols("z2")

# ╔═╡ 3dbf9b10-4192-11eb-04e5-1f29ccd9da2a
let
	dV2=derivada(V2,[x1 z1],[dx1 dz1])
	dV2=subs(dV2,x2=>z1+ϕ1) 
	dV2=subs(dV2,x3=>z2+ϕ2) #lo nuevo (si pudiese hacerlo claro)
	dV2=expand(dV2)
end

# ╔═╡ fb7102a0-28f8-11eb-131d-c1792c1aaa27
md"""
## Paso 3, controlar $x_1$, $z_1$ y $z_2$
Hemos llegado a dónde queríamos por fin va a salir $u$
## Ampliemos nuestro V otra vez"""

# ╔═╡ 28437790-28f9-11eb-18f1-11364b9029ee
V3=V2+z2^2/2

# ╔═╡ 9f114912-28f9-11eb-17f1-7fa461ff3371
begin
	dϕ2=symbols("dϕ2")
	dz2=dx3-dϕ2
end

# ╔═╡ 24332a00-28fa-11eb-2a88-8b5e6f124e6d
md"Aparece $u$ ¡Yupi!"

# ╔═╡ 85a6c900-28f9-11eb-229f-f78f9cd042e4
let
	dV3=derivada(V3,[x1 z1 z2],[dx1 dz1 dz2])
	dV3=subs(dV3,x2=>z1+ϕ1)
	dV3=subs(dV3,x3=>z2+ϕ2)
	dV3=expand(dV3)
end

# ╔═╡ c9b94770-28fa-11eb-25ba-2384e1b7672a
md"""queremos hacer que $uz_2 - z_2 \dot ϕ_2 + z_2z_1 = -z_2^2$

$u - \dot ϕ_2 + z_1 = -z_2$

$u =  \dot ϕ_2 - z_1 -z_2$

Ahora vamos a ver que funciona:"""

# ╔═╡ f8482150-4183-11eb-05f2-599741c800de
let
	dV3=derivada(V3,[x1 z1 z2],[dx1 dz1 dz2])
	dV3=subs(dV3,x2=>z1+ϕ1)
	dV3=subs(dV3,x3=>z2+ϕ2)
	dV3=subs(dV3,u=>dϕ2-z1-z2)
	dV3=expand(dV3)
end

# ╔═╡ f16d5e90-4183-11eb-16dc-8191ac698ce2
md"""

## Ya tenemos ley de control, ahora hay que implementarla

"""

# ╔═╡ d4991de0-4147-11eb-067b-4f9c143fb673
function control(x,p,t)
	#completemos esto interactivamente
    x1,x2,x3=x

	ϕ1=-x1 - x1^2	
	z1=x2-ϕ1
	dϕ1=-x2 - 2*x2*x1 - x1^2 - x1^3 + 2*x1^4
	ϕ2=dϕ1 -x1 -z1
	z2=x3-ϕ2
	dϕ2=dϕ1 - 2*x2 - x3 - 2*x2*x1 - 2*x1^2 + 2*x1^4
	
	u=dϕ2-z1-z2
end

# ╔═╡ c259a7d0-28fb-11eb-3cc6-ab3e44ea196c
function derivadas(x,p,t)
	x1,x2,x3=x
	u=control(x,p,t)
	dx1=x1^2-x1^3+x2
    dx2=x3
	dx3=u
	dx=[dx1 dx2 dx3]
end

# ╔═╡ f78bff70-4192-11eb-2f14-05fbee335f4f
md"""ha llegado el momento de hacer los cálculos"""

# ╔═╡ 451fc4fe-4193-11eb-1bcd-01e96df3c0e4
ϕ1

# ╔═╡ ed7974d0-4193-11eb-3302-d1aacff66696
derivada_ϕ1=expand(derivada(ϕ1,[x1],[dx1]))

# ╔═╡ 20b19900-4193-11eb-03bb-b3fe71fdcfb4
ϕ2

# ╔═╡ 61bbe080-4194-11eb-170d-3f827d77957a
derivada_ϕ2=expand(derivada(ϕ2,[x1 dϕ1 z1],[dx1 derivada_ϕ1 dz1]))

# ╔═╡ 628ad3d2-4195-11eb-1ae2-473f3862794c
md"""vamos a ver si nos hemos equivocado"""

# ╔═╡ 4031d6c0-4196-11eb-275f-cda9a32cce36
md"## En este caso si expandimos el control mejora (es sencillo) veamoslo"

# ╔═╡ 6ff23950-4195-11eb-250b-898515a5c031
expand(control([x1 x2 x3],[],0))

# ╔═╡ fc5e8b30-4147-11eb-30d3-4bef0bac9ca6
let
	x0 = [1 0 0]
	parametros=[1 1 1] #cambiar
	tspan = (0.0,10.0)
	prob = ODEProblem(derivadas,x0,tspan,parametros);
	sol = solve(prob,Tsit5());
	
	# Salida 
	fig1=plot(sol,vars=(0,1), xaxis=L"t",yaxis=L"y(t)",label=L"y(t)")
		
	#Estados
	fig2=plot(sol, xaxis="t",yaxis=L"x(t)", label=:none)
	    
	#Señal de control
	fig3=plot(sol, vars=((t,x1,x2,x3)->(t,control([x1,x2,x3],parametros,t)),0, 1, 2, 3),xaxis=L"t",yaxis=L"u(t)",legend=false)
	l = @layout [a ; b ; c]
	plot(fig1,fig2,fig3, layout=l)
end

# ╔═╡ Cell order:
# ╠═58b659de-28fb-11eb-3ad9-53a988ea9e56
# ╟─80213b30-28e7-11eb-388d-5501de9711e7
# ╠═ad471560-28e9-11eb-16b5-f7025684bc5f
# ╟─b5d53800-28ea-11eb-1858-4360a7d9fafc
# ╟─bd1a0730-28ea-11eb-2355-c1d9285733c5
# ╠═2890d210-28ea-11eb-2f16-23af078db24a
# ╠═e0c749a0-28e9-11eb-2a9d-3b5f6566605a
# ╟─4cce2a0e-28ea-11eb-26c6-672aeb7bfdf9
# ╠═646d6000-28ea-11eb-342d-b1c902896547
# ╠═224b2df0-28eb-11eb-2a45-c9eb740f149c
# ╟─4684dc70-28eb-11eb-1fde-2575c48354e8
# ╠═5159e1d0-28ec-11eb-14f4-99de193e35a1
# ╟─bbd3f930-4181-11eb-1007-7fecf56f39aa
# ╠═220f19e0-28ec-11eb-057e-0b369e37930c
# ╟─883fd650-28ec-11eb-0ba3-014589882c54
# ╠═c5680400-418b-11eb-32bb-ad1a2f2077ed
# ╟─95969500-418d-11eb-0132-fd7f0ed24652
# ╠═5f805e82-418b-11eb-05a4-dba1827925ec
# ╟─42851430-28ed-11eb-2823-a5014cb3eda0
# ╠═984ebfb0-28ed-11eb-2cee-fdbf3d493de8
# ╟─a5582280-28f4-11eb-2c05-51acb08f478b
# ╠═559eb340-28f3-11eb-3720-059fec2196ee
# ╟─03ef90e0-4190-11eb-1f86-9f2bf2582d96
# ╠═7d9e7930-28ed-11eb-20fd-970954cf18b7
# ╟─da519e60-28f6-11eb-3598-c12f14f24180
# ╠═2e63aed0-28f7-11eb-27a6-cbf8cb83c09c
# ╟─0dc9ed20-4192-11eb-1039-93d7bfd5d09f
# ╠═ebea590e-4182-11eb-2a0b-a9161adc7f97
# ╟─a58593c0-28f7-11eb-1674-51101547d596
# ╠═3094ca2e-28f8-11eb-2cc0-a1d2cebc1da6
# ╠═3dbf9b10-4192-11eb-04e5-1f29ccd9da2a
# ╟─fb7102a0-28f8-11eb-131d-c1792c1aaa27
# ╠═28437790-28f9-11eb-18f1-11364b9029ee
# ╠═9f114912-28f9-11eb-17f1-7fa461ff3371
# ╟─24332a00-28fa-11eb-2a88-8b5e6f124e6d
# ╠═85a6c900-28f9-11eb-229f-f78f9cd042e4
# ╟─c9b94770-28fa-11eb-25ba-2384e1b7672a
# ╠═f8482150-4183-11eb-05f2-599741c800de
# ╟─f16d5e90-4183-11eb-16dc-8191ac698ce2
# ╠═d4991de0-4147-11eb-067b-4f9c143fb673
# ╠═c259a7d0-28fb-11eb-3cc6-ab3e44ea196c
# ╠═f78bff70-4192-11eb-2f14-05fbee335f4f
# ╠═451fc4fe-4193-11eb-1bcd-01e96df3c0e4
# ╠═ed7974d0-4193-11eb-3302-d1aacff66696
# ╠═20b19900-4193-11eb-03bb-b3fe71fdcfb4
# ╠═61bbe080-4194-11eb-170d-3f827d77957a
# ╟─628ad3d2-4195-11eb-1ae2-473f3862794c
# ╟─4031d6c0-4196-11eb-275f-cda9a32cce36
# ╠═6ff23950-4195-11eb-250b-898515a5c031
# ╠═3f007f22-4148-11eb-1c8b-9b6e06ccc541
# ╠═42931bc0-4148-11eb-26b9-e5fee42c7580
# ╠═470b0f50-4148-11eb-3e5a-810c88e0c029
# ╠═fc5e8b30-4147-11eb-30d3-4bef0bac9ca6
