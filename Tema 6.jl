### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 42931bc0-4148-11eb-26b9-e5fee42c7580
begin
    import Pkg
	  Pkg.activate(@__DIR__)
	using OrdinaryDiffEq, Plots, LaTeXStrings, PlutoUI, SymEngine
end

# ╔═╡ 80213b30-28e7-11eb-388d-5501de9711e7
md"""# Tema 6 Backstepping

## Ejemplo final de la teoría, estabilizar:

$$\dot x_1=x_1^2-x_1^3+x_2$$
$$\dot x_2=x_3$$
$$\dot x_3 =u$$

# Plan de acción
* Con $u$ controlo $x_3$
* Con $x_3$ controlo $x_2$
* Con $x_2$ controlo $x_1$ 

Pero el plan se ejecuta al revés, por eso se llama back-stepping ☺

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
x1,x2,x3,u=symbols("x1,x2,x3,u");

# ╔═╡ e0c749a0-28e9-11eb-2a9d-3b5f6566605a
begin
	dx1=x1^2-x1^3+x2
    dx2=x3
	dx3=u
end;

# ╔═╡ 4cce2a0e-28ea-11eb-26c6-672aeb7bfdf9
md"## Paso 1, controlar $x_1$"

# ╔═╡ 646d6000-28ea-11eb-342d-b1c902896547
V1=x1^2/2

# ╔═╡ 224b2df0-28eb-11eb-2a45-c9eb740f149c
dV1=expand(derivada(V1,[x1],[dx1]))

# ╔═╡ 2d6d8fe0-c925-11eb-3a32-2fd273a45a10
md"Vamos a *tantear* un poco para entender cómo funciona esto, para ello vamos a jugar con $x_2$ a ver si conseguimos que V disminuya como nosotros queremos"

# ╔═╡ 4684dc70-28eb-11eb-1fde-2575c48354e8
md"""De forma directa $x_2$ como control poríamos eliminar el término **malo** y dejar uno cuadrático bueno bastaría con hacer que $x_2x_1 + x_1^3=-k_1x_1^2 \to x_2=-k_1x_1 - x_1^2$"""

# ╔═╡ d589c21e-5fd5-11eb-3c70-535e951c59cf
k1=symbols("k1");

# ╔═╡ 2c671c60-c925-11eb-303d-b398113bdb9b
expand(subs(dV1,x2=>-k1*x1-x1^2))

# ╔═╡ 01b23ad0-c926-11eb-2dc3-f143c1d2b6f4
md"Vamos a guardar lo que hemos encontrado en una funciónde control $\phi_1$"

# ╔═╡ 5159e1d0-28ec-11eb-14f4-99de193e35a1
ϕ1=-k1*x1-x1^2  #el  valor deseado de x2

# ╔═╡ bbd3f930-4181-11eb-1007-7fecf56f39aa
md"""Veamos que no nos hemos equivocado"""

# ╔═╡ 220f19e0-28ec-11eb-057e-0b369e37930c
expand(subs(dV1,x2=>ϕ1)) #si pudiera hacer esto realmente ya habría terminado

# ╔═╡ 883fd650-28ec-11eb-0ba3-014589882c54
md"Esto estaría guay, pero no puedo hacer que $x_2$ valga lo que yo quiera (todavía) así que defino un *error* $z_1=x_2-\phi_1$ y continuo, lo que tengo en realidad es que $x_2=\phi_1+z_1$"

# ╔═╡ c5680400-418b-11eb-32bb-ad1a2f2077ed
z1=symbols("z1");

# ╔═╡ 95969500-418d-11eb-0132-fd7f0ed24652
md"De modo que al derivar lo que realmente sale es"

# ╔═╡ 5f805e82-418b-11eb-05a4-dba1827925ec
expand(subs(dV1,x2=>ϕ1+z1)) 

# ╔═╡ 5eadc9c0-c926-11eb-2e4b-cb01e83f5e4f
md"""Resulta que el término $x_1z_1$ es **malo**. Para *cargármelo* necesito añadir a $V$ algo que al derivar produzca $z_1\cdot$"algo que pueda contolar" ... eso nos lleva al paso 2"""

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
md"Y derivemos a ver que sale, para ello primero calculo $\dot z=\dot x_2 - \dot \phi_1$"

# ╔═╡ 559eb340-28f3-11eb-3720-059fec2196ee
begin
	dϕ1=symbols("dϕ1")
	dz1=dx2-dϕ1
end;

# ╔═╡ 03ef90e0-4190-11eb-1f86-9f2bf2582d96
md"**Truco**: para no arrastrar *chorizos* creamos una variale que sea la derivada como haríamos **a mano**, luego ya calcularemos la verdadera derivada al final"

# ╔═╡ 7d9e7930-28ed-11eb-20fd-970954cf18b7
let
	dV2=derivada(V2,[x1 z1],[dx1 dz1])
	dV2=subs(dV2,x2=>z1+ϕ1) #La sustiucion es para hacer aparecer z1
	dV2=expand(dV2)
end

# ╔═╡ fb177d00-5fd7-11eb-0636-2f7f9dbd5598
k2=symbols("k2");

# ╔═╡ 2e63aed0-28f7-11eb-27a6-cbf8cb83c09c
ϕ2=-x1 -k2*z1 +dϕ1

# ╔═╡ da519e60-28f6-11eb-3598-c12f14f24180
md"""queremos hacer que $x_1z_1 + x_3z_1 - z_1 \dot ϕ_1 = -k_2z_1^2$, despejamos

$x_1 + x_3 - \dot ϕ_1 = -k_2z_1$

$x_3 = -k_2z_1 -x_1 + \dot ϕ_1$

Y hacemos $ϕ_2=-k_2z_1 -x_1 + \dot ϕ_1$

**Nota:** esto se podría despejar usando un paquete CAS más potente (Symbolics por ejemplo) pero no merece la pena, hacerlo a mano ayuda a *entender* lo que estamos haciendo. El CAS lo usaremos solamente para verificar los cálculos (que per se ya es muy útil)
"""

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
	dV2=subs(dV2,x3=>z2+ϕ2) #esto es lo que SI puedo hacer
	dV2=expand(dV2)
end

# ╔═╡ 8a020db0-c927-11eb-09b9-45afce72c1c8
md"""Al igual que antes queremos *cargarnos* $z_2z_1$ para ello hecesitamos un término adicional que de $z_2\cdot$"cosas" y a estas alturas os imaginaréis ya lo que vamos a poner ;)"""

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
end;

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
md"""queremos hacer que $uz_2 - z_2 \dot ϕ_2 + z_2z_1 = -k_3z_2^2$

$u - \dot ϕ_2 + z_1 = -k_3z_2$

$u =  \dot ϕ_2 - z_1 -k_3z_2$

Ahora vamos a ver que funciona:"""

# ╔═╡ db3aa380-5fd8-11eb-16fc-6b7ef56e8daa
k3=symbols("k3");

# ╔═╡ f8482150-4183-11eb-05f2-599741c800de
let
	dV3=derivada(V3,[x1 z1 z2],[dx1 dz1 dz2])
	dV3=subs(dV3,x2=>z1+ϕ1)
	dV3=subs(dV3,x3=>z2+ϕ2)
	dV3=subs(dV3,u=>dϕ2-z1-k3*z2)
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
	k1,k2,k3=p

	ϕ1=-k1*x1 - x1^2
	
	dϕ1=(-k1 - 2*x1)*(x2 + x1^2 - x1^3)
	
	z1=x2-ϕ1
	
	ϕ2= dϕ1 - x1 - k2*z1
	
	z2=x3-ϕ2
	
	dϕ2= -k2*(-dϕ1 + x3) + (-k1 - 2*x1)*x3 + (x2 + x1^2 - x1^3)*((-k1 - 2*x1)*(2*x1 - 3*x1^2) - 2*(x2 + x1^2 - x1^3)) - (x2 + x1^2 - x1^3)
	#si se hace a mano se puede sustituir todavía menos por ejemplo podemos darnos cuenta de que el término repetido x2 + x1^2 - x1^3 es dx2 y calcularlo una sola vez...
	
	u=dϕ2-z1-k3*z2
	
	#bonus para ver que funciona
	V=(x1^2 + z1^2 + z2^2)/2
	
	return (u,V)
end

# ╔═╡ c259a7d0-28fb-11eb-3cc6-ab3e44ea196c
function derivadas(x,p,t)
	x1,x2,x3=x
	u=control(x,p,t)[1]
	dx1=x1^2-x1^3+x2
    dx2=x3
	dx3=u
	dx=[dx1 dx2 dx3]
end

# ╔═╡ 4704bc10-c963-11eb-24d9-2dd79b9b5f33
#=prueba para ver que devuele números y por lo tanto hemos definido todo
si hemos dejado algo sin definir apareceran las variables simbólicas=#
control([1 2 3],[1 1 1],0)

# ╔═╡ f78bff70-4192-11eb-2f14-05fbee335f4f
md"""ha llegado el momento de hacer los cálculos"""

# ╔═╡ 19ce9560-c939-11eb-0678-6956501cddf9
ϕ1

# ╔═╡ 11015170-c939-11eb-3145-1ff24e53a7f8
dϕ1_b=derivada(ϕ1,[x1],[dx1])

# ╔═╡ b941ace2-c939-11eb-1414-a3befce0a650
ddϕ1=derivada(dϕ1_b,[x1 x2],[dx1 dx2])

# ╔═╡ 83fcb3a0-c938-11eb-1051-a93ef03e0639
ϕ2

# ╔═╡ 70c05310-c96e-11eb-26c8-334d2b13e99b


# ╔═╡ e2c93200-c938-11eb-2d02-752060d9d7d9
dϕ2_b=derivada(ϕ2,[x1 z1 dϕ1],[dx1 dz1 ddϕ1])

# ╔═╡ 4031d6c0-4196-11eb-275f-cda9a32cce36
md"## Vamos a ver cómo queda el *chorizo*"

# ╔═╡ 6ff23950-4195-11eb-250b-898515a5c031
chorizo_u=expand(control([x1 x2 x3],[k1,k2,k3],0)[1])

# ╔═╡ 9e20ff40-c963-11eb-16e2-cd955911eb41
md"## Bonus, depuración:
Para comprobar los cálculos los haremos de otra forma y veremos si concuerdan, para ello trabajaremos directamente con las x, primero ponemos $V$ en términos de $x_1$, $x_2$ y $x_3$"

# ╔═╡ ae7ec160-c963-11eb-07ae-2b1bb1cd78d4
begin
	V_x=V3
	V_x=subs(V_x,z1=>x2-ϕ1)
	V_x=subs(V_x,z2=>x3-ϕ2)
	V_x=subs(V_x,dϕ1=>dϕ1_b)
	V_x=subs(V_x,z1=>x2-ϕ1) 
	#El último parece redundante pero hay que insistir hasta que no queden zs
	#Algunas sustituciones hacen aparecer z1 de nuevo
end

# ╔═╡ 242cce30-c968-11eb-2c4b-4f2162e95f54
md"ahora derivamos y substituimos el conrtol quedando todo en términos de las x"

# ╔═╡ ae663050-c963-11eb-24d7-d3dcd3055b4b
begin
	dV_x=derivada(V_x,[x1 x2 x3],[dx1 dx2 dx3])
	dV_x=subs(dV_x,u=>chorizo_u)
end

# ╔═╡ 4ad29ab2-c968-11eb-26c6-f90433e890b1
md"la idea es comparar con lo que esperamos que salga"

# ╔═╡ ae4fe930-c963-11eb-30eb-f16f9e561190
begin
	dV_teórica=-k1*x1^2 - k2*z1^2 - k3*z2^2 - x1^4  #esto dice la teoría
	#para comparar lo ponemos en términos de las x solamente
	dV_teórica=subs(dV_teórica,z1=>x2-ϕ1)
	dV_teórica=subs(dV_teórica,z2=>x3-ϕ2)
	dV_teórica=subs(dV_teórica,dϕ1=>dϕ1_b)
	dV_teórica=subs(dV_teórica,z1=>x2-ϕ1)  #es necesario ya que vuelve a aparecer z1
end

# ╔═╡ cc2f0f30-c968-11eb-0814-bbe5107bdd3b
md"Si sale 0 es que lo hemos hecho bien"

# ╔═╡ ee6c43b0-c963-11eb-13f9-c51ef0a7f9a2
expand(dV_teórica-dV_x)

# ╔═╡ 9e7a974e-c96b-11eb-3a2b-759bf3f51058
md"# Y ahora a simular"

# ╔═╡ 49b79790-5fda-11eb-2917-4b055386f9db
k_1= @bind k_1 Slider(0:0.01:5,default=1,show_value=true)

# ╔═╡ 8ac8b250-5fda-11eb-2acd-4f317055f0cf
k_2= @bind k_2 Slider(0:0.01:5,default=1,show_value=true)

# ╔═╡ 8fbb0600-5fda-11eb-2f43-196140b76f69
k_3= @bind k_3 Slider(0:0.01:5,default=1,show_value=true)

# ╔═╡ fc5e8b30-4147-11eb-30d3-4bef0bac9ca6
begin
	x0 = [1 0 0]
	parametros=[k_1 k_2 k_3]
	tspan = (0.0,10.0)
	prob = ODEProblem(derivadas,x0,tspan,parametros);
	sol = solve(prob,Tsit5());
	
	Xs=sol.u
	ts=sol.t
	
	Vs=[]
	Us=[]
	
	for i=1:length(ts)
		u,V=control(Xs[i],parametros,ts[i])
		push!(Vs,V)
		push!(Us,u)
	end
	
	# Salida 
	fig1=plot(sol,idxs=(0,1), xaxis=L"t",yaxis=L"y(t)",label=L"y(t)")
		
	#Estados
	fig2=plot(sol, xaxis="t",yaxis=L"x(t)", label=:none)
	    
	#Señal de control
	fig3=plot(ts,Us,xaxis="t",yaxis=L"u(t)", label=:none)
	xlims!(tspan)
	
	#V
	fig4=plot(ts,Vs,xaxis="t",yaxis=L"V(t)", label=:none)
	xlims!(tspan)
	
	l = @layout [a ; b ; c; d]
	plot(fig1,fig2,fig3,fig4, layout=l)
end

# ╔═╡ Cell order:
# ╟─80213b30-28e7-11eb-388d-5501de9711e7
# ╠═42931bc0-4148-11eb-26b9-e5fee42c7580
# ╠═b5d53800-28ea-11eb-1858-4360a7d9fafc
# ╟─bd1a0730-28ea-11eb-2355-c1d9285733c5
# ╠═2890d210-28ea-11eb-2f16-23af078db24a
# ╠═e0c749a0-28e9-11eb-2a9d-3b5f6566605a
# ╟─4cce2a0e-28ea-11eb-26c6-672aeb7bfdf9
# ╠═646d6000-28ea-11eb-342d-b1c902896547
# ╠═224b2df0-28eb-11eb-2a45-c9eb740f149c
# ╟─2d6d8fe0-c925-11eb-3a32-2fd273a45a10
# ╠═2c671c60-c925-11eb-303d-b398113bdb9b
# ╟─4684dc70-28eb-11eb-1fde-2575c48354e8
# ╠═d589c21e-5fd5-11eb-3c70-535e951c59cf
# ╠═01b23ad0-c926-11eb-2dc3-f143c1d2b6f4
# ╠═5159e1d0-28ec-11eb-14f4-99de193e35a1
# ╟─bbd3f930-4181-11eb-1007-7fecf56f39aa
# ╠═220f19e0-28ec-11eb-057e-0b369e37930c
# ╟─883fd650-28ec-11eb-0ba3-014589882c54
# ╠═c5680400-418b-11eb-32bb-ad1a2f2077ed
# ╟─95969500-418d-11eb-0132-fd7f0ed24652
# ╠═5f805e82-418b-11eb-05a4-dba1827925ec
# ╟─5eadc9c0-c926-11eb-2e4b-cb01e83f5e4f
# ╟─42851430-28ed-11eb-2823-a5014cb3eda0
# ╠═984ebfb0-28ed-11eb-2cee-fdbf3d493de8
# ╟─a5582280-28f4-11eb-2c05-51acb08f478b
# ╠═559eb340-28f3-11eb-3720-059fec2196ee
# ╟─03ef90e0-4190-11eb-1f86-9f2bf2582d96
# ╠═7d9e7930-28ed-11eb-20fd-970954cf18b7
# ╠═fb177d00-5fd7-11eb-0636-2f7f9dbd5598
# ╠═2e63aed0-28f7-11eb-27a6-cbf8cb83c09c
# ╟─da519e60-28f6-11eb-3598-c12f14f24180
# ╟─0dc9ed20-4192-11eb-1039-93d7bfd5d09f
# ╠═ebea590e-4182-11eb-2a0b-a9161adc7f97
# ╟─a58593c0-28f7-11eb-1674-51101547d596
# ╠═3094ca2e-28f8-11eb-2cc0-a1d2cebc1da6
# ╠═3dbf9b10-4192-11eb-04e5-1f29ccd9da2a
# ╟─8a020db0-c927-11eb-09b9-45afce72c1c8
# ╟─fb7102a0-28f8-11eb-131d-c1792c1aaa27
# ╠═28437790-28f9-11eb-18f1-11364b9029ee
# ╠═9f114912-28f9-11eb-17f1-7fa461ff3371
# ╟─24332a00-28fa-11eb-2a88-8b5e6f124e6d
# ╠═85a6c900-28f9-11eb-229f-f78f9cd042e4
# ╟─c9b94770-28fa-11eb-25ba-2384e1b7672a
# ╠═db3aa380-5fd8-11eb-16fc-6b7ef56e8daa
# ╠═f8482150-4183-11eb-05f2-599741c800de
# ╟─f16d5e90-4183-11eb-16dc-8191ac698ce2
# ╠═c259a7d0-28fb-11eb-3cc6-ab3e44ea196c
# ╠═d4991de0-4147-11eb-067b-4f9c143fb673
# ╠═4704bc10-c963-11eb-24d9-2dd79b9b5f33
# ╟─f78bff70-4192-11eb-2f14-05fbee335f4f
# ╠═19ce9560-c939-11eb-0678-6956501cddf9
# ╠═11015170-c939-11eb-3145-1ff24e53a7f8
# ╠═b941ace2-c939-11eb-1414-a3befce0a650
# ╠═83fcb3a0-c938-11eb-1051-a93ef03e0639
# ╠═70c05310-c96e-11eb-26c8-334d2b13e99b
# ╠═e2c93200-c938-11eb-2d02-752060d9d7d9
# ╟─4031d6c0-4196-11eb-275f-cda9a32cce36
# ╠═6ff23950-4195-11eb-250b-898515a5c031
# ╟─9e20ff40-c963-11eb-16e2-cd955911eb41
# ╠═ae7ec160-c963-11eb-07ae-2b1bb1cd78d4
# ╟─242cce30-c968-11eb-2c4b-4f2162e95f54
# ╠═ae663050-c963-11eb-24d7-d3dcd3055b4b
# ╟─4ad29ab2-c968-11eb-26c6-f90433e890b1
# ╠═ae4fe930-c963-11eb-30eb-f16f9e561190
# ╟─cc2f0f30-c968-11eb-0814-bbe5107bdd3b
# ╠═ee6c43b0-c963-11eb-13f9-c51ef0a7f9a2
# ╟─9e7a974e-c96b-11eb-3a2b-759bf3f51058
# ╟─49b79790-5fda-11eb-2917-4b055386f9db
# ╟─8ac8b250-5fda-11eb-2acd-4f317055f0cf
# ╟─8fbb0600-5fda-11eb-2f43-196140b76f69
# ╟─fc5e8b30-4147-11eb-30d3-4bef0bac9ca6
