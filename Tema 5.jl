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

# ╔═╡ b5a8675e-0e2a-11eb-1787-55285590d7de
begin
    import Pkg
	  Pkg.activate(@__DIR__)
	using Plots, PlutoUI, SymEngine 
#hay muchos paquete simbólicos, el más completo (límites integrales, sustiuciones...) es SymPy, Symbolics es nuevo y en desarrollo y SymEngine es simple y suficiente para lo que vamos a hacer en este tema, usad el que más os convenga
end

# ╔═╡ bcd20d70-0e2a-11eb-34c1-5f49700df9cd
using LaTeXStrings

# ╔═╡ 5aded872-b8b4-11eb-0277-93d61bffe2fd
md"# Tema 5 Control en modo deslizante"

# ╔═╡ f8105520-250b-11eb-1a6b-0129a812bcee
md"""# Vamos a crear nuetro propio simulador con *Euler*
¿En serio? ¿Euler no era un método malo?
No es el mejor método pero... los otros tampoco funcionan bien si la función es discontinua y a nuestro *le podemos hacer perrerías* (usar funciones discontinuas, introducir ruido, etc...)
"""


# ╔═╡ 96dbe960-2509-11eb-1248-b7c75a8974dd
#uso euler por que en un sistema discontinuo el resto de solvers no ofrece mucha ventaja y así tenemos más control sobre lo que está pasando

function simular(x0,planta,control,parametros,tf,Δt)
	t=0:Δt:tf
	estado=copy(x0);
	x=[]
	push!(x,estado)
	u=Float64[control(estado,t[1],parametros)[1]]
	
	for n=1:(length(t)-1)
		señal_de_control=control(estado,t[n],parametros)[1] #el [1] se entenderá luego
		derivada=planta(estado,señal_de_control,t[n],parametros)
		estado=estado+derivada*Δt #euler
		push!(x,estado)
		push!(u,señal_de_control)
	end
	
	return (x,u,t)
end

# ╔═╡ ba989850-0a31-11eb-1f3a-e363a7e0207e
md"""# Primer ejemplo (primer orden):

$\dot x=sin(t)+acos(x)+u=f(x,t)+u$

$a \in [0.5,0.9]$

Luego sabemos que

$|f(x,t)|<|sin(t)|+|acos(x)| \le 1+|a| \le 1+0.9=1.9=\rho$

## Hay que aplicar un $u$ mayor que $\rho$ y con signo opuesto a $x$... luego haremos esto:

$u=-sign(x)(\rho+\epsilon)=-2sign(x)$

Donde hemos tomado $\epsilon=0.1$
"""

# ╔═╡ 94294500-0e2a-11eb-3a25-617242dfe932
function planta_ejemplo1(x,u,t,p)
	a=p[1]#parámetro desconocido
	#OJO no hagáis a=p aunque haya un paámetro no es lo mismo 1 que [1]
	dx=sin(t)+a*cos(x)+u
end

# ╔═╡ 1d402de0-0e2b-11eb-10e8-79bdc072ee9b
control_1a(x,t,p)=-2*sign(x)

# ╔═╡ 3b89eed0-0e2b-11eb-3cad-8d86b96accc5
@bind x1_0 Slider(-1:0.1:1)

# ╔═╡ 6b48fc60-0e2b-11eb-1ad1-075f865e4c78
@bind a1 Slider(0.5:0.01:0.9)

# ╔═╡ ff0eafe0-0e2a-11eb-3e1c-5dd70c43c3c1
let
	Δt=0.001
	tf=2
	x0=x1_0
	parametros=[a1]
	x,u,t=simular(x0,planta_ejemplo1,control_1a,parametros,tf,Δt)
	
	fig1=plot(t,x, xaxis=L"t",yaxis=L"x(t)",label=:none)
	title!(fig1,"a=$a1")
	fig2=plot(t,u, xaxis=L"t",yaxis=L"u(t)",legend=false)
	l = @layout [a ; b]
	plot(fig1,fig2, layout=l)
end

# ╔═╡ f77e1fc0-0e2c-11eb-0051-d7eee614b171
md"""# Hemos hecho *trampa*, ¿Dónde está el control equivalente?
El control equivalente es el que mantiene $s=x$ en cero es decir:

$x(t)=0 \to  \dot x = 0 =sin(t)+acos(x)+u=sin(t)+acos(0)+u$
$0=sin(t)+a+u \to u=-a-sin(t)$

No está ni hace falta, ¡Nunca se ejecuta! Vamos a verlo
"""

# ╔═╡ e3723e20-0efa-11eb-1331-e710fbcd289f
function control_1b(x,t,p)
	if x>0 
		u=-2
	elseif x<0
		u=2
	else
		u=1000000 #NO es el control equivalente es para ver que sucede ;)
	end	
	return u
end

# ╔═╡ 155f6de0-25c6-11eb-3a8e-851538372b3e
let
	Δt=0.001
	tf=2
	x0=x1_0
	parametros=[a1]
	x,u,t=simular(x0,planta_ejemplo1,control_1b,parametros,tf,Δt)
	
	fig1=plot(t,x, xaxis=L"t",yaxis=L"x(t)",label=:none)
	fig1=hline!([0],line=:dot,label=:none)
	title!(fig1,"a=$a1")
	fig2=plot(t,u, xaxis=L"t",yaxis=L"u(t)",legend=false)
	l = @layout [a ; b]
	plot(fig1,fig2, layout=l)
end

# ╔═╡ 7cc77a40-0efb-11eb-145d-4b535e701cf8
md"""no ha pasado nada, u nunca vale 1000000
... a no ser que x sea exactamente cero lo que es muy raro a no ser que lo forcemos
(poniendo $x_0=0$) pero ni así tiene mucha improtancia

## Significado del control equivalente

Aunque no lo veamos el control equivalente es "en promedio" igual que el control que estamos aplicando, vamos a verlo
""" 

# ╔═╡ 29de3480-25c7-11eb-1d21-2107f4712127
media_móvil(v,n) = [sum(@view v[i:(i+n-1)])/n for i in 1:(length(v)-(n-1))]

# ╔═╡ 1ce32ff0-25c8-11eb-1371-e360c5c3ed2b
@bind N Slider(1:30)

# ╔═╡ 674476e0-2810-11eb-000b-0742e1bce381
md"Ver real $(@bind ver CheckBox())"

# ╔═╡ e32619f0-0efa-11eb-020d-c9743701e335
let
	Δt=0.001
	tf=2
	x0=x1_0
	parametros=[a1]
	x,u,t=simular(x0,planta_ejemplo1,control_1a,parametros,tf,Δt)

	u_equivalente=-a1.-sin.(t)
	
	u_promedio=media_móvil(u,N)
	t_promedio=media_móvil(t,N)
	
	fig1=plot()
	if ver==true
	   fig1=plot!(fig1,t,u, xaxis=L"t",yaxis="ues",label="u_real")
	end
	fig1=plot!(fig1,t,u_equivalente,label="u_equivalente")
	fig1=plot!(fig1,t_promedio,u_promedio,label="u_promedio")
end

# ╔═╡ 14080630-0f07-11eb-068d-09c4055d081b
md"""# Regularización y control equivalente
En vez de usar el signo que es discontinuo podemos usar una versión continua del mismo
"""

# ╔═╡ f47ce110-24fe-11eb-036a-e5d8bad1955f
@bind ϵ Slider(0.01:0.01:0.2)   #\epsilon + TAB

# ╔═╡ 1f38d210-24ff-11eb-3cbe-9b1e70108439
let
	ordenadas=-1:0.01:1
	fig1=plot(ordenadas,sign.(ordenadas), xaxis=L"x",label="sign(x)")
	fig1=plot!(fig1,ordenadas,tanh.(ordenadas/ϵ), label="tahn(x/$ϵ)")
	
end

# ╔═╡ e71e9160-2500-11eb-1b6c-d59e43143b66
md"Veamos que efecto tiene camniar el signo por su vesion continua"

# ╔═╡ 24a65e50-2501-11eb-2363-67d7e5bcb405
control_1c(x,t,p)=-2*tanh(x/p[2])

# ╔═╡ e023ea92-2500-11eb-00ff-4f27947f46d4
let
	Δt=0.001
	tf=2
	x0=x1_0
	parametros=[a1,ϵ]
	x,u,t=simular(x0,planta_ejemplo1,control_1c,parametros,tf,Δt)
	
	u_equivalente=-a1.-sin.(t)
		
	fig1=plot(t,x, xaxis=L"t",yaxis=L"x(t)",label=:none)
	title!(fig1,"a=$a1 ϵ=$ϵ")
	fig1=hline!([0],line=:dot,label=:none)
	
	fig2=plot(t,u, xaxis=L"t",yaxis=L"u_{aplicado}",legend=false)
	fig2=plot!(fig2,t,u_equivalente,
		    xaxis=L"t",yaxis=L"u_{es}",label=L"u_{equivalente}")

	l = @layout [a ; b]
	plot(fig1,fig2, layout=l)
end

# ╔═╡ e01fcbe0-2500-11eb-3b94-a7662624c5ed
md"""## Entonces lo mejor es disminuir $\epsilon$ todo lo posible ¿No?

*No exactamente*, en la práctica tiene sus problemas, por ejemplo si la medida del estado es ruidosa...
Veámoslo:
"""

# ╔═╡ 16e8dac2-25cb-11eb-33d1-a95dc487b34d
function control_1d(x,t,p)
	_,ϵ,nivel_ruido=p
	ruido=nivel_ruido*2*(rand()-0.5)
	x_medido=x+ruido
	u=-2*tanh(x_medido/ϵ)
end

# ╔═╡ 612d43c0-2506-11eb-3a14-8145b79e349f
@bind ϵ2 Slider(0.01:0.01:0.2)  

# ╔═╡ d8e35250-25ca-11eb-3d59-b3c6396a4fb1
@bind nivel_ruido Slider(0:0.01:0.2)  

# ╔═╡ 5a31c47e-2504-11eb-37af-07b3b88360eb
let
	Δt=0.001
	tf=2
	x0=x1_0
	parametros=[a1,ϵ2,nivel_ruido]
	x,u,t=simular(x0,planta_ejemplo1,control_1d,parametros,tf,Δt)
	
	u_equivalente=-a1.-sin.(t)
		
	fig1=plot(t,x, xaxis=L"t",yaxis=L"x(t)",label=:none)
	title!(fig1,"a=$a1, ϵ=$ϵ, nivel del ruido=$nivel_ruido")
	fig1=hline!([0],line=:dot,label=:none)
	
	fig2=plot(t,u, xaxis=L"t",yaxis=L"u_{aplicado}",legend=false)
	fig2=plot!(fig2,t,u_equivalente,
		    xaxis=L"t",yaxis=L"u_{es}",label=L"u_{equivalente}")

	l = @layout [a ; b]
	plot(fig1,fig2, layout=l)
end

# ╔═╡ f694c060-2505-11eb-26ef-1fc536e2260d
md"""# Ejemplo 2 vamos al tracking directamente ☺		
La estabilización es como el tracking pero con $r=0$, como ya lo vimos en el tema anterior podemos ir más rápido en este.

$\dot x_1=sin(x_1)+x_2$

$\dot x_2=ax_1+bx_2^2+cu$

$y=x_1$

con  $a \in [-1,7]$,  $b \in [1,2]$ y $c \in [1.5,2.8]$

Queremos seguir la referenica $y_r=sin(t)$

"""

# ╔═╡ 595542f0-25cf-11eb-308d-bb524a624b04
md"Defimos las función con una interfaz igual que en el tema 4"

# ╔═╡ 41d2207e-25cf-11eb-1cbf-a599e3a5d2fc
#la diferencia es que Symengine tiene impmementada diff directamente
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

# ╔═╡ 32f2cfd0-2505-11eb-36a3-cdaec55a91cc
x1,x2,u,a,b,c=symbols("x1,x2,u,a,b,c") #o @vars x1... pero inicializando antes

# ╔═╡ b1c43a50-25ce-11eb-3478-f13e9174ccc9
begin
    dx1=sin(x1)+x2
    dx2=a*x1+b*x2^2+c*u
	dx=[dx1 dx2]
	y=x1
end

# ╔═╡ 992e79f0-25cf-11eb-0885-9f431d9e2ed6
dy=derivada(y,[x1 x2],[dx1 dx2],1)

# ╔═╡ 806542f0-25cf-11eb-39c8-db55ed5e86a7
ddy=derivada(y,[x1 x2],[dx1 dx2],2)

# ╔═╡ c5f21a70-25d2-11eb-2af8-395671e44537
md"Ok tiene el grado correcto podemos usar la técnica"

# ╔═╡ 2cd935f0-25d0-11eb-1b07-a7002584aa16
md"""## Defimos ahora el error y su dinámica

Recordemos $e_1 = y - y_r$  y ahora a calcular sus derivadas:
"""

# ╔═╡ 77f468f0-2816-11eb-2dc5-25eb17f3a834
yr,dyr,ddyr=symbols("yr,dyr,ddyr")

# ╔═╡ e1f2879e-25d2-11eb-04d4-b53a77d90c5f
e1=yr-y

# ╔═╡ 14e23520-2817-11eb-2f93-0779085f77a5
md"""La teoría es que $S=(\frac{d}{dt}+k)^{n-1}e_1$

Para el caso n=2

$S=\dot e_1 +ke_1$"""

# ╔═╡ 37fe65b0-2817-11eb-2cdb-2921b513978c
k=symbols("k")

# ╔═╡ 9ac55fb0-2816-11eb-12da-a9ddac1c913c
de1=derivada(e1,[x1 x2 yr dyr],[dx1 dx2 dyr ddyr])

# ╔═╡ 0ed1a02e-2817-11eb-0b92-e5feb0f17c63
S=de1+k*e1

# ╔═╡ a3caae20-2817-11eb-09e6-59d6464cef2a
dS=derivada(S,[x1 x2 yr dyr],[dx1 dx2 dyr ddyr])

# ╔═╡ aff35570-25cf-11eb-1627-43fe902b9db5
md"Identificamos términos $\dot S=-f-gu$"

# ╔═╡ a8a18df2-25cf-11eb-38b7-e128e53ba029
g=-diff(dS,u)

# ╔═╡ 081bb16e-25d0-11eb-28e6-c3e5a042fce2
f=-(dS+g*u) #el resto vamos.... OJO A LOS SIGNOS

# ╔═╡ f77c9930-25d2-11eb-287f-655c63ea721e
md"""Acotar $|f/g|$ un proceso **artesanal** se puede hacer **de muchas formas**
recordemos que $ a \in [-1,7] $,  $ b \in [1,2] $ y $ c \in [1.5,2.8] $
"""

# ╔═╡ cf8b68e2-b8ba-11eb-0177-8fa5c7f332bd
expand(f/g)

# ╔═╡ 097834c0-b8bb-11eb-3e52-c1d543ef9fb6
md"Y ahora vamos a acotar"

# ╔═╡ 9c2a7480-25d2-11eb-0479-751c0dc8b88e
rho=(abs(ddyr) + 7*abs(x1) + 2*x2^2 + k*abs(dyr)+ k*abs(x2) + k + abs(x2) + 1)/1.5

# ╔═╡ abc0604e-2823-11eb-1708-83268e22b7bc
md"## Todo junto"

# ╔═╡ 2b9f72b0-281b-11eb-1fbc-15b3e67b017d
referencia(t)=0.0.*(sin(t),cos(t),-sin(t))

# ╔═╡ 42884660-281a-11eb-015e-29ff1b21dfe0
function control_seguimiento(x,t,p)
	ϵ=0.01;
    x1,x2=x
	ε=0.1;#el epsilon bonito es varepsilon
    k=1;#una elección arbitraria
    yr,dyr,ddyr=referencia(t)
	
	y=x1; #mucho cuidado si no se define puede usar las globales y se cuelga
	e1=yr-y
	de1=dyr - (x2 + sin(x1))
	
	S=de1+k*e1
	
	#RELLENAR ESTO
	u=tanh(S / ϵ)*((abs(ddyr) + 7*abs(x1) + 2*x2^2 + k*abs(dyr)+ k*abs(x2) + k + abs(x2) + 1)/1.5 + ε)

	return (u,S,e1,de1)  #Saco todo esto para pintarlo luego fácilmente 
	#por eso la función simular se queda SOLO CON u
end

# ╔═╡ 9acf5880-281b-11eb-3aa5-9b6df9ac7501
function planta_ejemplo2(x,u,t,p)
	a,b,c=p
	x1,x2=x #OJO a las variables, si por ejemplo no las ponemos...
    dx1=sin(x1)+x2
    dx2=a*x1+b*x2^2+c*u
	dx=[dx1 dx2]
end

# ╔═╡ e8ee4712-281b-11eb-0fbf-7783935a0075
@bind a2 Slider(-1:0.1:7)

# ╔═╡ 0e04e810-281c-11eb-3e72-c58bf8d7ce0d
@bind b2 Slider(1:0.1:2)

# ╔═╡ 1e143c10-281c-11eb-1117-29f72834203e
@bind c2 Slider(1.5:0.1:2.8)

# ╔═╡ 756a3570-2833-11eb-2596-c37768e5a75b
@bind x1_0b Slider(-1:0.1:1)

# ╔═╡ be9d11e0-2833-11eb-1e9f-d3bfc428c893
@bind variables Select(["estados"=>"Ver estados", "errores"=>"Ver errores"])

# ╔═╡ 8bdbadb0-281b-11eb-1786-3bef9aee1d39
let
	Δt=0.001
	tf=10
	x0=[x1_0b 0]
	parametros=[a2,b2,c2]
	estados,u,t=simular(x0,planta_ejemplo2,control_seguimiento,parametros,tf,Δt)

	#hay que sacar los datos de los estados, es decir parsar de un vecror de vectores a vectores sencillos
	x1=[x[1] for x in estados]
	x2=[x[2] for x in estados]
	y=x1
	
	refs=referencia.(t)
	yr=[ref[1] for ref in refs]
	dyr=[ref[2] for ref in refs]
	
	#truco para pintar más cosas sin tener que reescribirlas
	Ss=[]
	e1s=[]
	de1s=[]
	for i=1:length(t)
		_,S,e1,de1=control_seguimiento(estados[i],t[i],parametros)
		push!(Ss,S)
		push!(e1s,e1)
        push!(de1s,de1)
	end
	
	if variables=="estados"
		plot(t,x1, xaxis=L"t",yaxis=L"x(t)",label=L"y=x_1")
	    plot!(t,x2, xaxis=L"t",label=L"x_2")
	    title!("a=$a2, b=$a2, c=$a2")
	    fig1=plot!(t,yr,line=:dot,label=L"y_r")
	else
		plot(t,e1s, xaxis=L"t",yaxis=L"e(t)",label=L"e_1")
	    plot!(t,de1s, xaxis=L"t",label=L"e_2=\dot e_1")
        title!("a=$a2, b=$a2, c=$a2")
		fig1=hline!([0],line=:dot,label=:none)
	end
	
	fig2=plot(t,Ss, xaxis=L"t",yaxis=L"s",legend=false)
	
	fig3=plot(t,u, xaxis=L"t",yaxis=L"u",legend=false)

    if variables=="estados"
		fig4=plot(x1,x2, xaxis=L"x_1", yaxis=L"x_2",label=:none)
	else
		fig4=plot(e1s,de1s, xaxis=L"e_1", yaxis=L"e_2",label=:none)
	end

	l = @layout [a ; b ;c ; d]
	plot(fig1,fig2,fig3, fig4, layout=l, size=(600,800))
	
end

# ╔═╡ Cell order:
# ╟─5aded872-b8b4-11eb-0277-93d61bffe2fd
# ╟─f8105520-250b-11eb-1a6b-0129a812bcee
# ╠═96dbe960-2509-11eb-1248-b7c75a8974dd
# ╟─ba989850-0a31-11eb-1f3a-e363a7e0207e
# ╠═b5a8675e-0e2a-11eb-1787-55285590d7de
# ╠═bcd20d70-0e2a-11eb-34c1-5f49700df9cd
# ╠═94294500-0e2a-11eb-3a25-617242dfe932
# ╠═1d402de0-0e2b-11eb-10e8-79bdc072ee9b
# ╠═3b89eed0-0e2b-11eb-3cad-8d86b96accc5
# ╠═6b48fc60-0e2b-11eb-1ad1-075f865e4c78
# ╠═ff0eafe0-0e2a-11eb-3e1c-5dd70c43c3c1
# ╟─f77e1fc0-0e2c-11eb-0051-d7eee614b171
# ╠═e3723e20-0efa-11eb-1331-e710fbcd289f
# ╠═155f6de0-25c6-11eb-3a8e-851538372b3e
# ╟─7cc77a40-0efb-11eb-145d-4b535e701cf8
# ╠═29de3480-25c7-11eb-1d21-2107f4712127
# ╠═1ce32ff0-25c8-11eb-1371-e360c5c3ed2b
# ╟─674476e0-2810-11eb-000b-0742e1bce381
# ╠═e32619f0-0efa-11eb-020d-c9743701e335
# ╟─14080630-0f07-11eb-068d-09c4055d081b
# ╟─1f38d210-24ff-11eb-3cbe-9b1e70108439
# ╠═f47ce110-24fe-11eb-036a-e5d8bad1955f
# ╠═e71e9160-2500-11eb-1b6c-d59e43143b66
# ╠═24a65e50-2501-11eb-2363-67d7e5bcb405
# ╟─e023ea92-2500-11eb-00ff-4f27947f46d4
# ╟─e01fcbe0-2500-11eb-3b94-a7662624c5ed
# ╠═16e8dac2-25cb-11eb-33d1-a95dc487b34d
# ╠═5a31c47e-2504-11eb-37af-07b3b88360eb
# ╠═612d43c0-2506-11eb-3a14-8145b79e349f
# ╠═d8e35250-25ca-11eb-3d59-b3c6396a4fb1
# ╟─f694c060-2505-11eb-26ef-1fc536e2260d
# ╠═595542f0-25cf-11eb-308d-bb524a624b04
# ╟─41d2207e-25cf-11eb-1cbf-a599e3a5d2fc
# ╠═32f2cfd0-2505-11eb-36a3-cdaec55a91cc
# ╠═b1c43a50-25ce-11eb-3478-f13e9174ccc9
# ╠═992e79f0-25cf-11eb-0885-9f431d9e2ed6
# ╠═806542f0-25cf-11eb-39c8-db55ed5e86a7
# ╟─c5f21a70-25d2-11eb-2af8-395671e44537
# ╟─2cd935f0-25d0-11eb-1b07-a7002584aa16
# ╠═77f468f0-2816-11eb-2dc5-25eb17f3a834
# ╠═e1f2879e-25d2-11eb-04d4-b53a77d90c5f
# ╟─14e23520-2817-11eb-2f93-0779085f77a5
# ╠═37fe65b0-2817-11eb-2cdb-2921b513978c
# ╠═9ac55fb0-2816-11eb-12da-a9ddac1c913c
# ╠═0ed1a02e-2817-11eb-0b92-e5feb0f17c63
# ╠═a3caae20-2817-11eb-09e6-59d6464cef2a
# ╟─aff35570-25cf-11eb-1627-43fe902b9db5
# ╠═a8a18df2-25cf-11eb-38b7-e128e53ba029
# ╠═081bb16e-25d0-11eb-28e6-c3e5a042fce2
# ╟─f77c9930-25d2-11eb-287f-655c63ea721e
# ╠═cf8b68e2-b8ba-11eb-0177-8fa5c7f332bd
# ╟─097834c0-b8bb-11eb-3e52-c1d543ef9fb6
# ╠═9c2a7480-25d2-11eb-0479-751c0dc8b88e
# ╟─abc0604e-2823-11eb-1708-83268e22b7bc
# ╠═2b9f72b0-281b-11eb-1fbc-15b3e67b017d
# ╠═42884660-281a-11eb-015e-29ff1b21dfe0
# ╠═9acf5880-281b-11eb-3aa5-9b6df9ac7501
# ╠═e8ee4712-281b-11eb-0fbf-7783935a0075
# ╠═0e04e810-281c-11eb-3e72-c58bf8d7ce0d
# ╠═1e143c10-281c-11eb-1117-29f72834203e
# ╠═756a3570-2833-11eb-2596-c37768e5a75b
# ╠═be9d11e0-2833-11eb-1e9f-d3bfc428c893
# ╠═8bdbadb0-281b-11eb-1786-3bef9aee1d39
