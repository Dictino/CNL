### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1be7c0c0-d4ba-11ea-1e43-ef66c35b0830
using Plots

# ╔═╡ 159f28e2-d4db-11ea-13e6-a97e2e7f7c60
#using DifferentialEquations
using OrdinaryDiffEq

# ╔═╡ 6cbb32f0-d665-11ea-3f24-41cbe1910b80
using LinearAlgebra

# ╔═╡ 6ba02e80-d4cf-11ea-19cb-33028d3dc67d
md"# _Cosas que ver en el primer tema_
+ ¿Qué es una ecuación diferencial y como se resuelve?
    * Euler
    * Métodos mejores
+ Sistema de segundo orden con dos polos
+ Explosión en tiempo finito
+ Linealización
+ Ejemplos
    * Barco
    * Generador
+ Ciclo límite
    * Van der Pool
+ Sistemas discretos
    * Ecuación logística (en otro libro)"

# ╔═╡ 7db3fdb0-d4d2-11ea-21fc-597c55dda572
md"# ¿Qué es una ecuación diferencial?

A nosotros nos interesan sistemas como este:

$$\frac{dx}{dt}=f(x,t)$$

$$x(0)=x_0$$

Donde x es un vector y la derovada respecto al tiempo se suele escribir como $$\dot x$$ "

# ╔═╡ 76079a70-d4cf-11ea-2645-8f43c78eae52
md" ## Euler
Si la función $$f$$ tiene buen comportamiento (lo veremos más adelante cuando repasemos la existencia y unicidad) hay una solución de la ecuación que es contínua y derivable y usando taylor:

$$x(\Delta t)=x(0) + f(X0)\Delta t + o(\Delta t^2)$$

Si se toman pasos pequeños y repetimos el procedimiento:

$$x(2\Delta t)=x(\Delta t) + f(x(\Delta t),\Delta t) \Delta t + o(\Delta t^2)$$
...

y llamando $$t_n$$ a $$n\Delta t$$ y $$x_n$$ a $$x(t_n)$$

$$x_{n+1} \approx x_n + f(x_n,t_n)\Delta t + o(\Delta t^2)$$

### Ejemplo

"

# ╔═╡ 66e9ffd0-d4d7-11ea-09c0-03ebaacab879
function Euler(f,x0,tf,Δt)
	t=0:Δt:10
	x=zeros(size(t))
	x[1]=x0
	for n=1:(length(t)-1)
		x[n+1]=x[n]+f(x[n],t[n])*Δt
	end
	return (x,t)
end

# ╔═╡ e885f5b0-d4d4-11ea-06d6-9ddf9619961f
begin
	f(x,t)=-x+sin(t)
	f(u,p,t)=f(u,t) #para luego
end

# ╔═╡ b99fc9b0-d4d9-11ea-0e18-cfd5af90188f
tf=10

# ╔═╡ 3be6b020-d4d8-11ea-3356-23afd148ebc9
@bind N html"<input type=range min=1 max=100>"

# ╔═╡ 220f6f12-d4d9-11ea-159e-c1e1551fd732
Δt=tf/N;

# ╔═╡ 23f5ece0-d4da-11ea-1330-c5d077bd3a89
md"N=$(N) $$\to \Delta t=$$ $(Δt)"

# ╔═╡ d8b5ed40-d4d7-11ea-27d4-bf594518802b
x,t =Euler(f,0,10,Δt);

# ╔═╡ 2b7df180-d4d8-11ea-118f-75e8709e36c3
plot(t,x,label="x")

# ╔═╡ 455e0c30-d4d7-11ea-22ce-0d04b4e8d860
md"## Hay métodos mejores:

-Ajustan los pasos automaticamente

-Tienen más precisión para los mismos pasos

-Está integrado con otras funciones como plots, permite interpolar los resultados automáticamente, etc..."

# ╔═╡ 04501380-d63c-11ea-2ed8-e11f02de80a4


# ╔═╡ c90980d0-d637-11ea-112d-c336b2080d87
begin
	u0 = 0
	tspan = (0.0,10.0)
	prob = ODEProblem(f,u0,tspan);
	sol = solve(prob,Tsit5(), reltol=1e-8, abstol=1e-8);
	sol.u
end

# ╔═╡ 4a7f5f72-0e0c-11eb-2b7c-4d7214c09c4b
sol(2.751) #¡se puede usar como una función continua de t!

# ╔═╡ 51059990-0e0c-11eb-2d55-89ac739619c7


# ╔═╡ a92de300-d63c-11ea-2446-5b099592d734
plot(sol, xaxis="t",yaxis="x(t)",label="x") # legend=false

# ╔═╡ 09980bd0-d63d-11ea-0752-1d033ee00cc9
md"## Sistema de segundo orden con dos polos"

# ╔═╡ e6afc70e-d63d-11ea-0d36-738b13f0de30
function f_lineal(x,p,t)
	dx=[0.0;0.0]
	dx[1]=x[2];
	dx[2]=+p[1]*x[1]+p[2]*x[2]
	return dx
end

# ╔═╡ 22398f60-d665-11ea-3352-714b467152b3
#eigen([0 1; -(real^2+imag^2) -2*real ])

# ╔═╡ 57ffd1f0-d637-11ea-1a81-7dd574291417
md""" `es real = ` $(@bind es_real html"<input type=checkbox >")"""

# ╔═╡ 5b26b050-d7eb-11ea-2488-b7da296fed24
@bind slider1 html"<input type=range min=-1 max=1 step=0.05>"

# ╔═╡ 84c3c3d0-d7eb-11ea-3f35-29559df926c3
@bind slider2 html"<input type=range min=-1 max=1 step=0.05>"

# ╔═╡ a0573fe0-d6a2-11ea-057e-5140c323c005
begin
	if es_real==true #el polinomio es (s-real1)*(s-real2)=s^2+(-real1-real2)s + real1*real2 --> k2= real1 + real2 k1=real1*real2
		k2= slider1 + slider2
		k1=-slider1*slider2
	else #el polinomio es (s-real-imag*im)*(s-real+imag*im)
		# es decir s^2 - 2*real*s + real^2 + imag^2 --> k2=-2*real k1=-real^2-imag^2
		k2=2*slider1
		k1=-slider1^2 - slider2^2
	end
	A=[0 1; k1 k2]
	autovalores,_=eigen(A)
	md"""
	k1= $(k1)
	k2= $(k2)
	autovalores= $(autovalores)"""
end

# ╔═╡ 5a802550-d665-11ea-3f9f-b103ce80822f
scatter(real(autovalores),imag(autovalores),grid=true,legend=false, framestyle = :origin )

# ╔═╡ 1790c760-0e0d-11eb-1c47-f5e902e158fb


# ╔═╡ f5ca89e0-0e0c-11eb-32a1-fb22511f42ae
function flechas!(figura,f,p,rangos;N=10)
	#modifica una figura añadiendo las flechas por eso tiene una exclamación !
	xmin,xmax,ymin,ymax=rangos
	long=max(abs(xmax-xmin),abs(ymax-ymin))/N
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax

    #barrido para las flechas
	X1=[];
	X2=[];
	U=[];
	V=[];
	for x1 in x1s #barrido para calcular las flechas
		for x2 in x2s
			push!(X1,x1);
			push!(X2,x2);
			F=f([x1,x2],p,0)
			X=0.1*F/(0.1+sqrt(F[1]^2+F[2]^2))
			push!(U,X[1]);
			push!(V,X[2]);
		end
	end
	quiver!(figura,X1,X2,quiver=(U,V), arrow=arrow(0.1*long, 0.1*long)) #añade
	return figura
end

# ╔═╡ 165464e2-0e0f-11eb-3851-3701a6c78502
#= Es un euler modificado que integra hacia delante y hacia atrás hasta que se sale o completa un número de pasos. 
	-hacia delante se detectan bien los PE y ciclos límite estables
	-Hacia atrás permite ver también los intestables
	Ajusta dt para que el paso tenga un tamaño que sea aproximadamente =#

function trayectoria(x0,f,p,rangos;N=10,pasosMax=100)
	xmin,xmax,ymin,ymax=rangos
	long=max(abs(xmax-xmin),abs(ymax-ymin))/N
	x,y=x0 #saco las componentes de un vector
	pasos=0
	
	X1=[x]
	Y1=[y]
	X2=[x] #parte de la trayectoria hacia atrás en el tiempo
	Y2=[y]
	
	while (xmin<x<xmax) & (xmin<x<xmax) & (pasos<pasosMax)
		derivada=f([x,y],p,0)
		dt=0.05*long/(0.1+sqrt(derivada[1]^2+derivada[2]^2))
		Xn=[x;y]+derivada*dt
		x,y=Xn
		pasos=pasos+1
		push!(X1,x)
		push!(Y1,y)
	end
	pasos=0
	x,y=x0
	while (xmin<x<xmax) & (xmin<x<xmax) & (pasos<pasosMax)
		derivada=f([x,y],p,0)
		dt=-0.05*long/(0.1+sqrt(derivada[1]^2+derivada[2]^2))
		Xn=[x;y]+derivada*dt
		x,y=Xn
		pasos=pasos+1
		push!(X2,x)
		push!(Y2,y)
	end
	X=append!(reverse(X2),X1) #junto ambas partes
	Y=append!(reverse(Y2),Y1)
	return (X,Y)
end

# ╔═╡ 915aa560-0e0e-11eb-18d4-e3a899ee2c1c
function trayectorias!(figura,f,p,rangos;N=10)
	#modifica una figura añadiendo las trayectorias por eso tiene una exclamación !
	    #barrido para las soluciones
	xmin,xmax,ymin,ymax=rangos
	N=round(N/2);
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax
	T=10;
	for x1 in x1s
		for x2 in x2s 
			X,Y=trayectoria([x1;x2],f,p,rangos,N=N,pasosMax=100)
			figura=plot!(figura,X,Y,legend=false)
		end
	end 
	return figura
end

# ╔═╡ 7115307e-d660-11ea-3f2f-4d4963a60e39
function fases(f,p;rangos=[-1,1,-1,1],N=10)
	xmin,xmax,ymin,ymax=rangos
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax
	
	p1=plot()
	
	p1=flechas!(p1,f,p,rangos,N=N) #argumento opcional
	p1=trayectorias!(p1,f,p,rangos,N=N)
	
	xlims!(p1,(xmin, xmax))
	ylims!(p1,(ymin, ymax))
	return p1
end

# ╔═╡ 443e1200-0e16-11eb-14d1-27b178d0fffe
append!([1 2],

# ╔═╡ 90208f82-d4ba-11ea-11e6-7f51e6ffc5d5
fases(f_lineal,[k1,k2])

# ╔═╡ b51d4d90-d7f4-11ea-07aa-adcb4f418c99
md"## Clasificación sistema de segundo orden"

# ╔═╡ 0b2775be-d7fb-11ea-0b82-0d2b4d544d31
md"### Estables autovalores en el semiplano izquierdo"

# ╔═╡ 4100c3f0-d7f5-11ea-10c3-593316528d50
md"Nodo **estable** (polos reales en -1 y -0.5)"


# ╔═╡ d846eed0-d7f8-11ea-0596-d56c65bbbb5e
fases(f_lineal,[-0.5,-1.5])

# ╔═╡ f8cc9600-d7f8-11ea-138d-5f40babd5dce
md"Foco **estable** (polos reales en $$-1\pm i$$)"

# ╔═╡ ee673070-d7f9-11ea-3cb9-330579db1f09
fases(f_lineal,[-2,-2])

# ╔═╡ 34bbafa0-d7fb-11ea-2ed3-8b300f45fbac
md"### En el limite de la estabilidad (autovalores imaginarios puros)"

# ╔═╡ 345ec290-d7fb-11ea-1f16-d919242fefc3
md"Centro, polos en $$\pm i$$"

# ╔═╡ 86cc4ac0-d7fb-11ea-2101-4b12fafd41f7
fases(f_lineal,[-1,0])

# ╔═╡ 8c6fa2e0-dd75-11ea-163b-c70c276c3f2a
md"### Caso *degeneredo*, autovalor nulo"

# ╔═╡ ad48c140-dd75-11ea-03fe-75a0b19f91bd
md"""No hay un punto de equilibrio sino una "*línea de equlibrio*"
autovalores -1 y 0"""

# ╔═╡ 8bc272f0-dd75-11ea-3f24-1766ca4f5a7c
fases(f_lineal,[0,-1])

# ╔═╡ 33f181d0-d7fb-11ea-0a51-d3a11d1613df
md"### Inestables, algún autovalor en el semiplano derecho"

# ╔═╡ fa4a8370-d7f8-11ea-12b0-fde507f287be
md"Nodo **inestable** (polos reales en +1 y +0.5)"

# ╔═╡ dcb4f750-d7f8-11ea-0027-47491be6b131
fases(f_lineal,[-0.5,1.5])

# ╔═╡ 08c0cb70-d7fa-11ea-0efd-ab3b7f76a8b8
md"Foco **inestable** (polos reales en $$+1\pm i$$)"

# ╔═╡ 7dacbd90-d7fa-11ea-1063-e5741bf3f618
fases(f_lineal,[-2,2])

# ╔═╡ b6db1b70-d7fa-11ea-2bea-43295b3179a6
md"punto de silla (**Inestable**) polos en $$\pm1$$"

# ╔═╡ df4f1840-d7fa-11ea-1e89-3d6d103458d1
fases(f_lineal,[1,0])

# ╔═╡ fd823200-d7fc-11ea-1af8-c91386a02ab5
md"""## Escape en tiempo finito
poner explicación"""

# ╔═╡ 2647d140-d7fd-11ea-201a-5744621fbe12
escape(x,p,t)=x^2

# ╔═╡ 7c13c6b0-d7fd-11ea-0cc5-55f5ab6106fe
begin
	prob2 = ODEProblem(escape,1,(0,10));
	sol2 = solve(prob2,Tsit5(), reltol=1e-5, abstol=1e-5);
	fig=plot(sol2.t,min.(sol2.u,10000),xlim=(0,1))
	ylims!(fig,(0.0,10))
end

# ╔═╡ 82023190-dd77-11ea-12d1-cfa1a2496348
md"# Linealización"

# ╔═╡ 16d7ee1e-dd7a-11ea-13a7-197edb480651
md""" Ejemplo del generador de corriente

$$\dot x_1 = x_2$$

$$\dot x_2 = \frac{P_m - P_e sin(x_1)}{H}$$  """

# ╔═╡ ab5bd370-dd77-11ea-034f-6153487b9658
function generador(x,p,t)
    dx=[0.0;0.0]
    H=p[1]
    Pm=p[2] #u
    Pe=p[3]
    dx[1]=x[2]
    dx[2]=(1/H)*(Pm-Pe*sin(x[1]))
	return dx
end

# ╔═╡ 1668b190-dd7a-11ea-3c80-c5f082ec8647
md"""## Puntos de Equilibrio
Los puntos *interesantes* es donde las felchas se desvanecen, es decir donde la derivada se anula:

$$0= x_2$$

$$0 = \frac{P_m - P_e sin(x_1)}{H} \to sin(x_1)=\frac{P_m}{P_e}$$ 

Las soluciones existen si ${|P_m| \le |P_e|}$

soluciones tipo 1: $$(arcsin\left(\frac{P_m}{P_e}\right)+2n\pi,0)$$

soluciones tipo 2: $$(\pi-arcsin\left(\frac{P_m}{P_e}\right)+2n\pi,0)$$

"""

# ╔═╡ 351f5332-dd7c-11ea-2fc4-c9dc71684f23
@bind sliderMotor html"<input type=range min=-1.1 max=1.1 step=0.05>"

# ╔═╡ 5fd238f0-dd7b-11ea-2a00-a587ee59c261
begin
	plot(-2*pi:0.01:2*pi,sin)
	plot!([-2*pi,2*pi],[sliderMotor,sliderMotor ])
end

# ╔═╡ f102ad80-dd78-11ea-08b5-2d95e6f1a678
begin
	H=1
	Pe=1
	Pm=Pe*sliderMotor
	dibujo=fases(generador,[H,Pm,Pe],rangos=[-2pi,2pi,-2,2],N=15)
	if abs(Pm/Pe)<=1
		x1a=asin(Pm/Pe)
	    x1b=pi-asin(Pm/Pe)
		scatter!([x1a,x1b],[0,0]) #puntos de equilibrio
	end
    dibujo
end

# ╔═╡ db58def0-dd78-11ea-335b-4140d22d9a78
md""" En el punto de equilibrio f=0 así que usando taylor:
...


"""

# ╔═╡ 8ce3af70-de3b-11ea-17e1-e318972faa86
md""" Calculando la matriz jacobiana:

$$A=\begin{pmatrix} 
0                     &  1 \\ 
\frac{P_ecos(x^*)}{H}  &  0
\end{pmatrix}$$

$$B=\begin{pmatrix} 
1 \\ 
\frac{1}{H}
\end{pmatrix}$$
"""

# ╔═╡ e0f77820-de3c-11ea-2124-2de08feef36a
begin
	if abs(Pm/Pe)<=1
	   A1=[0 1;
		  Pe*cos(x1a) 0]
	   autov1,_1=eigen(A1)
	end
	md" El primer PE es inestable (punto de silla), los autovalores son:
	$(autov1[1]) y $(autov1[2])"
end

# ╔═╡ d252c5be-de3e-11ea-2f1d-1dbc809c9980
fases((x,p,t)->A1*x,[])

# ╔═╡ 46dd7a80-de3e-11ea-0e9d-abe7496efef7
begin
	if abs(Pm/Pe)<=1
	   A2=[0 1;
		  Pe*cos(x1b) 0]
	   autov2,_2=eigen(A2)
	end
	md" El segundo PE es un centro en la aproximación lineal, los autovalores son:
	$(autov2[1]) y $(autov2[2])"
end

# ╔═╡ aa716bfe-de3e-11ea-0dc6-d59275895621
fases((x,p,t)->A2*x,[],N=10)

# ╔═╡ 33af7e82-d7f5-11ea-12bb-5982575f35d3
md"""## Ciclo límite
Ecuación de Van der Pol
"""

# ╔═╡ 3d70984e-d7f5-11ea-3243-9d3932ee155b
function van_der_pol(x,p,t)
    dx=[0.0;0.0]
    m=p[1]
    c=p[2]
    k=p[3]
    dx[1]=x[2]
    dx[2]=(1/m)*(-k*x[1]+c*(1-x[1]^2)*x[2])
	return dx
end

# ╔═╡ afa04610-d7f4-11ea-1cd0-e3a91d91c5c9
fases(van_der_pol,[1.0,1.0,1.0],rangos=[-3,3,-3,3])

# ╔═╡ Cell order:
# ╟─6ba02e80-d4cf-11ea-19cb-33028d3dc67d
# ╟─7db3fdb0-d4d2-11ea-21fc-597c55dda572
# ╟─76079a70-d4cf-11ea-2645-8f43c78eae52
# ╠═66e9ffd0-d4d7-11ea-09c0-03ebaacab879
# ╠═e885f5b0-d4d4-11ea-06d6-9ddf9619961f
# ╠═b99fc9b0-d4d9-11ea-0e18-cfd5af90188f
# ╠═3be6b020-d4d8-11ea-3356-23afd148ebc9
# ╠═220f6f12-d4d9-11ea-159e-c1e1551fd732
# ╟─23f5ece0-d4da-11ea-1330-c5d077bd3a89
# ╠═d8b5ed40-d4d7-11ea-27d4-bf594518802b
# ╠═1be7c0c0-d4ba-11ea-1e43-ef66c35b0830
# ╠═2b7df180-d4d8-11ea-118f-75e8709e36c3
# ╟─455e0c30-d4d7-11ea-22ce-0d04b4e8d860
# ╠═159f28e2-d4db-11ea-13e6-a97e2e7f7c60
# ╠═04501380-d63c-11ea-2ed8-e11f02de80a4
# ╠═c90980d0-d637-11ea-112d-c336b2080d87
# ╠═4a7f5f72-0e0c-11eb-2b7c-4d7214c09c4b
# ╠═51059990-0e0c-11eb-2d55-89ac739619c7
# ╠═a92de300-d63c-11ea-2446-5b099592d734
# ╠═09980bd0-d63d-11ea-0752-1d033ee00cc9
# ╠═e6afc70e-d63d-11ea-0d36-738b13f0de30
# ╠═6cbb32f0-d665-11ea-3f24-41cbe1910b80
# ╠═22398f60-d665-11ea-3352-714b467152b3
# ╠═5a802550-d665-11ea-3f9f-b103ce80822f
# ╟─57ffd1f0-d637-11ea-1a81-7dd574291417
# ╟─5b26b050-d7eb-11ea-2488-b7da296fed24
# ╟─84c3c3d0-d7eb-11ea-3f35-29559df926c3
# ╟─a0573fe0-d6a2-11ea-057e-5140c323c005
# ╠═1790c760-0e0d-11eb-1c47-f5e902e158fb
# ╠═f5ca89e0-0e0c-11eb-32a1-fb22511f42ae
# ╠═915aa560-0e0e-11eb-18d4-e3a899ee2c1c
# ╠═165464e2-0e0f-11eb-3851-3701a6c78502
# ╠═7115307e-d660-11ea-3f2f-4d4963a60e39
# ╠═443e1200-0e16-11eb-14d1-27b178d0fffe
# ╠═90208f82-d4ba-11ea-11e6-7f51e6ffc5d5
# ╟─b51d4d90-d7f4-11ea-07aa-adcb4f418c99
# ╟─0b2775be-d7fb-11ea-0b82-0d2b4d544d31
# ╟─4100c3f0-d7f5-11ea-10c3-593316528d50
# ╠═d846eed0-d7f8-11ea-0596-d56c65bbbb5e
# ╟─f8cc9600-d7f8-11ea-138d-5f40babd5dce
# ╠═ee673070-d7f9-11ea-3cb9-330579db1f09
# ╟─34bbafa0-d7fb-11ea-2ed3-8b300f45fbac
# ╟─345ec290-d7fb-11ea-1f16-d919242fefc3
# ╠═86cc4ac0-d7fb-11ea-2101-4b12fafd41f7
# ╟─8c6fa2e0-dd75-11ea-163b-c70c276c3f2a
# ╟─ad48c140-dd75-11ea-03fe-75a0b19f91bd
# ╠═8bc272f0-dd75-11ea-3f24-1766ca4f5a7c
# ╟─33f181d0-d7fb-11ea-0a51-d3a11d1613df
# ╟─fa4a8370-d7f8-11ea-12b0-fde507f287be
# ╠═dcb4f750-d7f8-11ea-0027-47491be6b131
# ╠═08c0cb70-d7fa-11ea-0efd-ab3b7f76a8b8
# ╠═7dacbd90-d7fa-11ea-1063-e5741bf3f618
# ╟─b6db1b70-d7fa-11ea-2bea-43295b3179a6
# ╠═df4f1840-d7fa-11ea-1e89-3d6d103458d1
# ╟─fd823200-d7fc-11ea-1af8-c91386a02ab5
# ╠═2647d140-d7fd-11ea-201a-5744621fbe12
# ╠═7c13c6b0-d7fd-11ea-0cc5-55f5ab6106fe
# ╟─82023190-dd77-11ea-12d1-cfa1a2496348
# ╟─16d7ee1e-dd7a-11ea-13a7-197edb480651
# ╠═ab5bd370-dd77-11ea-034f-6153487b9658
# ╟─1668b190-dd7a-11ea-3c80-c5f082ec8647
# ╠═5fd238f0-dd7b-11ea-2a00-a587ee59c261
# ╟─351f5332-dd7c-11ea-2fc4-c9dc71684f23
# ╠═f102ad80-dd78-11ea-08b5-2d95e6f1a678
# ╟─db58def0-dd78-11ea-335b-4140d22d9a78
# ╟─8ce3af70-de3b-11ea-17e1-e318972faa86
# ╟─e0f77820-de3c-11ea-2124-2de08feef36a
# ╠═d252c5be-de3e-11ea-2f1d-1dbc809c9980
# ╟─46dd7a80-de3e-11ea-0e9d-abe7496efef7
# ╠═aa716bfe-de3e-11ea-0dc6-d59275895621
# ╟─33af7e82-d7f5-11ea-12bb-5982575f35d3
# ╠═3d70984e-d7f5-11ea-3243-9d3932ee155b
# ╠═afa04610-d7f4-11ea-1cd0-e3a91d91c5c9
