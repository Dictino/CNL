### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ 88629870-3490-11eb-2215-d765d59a440c
using PlutoUI

# ╔═╡ 1be7c0c0-d4ba-11ea-1e43-ef66c35b0830
using Plots

# ╔═╡ 159f28e2-d4db-11ea-13e6-a97e2e7f7c60
#using DifferentialEquations
using OrdinaryDiffEq

# ╔═╡ 6cbb32f0-d665-11ea-3f24-41cbe1910b80
using LinearAlgebra

# ╔═╡ 6ba02e80-d4cf-11ea-19cb-33028d3dc67d
md"# _Primer tema_
+ ¿Qué es una ecuación diferencial y como se resuelve?
    * Euler
    * Métodos mejores
+ Sistema de segundo orden con dos polos
+ Explosión en tiempo finito
+ Linealización
+ Ejemplos
    * Generador
+ Ciclo límite
    * Van der Pool
+ Sistemas discretos
    * Ecuación logística, el camino hacia el Caos"

# ╔═╡ 7db3fdb0-d4d2-11ea-21fc-597c55dda572
md"# ¿Qué es una ecuación diferencial?

En la asignatura nos interesan sistemas como este:

$$\frac{dx}{dt}=f(x,t)$$

$$x(0)=x_0$$

Donde x es un vector y la derivada respecto al tiempo se suele escribir como $$\dot x$$ 

Las ecuaciones diferenciales son útiles cuando es fácil expresar cómo cambia una magnitud en el tiempo, pero no tanto cual es el estado $x(t)$ directamente, que es lo que ocurre en casi todos los sistemas físicos.
"

# ╔═╡ 76079a70-d4cf-11ea-2645-8f43c78eae52
md" ## Euler
Si la función $$f$$ tiene buen comportamiento (lo veremos más adelante cuando repasemos la existencia y unicidad) hay una solución de la ecuación que es continua y derivable y usando taylor:

$$x(\Delta t)=x(0) + f(x(0),0)\Delta t + o(\Delta t^2)$$

Si se toman pasos pequeños y repetimos el procedimiento:

$$x(2\Delta t)=x(\Delta t) + f(x(\Delta t),\Delta t) \Delta t + o(\Delta t^2)$$
...

y llamando $$t_n$$ a $$n\Delta t$$ y $$x_n$$ a $$x(t_n)$$

$$x_{n+1} \approx x_n + f(x_n,t_n)\Delta t + o(\Delta t^2)$$

### Ejemplo de implementación sencilla
Para el caso escalar

"

# ╔═╡ 66e9ffd0-d4d7-11ea-09c0-03ebaacab879
function euler(f,x0,tf,Δt)
	t=0:Δt:tf
	x=zeros(size(t))
	x[1]=x0
	for n=1:(length(t)-1)
		x[n+1]=x[n]+f(x[n],t[n])*Δt
	end
	return (x,t)
end

# ╔═╡ 0b94f1c0-7113-11eb-13da-c5c5b95ebbf5
md"Vamos a ver si funciona resolviendo
$\dot x=-x+sin(t)$ "

# ╔═╡ e885f5b0-d4d4-11ea-06d6-9ddf9619961f
begin
	f(x,t)=-x+sin(t)
		
	f(x,p,t)=f(x,t) #¡Despacho múltiple! (lo necesitamos para luego)
end

# ╔═╡ b99fc9b0-d4d9-11ea-0e18-cfd5af90188f
tf=10

# ╔═╡ 3be6b020-d4d8-11ea-3356-23afd148ebc9
@bind N Slider(1:100)

# ╔═╡ 220f6f12-d4d9-11ea-159e-c1e1551fd732
Δt=tf/N;

# ╔═╡ 23f5ece0-d4da-11ea-1330-c5d077bd3a89
md"N=$(N) $$\to \Delta t=$$ $(Δt)"

# ╔═╡ d8b5ed40-d4d7-11ea-27d4-bf594518802b
x,t =euler(f,0,10,Δt);

# ╔═╡ 2b7df180-d4d8-11ea-118f-75e8709e36c3
plot(t,x,label="x")

# ╔═╡ 455e0c30-d4d7-11ea-22ce-0d04b4e8d860
md"## Hay métodos mejores:

+ Ajustan los pasos automáticamente

+ Tienen más precisión para los mismos pasos

+ Están integrados con otras funciones como *plots*, permite interpolar los resultados automáticamente, etc...

Vamos a usar la librería OrdinaryDiffEq que es parte de DifferencialEquation.jl
"

# ╔═╡ c90980d0-d637-11ea-112d-c336b2080d87
begin
	x0 = 0
	tspan = (0.0,10.0)
	#ojo se necesita f(x,p,t) donde p son los parámetros (no tenemos)
	#de ahí la definición de antes
	prob = ODEProblem(f,x0,tspan);
	sol = solve(prob,Tsit5(), reltol=1e-8, abstol=1e-8);
	sol.u #su u es nuestro x
end

# ╔═╡ 1f44a88e-7114-11eb-0ac2-d12cd3fb40ef
sol.t

# ╔═╡ 4a7f5f72-0e0c-11eb-2b7c-4d7214c09c4b
sol(2.751) #¡se puede usar como una función continua de t!

# ╔═╡ 51059990-0e0c-11eb-2d55-89ac739619c7


# ╔═╡ a92de300-d63c-11ea-2446-5b099592d734
plot(sol, xaxis="t",yaxis="x(t)",label="x",legend=false)

# ╔═╡ 09980bd0-d63d-11ea-0752-1d033ee00cc9
md"## Sistema de segundo orden con dos polos

En un sistema lineal

$\dot x=Ax$

La solución es $x=x_0e^{At}$ y su carácter depende de los autovalores, es decir las soluciones de $|A- \lambda \mathbb{1} |=0$

La estabilidad depende de que la parte real de los $\lambda$ sea negativa (los autovalores estén en el semiplano izquierdo del plano complejo)

Nos interesa el sistema lineal segundo orden porque es simple y contiene todos los casos de interés:

$\dot x_1=x_2$
$\dot x_2=p_1x_1+p_2x_2$
"

# ╔═╡ e6afc70e-d63d-11ea-0d36-738b13f0de30
function f_lineal(x,p,t)
	dx=[0.0 ; 0.0]
	dx[1]=x[2];
	dx[2]=p[1]*x[1]+p[2]*x[2]
	return dx
end

# ╔═╡ 2ff44950-7116-11eb-3ad9-c3e145cdae4e
md"
Polos Reales? $(@bind son_reales CheckBox())
$(@bind slider1 Slider(-1:0.05:1))
$(@bind slider2 Slider(-1:0.05:1))
"

# ╔═╡ a0573fe0-d6a2-11ea-057e-5140c323c005
begin
	if son_reales==true 
		# el polinomio es (s-real1)*(s-real2)
		# s^2+(-real1-real2)s + real1*real2 --> k2= real1 + real2 k1=real1*real2
		k2= slider1 + slider2
		k1=-slider1*slider2
	else # el polinomio es (s-real-imag*im)*(s-real+imag*im)
		 # es decir s^2 - 2*real*s + real^2 + imag^2 --> k2=-2*real k1=-real^2-imag^2
		k2=2*slider1
		k1=-slider1^2 - slider2^2
	end
	A=[0 1; k1 k2]
	autovalores,_=eigen(A)  #para esto necesito LinearAlgebra
	md"""
	k1= $(k1), k2= $(k2)
	
	autovalores= $(autovalores)"""
end

# ╔═╡ 5a802550-d665-11ea-3f9f-b103ce80822f
scatter(real(autovalores),imag(autovalores),
	grid=true,legend=false, framestyle = :origin,
	xaxis=[-1,1], yaxis=[-1,1] )

# ╔═╡ 1790c760-0e0d-11eb-1c47-f5e902e158fb
md" # Muy bien pero ¡Enséñame el código!
... vosotros lo habéis querido ☺	"

# ╔═╡ f5ca89e0-0e0c-11eb-32a1-fb22511f42ae
function flechas!(figura,f,p,rangos;N=10)
	#modifica una figura añadiendo las flechas por eso tiene una exclamación !
	# N es un parámetro opcional, si no se lo paso el valor es 10
	# los otros son posicionales
	xmin,xmax,ymin,ymax=rangos #desestructurar un iterable (explicar)
	longitud=max(abs(xmax-xmin),abs(ymax-ymin))/N
	
    #hago una maya de x e ys y pongo las flechas en cada punto
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax

    #barrido para las flechas
	X1=[];
	X2=[];
	U=[];
	V=[];
	for x1 in x1s #barrido para calcular las flechas
		for x2 in x2s
			push!(X1,x1); #Push añade un componente al final de un vector
			push!(X2,x2);
			F=f([x1,x2],p,0)
			X=0.1*F/(0.1+sqrt(F[1]^2+F[2]^2)) #Saturo el tamaño de la flecha
			push!(U,X[1]);
			push!(V,X[2]);
		end
	end
	#Pinta las flechas
	quiver!(figura,X1,X2,quiver=(U,V), arrow=arrow(0.1*longitud, 0.1*longitud))
	return figura
end

# ╔═╡ 165464e2-0e0f-11eb-3851-3701a6c78502
#= Es un euler modificado que integra hacia delante y hacia atrás
    congtinúa integrando hasta que se sale o completa un número de pasos. 
	   -hacia delante se detectan bien los PE y ciclos límite estables
	   -Hacia atrás permite ver también los inestables
	Ajusta dt para que el paso tenga un tamaño que sea aproximadamente constante=#

function trayectoria(x0,f,p,rangos;N=10,pasosMax=100)
	xmin,xmax,ymin,ymax=rangos
	long=max(abs(xmax-xmin),abs(ymax-ymin))/N
	x,y=x0
	pasos=0
	
	X1=[x]
	Y1=[y]
	X2=[x] #parte de la trayectoria hacia atrás en el tiempo
	Y2=[y]
	
	while (xmin<x<xmax) & (ymin<y<ymax) & (pasos<pasosMax)
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
	while (xmin<x<xmax) & (ymin<y<ymax) & (pasos<pasosMax)
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

# ╔═╡ b51d4d90-d7f4-11ea-07aa-adcb4f418c99
md"## Clasificación sistema de segundo orden"

# ╔═╡ 0b2775be-d7fb-11ea-0b82-0d2b4d544d31
md"### Estables autovalores en el semiplano izquierdo"

# ╔═╡ 4100c3f0-d7f5-11ea-10c3-593316528d50
md"Nodo **estable** (polos reales en -1 y -0.5)"


# ╔═╡ f8cc9600-d7f8-11ea-138d-5f40babd5dce
md"Foco **estable** (polos reales en $$-1\pm i$$)"

# ╔═╡ 34bbafa0-d7fb-11ea-2ed3-8b300f45fbac
md"### En el limite de la estabilidad (autovalores imaginarios puros)"

# ╔═╡ 345ec290-d7fb-11ea-1f16-d919242fefc3
md"Centro, polos en $$\pm i$$"

# ╔═╡ 8c6fa2e0-dd75-11ea-163b-c70c276c3f2a
md"### Caso *degenerado*, autovalor nulo"

# ╔═╡ ad48c140-dd75-11ea-03fe-75a0b19f91bd
md"""No hay un punto de equilibrio sino una "*línea de equilibrio*"
autovalores -1 y 0"""

# ╔═╡ 33f181d0-d7fb-11ea-0a51-d3a11d1613df
md"### Inestables, algún autovalor en el semiplano derecho"

# ╔═╡ fa4a8370-d7f8-11ea-12b0-fde507f287be
md"Nodo **inestable** (polos reales en +1 y +0.5)"

# ╔═╡ 08c0cb70-d7fa-11ea-0efd-ab3b7f76a8b8
md"Foco **inestable** (polos reales en $$+1\pm i$$)"

# ╔═╡ b6db1b70-d7fa-11ea-2bea-43295b3179a6
md"punto de silla (**Inestable**) polos en $$\pm1$$"

# ╔═╡ fd823200-d7fc-11ea-1af8-c91386a02ab5
md"""## Escape en tiempo finito
Consideremos la ecuación aparentemente *inocente*:

$\dot x=x^2$

$x_0=1$

Es separable y como se vio en teoría la solución es obviamente (basta con derivar y sustituir)

$x(t)=\frac{1}{1-t}$

Que en t=1 se hace infinita, veámoslo

"""

# ╔═╡ 2647d140-d7fd-11ea-201a-5744621fbe12
escape(x,p,t)=x^2

# ╔═╡ 7c13c6b0-d7fd-11ea-0cc5-55f5ab6106fe
begin
	prob2 = ODEProblem(escape,1,(0,2));
	sol2 = solve(prob2,Tsit5(), reltol=1e-5, abstol=1e-5);
	fig=plot(sol2.t,min.(sol2.u,10000),xlim=(0,1))
	#saturo a 10000 para evitar que "pete"
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
	x1,x2=x #así es más bonito, es la forma RECOMENDADA para la legibilidad
	H,Pm,Pe=p #Pm es u de la teoría
    dx1=x2
    dx2=(1/H)*(Pm-Pe*sin(x1))
	#el return es opcional (si no se pone devuelve la última línea
	#pero me gusta ser explícito
	return [dx1; dx2] 
end

# ╔═╡ 1668b190-dd7a-11ea-3c80-c5f082ec8647
md"""## Puntos de Equilibrio
Los puntos *interesantes* es donde las flechas se **desvanecen**, es decir donde la derivada se anula:

$$0= x_2$$

$$0 = \frac{P_m - P_e sin(x_1)}{H} \to sin(x_1)=\frac{P_m}{P_e}$$ 

Las soluciones existen si ${|P_m| \le |P_e|}$

soluciones tipo 1: $$(arcsin\left(\frac{P_m}{P_e}\right)+2n\pi,0)$$

soluciones tipo 2: $$(\pi-arcsin\left(\frac{P_m}{P_e}\right)+2n\pi,0)$$

"""

# ╔═╡ 351f5332-dd7c-11ea-2fc4-c9dc71684f23
@bind sliderMotor Slider(-1.1:0.05:1.1)

# ╔═╡ 5fd238f0-dd7b-11ea-2a00-a587ee59c261
begin
	plot(-2*pi:0.01:2*pi,sin,legend=false)
	plot!([-2*pi,2*pi],[sliderMotor,sliderMotor ])
	title!("Pm/Pe=$sliderMotor")
end

# ╔═╡ db58def0-dd78-11ea-335b-4140d22d9a78
md""" En el punto de equilibrio $f(x^*)=0$ así que usando taylor:

definiendo $\delta x=x-x^*$

$\delta \dot x= \dot x=f(x^*+\delta x) \approx f(x^*)+\nabla f(x^*)\delta x=A\delta x$ 

donde A es la matriz jacobiana

$$A=\begin{pmatrix} 
\frac{df_1}{dx_1}  & ... & \frac{df1}{dx_n} \\
\vdots             & \ddots   &  \vdots \\
\frac{df_n}{dx_1}  & ... &\frac{df_n}{dx_n}
\end{pmatrix}$$

que en el caso del motor es:

$$A=\begin{pmatrix} 
0                     &  1 \\ 
\frac{P_ecos(x^*)}{H}  &  0
\end{pmatrix}$$

"""

# ╔═╡ 33af7e82-d7f5-11ea-12bb-5982575f35d3
md"""## Ciclo límite
Ecuación de Van der Pol

$\dot x_1=x_2$

$\dot x_2=\frac{-kx_1+c(1-x_1^2)x_2}{m}$

"""

# ╔═╡ 3d70984e-d7f5-11ea-3243-9d3932ee155b
function van_der_pol(x,p,t)
    m,c,k=p
	x1,x2=x
    dx1=x2
    dx2=(1/m)*(-k*x1+c*(1-x1^2)*x2)
	return [dx1;dx2]
end

# ╔═╡ 713de9b0-3490-11eb-041b-c7083c829ba7
md"""# El camino hacia el Caos
Presentado por el Dr. Chaos ☺

Ecuación logística

$y_{n+1}=ry_n(1-y_n)$
$r \in [0, 4]$
$y \in [0, 1]$

"""

# ╔═╡ 98165ef0-3490-11eb-0e57-259fcf7361c2
logistica(y,r)=r*y*(1-y)

# ╔═╡ 10625a80-7122-11eb-3601-273ea4390a23
md"Lo bueno de las ecuaciones en diferencias es que es trivial resolverlas, basta con iterar"

# ╔═╡ ea63b4f0-3490-11eb-26a8-6f5e1f2d3c05
function trayectoria(y_inicial,r,pasos)
	y=zeros(pasos)
	y[1]=y_inicial
	for i=1:(pasos-1)
		y[i+1]=logistica(y[i],r)
	end
	return y
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

# ╔═╡ 90208f82-d4ba-11ea-11e6-7f51e6ffc5d5
fases(f_lineal,[k1,k2])

# ╔═╡ d846eed0-d7f8-11ea-0596-d56c65bbbb5e
fases(f_lineal,[-0.5,-1.5])

# ╔═╡ ee673070-d7f9-11ea-3cb9-330579db1f09
fases(f_lineal,[-2,-2])

# ╔═╡ 86cc4ac0-d7fb-11ea-2101-4b12fafd41f7
fases(f_lineal,[-1,0])

# ╔═╡ 8bc272f0-dd75-11ea-3f24-1766ca4f5a7c
fases(f_lineal,[0,-1])

# ╔═╡ dcb4f750-d7f8-11ea-0027-47491be6b131
fases(f_lineal,[-0.5,1.5])

# ╔═╡ 7dacbd90-d7fa-11ea-1063-e5741bf3f618
fases(f_lineal,[-2,2])

# ╔═╡ df4f1840-d7fa-11ea-1e89-3d6d103458d1
fases(f_lineal,[1,0])

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

# ╔═╡ 46dd7a80-de3e-11ea-0e9d-abe7496efef7
begin
	if abs(Pm/Pe)<=1
	   A1=[0 1;
		  Pe*cos(x1a) 0]
	   autov1,_1=eigen(A1)
	   fases((x,p,t)->A1*x,[]) #otra forma e hacerlo con una función anónima
	   title!("""El primer PE es inestable (silla),
		λ1=$(autov1[1]) λ2=$(autov1[2])""")
	end
end

# ╔═╡ 44e6bce2-800b-11eb-3651-e3ba5df8617e
begin
	if abs(Pm/Pe)<=1
	   A2=[0 1;
		  Pe*cos(x1b) 0]
	   autov2,_2=eigen(A2)
	   fases((x,p,t)->A2*x,[]) #otra forma e hacerlo con una función anónima
       title!("""El primer PE es un centro (en la aproximación lineal),
	λ1=$(autov2[1]) λ2=$(autov2[2])""")
	end
end

# ╔═╡ afa04610-d7f4-11ea-1cd0-e3a91d91c5c9
fases(van_der_pol,[1.0,1.0,1.0],rangos=[-3,3,-3,3])

# ╔═╡ b8326ede-3490-11eb-1e0b-530752ed836e
@bind r s=Slider(0:0.05:4,default=2)

# ╔═╡ a07f3bc0-3490-11eb-11b3-6bedb92502cb
let  #introduce un "scope" nuevo no como begin-end lo que permite usar las mismas variables sin tener que inventar nombres nuevos
	x=0:0.01:1
	fig=plot(x,logistica.(x,r),xlim=(0,1),ylim=(0,1),
		label=:none, title="r=$r",ratio=1)
end

# ╔═╡ d43266e0-3490-11eb-0d9a-59ef7a663b50
@bind ver CheckBox(default=false)

# ╔═╡ 165343ee-3491-11eb-10a7-7b4d6fef56f4
function puntos(R)
	px=[]
	py=[]
	for r=0:0.02:R
		yt=trayectoria(0.5,r,100)
		k=length(yt)
		y_est=yt[(k÷2):k]
		append!(px,r*ones(length(y_est))) #append es como push pero añade un vector
		append!(py,y_est)
	end
	return (px,py)
end

# ╔═╡ e042eeee-3490-11eb-1406-058ee6753c1f
let
	y=trayectoria(0.5,r,50)
	p1=plot(y,title="r=$r",marker=:dot,markersize=:2,label=:none,ylimits=(0,1))
	k=length(y)
		p2=scatter()
	if ver
		px,py=puntos(4)
		p2=scatter!(p2,px,py,legend=:none,color=:grey, alpha=0.1,
		            linecolor=:grey, linealpha=0.1,marker=1)
	end
	y_estacionarios=y[(k÷2):k] #Nota ver cómo se escribe esto
	p2=scatter!(p2,r*ones(length(y_estacionarios)),y_estacionarios,
		        xlim=(0,4),ylim=(0,1),legend=:none,color=:red,marker=2)
	l = @layout [a ; b]  
	plot(p1,p2, layout=l, size=(650,500))
end

# ╔═╡ 315a35f0-3491-11eb-0ce1-636a82ef3dd9
md"""# ¿Cómo se llega al caos?

Haciendo lo mismo una y otra vez...

Primer paso, periodo 1:"""

# ╔═╡ dd340bd0-3491-11eb-31ab-5b1bf6286cf0
@bind Niter Slider(1:50,default=1,show_value=true)

# ╔═╡ 3765e8e2-3491-11eb-00a5-274f3caf6a15
@bind y0 Slider(0:0.05:1,default=0.5,show_value=true)

# ╔═╡ de04c952-3491-11eb-276b-b530d61da953
@bind r2 Slider(0:0.05:4,default=2,show_value=true)

# ╔═╡ 39b45040-3492-11eb-1e48-998b1e4b1911
function iterar(f,N)
	x=0:0.001:1
	fig2=plot(x,f.(x),xlim=(0,1),ylim=(0,1),label=:none, title="r=$r2")
	plot!(fig2,[0,1],[0,1],color=:red,label=:none, ratio=1)
	
	x=y0 #punto en el plano de la figura
	y=0
	for i=1:N
		#levanto y hasta f(x)
		xn=x
		yn=f(x)
		plot!(fig2,[x,xn],[y,yn],color=:green,label=:none)
		#desplazo x hasta y manteniendo yn
		x=yn
		y=yn 
		plot!(fig2,[xn,x],[yn,y],color=:green,label=:none)
	end
	fig2
end

# ╔═╡ 2e058370-729c-11eb-0822-793d3db17e74
begin
    f1(x)=logistica(x,r2)
	f2(x)=f1(f1(x))
	f4(x)=f2(f2(x))
	f8(x)=f4(f4(x))
	f16(x)=f8(f8(x))
	iterar(f1,Niter) # prueba con f1 f2 ...
end

# ╔═╡ 4c657840-3492-11eb-093d-59a4b3854b7a
md"""En un sistema discreto la estabilidad depende de la derivada, es parecido a los sistemas continuos.
Un punto de equilibrio $x^*$ es aquel en el que $x_{n+1}=x_n$, es decir

$x_{n+1}=f(x_n)=x_n=x^*$

Es decir que es un punto fijo de $f$, es decir $f(x^*)=x^*$

Si estamos cerca de un punto de equilibrio $x_n=x^*+\delta x_n$ luego haciendo taylor:

$f(x^*+\delta x_n)=f(x^*) + \nabla f(x^*)\delta x_n + o(\delta x_n^2)$


Pero entonces:

$\delta x_{n+1}=x_{n+1}-x^*=f(x^*)-x^* + \nabla f(x^*)\delta x_n + o(\delta x_n^2) \approx \nabla f(x^*) \delta x_n=A \delta x_n$

Cuya solución es $\delta x_n=A^n\delta x_0$ que es estable si los autovalores de $A$ están en el círculo unitario. Es parecido a los sistemas continuos cambiando el eje imaginario por el **círculo unidad**.

En el caso de una variable se reduce a que $|f'(x^*)|<1$ como se intuye gráficamente

"""

# ╔═╡ 55e2ab40-3492-11eb-2953-e9d8bd17735c
md"En cada paso el cambio de $r$ necesario para producir una bifurcación es más pequeño.

De hecho la bifurcación de orden 4 es como una versión pequeña de la de orden 2 como hemos visto, si llamamos $r_n$ al $r$ que produce una bifurcación tenemos que

$\frac{r_{n+2}-r_{n+1}}{r_{n+1}-r_{n}} \to 4.6992...$ que es la constante de Feigenbaum y que es **la misma** para muchos otros procesos esto se conoce como **Universalidad**.

Como el paso sigue aproximadamente una serie geométrica la suma de infinitas bifurcaciones se produce con un cambio finito de $r$ ya que:

$r_\infty-r_0 \approx \Delta r_0 +  \frac {\Delta r_0}{4.7} +\frac {\Delta r_0}{4.7^2}+\frac {\Delta r_0}{4.7^3}...=\Delta r_0 \frac {1}{1-1/4.7}\approx 1.27\Delta r_0$

" 

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
OrdinaryDiffEq = "~6.6.5"
Plots = "~1.25.7"
PlutoUI = "~0.7.32"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "5e98d6a6aa92e5758c4d58501b7bf23732699fa3"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.2"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "ba19d1c8ff6b9c680015033c66802dd817a9cf39"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.7"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f9982ef575e19b0e5c7a98c6e75ee496c0f73a93"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.0"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "7b8f09d58294dc8aa13d91a8544b37c8a1dcbc06"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.4"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "6b6f04f93710c71550ec7e16b650c1b9a612d0b6"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.16.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "31186e61936fbbccb41d809ad4338c9f7addf7ae"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.0"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "d75333ab19d6d01c53bb350a9aabb074ba768a9d"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.81.3"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "84083a5136b6abf426174a58325ffd159dd6d94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.1"

[[Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "24d26ca2197c158304ab2329af074fbe14c988e4"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.45"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "3e1289d9a6a54791c1ee60da0850f4fd71188da6"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.11.0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "0f8ef5dcb040dbb9edd98b1763ac10882ee1ff03"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.12"

[[FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "6eae72e9943d8992d14359c32aed5f892bda1569"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.10.0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "4a740db447aae0fbeb3ee730de1afbb14ac798a1"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.63.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "aa22e1ee9e722f1da183eb33370df4c1aeb6c2cd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "d727758173afef0af878b29ac364a0eca299fc6b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.5.1"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "3965a3216446a6b020f0d48f1ba94ef9ec01720d"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.6"

[[Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d8bccde6fc8300703673ef9e1383b11403ac1313"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.0+0"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "1ed18ccf6292d89abf85beba35b9399aeddff9b2"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.2.3"

[[Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "e60270d7871e7ffe66b3a90b477ecb5df037aa0c"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.7.11"

[[KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "0328ad9966ae29ccefb4e1b9bfd8c8867e4360df"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.3"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "3696fdc1d3ef6e4d19551c92626066702a5db91c"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.7.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "6dd77ee76188b0365f7d882d674b95796076fa2c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.5"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "a050cd5581a204eeda3ad13c1d2aabdc3c451b4e"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.11.1"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "67c0dfeae307972b50009ce220aae5684ea852d1"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.101"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[ManualMemory]]
git-tree-sha1 = "9cb207b18148b2199db259adfa923b45593fe08e"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "b61c51cd5b9d8b197dfcbbf2077a0a4e1505278d"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.14"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "1475c25a9dc4de848a5234543f2cb8601ada67d4"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.6.5"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "7e4920a7d4323b8ffc3db184580598450bde8a8e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.7"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "ae6145ca68947569058866e443df69587acc1806"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.32"

[[Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "55f5db122f19d8b5b26fe9576edc1ff819e499bb"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.3"

[[PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "0bc9e1a21ba066335a5207ac031ee41f72615181"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.3"

[[PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "e4cb8d4a2edf9b3804c1fb2c2de57d634ff3f36e"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.2.3"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "37c1631cb3cc36a535105e6d5557864c82cd8c2b"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.0"

[[RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "5144e1eafb2ecc75765888a4bdcd3a30a6a08b14"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.24.1"

[[RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "7ad4c2ef15b7aecd767b3921c0d255d39b3603ea"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.9"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "1410aad1c6b35862573c01b96cd1f6dbe3979994"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.28"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "f4862c0cb4e34ed182718221028ba1bf50742108"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.26.1"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "75c89362201983c500dd34923b015dbecdae7a90"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.0"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e6bf188613555c78062842777b116905a9f9dd49"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.0"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "a635a9333989a094bddc9f940c04c549cd66afcf"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.4"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f35e1879a71cca95f4826a14cdbf0b9e253ed918"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.15"

[[StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "12cf3253ebd8e2a3214ae171fbfe51e7e8d8ad28"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.9"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "884539ba8c4584a3a8173cb4ee7b61049955b79c"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.4.7"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "c3ab8b77b82fd92e2b6eea8a275a794d5a6e4011"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.9"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "e9a35d501b24c127af57ca5228bcfb806eda7507"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.24"

[[VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─6ba02e80-d4cf-11ea-19cb-33028d3dc67d
# ╟─7db3fdb0-d4d2-11ea-21fc-597c55dda572
# ╟─76079a70-d4cf-11ea-2645-8f43c78eae52
# ╠═66e9ffd0-d4d7-11ea-09c0-03ebaacab879
# ╟─0b94f1c0-7113-11eb-13da-c5c5b95ebbf5
# ╠═e885f5b0-d4d4-11ea-06d6-9ddf9619961f
# ╠═b99fc9b0-d4d9-11ea-0e18-cfd5af90188f
# ╠═88629870-3490-11eb-2215-d765d59a440c
# ╠═3be6b020-d4d8-11ea-3356-23afd148ebc9
# ╠═220f6f12-d4d9-11ea-159e-c1e1551fd732
# ╟─23f5ece0-d4da-11ea-1330-c5d077bd3a89
# ╠═d8b5ed40-d4d7-11ea-27d4-bf594518802b
# ╠═1be7c0c0-d4ba-11ea-1e43-ef66c35b0830
# ╠═2b7df180-d4d8-11ea-118f-75e8709e36c3
# ╟─455e0c30-d4d7-11ea-22ce-0d04b4e8d860
# ╠═159f28e2-d4db-11ea-13e6-a97e2e7f7c60
# ╠═c90980d0-d637-11ea-112d-c336b2080d87
# ╠═1f44a88e-7114-11eb-0ac2-d12cd3fb40ef
# ╠═4a7f5f72-0e0c-11eb-2b7c-4d7214c09c4b
# ╠═51059990-0e0c-11eb-2d55-89ac739619c7
# ╠═a92de300-d63c-11ea-2446-5b099592d734
# ╟─09980bd0-d63d-11ea-0752-1d033ee00cc9
# ╠═e6afc70e-d63d-11ea-0d36-738b13f0de30
# ╠═6cbb32f0-d665-11ea-3f24-41cbe1910b80
# ╟─5a802550-d665-11ea-3f9f-b103ce80822f
# ╟─2ff44950-7116-11eb-3ad9-c3e145cdae4e
# ╟─a0573fe0-d6a2-11ea-057e-5140c323c005
# ╠═90208f82-d4ba-11ea-11e6-7f51e6ffc5d5
# ╟─1790c760-0e0d-11eb-1c47-f5e902e158fb
# ╠═f5ca89e0-0e0c-11eb-32a1-fb22511f42ae
# ╠═165464e2-0e0f-11eb-3851-3701a6c78502
# ╠═915aa560-0e0e-11eb-18d4-e3a899ee2c1c
# ╠═7115307e-d660-11ea-3f2f-4d4963a60e39
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
# ╟─08c0cb70-d7fa-11ea-0efd-ab3b7f76a8b8
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
# ╟─351f5332-dd7c-11ea-2fc4-c9dc71684f23
# ╟─5fd238f0-dd7b-11ea-2a00-a587ee59c261
# ╠═f102ad80-dd78-11ea-08b5-2d95e6f1a678
# ╟─db58def0-dd78-11ea-335b-4140d22d9a78
# ╟─46dd7a80-de3e-11ea-0e9d-abe7496efef7
# ╟─44e6bce2-800b-11eb-3651-e3ba5df8617e
# ╟─33af7e82-d7f5-11ea-12bb-5982575f35d3
# ╠═3d70984e-d7f5-11ea-3243-9d3932ee155b
# ╠═afa04610-d7f4-11ea-1cd0-e3a91d91c5c9
# ╟─713de9b0-3490-11eb-041b-c7083c829ba7
# ╠═98165ef0-3490-11eb-0e57-259fcf7361c2
# ╟─10625a80-7122-11eb-3601-273ea4390a23
# ╠═ea63b4f0-3490-11eb-26a8-6f5e1f2d3c05
# ╠═a07f3bc0-3490-11eb-11b3-6bedb92502cb
# ╟─b8326ede-3490-11eb-1e0b-530752ed836e
# ╟─d43266e0-3490-11eb-0d9a-59ef7a663b50
# ╠═e042eeee-3490-11eb-1406-058ee6753c1f
# ╠═165343ee-3491-11eb-10a7-7b4d6fef56f4
# ╟─315a35f0-3491-11eb-0ce1-636a82ef3dd9
# ╠═dd340bd0-3491-11eb-31ab-5b1bf6286cf0
# ╠═3765e8e2-3491-11eb-00a5-274f3caf6a15
# ╠═de04c952-3491-11eb-276b-b530d61da953
# ╠═2e058370-729c-11eb-0822-793d3db17e74
# ╠═39b45040-3492-11eb-1e48-998b1e4b1911
# ╟─4c657840-3492-11eb-093d-59a4b3854b7a
# ╟─55e2ab40-3492-11eb-2953-e9d8bd17735c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
