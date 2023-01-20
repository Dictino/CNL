### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ ce44280b-1637-4866-b473-b9941c76c06d
begin
    import Pkg
	  Pkg.activate(@__DIR__)
	using Plots, PlutoUI
end

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
	return p1 #si lo quisieseis usar fuera de pluto hay que usar display(figura)
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

# ╔═╡ Cell order:
# ╟─6ba02e80-d4cf-11ea-19cb-33028d3dc67d
# ╟─7db3fdb0-d4d2-11ea-21fc-597c55dda572
# ╟─76079a70-d4cf-11ea-2645-8f43c78eae52
# ╠═66e9ffd0-d4d7-11ea-09c0-03ebaacab879
# ╟─0b94f1c0-7113-11eb-13da-c5c5b95ebbf5
# ╠═e885f5b0-d4d4-11ea-06d6-9ddf9619961f
# ╠═b99fc9b0-d4d9-11ea-0e18-cfd5af90188f
# ╠═ce44280b-1637-4866-b473-b9941c76c06d
# ╠═3be6b020-d4d8-11ea-3356-23afd148ebc9
# ╠═220f6f12-d4d9-11ea-159e-c1e1551fd732
# ╟─23f5ece0-d4da-11ea-1330-c5d077bd3a89
# ╠═d8b5ed40-d4d7-11ea-27d4-bf594518802b
# ╠═2b7df180-d4d8-11ea-118f-75e8709e36c3
# ╟─455e0c30-d4d7-11ea-22ce-0d04b4e8d860
# ╠═159f28e2-d4db-11ea-13e6-a97e2e7f7c60
# ╠═c90980d0-d637-11ea-112d-c336b2080d87
# ╠═1f44a88e-7114-11eb-0ac2-d12cd3fb40ef
# ╠═4a7f5f72-0e0c-11eb-2b7c-4d7214c09c4b
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
