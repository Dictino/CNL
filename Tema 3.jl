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

# ╔═╡ c09c0e20-d8be-11ea-11be-0966527e6f2d
begin
    import Pkg
	  Pkg.activate(@__DIR__)
	using Plots, PlutoUI
end

# ╔═╡ b711bf40-c912-11eb-120d-1bb53ca25d17
using ControlSystems

# ╔═╡ cfbeb720-90a8-11eb-1ac2-3bda60b2285c
using OrdinaryDiffEq

# ╔═╡ 889652b2-dbe8-11ea-0835-e966a39af471
md"""# Problema de estabilidad absoluta de Lur`e
## Partimos de una planta con función de transferencia
"""

# ╔═╡ bf0a5c82-d8bc-11ea-3b1c-b397e3b526fa
H1(s)=(s+3)/(s^2+7s+10)

# ╔═╡ 6b0b4680-9095-11eb-3afe-9b792355080b
md"lo primero que miramos es qué polos tiene, factorizando el denominador:

$H_1=\frac{s+3}{(s+5)(s+2)}$

No hay polos ni ceros en el semiplano derecho ni tampoco en el origen.
"

# ╔═╡ 9999f2c0-dbdd-11ea-04a5-dd56cf1eae5c
md"""# Conjeturas A-K
## Se trata de comprobar el rango para el cuál la realimentación negativa con ganancia K es estable:
+ Se calcula el número de polos $n$ en el semiplano derecho (inestables)
+ El criterio de Nyquist dice que la planta es estable si la respuesta en frecuecia rodea al punto $-1/K+0*i$ $n$ veces en sentido antihorario.
"""

# ╔═╡ 4628f590-dbe3-11ea-2a7a-83cabda49ee7
function conjeturas(H,K)	
	#respuesta en frecuencia
	wp=10.0.^(-50:0.01:50) #rango de 10^-50 a 10ymin=^50
	wn=-10.0.^(50:-0.01:-50) 
	w=vcat(wn,0,wp)#todas las frecuencias incluidas las negativas
	Hs=H.(w*im)
	
	#si hay polos en el eje imaginario habrá valores infinitos o NaN para pintar bien esos casos los "recorto", en este caso no se cumplen las condiciones de las conjeturas pero puede igualemnte ser interesante ver el rango de estabilidad. 
	
	recortar(x)=min(max(x,-10^3),10^3)
	
	#pinta la flecha en la frecuencia w0
if isfinite(maximum(norm.(Hs)))
	#pinto normalmente
	plot(real.(Hs),imag.(Hs), label=false)
	#pinta la flecha
	quiver!([real(H(0))],[imag(H(0))],
			quiver=([real(H(0.01*im)-H(0))],
				[imag(H(0.01*im)-H(0))]),label=false)
else
	#en este caso es cuando hay que recortar
	plot(recortar.(real.(Hs)),recortar.(imag.(Hs)), label=false)
	ylims!((-10,10))
	xlims!((-10,10))
end

#para que el movimiento sea intuitivo el slider ha de ir desde la derecha 

#pinta el punto
	scatter!([-1/K],[0],label=false)
	title!("""La respuesta en frecuencia ha de rodear al punto 
	en sentido antihorario tantas veces como polos inestables
	hay, K=$(K)""")
	#para usar fuera de pluto es necesario hacer display()
end

# ╔═╡ 490c8852-909f-11eb-2522-ad23744ce9d8
	#CSS que hace los Sliders más grandes
    html"""<style type="text/css">
    input[type=range] {
    width: 75%;
    } </style>
    """

# ╔═╡ 203ea980-9090-11eb-25f8-8ff76ab5262a
begin
	K_max_conjeturas=10
	@bind sliderAK Slider((-1:0.001:1)*K_max_conjeturas,default=1.0)
end

# ╔═╡ 93bfc0d0-cb03-11eb-283b-6d01e65bed3c
conjeturas(H1,sliderAK)

# ╔═╡ 29e11c20-90ae-11eb-1348-0f5a981d3361
md"Esto sugiere que es estable en el sector $[-\frac{10}{3}, \infty)$
Podemos ver analíticamente el -10/3 observando que el corte es el $w_c=0$ y por lo tanto $H(w_c)=\frac{0+3}{(0+5)(0+2)}=\frac{3}{10}
 \to H(0)k=-1 \to k=-10/3$

Veamos qué dicen los criterios *de verdad* y no las conjeturas
"

# ╔═╡ 66f0b770-d8c5-11ea-0ee3-7ba9b6dfb585
md"""# Criterio del círculo
## Antes de empezar comprueba que:
+ Como es una funcion de transferencia ya es contolable/observable (OK)
+ No hay integradores ni realimentación directa $$d=0$$
¿Es estable? Usa los casos A y B. Si no lo es usa los casos C y D.

**Cada caso** de los que se puede usar da **un rango** posible de ganancias
"""

# ╔═╡ 9ca66610-d985-11ea-0426-13ef5f785091
function circulo(H,K1,K2,caso)	
	"""hay cuatro casos se pueden determinar por los signos de k1 y k2
	------ Estables
	A) 0=k1<k2 recta a la izquierda
	B) k1<0<k2 dentro del circulo
	------ Inestables
	C) 0<k1<=k2 La respueta rodea al circulo en sentido antihorario n veces donde n es el número de polos inestables
	D) k1<k2<0 La respuesta de -H rodea en el sentido antihorario n veces el circulo D(-k1,-k2) es como C pero con todo negado"""

	#compruebo que las entradas son correctas
	if !isfinite(H(0))
		return "No puede haber polos en el origen"
	end
	if caso=='A'
		if (K1==0)&(0<K2)
			titulo="Recta a la izquierda"
		else
			return "Error: en el caso A 0=k1<k2"
		end
	elseif caso=='B'
		if (K1<0)&(0<K2)
			titulo="El circulo rodea la respuesta en frecuencia"
		else
			return "Error: en el caso B k1<0<k2"
		end
	elseif caso=='C'
		if (0<K1)&(K1<=K2)
			titulo="La respueta rodea al circulo en sentido antihorario n veces\n donde n es el número de polos inestables"
		else
			return "Error: en el caso C 0<k1<=k2 "
		end
	elseif caso=='D'
		if (K1<=K2)&(K2<0)
			H2(x)=-H(x)
            return circulo( H2,-K2,-K1,'C')
		else
			return "Error: en el caso D k1<k2<0"
		end
	else
		return "No es un caso valido, ha de ser 'A','B','C' o 'D'"
	end	
	#por fin los cálculos...
	#respuesta en frecuencia
	    wp=10.0.^(-50:0.01:50) #rango de 10^-50 a 10ymin=^50
	    wn=-10.0.^(50:-0.01:-50) 
		w=vcat(wn,wp)#todas las frecuencias incluidas las negativas
	   
	    Hs=H.(w*im)
	    plot(Hs,line=:arrow,label=false)
	    Hs=collect(Hs)
	    plot(real.(Hs),imag.(Hs), line=:arrow,label=false)
	    #pinta la flecha
	    quiver!([real(H(0))],[imag(H(0))],
		        quiver=([real(H(0.01*im)-H(0))],[imag(H(0.01*im)-H(0))])
		        ,label=false)
	
	#pinta el "circulo"
	if caso=='A'
		#el círculo es una recta 1/K1=inf
		ymin=minimum(imag.(Hs))
	    ymax=maximum(imag.(Hs))
	    y=[ymin,ymax]
		x=[-1/K2,-1/K2]
	    plot!(x,y,label=false)
	    title!(titulo)
	else
		#el círuclo es un círculo de vedad
		p1=-1/K1
		p2=-1/K2 #puntos de corte
	    centro=(p1+p2)/2
	    radio=abs(p2-p1)/2
	    theta=collect(0:0.01:2*pi)
	    circunferencia=centro.+radio*exp.(im*theta)
		scatter!([-1/K1 -1/K2],[0 0],label=false)
	    plot!(circunferencia,label="circulo")
		title!(titulo)
	end
	#idem aquí, para usar fuera de pluto hay que hacer display()
end

# ╔═╡ 00b900b0-d8c6-11ea-0e42-0f4bfaf17fd2
md"## Planta estable"

# ╔═╡ ed754f20-db56-11ea-02fa-3581c244dcbd
md"""### Caso A
$$0=k_1<k_2$$ 
La recta ha de quedar a la izquierda
"""

# ╔═╡ 9975e6e0-db57-11ea-30e1-3ffbd104839c
K_maxA=100

# ╔═╡ adee0ba0-d988-11ea-19c7-e324f1f86031
@bind k2a Slider((0:0.001:1)*K_maxA,default=1.0)

# ╔═╡ 839b4c72-9092-11eb-2237-1dc284e3041e
k2a

# ╔═╡ b642eb40-d988-11ea-26a4-374a74d5dc08
circulo(H1,0,k2a,'A')

# ╔═╡ 2b1196f0-90dc-11eb-2e20-89a61fcb1498
md"Este caso ha sido más conservador que las conjeturas $[0, \infty)$"

# ╔═╡ b47374ee-db55-11ea-2b8f-5b38e3a2c5a0
md"""### Caso B
$$k_1<0<k_2$$ 
El círculo ha de rodear la respuesta en frecuencia
"""

# ╔═╡ bd8ef6f0-dbb8-11ea-3fe0-3b17ccdb1fcd
begin
	K_maxB1=10
	K_maxB2=100
end

# ╔═╡ fe58cf80-dbb8-11ea-1c6c-93fef94ed79c
@bind k1b Slider((-1:0.001:0)*K_maxB1,default=-1.0)

# ╔═╡ d2de7080-dbb8-11ea-1d7e-11b40cd2ce8d
@bind k2b Slider((0:0.001:1)*K_maxB2,default=1.0)

# ╔═╡ d3077b80-9092-11eb-0f8b-85c6edde81f5
(k1=k1b,k2=k2b)

# ╔═╡ 9d833cf0-db03-11ea-3b28-33343306f4dc
circulo(H1,k1b,k2b,'B')

# ╔═╡ 65071560-90dc-11eb-24f3-3d03b8461507
md"Este caso es igual que las conjeturas AK es decir $[-\frac{10}{3}, \infty)$ "

# ╔═╡ 1ae51110-d98c-11ea-3b02-1ddccacadb56
md"## Planta inestable
En este caso **la planta de antes no nos sirve** ya que era estable, usamos estra otra

$H_2=\frac{s+3}{(s+5)(s-2)}$

Que tiene **un** polo inestable

"

# ╔═╡ 8054e0c0-90a2-11eb-20fb-ad26bded4db5
H2(s)=(s+3)/(s^2+3s-10)

# ╔═╡ 7697fac0-90a4-11eb-1687-5ff2a5cf415d
md"Primero veamos que dicen las conjeturas"

# ╔═╡ ab676b50-90a4-11eb-1632-f1f60fd49f07
begin
	K_max_conjeturas2=10
	@bind sliderAK2 Slider((-1:0.001:1)*K_max_conjeturas2,default=1.0)
end

# ╔═╡ 81144e40-90a4-11eb-2941-f53c46a9300c
conjeturas(H2,sliderAK2)

# ╔═╡ c2d18d10-90dc-11eb-1792-3779513d19cb
md"Las conjeturas dicen $[\frac{10}{3}, \infty)$"

# ╔═╡ f9938e60-d98b-11ea-1e6d-2fc1399f7bef
md"""### Caso C
$$0<k_1 \le k_2$$ 

Como es inestable habrá $n$ polos en el semiplano derecho.
El círculo ha de **estar rodeado** por la respuesta en frecuencia $n$ veces en sentido antihorario.
"""

# ╔═╡ 1be6ec20-dbba-11ea-2479-ed82409ebdd7
begin
	K_maxC1=10
	K_maxC2=100
end

# ╔═╡ 76a9f030-dbba-11ea-3326-95356d82a2c5
@bind k1c Slider((0:0.001:1)*K_maxC1,default=1.0)

# ╔═╡ 783969d0-dbba-11ea-10d0-b1cea6c3b727
@bind k2c Slider((0:0.001:1)*K_maxC2,default=2.0)

# ╔═╡ 5d374e70-9093-11eb-3126-bfa8080d6bc2
(k1=k1c,k2=k2c)

# ╔═╡ e3f469a0-dbb9-11ea-30bb-813e5688c466
circulo(H2,k1c,k2c,'C')

# ╔═╡ e748e0ce-90dc-11eb-14cf-5f0c6c4ae075
md"Lo mismo que las conjeturas   $[\frac{10}{3}, \infty)$"

# ╔═╡ c45e9020-d989-11ea-028f-5f36703d75da
md"""### Caso D
$$k_1<k_2<0$$ 

Es como el caso C pero con todos los signos cambiados y H cambiado por -H:
El círculo ha de **estar rodeado** por la respuesta en frecuencia (de -H) $n$ veces en sentido antihorario.
"""

# ╔═╡ 245ff590-dbba-11ea-1deb-1976fc8d165d
begin
	K_maxD1=100
	K_maxD2=100
end

# ╔═╡ 8d0b4eee-dbba-11ea-2fad-5d7ed6f8035b
@bind k1d Slider((-1:0.001:0)*K_maxD1,default=-2.0 )

# ╔═╡ 8f55db30-dbba-11ea-0370-ddadc3b11dfe
@bind k2d Slider((-1:0.001:0)*K_maxD1,default=-1.0)

# ╔═╡ 3f2dfc70-9094-11eb-340c-592edc764cd1
(k1d=k1d,k2d=k2d)

# ╔═╡ 11a30cf0-9094-11eb-2af0-6f424889df66
circulo(H2,k1d,k2d,'D')

# ╔═╡ 105992d0-90dd-11eb-1a3d-653111528c53
md"No es posible meter el círculo dentro de la respuesta en frecuencia, este caso **no nos dice nada**"

# ╔═╡ 8e2736c0-d8bb-11ea-2fde-43a428dd9281
md"""# Criterio de popov
## Antes de empezar comprueba que:
+ Como es una funcion de transferencia ya es contolable/observable
+ Hay a lo sumo un integrador y si lo hay la planta se puede poner como $$d/s+ H_{estable}$$ con d>0 es decir
    + La parte del ingregrador es positiva ($$H(s) \to \infty$$  si  $$s\to 0^+$$) 
    + Quitando el integrador la planta es estable"""

# ╔═╡ 7c27fd1e-c90f-11eb-0edf-29b2a8b2867b
md"Como ejemplo podemos usar la planta H1 añadiendo un polo en el origen:

$H_3=\frac{s+3}{s(s+5)(s+2)}$"

# ╔═╡ b7044bc0-9880-11eb-31b6-b1578086897a
H3(s)=H1(s)/s #OJO, cuando useis esto poner VUESTRA PLANTA no la misma del círculo dividida por s...

# ╔═╡ bc4b38a0-c913-11eb-2eb9-25b5e6468a61
md"""### Primero a ver que nos dicen las cojneturas AK
Siendo rigurosos no nos dicen **nada** ¿Por qué?

El problema es que no **son aplicables** si hay un polo con parte real igual a cero como es el caso...

Pero al fin y al cabo son conjeturas y lo que queremos es una idea *optimista* de dónde buscar y las conjeturas nos la dan ya que si una función cualquiera del sector es estable ha de serlo también el caso particular $f(x)=kx$ que es lo que calcularemos aquí.
"""

# ╔═╡ 31d44760-c90f-11eb-3a7a-ffb4d6c20b53
begin
	K_max_conjeturas3=10
	@bind sliderAK3 Slider((-1:0.001:1)*K_max_conjeturas3,default=1.0)
end

# ╔═╡ 439f2e10-c90f-11eb-1b76-eb23238ba418
begin
	conjeturas(H3,sliderAK3)
	ylims!(-1,1) #hay que hacer zoom para ver algo
	xlims!(-Inf,Inf) #esto quita el límite de x y lo pone automático
end

# ╔═╡ 710513f0-c910-11eb-2abb-3df94293cdb4
md"**OJO** Ahora no está claro cuando la respuesta *rodea* al punto ya que es abierta..

Lo que si es importante es ver cuando **cambiamos de lado** (la estabilidad cambia al pasar de un lado de la curva al otro lado).

En este caso se cambia de lado cuando $k=0$ y nos aproximamos a la curva cuando $k=\pm \infty$

Luego basta con tomar un ejemplo concreto de ganancia ($k=1$ por ejemplo) y ver si es estable para ver en qué lado es estable y en cual inestable)

Realimentando con $k=1$ tenemos

$H_{LazoCerrado}=\frac{H_3}{1+H_3}=\frac{\frac{s+3}{s(s+5)(s+2)}}{1+\frac{s+3}{s(s+5)(s+2)}}=\frac{s+3}{s(s+5)(s+2)+(s+3)}$

$H_{LazoCerrado}=\frac{s+3}{s^3+7s^2+11s+3}$

Que es estable por Routh–Hurwitz ya que 

$7>0$
$3>0$
$7 \cdot 11>3$

Luego el lado que contine $k=1$ es estable y el otro inestable por tanto el intervalo estable es $[0, \infty)$.

Si no queremos hacerlo a meno podemos usar la toolbox de control:"



# ╔═╡ d8394bbe-c912-11eb-250f-1b7385c0a536
S=tf([1.0, 0.0],[1.0])

# ╔═╡ 76c049a0-90ac-11eb-1892-a31cb390383d
md"Evaluando H(s) con S nos da la función de transferencia (esto funciona sólo por que H es racional)"

# ╔═╡ 33a205b0-c913-11eb-22a9-d51437c6252c
H3s=H3(S)

# ╔═╡ 420138b0-c913-11eb-2509-558e30b8048b
#el lazo cerrado con ganancia 1
H_LC=minreal(H3s/(1+H3s)) #minreal es para cancelar los ceros y polos sobrantes

# ╔═╡ 5c92fab0-c913-11eb-0542-2bd8ac7560de
#y por fin los ceros y los polos
zpk(H_LC)

# ╔═╡ 2cfce6fe-c915-11eb-1f34-5be3aed8045a
md"## Ahora a ver qué dice Popov"

# ╔═╡ 497f9250-d8bc-11ea-194d-77979b94814d
function popov(H,K,alpha)
	if K<0
		return "Error, en Popov K>=0"
	end
	W(w)=real(H(w*im))+w*im*imag(H(w*im)) #funcion asociada
	wp=10.0.^(-50:0.01:50) #rango de 10^-50 a 10ymin=^50
	wn=-10.0.^(-50:0.01:50) 
	w=hcat(wn,wp) #todas las frecuencias incluidas las negativas
	Ws=W.(w)
	p=plot(real(Ws),imag(Ws), line = :arrow,label=false)
	#recta
	xmin=minimum(real.(Ws))
	xmin=min(xmin,-1/K) #me aseguro de pintar el corte en el eje
	xmax=maximum(real.(Ws))
	x=[xmin,xmax]
	y=(x.+1/K)./alpha
	scatter!([-1/K],[0],label=false)
	plot!(x.+y.*im,label=false)
	title!("La recta ha de estar por encima de la función asociada")
end

# ╔═╡ fee20250-d8bb-11ea-299c-6d093b0edd20
K_max=1000

# ╔═╡ cb075200-d8bb-11ea-234b-6d4b6fca9083
@bind K Slider(0.0:0.01:K_max,default=10.0)

# ╔═╡ f2567110-dbe0-11ea-24a0-4d1a72151516
@bind pendiente Slider(0.0:0.1:90.0,default=45.0)

# ╔═╡ f8ad6aa0-d8bb-11ea-0dea-85d5051066fe
md"K=$(K)"

# ╔═╡ e443a3e0-d8c0-11ea-362f-75314f8b2cbb
alpha=1/tan(pendiente*pi/180)

# ╔═╡ dbcceee0-d985-11ea-2a3e-f9c3ec08cb75
popov(H3,K,alpha)

# ╔═╡ 367c9ac0-90dd-11eb-011b-0f7c4fca2c89
md"Es decir  $[0, \infty)$"

# ╔═╡ 3b303230-90a5-11eb-148b-d9c93a8c3772
md"# Simulando el sistema
Hasta ahora hemos usado funciones de tranferencia, pero ¿cómo lo simulamos?
Lo más fácil es pasar al espacio de estados

$\dot x =Ax+Bu$  

$y=Cx$

Asumiremos siempre en este tema que $D=0$
"

# ╔═╡ e7264d30-90a6-11eb-07fc-d1f2120385cd
H=H1(S)

# ╔═╡ 98a5d8b2-90dd-11eb-3061-1713c5934c8c
md"Obtenido esto podemos pasar al espacio de estados"

# ╔═╡ f2cdb330-90a6-11eb-3169-8fb43ab0c79c
estados=ss(H)

# ╔═╡ a4378a10-90ac-11eb-0041-a9145a07928c
md"Definimos ahora el sistema"

# ╔═╡ 17f33990-90ad-11eb-2fa5-b762b4f1138b
md"# Sector
Una función f(x) pertenece al sector $[k_1 k_2]$ si está entre las rectas $y=k_1x$ y $y=k_2x$, algebráicamente:

$k_1x^2 \le xf(x) \le k_2x^2$ 

u obviamdo la división por cero (suponiendo que $\frac{f(x)}{x}$ sea continua en $0$ definiendola como el $lim_{x\to 0}\frac{f(x)}{x}$)

$k_1\le \frac{f(x)}{x} \le k_2 \  \forall x\ne0$ 
"

# ╔═╡ b1a75e90-90ad-11eb-0f7d-1d2b98b9fce4
f(x)=x*abs(x)-10/3*sin(x)

# ╔═╡ a9dd2000-90df-11eb-0edb-e364fefafe32
@bind k1 Slider(-10:0.01:10,default=-1.0)

# ╔═╡ 3387bd60-90e0-11eb-1963-7dfc892de25c
@bind k2 Slider(0:0.01:100,default=1.0)

# ╔═╡ a134fad0-90e0-11eb-1145-df3b4d6da246
@bind maximo Slider(1:0.01:100,default=1.0)

# ╔═╡ 041d31e0-90b3-11eb-1d07-c3b6bb15e4ae
begin
	x=-maximo:0.01:maximo
	fx=f.(x)
	f_max=maximum(fx)
	f_min=minimum(fx)
	fig1=plot(x,fx)
	plot!(x,k1*x,label="$(k1)*x",color=:red)
	plot!(x,k2*x,label="$(k2)*x",color=:red)
	sector1 = Shape([(0.0, 0.0), (maximo, k1*maximo), 
			   (maximo, k2*maximo), (0.0, 0.0)])
	sector2 = Shape([(0.0, 0.0), (-maximo, -k1*maximo), 
			   (-maximo, -k2*maximo), (0.0, 0.0)])
    plot!(sector1, fillcolor = plot_color(:yellow, 0.3),label=false)
	plot!(sector2, fillcolor = plot_color(:yellow, 0.3),label=false)
	plot!(ylims = (f_min,f_max))
	title!("La función ha de estar dentro de la zona amarilla")
	
	fig2=plot(x,fx./x,label="f(x)/x")
	sector3 = Shape([(-maximo, k1), (maximo, k1), 
			   (maximo, k2), (-maximo, k2),(-maximo, k1)])
	plot!(sector3, fillcolor = plot_color(:yellow, 0.3),label=false)
	
	l = @layout [a; b]
	plot(fig1,fig2, layout=l)
end

# ╔═╡ d6fd4a4e-90e0-11eb-34bf-ab3bb29cc905
md"Jugando un poco es fácil ver que es  $[-\frac{10}{3}, \infty)$
Para verlo analíticamente habría que hacer

$\frac{f(x)}{x}=|x|-\frac{10}{3}\frac{sin(x)}{x}$

Y ver que toma cualquier valor entre $-\frac{10}{3}$ (cuando $x \to 0$) e $\infty$ (cuando $x \to \infty$) 
"

# ╔═╡ af179aae-90a7-11eb-2a7b-1b96a28c8ebc
function control(y)
	r=0
    u=r-f(y)
end

# ╔═╡ 21c1cfa0-90a7-11eb-3c54-23fa9b677080
function derivadas(x,p,t)
	A,B,C=p #extraemos los parámetros
	y=C*x
	y=y[1] #no quiero un vector sino un escalar [2.5] no es lo mismo que 2.5
	u=control(y) #calculamos el control
	dx=A*x+B*u
end

# ╔═╡ 1fc9b3b0-90a8-11eb-3b34-77d759a763d5
let
	x0 = zeros(size(estados.B))
	x0[1]=1
	parametros=(estados.A, estados.B, estados.C)
	tspan = (0.0, 30.0)
	prob = ODEProblem(derivadas,x0,tspan,parametros);
	sol = solve(prob,Tsit5());
	
	# Salida y referencia
	#la salida y el control no se guardan así que los reconstruyo
	ts=tspan[1]:0.01:tspan[end]
	X=sol.(ts)
	y=[(estados.C*x)[1] for x in X] #explicar paso a paso
	u=control.(y)
	
	plot(ts,y,xaxis="t",yaxis="y(t)",label="y(t)")
	#referencia
	fig1=plot!(ts,0.0*ts, linestyle=:dot, linewidth=2, label="referencia") 
		
	
	#Señal de control, u no se guarda, hay que recalcularlo
	fig2=plot(ts,u,xaxis="t",yaxis="u(t)",legend=false)
	l = @layout [a ; b]
	plot(fig1,fig2, layout=l)
end

# ╔═╡ Cell order:
# ╟─889652b2-dbe8-11ea-0835-e966a39af471
# ╠═bf0a5c82-d8bc-11ea-3b1c-b397e3b526fa
# ╟─6b0b4680-9095-11eb-3afe-9b792355080b
# ╟─9999f2c0-dbdd-11ea-04a5-dd56cf1eae5c
# ╠═4628f590-dbe3-11ea-2a7a-83cabda49ee7
# ╠═490c8852-909f-11eb-2522-ad23744ce9d8
# ╠═203ea980-9090-11eb-25f8-8ff76ab5262a
# ╠═c09c0e20-d8be-11ea-11be-0966527e6f2d
# ╠═93bfc0d0-cb03-11eb-283b-6d01e65bed3c
# ╟─29e11c20-90ae-11eb-1348-0f5a981d3361
# ╟─66f0b770-d8c5-11ea-0ee3-7ba9b6dfb585
# ╠═9ca66610-d985-11ea-0426-13ef5f785091
# ╟─00b900b0-d8c6-11ea-0e42-0f4bfaf17fd2
# ╟─ed754f20-db56-11ea-02fa-3581c244dcbd
# ╠═9975e6e0-db57-11ea-30e1-3ffbd104839c
# ╠═adee0ba0-d988-11ea-19c7-e324f1f86031
# ╠═839b4c72-9092-11eb-2237-1dc284e3041e
# ╠═b642eb40-d988-11ea-26a4-374a74d5dc08
# ╟─2b1196f0-90dc-11eb-2e20-89a61fcb1498
# ╟─b47374ee-db55-11ea-2b8f-5b38e3a2c5a0
# ╠═bd8ef6f0-dbb8-11ea-3fe0-3b17ccdb1fcd
# ╠═fe58cf80-dbb8-11ea-1c6c-93fef94ed79c
# ╠═d2de7080-dbb8-11ea-1d7e-11b40cd2ce8d
# ╟─d3077b80-9092-11eb-0f8b-85c6edde81f5
# ╠═9d833cf0-db03-11ea-3b28-33343306f4dc
# ╟─65071560-90dc-11eb-24f3-3d03b8461507
# ╟─1ae51110-d98c-11ea-3b02-1ddccacadb56
# ╠═8054e0c0-90a2-11eb-20fb-ad26bded4db5
# ╟─7697fac0-90a4-11eb-1687-5ff2a5cf415d
# ╟─ab676b50-90a4-11eb-1632-f1f60fd49f07
# ╠═81144e40-90a4-11eb-2941-f53c46a9300c
# ╟─c2d18d10-90dc-11eb-1792-3779513d19cb
# ╟─f9938e60-d98b-11ea-1e6d-2fc1399f7bef
# ╠═1be6ec20-dbba-11ea-2479-ed82409ebdd7
# ╠═76a9f030-dbba-11ea-3326-95356d82a2c5
# ╠═783969d0-dbba-11ea-10d0-b1cea6c3b727
# ╟─5d374e70-9093-11eb-3126-bfa8080d6bc2
# ╠═e3f469a0-dbb9-11ea-30bb-813e5688c466
# ╟─e748e0ce-90dc-11eb-14cf-5f0c6c4ae075
# ╟─c45e9020-d989-11ea-028f-5f36703d75da
# ╠═245ff590-dbba-11ea-1deb-1976fc8d165d
# ╠═8d0b4eee-dbba-11ea-2fad-5d7ed6f8035b
# ╠═8f55db30-dbba-11ea-0370-ddadc3b11dfe
# ╟─3f2dfc70-9094-11eb-340c-592edc764cd1
# ╠═11a30cf0-9094-11eb-2af0-6f424889df66
# ╟─105992d0-90dd-11eb-1a3d-653111528c53
# ╟─8e2736c0-d8bb-11ea-2fde-43a428dd9281
# ╟─7c27fd1e-c90f-11eb-0edf-29b2a8b2867b
# ╠═b7044bc0-9880-11eb-31b6-b1578086897a
# ╟─bc4b38a0-c913-11eb-2eb9-25b5e6468a61
# ╟─31d44760-c90f-11eb-3a7a-ffb4d6c20b53
# ╠═439f2e10-c90f-11eb-1b76-eb23238ba418
# ╟─710513f0-c910-11eb-2abb-3df94293cdb4
# ╠═b711bf40-c912-11eb-120d-1bb53ca25d17
# ╠═d8394bbe-c912-11eb-250f-1b7385c0a536
# ╟─76c049a0-90ac-11eb-1892-a31cb390383d
# ╠═33a205b0-c913-11eb-22a9-d51437c6252c
# ╠═420138b0-c913-11eb-2509-558e30b8048b
# ╠═5c92fab0-c913-11eb-0542-2bd8ac7560de
# ╟─2cfce6fe-c915-11eb-1f34-5be3aed8045a
# ╠═497f9250-d8bc-11ea-194d-77979b94814d
# ╠═fee20250-d8bb-11ea-299c-6d093b0edd20
# ╠═cb075200-d8bb-11ea-234b-6d4b6fca9083
# ╠═f2567110-dbe0-11ea-24a0-4d1a72151516
# ╟─f8ad6aa0-d8bb-11ea-0dea-85d5051066fe
# ╟─e443a3e0-d8c0-11ea-362f-75314f8b2cbb
# ╠═dbcceee0-d985-11ea-2a3e-f9c3ec08cb75
# ╟─367c9ac0-90dd-11eb-011b-0f7c4fca2c89
# ╟─3b303230-90a5-11eb-148b-d9c93a8c3772
# ╠═e7264d30-90a6-11eb-07fc-d1f2120385cd
# ╟─98a5d8b2-90dd-11eb-3061-1713c5934c8c
# ╠═f2cdb330-90a6-11eb-3169-8fb43ab0c79c
# ╠═a4378a10-90ac-11eb-0041-a9145a07928c
# ╠═21c1cfa0-90a7-11eb-3c54-23fa9b677080
# ╟─17f33990-90ad-11eb-2fa5-b762b4f1138b
# ╠═b1a75e90-90ad-11eb-0f7d-1d2b98b9fce4
# ╠═a9dd2000-90df-11eb-0edb-e364fefafe32
# ╠═3387bd60-90e0-11eb-1963-7dfc892de25c
# ╠═a134fad0-90e0-11eb-1145-df3b4d6da246
# ╟─041d31e0-90b3-11eb-1d07-c3b6bb15e4ae
# ╟─d6fd4a4e-90e0-11eb-34bf-ab3bb29cc905
# ╠═af179aae-90a7-11eb-2a7b-1b96a28c8ebc
# ╠═cfbeb720-90a8-11eb-1ac2-3bda60b2285c
# ╠═1fc9b3b0-90a8-11eb-3b34-77d759a763d5
