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

# ╔═╡ 02c15be0-defd-11ea-004b-45c4e09aa541
begin
    import Pkg
	  Pkg.activate(@__DIR__)
	using Plots, PlutoUI, SymEngine
end

# ╔═╡ 6cbb32f0-d665-11ea-3f24-41cbe1910b80
using LinearAlgebra

# ╔═╡ 28fccf20-d63c-11ea-3311-411f89a63431
md"# Curvas de nivel"

# ╔═╡ 309959f0-defd-11ea-34c3-71eeda7b523b
md"Curvas cerradas:

$$V(x)\to \infty$$ cuando $$x\to \infty$$"

# ╔═╡ cd88fade-defd-11ea-2097-25b2d1eac285
V1(x1,x2)=x1^2+x2^2

# ╔═╡ ce5fe800-deff-11ea-0b19-7fc23b14bd06
x=y=-1:0.1:1

# ╔═╡ c6c49460-deff-11ea-1c6d-673d863ea07b
begin
	plotly() #Back-end interactivo de Plots
	surface(x,y,V1)
end

# ╔═╡ 5b26b050-d7eb-11ea-2488-b7da296fed24
@bind L1 Slider(0:0.05:2)

# ╔═╡ 2e580002-defe-11ea-1b3d-d5ec461e0715
begin
	gr() #Back-End por defecto
	contour(x, y, V1, fill = true)
	contour!(x,y,V1,levels=[L1],linewidth=3,color=:red, ratio=1)
end

# ╔═╡ 2ec228f0-defd-11ea-329f-53576ab14004
md"Curvas abiertas

$$V(x) \nrightarrow \infty$$ cuando $$x\to \infty$$"

# ╔═╡ ca2fea72-3655-11eb-0984-a532e1f9efcc
md"""Por ejemplo 

$V_2(x_1,x_2)=\frac{x_1^2}{1+x_1^2}+x_2^2$

Ya que si $x_1\to \pm\infty$, $V_2 \to 1+x_2^2 \ne   \infty$"""

# ╔═╡ 204e6710-defe-11ea-1fb0-43c7ef863786
V2(x1,x2)=x1^2/(1+x1^2)+x2^2;

# ╔═╡ f3bcc450-df00-11ea-0080-03b54212b9ed
begin
	x2=-5:0.1:5
	y2=-1:0.1:1
	plotly()
	surface(x2,y2,V2,ratio=1)
end

# ╔═╡ ca99149e-df02-11ea-39a0-179e6a7d1344
@bind L2 Slider(0:0.05:2)

# ╔═╡ eb00be70-df00-11ea-0ce9-3b0a79a0f4a3
begin
	gr()
	contour(x2, y2, V2, fill = true)
	contour!(x2,y2,V2,levels=[L2],linewidth=3,color=:red) #, ratio=1)
end

# ╔═╡ d8678810-defc-11ea-05b2-51373f754e61
md"# Semi-Definida Vs Definida
Pensemos **antes** de pintarlo, para ver si hay *sorpresas*"

# ╔═╡ 34a9ca60-df03-11ea-018c-399409d7842e
begin
	f(x1,x2)=x1^2 + x2^2
	#f(x1,x2)=(x1+x2)^2
	#f(x1,x2)=-x1^2
	#f(x1,x2)=-2x1^2 - x2^2 - 2x1*x2
	#f(x1,x2)=-x1^2 + x2^2
end

# ╔═╡ f94dece0-defc-11ea-0546-39a094bbd2b6
begin
	gr() 
	contour(x, y, f, fill = true)
	contour!(x,y,f,levels=[0.0],linewidth=3,color=:red, ratio=1)
end

# ╔═╡ 440c2010-df04-11ea-2550-6f679efd6b07
md"""# La Salle y región de atracción
Partimos del ejemplo

$$\dot x_1=(x_2-x_1)(1+x_1^2)^2$$

$$\dot x_2=-x_1+x_1^2$$

"""

# ╔═╡ c3aa90e0-df04-11ea-0736-5176e619708a
function f_ejemplo1(x,p,t) 
	#p y t no se usan pero es para usar la función flechas del tema anterior
	x1,x2=x
	dx1=(x2-x1)*(1+x1^2)^2
    dx2=-x1+x1^2
	return [dx1;dx2]
end

# ╔═╡ 116e9d42-812f-11eb-327c-e107809f9603
md"## vamos a empezar a hacer cosas simbólicas (ahora que ya sabemos hacerlas a mano)"

# ╔═╡ d4547220-8121-11eb-1591-ebb8d7742fcc
#como usamos symengine definimos esto a mano, en próximos temas usaremos symbolics
function jacobiano(f,variables)
	N=length(variables)
	fun=f(variables,[],0)
	M=length(fun)
	J=[diff(fun[i],variables[j]) for i in 1:M, j in 1:N]
end

# ╔═╡ 5eebc400-8123-11eb-029e-8bd1bef4612e
x_1,x_2=symbols("x1,x2") #permite operar simbólicamente

# ╔═╡ f9dcadc0-812e-11eb-005f-d751ee820296
J=jacobiano(f_ejemplo1,[x_1,x_2])

# ╔═╡ f13ed140-81bd-11eb-200f-9fc7aa549363
f_ejemplo1([x_1,x_2],[],0)

# ╔═╡ 6eb00ba0-812b-11eb-06b8-9d782689846a
function linealizar(f,variables,puntos)
	#evalua el jacobiano en un punto
	J=jacobiano(f,variables)
	for i in 1:length(variables)
		for j in 1:length(J)
		   J[j]=subs(J[j],variables[i],puntos[i])
		end
	end
	return J
end

# ╔═╡ 83941610-df41-11ea-0f0b-636bcab6c386
A=linealizar(f_ejemplo1,[x_1,x_2],[0,0])

# ╔═╡ c4a86fb0-df42-11ea-2ee5-63a5bdd03376
begin
	autovalores=eigvals(Float64.(A))
	md"Los autovalores son $autovalores"
end

# ╔═╡ 3cbc07a0-df43-11ea-0df8-5ba5311f2f43
md"""Sabemos que es estable (**localmente**)
para saber más derivamos $V_2$ del ejemplo anterior:


$$V_2=\frac{x_1^2}{1+x_1^2}+x_2^2$$

$$\dot V_2=\frac{2x_1}{(1+x_1^2)^2}\dot x_1+2x_2\dot x_2


=2x_1(x_2-x_1)+2x_2(-x_1+x_1^2)$$

$$=-2x_1^2(1-x_2)$$

De momento a mano en el tema 4 veremos cómo decirle a la máquina que lo haga

"""

# ╔═╡ c82134f0-8123-11eb-279a-2b03ca565783
md"V no aumenta luego es **Estable**, pero ¿Podemos saber algo más?"

# ╔═╡ a7b42e60-dfcb-11ea-2ae9-87311393dbef
md"""# Usemos LaSalle para ver la región de atracción
## Paso 1 Elegir L
Si elijo bien L V no aumenta en la región V(x)<L (las flechas no apuntan hacia fuera)"""

# ╔═╡ dd573c40-dfc3-11ea-3ec8-7d0b939123c1
dV2(x1,x2)=-2x1^2*(1-x2)

# ╔═╡ 9b580610-dfc5-11ea-2cfb-bb5db8322e86
@bind L3 Slider(0:0.05:2)

# ╔═╡ 5b862210-dfc5-11ea-3958-d5f618a37108
begin
	x3=-2:0.01:2
	y3=x3
	gr()
	contour(x3, y3, dV2, fill = false, linestyle=:dot)
	contour!(x3, y3, V2, fill = false, colorbar=false)
	contour!(x3,y3,V2,levels=[L3],linewidth=3,color=:red)
	contour!(x3,y3,dV2,levels=[0],linewidth=3,color=:green)
	plot!(title="L=$(L3)")
end

# ╔═╡ d64c1460-f9d2-11ea-307b-c7d59f18b084
md"""**Ojo no siempre se detectan bien las curvas de nivel**

¿Echais algo en falta?"""

# ╔═╡ 94fc4d7e-f9d3-11ea-0b41-a1d22cd9376b
begin
	plotly()
	surface(x3,y3,dV2)
end

# ╔═╡ 698f68b0-dfcc-11ea-2a38-1fd314c52157
md"""Si elijo L=1 V no aumenta en V(x)<L.

Luego si empiezo en $V(x_0)<1$ entonces $V(x(t))<1$ para todo $t$
estamos atrapados en V<L.

V No aumenta y de hecho disminuye en algunos puntos ¿Dónde deja de disminuir? donde $$\dot V=0)$$"""

# ╔═╡ f2777960-dfcc-11ea-089a-af9cb98c26e7
md"""## Paso 2 calcular la región R
R es simplemente la región en la que $V<L$ y además $\dot V=0$
es decir:
$$\dot V_1=-2x_1^2(1-x_2)=0$$
de modo que:

A) $x_1=0$

B) $1-x_2=0$ esta **NO** puede ser ya que si $x_2=1$ entonces $V=\frac{x_1^2}{1+x_1^2}+1^2 \ge 1$ que está **fuera** de $V<1$

Sabemos entonces que $x_1 \to 0$
"""

# ╔═╡ 9370bbb0-dfcd-11ea-3079-d3d74eee20ea
md"""## Paso 3 calcular la región M
M es el conjunto de todas las solucines que viven en R.
Es decir tenemos que hacer $x_1(t)=0$ (y por tanto $\dot x_1=0$ )y ver que pasa al meterlo en las ecuaciones:

$$\dot x_1=$$

$$0=(x_2-0)(1+0^2)^2 \to x_2=0$$

$$\dot x_2=-0+0^2$$

Solo hay un punto en el que las soluciones *se quedan* en la región R que es (0,0) las otras soluciones "pasan de largo" (ver el dibujo)

"""

# ╔═╡ f9362140-dfcf-11ea-2882-77e3c638d875
md"""## Hemos demostrado que si empezamos en $V<1$ entonces $x \to 0$

O lo que es lo mismo que $V<1$ es una región de atración del origen"""

# ╔═╡ 503ad610-dfd6-11ea-1021-5fb2a45bc262
md"""# ¿Por qué ese V?
Podríamos haber usado $V_1$, veamos que sale:
$$V_1=x_1^2+x_2^2$$

$$\dot V_1=2x_1\dot x_1+2x_2\dot x_2


=2x_1(x_2-x_1)(1+x_1^2)^2+2x_2(-x_1+x_1^2)$$


V no aumenta luego es **Estable**, pero ¿Podemos saber algo más?

"""

# ╔═╡ bf615090-dfd7-11ea-12bb-11d922c37139
dV1(x1,x2)=2x1*(x2-x1)*(1+x1^2)^2+2x2*(-x1+x1^2)

# ╔═╡ f0dd5060-dfd7-11ea-03d8-63269e0c3394
@bind L4 Slider(0:0.05:2)

# ╔═╡ a4399020-dfd7-11ea-2e36-ddc90140e91d
begin
	x4=-2:0.01:2
	y4=x3
	gr()
	contour(x4, y4, dV1, fill = false,colorbar=false)
	contour!(x4, y4, V1, fill = false)
	contour!(x4,y4,V1,levels=[L4],linewidth=3,color=:red)
	contour!(x4,y4,dV1,levels=[0],linewidth=3,color=:green)
	plot!(title="L=$(L4)",ratio=1)
end

# ╔═╡ 0edbfe90-dfd8-11ea-038a-97c39cddb672
md"""con L=0.6, V no aumenta

R es la región en la que \dot v=0 que al elegir el l adecuado se reduce a $x_1=0$ el resto es igual que el caso anterior.
"""

# ╔═╡ a96a7e30-d639-11ea-08df-5d727910849c
function flechas!(figura,f,p,rangos;N=10)
	xmin,xmax,ymin,ymax=rangos 
	longitud=max(abs(xmax-xmin),abs(ymax-ymin))/N
	
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax

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
	quiver!(figura,X1,X2,quiver=(U,V), arrow=arrow(0.1*longitud, 0.1*longitud))
	return figura
end

# ╔═╡ 93e7ebb0-dfc6-11ea-26f1-f7e281edf001
begin
	gr()
	p1=plot()
	flechas!(p1,f_ejemplo1,[],[-2,2,-2,2],N=20)
	contour!(x3,y3,V2,levels=[L3],linewidth=3,color=:red,colorbar=false)
	contour!(x3,y3,dV2,levels=[0],linewidth=3,color=:green)
end

# ╔═╡ 8f8b3390-8130-11eb-38df-99e46c4051e0
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

# ╔═╡ e4e5a27e-8130-11eb-0c8a-91d509d7e897
function trayectorias!(figura,f,p,rangos;N=10)
	
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

# ╔═╡ 42ab643e-dfd8-11ea-0ed7-392d94c7904e
begin
	gr()
	p=plot()
	flechas!(p,f_ejemplo1,[],[-2,2,-2,2],N=20)
	trayectorias!(p,f_ejemplo1,[],[-2,2,-2,2],N=30)
	contour!(x3,y3,V2,levels=[L3],linewidth=3,color=:red,colorbar=false)
	contour!(x3,y3,V1,levels=[L4],linewidth=3,color=:red)
	plot!(title="Comparación de regiones",ratio=1)
end

# ╔═╡ Cell order:
# ╠═02c15be0-defd-11ea-004b-45c4e09aa541
# ╟─28fccf20-d63c-11ea-3311-411f89a63431
# ╟─309959f0-defd-11ea-34c3-71eeda7b523b
# ╠═cd88fade-defd-11ea-2097-25b2d1eac285
# ╠═ce5fe800-deff-11ea-0b19-7fc23b14bd06
# ╠═c6c49460-deff-11ea-1c6d-673d863ea07b
# ╠═5b26b050-d7eb-11ea-2488-b7da296fed24
# ╠═2e580002-defe-11ea-1b3d-d5ec461e0715
# ╟─2ec228f0-defd-11ea-329f-53576ab14004
# ╟─ca2fea72-3655-11eb-0984-a532e1f9efcc
# ╠═204e6710-defe-11ea-1fb0-43c7ef863786
# ╠═f3bcc450-df00-11ea-0080-03b54212b9ed
# ╠═ca99149e-df02-11ea-39a0-179e6a7d1344
# ╠═eb00be70-df00-11ea-0ce9-3b0a79a0f4a3
# ╟─d8678810-defc-11ea-05b2-51373f754e61
# ╠═34a9ca60-df03-11ea-018c-399409d7842e
# ╠═f94dece0-defc-11ea-0546-39a094bbd2b6
# ╟─440c2010-df04-11ea-2550-6f679efd6b07
# ╠═c3aa90e0-df04-11ea-0736-5176e619708a
# ╟─116e9d42-812f-11eb-327c-e107809f9603
# ╠═d4547220-8121-11eb-1591-ebb8d7742fcc
# ╠═5eebc400-8123-11eb-029e-8bd1bef4612e
# ╠═f9dcadc0-812e-11eb-005f-d751ee820296
# ╠═f13ed140-81bd-11eb-200f-9fc7aa549363
# ╠═6eb00ba0-812b-11eb-06b8-9d782689846a
# ╠═83941610-df41-11ea-0f0b-636bcab6c386
# ╠═6cbb32f0-d665-11ea-3f24-41cbe1910b80
# ╠═c4a86fb0-df42-11ea-2ee5-63a5bdd03376
# ╟─3cbc07a0-df43-11ea-0df8-5ba5311f2f43
# ╟─c82134f0-8123-11eb-279a-2b03ca565783
# ╟─a7b42e60-dfcb-11ea-2ae9-87311393dbef
# ╠═dd573c40-dfc3-11ea-3ec8-7d0b939123c1
# ╠═5b862210-dfc5-11ea-3958-d5f618a37108
# ╠═9b580610-dfc5-11ea-2cfb-bb5db8322e86
# ╟─d64c1460-f9d2-11ea-307b-c7d59f18b084
# ╠═94fc4d7e-f9d3-11ea-0b41-a1d22cd9376b
# ╠═93e7ebb0-dfc6-11ea-26f1-f7e281edf001
# ╟─698f68b0-dfcc-11ea-2a38-1fd314c52157
# ╟─f2777960-dfcc-11ea-089a-af9cb98c26e7
# ╟─9370bbb0-dfcd-11ea-3079-d3d74eee20ea
# ╟─f9362140-dfcf-11ea-2882-77e3c638d875
# ╟─503ad610-dfd6-11ea-1021-5fb2a45bc262
# ╠═bf615090-dfd7-11ea-12bb-11d922c37139
# ╠═f0dd5060-dfd7-11ea-03d8-63269e0c3394
# ╠═a4399020-dfd7-11ea-2e36-ddc90140e91d
# ╟─0edbfe90-dfd8-11ea-038a-97c39cddb672
# ╟─a96a7e30-d639-11ea-08df-5d727910849c
# ╟─8f8b3390-8130-11eb-38df-99e46c4051e0
# ╟─e4e5a27e-8130-11eb-0c8a-91d509d7e897
# ╠═42ab643e-dfd8-11ea-0ed7-392d94c7904e
