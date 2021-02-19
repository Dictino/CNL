### A Pluto.jl notebook ###
# v0.12.20

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

# ╔═╡ 02c15be0-defd-11ea-004b-45c4e09aa541
using Plots

# ╔═╡ 50cdc34e-df3e-11ea-0965-793860f72c45
using ModelingToolkit

# ╔═╡ 6cbb32f0-d665-11ea-3f24-41cbe1910b80
using LinearAlgebra

# ╔═╡ a96a7e30-d639-11ea-08df-5d727910849c
function flechas(f,p,rangos;N=10)
	#version no modificante del tema anterior
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
			X=f([x1,x2],p,0)
			X=0.1*X/(0.1+sqrt(X[1]^2+X[2]^2))
			push!(U,X[1]);
			push!(V,X[2]);
		end
	end
	p1=quiver(X1,X2,quiver=(U,V), arrow=arrow(0.1*long, 0.1*long))
	return p1
end

# ╔═╡ 28fccf20-d63c-11ea-3311-411f89a63431
md"# Curvas de nivel"

# ╔═╡ 309959f0-defd-11ea-34c3-71eeda7b523b
md"Curvas cerradas:

$$V(x)\to \infty$$ cuando $$x\to \infty$$"

# ╔═╡ cd88fade-defd-11ea-2097-25b2d1eac285
V1(x1,x2)=x1^2+x2^2

# ╔═╡ ce5fe800-deff-11ea-0b19-7fc23b14bd06
begin
		x=-1:0.1:1
		y=x;
end

# ╔═╡ c6c49460-deff-11ea-1c6d-673d863ea07b
begin
	plotly()
	surface(x,y,V1)
end

# ╔═╡ 5b26b050-d7eb-11ea-2488-b7da296fed24
@bind L1 html"<input type=range min=0 max=2 step=0.05>"

# ╔═╡ 2e580002-defe-11ea-1b3d-d5ec461e0715
begin
	gr()
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

# ╔═╡ ca99149e-df02-11ea-39a0-179e6a7d1344
@bind L2 html"<input type=range min=0 max=2 step=0.05>"

# ╔═╡ fe41070e-df00-11ea-3be2-b72adee0837d
begin
		x2=-5:0.1:5
		y2=-1:0.1:1
end

# ╔═╡ f3bcc450-df00-11ea-0080-03b54212b9ed
begin
	plotly()
	surface(x2,y2,V2,ratio=1)
end

# ╔═╡ eb00be70-df00-11ea-0ce9-3b0a79a0f4a3
begin
	gr()
	contour(x2, y2, V2, fill = true)
	contour!(x2,y2,V2,levels=[L2],linewidth=3,color=:red) #, ratio=1)
end

# ╔═╡ d8678810-defc-11ea-05b2-51373f754e61
md"# Semi-Definida Vs Definida"

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
	contour!(x,y,f,levels=[0],linewidth=3,color=:red, ratio=1)
end

# ╔═╡ 440c2010-df04-11ea-2550-6f679efd6b07
md"""# La Salle y región de atracción
Partimos de ejemplo

$$\dot x_1=(x_2-x_1)(1+x_1^2)^2$$

$$\dot x_2=-x_1+x_1^2$$

"""

# ╔═╡ e062be10-df40-11ea-1565-3d684f59321c
function linealizacion(f,x0,p,t)
	n=length(x0)
	@variables X[1:n]
	op=f(X,p,t)
	jac=ModelingToolkit.jacobian(op,X)
	for i in 1:n
	    jac=substitute.(jac,[ X[i]=>x0[i] ])
	end
	map(x->x.value,jac)
end

# ╔═╡ c3aa90e0-df04-11ea-0736-5176e619708a
function f_ejemplo1(x,p,t)
	dx1=(+x[2]-x[1])*(1+x[1]^2)^2
    dx2=-x[1]+x[1]^2 #busca otro ejemplo
	return [dx1;dx2]
end

# ╔═╡ 83941610-df41-11ea-0f0b-636bcab6c386
A=linealizacion(f_ejemplo1,[0;0],[],0)

# ╔═╡ c4a86fb0-df42-11ea-2ee5-63a5bdd03376
begin
	autovalores,_=eigen(A)
	md"Los autovalores son $autovalores"
end

# ╔═╡ 3cbc07a0-df43-11ea-0df8-5ba5311f2f43
md"""Sabemos que es estable (localmente)
para saber mas derivamos $V_2$ del ejemplo anterior:


$$V_2=\frac{x_1^2}{1+x_1^2}+x_2^2$$

$$\dot V_2=\frac{2x_1}{(1+x_1^2)^2}\dot x_1+2x_2\dot x_2


=2x_1(x_2-x_1)+2x_2(-x_1+x_1^2)$$

$$=-2x_1^2(1-x_2)$$

V no aumenta luego es **Estable**, pero ¿Podemos saber algo más?
"""

# ╔═╡ a7b42e60-dfcb-11ea-2ae9-87311393dbef
md"""# Usemos LaSalle para ver la región de atracción
## Paso 1 Elegir L
Si elijo bien L V no aumenta en la región V(x)<L (las flechas no apuntan hacia fuera)"""

# ╔═╡ dd573c40-dfc3-11ea-3ec8-7d0b939123c1
dV2(x1,x2)=-2x1^2*(1-x2)

# ╔═╡ 9b580610-dfc5-11ea-2cfb-bb5db8322e86
@bind L3 html"<input type=range min=0 max=2 step=0.05>"

# ╔═╡ 5b862210-dfc5-11ea-3958-d5f618a37108
begin
	x3=-2:0.01:2
	y3=x3
	gr()
	contour(x3, y3, dV2, fill = false, linestyle=:dot)
	contour!(x3, y3, V2, fill = false)
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

# ╔═╡ 93e7ebb0-dfc6-11ea-26f1-f7e281edf001
begin
	gr()
	flechas(f_ejemplo1,[],[-2,2,-2,2];N=20)
	contour!(x3,y3,V2,levels=[L3],linewidth=3,color=:red)
	contour!(x3,y3,dV2,levels=[0],linewidth=3,color=:green)
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
@bind L4 html"<input type=range min=0 max=2 step=0.05>"

# ╔═╡ a4399020-dfd7-11ea-2e36-ddc90140e91d
begin
	x4=-2:0.01:2
	y4=x3
	gr()
	contour(x4, y4, dV1, fill = false)
	contour!(x4, y4, V1, fill = false)
	contour!(x4,y4,V1,levels=[L4],linewidth=3,color=:red)
	contour!(x4,y4,dV1,levels=[0],linewidth=3,color=:green)
	plot!(title="L=$(L4)",ratio=1)
end

# ╔═╡ 0edbfe90-dfd8-11ea-038a-97c39cddb672
md"""con L=0.6, V no aumenta

R es la región en la que \dot v=0 que al elegir el l adecuado se reduce a $x_1=0$ el resto es igual que el caso anterior.
"""

# ╔═╡ 42ab643e-dfd8-11ea-0ed7-392d94c7904e
begin
	gr()
	flechas(f_ejemplo1,[],[-2,2,-2,2],N=20)
	contour!(x3,y3,V2,levels=[L3],linewidth=3,color=:red)
	contour!(x3,y3,V1,levels=[L4],linewidth=3,color=:red)
	plot!(title="Comparación de regiones",ratio=1)
end

# ╔═╡ Cell order:
# ╠═02c15be0-defd-11ea-004b-45c4e09aa541
# ╟─a96a7e30-d639-11ea-08df-5d727910849c
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
# ╠═fe41070e-df00-11ea-3be2-b72adee0837d
# ╠═eb00be70-df00-11ea-0ce9-3b0a79a0f4a3
# ╠═d8678810-defc-11ea-05b2-51373f754e61
# ╠═34a9ca60-df03-11ea-018c-399409d7842e
# ╠═f94dece0-defc-11ea-0546-39a094bbd2b6
# ╟─440c2010-df04-11ea-2550-6f679efd6b07
# ╠═50cdc34e-df3e-11ea-0965-793860f72c45
# ╟─e062be10-df40-11ea-1565-3d684f59321c
# ╠═c3aa90e0-df04-11ea-0736-5176e619708a
# ╠═83941610-df41-11ea-0f0b-636bcab6c386
# ╠═6cbb32f0-d665-11ea-3f24-41cbe1910b80
# ╟─c4a86fb0-df42-11ea-2ee5-63a5bdd03376
# ╟─3cbc07a0-df43-11ea-0df8-5ba5311f2f43
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
# ╠═0edbfe90-dfd8-11ea-038a-97c39cddb672
# ╠═42ab643e-dfd8-11ea-0ed7-392d94c7904e
