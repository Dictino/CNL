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

# ╔═╡ 8659d240-0be2-11eb-23bf-a3d1f4c1b2fc
using Plots

# ╔═╡ c5f84a30-0be2-11eb-1678-d3e648fb3856
using PlutoUI

# ╔═╡ 01f3a030-0be2-11eb-2cee-85e45d0130a2
md"""# El camino hacia el caos
Presentado por el Dr. Chaos ☺

Ecuación logística

$y_{n+1}=ry_n(1-y_n)$
$r \in [0, 4]$
$y \in [0, 1]$

"""

# ╔═╡ 264ff0e0-0bed-11eb-38f3-912127308019
logistica(y,r)=r*y*(1-y)

# ╔═╡ ae42b0a0-0be3-11eb-1055-e13728ef27b3
function trayectoria(y_inicial,r,pasos)
	y=zeros(pasos)
	y[1]=y_inicial
	for i=1:(pasos-1)
		y[i+1]=logistica(y[i],r)
	end
	return y
end

# ╔═╡ cdfb1230-0be2-11eb-0d7a-19d47f6c077e
@bind r Slider(0:0.05:4,default=2)

# ╔═╡ b40bd1c0-0be2-11eb-2659-ab71da779ed4
begin
	x=0:0.01:1
	fig=plot(x,logistica.(x,r),ylim=(0,1),label=:none, title="r=$r")
end

# ╔═╡ 3e8d8d10-0c16-11eb-3e59-6d658fca5049
md"Ver rastro"

# ╔═╡ 72621fae-0c12-11eb-2078-9d12085b8980
begin
	@bind ver CheckBox(default=false)
end

# ╔═╡ aadc88f0-0be8-11eb-1a60-33a305646fb0
function puntos(R)
	px=[]
	py=[]
	for r=0:0.02:R
		yt=trayectoria(0.5,r,100)
		k=length(yt)
		y_est=yt[(k÷2):k]
		append!(px,r*ones(length(y_est)))
		append!(py,y_est)
	end
	return (px,py)
end

# ╔═╡ 6d4ca280-0be4-11eb-00b6-1f5ce746261f
let
	y=trayectoria(0.5,r,50)
	p1=plot(y,title="r=$r",marker=:dot,markersize=:2,label=:none,ylimits=(0,1))
	k=length(y)
	px,py=puntos(4)
	p2=scatter()
	if ver
		p2=scatter!(p2,px,py,legend=:none,color=:grey, alpha=0.1,
		            linecolor=:grey, linealpha=0.1,marker=1)
	end
	y_estacionarios=y[(k÷2):k] #Nota ver cómo se escribe esto
	p2=scatter!(p2,r*ones(length(y_estacionarios)),y_estacionarios,
		        xlim=(0,4),ylim=(0,1),legend=:none,color=:red,marker=2)
	l = @layout [a ; b]
	plot(p1,p2, layout=l, size=(650,500))
end

# ╔═╡ d8c3acc0-0be9-11eb-3988-7138e6bd331a
md"""# ¿Cómo se llega al caos?

Haciendo lo mismo una y otra vez...

Primer paso, periodo 1:"""

# ╔═╡ 9c1d3b20-0c15-11eb-1eea-158897d1327f
@bind y0 Slider(0:0.05:1,default=0.5)

# ╔═╡ 7f5d8020-0c16-11eb-0305-01c68de19429
@bind N Slider(1:50,default=1)

# ╔═╡ 8d99de50-0c15-11eb-3d1d-1747eb4eb117
@bind r2 Slider(0:0.05:4,default=2)

# ╔═╡ dd67b2a0-0e08-11eb-2f9e-bb83de99d049
@bind n Slider(2:2:20,default=2)

# ╔═╡ afe95392-0c16-11eb-10dd-b7d54a24c445
begin
	f(x)=logistica(x,r2)
	f2(x)=f(f(x))
	f3(x)=f(f2(x))
	f4(x)=f(f3(x))
	function fn(x)
		y=x
		for i=1:n
			y=f(y)
		end
		return y
	end
end

# ╔═╡ 771e8900-0e04-11eb-3b0f-3552d2de3910
md"""En un sistema discreto la estabilidad depende de la derivada, es parecido a los sistemas continuos.
Un punto de equilibiro $x^*$ es aquel en el que $x_{n+1}=x_n$, es decir

$x_{n+1}=f(x_n)=x_n=x^*$

Es decir que es un punto fijo de $f$, es decir $f(x^*)=x^*$

Si estámos cerca de un punto de equilibrio $x_n=x^*+\delta x_n$ luego haciendo taylor:

$f(x^*+\delta x_n)=f(x^*) + \nabla f(x^*)\delta x_n + o(\delta x_n^2)$


Pero entonces:

$\delta x_{n+1}=x_{n+1}-x^*=f(x^*)-x^* + \nabla f(x^*)\delta x_n + o(\delta x_n^2) \approx \nabla f(x^*) \delta x_n=A \delta x_n$

Cuya solución es $\delta x_n=A^n\delta x_0$ que es estable si los autovalores de $A$ están en el círculo unitario. Es parecido a los sistemas continuos cambiando el eje imaginario por el **círculo unidad**.

En el caso de una variable se reduce a que $|f'(x^*)|<1 como se intuye gráficamente$

"""

# ╔═╡ 754e14b0-0c15-11eb-1230-01eaa309edac
function iterar(f,N)
	x=0:0.001:1
	fig2=plot(x,f.(x),ylim=(0,1),label=:none, title="r=$r2")
	plot!(fig2,[0,1],[0,1],color=:red,label=:none)
	
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

# ╔═╡ 8c0c03a0-0c1b-11eb-35a5-37dbeb561d42
iterar(f,N) # f1 f2 y luego pasarse a fn

# ╔═╡ 20039fc0-0c1f-11eb-3007-7dbf66a6c5b5
md"En cada paso el cambio de $r$ necesario para producir una bifurcación es más pequeño.

De hecho la bifurcaciñon de orden 4 es como una versión pequeña de la de orden 2 como hemos visto, si llamamos $r_n$ al $r$ que produce una bifurcación tenemos que

$\frac{r_{n+2}-r_{n+1}}{r_{n+2}-r_{n+1}} \to 4.6992...$ que es la constante de Feigenbaum y que es **la misma** para muchos otros procesos esto se conoce como **Universalidad**.

Como el paso sigue aproximadamente una serie geométrica la suma de infintas bifurcaciones se produce con un cambio finito de $r$ ya que:

$r_\infty-r_0 \approx \Delta r_0 +  \frac {\Delta r_0}{4.7} +\frac {\Delta r_0}{4.7^2}+\frac {\Delta r_0}{4.7^3}...=\Delta r_0 \frac {1}{1-1/4.7}\approx 1.27\Delta r_0$

" 

# ╔═╡ Cell order:
# ╟─01f3a030-0be2-11eb-2cee-85e45d0130a2
# ╠═8659d240-0be2-11eb-23bf-a3d1f4c1b2fc
# ╠═c5f84a30-0be2-11eb-1678-d3e648fb3856
# ╠═264ff0e0-0bed-11eb-38f3-912127308019
# ╠═b40bd1c0-0be2-11eb-2659-ab71da779ed4
# ╠═ae42b0a0-0be3-11eb-1055-e13728ef27b3
# ╠═cdfb1230-0be2-11eb-0d7a-19d47f6c077e
# ╟─3e8d8d10-0c16-11eb-3e59-6d658fca5049
# ╟─72621fae-0c12-11eb-2078-9d12085b8980
# ╠═6d4ca280-0be4-11eb-00b6-1f5ce746261f
# ╠═aadc88f0-0be8-11eb-1a60-33a305646fb0
# ╟─d8c3acc0-0be9-11eb-3988-7138e6bd331a
# ╠═9c1d3b20-0c15-11eb-1eea-158897d1327f
# ╠═7f5d8020-0c16-11eb-0305-01c68de19429
# ╠═8d99de50-0c15-11eb-3d1d-1747eb4eb117
# ╠═8c0c03a0-0c1b-11eb-35a5-37dbeb561d42
# ╟─afe95392-0c16-11eb-10dd-b7d54a24c445
# ╠═dd67b2a0-0e08-11eb-2f9e-bb83de99d049
# ╟─771e8900-0e04-11eb-3b0f-3552d2de3910
# ╠═754e14b0-0c15-11eb-1230-01eaa309edac
# ╟─20039fc0-0c1f-11eb-3007-7dbf66a6c5b5
