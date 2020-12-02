### A Pluto.jl notebook ###
# v0.12.10

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

# ╔═╡ c09c0e20-d8be-11ea-11be-0966527e6f2d
using Plots

# ╔═╡ 889652b2-dbe8-11ea-0835-e966a39af471
md"""# Problema de estabilidad absoluta de Lur`e
## Partimos de una planta con función de transferencia
"""

# ╔═╡ bf0a5c82-d8bc-11ea-3b1c-b397e3b526fa
H(s)=(s+3)/(s^2+7s+10)

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
		w=vcat(wn,wp)#todas las frecuencias incluidas las negativas
	    Hs=H.(w*im)
	    plot(real.(Hs),imag.(Hs), line = :arrow,label=false)
	#para que el movimiento sea intuitivo el slider ha de ir desde la derecha 
	
	#pinta el punto
		scatter!([-1/K],[0],label=false)
		title!("""La respuesta en frecuencia ha de rodear al punto 
		en sentido antihorario tantas veces como polos inestables hay
		K=$(K)""")
end

# ╔═╡ 7235ff00-dbe0-11ea-3285-37398b545fc4
@bind sliderAK html""" <input type=range min=-1 max=1 step=0.01>
<style type="text/css">
input[type=range] {
  width: 75%;
} </style>
"""

# ╔═╡ 05d824f0-dbe0-11ea-3cd7-a7a012a706f3
K_max_conjeturas=10

# ╔═╡ 235f7e60-dcdf-11ea-2af4-73fc564d9c84
#plotly()

# ╔═╡ c87918e0-dbe3-11ea-1d12-fb31f8ea5bfb
conjeturas(H,sliderAK*K_max_conjeturas)

# ╔═╡ 66f0b770-d8c5-11ea-0ee3-7ba9b6dfb585
md"""# Criterio del circulo
## Antes de empezar comprueba que:
+ Como es una funcion de transferencia ya es contolable/observable (OK)
+ No hay integradores ni realimentación directa $$d=0$$
¿Es estable? Usa los casos A y B. Si no lo es usa los casos C y D.

Cada caso de los que se puede usar da un rango posible de ganancias
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
	#respuesta en frecuencia
	    wp=10.0.^(-50:0.01:50) #rango de 10^-50 a 10ymin=^50
	    wn=-10.0.^(50:-0.01:-50) 
		w=vcat(wn,wp)#todas las frecuencias incluidas las negativas
	   
	    Hs=H.(w*im)
	     plot(Hs,line=:arrow,label=false)
	   Hs=collect(Hs)
	    plot(real.(Hs),imag.(Hs), line=:arrow,label=false)
	
	#pinta el circulo
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
@bind sliderA html"<input type=range min=0 max=1 step=0.01>"

# ╔═╡ b7b101d0-db57-11ea-3a04-c7d11f15bf0c
k1a=K_maxA*sliderA

# ╔═╡ b642eb40-d988-11ea-26a4-374a74d5dc08
circulo(H,0,k1a,'A')

# ╔═╡ b47374ee-db55-11ea-2b8f-5b38e3a2c5a0
md"""### Caso B
$$k_1<0<k_2$$ 
El círculo ha de rodear la respuesta en frecuencia
"""

# ╔═╡ bd8ef6f0-dbb8-11ea-3fe0-3b17ccdb1fcd
begin
	K_maxB1=100
	K_maxB2=100
end

# ╔═╡ fe58cf80-dbb8-11ea-1c6c-93fef94ed79c
@bind sliderB1 html"<input type=range min=-1 max=0 step=0.01>"

# ╔═╡ d2de7080-dbb8-11ea-1d7e-11b40cd2ce8d
@bind sliderB2 html"<input type=range min=0 max=1 step=0.01>"

# ╔═╡ ef6b8d50-dbb8-11ea-3d5d-551e020edd24
begin
	k1b=K_maxB1*sliderB1
	k2b=K_maxB2*sliderB2
	md"""k1=$(k1b),
	
	k2=$(k2b)"""
end

# ╔═╡ 9d833cf0-db03-11ea-3b28-33343306f4dc
circulo(H,k1b,k2b,'B')

# ╔═╡ 1ae51110-d98c-11ea-3b02-1ddccacadb56
md"## Planta inestable"

# ╔═╡ f9938e60-d98b-11ea-1e6d-2fc1399f7bef
md"""### Caso C
$$0<k_1 \le k_2$$ 

Como es inestable habrá $n$ polos en el semiplano derecho.
El círculo ha de **estar rodeado** por la respuesta en frecuencia $n$ veces en sentido antihorario.
"""

# ╔═╡ 1be6ec20-dbba-11ea-2479-ed82409ebdd7
begin
	K_maxC1=50
	K_maxC2=100
end

# ╔═╡ 76a9f030-dbba-11ea-3326-95356d82a2c5
@bind sliderC1 html"<input type=range min=0 max=1 step=0.01>"

# ╔═╡ 783969d0-dbba-11ea-10d0-b1cea6c3b727
@bind sliderC2 html"<input type=range min=0 max=1 step=0.01>"

# ╔═╡ 3fa5ec60-dbba-11ea-0222-054f30eec53d
begin
	k1c=K_maxC1*sliderC1
	k2c=K_maxC2*sliderC2
	md"""k1=$(k1c),
	
	k2=$(k2c)"""
end

# ╔═╡ e3f469a0-dbb9-11ea-30bb-813e5688c466
circulo(H,k1c,k2c,'C')

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
@bind sliderD1 html"<input type=range min=0 max=1 step=0.01>"

# ╔═╡ 8f55db30-dbba-11ea-0370-ddadc3b11dfe
@bind sliderD2 html"<input type=range min=0 max=1 step=0.01>"

# ╔═╡ 4210f8f0-dbba-11ea-071d-9743b745ebbe
begin
	k1d=-K_maxD1*sliderD1
	k2d=-K_maxD2*sliderD2
	md"""k1=$(k1d),
	
	k2=$(k2d)"""
end

# ╔═╡ e660d5c0-dbb9-11ea-2b91-5d86b094cb03
circulo(H,k1d,k2d,'D')

# ╔═╡ 8e2736c0-d8bb-11ea-2fde-43a428dd9281
md"""# Criterio de popov
## Antes de empezar comprueba que:
+ Como es una funcion de transferencia ya es contolable/observable
+ Hay a lo sumo un integrador y si lo hay la planta se puede poner como $$d/s+ H_{estable}$$ con d>0 es decir
    + La parte del ingregrador es positiva ($$H(s) \to \infty$$  si  $$s\to 0^+$$) 
    + Quitando el integrador la planta es estable"""

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
	p=plot(Ws, line = :arrow,label=false)
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
K_max=100

# ╔═╡ cb075200-d8bb-11ea-234b-6d4b6fca9083
@bind slider html"<input type=range min=0 max=1 step=0.01>"

# ╔═╡ f2567110-dbe0-11ea-24a0-4d1a72151516
@bind pendiente html"<input type=range min=0 max=90 step=1>"

# ╔═╡ f8ad6aa0-d8bb-11ea-0dea-85d5051066fe
K=K_max*slider

# ╔═╡ e443a3e0-d8c0-11ea-362f-75314f8b2cbb
alpha=tan(pendiente*pi/180)

# ╔═╡ dbcceee0-d985-11ea-2a3e-f9c3ec08cb75
popov(H,K,alpha)

# ╔═╡ Cell order:
# ╟─889652b2-dbe8-11ea-0835-e966a39af471
# ╠═bf0a5c82-d8bc-11ea-3b1c-b397e3b526fa
# ╟─9999f2c0-dbdd-11ea-04a5-dd56cf1eae5c
# ╠═4628f590-dbe3-11ea-2a7a-83cabda49ee7
# ╟─7235ff00-dbe0-11ea-3285-37398b545fc4
# ╠═05d824f0-dbe0-11ea-3cd7-a7a012a706f3
# ╠═c09c0e20-d8be-11ea-11be-0966527e6f2d
# ╠═235f7e60-dcdf-11ea-2af4-73fc564d9c84
# ╠═c87918e0-dbe3-11ea-1d12-fb31f8ea5bfb
# ╟─66f0b770-d8c5-11ea-0ee3-7ba9b6dfb585
# ╠═9ca66610-d985-11ea-0426-13ef5f785091
# ╟─00b900b0-d8c6-11ea-0e42-0f4bfaf17fd2
# ╟─ed754f20-db56-11ea-02fa-3581c244dcbd
# ╠═9975e6e0-db57-11ea-30e1-3ffbd104839c
# ╟─adee0ba0-d988-11ea-19c7-e324f1f86031
# ╟─b7b101d0-db57-11ea-3a04-c7d11f15bf0c
# ╠═b642eb40-d988-11ea-26a4-374a74d5dc08
# ╟─b47374ee-db55-11ea-2b8f-5b38e3a2c5a0
# ╠═bd8ef6f0-dbb8-11ea-3fe0-3b17ccdb1fcd
# ╟─fe58cf80-dbb8-11ea-1c6c-93fef94ed79c
# ╠═d2de7080-dbb8-11ea-1d7e-11b40cd2ce8d
# ╟─ef6b8d50-dbb8-11ea-3d5d-551e020edd24
# ╠═9d833cf0-db03-11ea-3b28-33343306f4dc
# ╟─1ae51110-d98c-11ea-3b02-1ddccacadb56
# ╟─f9938e60-d98b-11ea-1e6d-2fc1399f7bef
# ╠═1be6ec20-dbba-11ea-2479-ed82409ebdd7
# ╟─76a9f030-dbba-11ea-3326-95356d82a2c5
# ╟─783969d0-dbba-11ea-10d0-b1cea6c3b727
# ╟─3fa5ec60-dbba-11ea-0222-054f30eec53d
# ╠═e3f469a0-dbb9-11ea-30bb-813e5688c466
# ╟─c45e9020-d989-11ea-028f-5f36703d75da
# ╠═245ff590-dbba-11ea-1deb-1976fc8d165d
# ╟─8d0b4eee-dbba-11ea-2fad-5d7ed6f8035b
# ╟─8f55db30-dbba-11ea-0370-ddadc3b11dfe
# ╟─4210f8f0-dbba-11ea-071d-9743b745ebbe
# ╠═e660d5c0-dbb9-11ea-2b91-5d86b094cb03
# ╟─8e2736c0-d8bb-11ea-2fde-43a428dd9281
# ╟─497f9250-d8bc-11ea-194d-77979b94814d
# ╠═fee20250-d8bb-11ea-299c-6d093b0edd20
# ╟─cb075200-d8bb-11ea-234b-6d4b6fca9083
# ╟─f2567110-dbe0-11ea-24a0-4d1a72151516
# ╟─f8ad6aa0-d8bb-11ea-0dea-85d5051066fe
# ╟─e443a3e0-d8c0-11ea-362f-75314f8b2cbb
# ╠═dbcceee0-d985-11ea-2a3e-f9c3ec08cb75
